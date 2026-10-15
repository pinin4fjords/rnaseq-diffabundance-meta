nextflow.enable.types = true

process PARABRICKS_RNAFQ2BAM {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'
    // needed by the module to work properly can be removed when fixed upstream - see: https://github.com/nf-core/modules/issues/7226
    stageInMode 'copy'

    container "nvcr.io/nvidia/clara/clara-parabricks:4.7.1-1"

    input:
    tuple(
        meta:            Map,
        reads:           List<Path>,
        fasta:           Path,
        index:           Path,
        qc_metrics:      Boolean,
        mark_duplicates: Boolean,
        args:            String?
    )

    output:
    // (meta, raw_bams, bam_sorted, bam_sorted_aligned, bam_unsorted, transcriptome_bam, unmapped, sam, junction, spl_junc_tab, read_per_gene_tab, wig, bedgraph, orig_bai, qc_metrics, duplicate_metrics, star(log_final, log_out, log_progress, tab))
    tuple(
        meta,
        files("${prefix}.bam", optional: true).toSorted { f -> f.name },
        file("${prefix}.sortedByCoord.out.bam", optional: true),
        file("${prefix}.Aligned.sortedByCoord.out.bam", optional: true),
        file('*Aligned.unsort.out.bam', optional: true),
        file('*toTranscriptome.out.bam', optional: true),
        files('*fastq.gz', optional: true).toSorted { f -> f.name },
        file('*.out.sam', optional: true),
        file('*.out.junction', optional: true),
        file('*.SJ.out.tab', optional: true),
        file('*.ReadsPerGene.out.tab', optional: true),
        files('*.wig', optional: true).toSorted { f -> f.name },
        files('*.bg', optional: true).toSorted { f -> f.name },
        file("${prefix}.bam.bai", optional: true),
        file("${prefix}_qc_metrics", optional: true),
        file("${prefix}.duplicate-metrics.txt", optional: true),
        tuple(file("${prefix}.Log.final.out"), file("${prefix}.Log.out"), file("${prefix}.Log.progress.out"), files('*.tab', optional: true).toSorted { f -> f.name })
    )

    topic:
    tuple(task.process, 'parabricks', eval("pbrun version 2>&1 | grep -Po '(?<=^pbrun: ).*'")) >> 'versions'

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    def tool_args = task.ext.args ?: args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def in_fq_command = meta.single_end ? "--in-se-fq ${reads.join(' ')}" : "--in-fq ${reads.join(' ')}"
    def num_gpus = task.accelerator ? "--num-gpus ${task.accelerator.request}" : ''

    def qc_metrics_command = qc_metrics ? "--out-qc-metrics-dir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_command = mark_duplicates ? "--out-duplicate-metrics ${prefix}.duplicate-metrics.txt" : "--no-markdups"

    """
    pbrun \\
        rna_fq2bam  \\
        --ref ${fasta} \\
        ${in_fq_command} \\
        --output-dir . \\
        --genome-lib-dir ${index} \\
        --out-bam ${prefix}.bam \\
        --logfile ${prefix}.Log.final.out \\
        --out-prefix ${prefix}. \\
        ${num_gpus} \\
        ${qc_metrics_command} \\
        ${duplicate_metrics_command} \\
        ${tool_args}
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Parabricks module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }
    prefix = task.ext.prefix ?: "${meta.id}"
    def qc_metrics_output = qc_metrics ? "mkdir ${prefix}_qc_metrics" : ""
    def duplicate_metrics_output = mark_duplicates ? "touch ${prefix}.duplicate-metrics.txt" : ""
    """
    echo "" | gzip > ${prefix}.unmapped_1.fastq.gz
    echo "" | gzip > ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}.bam
    touch ${prefix}.bam.bai
    touch ${prefix}.Log.final.out
    touch ${prefix}.Log.out
    touch ${prefix}.Log.progress.out
    touch ${prefix}.sortedByCoord.out.bam
    touch ${prefix}.toTranscriptome.out.bam
    touch ${prefix}.Aligned.unsort.out.bam
    touch ${prefix}.Aligned.sortedByCoord.out.bam
    touch ${prefix}.tab
    touch ${prefix}.SJ.out.tab
    touch ${prefix}.ReadsPerGene.out.tab
    touch ${prefix}.Chimeric.out.junction
    touch ${prefix}.out.sam
    touch ${prefix}.Signal.UniqueMultiple.str1.out.wig
    touch ${prefix}.Signal.UniqueMultiple.str1.out.bg
    ${qc_metrics_output}
    ${duplicate_metrics_output}
    """
}
