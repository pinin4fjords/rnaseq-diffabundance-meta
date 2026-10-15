nextflow.enable.types = true

process SENTIEON_STARALIGN {
    tag "$meta.id"
    label 'process_high'
    label 'sentieon'


    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73e9111552beb76e2ad3ad89eb75bed162d7c5b85b2433723ecb4fc96a02674a/data'
        : 'community.wave.seqera.io/library/sentieon:202503.02--def60555294d04fa'}"

    input:
    tuple(
        meta:                Map,
        reads:               List<Path>,
        args:                String?,
        index:               Path,
        gtf:                 Path?,
        star_ignore_sjdbgtf: Boolean
    )

    stage:
    stageAs reads, 'input*/*'

    output:
    // (meta, raw_bams, bam_sorted, bam_sorted_aligned, bam_unsorted, transcriptome_bam, unmapped, sam, junction, spl_junc_tab, read_per_gene_tab, wig, bedgraph, orig_bai, qc_metrics, duplicate_metrics, star(log_final, log_out, log_progress, tab))
    tuple(
        meta,
        files('*d.out.bam', optional: true).toSorted { f -> f.name },
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
        null,
        null,
        null,
        tuple(file('*Log.final.out'), file('*Log.out'), file('*Log.progress.out'), files('*.tab', optional: true).toSorted { f -> f.name })
    )

    topic:
    tuple(task.process, 'star', eval('sentieon STAR --version | sed -e "s/STAR_//g"')) >> 'versions'
    tuple(task.process, 'sentieon', eval('sentieon driver --version 2>&1 | sed -e "s/sentieon-genomics-//g"')) >> 'versions'

    script:
    def tool_args = task.ext.args ?: args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def read_pairs = reads.collate(2)
    def reads1 = meta.single_end ? reads : read_pairs.collect { pair -> pair[0] }.toList()
    def reads2 = meta.single_end ? [] : read_pairs.collect { pair -> pair[1] }.toList()
    def ignore_gtf = star_ignore_sjdbgtf ? '' : "--sjdbGTFfile ${gtf}"
    attrRG = tool_args.contains("--outSAMattrRGline") ? "" : "--outSAMattrRGline 'ID:${prefix}' 'SM:${prefix}'"
    def out_sam_type = tool_args.contains('--outSAMtype') ? '' : '--outSAMtype BAM Unsorted'
    mv_unsorted_bam = tool_args.contains('--outSAMtype BAM Unsorted SortedByCoordinate') ? "mv ${prefix}.Aligned.out.bam ${prefix}.Aligned.unsort.out.bam" : ''

    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    $sentieonLicense

    sentieon STAR \\
        --genomeDir ${index} \\
        --readFilesIn ${reads1.join(",")} ${reads2.join(",")} \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${prefix}. \\
        ${out_sam_type} \\
        ${ignore_gtf} \\
        ${attrRG} \\
        ${tool_args}

    ${mv_unsorted_bam}

    if [ -f ${prefix}.Unmapped.out.mate1 ]; then
        mv ${prefix}.Unmapped.out.mate1 ${prefix}.unmapped_1.fastq
        gzip ${prefix}.unmapped_1.fastq
    fi
    if [ -f ${prefix}.Unmapped.out.mate2 ]; then
        mv ${prefix}.Unmapped.out.mate2 ${prefix}.unmapped_2.fastq
        gzip ${prefix}.unmapped_2.fastq
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.unmapped_1.fastq.gz
    echo "" | gzip > ${prefix}.unmapped_2.fastq.gz
    touch ${prefix}Xd.out.bam
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
    """
}
