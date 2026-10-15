nextflow.enable.types = true

// Flatten the tool subdirectories into the task root so every record field is a plain file,
// keeping qualimap's own directory structure under the sample-named directory (the tool's
// --outdir stays ${prefix} because featureCounts records it in its output header).
def flattenOutputs(prefix: String) {
    """
    dups=\$(find ${prefix} -type f -not -path '${prefix}/qualimap/*' | sed 's#.*/##' | sort | uniq -d)
    if [ -n "\$dups" ]; then echo "RustQC output basenames collide: \$dups" >&2; exit 1; fi
    find ${prefix} -type f -not -path '${prefix}/qualimap/*' -exec mv {} . \\;
    mv ${prefix}/qualimap qualimap_tmp
    rm -rf ${prefix}
    mv qualimap_tmp ${prefix}
    """
}

process RUSTQC {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2a/2a8a0514855c54307399fd0f664c2685e76c8cc07631e767c1e37c575b18d59f/data'
        : 'community.wave.seqera.io/library/rustqc:0.2.1--00df1502b490e005'}"

    input:
    tuple(
        meta: Map,
        bam:  Path,
        bai:  Path,
        gtf:  Path,
        args: String?
    )

    output:
    // (meta, samtools(stats, flagstat, idxstats), preseq_lc_extrap, dupradar(scatter2d, boxplot, hist, dupmatrix, intercept_slope, multiqc), featurecounts(counts, summary), biotype(tsv, mqc, rrna), rseqc(bamstat, inferexperiment, readdistribution, tin(txt, xls), innerdistance(distance, freq, mean, summary, plot, rscript), junctionannotation(bed, interact_bed, xls, log, plot, rscript), junctionsaturation(summary, plot, rscript), readduplication(seq_xls, pos_xls, plot, rscript)), qualimap, all_files)
    tuple(
        meta,
        tuple(file("*.stats", optional: true), file("*.flagstat", optional: true), file("*.idxstats", optional: true)),
        file("*.lc_extrap.txt", optional: true),
        tuple(files("*duprateExpDens.*", optional: true), files("*duprateExpBoxplot.*", optional: true), files("*expressionHist.*", optional: true), file("*dupMatrix.*", optional: true), file("*intercept_slope.*", optional: true), files("*{dup_intercept,duprateExpDensCurve}_mqc.txt", optional: true)),
        tuple(file("*.featureCounts.tsv", optional: true), file("*.featureCounts.biotype.tsv.summary", optional: true)),
        tuple(file("*.biotype_counts.tsv", optional: true), file("*.biotype_counts_mqc.tsv", optional: true), file("*.biotype_counts_rrna_mqc.tsv", optional: true)),
        tuple(file("*.bam_stat.txt", optional: true), file("*.infer_experiment.txt", optional: true), file("*.read_distribution.txt", optional: true), tuple(file("*.summary.txt", optional: true), file("*.tin.xls", optional: true)), tuple(file("*.inner_distance.txt", optional: true), file("*.inner_distance_freq.txt", optional: true), file("*.inner_distance_mean.txt", optional: true), file("*.inner_distance_summary.txt", optional: true), files("*.inner_distance_plot.{png,svg}", optional: true), file("*.inner_distance_plot.r", optional: true)), tuple(file("*.junction.bed", optional: true), file("*.junction.Interact.bed", optional: true), file("*.junction.xls", optional: true), file("*.junction_annotation.log", optional: true), files("*.splice_*.{png,svg}", optional: true), file("*.junction_plot.r", optional: true)), tuple(file("*.junctionSaturation_summary.txt", optional: true), files("*.junctionSaturation_plot.{png,svg}", optional: true), file("*.junctionSaturation_plot.r", optional: true)), tuple(file("*.seq.DupRate.xls", optional: true), file("*.pos.DupRate.xls", optional: true), files("*.DupRate_plot.{png,svg}", optional: true), file("*.DupRate_plot.r", optional: true))),
        file("${prefix}", optional: true),
        files("{*.txt,*.tsv,*.xls,*.log,*.stats,*.flagstat,*.idxstats,*.html,*_mqc.*,${prefix}/**}", optional: true)
    )

    topic:
    tuple(task.process, 'rustqc', eval("rustqc --version 2>&1 | sed -n '1s/rustqc //; 1s/ .*//p'")) >> 'versions'

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tool_args = task.ext.args ?: args ?: ''
    def paired = meta.single_end ? '' : '--paired'
    """
    rustqc rna \\
        ${bam} \\
        --gtf ${gtf} \\
        ${paired} \\
        --threads ${task.cpus} \\
        --outdir ${prefix} \\
        --sample-name ${prefix} \\
        ${tool_args}
    ${flattenOutputs(prefix)}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/{dupradar,featurecounts,preseq,samtools} \\
            ${prefix}/rseqc/{bam_stat,infer_experiment,read_duplication,read_distribution,junction_annotation,junction_saturation,inner_distance,tin} \\
            ${prefix}/qualimap/{raw_data_qualimapReport,images_qualimapReport}

    touch ${prefix}/dupradar/${prefix}_{duprateExpDens,duprateExpBoxplot,expressionHist}.png \\
          ${prefix}/dupradar/${prefix}_{dupMatrix,intercept_slope,dup_intercept_mqc,duprateExpDensCurve_mqc}.txt

    touch ${prefix}/featurecounts/${prefix}.featureCounts.tsv \\
          ${prefix}/featurecounts/${prefix}.featureCounts.tsv.summary \\
          ${prefix}/featurecounts/${prefix}.featureCounts.biotype.tsv.summary \\
          ${prefix}/featurecounts/${prefix}.{biotype_counts,biotype_counts_mqc,biotype_counts_rrna_mqc}.tsv

    touch ${prefix}/rseqc/bam_stat/${prefix}.bam_stat.txt \\
          ${prefix}/rseqc/infer_experiment/${prefix}.infer_experiment.txt \\
          ${prefix}/rseqc/read_duplication/${prefix}.{pos.DupRate,seq.DupRate}.xls \\
          ${prefix}/rseqc/read_duplication/${prefix}.DupRate_plot.{r,png} \\
          ${prefix}/rseqc/read_distribution/${prefix}.read_distribution.txt \\
          ${prefix}/rseqc/junction_annotation/${prefix}.{junction.xls,junction.bed,junction_plot.r,junction_annotation.log,splice_events.png,splice_junction.png} \\
          ${prefix}/rseqc/junction_saturation/${prefix}.junctionSaturation_{plot.r,plot.png,summary.txt} \\
          ${prefix}/rseqc/inner_distance/${prefix}.inner_distance{.txt,_freq.txt,_plot.r,_plot.png,_summary.txt,_mean.txt} \\
          ${prefix}/rseqc/tin/${prefix}.{tin.xls,summary.txt}

    touch ${prefix}/preseq/${prefix}.lc_extrap.txt \\
          ${prefix}/samtools/${prefix}.{flagstat,idxstats,stats} \\
          ${prefix}/qualimap/rnaseq_qc_results.txt \\
          ${prefix}/qualimap/qualimapReport.html
    touch "${prefix}/qualimap/raw_data_qualimapReport/coverage_profile_along_genes_(total).txt"
    ${flattenOutputs(prefix)}
    """
}
