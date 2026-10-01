nextflow.enable.types = true

include { params as RnaseqParams ; workflow as NFCORE_RNASEQ } from './pipelines/nf-core/rnaseq'
include { SampleRow                                          } from './pipelines/nf-core/rnaseq/modules/nf-core/types'
include { params as DiffabParams ; workflow as NFCORE_DIFFERENTIALABUNDANCE } from './pipelines/nf-core/differentialabundance'

params {
    samples:         Channel<SampleRow>   // rnaseq samplesheet, one row per sequencing run
    sample_metadata: Path                 // differentialabundance observations: a sample column and the experimental variables
    contrasts:       Path                 // differentialabundance contrasts
    rnaseq:          RnaseqParams         // references and options of nf-core/rnaseq
    diffab:          DiffabParams         // options of nf-core/differentialabundance

    // rnaseq params that its process config reads as top-level params as well as the pipeline itself.
    // Set them here, not in params.rnaseq, so that both see the same value.
    aligner:                     String  = 'star_salmon'
    bam_csi_index:               Boolean = false
    contaminant_screening:       String?
    contaminant_screening_input: String  = 'unmapped'
    featurecounts_group_type:    String  = 'gene_biotype'
    gencode:                     Boolean = false
    prokaryotic:                 Boolean = false
    pseudo_aligner:              String?
    save_unaligned:              Boolean = false
    seq_center:                  String?
    seq_platform:                String?
    with_umi:                    Boolean = false
}

workflow {
    main:
    // Quantify the samples. rnaseq starts differential abundance as soon as the merged matrices exist.
    rnaseq = NFCORE_RNASEQ(
        params.rnaseq + record(
            input:                       params.samples,
            aligner:                     params.aligner,
            bam_csi_index:               params.bam_csi_index,
            contaminant_screening:       params.contaminant_screening,
            contaminant_screening_input: params.contaminant_screening_input,
            featurecounts_group_type:    params.featurecounts_group_type,
            gencode:                     params.gencode,
            prokaryotic:                 params.prokaryotic,
            pseudo_aligner:              params.pseudo_aligner,
            save_unaligned:              params.save_unaligned,
            seq_center:                  params.seq_center,
            seq_platform:                params.seq_platform,
            with_umi:                    params.with_umi
        )
    )

    // differentialabundance takes rnaseq's merged gene-level outputs as its matrix and feature lengths
    def ch_quant = rnaseq.gene_quant.collect().map { quants -> quants.toList().first() }

    abundance = NFCORE_DIFFERENTIALABUNDANCE(
        params.diffab + record(
            input:                 channel.value(params.sample_metadata),
            contrasts:             channel.value(params.contrasts),
            matrix:                ch_quant.map { quant -> quant.counts_gene },
            feature_length_matrix: ch_quant.map { quant -> quant.lengths_gene },
            gtf:                   rnaseq.gtf
        )
    )

    publish:
    multiqc      = rnaseq.multiqc.map { r -> r.report }
    gene_counts  = rnaseq.gene_quant.map { r -> r.counts_gene }
    gene_lengths = rnaseq.gene_quant.map { r -> r.lengths_gene }
    gene_tpm     = rnaseq.gene_quant.map { r -> r.tpm_gene }
    report       = abundance.report.filter { r -> r.name == 'report_html' }.flatMap { r -> r.files }
}

output {
    multiqc: Channel<Path> {
        path 'rnaseq/multiqc'
    }
    gene_counts: Channel<Path> {
        path 'rnaseq/quant'
    }
    gene_lengths: Channel<Path> {
        path 'rnaseq/quant'
    }
    gene_tpm: Channel<Path> {
        path 'rnaseq/quant'
    }
    report: Channel<Path> {
        path 'differentialabundance/report'
    }
}
