nextflow.enable.types = true

include { params as RnaseqParams ; workflow as NFCORE_RNASEQ } from './pipelines/nf-core/rnaseq'
include { params as DiffabParams ; workflow as NFCORE_DIFFERENTIALABUNDANCE } from './pipelines/nf-core/differentialabundance'

params {
    rnaseq:          RnaseqParams         // references and options of nf-core/rnaseq
    diffab:          DiffabParams         // options of nf-core/differentialabundance
}

workflow {
    main:
    // Quantify the samples. rnaseq starts differential abundance as soon as the merged matrices exist.
    rnaseq = NFCORE_RNASEQ(params.rnaseq)

    // differentialabundance takes rnaseq's merged gene-level outputs as its matrix and feature lengths
    def ch_quant = rnaseq.gene_quant.collect().map { quants -> quants.toList().first() }

    abundance = NFCORE_DIFFERENTIALABUNDANCE(
        params.diffab + record(
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
