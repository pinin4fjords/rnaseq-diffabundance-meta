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

    // rnaseq's merged gene-level matrices: (tpm_gene, counts_gene, lengths_gene), from its gene_quant
    // (meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, tpm_transcript,
    //  counts_transcript, lengths_transcript, tx2gene_augmented, merged_gene_rds, merged_transcript_rds, tx2gene)
    def ch_gene_matrices = rnaseq.gene_quant.map { _meta, tpm_gene, counts_gene, lengths_gene, _counts_gene_length_scaled, _counts_gene_scaled, _tpm_transcript,
                                                   _counts_transcript, _lengths_transcript, _tx2gene_augmented, _merged_gene_rds, _merged_transcript_rds, _tx2gene ->
        tuple(tpm_gene, counts_gene, lengths_gene)
    }

    // differentialabundance takes rnaseq's merged gene-level outputs as its matrix and feature lengths
    def ch_quant = ch_gene_matrices.collect().map { quants -> quants.toList().first() }

    abundance = NFCORE_DIFFERENTIALABUNDANCE(
        params.diffab + record(
            matrix:                ch_quant.map { _tpm_gene, counts_gene, _lengths_gene -> counts_gene },
            feature_length_matrix: ch_quant.map { _tpm_gene, _counts_gene, lengths_gene -> lengths_gene },
            gtf:                   rnaseq.gtf
        )
    )

    publish:
    multiqc      = rnaseq.multiqc.map { _meta, report, _data, _plots -> report }
    gene_counts  = ch_gene_matrices.map { _tpm_gene, counts_gene, _lengths_gene -> counts_gene }
    gene_lengths = ch_gene_matrices.map { _tpm_gene, _counts_gene, lengths_gene -> lengths_gene }
    gene_tpm     = ch_gene_matrices.map { tpm_gene, _counts_gene, _lengths_gene -> tpm_gene }
    report       = abundance.report.filter { name, _meta, _files -> name == 'report_html' }.flatMap { _name, _meta, files -> files }
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
