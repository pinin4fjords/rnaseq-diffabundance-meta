nextflow.enable.types = true

include { params as RnaseqParams ; workflow as NFCORE_RNASEQ } from './pipelines/nf-core/rnaseq'
include { SampleRow                                          } from './pipelines/nf-core/rnaseq/modules/nf-core/types'
include { DIFFERENTIALABUNDANCE                              } from './pipelines/nf-core/differentialabundance/workflows/differentialabundance'
include { diffabParamset                                     } from './lib/diffab'

params {
    samples:         Channel<SampleRow>   // rnaseq samplesheet, one row per sequencing run
    sample_metadata: Path                 // differentialabundance observations: a sample column and the experimental variables
    contrasts:       Path                 // differentialabundance contrasts
    rnaseq:          RnaseqParams         // references and options of nf-core/rnaseq
}

workflow {
    main:
    // Quantify the samples. rnaseq starts differential abundance as soon as the merged matrices exist.
    rnaseq = NFCORE_RNASEQ( params.rnaseq + record(input: params.samples) )

    // One differentialabundance paramset built from rnaseq's merged gene-level outputs
    ch_paramsets = rnaseq.gene_quant
        .combine(gtf: rnaseq.gtf)
        .map { quant ->
            diffabParamset(quant.counts_gene, quant.lengths_gene, quant.gtf, params.sample_metadata, params.contrasts)
        }

    abundance = DIFFERENTIALABUNDANCE( ch_paramsets )

    publish:
    multiqc      = rnaseq.multiqc.map { r -> r.report }
    gene_counts  = rnaseq.gene_quant.map { r -> r.counts_gene }
    gene_lengths = rnaseq.gene_quant.map { r -> r.lengths_gene }
    gene_tpm     = rnaseq.gene_quant.map { r -> r.tpm_gene }
    report       = abundance.report_html.map { r -> r[1] }
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
