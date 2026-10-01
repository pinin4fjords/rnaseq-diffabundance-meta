nextflow.enable.types = true

include { params as RnaseqParams ; workflow as NFCORE_RNASEQ } from './pipelines/nf-core/rnaseq'
include { SampleRow                                          } from './pipelines/nf-core/rnaseq/modules/nf-core/types'
include { DIFFERENTIALABUNDANCE                              } from './pipelines/nf-core/differentialabundance/workflows/differentialabundance'
include { diffabParamset                                     } from './lib/diffab'

params {
    samples:         Channel<SampleRow>   // rnaseq samplesheet, one row per sequencing run
    sample_metadata: Path                 // differentialabundance observations: a sample column and the experimental variables
    contrasts:       Path                 // differentialabundance contrasts
    quantification:  String = 'pseudo'    // which merged rnaseq matrices to analyse: 'pseudo' (salmon pseudo-alignment) or 'aligned'
    rnaseq:          RnaseqParams         // references and options of nf-core/rnaseq
}

workflow {
    main:
    // Quantify the samples. rnaseq starts differential abundance as soon as the merged matrices exist.
    rnaseq = NFCORE_RNASEQ( params.rnaseq + record(input: params.samples) )

    ch_quant = params.quantification == 'aligned' ? rnaseq.quant_merged : rnaseq.quant_merged_pseudo

    ch_gtf = rnaseq.genome_references
        .filter { r -> r.kind == 'gtf' }
        .map { r -> r.file }

    // One differentialabundance paramset built from rnaseq's merged gene-level outputs
    ch_paramsets = ch_quant
        .combine(ch_gtf)
        .map { quant, gtf ->
            diffabParamset(quant.counts_gene, quant.lengths_gene, gtf, params.sample_metadata, params.contrasts)
        }

    abundance = DIFFERENTIALABUNDANCE( ch_paramsets )

    publish:
    multiqc      = rnaseq.multiqc.map { r -> r.report }
    gene_counts  = ch_quant.map { r -> r.counts_gene }
    gene_lengths = ch_quant.map { r -> r.lengths_gene }
    gene_tpm     = ch_quant.map { r -> r.tpm_gene }
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
