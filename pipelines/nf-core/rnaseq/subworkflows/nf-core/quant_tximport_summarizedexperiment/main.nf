nextflow.enable.types = true

//
// Quantification post-processing with tximport and SummarizedExperiment
//

include { CUSTOM_TX2GENE   } from '../../../modules/nf-core/custom/tx2gene'
include { TXIMETA_TXIMPORT } from '../../../modules/nf-core/tximeta/tximport'

include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_GENE_UNIFIED       } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'
include { SUMMARIZEDEXPERIMENT_SUMMARIZEDEXPERIMENT as SE_TRANSCRIPT_UNIFIED } from '../../../modules/nf-core/summarizedexperiment/summarizedexperiment'

workflow QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT {
    take:
    samplesheet: Value<Path>
    ch_quants: Channel<Tuple> // (meta, quants): per-sample quantification results
    gtf: Value<Path>
    gtf_id_attribute: String // GTF gene ID attribute
    gtf_extra_attribute: String // GTF alternative gene attribute (e.g. gene_name)
    quant_type: String // 'salmon', 'kallisto', or 'rsem'
    skip_merge: Boolean // skip cross-sample merging, run tximport per-sample
    tool_args: Record // se_prefix: prefix of the SummarizedExperiment files (default: the sample id)

    main:

    //
    // Create tx2gene mapping from GTF + a single sample's quantification files.
    // tx2gene only uses quant files to discover which GTF attribute corresponds
    // to transcript IDs, so a single sample suffices. This assumes all samples
    // were quantified against the same transcriptome, avoiding the need to run
    // tx2gene independently per sample. If a use case arises requiring mixed
    // transcriptomes, this will need to be revisited.
    //
    // The sample is selected by sorting on the staged results name so the same
    // one is chosen on every run, keeping the CUSTOM_TX2GENE cache key stable
    // across -resume. Picking by arrival order would vary between runs and
    // invalidate the cache for tx2gene and everything downstream of it. The
    // empty-list guard keeps tx2gene from running when no quant results are
    // supplied, since collect still emits an empty list in that case.
    //
    // (meta, quants)
    ch_tx2gene_quants = ch_quants
        .collect()
        .flatMap { rs -> rs.isEmpty() ? [] : [ tuple([id: null], [ rs.collect { _meta, quants -> quants[0] }.toSorted { f -> f.name }.first() ]) ] }

    // (meta, tx2gene)
    ch_tx2gene = CUSTOM_TX2GENE(
        ch_tx2gene_quants
            .combine(gtf)
            .map { meta, quants, gtf_file -> tuple(meta, quants, gtf_file, quant_type, gtf_id_attribute, gtf_extra_attribute) }
    )
    // toSorted treats a two-parameter closure as a comparator, so the key closure indexes the tuple
    def ch_tx2gene_file: Value<Path?> = ch_tx2gene
        .collect()
        .map { rs ->
            if( rs.isEmpty() ) {
                return null
            }
            def (_meta, tx2gene) = rs.toSorted { t -> t[0].id }.first()
            return tx2gene
        }

    //
    // Import and summarize quantifications with tximport
    // In per-sample mode, run once per sample instead of collecting all
    //
    // Sorted by name for a stable cache key; the R script derives sample
    // identity from staged file names, not list position. Filtered to
    // skip TXIMETA_TXIMPORT when there are no samples, since
    // collect emits [] rather than nothing on an empty channel.
    //
    if (skip_merge) {
        ch_tximport_input = ch_quants
    } else {
        ch_tximport_input = ch_quants
            .collect()
            .flatMap { rs -> rs.isEmpty() ? [] : [ tuple([id: 'all_samples'], rs.collect { _meta, quants -> quants[0] }.toSorted { f -> f.name }) ] }
    }

    // (meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, tpm_transcript, counts_transcript, lengths_transcript, tx2gene_augmented)
    ch_tximport = TXIMETA_TXIMPORT(
        ch_tximport_input
            .combine(ch_tx2gene_file)
            .map { meta, quants, tx2gene -> tuple(meta, quants, tx2gene, quant_type) }
    )

    //
    // Build SummarizedExperiment objects (only when merging)
    //
    // ch_se_gene: (meta, merged_gene_rds), ch_se_transcript: (meta, merged_transcript_rds)
    if (skip_merge) {
        ch_se_gene       = ch_tximport.map { meta, _tpm_gene, _counts_gene, _lengths_gene, _counts_gene_length_scaled, _counts_gene_scaled, _tpm_transcript, _counts_transcript, _lengths_transcript, _tx2gene_augmented -> tuple(meta, null) }
        ch_se_transcript = ch_tximport.map { meta, _tpm_gene, _counts_gene, _lengths_gene, _counts_gene_length_scaled, _counts_gene_scaled, _tpm_transcript, _counts_transcript, _lengths_transcript, _tx2gene_augmented -> tuple(meta, null) }
    } else {
        //
        // Build gene-level SummarizedExperiment
        //
        ch_se_gene = SE_GENE_UNIFIED(
            ch_tximport.map { meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, _tpm_transcript, _counts_transcript, _lengths_transcript, _tx2gene_augmented ->
                tuple(meta, [ counts_gene, counts_gene_length_scaled, counts_gene_scaled, lengths_gene, tpm_gene ], tool_args.se_prefix ? "${tool_args.se_prefix}.gene" : "${meta.id}_gene")
            }
                .combine(ch_tx2gene_file)
                .combine(samplesheet)
        ).map { meta, rds, _log -> tuple(meta, rds) }

        //
        // Build transcript-level SummarizedExperiment
        //
        def ch_tx2gene_augmented: Value<Path?> = ch_tximport
            .collect()
            .map { rs ->
                if( rs.isEmpty() ) {
                    return null
                }
                def (_meta, _tpm_gene, _counts_gene, _lengths_gene, _counts_gene_length_scaled, _counts_gene_scaled, _tpm_transcript, _counts_transcript, _lengths_transcript, tx2gene_augmented) = rs.toSorted { t -> t[0].id }.first()
                return tx2gene_augmented
            }
        ch_se_transcript = SE_TRANSCRIPT_UNIFIED(
            ch_tximport.map { meta, _tpm_gene, _counts_gene, _lengths_gene, _counts_gene_length_scaled, _counts_gene_scaled, tpm_transcript, counts_transcript, lengths_transcript, _tx2gene_augmented ->
                tuple(meta, [ counts_transcript, lengths_transcript, tpm_transcript ], tool_args.se_prefix ? "${tool_args.se_prefix}.transcript" : "${meta.id}_transcript")
            }
                .combine(ch_tx2gene_augmented)
                .combine(samplesheet)
        ).map { meta, rds, _log -> tuple(meta, rds) }
    }

    //
    // One tuple per TXIMETA_TXIMPORT row: a single 'all_samples' row when
    // merging, or one per sample under skip_merge. The SE outputs only exist
    // when merging and are null otherwise.
    //
    // (meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, tpm_transcript, counts_transcript, lengths_transcript, tx2gene_augmented, merged_gene_rds, merged_transcript_rds, tx2gene)
    ch_results = ch_tximport
        .join(ch_se_gene, by: 0)
        .join(ch_se_transcript, by: 0)
        .combine(ch_tx2gene_file)

    emit:
    ch_results
}
