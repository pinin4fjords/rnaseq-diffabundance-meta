nextflow.enable.types = true

//
// Gene/transcript quantification with RSEM
//

include { RSEM_CALCULATEEXPRESSION           } from '../../../modules/nf-core/rsem/calculateexpression'
include { CUSTOM_RSEMMERGECOUNTS             } from '../../../modules/nf-core/custom/rsemmergecounts'
include { SENTIEON_RSEMCALCULATEEXPRESSION   } from '../../../modules/nf-core/sentieon/rsemcalculateexpression'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'

workflow QUANTIFY_RSEM {
    take:
    samplesheet: Value<Path>
    ch_samples: Channel<Tuple> // (meta, reads): FASTQ or BAM files
    index: Value<Path> // RSEM index
    gtf: Value<Path>
    gtf_id_attribute: String // GTF gene ID attribute
    gtf_extra_attribute: String // GTF alternative gene attribute (e.g. gene_name)
    use_sentieon_star: Boolean // use Sentieon-accelerated STAR (FASTQ mode only)
    skip_merge: Boolean // skip cross-sample merging, run tximport per-sample

    main:

    //
    // Quantify reads with RSEM
    //
    // (meta, counts_gene, counts_transcript, stat, log, bam_star, bam_genome, bam_transcript)
    def ch_rsem: Channel<Tuple>
    if (use_sentieon_star) {
        // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
        ch_rsem = SENTIEON_RSEMCALCULATEEXPRESSION(ch_samples.map { t -> tuple(t) }.combine(index))
    } else {
        ch_rsem = RSEM_CALCULATEEXPRESSION(ch_samples.map { t -> tuple(t) }.combine(index))
    }

    //
    // Merge counts across samples (only when not skipping merge)
    //
    // Sorted by name for a stable cache key; the script globs the
    // staged directory, so order doesn't affect output. The emptiness
    // check skips CUSTOM_RSEMMERGECOUNTS when there are no samples, since
    // collect emits [] rather than nothing on an empty channel.
    //
    // (meta, rsem_merge(counts_gene, tpm_gene, counts_transcript, tpm_transcript, genes_long, isoforms_long))
    def ch_rsem_merge: Channel<Tuple>
    if (skip_merge) {
        ch_rsem_merge = channel.empty()
    } else {
        ch_rsem_merge = CUSTOM_RSEMMERGECOUNTS(
            ch_rsem
                .collect()
                .flatMap { rs ->
                    rs.isEmpty()
                        ? []
                        : [ tuple(
                            [id: 'all_samples'],
                            rs.collect { _meta, counts_gene, _counts_transcript, _stat, _log, _bam_star, _bam_genome, _bam_transcript -> counts_gene }.toSorted { f -> f.name },
                            rs.collect { _meta, _counts_gene, counts_transcript, _stat, _log, _bam_star, _bam_genome, _bam_transcript -> counts_transcript }.toSorted { f -> f.name }
                        ) ]
                }
        )
    }

    //
    // Post-process quantifications with tximport and SummarizedExperiment
    //
    ch_quant_merged = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT(
        samplesheet,
        ch_rsem.map { meta, _counts_gene, counts_transcript, _stat, _log, _bam_star, _bam_genome, _bam_transcript -> tuple(meta, [ counts_transcript ]) },
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        'rsem',
        skip_merge,
        record(se_prefix: null)
    )

    emit:
    // per sample: (meta, counts_gene, counts_transcript, stat, log, bam_star, bam_genome, bam_transcript)
    samples:    Channel<Tuple> = ch_rsem
    // one row per sample under skip_merge, a single 'all_samples' row otherwise:
    // (meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, tpm_transcript, counts_transcript, lengths_transcript, tx2gene_augmented, merged_gene_rds, merged_transcript_rds, tx2gene)
    merged:     Channel<Tuple> = ch_quant_merged
    // a single 'all_samples' row; empty under skip_merge:
    // (meta, rsem_merge(counts_gene, tpm_gene, counts_transcript, tpm_transcript, genes_long, isoforms_long))
    rsem_merge: Channel<Tuple> = ch_rsem_merge
}
