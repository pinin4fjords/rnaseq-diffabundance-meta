nextflow.enable.types = true

//
// BAM deduplication with UMI processing
//

include { BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE as BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME } from '../bam_dedup_stats_samtools_umicollapse'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME       } from '../bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE as BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME        } from '../bam_dedup_stats_samtools_umicollapse'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME              } from '../bam_dedup_stats_samtools_umitools'
include { BAM_SORT_STATS_SAMTOOLS                                                                    } from '../bam_sort_stats_samtools'

include { UMITOOLS_PREPAREFORRSEM                                                                    } from '../../../modules/nf-core/umitools/prepareforrsem'
include { SAMTOOLS_SORT                                                                              } from '../../../modules/nf-core/samtools/sort/main'

workflow BAM_DEDUP_UMI {
    take:
    ch_genome_bam: Channel<Tuple> // (meta, bam, bai)
    fasta: Value<Path?>
    fai: Value<Path?>
    umi_dedup_tool: String // 'umicollapse' or 'umitools'
    umitools_dedup_stats: Boolean // whether to generate UMI-tools dedup stats
    ch_transcriptome_bam: Channel<Tuple> // (meta, transcriptome_bam)
    transcript_fasta: Value<Path?>
    umitools_dedup_primary_only: Boolean // whether to filter to primary alignments before dedup
    tool_args: Record // samtools_index: the samtools index options, forwarded
    umi_grouping_method: String? // UMI grouping method
    umi_separator: String? // UMI separator

    main:
    if (umi_dedup_tool != "umicollapse" && umi_dedup_tool != "umitools") {
        error("Unknown umi_dedup_tool '${umi_dedup_tool}'")
    }

    // Both tools are reduced to the same tuple layout; only umitools produces stats tables (tsv).
    // ch_genome_dedup and ch_transcriptome_dedup:
    // (meta, bam, bai, log, samtools(stats, flagstat, idxstats), tsv(edit_distance, per_umi, umi_per_position))

    // Genome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        ch_genome_dedup = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_GENOME(
            ch_genome_bam.map { meta, bam, _bai -> tuple(meta, bam) },
            tool_args,
            umi_grouping_method,
            umi_separator
        )
            .map { meta, bam, bai, log, samtools -> tuple(meta, bam, bai, log, samtools, null) }
    }
    else {
        ch_genome_dedup = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_GENOME(
            ch_genome_bam,
            umitools_dedup_stats,
            umitools_dedup_primary_only,
            tool_args,
            umi_grouping_method,
            umi_separator,
        )
            .map { meta, bam, bai, log, edit_distance, per_umi, umi_per_position, samtools ->
                tuple(
                    meta, bam, bai, log, samtools,
                    edit_distance != null ? tuple(edit_distance, per_umi, umi_per_position) : null
                )
            }
    }

    // Co-ordinate sort, index and run stats on transcriptome BAM. This takes
    // some preparation- we have to coordinate sort the BAM, run the
    // deduplication, then restore name sorting and run a script from umitools
    // to prepare for rsem or salmon

    // 1. Coordinate sort
    // (meta, bam, bai, samtools(stats, flagstat, idxstats))
    ch_coord_sorted = BAM_SORT_STATS_SAMTOOLS(
        ch_transcriptome_bam.map { meta, transcriptome_bam -> tuple(meta, [transcriptome_bam]) },
        transcript_fasta,
        null,
        tool_args
    )

    // 2. Transcriptome BAM deduplication
    if (umi_dedup_tool == "umicollapse") {
        ch_transcriptome_dedup = BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE_TRANSCRIPTOME(
            ch_coord_sorted.map { meta, bam, _bai, _samtools -> tuple(meta, bam) },
            tool_args,
            umi_grouping_method,
            umi_separator
        )
            .map { meta, bam, bai, log, samtools -> tuple(meta, bam, bai, log, samtools, null) }
    }
    else {
        ch_transcriptome_dedup = BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS_TRANSCRIPTOME(
            ch_coord_sorted.map { meta, bam, bai, _samtools -> tuple(meta, bam, bai) },
            umitools_dedup_stats,
            umitools_dedup_primary_only,
            tool_args,
            umi_grouping_method,
            umi_separator,
        )
            .map { meta, bam, bai, log, edit_distance, per_umi, umi_per_position, samtools ->
                tuple(
                    meta, bam, bai, log, samtools,
                    edit_distance != null ? tuple(edit_distance, per_umi, umi_per_position) : null
                )
            }
    }

    // 3. Restore name sorting
    // (meta, bam)
    ch_name_sorted = SAMTOOLS_SORT(
        ch_transcriptome_dedup
            .map { meta, bam, _bai, _log, _samtools, _tsv -> tuple(meta, [bam]) }
            .combine(fasta)
            .combine(fai)
            .combine(channel.value(''))
    )
        .filter { _meta, bam, _cram, _sam, _bai, _csi, _crai -> bam != null }
        .map { meta, bam, _cram, _sam, _bai, _csi, _crai -> tuple(meta, bam) }

    // 4. Run prepare_for_rsem.py on paired-end BAM files
    // This fixes paired-end reads in name sorted BAM files
    // See: https://github.com/nf-core/rnaseq/issues/828
    // (meta, bam, log)
    ch_prepared = UMITOOLS_PREPAREFORRSEM(ch_name_sorted.filter { meta, _bam -> !meta.single_end })

    // Only paired-end samples pass through UMITOOLS_PREPAREFORRSEM, so the remainder
    // join leaves `prepared` null for single-end samples. That side is nested so the
    // single null padding keeps the tuple width fixed.
    // (meta, transcriptomic_dedup_log, prepare_for_rsem_log, transcriptome_bam,
    //  transcriptome(transcriptome_bam, dedup_bam, sorted_bam, sorted_bam_index, filtered_bam, samtools, tsv,
    //                coord_sorted_bam, coord_sorted_bam_index, coord_sorted_samtools))
    ch_transcriptome_part = ch_transcriptome_bam
        .join(ch_coord_sorted, by: 0)
        .join(ch_transcriptome_dedup, by: 0)
        .join(ch_name_sorted, by: 0)
        .join(ch_prepared.map { meta, bam, log -> tuple(meta, tuple(bam, log)) }, by: 0, remainder: true)
        .map { meta, transcriptome_bam,
               coord_sorted_bam, coord_sorted_bai, coord_sorted_samtools,
               dedup_bam, dedup_bai, dedup_log, dedup_samtools, dedup_tsv,
               name_sorted_bam, prepared ->
            def (prepared_bam, prepared_log) = prepared ?: [null, null]
            tuple(
                meta,
                dedup_log,
                prepared_log,
                prepared_bam ?: name_sorted_bam,
                tuple(
                    transcriptome_bam,
                    dedup_bam,
                    name_sorted_bam,
                    dedup_bai,
                    prepared_bam,
                    dedup_samtools,
                    dedup_tsv,
                    coord_sorted_bam,
                    coord_sorted_bai,
                    coord_sorted_samtools
                )
            )
        }

    // The transcriptome side is absent when there is no transcriptome BAM (e.g. HISAT2). Both
    // sides of the remainder join are nested so the single null padding keeps the tuple width fixed.
    // (meta, bam, bai, samtools(stats, flagstat, idxstats), genomic_dedup_log,
    //  tsv(edit_distance, per_umi, umi_per_position),
    //  transcriptomic_dedup_log, prepare_for_rsem_log, transcriptome_bam,
    //  transcriptome(transcriptome_bam, dedup_bam, sorted_bam, sorted_bam_index, filtered_bam, samtools, tsv,
    //                coord_sorted_bam, coord_sorted_bam_index, coord_sorted_samtools))
    ch_results = ch_genome_dedup
        .map { meta, bam, bai, log, samtools, tsv -> tuple(meta, tuple(bam, bai, samtools, log, tsv)) }
        .join(
            ch_transcriptome_part.map { meta, dedup_log, prepare_log, transcriptome_bam, transcriptome ->
                tuple(meta, tuple(dedup_log, prepare_log, transcriptome_bam, transcriptome))
            },
            by: 0,
            remainder: true
        )
        .map { meta, genome, transcriptome_part ->
            def (bam, bai, samtools, genomic_dedup_log, tsv) = genome ?: [null, null, null, null, null]
            def (transcriptomic_dedup_log, prepare_for_rsem_log, transcriptome_bam, transcriptome) = transcriptome_part ?: [null, null, null, null]
            tuple(
                meta, bam, bai, samtools, genomic_dedup_log, tsv,
                transcriptomic_dedup_log, prepare_for_rsem_log, transcriptome_bam, transcriptome
            )
        }

    emit:
    ch_results
}
