nextflow.enable.types = true

//
// UMI-tools dedup, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMITOOLS_DEDUP                           } from '../../../modules/nf-core/umitools/dedup/main'
include { SAMTOOLS_INDEX                           } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_PRIMARY   } from '../../../modules/nf-core/samtools/view'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_PRIMARY } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS                       } from '../bam_stats_samtools/main'

//
// UMI-tools dedup options
//
def umitoolsDedupArgs(grouping_method, umi_separator, meta) {
    return [
        meta.single_end ? '' : '--unpaired-reads=discard --chimeric-pairs=discard',
        grouping_method ? "--method='${grouping_method}'" : '',
        umi_separator   ? "--umi-separator='${umi_separator}'" : ''
    ].join(' ').trim()
}

workflow BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS {
    take:
    ch_bam_bai: Channel<Tuple> // (meta, bam, bai)
    val_get_dedup_stats: Boolean
    val_primary_only: Boolean
    tool_args: Record // samtools_index: samtools index options
    grouping_method: String? // UMI grouping method
    umi_separator: String? // UMI separator

    main:

    //
    // Optionally filter to primary alignments before deduplication
    //
    // ch_dedup_input: (meta, bam, bai)
    if (val_primary_only) {
        // (meta, bam)
        ch_primary = SAMTOOLS_VIEW_PRIMARY(
            ch_bam_bai.map { meta, bam, bai -> tuple(meta, bam, bai, null, null, null, null, '') }
        )
            .filter { _meta, bam, _bai, _unselected, _unselected_index -> bam != null }
            .map { meta, bam, _bai, _unselected, _unselected_index -> tuple(meta, bam) }

        ch_dedup_input = ch_primary
            .join(SAMTOOLS_INDEX_PRIMARY(ch_primary.map { meta, bam -> tuple(meta, bam, tool_args.samtools_index) }), by: 0)
    }
    else {
        ch_dedup_input = ch_bam_bai
    }

    //
    // UMI-tools dedup
    //
    // (meta, bam, log, tsv_edit_distance, tsv_per_umi, tsv_umi_per_position)
    ch_dedup = UMITOOLS_DEDUP(
        ch_dedup_input.map { meta, bam, bai -> tuple(meta, bam, bai, umitoolsDedupArgs(grouping_method, umi_separator, meta), val_get_dedup_stats) }
    )

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_index = SAMTOOLS_INDEX(ch_dedup.map { meta, bam, _log, _ed, _pu, _upp -> tuple(meta, bam, tool_args.samtools_index) })
    // (meta, bam, bai, log, tsv_edit_distance, tsv_per_umi, tsv_umi_per_position)
    ch_indexed = ch_dedup
        .join(ch_index, by: 0)
        .map { meta, bam, log, ed, pu, upp, bai -> tuple(meta, bam, bai, log, ed, pu, upp) }

    // (meta, bam, bai, log, tsv_edit_distance, tsv_per_umi, tsv_umi_per_position, samtools(stats, flagstat, idxstats))
    ch_results = ch_indexed.join(
        BAM_STATS_SAMTOOLS(ch_indexed.map { meta, bam, bai, _log, _ed, _pu, _upp -> tuple(meta, bam, bai) }, channel.value(null as Path?), channel.value(null as Path?)),
        by: 0
    )

    emit:
    ch_results
}
