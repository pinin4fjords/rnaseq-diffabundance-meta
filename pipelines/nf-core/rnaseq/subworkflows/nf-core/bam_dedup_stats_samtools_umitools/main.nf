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
    ch_bam_bai: Channel<Record>
    val_get_dedup_stats: Boolean
    val_primary_only: Boolean
    tool_args: Record // samtools_index: samtools index options
    grouping_method: String? // UMI grouping method
    umi_separator: String? // UMI separator

    main:

    //
    // Optionally filter to primary alignments before deduplication
    //
    if (val_primary_only) {
        ch_primary = SAMTOOLS_VIEW_PRIMARY(
            ch_bam_bai.map { r -> record(meta: r.meta, bam: r.bam, bai: r.bai, fasta: null, fai: null, qname: null, bed: null, index_format: '') }
        ).filter { r -> r.bam != null }

        ch_dedup_input = ch_primary
            .map { r -> record(meta: r.meta, bam: r.bam) }
            .join(SAMTOOLS_INDEX_PRIMARY(ch_primary.map { r -> record(meta: r.meta, bam: r.bam, args: tool_args.samtools_index) }), by: 'meta')
    }
    else {
        ch_dedup_input = ch_bam_bai.map { r -> record(meta: r.meta, bam: r.bam, bai: r.bai) }
    }

    //
    // UMI-tools dedup
    //
    ch_dedup = UMITOOLS_DEDUP(
        ch_dedup_input.map { r -> record(meta: r.meta, bam: r.bam, bai: r.bai, args: umitoolsDedupArgs(grouping_method, umi_separator, r.meta), get_output_stats: val_get_dedup_stats) }
    )

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_index = SAMTOOLS_INDEX(ch_dedup.map { r -> record(meta: r.meta, bam: r.bam, args: tool_args.samtools_index) })
    ch_indexed = ch_dedup.join(ch_index, by: 'meta')

    ch_results = ch_indexed.join(BAM_STATS_SAMTOOLS(ch_indexed, null, null), by: 'meta')

    emit:
    ch_results
}
