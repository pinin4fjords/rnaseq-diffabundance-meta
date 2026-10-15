nextflow.enable.types = true

//
// umicollapse, index BAM file and run samtools stats, flagstat and idxstats
//

include { UMICOLLAPSE        } from '../../../modules/nf-core/umicollapse/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

//
// UMICollapse options
//
def umicollapseArgs(grouping_method: String?, umi_separator: String?, single_end: Boolean) -> String {
    def algos = ['directional': 'dir', 'adjacency': 'adj', 'cluster': 'cc']
    def paired_args = single_end ? '' : '--paired --remove-unpaired --remove-chimeric'
    def algo_arg = grouping_method ? "--algo '${algos[grouping_method] ?: ''}'" : ''
    def separator_arg = umi_separator ? "--umi-sep '${umi_separator}'" : ''
    return ['--two-pass', paired_args, algo_arg, separator_arg].join(' ').trim()
}

workflow BAM_DEDUP_STATS_SAMTOOLS_UMICOLLAPSE {
    take:
    ch_bam: Channel<Tuple> // (meta, bam)
    tool_args: Record // samtools_index: samtools index options
    grouping_method: String? // UMI grouping method
    umi_separator: String? // UMI separator

    main:
    //
    // umicollapse in bam mode (thus hardcode mode input to 'bam')
    //
    // (meta, bam, log)
    ch_dedup = UMICOLLAPSE(ch_bam.map { meta, bam -> tuple(meta, bam, umicollapseArgs(grouping_method, umi_separator, meta.single_end), 'bam') })
        .filter { _meta, bam, _fastq, _log -> bam != null }
        .map { meta, bam, _fastq, log -> tuple(meta, bam, log) }

    //
    // Index BAM file and run samtools stats, flagstat and idxstats
    //
    ch_index = SAMTOOLS_INDEX(ch_dedup.map { meta, bam, _log -> tuple(meta, bam, tool_args.samtools_index) })
    // (meta, bam, bai, log)
    ch_indexed = ch_dedup
        .join(ch_index, by: 0)
        .map { meta, bam, log, bai -> tuple(meta, bam, bai, log) }

    // (meta, bam, bai, log, samtools(stats, flagstat, idxstats))
    ch_results = ch_indexed.join(
        BAM_STATS_SAMTOOLS(ch_indexed.map { meta, bam, bai, _log -> tuple(meta, bam, bai) }, channel.value(null as Path?), channel.value(null as Path?)),
        by: 0
    )

    emit:
    ch_results
}
