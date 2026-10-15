nextflow.enable.types = true

include { HISAT2_ALIGN            } from '../../../modules/nf-core/hisat2/align/main'
include { BAM_SORT_STATS_SAMTOOLS } from '../bam_sort_stats_samtools/main'

workflow FASTQ_ALIGN_HISAT2 {
    take:
    ch_samples: Channel<Tuple> // (meta, reads, args)
    index: Value<Path>
    splicesites: Value<Path?>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>
    save_unaligned: Boolean
    tool_args: Record // samtools_index: samtools index options, forwarded

    main:
    //
    // Map reads with HISAT2
    //
    // (meta, raw_bams, unmapped, hisat2_summary)
    ch_hisat2 = HISAT2_ALIGN(
        ch_samples
            .combine(index)
            .combine(splicesites)
            .combine(channel.value(save_unaligned))
    )

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    // (meta, bam, bai, samtools(stats, flagstat, idxstats))
    ch_sorted = BAM_SORT_STATS_SAMTOOLS(
        ch_hisat2.map { meta, raw_bams, _unmapped, _hisat2_summary -> tuple(meta, raw_bams) },
        ch_fasta,
        ch_fai,
        tool_args
    )

    // (meta, raw_bams, unmapped, hisat2_summary, aligner, bam, bai, samtools(stats, flagstat, idxstats))
    ch_results = ch_hisat2
        .map { meta, raw_bams, unmapped, hisat2_summary -> tuple(meta, raw_bams, unmapped, hisat2_summary, 'hisat2') }
        .join(ch_sorted, by: 0)

    emit:
    ch_results
}
