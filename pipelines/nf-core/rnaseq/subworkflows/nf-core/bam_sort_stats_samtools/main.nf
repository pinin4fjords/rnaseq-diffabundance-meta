nextflow.enable.types = true

//
// Sort, index BAM file and run samtools stats, flagstat and idxstats
//

include { SAMTOOLS_SORT      } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS } from '../bam_stats_samtools/main'

workflow BAM_SORT_STATS_SAMTOOLS {
    take:
    ch_bam: Channel<Tuple> // (meta, raw_bams)
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>
    tool_args: Record // samtools_index: samtools index options

    main:
    // (meta, bam, cram, sam, bai, csi, crai)
    // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
    ch_sorted = SAMTOOLS_SORT(
        ch_bam
            .map { t -> tuple(t) }
            .combine(ch_fasta)
            .combine(ch_fai)
            .combine(channel.value(''))
    )
        .filter { _meta, bam, _cram, _sam, _bai, _csi, _crai -> bam != null }
        .map { meta, bam, _cram, _sam, _bai, _csi, _crai -> tuple(meta, bam) }

    ch_index = SAMTOOLS_INDEX(ch_sorted.map { meta, bam -> tuple(meta, bam, tool_args.samtools_index) })
    ch_sorted_indexed = ch_sorted.join(ch_index, by: 0, remainder: true)
    ch_sorted_indexed.subscribe { meta, bam, bai ->
        if( bam == null || bai == null ) {
            error "Sample '${meta.id}' is missing its samtools index result"
        }
    }
    ch_indexed = ch_sorted_indexed.filter { _meta, bam, bai -> bam != null && bai != null }

    // (meta, bam, bai, samtools(stats, flagstat, idxstats))
    ch_results = ch_indexed
        .join(BAM_STATS_SAMTOOLS(ch_indexed, ch_fasta, ch_fai), by: 0)

    emit:
    ch_results
}
