nextflow.enable.types = true

//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'

workflow BAM_MARKDUPLICATES_PICARD {
    take:
    ch_bam: Channel<Record>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>
    run_stats: Boolean
    tool_args: Record // samtools_index: samtools index options

    main:
    ch_markdup = PICARD_MARKDUPLICATES(
        ch_bam
            .map { r -> record(meta: r.meta, bam: r.bam) }
            .combine(fasta: ch_fasta, fai: ch_fai)
    )

    // Picard writes exactly one of bam/cram per sample, matching the input format.
    ch_marked = ch_markdup.map { r -> record(meta: r.meta, bam: r.bam ?: r.cram) }

    ch_index = SAMTOOLS_INDEX(ch_marked.map { r -> record(meta: r.meta, bam: r.bam, args: tool_args.samtools_index) })

    def ch_stats = ch_marked.map { r -> record(meta: r.meta, samtools: null as Record?) }
    if (run_stats) {
        ch_stats = BAM_STATS_SAMTOOLS(ch_marked.join(ch_index, by: 'meta'), ch_fasta, ch_fai)
    }

    ch_markdup_indexed = ch_markdup.join(ch_index, by: 'meta', remainder: true)
    ch_markdup_indexed.subscribe { r ->
        if( r.metrics == null || r.bai == null ) {
            error "Sample '${r.meta.id}' is missing its samtools index result"
        }
    }
    ch_results = ch_markdup_indexed
        .filter { r -> r.metrics != null && r.bai != null }
        .join(ch_stats, by: 'meta')

    emit:
    ch_results
}
