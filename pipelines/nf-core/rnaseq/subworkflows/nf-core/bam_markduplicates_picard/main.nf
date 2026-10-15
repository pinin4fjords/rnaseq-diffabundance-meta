nextflow.enable.types = true

//
// Picard MarkDuplicates, index BAM file and run samtools stats, flagstat and idxstats
//

include { PICARD_MARKDUPLICATES } from '../../../modules/nf-core/picard/markduplicates/main'
include { SAMTOOLS_INDEX        } from '../../../modules/nf-core/samtools/index/main'
include { BAM_STATS_SAMTOOLS    } from '../bam_stats_samtools/main'

workflow BAM_MARKDUPLICATES_PICARD {
    take:
    ch_bam: Channel<Tuple> // (meta, bam)
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>
    run_stats: Boolean
    tool_args: Record // samtools_index: samtools index options

    main:
    // (meta, bam, bai, cram, metrics)
    // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
    ch_markdup = PICARD_MARKDUPLICATES(
        ch_bam
            .map { t -> tuple(t) }
            .combine(ch_fasta)
            .combine(ch_fai)
    )

    // Picard writes exactly one of bam/cram per sample, matching the input format.
    // (meta, bam_or_cram)
    ch_marked = ch_markdup.map { meta, bam, _bai, cram, _metrics -> tuple(meta, bam ?: cram) }

    // (meta, bai)
    ch_index = SAMTOOLS_INDEX(ch_marked.map { meta, bam -> tuple(meta, bam, tool_args.samtools_index) })

    // (meta, samtools(stats, flagstat, idxstats)), samtools is null when stats are skipped
    def ch_stats = ch_marked.map { meta, _bam -> tuple(meta, null) }
    if (run_stats) {
        ch_stats = BAM_STATS_SAMTOOLS(ch_marked.join(ch_index, by: 0), ch_fasta, ch_fai)
    }

    // The samtools index replaces the Picard bai. A remainder join pads a missing side with a
    // single null, so the multi-element Picard side is nested to keep the tuple width fixed.
    // (meta, markdup(bam, cram, metrics), bai)
    ch_markdup_indexed = ch_markdup
        .map { meta, bam, _bai, cram, metrics -> tuple(meta, tuple(bam, cram, metrics)) }
        .join(ch_index, by: 0, remainder: true)
    ch_markdup_indexed.subscribe { meta, markdup, bai ->
        if( markdup == null || bai == null ) {
            error "Sample '${meta.id}' is missing its samtools index result"
        }
    }
    // (meta, bam, bai, cram, metrics, samtools(stats, flagstat, idxstats))
    ch_results = ch_markdup_indexed
        .filter { _meta, markdup, bai -> markdup != null && bai != null }
        .map { meta, markdup, bai -> tuple(meta, markdup[0], bai, markdup[1], markdup[2]) }
        .join(ch_stats, by: 0)

    emit:
    ch_results
}
