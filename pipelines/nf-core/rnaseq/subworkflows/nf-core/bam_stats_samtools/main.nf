nextflow.enable.types = true

//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai: Channel<Tuple> // (meta, bam, bai)
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>

    main:
    // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
    ch_stats = SAMTOOLS_STATS(ch_bam_bai.map { t -> tuple(t) }.combine(ch_fasta).combine(ch_fai))
    ch_flagstat = SAMTOOLS_FLAGSTAT(ch_bam_bai)
    ch_idxstats = SAMTOOLS_IDXSTATS(ch_bam_bai)

    ch_stats_flagstat = ch_stats.join(ch_flagstat, by: 0, remainder: true)
    ch_stats_flagstat.subscribe { meta, stats, flagstat ->
        if( stats == null || flagstat == null ) {
            error "Sample '${meta.id}' is missing its samtools flagstat result"
        }
    }
    // A remainder join pads a missing side with a single null, so the two-element side is nested
    // to keep the tuple width fixed.
    ch_stats_idxstats = ch_stats_flagstat
        .filter { _meta, stats, flagstat -> stats != null && flagstat != null }
        .map { meta, stats, flagstat -> tuple(meta, tuple(stats, flagstat)) }
        .join(ch_idxstats, by: 0, remainder: true)
    ch_stats_idxstats.subscribe { meta, stats_flagstat, idxstats ->
        if( stats_flagstat == null || idxstats == null ) {
            error "Sample '${meta.id}' is missing its samtools idxstats result"
        }
    }
    // (meta, samtools(stats, flagstat, idxstats))
    ch_results = ch_stats_idxstats
        .filter { _meta, stats_flagstat, idxstats -> stats_flagstat != null && idxstats != null }
        .map { meta, stats_flagstat, idxstats -> tuple(meta, tuple(stats_flagstat[0], stats_flagstat[1], idxstats)) }

    emit:
    ch_results
}
