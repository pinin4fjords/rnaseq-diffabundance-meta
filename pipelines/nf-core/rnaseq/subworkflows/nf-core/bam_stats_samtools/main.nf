nextflow.enable.types = true

//
// Run SAMtools stats, flagstat and idxstats
//

include { SAMTOOLS_STATS    } from '../../../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_IDXSTATS } from '../../../modules/nf-core/samtools/idxstats/main'
include { SAMTOOLS_FLAGSTAT } from '../../../modules/nf-core/samtools/flagstat/main'

workflow BAM_STATS_SAMTOOLS {
    take:
    ch_bam_bai: Channel<Record>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>

    main:
    ch_bam_bai_only = ch_bam_bai.map { r -> record(meta: r.meta, bam: r.bam, bai: r.bai) }
    ch_stats = SAMTOOLS_STATS(ch_bam_bai_only.combine(fasta: ch_fasta, fai: ch_fai))
    ch_flagstat = SAMTOOLS_FLAGSTAT(ch_bam_bai_only)
    ch_idxstats = SAMTOOLS_IDXSTATS(ch_bam_bai_only)

    ch_stats_flagstat = ch_stats.join(ch_flagstat, by: 'meta', remainder: true)
    ch_stats_flagstat.subscribe { r ->
        if( r.stats == null || r.flagstat == null ) {
            error "Sample '${r.meta.id}' is missing its samtools flagstat result"
        }
    }
    ch_stats_idxstats = ch_stats_flagstat
        .filter { r -> r.stats != null && r.flagstat != null }
        .join(ch_idxstats, by: 'meta', remainder: true)
    ch_stats_idxstats.subscribe { r ->
        if( r.flagstat == null || r.idxstats == null ) {
            error "Sample '${r.meta.id}' is missing its samtools idxstats result"
        }
    }
    ch_results = ch_stats_idxstats
        .filter { r -> r.flagstat != null && r.idxstats != null }
        .map { r -> record(meta: r.meta, samtools: record(stats: r.stats, flagstat: r.flagstat, idxstats: r.idxstats)) }

    emit:
    ch_results
}
