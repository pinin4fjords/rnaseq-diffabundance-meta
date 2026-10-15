nextflow.enable.types = true

//
// Run RSeQC modules
//

include { RSEQC_BAMSTAT            } from '../../../modules/nf-core/rseqc/bamstat/main'
include { RSEQC_INNERDISTANCE      } from '../../../modules/nf-core/rseqc/innerdistance/main'
include { RSEQC_INFEREXPERIMENT    } from '../../../modules/nf-core/rseqc/inferexperiment/main'
include { RSEQC_JUNCTIONANNOTATION } from '../../../modules/nf-core/rseqc/junctionannotation/main'
include { RSEQC_JUNCTIONSATURATION } from '../../../modules/nf-core/rseqc/junctionsaturation/main'
include { RSEQC_READDISTRIBUTION   } from '../../../modules/nf-core/rseqc/readdistribution/main'
include { RSEQC_READDUPLICATION    } from '../../../modules/nf-core/rseqc/readduplication/main'
include { RSEQC_TIN                } from '../../../modules/nf-core/rseqc/tin/main'

workflow BAM_RSEQC {
    take:
    ch_bam_bai: Channel<Record>
    bed: Value<Path>
    rseqc_modules: List<String>

    main:
    ch_bam_bai_only = ch_bam_bai.map { r -> record(meta: r.meta, bam: r.bam, bai: r.bai) }
    ch_bam_bai_bed = ch_bam_bai_only.combine(bed: bed)

    // Every field starts null and is overwritten by the join of the tool that ran, so
    // skipped tools leave a null field and no join is made against an empty channel.
    ch_results = ch_bam_bai.map { r ->
        record(
            meta:               r.meta,
            bamstat:            null,
            inferexperiment:    null,
            innerdistance:      null,
            junctionannotation: null,
            junctionsaturation: null,
            readdistribution:   null,
            readduplication:    null,
            tin:                null
        )
    }

    //
    // Run RSeQC bam_stat.py
    //
    if ('bam_stat' in rseqc_modules) {
        ch_results = ch_results.join(RSEQC_BAMSTAT(ch_bam_bai_only), by: 'meta')
    }

    //
    // Run RSeQC inner_distance.py
    //
    if ('inner_distance' in rseqc_modules) {
        ch_results = ch_results.join(
            RSEQC_INNERDISTANCE(ch_bam_bai_bed).map { r -> record(meta: r.meta, innerdistance: r) },
            by: 'meta'
        )
    }

    //
    // Run RSeQC infer_experiment.py
    //
    if ('infer_experiment' in rseqc_modules) {
        ch_results = ch_results.join(RSEQC_INFEREXPERIMENT(ch_bam_bai_bed), by: 'meta')
    }

    //
    // Run RSeQC junction_annotation.py
    //
    if ('junction_annotation' in rseqc_modules) {
        ch_results = ch_results.join(
            RSEQC_JUNCTIONANNOTATION(ch_bam_bai_bed).map { r -> record(meta: r.meta, junctionannotation: r) },
            by: 'meta'
        )
    }

    //
    // Run RSeQC junction_saturation.py
    //
    if ('junction_saturation' in rseqc_modules) {
        ch_results = ch_results.join(
            RSEQC_JUNCTIONSATURATION(ch_bam_bai_bed).map { r -> record(meta: r.meta, junctionsaturation: r) },
            by: 'meta'
        )
    }

    //
    // Run RSeQC read_distribution.py
    //
    if ('read_distribution' in rseqc_modules) {
        ch_results = ch_results.join(RSEQC_READDISTRIBUTION(ch_bam_bai_bed), by: 'meta')
    }

    //
    // Run RSeQC read_duplication.py
    //
    if ('read_duplication' in rseqc_modules) {
        ch_results = ch_results.join(
            RSEQC_READDUPLICATION(ch_bam_bai_only).map { r -> record(meta: r.meta, readduplication: r) },
            by: 'meta'
        )
    }

    //
    // Run RSeQC tin.py
    //
    if ('tin' in rseqc_modules) {
        ch_results = ch_results.join(
            RSEQC_TIN(ch_bam_bai_bed).map { r -> record(meta: r.meta, tin: r) },
            by: 'meta'
        )
    }

    emit:
    ch_results
}
