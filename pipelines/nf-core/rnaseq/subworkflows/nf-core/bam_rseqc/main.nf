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
    ch_bam_bai: Channel<Tuple> // (meta, bam, bai)
    bed: Value<Path>
    rseqc_modules: List<String>

    main:
    // (meta, bam, bai, bed)
    ch_bam_bai_bed = ch_bam_bai.combine(bed)

    // Each tool contributes one (meta, result) channel. A skipped tool contributes a null result
    // for every sample, so the final inner join keeps each sample and leaves that slot null.
    // (meta, null)
    def ch_skipped: Channel<Tuple> = ch_bam_bai.map { meta, _bam, _bai -> tuple(meta, null) }

    //
    // Run RSeQC bam_stat.py
    //
    // (meta, bamstat)
    def ch_bamstat: Channel<Tuple> = ch_skipped
    if ('bam_stat' in rseqc_modules) {
        ch_bamstat = RSEQC_BAMSTAT(ch_bam_bai)
    }

    //
    // Run RSeQC inner_distance.py
    //
    // (meta, innerdistance(distance, freq, mean, pdf, rscript))
    def ch_innerdistance: Channel<Tuple> = ch_skipped
    if ('inner_distance' in rseqc_modules) {
        ch_innerdistance = RSEQC_INNERDISTANCE(ch_bam_bai_bed)
            .map { meta, distance, freq, mean, pdf, rscript -> tuple(meta, tuple(distance, freq, mean, pdf, rscript)) }
    }

    //
    // Run RSeQC infer_experiment.py
    //
    // (meta, inferexperiment)
    def ch_inferexperiment: Channel<Tuple> = ch_skipped
    if ('infer_experiment' in rseqc_modules) {
        ch_inferexperiment = RSEQC_INFEREXPERIMENT(ch_bam_bai_bed)
    }

    //
    // Run RSeQC junction_annotation.py
    //
    // (meta, junctionannotation(bed, interact_bed, xls, pdf, events_pdf, rscript, log))
    def ch_junctionannotation: Channel<Tuple> = ch_skipped
    if ('junction_annotation' in rseqc_modules) {
        ch_junctionannotation = RSEQC_JUNCTIONANNOTATION(ch_bam_bai_bed)
            .map { meta, bed_file, interact_bed, xls, pdf, events_pdf, rscript, log ->
                tuple(meta, tuple(bed_file, interact_bed, xls, pdf, events_pdf, rscript, log))
            }
    }

    //
    // Run RSeQC junction_saturation.py
    //
    // (meta, junctionsaturation(pdf, rscript))
    def ch_junctionsaturation: Channel<Tuple> = ch_skipped
    if ('junction_saturation' in rseqc_modules) {
        ch_junctionsaturation = RSEQC_JUNCTIONSATURATION(ch_bam_bai_bed)
            .map { meta, pdf, rscript -> tuple(meta, tuple(pdf, rscript)) }
    }

    //
    // Run RSeQC read_distribution.py
    //
    // (meta, readdistribution)
    def ch_readdistribution: Channel<Tuple> = ch_skipped
    if ('read_distribution' in rseqc_modules) {
        ch_readdistribution = RSEQC_READDISTRIBUTION(ch_bam_bai_bed)
    }

    //
    // Run RSeQC read_duplication.py
    //
    // (meta, readduplication(seq_xls, pos_xls, pdf, rscript))
    def ch_readduplication: Channel<Tuple> = ch_skipped
    if ('read_duplication' in rseqc_modules) {
        ch_readduplication = RSEQC_READDUPLICATION(ch_bam_bai)
            .map { meta, seq_xls, pos_xls, pdf, rscript -> tuple(meta, tuple(seq_xls, pos_xls, pdf, rscript)) }
    }

    //
    // Run RSeQC tin.py
    //
    // (meta, tin(txt, xls))
    def ch_tin: Channel<Tuple> = ch_skipped
    if ('tin' in rseqc_modules) {
        ch_tin = RSEQC_TIN(ch_bam_bai_bed)
            .map { meta, txt, xls -> tuple(meta, tuple(txt, xls)) }
    }

    // (meta, bamstat, inferexperiment, innerdistance(distance, freq, mean, pdf, rscript),
    //  junctionannotation(bed, interact_bed, xls, pdf, events_pdf, rscript, log), junctionsaturation(pdf, rscript),
    //  readdistribution, readduplication(seq_xls, pos_xls, pdf, rscript), tin(txt, xls))
    ch_results = ch_bamstat
        .join(ch_inferexperiment, by: 0)
        .join(ch_innerdistance, by: 0)
        .join(ch_junctionannotation, by: 0)
        .join(ch_junctionsaturation, by: 0)
        .join(ch_readdistribution, by: 0)
        .join(ch_readduplication, by: 0)
        .join(ch_tin, by: 0)

    emit:
    ch_results
}
