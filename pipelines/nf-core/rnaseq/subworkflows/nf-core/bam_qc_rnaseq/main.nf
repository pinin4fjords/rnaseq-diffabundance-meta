nextflow.enable.types = true

//
// Run post-alignment QC tools on RNA-seq BAM files
//

include { DUPRADAR                        } from '../../../modules/nf-core/dupradar/main'
include { PRESEQ_LCEXTRAP                 } from '../../../modules/nf-core/preseq/lcextrap/main'
include { QUALIMAP_RNASEQ                 } from '../../../modules/nf-core/qualimap/rnaseq/main'
include { SUBREAD_FEATURECOUNTS           } from '../../../modules/nf-core/subread/featurecounts/main'
include { CUSTOM_MULTIQCCUSTOMBIOTYPE     } from '../../../modules/nf-core/custom/multiqccustombiotype/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_QUALIMAP } from '../../../modules/nf-core/samtools/sort/main'
include { BAM_RSEQC                       } from '../bam_rseqc/main'

workflow BAM_QC_RNASEQ {

    take:
    ch_bam_bai: Channel<Record>
    ch_gtf: Value<Path>
    ch_gene_bed: Value<Path>
    ch_fasta: Value<Path?>
    ch_fai: Value<Path?>
    biotypes_header: Path
    tools: List<String>  // e.g. ['preseq', 'biotype_qc', 'qualimap', 'dupradar', 'rseqc_bam_stat', 'rseqc_infer_experiment', ...]
    biotype: String      // e.g. "gene_type" or "gene_biotype"
    featurecounts_feature_type: String // feature type counted by featureCounts, e.g. "exon"

    main:
    def rseqc_modules = tools.findAll { tool -> tool.startsWith('rseqc_') }.collect { tool -> tool.replace('rseqc_', '') }.toList()

    // Every field starts null and is overwritten by the join of the tool group that ran,
    // so skipped groups leave a null field and no join is made against an empty channel.
    ch_qc = ch_bam_bai.map { r ->
        record(
            meta:          r.meta,
            preseq:        null,
            featurecounts: null,
            biotype:       null,
            qualimap:      null,
            dupradar:      null,
            rseqc:         null
        )
    }

    //
    // MODULE: Preseq library complexity
    //
    if ('preseq' in tools) {
        // Remainder join: preseq can fail on low-duplication BAMs, and callers
        // may set errorStrategy 'ignore' rather than lose the sample.
        ch_qc = ch_qc.join(
            PRESEQ_LCEXTRAP(ch_bam_bai.map { r -> record(meta: r.meta, bam: r.bam) }).map { r -> record(meta: r.meta, preseq: r) },
            by: 'meta',
            remainder: true
        )
    }

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    if ('biotype_qc' in tools && biotype) {
        ch_featurecounts = SUBREAD_FEATURECOUNTS(
            ch_bam_bai
                .map { r -> record(meta: r.meta, bam: r.bam, args: ['-B -C', "-g ${biotype}", "-t ${featurecounts_feature_type}"].join(' ')) }
                .combine(annotation: ch_gtf)
        )
        ch_biotype = CUSTOM_MULTIQCCUSTOMBIOTYPE(ch_featurecounts.map { r -> record(meta: r.meta, counts: r.counts, header: biotypes_header) })
        ch_qc = ch_qc
            .join(ch_featurecounts.map { r -> record(meta: r.meta, featurecounts: r) }, by: 'meta')
            .join(ch_biotype.map { r -> record(meta: r.meta, biotype: record(tsv: r.tsv, rrna: r.rrna)) }, by: 'meta')
    }

    //
    // MODULE: Qualimap
    //
    if ('qualimap' in tools) {
        ch_sort_in = ch_bam_bai
            .map { r -> record(meta: r.meta, raw_bams: [ r.bam ]) }
            .combine(fasta: ch_fasta, fai: ch_fai, index_format: '')

        // Name-sorted BAM via samtools sort; requires ext.args = '-n' to be set by the caller for SAMTOOLS_SORT_QUALIMAP
        ch_name_sorted = SAMTOOLS_SORT_QUALIMAP(ch_sort_in)
        ch_qc = ch_qc.join(QUALIMAP_RNASEQ(ch_name_sorted.map { r -> record(meta: r.meta, bam: r.bam) }.combine(gtf: ch_gtf)), by: 'meta')
    }

    //
    // MODULE: dupRadar
    //
    if ('dupradar' in tools) {
        ch_qc = ch_qc.join(
            DUPRADAR(ch_bam_bai.map { r -> record(meta: r.meta, bam: r.bam) }.combine(gtf: ch_gtf)).map { r -> record(meta: r.meta, dupradar: r) },
            by: 'meta'
        )
    }

    //
    // SUBWORKFLOW: RSeQC
    //
    if (rseqc_modules.size() > 0) {
        ch_qc = ch_qc.join(
            BAM_RSEQC(ch_bam_bai, ch_gene_bed, rseqc_modules).map { r -> record(meta: r.meta, rseqc: r) },
            by: 'meta'
        )
    }

    // Files MultiQC reads for each sample, in the order the tools report them.
    ch_results = ch_qc.map { r ->
        r + record(
            mqc_files: [
                r.preseq?.lc_extrap,
                r.biotype?.tsv,
                r.qualimap,
                r.rseqc?.bamstat,
                r.rseqc?.inferexperiment,
                r.rseqc?.innerdistance?.freq,
                r.rseqc?.junctionannotation?.log,
                r.rseqc?.junctionsaturation?.rscript,
                r.rseqc?.readdistribution,
                r.rseqc?.readduplication?.pos_xls,
                r.rseqc?.tin?.txt
            ].findAll { f -> f != null }.toList() + ((r.dupradar?.multiqc ?: []) as List<Path>)
        )
    }

    emit:
    ch_results
}
