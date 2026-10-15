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
    ch_bam_bai: Channel<Tuple> // (meta, bam, bai)
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

    // (meta, bam)
    ch_bam = ch_bam_bai.map { meta, bam, _bai -> tuple(meta, bam) }

    // Each tool group contributes one (meta, result) channel. A skipped group contributes a null
    // result for every sample, so the final inner join keeps each sample and leaves that slot null.
    // (meta, null)
    def ch_skipped: Channel<Tuple> = ch_bam_bai.map { meta, _bam, _bai -> tuple(meta, null) }

    //
    // MODULE: Preseq library complexity
    //
    // (meta, preseq(lc_extrap, log))
    def ch_preseq: Channel<Tuple> = ch_skipped
    if ('preseq' in tools) {
        // Remainder join: preseq can fail on low-duplication BAMs, and callers
        // may set errorStrategy 'ignore' rather than lose the sample.
        ch_preseq = ch_skipped
            .join(
                PRESEQ_LCEXTRAP(ch_bam).map { meta, lc_extrap, log -> tuple(meta, tuple(lc_extrap, log)) },
                by: 0,
                remainder: true
            )
            .map { meta, _skipped, preseq -> tuple(meta, preseq) }
    }

    //
    // MODULE: Feature biotype QC using featureCounts
    //
    // (meta, featurecounts(counts, summary))
    def ch_featurecounts: Channel<Tuple> = ch_skipped
    // (meta, biotype(tsv, rrna))
    def ch_biotype: Channel<Tuple> = ch_skipped
    if ('biotype_qc' in tools && biotype) {
        // (meta, counts, summary)
        ch_featurecounts_out = SUBREAD_FEATURECOUNTS(
            ch_bam
                .map { meta, bam -> tuple(meta, bam, ['-B -C', "-g ${biotype}", "-t ${featurecounts_feature_type}"].join(' ')) }
                .combine(ch_gtf)
        )
        ch_featurecounts = ch_featurecounts_out.map { meta, counts, summary -> tuple(meta, tuple(counts, summary)) }
        ch_biotype = CUSTOM_MULTIQCCUSTOMBIOTYPE(ch_featurecounts_out.map { meta, counts, _summary -> tuple(meta, counts, biotypes_header) })
            .map { meta, tsv, rrna -> tuple(meta, tuple(tsv, rrna)) }
    }

    //
    // MODULE: Qualimap
    //
    // (meta, qualimap)
    def ch_qualimap: Channel<Tuple> = ch_skipped
    if ('qualimap' in tools) {
        // (meta, raw_bams, fasta, fai, index_format)
        ch_sort_in = ch_bam
            .map { meta, bam -> tuple(meta, [ bam ]) }
            .combine(ch_fasta)
            .combine(ch_fai)
            .combine(channel.value(''))

        // Name-sorted BAM via samtools sort; requires ext.args = '-n' to be set by the caller for SAMTOOLS_SORT_QUALIMAP
        ch_name_sorted = SAMTOOLS_SORT_QUALIMAP(ch_sort_in)
        ch_qualimap = QUALIMAP_RNASEQ(
            ch_name_sorted
                .map { meta, bam, _cram, _sam, _bai, _csi, _crai -> tuple(meta, bam) }
                .combine(ch_gtf)
        )
    }

    //
    // MODULE: dupRadar
    //
    // (meta, dupradar(scatter2d, boxplot, hist, dupmatrix, intercept_slope, multiqc, session_info))
    def ch_dupradar: Channel<Tuple> = ch_skipped
    if ('dupradar' in tools) {
        // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
        ch_dupradar = DUPRADAR(ch_bam.map { t -> tuple(t) }.combine(ch_gtf))
            .map { meta, scatter2d, boxplot, hist, dupmatrix, intercept_slope, multiqc, session_info ->
                tuple(meta, tuple(scatter2d, boxplot, hist, dupmatrix, intercept_slope, multiqc, session_info))
            }
    }

    //
    // SUBWORKFLOW: RSeQC
    //
    // (meta, rseqc(bamstat, inferexperiment, innerdistance(...), junctionannotation(...), junctionsaturation(...),
    //  readdistribution, readduplication(...), tin(...))), nested layouts as emitted by BAM_RSEQC
    def ch_rseqc: Channel<Tuple> = ch_skipped
    if (rseqc_modules.size() > 0) {
        ch_rseqc = BAM_RSEQC(ch_bam_bai, ch_gene_bed, rseqc_modules)
            .map { meta, bamstat, inferexperiment, innerdistance, junctionannotation, junctionsaturation, readdistribution, readduplication, tin ->
                tuple(meta, tuple(bamstat, inferexperiment, innerdistance, junctionannotation, junctionsaturation, readdistribution, readduplication, tin))
            }
    }

    // (meta, preseq, featurecounts, biotype, qualimap, dupradar, rseqc)
    ch_qc = ch_preseq
        .join(ch_featurecounts, by: 0)
        .join(ch_biotype, by: 0)
        .join(ch_qualimap, by: 0)
        .join(ch_dupradar, by: 0)
        .join(ch_rseqc, by: 0)

    // Files MultiQC reads for each sample, in the order the tools report them.
    // Nested results are indexed by the positions documented on their channels above and in BAM_RSEQC.
    // (meta, preseq(lc_extrap, log), featurecounts(counts, summary), biotype(tsv, rrna), qualimap,
    //  dupradar(scatter2d, boxplot, hist, dupmatrix, intercept_slope, multiqc, session_info),
    //  rseqc(bamstat, inferexperiment, innerdistance(distance, freq, mean, pdf, rscript),
    //        junctionannotation(bed, interact_bed, xls, pdf, events_pdf, rscript, log), junctionsaturation(pdf, rscript),
    //        readdistribution, readduplication(seq_xls, pos_xls, pdf, rscript), tin(txt, xls)),
    //  mqc_files)
    ch_results = ch_qc.map { meta, preseq, featurecounts, biotype_qc, qualimap, dupradar, rseqc ->
        def innerdistance      = rseqc ? rseqc[2] : null
        def junctionannotation = rseqc ? rseqc[3] : null
        def junctionsaturation = rseqc ? rseqc[4] : null
        def readduplication    = rseqc ? rseqc[6] : null
        def tin                = rseqc ? rseqc[7] : null
        def mqc_files = [
            preseq             ? preseq[0]             : null,
            biotype_qc         ? biotype_qc[0]         : null,
            qualimap,
            rseqc              ? rseqc[0]              : null,
            rseqc              ? rseqc[1]              : null,
            innerdistance      ? innerdistance[1]      : null,
            junctionannotation ? junctionannotation[6] : null,
            junctionsaturation ? junctionsaturation[1] : null,
            rseqc              ? rseqc[5]              : null,
            readduplication    ? readduplication[1]    : null,
            tin                ? tin[0]                : null
        ].findAll { f -> f != null }.toList() + ((dupradar ? dupradar[5] : []) as List<Path>)
        tuple(meta, preseq, featurecounts, biotype_qc, qualimap, dupradar, rseqc, mqc_files)
    }

    emit:
    ch_results
}
