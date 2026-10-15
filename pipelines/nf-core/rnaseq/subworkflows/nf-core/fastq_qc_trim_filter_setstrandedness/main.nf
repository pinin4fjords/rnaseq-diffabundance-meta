nextflow.enable.types = true

include { BBMAP_BBSPLIT                         } from '../../../modules/nf-core/bbmap/bbsplit'
include { FASTQC as FASTQC_FILTERED             } from '../../../modules/nf-core/fastqc'
include { CAT_FASTQ                             } from '../../../modules/nf-core/cat/fastq/main'
include { FQ_LINT                               } from '../../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_TRIMMING     } from '../../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_BBSPLIT      } from '../../../modules/nf-core/fq/lint/main'
include { FQ_LINT as FQ_LINT_AFTER_RIBO_REMOVAL } from '../../../modules/nf-core/fq/lint/main'
include { FASTQ_REMOVE_RRNA                     } from '../fastq_remove_rrna'
include { FASTQ_SUBSAMPLE_FQ_SALMON             } from '../fastq_subsample_fq_salmon'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE      } from '../fastq_fastqc_umitools_trimgalore'
include { FASTQ_FASTQC_UMITOOLS_FASTP           } from '../fastq_fastqc_umitools_fastp'

//
// Function to determine library type by comparing type counts.
//

//
def calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    def totalFragments = forwardFragments + reverseFragments + unstrandedFragments
    def totalStrandedFragments = forwardFragments + reverseFragments

    def strandedness = 'undetermined'
    if (totalStrandedFragments > 0) {
        def forwardProportion = forwardFragments / (totalStrandedFragments as double)
        def reverseProportion = reverseFragments / (totalStrandedFragments as double)
        def proportionDifference = Math.abs(forwardProportion - reverseProportion)

        if (forwardProportion >= stranded_threshold) {
            strandedness = 'forward'
        }
        else if (reverseProportion >= stranded_threshold) {
            strandedness = 'reverse'
        }
        else if (proportionDifference <= unstranded_threshold) {
            strandedness = 'unstranded'
        }
    }

    return [
        inferred_strandedness: strandedness,
        forwardFragments: (forwardFragments / (totalFragments as double)) * 100,
        reverseFragments: (reverseFragments / (totalFragments as double)) * 100,
        unstrandedFragments: (unstrandedFragments / (totalFragments as double)) * 100,
    ]
}

//
// Function that parses Salmon quant 'lib_format_counts.json' output file to get inferred strandedness
//
def getSalmonInferredStrandedness(json_file, stranded_threshold = 0.8, unstranded_threshold = 0.1) {
    // Parse the JSON content of the file
    def libCounts = new groovy.json.JsonSlurper().parseText(json_file.text) as Map<String, Double>

    // Calculate unstranded fragments (IU and U)
    // NOTE: this is here for completeness, but actually all fragments have a
    // strandedness (even if the overall library does not), so all these values
    // will be '0'. See
    // https://groups.google.com/g/sailfish-users/c/yxzBDv6NB6I
    def unstrandedKeys = ['IU', 'U', 'MU']
    def unstrandedFragments = unstrandedKeys.collect { key -> libCounts[key] ?: 0.0 }.sum()

    // The per-format keys (SF, ISF, MSF, OSF, ...) are not fragment-exclusive:
    // a fragment's candidate placements can increment both the sense and
    // antisense variant. `strand_mapping_bias` is salmon's fragment-exclusive
    // measure of the sense share of the resolved orientation instead.
    def numAssignedFragments = (libCounts['num_assigned_fragments'] ?: 0.0) as double
    def strandMappingBias = (libCounts['strand_mapping_bias'] ?: 0.0) as double
    def forwardFragments = strandMappingBias * numAssignedFragments
    def reverseFragments = (1 - strandMappingBias) * numAssignedFragments

    // Use shared calculation function to determine strandedness
    return calculateStrandedness(forwardFragments, reverseFragments, unstrandedFragments, stranded_threshold, unstranded_threshold)
}

// The emitted record of one sample, from the record joined across the stages
def preprocessedSample(r, meta, skip_linting) {
    def has_fastqc = [
        r.fastqc_raw_html, r.fastqc_raw_zip,
        r.fastqc_trim_html, r.fastqc_trim_zip,
        r.fastqc_filtered_html, r.fastqc_filtered_zip
    ].any { f -> f != null }
    return record(
        meta:              meta,
        reads:             r.reads,
        reads_cat:         r.reads_cat,
        reads_trimmed:     r.reads_trimmed,
        num_trimmed_reads: r.num_trimmed_reads,
        fastqc:            has_fastqc
            ? record(
                raw_html:      r.fastqc_raw_html,
                raw_zip:       r.fastqc_raw_zip,
                trim_html:     r.fastqc_trim_html,
                trim_zip:      r.fastqc_trim_zip,
                filtered_html: r.fastqc_filtered_html,
                filtered_zip:  r.fastqc_filtered_zip
            )
            : null,
        trim:              r.trim,
        umi:               r.umi,
        bbsplit:           r.bbsplit,
        lint:              skip_linting
            ? null
            : record(
                raw:     r.lint_raw,
                trimmed: r.lint_trimmed,
                bbsplit: r.lint_bbsplit,
                ribo:    r.lint_ribo
            ),
        rrna:              r.rrna
    )
}

workflow FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS {
    take:
    // Input channels
    ch_reads: Channel<Record> // meta, reads
    ch_fasta: Value<Path> // genome.fasta
    ch_transcript_fasta: Value<Path> // transcript.fasta
    ch_gtf: Value<Path> // genome.gtf
    ch_salmon_index: Value<Path> // salmon/index/ (optional)
    ch_sortmerna_index: Value<Path> // sortmerna/index/ (optional)
    ch_bowtie2_index: Value<Path> // bowtie2/index/ (optional)
    ch_bbsplit_index: Value<Path> // bbsplit/index/ (optional)
    ch_rrna_fastas: Channel<Path> // one or more fasta files containing rrna sequences to be passed to SortMeRNA/Bowtie2 (optional)

    // Skip options
    skip_bbsplit: Boolean // Skip BBSplit for removal of non-reference genome reads.
    skip_fastqc: Boolean // true/false
    skip_trimming: Boolean // true/false
    skip_umi_extract: Boolean // true/false
    skip_linting: Boolean // true/false

    // Index generation
    make_salmon_index: Boolean // Whether to create salmon index before running salmon quant
    make_sortmerna_index: Boolean // Whether to create a sortmerna index before running sortmerna
    make_bowtie2_index: Boolean // Whether to create a bowtie2 index before running bowtie2

    // Trimming options
    trimmer: String // 'fastp' or 'trimgalore'
    min_trimmed_reads: Integer // > 0
    save_trimmed: Boolean // true/false
    fastp_merge: Boolean // true/false: whether to stitch paired end reads together in FASTP output

    // rRNA removal options
    remove_ribo_rna: Boolean // true/false: whether to remove rRNA
    ribo_removal_tool: String // 'sortmerna', 'ribodetector', or 'bowtie2'

    // UMI options
    with_umi: Boolean // true/false: Enable UMI-based read deduplication.
    umi_discard_read: Integer // 0, 1 or 2

    // Merging options
    save_merged_fastq: Boolean // true/false: Save merged FastQ files even for single-library samples

    // Strandedness thresholds
    stranded_threshold: Float // The fraction of stranded reads that must be assigned to a strandedness for confident assignment. Must be at least 0.5
    unstranded_threshold: Float // The difference in fraction of stranded reads assigned to 'forward' and 'reverse' below which a sample is classified as 'unstranded'

    // Tool options
    tool_args: Record // fq_lint, plus umi_extract, fastp, trimgalore, use_gpu_ribodetector and salmon_index_strandedness for the subworkflows it runs

    main:

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    // `reads` is flat, so a paired-end sample holds two files per run: one run means one file (single-end) or two (paired-end)
    ch_fastq_single = ch_reads.filter { r ->
        (r.meta.single_end ? r.reads.size() == 1 : r.reads.size() == 2) && r.reads[0].name.endsWith('.gz') && !save_merged_fastq
    }
    ch_fastq_multiple = ch_reads.filter { r ->
        !((r.meta.single_end ? r.reads.size() == 1 : r.reads.size() == 2) && r.reads[0].name.endsWith('.gz') && !save_merged_fastq)
    }
    ch_cat = CAT_FASTQ(ch_fastq_multiple.map { r -> record(meta: r.meta, reads: r.reads) })
        .mix(ch_fastq_single.map { r -> record(meta: r.meta, reads: r.reads) })

    // Each stage that runs joins its outputs onto this per-sample record, overwriting the
    // null placeholders of the fields it owns. `reads` always holds the reads handed to the
    // next stage and is null for samples that failed a filter. Stages after the
    // min_trimmed_reads filter run on samples with reads and join back with `remainder` so
    // that samples failing it keep their record. Joins key on meta, which the trimming
    // subworkflows can change (UMI read discard sets single_end), so their output replaces
    // the record instead of being joined onto it.
    ch_samples = ch_cat.map { r ->
        record(
            meta:                 r.meta,
            reads:                r.reads,
            reads_cat:            r.reads,
            reads_trimmed:        null,
            num_trimmed_reads:    null,
            lint_raw:             null,
            lint_trimmed:         null,
            lint_bbsplit:         null,
            lint_ribo:            null,
            fastqc_raw_html:      null,
            fastqc_raw_zip:       null,
            fastqc_trim_html:     null,
            fastqc_trim_zip:      null,
            fastqc_filtered_html: null,
            fastqc_filtered_zip:  null,
            trim:                 null,
            umi:                  null,
            bbsplit:              null,
            rrna:                 null
        )
    }

    //
    // MODULE: Lint FastQ files
    //

    if (!skip_linting) {
        ch_lint_raw = FQ_LINT(ch_cat.map { r -> record(meta: r.meta, reads: r.reads, args: tool_args.fq_lint) })
        ch_samples = ch_samples.join(ch_lint_raw.map { r -> record(meta: r.meta, lint_raw: r.lint) }, by: 'meta')
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    //
    if (trimmer == 'trimgalore') {
        ch_trimgalore = FASTQ_FASTQC_UMITOOLS_TRIMGALORE(
            ch_samples,
            skip_fastqc,
            with_umi,
            skip_umi_extract,
            skip_trimming,
            umi_discard_read,
            min_trimmed_reads,
            tool_args,
        )

        // TrimGalore's own html and zip are FastQC reports on the trimmed reads
        ch_samples = ch_trimgalore.map { r ->
            r + record(
                fastqc_trim_html: r.trim != null ? r.trim.html : null,
                fastqc_trim_zip:  r.trim != null ? r.trim.zip : null,
                trim: r.trim != null
                    ? record(html: null, log: r.trim.log, json: r.trim.json, unpaired: r.trim.unpaired, reads_fail: null, reads_merged: null)
                    : null
            )
        }
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (trimmer == 'fastp') {
        ch_fastp = FASTQ_FASTQC_UMITOOLS_FASTP(
            ch_samples.map { r -> r + record(adapter_fasta: null) }, // No adapter fasta
            skip_fastqc,
            with_umi,
            skip_umi_extract,
            umi_discard_read,
            skip_trimming,
            save_trimmed,
            fastp_merge,
            min_trimmed_reads,
            tool_args,
        )

        ch_samples = ch_fastp.map { r ->
            r + record(
                trim: r.trim != null
                    ? record(html: [r.trim.html], log: [r.trim.log], json: [r.trim.json], unpaired: null, reads_fail: r.trim.reads_fail, reads_merged: r.trim.reads_merged)
                    : null
            )
        }
    }

    if (!skip_linting && !skip_trimming) {
        ch_lint_trimmed = FQ_LINT_AFTER_TRIMMING(ch_samples.filter { r -> r.reads != null }.map { r -> record(meta: r.meta, reads: r.reads, args: tool_args.fq_lint) })
        ch_samples = ch_samples.join(
            ch_lint_trimmed.map { r -> record(meta: r.meta, lint_trimmed: r.lint) },
            by: 'meta',
            remainder: true
        )
    }

    ch_samples = ch_samples.map { r -> r + record(reads_trimmed: r.reads) }

    //
    // MODULE: Remove genome contaminant reads
    //
    if (!skip_bbsplit) {
        ch_bbsplit = BBMAP_BBSPLIT(
            ch_samples
                .filter { r -> r.reads != null }
                .map { r -> record(meta: r.meta, reads: r.reads, primary_ref: null, other_ref_names: [], other_ref_paths: [], only_build_index: false) }
                .combine(index: ch_bbsplit_index)
        )

        ch_samples = ch_samples.join(
            ch_bbsplit.map { r ->
                record(
                    meta:    r.meta,
                    reads:   r.reads.isEmpty() ? null : r.reads,
                    bbsplit: r.stats != null
                        ? record(
                            stats:              r.stats,
                            primary_reads:      r.reads.isEmpty() ? null : r.reads,
                            other_genome_reads: r.reads.isEmpty() ? null : r.other_genome_reads
                        )
                        : null
                )
            },
            by: 'meta',
            remainder: true
        )

        if (!skip_linting) {
            ch_lint_bbsplit = FQ_LINT_AFTER_BBSPLIT(ch_samples.filter { r -> r.reads != null }.map { r -> record(meta: r.meta, reads: r.reads, args: tool_args.fq_lint) })
            ch_samples = ch_samples.join(
                ch_lint_bbsplit.map { r -> record(meta: r.meta, lint_bbsplit: r.lint) },
                by: 'meta',
                remainder: true
            )
        }
    }

    //
    // SUBWORKFLOW: Remove ribosomal RNA reads
    //
    def val_rrna_references: Value<Record> = channel.value(
        record(sortmerna_index: null, bowtie2_index: null, seqkit_prefixed: null, seqkit_converted: null)
    )
    if (remove_ribo_rna) {
        ch_rrna_removed = FASTQ_REMOVE_RRNA(
            ch_samples.filter { r -> r.reads != null },
            ch_rrna_fastas,
            ch_sortmerna_index,
            ch_bowtie2_index,
            ribo_removal_tool,
            make_sortmerna_index,
            make_bowtie2_index,
            tool_args,
        )

        val_rrna_references = ch_rrna_removed.references

        ch_samples = ch_samples.join(
            ch_rrna_removed.samples.map { r ->
                record(
                    meta:  r.meta,
                    reads: r.reads,
                    rrna:  record(
                        sortmerna_log:    r.sortmerna_log,
                        ribodetector_log: r.ribodetector_log,
                        seqkit_stats:     r.seqkit_stats,
                        bowtie2_log:      r.bowtie2_log
                    )
                )
            },
            by: 'meta',
            remainder: true
        )

        if (!skip_linting) {
            ch_lint_ribo = FQ_LINT_AFTER_RIBO_REMOVAL(ch_samples.filter { r -> r.reads != null }.map { r -> record(meta: r.meta, reads: r.reads, args: tool_args.fq_lint) })
            ch_samples = ch_samples.join(
                ch_lint_ribo.map { r -> record(meta: r.meta, lint_ribo: r.lint) },
                by: 'meta',
                remainder: true
            )
        }
    }

    //
    // MODULE: Run FastQC on filtered reads (after BBSplit and/or rRNA removal)
    //
    if (!skip_fastqc && (!skip_bbsplit || remove_ribo_rna)) {
        ch_fastqc_filtered = FASTQC_FILTERED(ch_samples.filter { r -> r.reads != null })
        ch_samples = ch_samples.join(
            ch_fastqc_filtered.map { r -> record(meta: r.meta, fastqc_filtered_html: r.html, fastqc_filtered_zip: r.zip) },
            by: 'meta',
            remainder: true
        )
    }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudoalign with Salmon to auto-infer strandedness
    //
    ch_auto_strand = ch_samples.filter { r -> r.reads != null && r.meta.strandedness == 'auto' }

    ch_salmon = FASTQ_SUBSAMPLE_FQ_SALMON(
        ch_auto_strand,
        ch_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        make_salmon_index,
        tool_args,
    )
    ch_lib_format_counts = ch_salmon.samples.map { r -> record(meta: r.meta, lib_format_counts: r.lib_format_counts) }
    ch_salmon_index_built = ch_salmon.index_built

    // Joined on meta before the inferred strandedness is added to it
    ch_inferred = ch_auto_strand
        .join(ch_lib_format_counts, by: 'meta', remainder: true)
        .map { r ->
            if (r.lib_format_counts == null) {
                error("Salmon failed to produce lib_format_counts for sample '${r.meta.id}' " +
                    "which was set to 'auto' strandedness. Check that the Salmon " +
                    "index matches your input reads, or set strandedness explicitly in the samplesheet.")
            }
            def salmon_strand_analysis = getSalmonInferredStrandedness(r.lib_format_counts, stranded_threshold, unstranded_threshold)
            def strandedness = salmon_strand_analysis.inferred_strandedness
            if (strandedness == 'undetermined') {
                strandedness = 'unstranded'
            }
            preprocessedSample(r, r.meta + record(strandedness: strandedness, salmon_strand_analysis: salmon_strand_analysis), skip_linting)
        }

    ch_known_strand = ch_samples
        .filter { r -> r.reads == null || r.meta.strandedness != 'auto' }
        .map { r -> preprocessedSample(r, r.meta, skip_linting) }

    ch_results = ch_inferred.mix(ch_known_strand)

    emit:
    samples: Channel<Record> = ch_results
    rrna_references: Value<Record> = val_rrna_references
    salmon_index_built: Value<Path?> = ch_salmon_index_built
}
