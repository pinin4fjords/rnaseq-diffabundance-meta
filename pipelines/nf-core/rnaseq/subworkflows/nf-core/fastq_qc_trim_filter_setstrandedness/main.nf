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

workflow FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS {
    take:
    // Input channels
    ch_reads: Channel<Tuple> // (meta, reads)
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
    ch_fastq_single = ch_reads.filter { meta, reads ->
        (meta.single_end ? reads.size() == 1 : reads.size() == 2) && reads[0].name.endsWith('.gz') && !save_merged_fastq
    }
    ch_fastq_multiple = ch_reads.filter { meta, reads ->
        !((meta.single_end ? reads.size() == 1 : reads.size() == 2) && reads[0].name.endsWith('.gz') && !save_merged_fastq)
    }
    // (meta, reads)
    ch_cat = CAT_FASTQ(ch_fastq_multiple).mix(ch_fastq_single)

    //
    // MODULE: Lint FastQ files
    //
    // (meta, reads_cat, lint_raw)
    ch_raw = ch_cat.map { meta, reads -> tuple(meta, reads, null) }
    if (!skip_linting) {
        ch_lint_raw = FQ_LINT(ch_cat.map { meta, reads -> tuple(meta, reads, tool_args.fq_lint) })
        ch_raw = ch_cat.join(ch_lint_raw, by: 0)
    }

    // Both trimmers are brought to one layout. Elements of stages that do not run are null.
    // (meta, reads, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi(log, reads),
    //  trim(html, log, json, unpaired, reads_fail, reads_merged), num_trimmed_reads)
    ch_trimmed = ch_cat.map { meta, reads -> tuple(meta, reads, null, null, null, null, null, null, null) }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with TrimGalore!
    //
    if (trimmer == 'trimgalore') {
        // (meta, reads, fastqc_raw_html, fastqc_raw_zip, umi(log, reads), trim(html, zip, log, json, unpaired), num_trimmed_reads)
        ch_trimgalore = FASTQ_FASTQC_UMITOOLS_TRIMGALORE(
            ch_cat,
            skip_fastqc,
            with_umi,
            skip_umi_extract,
            skip_trimming,
            umi_discard_read,
            min_trimmed_reads,
            tool_args,
        )

        // TrimGalore's own html and zip are FastQC reports on the trimmed reads
        ch_trimmed = ch_trimgalore.map { meta, reads, fastqc_raw_html, fastqc_raw_zip, umi, trim, num_trimmed_reads ->
            tuple(
                meta,
                reads,
                fastqc_raw_html,
                fastqc_raw_zip,
                trim != null ? trim[0] : null,
                trim != null ? trim[1] : null,
                umi,
                trim != null ? tuple(null, trim[2], trim[3], trim[4], null, null) : null,
                num_trimmed_reads
            )
        }
    }

    //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters with fastp
    //
    if (trimmer == 'fastp') {
        // (meta, reads, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi(log, reads),
        //  trim(html, json, log, reads_fail, reads_merged), adapter_seq, num_trimmed_reads)
        ch_fastp = FASTQ_FASTQC_UMITOOLS_FASTP(
            ch_cat.map { meta, reads -> tuple(meta, reads, null) }, // No adapter fasta
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

        ch_trimmed = ch_fastp.map { meta, reads, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi, trim, _adapter_seq, num_trimmed_reads ->
            tuple(
                meta,
                reads,
                fastqc_raw_html,
                fastqc_raw_zip,
                fastqc_trim_html,
                fastqc_trim_zip,
                umi,
                trim != null ? tuple([trim[0]], [trim[2]], [trim[1]], null, trim[3], trim[4]) : null,
                num_trimmed_reads
            )
        }
    }

    // The trimming subworkflows can change meta (UMI read discard sets single_end), so the
    // pre-trimming elements are joined back on the sample id rather than on meta.
    // (meta, reads_trimmed, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi, trim, num_trimmed_reads,
    //  reads_cat, lint_raw)
    ch_samples = ch_trimmed
        .map { meta, reads, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi, trim, num_trimmed_reads ->
            tuple(meta.id, meta, reads, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi, trim, num_trimmed_reads)
        }
        .join(ch_raw.map { meta, reads_cat, lint_raw -> tuple(meta.id, reads_cat, lint_raw) }, by: 0)
        .map { _id, meta, reads, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi, trim, num_trimmed_reads, reads_cat, lint_raw ->
            tuple(meta, reads, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi, trim, num_trimmed_reads, reads_cat, lint_raw)
        }

    // Each later stage runs on the samples that still have reads and emits its own channel, which
    // is joined onto ch_samples at the end with `remainder` so that samples it did not run on keep
    // null elements. ch_reads_current holds every sample with the reads handed to the next stage,
    // null for samples that failed a filter.
    // (meta, reads)
    ch_reads_current = ch_trimmed.map { meta, reads, _fastqc_raw_html, _fastqc_raw_zip, _fastqc_trim_html, _fastqc_trim_zip, _umi, _trim, _num_trimmed_reads ->
        tuple(meta, reads)
    }

    // (meta, lint)
    def ch_lint_trimmed: Channel<Tuple> = channel.empty()
    if (!skip_linting && !skip_trimming) {
        ch_lint_trimmed = FQ_LINT_AFTER_TRIMMING(
            ch_reads_current
                .filter { _meta, reads -> reads != null }
                .map { meta, reads -> tuple(meta, reads, tool_args.fq_lint) }
        )
    }

    //
    // MODULE: Remove genome contaminant reads
    //
    // (meta, reads, bbsplit(stats, primary_reads, other_genome_reads))
    def ch_bbsplit: Channel<Tuple> = channel.empty()
    // (meta, lint)
    def ch_lint_bbsplit: Channel<Tuple> = channel.empty()
    if (!skip_bbsplit) {
        // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
        ch_bbsplit = BBMAP_BBSPLIT(
            ch_reads_current
                .filter { _meta, reads -> reads != null }
                .map { t -> tuple(t) }
                .combine(ch_bbsplit_index)
                .map { meta, reads, index -> tuple(meta, reads, index, null, [] as List<String>, [] as List<Path>, false) }
        ).map { meta, _index, reads, other_genome_reads, stats, _log ->
            tuple(
                meta,
                reads.isEmpty() ? null : reads,
                stats != null
                    ? tuple(stats, reads.isEmpty() ? null : reads, reads.isEmpty() ? null : other_genome_reads)
                    : null
            )
        }

        // BBSplit replaces the reads of the samples it ran on. Its two elements are nested so that
        // a sample it did not run on gets one null, told apart from a run that produced no reads.
        ch_reads_current = ch_reads_current
            .join(ch_bbsplit.map { meta, reads, bbsplit -> tuple(meta, tuple(reads, bbsplit)) }, by: 0, remainder: true)
            .map { meta, reads, bbsplit_out -> tuple(meta, bbsplit_out != null ? bbsplit_out[0] : reads) }

        if (!skip_linting) {
            ch_lint_bbsplit = FQ_LINT_AFTER_BBSPLIT(
                ch_reads_current
                    .filter { _meta, reads -> reads != null }
                    .map { meta, reads -> tuple(meta, reads, tool_args.fq_lint) }
            )
        }
    }

    //
    // SUBWORKFLOW: Remove ribosomal RNA reads
    //
    // (sortmerna_index, bowtie2_index, seqkit_prefixed, seqkit_converted)
    def val_rrna_references: Value<Tuple> = channel.value(tuple(null, null, null, null))
    // (meta, reads, rrna(sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log))
    def ch_rrna: Channel<Tuple> = channel.empty()
    // (meta, lint)
    def ch_lint_ribo: Channel<Tuple> = channel.empty()
    if (remove_ribo_rna) {
        ch_rrna_removed = FASTQ_REMOVE_RRNA(
            ch_reads_current.filter { _meta, reads -> reads != null },
            ch_rrna_fastas,
            ch_sortmerna_index,
            ch_bowtie2_index,
            ribo_removal_tool,
            make_sortmerna_index,
            make_bowtie2_index,
            tool_args,
        )

        val_rrna_references = ch_rrna_removed.references

        ch_rrna = ch_rrna_removed.samples.map { meta, reads, sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log ->
            tuple(meta, reads, tuple(sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log))
        }

        ch_reads_current = ch_reads_current
            .join(ch_rrna.map { meta, reads, rrna -> tuple(meta, tuple(reads, rrna)) }, by: 0, remainder: true)
            .map { meta, reads, rrna_out -> tuple(meta, rrna_out != null ? rrna_out[0] : reads) }

        if (!skip_linting) {
            ch_lint_ribo = FQ_LINT_AFTER_RIBO_REMOVAL(
                ch_reads_current
                    .filter { _meta, reads -> reads != null }
                    .map { meta, reads -> tuple(meta, reads, tool_args.fq_lint) }
            )
        }
    }

    //
    // MODULE: Run FastQC on filtered reads (after BBSplit and/or rRNA removal)
    //
    // (meta, fastqc_filtered(html, zip))
    def ch_fastqc_filtered: Channel<Tuple> = channel.empty()
    if (!skip_fastqc && (!skip_bbsplit || remove_ribo_rna)) {
        ch_fastqc_filtered = FASTQC_FILTERED(ch_reads_current.filter { _meta, reads -> reads != null })
            .map { meta, html, zip -> tuple(meta, tuple(html, zip)) }
    }

    //
    // SUBWORKFLOW: Sub-sample FastQ files and pseudoalign with Salmon to auto-infer strandedness
    //
    // (meta, reads)
    ch_auto_strand = ch_reads_current.filter { meta, reads -> reads != null && meta.strandedness == 'auto' }

    ch_salmon = FASTQ_SUBSAMPLE_FQ_SALMON(
        ch_auto_strand,
        ch_fasta,
        ch_transcript_fasta,
        ch_gtf,
        ch_salmon_index,
        make_salmon_index,
        tool_args,
    )
    // (meta, lib_format_counts)
    ch_lib_format_counts = ch_salmon.samples.map { meta, _reads, _quant_dir, _json_info, lib_format_counts -> tuple(meta, lib_format_counts) }
    ch_salmon_index_built = ch_salmon.index_built

    // Joined on meta before the inferred strandedness is added to it
    // (meta, inferred_meta)
    ch_inferred_meta = ch_auto_strand
        .join(ch_lib_format_counts, by: 0, remainder: true)
        .map { meta, _reads, lib_format_counts ->
            if (lib_format_counts == null) {
                error("Salmon failed to produce lib_format_counts for sample '${meta.id}' " +
                    "which was set to 'auto' strandedness. Check that the Salmon " +
                    "index matches your input reads, or set strandedness explicitly in the samplesheet.")
            }
            def salmon_strand_analysis = getSalmonInferredStrandedness(lib_format_counts, stranded_threshold, unstranded_threshold)
            def strandedness = salmon_strand_analysis.inferred_strandedness
            if (strandedness == 'undetermined') {
                strandedness = 'unstranded'
            }
            tuple(meta, meta + [strandedness: strandedness, salmon_strand_analysis: salmon_strand_analysis])
        }

    // Every stage channel holds a subset of the samples in ch_samples under the same meta, so only
    // the stage side of each remainder join can be missing, and each stage side is one element.
    // (meta, reads_trimmed, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi, trim, num_trimmed_reads,
    //  reads_cat, lint_raw, reads, lint_trimmed, bbsplit, lint_bbsplit, rrna, lint_ribo, fastqc_filtered(html, zip), inferred_meta)
    ch_joined = ch_samples
        .join(ch_reads_current, by: 0)
        .join(ch_lint_trimmed, by: 0, remainder: true)
        .join(ch_bbsplit.map { meta, _reads, bbsplit -> tuple(meta, bbsplit) }, by: 0, remainder: true)
        .join(ch_lint_bbsplit, by: 0, remainder: true)
        .join(ch_rrna.map { meta, _reads, rrna -> tuple(meta, rrna) }, by: 0, remainder: true)
        .join(ch_lint_ribo, by: 0, remainder: true)
        .join(ch_fastqc_filtered, by: 0, remainder: true)
        .join(ch_inferred_meta, by: 0, remainder: true)

    // (meta, reads, reads_cat, reads_trimmed, num_trimmed_reads,
    //  fastqc(raw_html, raw_zip, trim_html, trim_zip, filtered_html, filtered_zip),
    //  trim(html, log, json, unpaired, reads_fail, reads_merged), umi(log, reads),
    //  bbsplit(stats, primary_reads, other_genome_reads), lint(raw, trimmed, bbsplit, ribo),
    //  rrna(sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log))
    ch_results = ch_joined.map { meta, reads_trimmed, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi, trim, num_trimmed_reads,
                                 reads_cat, lint_raw, reads, lint_trimmed, bbsplit, lint_bbsplit, rrna, lint_ribo, fastqc_filtered, inferred_meta ->
        def fastqc_filtered_html = fastqc_filtered != null ? fastqc_filtered[0] : null
        def fastqc_filtered_zip = fastqc_filtered != null ? fastqc_filtered[1] : null
        def has_fastqc = [
            fastqc_raw_html, fastqc_raw_zip,
            fastqc_trim_html, fastqc_trim_zip,
            fastqc_filtered_html, fastqc_filtered_zip
        ].any { f -> f != null }
        tuple(
            inferred_meta ?: meta,
            reads,
            reads_cat,
            reads_trimmed,
            num_trimmed_reads,
            has_fastqc
                ? tuple(fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, fastqc_filtered_html, fastqc_filtered_zip)
                : null,
            trim,
            umi,
            bbsplit,
            skip_linting ? null : tuple(lint_raw, lint_trimmed, lint_bbsplit, lint_ribo),
            rrna
        )
    }

    emit:
    // (meta, reads, reads_cat, reads_trimmed, num_trimmed_reads,
    //  fastqc(raw_html, raw_zip, trim_html, trim_zip, filtered_html, filtered_zip),
    //  trim(html, log, json, unpaired, reads_fail, reads_merged), umi(log, reads),
    //  bbsplit(stats, primary_reads, other_genome_reads), lint(raw, trimmed, bbsplit, ribo),
    //  rrna(sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log))
    samples: Channel<Tuple> = ch_results
    rrna_references: Value<Tuple> = val_rrna_references // (sortmerna_index, bowtie2_index, seqkit_prefixed, seqkit_converted)
    salmon_index_built: Value<Path?> = ch_salmon_index_built
}
