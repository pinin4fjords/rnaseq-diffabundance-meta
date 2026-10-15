#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/rnaseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/rnaseq
    Website: https://nf-co.re/rnaseq
    Slack  : https://nfcore.slack.com/channels/rnaseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PIPELINE PARAMETERS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Parameters that nextflow.config or the conf/ files read take their default from
    nextflow.config, because the config is resolved before this block; the block declares
    only their type.
*/

// One row of the input samplesheet. The params block needs a record type to load the
// samplesheet with typed fields: a plain Record leaves every value a String.
record SampleRow {
    sample:            String
    fastq_1:           Path
    fastq_2:           Path?
    strandedness:      String
    seq_platform:      String?
    seq_center:        String?
    genome_bam:        Path?
    transcriptome_bam: Path?
    percent_mapped:    Float?
}

params {
    // Input/output options
    input:                       Channel<SampleRow>
    email:                       String? = null

    // Reference genome options
    genome:                      String? = null
    fasta:                       String? = getGenomeAttribute('fasta')
    gtf:                         String? = getGenomeAttribute('gtf')
    gff:                         String? = getGenomeAttribute('gff')
    gene_bed:                    String? = getGenomeAttribute('bed12')
    transcript_fasta:            String? = getGenomeAttribute('transcript_fasta')
    additional_fasta:            String? = getGenomeAttribute('additional_fasta')
    splicesites:                 String? = null
    star_index:                  String? = getGenomeAttribute('star')
    hisat2_index:                String? = getGenomeAttribute('hisat2')
    rsem_index:                  String? = getGenomeAttribute('rsem')
    salmon_index:                String? = getGenomeAttribute('salmon')
    kallisto_index:              String? = getGenomeAttribute('kallisto')
    bowtie2_index:               String? = getGenomeAttribute('bowtie2')
    hisat2_build_memory:         String = '200.GB'
    gencode:                     Boolean = false
    prokaryotic:                 Boolean = false
    gffread_transcript_fasta:    Boolean = false
    gtf_extra_attributes:        String = 'gene_name'
    gtf_group_features:          String = 'gene_id'
    featurecounts_group_type:    String = 'gene_biotype'
    featurecounts_feature_type:  String = 'exon'
    skip_gtf_transcript_filter:  Boolean = false
    arm:                         Boolean = false
    // Read trimming options
    trimmer:                     String = 'trimgalore'
    min_trimmed_reads:           Integer = 10000
    extra_trimgalore_args:       String? = null
    extra_fastp_args:            String? = null
    extra_fqlint_args:           String = '--disable-validator P001'

    // Read filtering options
    bbsplit_fasta_list:          String? = null
    bbsplit_index:               String? = getGenomeAttribute('bbsplit')
    sortmerna_index:             String? = getGenomeAttribute('sortmerna')
    remove_ribo_rna:             Boolean = false
    ribo_removal_tool:           String = 'sortmerna'
    bowtie2_rrna_index:          String? = null
    ribo_database_manifest:      String = "${moduleDir}/assets/rrna-db-defaults.txt"
    use_gpu_ribodetector:        Boolean = false

    // UMI options
    with_umi:                    Boolean = false
    umitools_extract_method:     String = 'string'
    umitools_grouping_method:    String = 'directional'
    umitools_bc_pattern:         String?
    umitools_bc_pattern2:        String?
    umitools_umi_separator:      String?
    umi_dedup_tool:              String = 'umitools'
    umi_discard_read:            Integer = 0
    umitools_dedup_stats:        Boolean = false
    umitools_dedup_primary_only: Boolean = false

    // Alignment options
    aligner:                     String = 'star_salmon'
    use_sentieon_star:           Boolean = false
    use_parabricks_star:         Boolean = false
    pseudo_aligner:              String?
    bam_csi_index:               Boolean = false
    star_ignore_sjdbgtf:         Boolean = false
    pseudo_aligner_kmer_size:    Integer = 31
    salmon_quant_libtype:        String?
    extra_star_align_args:       String?
    extra_hisat2_align_args:     String?
    extra_bowtie2_align_args:    String?
    extra_salmon_index_args:     String?
    extra_salmon_quant_args:     String?
    extra_kallisto_quant_args:   String?
    min_mapped_reads:            Float = 5.0
    seq_center:                  String?
    seq_platform:                String?
    stringtie_ignore_gtf:        Boolean = false
    kallisto_quant_fraglen:      Integer = 200
    kallisto_quant_fraglen_sd:   Integer = 200
    stranded_threshold:          Float = 0.8
    unstranded_threshold:        Float = 0.1

    // Optional outputs
    save_merged_fastq:           Boolean = false
    save_umi_intermeds:          Boolean = false
    save_non_ribo_reads:         Boolean = false
    save_bbsplit_reads:          Boolean = false
    save_reference:              Boolean = false
    save_trimmed:                Boolean = false
    save_align_intermeds:        Boolean = false
    save_unaligned:              Boolean = false
    save_kraken_assignments:     Boolean = false
    save_kraken_unassigned:      Boolean = false

    // Quality Control
    rseqc_modules:               String = 'bam_stat,inner_distance,infer_experiment,junction_annotation,junction_saturation,read_distribution,read_duplication'
    contaminant_screening:       String?
    contaminant_screening_input: String = 'unmapped'
    bracken_precision:           String = 'S'
    deseq2_vst:                  Boolean = true
    kraken_db:                   String? = null
    sylph_db:                    Path? = null
    sylph_taxonomy:              Path? = null

    // Process skipping options
    skip_gtf_filter:             Boolean = false
    skip_bbsplit:                Boolean = true
    skip_umi_extract:            Boolean = false
    skip_linting:                Boolean = false
    skip_trimming:               Boolean = false
    skip_alignment:              Boolean = false
    skip_pseudo_alignment:       Boolean = false
    skip_quantification_merge:   Boolean = false
    skip_markduplicates:         Boolean = false
    skip_bigwig:                 Boolean = false
    skip_stringtie:              Boolean = false
    skip_fastqc:                 Boolean = false
    use_rustqc:                  Boolean = false
    skip_preseq:                 Boolean = true
    skip_dupradar:               Boolean = false
    skip_qualimap:               Boolean = false
    skip_rseqc:                  Boolean = false
    skip_biotype_qc:             Boolean = false
    skip_deseq2_qc:              Boolean = false
    skip_multiqc:                Boolean = false
    skip_qc:                     Boolean = false

    // Generic options
    version:                     Boolean = false
    email_on_fail:               String? = null
    plaintext_email:             Boolean = false
    monochrome_logs:             Boolean = false
    multiqc_config:              Path? = null
    multiqc_logo:                Path? = null
    multiqc_methods_description: Path? = null
    multiqc_title:               String?
    validate_params:             Boolean = true
    help:                        String? = null
    help_full:                   Boolean = false
    show_hidden:                 Boolean = false
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RNASEQ                     } from './workflows/rnaseq'
include { PREPARE_GENOME_REFERENCES  } from './subworkflows/local/prepare_genome_references'
include { PREPARE_GENOME_INDICES     } from './subworkflows/local/prepare_genome_indices'
include { PIPELINE_INITIALISATION    } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { PIPELINE_COMPLETION        } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { checkMaxContigSize         } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { defineQcTools              } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getGenomeAttribute         } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { isStarIndexLegacy          } from './subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { toolArgs } from './subworkflows/local/utils_nfcore_rnaseq_pipeline/tool_args'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline
//
workflow NFCORE_RNASEQ {

    take:
    params: Record   // the pipeline's params, passed in because an including pipeline has its own

    main:

    //
    // SUBWORKFLOW: Prepare reference genome files (FASTA, GTF, BED, transcript FASTA, chrom.sizes, rRNA FASTAs, Kraken DB)
    //
    def references_ribo_tool: String? = params.remove_ribo_rna && !(params.ribo_removal_tool == "bowtie2" && params.bowtie2_rrna_index) ? params.ribo_removal_tool as String : null
    def indices_ribo_tool: String? = params.remove_ribo_rna ? params.ribo_removal_tool as String : null

    references = PREPARE_GENOME_REFERENCES (
        params.fasta,
        params.gtf,
        params.gff,
        params.additional_fasta,
        params.transcript_fasta,
        params.gene_bed,
        params.ribo_database_manifest,
        params.kraken_db,
        params.gencode,
        params.gffread_transcript_fasta,
        params.featurecounts_group_type,
        params.aligner,
        params.pseudo_aligner,
        params.skip_gtf_filter,
        references_ribo_tool,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.use_sentieon_star,
        params.contaminant_screening,
        params.prokaryotic,
        params.skip_gtf_transcript_filter,
        params.genome
    )

    //
    // SUBWORKFLOW: Build or load aligner / pseudo-aligner / filtering indices
    //
    def tool_args: Record = toolArgs(params)

    indices = PREPARE_GENOME_INDICES (
        references.fasta,
        references.gtf,
        references.transcript_fasta,
        references.rrna_fastas,
        params.fasta ? true : false,
        params.splicesites,
        params.bbsplit_fasta_list,
        params.star_index,
        params.rsem_index,
        params.salmon_index,
        params.kallisto_index,
        params.hisat2_index,
        params.bowtie2_index,
        params.bbsplit_index,
        params.sortmerna_index,
        params.bowtie2_rrna_index,
        params.aligner,
        params.pseudo_aligner,
        params.skip_bbsplit,
        indices_ribo_tool,
        params.skip_alignment,
        params.skip_pseudo_alignment,
        params.use_sentieon_star,
        params.use_parabricks_star,
        isStarIndexLegacy(params) ? true : false,
        params.hisat2_build_memory,
        params.input.collect().map { rows -> rows.any { row -> row.strandedness == 'auto' } },
        params.prokaryotic,
        tool_args
    )

    // Check if contigs in genome fasta file > 512 Mbp
    if (!params.skip_alignment && !params.bam_csi_index) {
        references.fai.map { fai -> fai != null ? checkMaxContigSize(fai) : null }
    }

    //
    // WORKFLOW: Run nf-core/rnaseq workflow
    //
    qc_tools = defineQcTools(params)

    results = RNASEQ (
        params,
        tool_args,
        params.input,
        references.fasta,
        references.fai,
        references.gtf,
        references.chrom_sizes,
        references.gene_bed,
        references.transcript_fasta,
        indices.star_index,
        indices.rsem_index,
        indices.hisat2_index,
        indices.bowtie2_index,
        indices.salmon_index,
        indices.kallisto_index,
        indices.bbsplit_index,
        references.rrna_fastas,
        indices.sortmerna_index,
        indices.bowtie2_rrna_index,
        indices.splicesites,
        references.kraken_db,
        qc_tools
    )

    // No task-workdir guard: a user-supplied --bowtie2_rrna_index is republished too.
    def ch_rrna_bowtie2_index: Channel<Path> = results.rrna_references
        .map { _sortmerna_index, bowtie2_index, _seqkit_prefixed, _seqkit_converted -> bowtie2_index }
        .flatMap { index -> index == null ? [] : [index] }

    // Same-basename files split out per stage to avoid a >> rename-key collision (nextflow-io/nextflow#6617).
    // Each is (meta, lint_report), from lint(raw, trimmed, bbsplit, ribo) of the preprocessed tuple.
    def ch_lint_raw     = results.preprocessed
        .map { meta, _reads, _reads_cat, _reads_trimmed, _num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, lint, _rrna -> tuple(meta, nth(lint, 0)) }
        .filter { _meta, lint_report -> lint_report != null }
    def ch_lint_trimmed = results.preprocessed
        .map { meta, _reads, _reads_cat, _reads_trimmed, _num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, lint, _rrna -> tuple(meta, nth(lint, 1)) }
        .filter { _meta, lint_report -> lint_report != null }
    def ch_lint_bbsplit = results.preprocessed
        .map { meta, _reads, _reads_cat, _reads_trimmed, _num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, lint, _rrna -> tuple(meta, nth(lint, 2)) }
        .filter { _meta, lint_report -> lint_report != null }
    def ch_lint_ribo    = results.preprocessed
        .map { meta, _reads, _reads_cat, _reads_trimmed, _num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, lint, _rrna -> tuple(meta, nth(lint, 3)) }
        .filter { _meta, lint_report -> lint_report != null }

    // The prepared rRNA FASTAs publish whether or not --save_reference is set, so they cannot ride the genome target.
    // (sortmerna_index, bowtie2_index, seqkit_prefixed, seqkit_converted)
    def ch_rrna_seqkit = results.rrna_references.flatMap { sortmerna_index, bowtie2_index, seqkit_prefixed, seqkit_converted ->
        ((seqkit_prefixed ?: []) as List<Path>).isEmpty() && ((seqkit_converted ?: []) as List<Path>).isEmpty()
            ? []
            : [tuple(sortmerna_index, bowtie2_index, seqkit_prefixed, seqkit_converted)]
    }

    // samplesheet_with_bams.csv rows: one per sequencing run of a sample, with the aligned
    // meta so an inferred strandedness replaces 'auto'. genome_bam is the
    // coordinate-sorted BAM; bowtie2_salmon aligns to the transcriptome, so its unsorted bowtie2 BAM is the transcriptome_bam.
    // The BAMs are published by the `aligned` output, so the row carries their published location
    // as a string; routing the Paths through `samplesheet` too would copy each file twice.
    // The runs carry the samplesheet meta, before strandedness inference, so they join on the sample name.
    // (id, runs)
    ch_runs = results.reads.map { meta, runs -> tuple(meta.id, runs) }

    // (sample, fastq_1, fastq_2, strandedness, seq_platform, seq_center, genome_bam, percent_mapped, transcriptome_bam):
    // the samplesheet output's `header` lists these column names in this order
    def ch_samplesheet_rows = results.aligned
        .map { meta, bam, _bai, _samtools, percent_mapped, _raw_bams, transcriptome_bam, _unmapped, _logs ->
            tuple(meta.id, meta, bam, percent_mapped, transcriptome_bam)
        }
        .join(ch_runs, by: 0)
        .flatMap { _id, meta, bam, percent_mapped, transcriptome_bam, runs ->
            def aligned_dir = "${params.skip_quantification_merge ? "${meta.id}/" : ''}${params.aligner}/"
            runs.collect { run ->
                tuple(
                    meta.id,
                    run[0],
                    run.size() > 1 ? run[1] : null,
                    meta.strandedness,
                    meta.seq_platform ?: params.seq_platform,
                    meta.seq_center ?: params.seq_center,
                    bam != null ? "${workflow.outputDir}/${aligned_dir}${bam.name}" : null,
                    percent_mapped,
                    transcriptome_bam != null ? "${workflow.outputDir}/${aligned_dir}${transcriptome_bam.name}" : null
                )
            }
        }

    // Merged gene matrices of the primary quantifier: the alignment-based one unless alignment is skipped
    def ch_gene_quant = params.skip_alignment ? results.quant_merged_pseudo : results.quant_merged

    emit:
    trim_status:           Channel<Tuple>   = results.trim_status     // (meta, pass)
    map_status:            Channel<Tuple>   = results.map_status      // (meta, pass)
    strand_status:         Channel<Tuple>   = results.strand_status   // (meta, pass)
    multiqc_report:        Channel<Path>    = results.multiqc_report
    gtf:                   Value<Path?>     = references.gtf
    // (meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, tpm_transcript,
    //  counts_transcript, lengths_transcript, tx2gene_augmented, merged_gene_rds, merged_transcript_rds, tx2gene)
    gene_quant:            Channel<Tuple>   = ch_gene_quant
    genome_references:     Channel<Tuple>   = references.references     // (kind, file)
    genome_intermediates:  Channel<Tuple>   = references.intermediates  // (kind, file)
    genome_indices:        Channel<Tuple>   = indices.indices           // (kind, file)
    rrna_bowtie2_index:    Channel<Path>    = ch_rrna_bowtie2_index
    rrna_seqkit:           Channel<Tuple>   = ch_rrna_seqkit

    // Per-sample stage results, keyed on meta. The RNASEQ workflow documents each layout.
    preprocessed:          Channel<Tuple>   = results.preprocessed
    lint_raw:              Channel<Tuple>   = ch_lint_raw
    lint_trimmed:          Channel<Tuple>   = ch_lint_trimmed
    lint_bbsplit:          Channel<Tuple>   = ch_lint_bbsplit
    lint_ribo:             Channel<Tuple>   = ch_lint_ribo
    aligned:               Channel<Tuple>   = results.aligned
    umi_dedup:             Channel<Tuple>   = results.umi_dedup
    markdup:               Channel<Tuple>   = results.markdup
    bam_qc:                Channel<Tuple>   = results.bam_qc
    bam_qc_rustqc:         Channel<Tuple>   = results.bam_qc_rustqc
    samplesheet:           Channel<Tuple>   = ch_samplesheet_rows
    quant:                 Channel<Tuple>   = results.quant
    quant_salmon:          Channel<Tuple>   = results.quant_salmon
    quant_merged:          Channel<Tuple>   = results.quant_merged
    quant_rsem_merge:      Channel<Tuple>   = results.quant_rsem_merge
    quant_pseudo:          Channel<Tuple>   = results.quant_pseudo
    quant_pseudo_kallisto: Channel<Tuple>   = results.quant_pseudo_kallisto
    quant_merged_pseudo:   Channel<Tuple>   = results.quant_merged_pseudo
    contaminants:          Channel<Tuple>   = results.contaminants
    stringtie:             Channel<Tuple>   = results.stringtie
    bigwig:                Channel<Tuple>   = results.bigwig

    // Run-level results
    stringtie_merged:      Channel<Tuple>   = results.stringtie_merged
    deseq2:                Channel<Tuple>   = results.deseq2
    deseq2_pseudo:         Channel<Tuple>   = results.deseq2_pseudo
    multiqc:               Channel<Tuple>   = results.multiqc          // (meta, report, data, plots)
    pipeline_info:         Channel<Path>    = results.pipeline_info
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden,
        params
    )

    //
    // WORKFLOW: Run main workflow
    //
    results = NFCORE_RNASEQ (params)

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.monochrome_logs,
        results.multiqc_report,
        results.trim_status,
        results.map_status,
        results.strand_status,
        params
    )

    publish:
    contaminants         = results.contaminants
    stringtie            = results.stringtie
    stringtie_merged     = results.stringtie_merged
    bigwig               = results.bigwig
    gtf                  = results.gtf
    gene_quant           = results.gene_quant
    genome_references    = results.genome_references
    genome_intermediates = results.genome_intermediates
    genome_indices       = results.genome_indices
    rrna_bowtie2_index   = results.rrna_bowtie2_index
    rrna_seqkit          = results.rrna_seqkit
    preprocessed         = results.preprocessed
    lint_raw             = results.lint_raw
    lint_trimmed         = results.lint_trimmed
    lint_bbsplit         = results.lint_bbsplit
    lint_ribo            = results.lint_ribo
    aligned              = results.aligned
    umi_dedup            = results.umi_dedup
    markdup              = results.markdup
    samplesheet          = results.samplesheet
    quant                = results.quant
    quant_salmon         = results.quant_salmon
    quant_merged         = results.quant_merged
    quant_rsem_merge     = results.quant_rsem_merge
    quant_pseudo         = results.quant_pseudo
    quant_pseudo_kallisto = results.quant_pseudo_kallisto
    quant_merged_pseudo  = results.quant_merged_pseudo
    deseq2               = results.deseq2
    deseq2_pseudo        = results.deseq2_pseudo
    bam_qc               = results.bam_qc
    bam_qc_rustqc        = results.bam_qc_rustqc
    multiqc              = results.multiqc
    pipeline_info        = results.pipeline_info
}

// Per-sample directory prefix for per-sample outputs under --skip_quantification_merge.
// Run-level results (e.g. a cross-sample merged file) are never sample-prefixed.
def samplePrefix(id: String) -> String { params.skip_quantification_merge ? "${id}/" : '' }

// Directory the aligned BAM/BAI/transcriptome BAM publish to. Also used to rebuild the
// samplesheet's genome_bam/transcriptome_bam paths in the entry workflow, so the two can't drift.
def alignedDir(id: String) -> String { "${samplePrefix(id)}${params.aligner}/" }

// Element `i` of a nested result tuple, or null when the tuple is null because its stage did not run.
def nth(t, i: Integer) { t != null ? (t as List)[i] : null }

// Each `path` closure destructures the tuple layout documented on its emit, and returns a map of
// files to target directories. When a file is a key more than once (two elements can hold the same
// file), the last entry wins, including a null target that suppresses it, so the order of entries
// is part of the routing. Where elements can overlap, the closure says which entry has to come last.
output {
    contaminants: Channel<Tuple> {   // exactly one tool branch is populated per run
        enabled !params.skip_qc && params.contaminant_screening
        path { meta, kraken2, bracken, sylph_profile, sylphtax_taxprof ->
            // kraken2: (report, classified_reads_fastq, unclassified_reads_fastq, classified_reads_assignment); bracken: (abundance, report)
            def dir     = "${alignedDir(meta.id)}contaminants/"
            def kraken  = "${dir}kraken2/kraken_reports/"
            def brackenDir = "${dir}bracken/"
            def sylph   = "${dir}sylph/"
            def saveAssignments = params.save_kraken_assignments ? kraken : null
            [
                (nth(kraken2, 0)):                           kraken,
                (nth(kraken2, 1)):                           saveAssignments,
                (nth(kraken2, 2)):                           saveAssignments,
                (nth(kraken2, 3)):                           params.save_kraken_unassigned ? kraken : null,
                ([nth(bracken, 0), nth(bracken, 1)]):        brackenDir,
                ([sylph_profile, sylphtax_taxprof]):         sylph,
            ]
        }
    }

    stringtie: Channel<Tuple> {
        enabled !params.skip_stringtie
        path { meta, transcript_gtf, abundance, coverage_gtf, ballgown, denovo ->
            // denovo: (transcript_gtf, abundance, coverage_gtf, ballgown)
            def dir = "${alignedDir(meta.id)}stringtie/"
            [
                ([transcript_gtf, abundance, coverage_gtf, nth(denovo, 0)]): dir,
                (ballgown):                                                  dir,
            ]
        }
    }

    stringtie_merged: Channel<Tuple> {   // cross-sample, never sample-prefixed
        enabled !params.skip_stringtie && params.stringtie_ignore_gtf
        path { _meta, _gtf, _assemblies, merged_gtf -> merged_gtf >> "${params.aligner}/stringtie/" }
    }

    bigwig: Channel<Tuple> {   // bedgraph stays unrouted
        enabled !params.skip_bigwig
        path { meta, combined, forward, reverse ->
            // each track: (bedgraph, bigwig)
            [([nth(combined, 1), nth(forward, 1), nth(reverse, 1)]): "${alignedDir(meta.id)}bigwig/"]
        }
    }

    // The genome targets route on the kind element: 'kraken_db' (genome_references) and 'sortmerna' (genome_indices).
    // gtf and gene_quant are only returned to an including pipeline, never published
    gtf: Path? {
        enabled false
        path { f -> f >> "genome/" }
    }

    gene_quant: Channel<Tuple> {
        enabled false
        path { _meta, _tpm_gene, counts_gene, _lengths_gene, _counts_gene_length_scaled, _counts_gene_scaled, _tpm_transcript,
               _counts_transcript, _lengths_transcript, _tx2gene_augmented, _merged_gene_rds, _merged_transcript_rds, _tx2gene ->
            counts_gene >> "quant/"
        }
    }

    genome_references: Channel<Tuple> {   // never sample-prefixed
        enabled params.save_reference
        path { kind, ref_file -> ref_file >> (kind == 'kraken_db' ? 'genome/index/' : 'genome/') }
    }

    genome_intermediates: Channel<Tuple> {   // never sample-prefixed
        enabled params.save_reference
        path { _kind, ref_file -> ref_file >> 'genome/' }
    }

    genome_indices: Channel<Tuple> {   // never sample-prefixed
        enabled params.save_reference
        path { kind, index_file -> index_file >> (kind == 'sortmerna' ? 'genome/sortmerna/' : 'genome/index/') }
    }

    rrna_bowtie2_index: Channel<Path> {   // never sample-prefixed
        enabled params.save_reference
        path { p -> p >> 'bowtie2_rrna/index/' }
    }

    rrna_seqkit: Channel<Tuple> {   // ch_rrna_seqkit guarantees at least one FASTA list is non-empty
        path { _sortmerna_index, _bowtie2_index, seqkit_prefixed, seqkit_converted ->
            seqkit_prefixed >> 'seqkit/'
            seqkit_converted >> 'seqkit/'
        }
    }

    preprocessed: Channel<Tuple> {   // no single anchor element survives every skip combination
        path { meta, reads, reads_cat, reads_trimmed, _num_trimmed_reads, fastqc, trim, umi, bbsplit, _lint, rrna ->
            // fastqc: (raw_html, raw_zip, trim_html, trim_zip, filtered_html, filtered_zip)
            // trim: (html, log, json, unpaired, reads_fail, reads_merged); umi: (log, reads)
            // bbsplit: (stats, primary_reads, other_genome_reads); rrna: (sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log)
            def sp       = samplePrefix(meta.id)
            def fastqcDir = "${sp}fastqc/"
            def trimDir  = "${sp}${params.trimmer}/"
            def umitools = "${sp}umitools/"
            def bbsplitDir = "${sp}bbsplit/"
            def riboDir  = "${sp}${params.ribo_removal_tool == 'bowtie2' ? 'bowtie2_rrna' : params.ribo_removal_tool}/"
            def isFastp  = params.trimmer == 'fastp'
            def saveTrimmed = params.save_trimmed ? trimDir : null
            def saveBbsplit = params.save_bbsplit_reads ? bbsplitDir : null
            // Only the single-end --un-gz FASTQs are published; paired-end reads rebuilt by SAMTOOLS_FASTQ_BOWTIE2 are not.
            def saveNonRibo = params.remove_ribo_rna && params.save_non_ribo_reads && (params.ribo_removal_tool != 'bowtie2' || meta.single_end)
            [
                // `reads` can be the same files as `reads_cat` or `reads_trimmed`. A later entry for the same files
                // replaces an earlier one even when its target is null, so `reads` must stay first.
                (reads):                  saveNonRibo ? riboDir : null,
                (reads_cat):              params.save_merged_fastq ? "${sp}fastq/" : null,
                (reads_trimmed):          params.skip_trimming ? null : saveTrimmed,
                (nth(fastqc, 0)):         "${fastqcDir}raw/",
                (nth(fastqc, 1)):         "${fastqcDir}raw/",
                (nth(fastqc, 2)):         "${fastqcDir}trim/",
                (nth(fastqc, 3)):         "${fastqcDir}trim/",
                (nth(fastqc, 4)):         "${fastqcDir}filtered/",
                (nth(fastqc, 5)):         "${fastqcDir}filtered/",
                (nth(trim, 0)):           trimDir,
                (nth(trim, 1)):           isFastp ? "${trimDir}log/" : trimDir,
                (nth(trim, 2)):           isFastp ? trimDir : null,
                (nth(trim, 3)):           saveTrimmed,
                (nth(trim, 4)):           saveTrimmed,
                (nth(trim, 5)):           saveTrimmed,
                (nth(umi, 0)):            umitools,
                (nth(umi, 1)):            params.save_umi_intermeds ? umitools : null,
                (nth(bbsplit, 0)):        bbsplitDir,
                (nth(bbsplit, 1)):        saveBbsplit,
                (nth(bbsplit, 2)):        saveBbsplit,
                (nth(rrna, 0)):           "${sp}sortmerna/",
                ([nth(rrna, 1), nth(rrna, 2)]): "${sp}ribodetector/",
                (nth(rrna, 3)):           "${sp}bowtie2_rrna/",
            ]
        }
    }

    // One target per FQ_LINT stage: same basename, avoids the >> collision (nextflow-io/nextflow#6617).
    lint_raw: Channel<Tuple> {   // lint_report is non-null: ch_lint_raw filters nulls
        path { meta, lint_report -> lint_report >> "${samplePrefix(meta.id)}fq_lint/raw/" }
    }

    lint_trimmed: Channel<Tuple> {   // lint_report is non-null: ch_lint_trimmed filters nulls
        path { meta, lint_report -> lint_report >> "${samplePrefix(meta.id)}fq_lint/trimmed/" }
    }

    lint_bbsplit: Channel<Tuple> {   // lint_report is non-null: ch_lint_bbsplit filters nulls
        path { meta, lint_report -> lint_report >> "${samplePrefix(meta.id)}fq_lint/bbsplit/" }
    }

    lint_ribo: Channel<Tuple> {   // lint_report is non-null: ch_lint_ribo filters nulls
        path { meta, lint_report -> lint_report >> "${samplePrefix(meta.id)}fq_lint/${params.ribo_removal_tool ?: 'sortmerna'}/" }
    }

    aligned: Channel<Tuple> {
        path { meta, bam, bai, samtools, _percent_mapped, raw_bams, transcriptome_bam, unmapped, logs ->
            // samtools: (stats, flagstat, idxstats)
            def dir         = alignedDir(meta.id)
            def bamDir      = (params.save_align_intermeds || params.skip_markduplicates) ? dir : null
            def intermedDir = params.save_align_intermeds ? dir : null
            [
                ([nth(samtools, 0), nth(samtools, 1), nth(samtools, 2)]): "${dir}samtools_stats/",
                ([bam, bai]):                                             bamDir,
                ((raw_bams + [transcriptome_bam]).toSet()):               intermedDir,
                (unmapped):                                               params.save_unaligned ? "${dir}unmapped/" : null,
                (logs):                                                   "${dir}log/",
            ]
        }
    }

    samplesheet: Channel<Tuple> {
        enabled params.save_align_intermeds && !params.skip_alignment
        index {
            path 'samplesheets/samplesheet_with_bams.csv'
            header(['sample', 'fastq_1', 'fastq_2', 'strandedness', 'seq_platform', 'seq_center', 'genome_bam', 'percent_mapped', 'transcriptome_bam'])
        }
    }

    umi_dedup: Channel<Tuple> {
        path { meta, bam, bai, samtools, genomic_dedup_log, tsv, transcriptomic_dedup_log, prepare_for_rsem_log, _transcriptome_bam, t ->
            // samtools: (stats, flagstat, idxstats); tsv: (edit_distance, per_umi, umi_per_position)
            // t: (transcriptome_bam, dedup_bam, sorted_bam, sorted_bam_index, filtered_bam, samtools, tsv,
            //     coord_sorted_bam, coord_sorted_bam_index, coord_sorted_samtools)
            def dir      = alignedDir(meta.id)
            def stats    = "${dir}samtools_stats/"
            def umitools = "${dir}umitools/"
            def tool     = "${dir}${params.umi_dedup_tool}/"
            def saveBam  = params.save_align_intermeds || params.save_umi_intermeds
            def t_samtools = nth(t, 5)
            def t_tsv      = nth(t, 6)
            def t_coord_samtools = nth(t, 9)
            [
                ([nth(samtools, 0), nth(samtools, 1), nth(samtools, 2),
                  nth(t_samtools, 0), nth(t_samtools, 1), nth(t_samtools, 2)]): stats,
                ([bam, bai, nth(t, 7), nth(t, 8),
                  nth(t, 2), nth(t, 3), nth(t, 4), nth(t, 1)]): saveBam ? dir : null,
                ([nth(t_coord_samtools, 0), nth(t_coord_samtools, 1),
                  nth(t_coord_samtools, 2)]): saveBam ? stats : null,
                ([nth(tsv, 0), nth(tsv, 1), nth(tsv, 2),
                  nth(t_tsv, 0), nth(t_tsv, 1), nth(t_tsv, 2)]): umitools,
                (prepare_for_rsem_log):        "${umitools}prepare_for_quantification_log/",
                (genomic_dedup_log):           "${tool}genomic_dedup_log/",
                (transcriptomic_dedup_log):    "${tool}transcriptomic_dedup_log/",
            ]
        }
    }

    markdup: Channel<Tuple> {
        path { meta, bam, bai, cram, metrics, samtools ->
            // samtools: (stats, flagstat, idxstats), null when RustQC replaces the stats
            def dir = alignedDir(meta.id)
            [
                ([bam, cram, bai]):                                        dir,
                (metrics):                                                 "${dir}picard_metrics/",
                ([nth(samtools, 0), nth(samtools, 1), nth(samtools, 2)]): "${dir}samtools_stats/",
            ]
        }
    }

    quant: Channel<Tuple> {
        path { meta, counts_gene, counts_transcript, stat, log, _bam_star, _bam_genome, _bam_transcript ->
            def dir = alignedDir(meta.id)
            [
                ([counts_gene, counts_transcript, stat]): dir,
                (log):                                    "${dir}log/",
            ]
        }
    }

    quant_salmon: Channel<Tuple> {
        path { meta, quant_dir, _json_info, _lib_format_counts -> [(quant_dir): alignedDir(meta.id)] }
    }

    quant_merged: Channel<Tuple> {   // never sample-prefixed except under --skip_quantification_merge
        path { meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, tpm_transcript,
               counts_transcript, lengths_transcript, tx2gene_augmented, merged_gene_rds, merged_transcript_rds, tx2gene ->
            def top = "${params.aligner}/"
            def dir = "${samplePrefix(meta.id)}${top}"
            [
                ([tpm_gene, counts_gene, lengths_gene, counts_gene_scaled,
                  tpm_transcript, counts_transcript, lengths_transcript,
                  tx2gene_augmented, merged_gene_rds, merged_transcript_rds]): dir,
                (counts_gene_length_scaled): params.skip_quantification_merge ? null : top,
                (tx2gene):                   top,
            ]
        }
    }

    // CUSTOM_RSEMMERGECOUNTS and tximport both write rsem.merged.* basenames, so they need separate targets (nextflow-io/nextflow#6617).
    quant_rsem_merge: Channel<Tuple> {   // rsem_merge is non-null in every tuple
        path { meta, rsem_merge ->
            // rsem_merge: (counts_gene, tpm_gene, counts_transcript, tpm_transcript, genes_long, isoforms_long)
            [
                ([rsem_merge[0], rsem_merge[1], rsem_merge[2], rsem_merge[3], rsem_merge[4], rsem_merge[5]]): "${alignedDir(meta.id)}rsem_merge_counts/",
            ]
        }
    }

    quant_pseudo: Channel<Tuple> {
        path { meta, quant_dir, _json_info, _lib_format_counts -> quant_dir >> "${samplePrefix(meta.id)}${params.pseudo_aligner}/" }
    }

    // The Kallisto log lives inside quant_dir already and is not routed separately.
    quant_pseudo_kallisto: Channel<Tuple> {
        path { meta, quant_dir, _json_info, _log -> quant_dir >> "${samplePrefix(meta.id)}${params.pseudo_aligner}/" }
    }

    quant_merged_pseudo: Channel<Tuple> {
        path { meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, tpm_transcript,
               counts_transcript, lengths_transcript, tx2gene_augmented, merged_gene_rds, merged_transcript_rds, tx2gene ->
            def top = "${params.pseudo_aligner}/"
            def dir = "${samplePrefix(meta.id)}${top}"
            [
                ([tpm_gene, counts_gene, lengths_gene, counts_gene_scaled,
                  tpm_transcript, counts_transcript, lengths_transcript,
                  tx2gene_augmented, merged_gene_rds, merged_transcript_rds]): dir,
                (counts_gene_length_scaled): params.skip_quantification_merge ? null : top,
                (tx2gene):                   top,
            ]
        }
    }

    deseq2: Channel<Tuple> {   // never sample-prefixed
        path { _meta, rdata, pca_vals, plots_pdf, sample_dists, size_factors, log, _pca_multiqc, _dists_multiqc ->
            [([rdata, pca_vals, plots_pdf, sample_dists, size_factors, log]): "${params.aligner}/deseq2_qc/"]
        }
    }

    deseq2_pseudo: Channel<Tuple> {   // never sample-prefixed
        path { _meta, rdata, pca_vals, plots_pdf, sample_dists, size_factors, log, _pca_multiqc, _dists_multiqc ->
            [([rdata, pca_vals, plots_pdf, sample_dists, size_factors, log]): "${params.pseudo_aligner}/deseq2_qc/"]
        }
    }

    bam_qc: Channel<Tuple> {
        enabled defineQcTools(params).size() > 0   // same check that decides whether any of these tools ran
        path { meta, preseq, featurecounts, biotype, qualimap, dupradar, rseqc, _mqc_files ->
            // preseq: (lc_extrap, log); featurecounts: (counts, summary); biotype: (tsv, rrna)
            // dupradar: (scatter2d, boxplot, hist, dupmatrix, intercept_slope, multiqc, session_info)
            // rseqc: (bamstat, inferexperiment, innerdistance, junctionannotation, junctionsaturation,
            //         readdistribution, readduplication, tin), its nested groups laid out below
            def dir      = alignedDir(meta.id)
            def dupDir   = "${dir}dupradar/"
            def rseqcDir = "${dir}rseqc/"
            def annDir   = "${rseqcDir}junction_annotation/"
            def satDir   = "${rseqcDir}junction_saturation/"
            def rdupDir  = "${rseqcDir}read_duplication/"
            def innerDir = "${rseqcDir}inner_distance/"
            def inner    = nth(rseqc, 2)   // (distance, freq, mean, pdf, rscript)
            def ann      = nth(rseqc, 3)   // (bed, interact_bed, xls, pdf, events_pdf, rscript, log)
            def sat      = nth(rseqc, 4)   // (pdf, rscript)
            def rdup     = nth(rseqc, 6)   // (seq_xls, pos_xls, pdf, rscript)
            def tin      = nth(rseqc, 7)   // (txt, xls)
            [
                (nth(preseq, 0)):               "${dir}preseq/",
                (nth(preseq, 1)):               "${dir}preseq/log/",
                ([nth(featurecounts, 0), nth(featurecounts, 1), nth(biotype, 0), nth(biotype, 1)]): "${dir}featurecounts/",
                (qualimap):                     "${dir}qualimap/",
                (nth(dupradar, 0)):             "${dupDir}scatter_plot/",
                (nth(dupradar, 1)):             "${dupDir}box_plot/",
                (nth(dupradar, 2)):             "${dupDir}histogram/",
                (nth(dupradar, 3)):             "${dupDir}gene_data/",
                (nth(dupradar, 4)):             "${dupDir}intercepts_slope/",
                (nth(rseqc, 0)):                "${rseqcDir}bam_stat/",
                (nth(rseqc, 1)):                "${rseqcDir}infer_experiment/",
                ([nth(ann, 3), nth(ann, 4)]):   "${annDir}pdf/",
                ([nth(ann, 0), nth(ann, 1)]):   "${annDir}bed/",
                (nth(ann, 2)):                  "${annDir}xls/",
                (nth(ann, 6)):                  "${annDir}log/",
                (nth(ann, 5)):                  "${annDir}rscript/",
                (nth(sat, 0)):                  "${satDir}pdf/",
                (nth(sat, 1)):                  "${satDir}rscript/",
                (nth(rseqc, 5)):                "${rseqcDir}read_distribution/",
                (nth(rdup, 2)):                 "${rdupDir}pdf/",
                ([nth(rdup, 0), nth(rdup, 1)]): "${rdupDir}xls/",
                (nth(rdup, 3)):                 "${rdupDir}rscript/",
                ([nth(inner, 0), nth(inner, 1), nth(inner, 2)]): "${innerDir}txt/",
                (nth(inner, 3)):                "${innerDir}pdf/",
                (nth(inner, 4)):                "${innerDir}rscript/",
                ([nth(tin, 0), nth(tin, 1)]):   "${rseqcDir}tin/",
            ]
        }
    }

    bam_qc_rustqc: Channel<Tuple> {   // never sample-prefixed
        path { _meta, samtools, preseq_lc_extrap, dupradar, featurecounts, biotype, rseqc, qualimap, _all_files ->
            // samtools: (stats, flagstat, idxstats); dupradar: (scatter2d, boxplot, hist, dupmatrix, intercept_slope, multiqc)
            // featurecounts: (counts, summary); biotype: (tsv, mqc, rrna)
            // rseqc: (bamstat, inferexperiment, readdistribution, tin, innerdistance, junctionannotation,
            //         junctionsaturation, readduplication), its nested groups laid out below
            def dir      = "${params.aligner}/rustqc/"
            def dupDir   = "${dir}dupradar/"
            def fcDir    = "${dir}featurecounts/"
            def rseqcDir = "${dir}rseqc/"
            def annDir   = "${rseqcDir}junction_annotation/"
            def satDir   = "${rseqcDir}junction_saturation/"
            def rdupDir  = "${rseqcDir}read_duplication/"
            def innerDir = "${rseqcDir}inner_distance/"
            def tin      = nth(rseqc, 3)   // (txt, xls)
            def inner    = nth(rseqc, 4)   // (distance, freq, mean, summary, plot, rscript)
            def ann      = nth(rseqc, 5)   // (bed, interact_bed, xls, log, plot, rscript)
            def sat      = nth(rseqc, 6)   // (summary, plot, rscript)
            def rdup     = nth(rseqc, 7)   // (seq_xls, pos_xls, plot, rscript)
            def fcSummary = nth(featurecounts, 1) as Path?
            [
                ([nth(samtools, 0), nth(samtools, 1), nth(samtools, 2)]): "${dir}samtools_stats/",
                (preseq_lc_extrap):             "${dir}preseq/",
                (nth(dupradar, 0)):             "${dupDir}scatter_plot/",
                (nth(dupradar, 1)):             "${dupDir}box_plot/",
                (nth(dupradar, 2)):             "${dupDir}histogram/",
                (nth(dupradar, 3)):             "${dupDir}gene_data/",
                (nth(dupradar, 4)):             "${dupDir}intercepts_slope/",
                (nth(dupradar, 5)):             dupDir,
                ([nth(featurecounts, 0), nth(biotype, 0), nth(biotype, 1), nth(biotype, 2)]): fcDir,
                // Published under the name featureCounts gives its non-biotype summary
                (fcSummary):                    "${fcDir}${fcSummary?.name?.replace('.biotype.tsv.summary', '.tsv.summary')}",
                (qualimap):                     "${dir}qualimap/",
                (nth(rseqc, 0)):                "${rseqcDir}bam_stat/",
                (nth(rseqc, 1)):                "${rseqcDir}infer_experiment/",
                (nth(rseqc, 2)):                "${rseqcDir}read_distribution/",
                ([nth(tin, 0), nth(tin, 1)]):   "${rseqcDir}tin/",
                ([nth(ann, 0), nth(ann, 1)]):   "${annDir}bed/",
                (nth(ann, 2)):                  "${annDir}xls/",
                (nth(ann, 3)):                  "${annDir}log/",
                (nth(ann, 4)):                  "${annDir}plot/",
                (nth(ann, 5)):                  "${annDir}rscript/",
                (nth(sat, 0)):                  "${satDir}txt/",
                (nth(sat, 1)):                  "${satDir}plot/",
                (nth(sat, 2)):                  "${satDir}rscript/",
                ([nth(rdup, 0), nth(rdup, 1)]): "${rdupDir}xls/",
                (nth(rdup, 2)):                 "${rdupDir}plot/",
                (nth(rdup, 3)):                 "${rdupDir}rscript/",
                ([nth(inner, 0), nth(inner, 1), nth(inner, 2), nth(inner, 3)]): "${innerDir}txt/",
                (nth(inner, 4)):                "${innerDir}plot/",
                (nth(inner, 5)):                "${innerDir}rscript/",
            ]
        }
    }

    multiqc: Channel<Tuple> {
        path { meta, report, data, plots ->
            [([report, data, plots]): "${samplePrefix(meta.id)}multiqc${params.skip_alignment ? '' : "/${params.aligner}"}/"]
        }
    }

    pipeline_info: Channel<Path> {
        path { versions -> versions >> 'pipeline_info/' }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
