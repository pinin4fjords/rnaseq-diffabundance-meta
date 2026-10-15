nextflow.enable.types = true

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { DESEQ2_QC as DESEQ2_QC_BAM_SALMON } from '../../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_RSEM        } from '../../modules/local/deseq2_qc'
include { DESEQ2_QC as DESEQ2_QC_PSEUDO      } from '../../modules/local/deseq2_qc'
include { RUSTQC                              } from '../../modules/nf-core/rustqc/main'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { ALIGN_STAR                            } from '../../subworkflows/local/align_star'
include { ALIGN_BOWTIE2                         } from '../../subworkflows/local/align_bowtie2'
include { MULTIQC_RNASEQ                        } from '../../subworkflows/local/multiqc_rnaseq'
include { BAM_QC_RNASEQ                         } from '../../subworkflows/nf-core/bam_qc_rnaseq'
include { QUANTIFY_RSEM                         } from '../../subworkflows/nf-core/quantify_rsem'
include { BAM_DEDUP_UMI                         } from '../../subworkflows/nf-core/bam_dedup_umi'


include { deseq2QcArgs; rustqcArgs; starAlignArgs; hisat2AlignArgs; bowtie2AlignArgs } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline/tool_args'
include { readSamplesheet; samplesheetRowsToCsv } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline/samplesheet'
include { classifyStrand                 } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getHisat2PercentMapped         } from '../../subworkflows/local/utils_nfcore_rnaseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { STRINGTIE_STRINGTIE        } from '../../modules/nf-core/stringtie/stringtie'
include { KRAKEN2_KRAKEN2 as KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN as BRACKEN } from '../../modules/nf-core/bracken/bracken/main'
include { SYLPH_PROFILE              } from '../../modules/nf-core/sylph/profile/main'
include { SYLPHTAX_TAXPROF           } from '../../modules/nf-core/sylphtax/taxprof/main'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_FW          } from '../../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_REV         } from '../../modules/nf-core/bedtools/genomecov'
include { BEDTOOLS_GENOMECOV as BEDTOOLS_GENOMECOV_COMBINED    } from '../../modules/nf-core/bedtools/genomecov'
include { SAMTOOLS_INDEX                                       } from '../../modules/nf-core/samtools/index'

//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { softwareVersionsToYAML           } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { FASTQ_ALIGN_HISAT2               } from '../../subworkflows/nf-core/fastq_align_hisat2'
include { BAM_MARKDUPLICATES_PICARD        } from '../../subworkflows/nf-core/bam_markduplicates_picard'
include { BAM_STRINGTIE_MERGE              } from '../../subworkflows/nf-core/bam_stringtie_merge/main'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG as BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_COMBINED } from '../../subworkflows/nf-core/bedgraph_bedclip_bedgraphtobigwig'
include { QUANTIFY_PSEUDO_ALIGNMENT as QUANTIFY_BAM_SALMON } from '../../subworkflows/nf-core/quantify_pseudo_alignment'
include { QUANTIFY_PSEUDO_ALIGNMENT                         } from '../../subworkflows/nf-core/quantify_pseudo_alignment'
include { FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS              } from '../../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow RNASEQ {

    take:
    params: Record                                     // the pipeline's params
    tool_args: Record                                  // the tool arguments that follow from the params
    ch_sample_rows: Channel<Record>                 // one row per sequencing run of each sample
    ch_fasta: Value<Path?>                             // genome.fasta
    ch_fai: Value<Path?>                               // genome.fai
    ch_gtf: Value<Path?>                               // genome.gtf
    ch_chrom_sizes: Value<Path?>                       // genome.sizes
    ch_gene_bed: Value<Path?>                          // gene.bed
    ch_transcript_fasta: Value<Path?>                  // transcript.fasta
    ch_star_index: Value<Path?>                        // star/index/
    ch_rsem_index: Value<Path?>                        // rsem/index/
    ch_hisat2_index: Value<Path?>                      // hisat2/index/
    ch_bowtie2_index: Value<Path?>                     // bowtie2/index/ for alignment
    ch_salmon_index: Value<Path?>                      // salmon/index/
    ch_kallisto_index: Value<Path?>                    // kallisto/index/
    ch_bbsplit_index: Value<Path?>                     // bbsplit/index/
    ch_ribo_db: Channel<Path>                          // sortmerna_fasta_list
    ch_sortmerna_index: Value<Path?>                   // sortmerna/index/
    ch_bowtie2_rrna_index: Value<Path?>                // bowtie2/index/ for rRNA removal
    ch_splicesites: Value<Path?>                       // genome.splicesites.txt
    ch_kraken_db: Value<Path?>                         // kraken2/db/
    qc_tools: List<String>                             // QC tools to run, e.g. ['preseq', 'qualimap', 'rseqc_bam_stat', ...]

    main:

    def ch_pca_header_multiqc        = file("${moduleDir}/../../assets/deseq2_pca_header.txt", checkIfExists: true)
    def sample_status_header_multiqc = file("${moduleDir}/../../assets/sample_status_header.txt", checkIfExists: true)
    def ch_clustering_header_multiqc = file("${moduleDir}/../../assets/deseq2_clustering_header.txt", checkIfExists: true)
    def ch_biotypes_header_multiqc   = file("${moduleDir}/../../assets/biotypes_header.txt", checkIfExists: true)

    // Match the General Statistics column the active aligner emits so the
    // MultiQC fail_mapped row reads consistently with the rest of the report.
    def aligner_display_name = [
        'star_salmon'    : 'STAR uniquely mapped reads',
        'star_rsem'      : 'STAR uniquely mapped reads',
        'hisat2'         : 'HISAT2 overall alignment rate',
        'bowtie2_salmon' : 'Bowtie2 overall alignment rate',
    ][params.aligner as String] ?: 'Aligned reads'

    // Files each stage contributes to MultiQC, per sample, as (meta, files). The report-only and
    // sample-only channels hold the files that go to just the merged report or just the per-sample
    // reports (--skip_quantification_merge). fail_* rows are appended inside MULTIQC_RNASEQ.
    def ch_mqc_files: Channel<Tuple> = channel.empty()
    def ch_mqc_sample_only: Channel<Tuple> = channel.empty()
    def ch_mqc_report_only: Channel<Path> = channel.empty()

    // The rows are validated and merged per sample as one batch, since a sample's runs can be
    // spread over any rows.
    // (meta, reads, runs, bam, transcriptome_bam, percent_mapped, prealigned)
    def ch_input = ch_sample_rows
        .collect()
        .flatMap { rows -> readSamplesheet(rows, "${moduleDir}/../../assets/schema_input.json", params.skip_alignment) as List<Tuple> }

    // Samplesheet re-assembled from the validated rows, for the SummarizedExperiment sample metadata
    def ch_samplesheet: Value<Path> = ch_sample_rows
        .collect()
        .flatMap { rows -> [ samplesheetRowsToCsv(rows) ] }
        .collectFile(name: 'samplesheet.csv')
        .collect()
        .map { files -> files.toList().first() as Path }

    // Samples that go through FASTQ preprocessing and samples supplied as pre-aligned BAM files
    ch_fastq_samples = ch_input.filter { _meta, _reads, _runs, _bam, _transcriptome_bam, _percent_mapped, prealigned -> !prealigned }
    ch_bam_samples   = ch_input.filter { _meta, _reads, _runs, _bam, _transcriptome_bam, _percent_mapped, prealigned -> prealigned }

    // One entry per sample, holding one list of FASTQ files per sequencing run: (meta, runs)
    def ch_fastq = ch_fastq_samples.map { meta, _reads, runs, _bam, _transcriptome_bam, _percent_mapped, _prealigned -> tuple(meta, runs) }

    // Index pre-aligned genome BAM files; a sample may supply only a transcriptome BAM
    // (meta, bam, percent_mapped)
    ch_prealigned_genome = ch_bam_samples
        .filter { _meta, _reads, _runs, bam, _transcriptome_bam, _percent_mapped, _prealigned -> bam != null }
        .map { meta, _reads, _runs, bam, _transcriptome_bam, percent_mapped, _prealigned -> tuple(meta, bam, percent_mapped) }
    // (meta, bai)
    def ch_bam_index = SAMTOOLS_INDEX(ch_prealigned_genome.map { meta, bam, _percent_mapped -> tuple(meta, bam, null) })
    // The BAM side is nested so the single null a remainder join pads with keeps the width fixed:
    // (meta, genome(bam, percent_mapped), bai)
    ch_prealigned_indexed = ch_prealigned_genome
        .map { meta, bam, percent_mapped -> tuple(meta, tuple(bam, percent_mapped)) }
        .join(ch_bam_index, by: 0, remainder: true)
    ch_prealigned_indexed.subscribe { meta, genome, bai ->
        if( genome == null || bai == null ) {
            error "Sample '${meta.id}' is missing its pre-aligned BAM index result"
        }
    }
    // (meta, bam, bai, percent_mapped)
    ch_prealigned = ch_prealigned_indexed
        .filter { _meta, genome, bai -> genome != null && bai != null }
        .map { meta, genome, bai -> tuple(meta, genome[0], bai, genome[1]) }

    //
    // Run RNA-seq FASTQ preprocessing subworkflow
    //

    // The bowtie2 rRNA index is built inside the FASTQ subworkflow, not in PREPARE_GENOME_INDICES.
    def make_bowtie2_index = !params.bowtie2_rrna_index && params.remove_ribo_rna && params.ribo_removal_tool == 'bowtie2'

    fastq_preprocessed = FASTQ_QC_TRIM_FILTER_SETSTRANDEDNESS (
        ch_fastq_samples.map { meta, reads, _runs, _bam, _transcriptome_bam, _percent_mapped, _prealigned -> tuple(meta, reads) }, // ch_reads
        ch_fasta,                                   // ch_fasta
        ch_transcript_fasta,                        // ch_transcript_fasta
        ch_gtf,                                     // ch_gtf
        ch_salmon_index,                            // ch_salmon_index
        ch_sortmerna_index,                         // ch_sortmerna_index
        ch_bowtie2_rrna_index,                      // ch_bowtie2_index (for rRNA removal)
        ch_bbsplit_index,                           // ch_bbsplit_index
        ch_ribo_db,                                 // ch_rrna_fastas
        params.skip_bbsplit || !params.fasta,       // skip_bbsplit
        params.skip_fastqc || params.skip_qc,       // skip_fastqc
        params.skip_trimming,                       // skip_trimming
        params.skip_umi_extract,                    // skip_umi_extract
        params.skip_linting,                        // skip_linting
        false,                                      // make_salmon_index (PREPARE_GENOME_INDICES already builds/loads this)
        false,                                      // make_sortmerna_index (PREPARE_GENOME_INDICES already builds/loads this)
        make_bowtie2_index,                         // make_bowtie2_index
        params.trimmer,                             // trimmer
        params.min_trimmed_reads,                   // min_trimmed_reads
        params.save_trimmed,                        // save_trimmed
        false,                                      // fastp_merge
        params.remove_ribo_rna,                     // remove_ribo_rna
        params.ribo_removal_tool,                   // ribo_removal_tool
        params.with_umi,                            // with_umi
        params.umi_discard_read,                    // umi_discard_read
        params.save_merged_fastq,                   // save_merged_fastq
        params.stranded_threshold,                  // stranded_threshold
        params.unstranded_threshold,                // unstranded_threshold
        tool_args                                   // tool_args
    )

    // (meta, reads, reads_cat, reads_trimmed, num_trimmed_reads,
    //  fastqc(raw_html, raw_zip, trim_html, trim_zip, filtered_html, filtered_zip),
    //  trim(html, log, json, unpaired, reads_fail, reads_merged), umi(log, reads),
    //  bbsplit(stats, primary_reads, other_genome_reads), lint(raw, trimmed, bbsplit, ribo),
    //  rrna(sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log))
    def ch_preprocessed = fastq_preprocessed.samples

    // Run-level rRNA references, built by FASTQ_REMOVE_RRNA from the rRNA
    // FASTAs only when no bowtie2 rRNA index was supplied:
    // (sortmerna_index, bowtie2_index, seqkit_prefixed, seqkit_converted)
    def ch_rrna_references = fastq_preprocessed.rrna_references

    // Samples that fail min_trimmed_reads have no filtered reads and go no further: (meta, reads)
    ch_reads_ok = ch_preprocessed
        .filter { _meta, reads, _reads_cat, _reads_trimmed, _num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, _lint, _rrna -> reads != null }
        .map { meta, reads, _reads_cat, _reads_trimmed, _num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, _lint, _rrna -> tuple(meta, reads) }

    // MultiQC files from FASTQ preprocessing. The trimmer log only reaches the merged
    // report when the trimmer is not fastp. A stage that did not run leaves its group null.
    ch_fastq_qc = ch_preprocessed.map { meta, _reads, _reads_cat, _reads_trimmed, _num_trimmed_reads, fastqc, trim, umi, bbsplit, _lint, rrna ->
        def fastqc_head = ((fastqc != null ? fastqc[1] : null) ?: []) as List<Path>
        fastqc_head     = fastqc_head + (((fastqc != null ? fastqc[3] : null) ?: []) as List<Path>)
        def trim_json   = ((trim != null ? trim[2] : null) ?: []) as List<Path>
        def trim_log    = ((trim != null ? trim[1] : null) ?: []) as List<Path>
        def other       = [
            umi     != null ? umi[0]     : null,
            bbsplit != null ? bbsplit[0] : null,
            rrna    != null ? rrna[0]    : null,
            rrna    != null ? rrna[1]    : null,
            rrna    != null ? rrna[2]    : null,
            rrna    != null ? rrna[3]    : null
        ] + (((fastqc != null ? fastqc[5] : null) ?: []) as List<Path>)
        tuple(meta, (fastqc_head + (params.trimmer == 'fastp' ? [] : trim_log) + trim_json + other).findAll { f -> f != null })
    }
    ch_mqc_files = ch_mqc_files.mix(ch_fastq_qc)
    if (params.trimmer == 'fastp') {
        ch_mqc_sample_only = ch_mqc_sample_only.mix(
            ch_preprocessed.map { meta, _reads, _reads_cat, _reads_trimmed, _num_trimmed_reads, _fastqc, trim, _umi, _bbsplit, _lint, _rrna ->
                tuple(meta, (trim != null ? trim[1] : null) ?: [])
            }
        )
    }

    // (meta, num_reads)
    def ch_trim_read_count = ch_preprocessed
        .filter { _meta, _reads, _reads_cat, _reads_trimmed, num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, _lint, _rrna -> num_trimmed_reads != null }
        .map { meta, _reads, _reads_cat, _reads_trimmed, num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, _lint, _rrna -> tuple(meta, num_trimmed_reads) }

    // (meta, pass)
    def ch_trim_status = ch_trim_read_count
        .map { meta, num_reads -> tuple(meta, num_reads > params.min_trimmed_reads.toFloat()) }

    //
    // The three aligners emit different layouts. Each is projected onto one common layout before
    // they are mixed; `logs` starts with the aligner log MultiQC reads (STAR Log.final.out, the
    // HISAT2 summary or the Bowtie2 log), and holds every aligner log that is published:
    // (meta, bam, bai, samtools(stats, flagstat, idxstats), percent_mapped, raw_bams, transcriptome_bam, unmapped, logs)
    //

    //
    // SUBWORKFLOW: Alignment with STAR and gene/transcript quantification with Salmon
    //
    def ch_star: Channel<Tuple> = channel.empty()
    if (!params.skip_alignment && (params.aligner == 'star_salmon' || params.aligner == 'star_rsem')) {
        def star_tool = params.use_sentieon_star ? 'sentieon' : (params.use_parabricks_star ? 'parabricks' : 'star')
        ch_star = ALIGN_STAR (
            ch_reads_ok.map { meta, reads -> tuple(meta, reads, starAlignArgs(params, star_tool, meta)) },
            ch_star_index,
            ch_gtf,
            params.star_ignore_sjdbgtf,
            ch_fasta,
            ch_fai,
            params.use_sentieon_star,
            params.use_parabricks_star,
            params.skip_markduplicates,
            tool_args
        ).map { meta, raw_bams, _bam_sorted, _bam_sorted_aligned, _bam_unsorted, transcriptome_bam, unmapped, _sam, _junction, _spl_junc_tab, _read_per_gene_tab, _wig, _bedgraph, _orig_bai, _qc_metrics, _duplicate_metrics, star, bam, bai, samtools, _aligner, percent_mapped ->
            // star: (log_final, log_out, log_progress, tab)
            tuple(meta, bam, bai, samtools, percent_mapped, raw_bams, transcriptome_bam, unmapped, [star[0], star[1], star[2]] + (star[3] as List<Path>))
        }

        ch_mqc_files = ch_mqc_files.mix(
            ch_star.map { meta, _bam, _bai, _samtools, _percent_mapped, _raw_bams, _transcriptome_bam, _unmapped, logs -> tuple(meta, [logs[0]]) }
        )

        if (!params.with_umi && (params.skip_markduplicates || params.use_parabricks_star)) {
            // When Picard markduplicates runs, its stats are added later; adding these too would
            // duplicate flagstat files in MultiQC. Parabricks skips Picard, so add them here.
            ch_mqc_files = ch_mqc_files.mix(
                ch_star.map { meta, _bam, _bai, samtools, _percent_mapped, _raw_bams, _transcriptome_bam, _unmapped, _logs -> tuple(meta, [samtools[0], samtools[1], samtools[2]]) }
            )
        }
    }

    //
    // SUBWORKFLOW: Alignment with Bowtie2
    //
    def ch_bowtie2: Channel<Tuple> = channel.empty()
    if (!params.skip_alignment && params.aligner == 'bowtie2_salmon') {
        ch_bowtie2 = ALIGN_BOWTIE2 (
            ch_reads_ok.map { meta, reads -> tuple(meta, reads, bowtie2AlignArgs(params, meta)) },
            ch_bowtie2_index,
            ch_fasta,
            ch_fai,
            params.save_unaligned,
            tool_args
        ).map { meta, _sam, raw_bams, _cram, _csi, _crai, unmapped, bowtie2_log, _aligner, percent_mapped, transcriptome_bam, bam, bai, samtools ->
            tuple(meta, bam, bai, samtools, percent_mapped, raw_bams, transcriptome_bam, unmapped, [bowtie2_log])
        }

        ch_mqc_files = ch_mqc_files.mix(
            ch_bowtie2.map { meta, _bam, _bai, _samtools, _percent_mapped, _raw_bams, _transcriptome_bam, _unmapped, logs -> tuple(meta, [logs[0]]) }
        )

        if (!params.with_umi && params.skip_markduplicates) {
            ch_mqc_files = ch_mqc_files.mix(
                ch_bowtie2.map { meta, _bam, _bai, samtools, _percent_mapped, _raw_bams, _transcriptome_bam, _unmapped, _logs -> tuple(meta, [samtools[0], samtools[1], samtools[2]]) }
            )
        }
    }

    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    def ch_hisat2: Channel<Tuple> = channel.empty()
    if (!params.skip_alignment && params.aligner == 'hisat2') {
        ch_hisat2 = FASTQ_ALIGN_HISAT2 (
            ch_reads_ok.map { meta, reads -> tuple(meta, reads, hisat2AlignArgs(params, meta)) },
            ch_hisat2_index,
            ch_splicesites,
            ch_fasta,
            ch_fai,
            params.save_unaligned || (params.contaminant_screening && params.contaminant_screening_input == 'unmapped'),
            tool_args
        ).map { meta, raw_bams, unmapped, hisat2_summary, _aligner, bam, bai, samtools ->
            tuple(meta, bam, bai, samtools, getHisat2PercentMapped(hisat2_summary), raw_bams, null, unmapped, [hisat2_summary])
        }

        ch_mqc_files = ch_mqc_files.mix(
            ch_hisat2.map { meta, _bam, _bai, _samtools, _percent_mapped, _raw_bams, _transcriptome_bam, _unmapped, logs -> tuple(meta, [logs[0]]) }
        )

        if (!params.with_umi && params.skip_markduplicates) {
            ch_mqc_files = ch_mqc_files.mix(
                ch_hisat2.map { meta, _bam, _bai, samtools, _percent_mapped, _raw_bams, _transcriptome_bam, _unmapped, _logs -> tuple(meta, [samtools[0], samtools[1], samtools[2]]) }
            )
        }
    }

    // (meta, bam, bai, samtools(stats, flagstat, idxstats), percent_mapped, raw_bams, transcriptome_bam, unmapped, logs)
    def ch_aligned: Channel<Tuple> = channel.empty().mix(ch_star).mix(ch_bowtie2).mix(ch_hisat2)

    //
    // Genome-aligned BAMs with their index and mapping percentage, from the aligner or the samplesheet:
    // (meta, bam, bai, percent_mapped)
    //
    def ch_genome_bam: Channel<Tuple> = ch_prealigned
        .mix(ch_aligned.map { meta, bam, bai, _samtools, percent_mapped, _raw_bams, _transcriptome_bam, _unmapped, _logs -> tuple(meta, bam, bai, percent_mapped) })

    // (meta, transcriptome_bam)
    def ch_transcriptome_bam: Channel<Tuple> = ch_aligned
        .map { meta, _bam, _bai, _samtools, _percent_mapped, _raw_bams, transcriptome_bam, _unmapped, _logs -> tuple(meta, transcriptome_bam) }
        .mix(ch_bam_samples.map { meta, _reads, _runs, _bam, transcriptome_bam, _percent_mapped, _prealigned -> tuple(meta, transcriptome_bam) })
        .filter { _meta, transcriptome_bam -> transcriptome_bam != null }

    //
    // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
    //
    // (meta, bam, bai, samtools(stats, flagstat, idxstats), genomic_dedup_log,
    //  tsv(edit_distance, per_umi, umi_per_position),
    //  transcriptomic_dedup_log, prepare_for_rsem_log, transcriptome_bam,
    //  transcriptome(transcriptome_bam, dedup_bam, sorted_bam, sorted_bam_index, filtered_bam, samtools, tsv,
    //                coord_sorted_bam, coord_sorted_bam_index, coord_sorted_samtools))
    def ch_umi_dedup: Channel<Tuple> = channel.empty()
    if (!params.skip_alignment && params.with_umi) {
        ch_umi_dedup = BAM_DEDUP_UMI(
            ch_genome_bam.map { meta, bam, bai, _percent_mapped -> tuple(meta, bam, bai) },
            ch_fasta,
            ch_fai,
            params.umi_dedup_tool,
            params.umitools_dedup_stats,
            ch_transcriptome_bam,
            ch_transcript_fasta,
            params.umitools_dedup_primary_only,
            tool_args,
            params.umitools_grouping_method,
            params.umitools_umi_separator
        )

        // The deduplicated bam and bai replace the aligner's; the mapping percentage is kept.
        // The dedup side is nested so the single null a remainder join pads with keeps the width fixed:
        // (meta, percent_mapped, dedup(bam, bai, genomic_dedup_log))
        ch_genome_deduped = ch_genome_bam
            .map { meta, _bam, _bai, percent_mapped -> tuple(meta, percent_mapped) }
            .join(
                ch_umi_dedup.map { meta, bam, bai, _samtools, genomic_dedup_log, _tsv, _transcriptomic_dedup_log, _prepare_for_rsem_log, _transcriptome_bam, _transcriptome ->
                    tuple(meta, tuple(bam, bai, genomic_dedup_log))
                },
                by: 0,
                remainder: true
            )
        ch_genome_deduped.subscribe { meta, _percent_mapped, dedup ->
            if( dedup == null || dedup[2] == null ) {
                error "Sample '${meta.id}' is missing its UMI deduplication result"
            }
        }
        ch_genome_bam = ch_genome_deduped
            .filter { _meta, _percent_mapped, dedup -> dedup != null && dedup[2] != null }
            .map { meta, percent_mapped, dedup -> tuple(meta, dedup[0], dedup[1], percent_mapped) }
        ch_transcriptome_bam = ch_umi_dedup
            .filter { _meta, _bam, _bai, _samtools, _genomic_dedup_log, _tsv, _transcriptomic_dedup_log, _prepare_for_rsem_log, transcriptome_bam, _transcriptome -> transcriptome_bam != null }
            .map { meta, _bam, _bai, _samtools, _genomic_dedup_log, _tsv, _transcriptomic_dedup_log, _prepare_for_rsem_log, transcriptome_bam, _transcriptome -> tuple(meta, transcriptome_bam) }

        // Genome-side files only; MultiQC cannot tell transcriptome stats apart from genome stats
        ch_mqc_files = ch_mqc_files.mix(
            ch_umi_dedup.map { meta, _bam, _bai, samtools, genomic_dedup_log, _tsv, _transcriptomic_dedup_log, _prepare_for_rsem_log, _transcriptome_bam, _transcriptome ->
                tuple(meta, [genomic_dedup_log, samtools[0], samtools[1], samtools[2]])
            }
        )
    }

    //
    // Quantification
    //
    // (meta, reads)
    def ch_transcriptome_reads = ch_transcriptome_bam.map { meta, transcriptome_bam -> tuple(meta, [transcriptome_bam]) }

    def run_deseq2_qc = !params.skip_qc && !params.skip_deseq2_qc && !params.skip_quantification_merge

    // The merged quantification channels hold one row per sample under skip_quantification_merge
    // and a single 'all_samples' row otherwise:
    // (meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, tpm_transcript,
    //  counts_transcript, lengths_transcript, tx2gene_augmented, merged_gene_rds, merged_transcript_rds, tx2gene)
    // DESEQ2_QC: (meta, rdata, pca_vals, plots_pdf, sample_dists, size_factors, log, pca_multiqc, dists_multiqc)
    def ch_quant: Channel<Tuple> = channel.empty()
    def ch_quant_salmon: Channel<Tuple> = channel.empty()
    def ch_quant_merged: Channel<Tuple> = channel.empty()
    def ch_quant_rsem_merge: Channel<Tuple> = channel.empty()
    def ch_deseq2: Channel<Tuple> = channel.empty()
    if (params.aligner == 'star_rsem') {
        rsem = QUANTIFY_RSEM (
            ch_samplesheet,
            ch_transcriptome_reads,
            ch_rsem_index,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            params.use_sentieon_star,
            params.skip_quantification_merge
        )

        ch_quant            = rsem.samples
        ch_quant_merged     = rsem.merged
        ch_quant_rsem_merge = rsem.rsem_merge

        ch_mqc_files = ch_mqc_files.mix(
            rsem.samples.map { meta, _counts_gene, _counts_transcript, stat, _log, _bam_star, _bam_genome, _bam_transcript -> tuple(meta, [stat]) }
        )

        if (run_deseq2_qc) {
            ch_deseq2 = DESEQ2_QC_RSEM (
                ch_quant_merged.map { meta, _tpm_gene, _counts_gene, _lengths_gene, counts_gene_length_scaled, _counts_gene_scaled, _tpm_transcript, _counts_transcript, _lengths_transcript, _tx2gene_augmented, _merged_gene_rds, _merged_transcript_rds, _tx2gene ->
                    tuple(meta, counts_gene_length_scaled, deseq2QcArgs(params), params.aligner, ch_pca_header_multiqc, ch_clustering_header_multiqc)
                }
            )
        }
    } else if (params.aligner in ['star_salmon', 'bowtie2_salmon']) {

        //
        // SUBWORKFLOW: Count reads from BAM alignments using Salmon
        //
        bam_salmon = QUANTIFY_BAM_SALMON (
            ch_samplesheet,
            ch_transcriptome_reads,
            channel.value(null as Path?),
            ch_transcript_fasta,
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            'salmon',
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd,
            params.skip_quantification_merge,
            params.salmon_quant_libtype,
            tool_args
        )

        ch_quant_salmon = bam_salmon.salmon
        ch_quant_merged = bam_salmon.merged

        if (run_deseq2_qc) {
            ch_deseq2 = DESEQ2_QC_BAM_SALMON (
                ch_quant_merged.map { meta, _tpm_gene, _counts_gene, _lengths_gene, counts_gene_length_scaled, _counts_gene_scaled, _tpm_transcript, _counts_transcript, _lengths_transcript, _tx2gene_augmented, _merged_gene_rds, _merged_transcript_rds, _tx2gene ->
                    tuple(meta, counts_gene_length_scaled, deseq2QcArgs(params), params.aligner, ch_pca_header_multiqc, ch_clustering_header_multiqc)
                }
            )
        }
    }

    // (meta, bam, bai, percent_mapped, pass)
    ch_mapped = ch_genome_bam.map { meta, bam, bai, percent_mapped ->
        tuple(meta, bam, bai, percent_mapped, percent_mapped != null ? percent_mapped >= params.min_mapped_reads.toFloat() : null)
    }

    // (meta, percent_mapped)
    def ch_percent_mapped = ch_mapped.map { meta, _bam, _bai, percent_mapped, _pass -> tuple(meta, percent_mapped) }

    // (meta, pass)
    def ch_map_status = ch_mapped
        .filter { _meta, _bam, _bai, _percent_mapped, pass -> pass != null }
        .map { meta, _bam, _bai, _percent_mapped, pass -> tuple(meta, pass as Boolean) }

    // (meta, percent_mapped, pass)
    def ch_percent_mapped_pass = ch_mapped.map { meta, _bam, _bai, percent_mapped, pass -> tuple(meta, percent_mapped, pass) }

    // Samples without a mapping percentage are never filtered
    ch_genome_bam = ch_mapped
        .filter { _meta, _bam, _bai, _percent_mapped, pass -> pass == null || pass }
        .map { meta, bam, bai, percent_mapped, _pass -> tuple(meta, bam, bai, percent_mapped) }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //

    // Some tools (Ex. Parabricks) may have already run marked duplicates during alignment
    def markdups_done = !params.skip_markduplicates && params.use_parabricks_star
    // (meta, bam, bai, cram, metrics, samtools(stats, flagstat, idxstats)), samtools null when RustQC replaces the stats
    def ch_markdup: Channel<Tuple> = channel.empty()
    if (!params.skip_markduplicates && !params.with_umi && !markdups_done) {
        ch_markdup = BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam.map { meta, bam, _bai, _percent_mapped -> tuple(meta, bam) },
            ch_fasta,
            ch_fai,
            !params.use_rustqc,
            tool_args
        )

        // The marked bam and bai replace the aligner's; the mapping percentage is kept.
        // The markdup side is nested so the single null a remainder join pads with keeps the width fixed:
        // (meta, percent_mapped, marked(bam, bai, metrics))
        ch_genome_marked = ch_genome_bam
            .map { meta, _bam, _bai, percent_mapped -> tuple(meta, percent_mapped) }
            .join(
                ch_markdup.map { meta, bam, bai, _cram, metrics, _samtools -> tuple(meta, tuple(bam, bai, metrics)) },
                by: 0,
                remainder: true
            )
        ch_genome_marked.subscribe { meta, _percent_mapped, marked ->
            if( marked == null || marked[0] == null || marked[2] == null ) {
                error "Sample '${meta.id}' is missing its BAM from mark duplicates"
            }
        }
        ch_genome_bam = ch_genome_marked
            .filter { _meta, _percent_mapped, marked -> marked != null && marked[0] != null && marked[2] != null }
            .map { meta, percent_mapped, marked -> tuple(meta, marked[0], marked[1], percent_mapped) }

        ch_mqc_files = ch_mqc_files.mix(
            ch_markdup
                .filter { _meta, _bam, _bai, _cram, _metrics, samtools -> samtools != null }
                .map { meta, _bam, _bai, _cram, metrics, samtools -> tuple(meta, [samtools[0], samtools[1], samtools[2], metrics]) }
        )
        ch_mqc_report_only = ch_mqc_report_only.mix(
            ch_markdup
                .filter { _meta, _bam, _bai, _cram, _metrics, samtools -> samtools == null }
                .map { _meta, _bam, _bai, _cram, metrics, _samtools -> metrics }
        )
    }

    //
    // MODULE: StringTie assembly and quantification
    //
    // (meta, transcript_gtf, abundance, coverage_gtf, ballgown, denovo(transcript_gtf, abundance, coverage_gtf, ballgown)),
    // denovo being the sample's de novo assembly that fed the merged GTF, null unless --stringtie_ignore_gtf
    def ch_stringtie: Channel<Tuple> = channel.empty()
    // (meta, gtf, assemblies, merged_gtf), assemblies holding every (meta, transcript_gtf, abundance, coverage_gtf, ballgown)
    def ch_stringtie_merged: Channel<Tuple> = channel.empty()
    if (!params.skip_stringtie) {
        // (meta, bam, lrbam)
        def ch_stringtie_input = ch_genome_bam.map { meta, bam, _bai, _percent_mapped -> tuple(meta, bam, null) }

        if (params.stringtie_ignore_gtf) {
            stringtie_merge = BAM_STRINGTIE_MERGE(
                ch_stringtie_input,
                channel.value([]),
                ch_gtf
            )
            ch_stringtie_merged = stringtie_merge
            ch_stringtie_gtf = stringtie_merge
                .collect()
                .map { merged -> merged.isEmpty() ? null : merged.toList().collect { _meta, _gtf, _assemblies, merged_gtf -> merged_gtf }.first() }
        } else {
            ch_stringtie_gtf = ch_gtf
        }
        // (meta, transcript_gtf, abundance, coverage_gtf, ballgown)
        def ch_stringtie_samples = STRINGTIE_STRINGTIE(
            ch_stringtie_input
                .map { meta, bam, lrbam -> tuple(meta, bam, lrbam, ['expression-estimation']) }
                .combine(ch_stringtie_gtf)
        )

        if (params.stringtie_ignore_gtf) {
            ch_stringtie = ch_stringtie_samples
                .join(
                    stringtie_merge
                        .flatMap { _meta, _gtf, assemblies, _merged_gtf -> assemblies }
                        .map { meta, transcript_gtf, abundance, coverage_gtf, ballgown -> tuple(meta, tuple(transcript_gtf, abundance, coverage_gtf, ballgown)) },
                    by: 0
                )
        } else {
            ch_stringtie = ch_stringtie_samples.map { meta, transcript_gtf, abundance, coverage_gtf, ballgown ->
                tuple(meta, transcript_gtf, abundance, coverage_gtf, ballgown, null)
            }
        }
    }

    //
    // Pre-compute param-derived values for QC subworkflow
    //
    def biotype = (params.gencode ? "gene_type" : params.featurecounts_group_type) as String
    def rseqc_modules = qc_tools.findAll { tool -> tool.startsWith('rseqc_') }.collect { tool -> tool.replace('rseqc_', '') }

    // BAM_QC_RNASEQ, skipped groups null:
    // (meta, preseq(lc_extrap, log), featurecounts(counts, summary), biotype(tsv, rrna), qualimap,
    //  dupradar(scatter2d, boxplot, hist, dupmatrix, intercept_slope, multiqc, session_info),
    //  rseqc(bamstat, inferexperiment, innerdistance(distance, freq, mean, pdf, rscript),
    //        junctionannotation(bed, interact_bed, xls, pdf, events_pdf, rscript, log), junctionsaturation(pdf, rscript),
    //        readdistribution, readduplication(seq_xls, pos_xls, pdf, rscript), tin(txt, xls)),
    //  mqc_files)
    def ch_bam_qc: Channel<Tuple> = channel.empty()
    // RUSTQC:
    // (meta, samtools(stats, flagstat, idxstats), preseq_lc_extrap,
    //  dupradar(scatter2d, boxplot, hist, dupmatrix, intercept_slope, multiqc), featurecounts(counts, summary),
    //  biotype(tsv, mqc, rrna),
    //  rseqc(bamstat, inferexperiment, readdistribution, tin(txt, xls),
    //        innerdistance(distance, freq, mean, summary, plot, rscript),
    //        junctionannotation(bed, interact_bed, xls, log, plot, rscript), junctionsaturation(summary, plot, rscript),
    //        readduplication(seq_xls, pos_xls, plot, rscript)),
    //  qualimap, all_files)
    def ch_bam_qc_rustqc: Channel<Tuple> = channel.empty()
    // (meta, inferexperiment)
    def ch_inferexperiment: Channel<Tuple> = channel.empty()

    if (!params.skip_qc) {
        if (params.use_rustqc) {
            //
            // MODULE: RustQC - single-pass replacement for multiple QC tools
            //
            ch_bam_qc_rustqc = RUSTQC (
                ch_genome_bam
                    .map { meta, bam, bai, _percent_mapped -> tuple(meta, bam, bai) }
                    .combine(ch_gtf)
                    .map { meta, bam, bai, gtf -> tuple(meta, bam, bai, gtf, rustqcArgs(params, meta)) }
            )

            // Drop non-MultiQC files. Excluding `*.featureCounts.tsv.summary`
            // keeps only the biotype summary, matching the default pipeline's
            // `featureCounts -g gene_biotype` output.
            ch_mqc_files = ch_mqc_files.mix(
                ch_bam_qc_rustqc.map { meta, _samtools, _preseq_lc_extrap, _dupradar, _featurecounts, _biotype, _rseqc, _qualimap, all_files ->
                    tuple(
                        meta,
                        all_files.findAll { f ->
                            !f.name.endsWith('.featureCounts.tsv.summary') &&
                                ((f.name =~ /(?i)\.(txt|tsv|xls|log|stats|flagstat|idxstats|html)$/).find() || f.name.contains('_mqc.'))
                        }.toList()
                    )
                }
            )

            // Extract infer_experiment from rseqc channel
            ch_inferexperiment = ch_bam_qc_rustqc
                .filter { _meta, _samtools, _preseq_lc_extrap, _dupradar, _featurecounts, _biotype, rseqc, _qualimap, _all_files -> rseqc[1] != null }
                .map { meta, _samtools, _preseq_lc_extrap, _dupradar, _featurecounts, _biotype, rseqc, _qualimap, _all_files -> tuple(meta, rseqc[1]) }
        } else {
            //
            // SUBWORKFLOW: Post-alignment QC
            //
            ch_bam_qc = BAM_QC_RNASEQ (
                ch_genome_bam.map { meta, bam, bai, _percent_mapped -> tuple(meta, bam, bai) },
                ch_gtf,
                ch_gene_bed,
                ch_fasta,
                ch_fai,
                ch_biotypes_header_multiqc,
                qc_tools,
                biotype,
                params.featurecounts_feature_type
            )

            ch_mqc_files = ch_mqc_files.mix(
                ch_bam_qc.map { meta, _preseq, _featurecounts, _biotype, _qualimap, _dupradar, _rseqc, mqc_files -> tuple(meta, mqc_files) }
            )
            ch_inferexperiment = ch_bam_qc
                .filter { _meta, _preseq, _featurecounts, _biotype, _qualimap, _dupradar, rseqc, _mqc_files -> rseqc != null && rseqc[1] != null }
                .map { meta, _preseq, _featurecounts, _biotype, _qualimap, _dupradar, rseqc, _mqc_files -> tuple(meta, rseqc[1]) }
        }
    }

    //
    // Build the per-sample strand-classification tuple consumed by the
    // MultiQC Strandedness checks section. When RSeQC / RustQC ran we
    // classify via `classifyStrand`; otherwise we surface Salmon's
    // auto-inference so --skip_rseqc / --skip_qc users still see the
    // available signal. RustQC always emits infer_experiment regardless
    // of rseqc_modules.
    //
    def run_infer_experiment = !params.skip_qc && (params.use_rustqc || rseqc_modules.contains('infer_experiment'))
    // (meta, provided, status, salmon, rseqc)
    def ch_strand_data: Channel<Tuple> = channel.empty()
    // (meta, pass)
    def ch_strand_status: Channel<Tuple> = channel.empty()
    if (run_infer_experiment) {
        ch_strand_data = ch_inferexperiment.map { meta, inferexperiment ->
            classifyStrand(meta, inferexperiment, params.stranded_threshold, params.unstranded_threshold)
        }
        ch_strand_status = ch_strand_data.map { meta, _provided, status, _salmon, _rseqc -> tuple(meta, status == 'pass') }
    }
    else {
        ch_strand_data = ch_reads_ok
            .filter { meta, _reads -> meta.salmon_strand_analysis != null }
            .map { meta, _reads -> tuple(meta, 'auto', '-', meta.salmon_strand_analysis, null) }
    }

    //
    // MODULE: Genome-wide coverage with BEDTools
    // Stranded libraries get per-strand + combined bigWigs; unstranded libraries get only the combined one.
    //
    // (meta, combined(bedgraph, bigwig), forward(bedgraph, bigwig), reverse(bedgraph, bigwig)), forward and reverse null when unstranded
    def ch_bigwig: Channel<Tuple> = channel.empty()
    if (!params.skip_bigwig) {
        // (meta, intervals, scale, sizes, extension, sort)
        ch_genomecov_input = ch_genome_bam.map { meta, bam, _bai, _percent_mapped ->
            tuple(meta, bam, 1 as Float, null, 'bedGraph', true)
        }
        ch_genomecov_input_stranded = ch_genomecov_input.filter { meta, _intervals, _scale, _sizes, _extension, _sort -> meta.strandedness in ['forward', 'reverse'] }

        ch_bedgraph_fw = BEDTOOLS_GENOMECOV_FW (
            ch_genomecov_input_stranded
        )
        ch_bedgraph_rev = BEDTOOLS_GENOMECOV_REV (
            ch_genomecov_input_stranded
        )
        ch_bedgraph_combined = BEDTOOLS_GENOMECOV_COMBINED (
            ch_genomecov_input
        )

        //
        // SUBWORKFLOW: Convert bedGraph to bigWig
        //
        // (meta, bedgraph, bigwig)
        ch_bigwig_fw = BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_FORWARD (
            ch_bedgraph_fw,
            ch_chrom_sizes
        )

        ch_bigwig_rev = BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_REVERSE (
            ch_bedgraph_rev,
            ch_chrom_sizes
        )

        ch_bigwig_combined = BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG_COMBINED (
            ch_bedgraph_combined,
            ch_chrom_sizes
        )

        // Every sample gets a combined track; only stranded ones get forward/reverse. Each track is
        // nested so a missing one is the single null a remainder join pads with. The combined side
        // is never missing, so the second join's left side needs no nesting.
        ch_bigwig = ch_bigwig_combined
            .map { meta, bedgraph, bigwig -> tuple(meta, tuple(bedgraph, bigwig)) }
            .join(ch_bigwig_fw.map { meta, bedgraph, bigwig -> tuple(meta, tuple(bedgraph, bigwig)) }, by: 0, remainder: true)
            .join(ch_bigwig_rev.map { meta, bedgraph, bigwig -> tuple(meta, tuple(bedgraph, bigwig)) }, by: 0, remainder: true)
    }

    // (meta, kraken2(report, classified_reads_fastq, unclassified_reads_fastq, classified_reads_assignment),
    //  bracken(abundance, report), sylph_profile, sylphtax_taxprof); only the screening tool's elements are set
    def ch_contaminants: Channel<Tuple> = channel.empty()
    if (!params.skip_qc) {
        //
        // Contaminant screening (Kraken2/Bracken/Sylph)
        //
        // (meta, reads)
        if (params.contaminant_screening_input == 'trim_only') {
            ch_contaminant_reads = ch_preprocessed
                .filter { _meta, _reads, _reads_cat, reads_trimmed, _num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, _lint, _rrna -> reads_trimmed != null }
                .map { meta, _reads, _reads_cat, reads_trimmed, _num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, _lint, _rrna -> tuple(meta, reads_trimmed) }
        } else if (params.contaminant_screening_input == 'trimmed') {
            ch_contaminant_reads = ch_reads_ok
        } else if (params.contaminant_screening_input == 'raw') {
            ch_contaminant_reads = ch_preprocessed
                .map { meta, _reads, reads_cat, _reads_trimmed, _num_trimmed_reads, _fastqc, _trim, _umi, _bbsplit, _lint, _rrna -> tuple(meta, reads_cat) }
        } else {
            ch_contaminant_reads = ch_star.mix(ch_hisat2)
                .filter { _meta, _bam, _bai, _samtools, _percent_mapped, _raw_bams, _transcriptome_bam, unmapped, _logs -> !unmapped.isEmpty() }
                .map { meta, _bam, _bai, _samtools, _percent_mapped, _raw_bams, _transcriptome_bam, unmapped, _logs -> tuple(meta, unmapped) }
        }

        if (params.contaminant_screening in ['kraken2', 'kraken2_bracken'] ) {
            // (meta, report, classified_reads_fastq, unclassified_reads_fastq, classified_reads_assignment)
            // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
            ch_kraken2 = KRAKEN2 (
                ch_contaminant_reads
                    .map { t -> tuple(t) }
                    .combine(ch_kraken_db)
                    .map { meta, reads, db -> tuple(meta, reads, db, params.save_kraken_assignments, params.save_kraken_unassigned) }
            )
            // (meta, kraken2(report, classified_reads_fastq, unclassified_reads_fastq, classified_reads_assignment))
            ch_kraken2_nested = ch_kraken2.map { meta, report, classified_reads_fastq, unclassified_reads_fastq, classified_reads_assignment ->
                tuple(meta, tuple(report, classified_reads_fastq, unclassified_reads_fastq, classified_reads_assignment))
            }

            if (params.contaminant_screening == 'kraken2') {
                ch_mqc_files = ch_mqc_files.mix(
                    ch_kraken2.map { meta, report, _classified_reads_fastq, _unclassified_reads_fastq, _classified_reads_assignment -> tuple(meta, [report]) }
                )
                ch_contaminants = ch_kraken2_nested.map { meta, kraken2 -> tuple(meta, kraken2, null, null, null) }
            } else if (params.contaminant_screening == 'kraken2_bracken') {
                // (meta, abundance, report)
                ch_bracken = BRACKEN (
                    ch_kraken2
                        .map { meta, report, _classified_reads_fastq, _unclassified_reads_fastq, _classified_reads_assignment -> tuple(meta, report, "-l ${params.bracken_precision}") }
                        .combine(ch_kraken_db)
                )
                ch_mqc_files = ch_mqc_files.mix(ch_bracken.map { meta, _abundance, report -> tuple(meta, [report]) })
                ch_contaminants = ch_kraken2_nested
                    .join(ch_bracken.map { meta, abundance, report -> tuple(meta, tuple(abundance, report)) }, by: 0)
                    .map { meta, kraken2, bracken -> tuple(meta, kraken2, bracken, null, null) }
            }
        } else if (params.contaminant_screening == 'sylph') {
            def sylph_databases = (params.sylph_db ? params.sylph_db.split(',').collect{ path -> file(path.trim()) } : []) as List<Path>
            // (meta, profile_out)
            ch_sylph = SYLPH_PROFILE (
                ch_contaminant_reads.map { meta, reads -> tuple(meta, reads, sylph_databases) }
            )
            ch_sylph_profile = ch_sylph.filter { _meta, profile_out -> !profile_out.isEmpty() }

            def sylph_taxonomies = (params.sylph_taxonomy ? params.sylph_taxonomy.split(',').collect{ path -> file(path.trim()) } : []) as List<Path>
            // (meta, taxprof_output)
            ch_sylphtax = SYLPHTAX_TAXPROF (
                ch_sylph_profile.map { meta, profile_out -> tuple(meta, profile_out, sylph_taxonomies) }
            )
            ch_mqc_files = ch_mqc_files.mix(ch_sylphtax.map { meta, taxprof_output -> tuple(meta, [taxprof_output]) })

            // Every profile is published, but empty ones never reach SYLPHTAX_TAXPROF
            ch_contaminants = ch_sylph
                .join(ch_sylphtax, by: 0, remainder: true)
                .map { meta, profile_out, taxprof_output -> tuple(meta, null, null, profile_out, taxprof_output) }
        }
    }

    //
    // SUBWORKFLOW: Pseudoalignment and quantification with Salmon
    //
    // (meta, quant_dir, json_info, lib_format_counts)
    def ch_quant_pseudo: Channel<Tuple> = channel.empty()
    // (meta, quant_dir, json_info, log)
    def ch_quant_pseudo_kallisto: Channel<Tuple> = channel.empty()
    def ch_quant_merged_pseudo: Channel<Tuple> = channel.empty()
    def ch_deseq2_pseudo: Channel<Tuple> = channel.empty()
    if (!params.skip_pseudo_alignment && params.pseudo_aligner) {

        if (params.pseudo_aligner == 'salmon') {
            ch_pseudo_index = ch_salmon_index
        } else {
            ch_pseudo_index = ch_kallisto_index
        }

        pseudo = QUANTIFY_PSEUDO_ALIGNMENT (
            ch_samplesheet,
            ch_reads_ok,
            ch_pseudo_index,
            channel.value(null as Path?),
            ch_gtf,
            params.gtf_group_features,
            params.gtf_extra_attributes,
            params.pseudo_aligner,
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd,
            params.skip_quantification_merge,
            params.salmon_quant_libtype,
            tool_args + record(se_prefix: "${params.pseudo_aligner}.merged")
        )

        ch_quant_pseudo          = pseudo.salmon
        ch_quant_pseudo_kallisto = pseudo.kallisto
        ch_quant_merged_pseudo   = pseudo.merged

        // MultiQC parses the Salmon quant directory and the Kallisto log
        def ch_pseudo_mqc = ch_quant_pseudo
            .map { meta, quant_dir, _json_info, _lib_format_counts -> tuple(meta, [quant_dir]) }
            .mix(ch_quant_pseudo_kallisto.map { meta, _quant_dir, _json_info, log -> tuple(meta, [log]) })
        ch_mqc_files = ch_mqc_files.mix(ch_pseudo_mqc)

        if (run_deseq2_qc) {
            ch_deseq2_pseudo = DESEQ2_QC_PSEUDO (
                ch_quant_merged_pseudo.map { meta, _tpm_gene, _counts_gene, _lengths_gene, counts_gene_length_scaled, _counts_gene_scaled, _tpm_transcript, _counts_transcript, _lengths_transcript, _tx2gene_augmented, _merged_gene_rds, _merged_transcript_rds, _tx2gene ->
                    tuple(meta, counts_gene_length_scaled, deseq2QcArgs(params), params.pseudo_aligner, ch_pca_header_multiqc, ch_clustering_header_multiqc)
                }
            )
        }
    }

    ch_mqc_report_only = ch_mqc_report_only.mix(
        ch_deseq2.mix(ch_deseq2_pseudo).flatMap { _meta, _rdata, _pca_vals, _plots_pdf, _sample_dists, _size_factors, _log, pca_multiqc, dists_multiqc ->
            [pca_multiqc, dists_multiqc].findAll { f -> f != null }
        }
    )

    //
    // Collate and save software versions from the `versions` topic. Entries are either
    // `path(versions.yml)` (legacy file-emit style) or `(task.process, tool, version)`
    // tuples (inline `eval` style).
    //
    ch_topic_versions = channel.topic('versions').unique()

    ch_versions_string = ch_topic_versions
        .filter { entry -> !(entry instanceof Path) }
        .map { entry -> tuple(entry[0].tokenize(':').last(), "  ${entry[1]}: ${entry[2]}") }
        .groupBy()
        .map { process, tool_versions ->
            "${process}:\n${tool_versions.toSet().toSorted().join('\n')}"
        }

    ch_collated_versions = softwareVersionsToYAML(ch_topic_versions.filter { entry -> entry instanceof Path })
        .mix(ch_versions_string)
        .collectFile(name: 'nf_core_rnaseq_software_mqc_versions.yml', sort: true, newLine: true)
        .map { p -> p as Path }

    //
    // SUBWORKFLOW: MultiQC
    //
    // (meta, report, data, plots)
    def ch_multiqc: Channel<Tuple> = channel.empty()
    def ch_multiqc_report: Channel<Path> = channel.empty()

    if (!params.skip_multiqc) {
        ch_multiqc = MULTIQC_RNASEQ(
            ch_input.map { meta, _reads, _runs, _bam, _transcriptome_bam, _percent_mapped, _prealigned -> meta.id as String },
            ch_mqc_files,
            ch_mqc_sample_only,
            ch_mqc_report_only,
            ch_strand_data,
            ch_trim_read_count,
            ch_percent_mapped_pass,
            aligner_display_name,
            ch_fastq,
            ch_collated_versions,
            file("${moduleDir}/../../assets/multiqc_config.yml", checkIfExists: true),
            params.multiqc_config ? file(params.multiqc_config, checkIfExists: true) : null,
            params.multiqc_logo   ? file(params.multiqc_logo,   checkIfExists: true) : null,
            params.multiqc_methods_description
                ? file(params.multiqc_methods_description)
                : file("${moduleDir}/../../assets/methods_description_template.yml", checkIfExists: true),
            file("${moduleDir}/../../assets/strand_check_summary.yaml",     checkIfExists: true),
            file("${moduleDir}/../../assets/strand_check_composition.yaml", checkIfExists: true),
            sample_status_header_multiqc,
            params.min_trimmed_reads,
            params.skip_quantification_merge,
            tool_args
        )
        ch_multiqc_report = ch_multiqc.map { _meta, report, _data, _plots -> report }
    }

    emit:
    trim_status:           Channel<Tuple> = ch_trim_status           // (meta, pass): meets min_trimmed_reads
    map_status:            Channel<Tuple> = ch_map_status            // (meta, pass): meets min_mapped_reads; samples with a mapping percentage only
    strand_status:         Channel<Tuple> = ch_strand_status         // (meta, pass): strandedness check passed
    multiqc_report:        Channel<Path>  = ch_multiqc_report
    reads:                 Channel<Tuple> = ch_fastq                 // (meta, runs)
    percent_mapped:        Channel<Tuple> = ch_percent_mapped        // (meta, percent_mapped)

    // Per-sample stage results, keyed on meta. Layouts are documented on the channels above.
    preprocessed:          Channel<Tuple> = ch_preprocessed
    aligned:               Channel<Tuple> = ch_aligned               // (meta, bam, bai, samtools(stats, flagstat, idxstats), percent_mapped, raw_bams, transcriptome_bam, unmapped, logs)
    umi_dedup:             Channel<Tuple> = ch_umi_dedup
    markdup:               Channel<Tuple> = ch_markdup
    bam_qc:                Channel<Tuple> = ch_bam_qc
    bam_qc_rustqc:         Channel<Tuple> = ch_bam_qc_rustqc
    quant:                 Channel<Tuple> = ch_quant                 // (meta, counts_gene, counts_transcript, stat, log, bam_star, bam_genome, bam_transcript)
    quant_salmon:          Channel<Tuple> = ch_quant_salmon          // Salmon on the transcriptome BAM: (meta, quant_dir, json_info, lib_format_counts)
    quant_merged:          Channel<Tuple> = ch_quant_merged          // alignment-based quantifier, tximport layout
    quant_rsem_merge:      Channel<Tuple> = ch_quant_rsem_merge      // empty unless --aligner star_rsem: (meta, rsem_merge(counts_gene, tpm_gene, counts_transcript, tpm_transcript, genes_long, isoforms_long))
    quant_pseudo:          Channel<Tuple> = ch_quant_pseudo          // Salmon pseudo-alignment
    quant_pseudo_kallisto: Channel<Tuple> = ch_quant_pseudo_kallisto // Kallisto pseudo-alignment
    quant_merged_pseudo:   Channel<Tuple> = ch_quant_merged_pseudo   // pseudo-aligner, tximport layout
    contaminants:          Channel<Tuple> = ch_contaminants
    stringtie:             Channel<Tuple> = ch_stringtie
    bigwig:                Channel<Tuple> = ch_bigwig

    // Run-level results
    stringtie_merged:      Channel<Tuple> = ch_stringtie_merged
    deseq2:                Channel<Tuple> = ch_deseq2                // alignment-based quantifier
    deseq2_pseudo:         Channel<Tuple> = ch_deseq2_pseudo         // pseudo-aligner
    rrna_references:       Value<Tuple>   = ch_rrna_references       // (sortmerna_index, bowtie2_index, seqkit_prefixed, seqkit_converted)
    multiqc:               Channel<Tuple> = ch_multiqc               // (meta, report, data, plots); per sample under skip_quantification_merge
    pipeline_info:         Channel<Path>  = ch_collated_versions     // collated software versions YAML
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
