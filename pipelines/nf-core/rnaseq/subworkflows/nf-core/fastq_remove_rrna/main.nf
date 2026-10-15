nextflow.enable.types = true

include { BOWTIE2_ALIGN                            } from '../../../modules/nf-core/bowtie2/align'
include { BOWTIE2_ALIGN as BOWTIE2_ALIGN_PE        } from '../../../modules/nf-core/bowtie2/align'
include { BOWTIE2_BUILD                            } from '../../../modules/nf-core/bowtie2/build'
include { RIBODETECTOR                             } from '../../../modules/nf-core/ribodetector'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_BOWTIE2 } from '../../../modules/nf-core/samtools/fastq'
include { SAMTOOLS_VIEW as SAMTOOLS_VIEW_BOWTIE2   } from '../../../modules/nf-core/samtools/view'
include { SEQKIT_REPLACE                           } from '../../../modules/nf-core/seqkit/replace'
include { SEQKIT_REPLACE as SEQKIT_REPLACE_U2T     } from '../../../modules/nf-core/seqkit/replace'
include { SEQKIT_STATS                             } from '../../../modules/nf-core/seqkit/stats'
include { SORTMERNA                                } from '../../../modules/nf-core/sortmerna'
include { SORTMERNA as SORTMERNA_INDEX             } from '../../../modules/nf-core/sortmerna'

//
// Function that parses seqkit stats TSV output to extract the mean read length
// for use with RiboDetector's -l parameter
//
def getReadLengthFromSeqkitStats(stats_file) {
    def lines = stats_file.text.readLines()
    if (lines.size() < 2) {
        return 100 // Default fallback
    }

    def header = lines[0].split('\t')
    def avgLenIdx = header.findIndexOf { col -> col == 'avg_len' }
    if (avgLenIdx < 0) {
        return 100 // Default fallback if column not found
    }

    // Calculate mean avg_len across all files in the stats output
    def avgLens = lines[1..-1].collect { line -> line.split('\t')[avgLenIdx] as float }
    def meanAvgLen = avgLens.sum() / avgLens.size()

    return Math.round(meanAvgLen) as int
}

workflow FASTQ_REMOVE_RRNA {
    take:
    ch_reads: Channel<Tuple> // (meta, reads)
    ch_rrna_fastas: Channel<Path> // one or more fasta files containing rrna sequences
    ch_sortmerna_index: Value<Path> // sortmerna index (optional)
    ch_bowtie2_index: Value<Path> // bowtie2 index (optional)
    ribo_removal_tool: String // 'sortmerna', 'ribodetector', or 'bowtie2'
    make_sortmerna_index: Boolean // Whether to create a sortmerna index before running sortmerna
    make_bowtie2_index: Boolean // Whether to create a bowtie2 index before running bowtie2
    tool_args: Record // use_gpu_ribodetector: whether RiboDetector runs on a GPU

    main:

    // Run-level references built here, emitted separately from the per-sample results
    // (sortmerna_index, bowtie2_index, seqkit_prefixed, seqkit_converted)
    val_refs = channel.value(tuple(null, null, null, null))

    if (ribo_removal_tool == 'sortmerna') {
        ch_sortmerna_fastas = ch_rrna_fastas.collect().map { refs -> refs.toList() }

        ch_sortmerna_idx = ch_sortmerna_index
        if (make_sortmerna_index) {
            // (meta, reads, log, index)
            ch_sortmerna_built = SORTMERNA_INDEX(
                ch_sortmerna_fastas.map { fastas -> tuple([id: 'rrna_refs'], [], fastas, null) }
            )
            ch_sortmerna_idx = ch_sortmerna_built.map { _meta, _reads, _log, index -> index }
            val_refs = ch_sortmerna_built.map { _meta, _reads, _log, index -> tuple(index, null, null, null) }
        }

        // (meta, reads, log, index)
        // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
        ch_sortmerna = SORTMERNA(
            ch_reads
                .map { t -> tuple(t) }
                .combine(ch_sortmerna_fastas)
                .combine(ch_sortmerna_idx)
        )

        // (meta, reads, sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log)
        ch_results = ch_sortmerna
            .map { meta, reads, log, _index -> tuple(meta, reads.isEmpty() ? null : reads, log, null, null, null) }
    }
    else if (ribo_removal_tool == 'ribodetector') {
        // Run seqkit stats to determine average read length
        // (meta, stats)
        ch_seqkit_stats = SEQKIT_STATS(ch_reads)

        // Join stats with reads and calculate read length for RiboDetector
        // (meta, reads, stats)
        ch_reads_with_stats = ch_reads.join(ch_seqkit_stats, by: 0)
        // (meta, reads, log)
        ch_ribodetector = RIBODETECTOR(
            ch_reads_with_stats.map { meta, reads, stats ->
                tuple(meta, reads, getReadLengthFromSeqkitStats(stats), (tool_args.use_gpu_ribodetector ?: false) as Boolean)
            }
        )

        // (meta, reads, sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log)
        ch_results = ch_ribodetector
            .join(ch_seqkit_stats, by: 0)
            .map { meta, reads, log, stats -> tuple(meta, reads, null, log, stats, null) }
    }
    else {
        ch_bowtie2_idx = ch_bowtie2_index
        if (make_bowtie2_index) {
            // Process each rRNA file to add unique prefixes and convert U to T
            // This prevents duplicate sequence IDs in SAM header when combining databases
            // (meta, fastx, out_ext)
            ch_rrna_with_meta = ch_rrna_fastas.map { fasta_file -> tuple([id: fasta_file.baseName], fasta_file, '') }

            // Step 1: Add filename prefixes to sequence headers
            // (meta, fastx)
            ch_seqkit_prefixed = SEQKIT_REPLACE(ch_rrna_with_meta)

            // Step 2: Convert U to T in sequences (RNA to DNA)
            ch_prefixed_fastas = ch_seqkit_prefixed.map { meta, fastx -> tuple([id: "${meta.id}_dna"], fastx, '') }
            // (meta, fastx)
            ch_seqkit_converted = SEQKIT_REPLACE_U2T(ch_prefixed_fastas)

            // Collect processed files (already prefixed and U->T converted)
            // (meta, fasta)
            def ch_combined_fasta = ch_seqkit_converted
                .map { _meta, fastx -> fastx }
                .collectFile(name: 'rrna_combined_dna.fasta', sort: { a, b -> a.name <=> b.name }, newLine: true)
                .collect()
                .map { fastas -> tuple([id: 'rrna_refs'], fastas.toList().first() as Path) }

            // (meta, index)
            ch_bowtie2_built = BOWTIE2_BUILD(ch_combined_fasta)
            ch_bowtie2_idx = ch_bowtie2_built.map { _meta, index -> index }
            val_seqkit_prefixed = ch_seqkit_prefixed
                .map { _meta, fastx -> fastx }
                .collect()
                .map { fastxs -> fastxs.toList().toSorted { f -> f.name } }
            val_seqkit_converted = ch_seqkit_converted
                .map { _meta, fastx -> fastx }
                .collect()
                .map { fastxs -> fastxs.toList().toSorted { f -> f.name } }
            val_refs = ch_bowtie2_idx
                .combine(val_seqkit_prefixed)
                .combine(val_seqkit_converted)
                .map { index, prefixed, converted -> tuple(null, index, prefixed, converted) }
        }

        // For single-end reads: bowtie2's --un-gz works correctly
        // save_unaligned=true outputs unmapped reads directly
        // (meta, sam, raw_bams, cram, csi, crai, unmapped, bowtie2_log)
        ch_bowtie2_se = BOWTIE2_ALIGN(
            ch_reads
                .filter { meta, _reads -> meta.single_end }
                .map { t -> tuple(t) }
                .combine(ch_bowtie2_idx)
                .map { meta, reads, index -> tuple(meta, reads, null, index, null, true, false) }
        )

        // For paired-end reads: bowtie2's --un-conc-gz outputs pairs that didn't
        // align concordantly, which INCLUDES pairs where one mate aligned.
        // We need to filter via samtools to get pairs where BOTH mates are unmapped.
        // (meta, sam, raw_bams, cram, csi, crai, unmapped, bowtie2_log)
        ch_bowtie2_pe = BOWTIE2_ALIGN_PE(
            ch_reads
                .filter { meta, _reads -> !meta.single_end }
                .map { t -> tuple(t) }
                .combine(ch_bowtie2_idx)
                .map { meta, reads, index -> tuple(meta, reads, null, index, null, false, false) }
        )

        // Filter BAM for read pairs where BOTH mates are unmapped (flag 12 = 4 + 8)
        // This removes any pair where at least one mate aligned to rRNA
        // (meta, bam, bai, unselected, unselected_index)
        ch_view = SAMTOOLS_VIEW_BOWTIE2(
            ch_bowtie2_pe
                .filter { _meta, _sam, raw_bams, _cram, _csi, _crai, _unmapped, _bowtie2_log -> !raw_bams.isEmpty() }
                .map { meta, _sam, raw_bams, _cram, _csi, _crai, _unmapped, _bowtie2_log ->
                    tuple(meta, raw_bams[0], null, null, null, null, null, '')
                }
        )
        // Note: samtools/view versions collected via topic
        ch_view_bam = ch_view.filter { _meta, bam, _bai, _unselected, _unselected_index -> bam != null }

        // Convert filtered BAM back to paired FASTQ
        // (meta, reads, interleaved, singleton, other)
        ch_fastq_pe = SAMTOOLS_FASTQ_BOWTIE2(
            ch_view_bam.map { meta, bam, _bai, _unselected, _unselected_index -> tuple(meta, bam, false) }
        )

        // Combine single-end and paired-end results
        // (meta, reads)
        ch_filtered_reads = ch_bowtie2_se
            .filter { _meta, _sam, _raw_bams, _cram, _csi, _crai, unmapped, _bowtie2_log -> !unmapped.isEmpty() }
            .map { meta, _sam, _raw_bams, _cram, _csi, _crai, unmapped, _bowtie2_log -> tuple(meta, unmapped) }
            .mix(
                ch_fastq_pe
                    .filter { _meta, reads, _interleaved, _singleton, _other -> !reads.isEmpty() }
                    .map { meta, reads, _interleaved, _singleton, _other -> tuple(meta, reads) }
            )

        // (meta, bowtie2_log)
        ch_bowtie2_logs = ch_bowtie2_se
            .mix(ch_bowtie2_pe)
            .map { meta, _sam, _raw_bams, _cram, _csi, _crai, _unmapped, bowtie2_log -> tuple(meta, bowtie2_log) }

        // (meta, reads, sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log)
        ch_results = ch_bowtie2_logs
            .join(ch_filtered_reads, by: 0, remainder: true)
            .map { meta, bowtie2_log, reads -> tuple(meta, reads, null, null, null, bowtie2_log) }
    }

    emit:
    samples: Channel<Tuple> = ch_results // (meta, reads, sortmerna_log, ribodetector_log, seqkit_stats, bowtie2_log)
    references: Value<Tuple> = val_refs // (sortmerna_index, bowtie2_index, seqkit_prefixed, seqkit_converted)
}
