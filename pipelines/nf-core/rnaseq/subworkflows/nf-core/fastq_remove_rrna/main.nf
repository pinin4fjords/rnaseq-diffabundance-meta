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
    ch_reads: Channel<Record> // meta, reads
    ch_rrna_fastas: Channel<Path> // one or more fasta files containing rrna sequences
    ch_sortmerna_index: Value<Path> // sortmerna index (optional)
    ch_bowtie2_index: Value<Path> // bowtie2 index (optional)
    ribo_removal_tool: String // 'sortmerna', 'ribodetector', or 'bowtie2'
    make_sortmerna_index: Boolean // Whether to create a sortmerna index before running sortmerna
    make_bowtie2_index: Boolean // Whether to create a bowtie2 index before running bowtie2
    tool_args: Record // use_gpu_ribodetector: whether RiboDetector runs on a GPU

    main:

    // Run-level references built here, emitted separately from the per-sample results
    val_refs = channel.value(
        record(sortmerna_index: null, bowtie2_index: null, seqkit_prefixed: null, seqkit_converted: null)
    )

    if (ribo_removal_tool == 'sortmerna') {
        ch_sortmerna_fastas = ch_rrna_fastas.collect().map { refs -> refs.toList() }

        ch_sortmerna_idx = ch_sortmerna_index
        if (make_sortmerna_index) {
            ch_sortmerna_built = SORTMERNA_INDEX(
                channel.value(record(meta: record(id: 'rrna_refs'), reads: [], index: null))
                    .combine(fastas: ch_sortmerna_fastas)
            )
            ch_sortmerna_idx = ch_sortmerna_built.map { r -> r.index }
            val_refs = ch_sortmerna_built.map { r ->
                record(sortmerna_index: r.index, bowtie2_index: null, seqkit_prefixed: null, seqkit_converted: null)
            }
        }

        ch_sortmerna = SORTMERNA(
            ch_reads
                .map { r -> record(meta: r.meta, reads: r.reads) }
                .combine(fastas: ch_sortmerna_fastas, index: ch_sortmerna_idx)
        )

        ch_results = ch_sortmerna
            .map { r ->
                record(
                    meta:             r.meta,
                    reads:            r.reads.isEmpty() ? null : r.reads,
                    sortmerna_log:    r.log,
                    ribodetector_log: null,
                    seqkit_stats:     null,
                    bowtie2_log:      null
                )
            }
    }
    else if (ribo_removal_tool == 'ribodetector') {
        // Run seqkit stats to determine average read length
        ch_seqkit_stats = SEQKIT_STATS(ch_reads.map { r -> record(meta: r.meta, reads: r.reads) })

        // Join stats with reads and calculate read length for RiboDetector
        ch_reads_with_stats = ch_reads.map { r -> record(meta: r.meta, reads: r.reads) }.join(ch_seqkit_stats, by: 'meta')
        ch_ribodetector = RIBODETECTOR(
            ch_reads_with_stats.map { r -> record(meta: r.meta, reads: r.reads, length: getReadLengthFromSeqkitStats(r.stats), gpu: (tool_args.use_gpu_ribodetector ?: false) as Boolean) }
        )

        ch_results = ch_ribodetector
            .join(ch_seqkit_stats.map { r -> record(meta: r.meta, seqkit_stats: r.stats) }, by: 'meta')
            .map { r ->
                record(
                    meta:             r.meta,
                    reads:            r.reads,
                    sortmerna_log:    null,
                    ribodetector_log: r.log,
                    seqkit_stats:     r.seqkit_stats,
                    bowtie2_log:      null
                )
            }
    }
    else {
        ch_bowtie2_idx = ch_bowtie2_index
        if (make_bowtie2_index) {
            // Process each rRNA file to add unique prefixes and convert U to T
            // This prevents duplicate sequence IDs in SAM header when combining databases
            ch_rrna_with_meta = ch_rrna_fastas.map { fasta_file ->
                record(meta: record(id: fasta_file.baseName), fastx: fasta_file, out_ext: '')
            }

            // Step 1: Add filename prefixes to sequence headers
            ch_seqkit_prefixed = SEQKIT_REPLACE(ch_rrna_with_meta)

            // Step 2: Convert U to T in sequences (RNA to DNA)
            ch_prefixed_fastas = ch_seqkit_prefixed.map { r ->
                record(meta: record(id: "${r.meta.id}_dna"), fastx: r.fastx, out_ext: '')
            }
            ch_seqkit_converted = SEQKIT_REPLACE_U2T(ch_prefixed_fastas)

            // Collect processed files (already prefixed and U->T converted)
            def ch_combined_fasta = ch_seqkit_converted
                .map { r -> r.fastx }
                .collectFile(name: 'rrna_combined_dna.fasta', sort: { a, b -> a.name <=> b.name }, newLine: true)
                .collect()
                .map { fastas -> record(meta: record(id: 'rrna_refs'), fasta: fastas.toList().first() as Path) }

            ch_bowtie2_built = BOWTIE2_BUILD(ch_combined_fasta)
            ch_bowtie2_idx = ch_bowtie2_built.map { built -> built.index }
            val_seqkit_prefixed = ch_seqkit_prefixed
                .collect()
                .map { built -> built.collect { r -> r.fastx }.toSorted { f -> f.name } }
            val_seqkit_converted = ch_seqkit_converted
                .collect()
                .map { built -> built.collect { r -> r.fastx }.toSorted { f -> f.name } }
            val_refs = ch_bowtie2_built
                .combine(val_seqkit_prefixed)
                .combine(val_seqkit_converted)
                .map { built, prefixed, converted ->
                    record(sortmerna_index: null, bowtie2_index: built.index, seqkit_prefixed: prefixed, seqkit_converted: converted)
                }
        }

        // For single-end reads: bowtie2's --un-gz works correctly
        // save_unaligned=true outputs unmapped reads directly
        ch_bowtie2_se = BOWTIE2_ALIGN(
            ch_reads
                .filter { r -> r.meta.single_end }
                .map { r -> record(meta: r.meta, reads: r.reads, fasta: null, save_unaligned: true, sort_bam: false) }
                .combine(index: ch_bowtie2_idx)
        )

        // For paired-end reads: bowtie2's --un-conc-gz outputs pairs that didn't
        // align concordantly, which INCLUDES pairs where one mate aligned.
        // We need to filter via samtools to get pairs where BOTH mates are unmapped.
        ch_bowtie2_pe = BOWTIE2_ALIGN_PE(
            ch_reads
                .filter { r -> !r.meta.single_end }
                .map { r -> record(meta: r.meta, reads: r.reads, fasta: null, save_unaligned: false, sort_bam: false) }
                .combine(index: ch_bowtie2_idx)
        )

        // Filter BAM for read pairs where BOTH mates are unmapped (flag 12 = 4 + 8)
        // This removes any pair where at least one mate aligned to rRNA
        ch_view = SAMTOOLS_VIEW_BOWTIE2(
            ch_bowtie2_pe
                .filter { r -> !r.raw_bams.isEmpty() }
                .map { r -> record(meta: r.meta, bam: r.raw_bams[0], bai: null, fasta: null, fai: null, qname: null, bed: null, index_format: '') }
        )
        // Note: samtools/view versions collected via topic
        ch_view_bam = ch_view.filter { r -> r.bam != null }

        // Convert filtered BAM back to paired FASTQ
        ch_fastq_pe = SAMTOOLS_FASTQ_BOWTIE2(ch_view_bam.map { r -> record(meta: r.meta, bam: r.bam, interleave: false) })

        // Combine single-end and paired-end results
        ch_filtered_reads = ch_bowtie2_se
            .filter { r -> !r.unmapped.isEmpty() }
            .map { r -> record(meta: r.meta, reads: r.unmapped) }
            .mix(ch_fastq_pe.filter { r -> !r.reads.isEmpty() })

        ch_bowtie2_logs = ch_bowtie2_se
            .mix(ch_bowtie2_pe)
            .map { r ->
                record(
                    meta:             r.meta,
                    reads:            null,
                    sortmerna_log:    null,
                    ribodetector_log: null,
                    seqkit_stats:     null,
                    bowtie2_log:      r.bowtie2.log
                )
            }

        ch_results = ch_bowtie2_logs
            .join(ch_filtered_reads, by: 'meta', remainder: true)
    }

    emit:
    samples: Channel<Record> = ch_results
    references: Value<Record> = val_refs
}
