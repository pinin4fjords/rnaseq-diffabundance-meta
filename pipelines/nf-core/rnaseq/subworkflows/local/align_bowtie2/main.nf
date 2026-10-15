nextflow.enable.types = true

//
// Alignment with Bowtie2
//
// Bowtie2 is a splice-unaware aligner, suitable for prokaryotic RNA-seq data.
// It aligns reads to the transcriptome (not genome) for Salmon quantification.
// The Bowtie2 index is built from transcript FASTA, enabling alignment-based
// Salmon quantification similar to the STAR transcriptome BAM workflow.
//

include { BOWTIE2_ALIGN           } from '../../../modules/nf-core/bowtie2/align'
include { BAM_SORT_STATS_SAMTOOLS } from '../../nf-core/bam_sort_stats_samtools'

//
// Function that parses and returns the alignment rate from the Bowtie2 log output
//
def getBowtie2PercentMapped(align_log) {
    def percent_aligned = 0
    def pattern = /(\d+\.\d+)% overall alignment rate/
    align_log.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_aligned = matcher[0][1].toFloat()
        }
    }
    return percent_aligned
}

workflow ALIGN_BOWTIE2 {
    take:
    ch_samples: Channel<Tuple> // (meta, reads, args)
    index: Value<Path> // /path/to/bowtie2/index/
    fasta: Value<Path?>
    fai: Value<Path?>
    save_unaligned: Boolean
    tool_args: Record // the samtools index options, forwarded

    main:

    //
    // Map reads with Bowtie2. No fasta is needed for BAM output, and the BAM is sorted with samtools for consistency.
    //
    // (meta, sam, raw_bams, cram, csi, crai, unmapped, bowtie2_log)
    // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
    def ch_bowtie2 = BOWTIE2_ALIGN(
        ch_samples
            .map { t -> tuple(t) }
            .combine(index)
            .map { meta, reads, args, index_dir -> tuple(meta, reads, args, index_dir, null, save_unaligned, false) }
    ).filter { _meta, _sam, raw_bams, _cram, _csi, _crai, _unmapped, _bowtie2_log -> !raw_bams.isEmpty() }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    // (meta, bam, bai, samtools(stats, flagstat, idxstats))
    def ch_sorted = BAM_SORT_STATS_SAMTOOLS(
        ch_bowtie2.map { meta, _sam, raw_bams, _cram, _csi, _crai, _unmapped, _bowtie2_log -> tuple(meta, raw_bams) },
        fasta,
        fai,
        tool_args
    )

    // The BAM is aligned to the transcriptome, and its unsorted form is what Salmon quantifies:
    // a coordinate-sorted BAM breaks paired-end quantification.
    // A remainder join pads a missing side with a single null, so each side's payload is nested
    // to keep the tuple width fixed: (meta, bowtie2_out, sorted)
    ch_bowtie2_sorted = ch_bowtie2
        .map { meta, sam, raw_bams, cram, csi, crai, unmapped, bowtie2_log ->
            tuple(meta, tuple(sam, raw_bams, cram, csi, crai, unmapped, bowtie2_log, 'bowtie2', getBowtie2PercentMapped(bowtie2_log), raw_bams[0]))
        }
        .join(ch_sorted.map { meta, bam, bai, samtools -> tuple(meta, tuple(bam, bai, samtools)) }, by: 0, remainder: true)
    ch_bowtie2_sorted.subscribe { meta, bowtie2_out, sorted ->
        if( bowtie2_out == null || sorted == null ) {
            error "Sample '${meta.id}' is missing its Bowtie2 sorted BAM result"
        }
    }

    // (meta, sam, raw_bams, cram, csi, crai, unmapped, bowtie2_log, aligner, percent_mapped, transcriptome_bam, bam, bai, samtools(stats, flagstat, idxstats))
    ch_results = ch_bowtie2_sorted
        .filter { _meta, bowtie2_out, sorted -> bowtie2_out != null && sorted != null }
        .map { meta, bowtie2_out, sorted ->
            def (sam, raw_bams, cram, csi, crai, unmapped, bowtie2_log, aligner, percent_mapped, transcriptome_bam) = bowtie2_out
            def (bam, bai, samtools) = sorted
            tuple(meta, sam, raw_bams, cram, csi, crai, unmapped, bowtie2_log, aligner, percent_mapped, transcriptome_bam, bam, bai, samtools)
        }

    emit:
    ch_results
}
