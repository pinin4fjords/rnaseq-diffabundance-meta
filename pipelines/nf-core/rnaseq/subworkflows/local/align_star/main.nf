nextflow.enable.types = true

//
// Alignment with STAR
//
include { SENTIEON_STARALIGN as SENTIEON_STAR_ALIGN } from '../../../modules/nf-core/sentieon/staralign/main'
include { PARABRICKS_RNAFQ2BAM as PARABRICKS_RNA_FQ2BAM } from '../../../modules/nf-core/parabricks/rnafq2bam/main'
include { STAR_ALIGN                                } from '../../../modules/nf-core/star/align'
include { BAM_SORT_STATS_SAMTOOLS                   } from '../../nf-core/bam_sort_stats_samtools'


//
// Function that parses and returns the alignment rate from the STAR log output
//
def getStarPercentMapped(_params, align_log) {
    def percent_aligned = 0
    def pattern = /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/
    align_log.eachLine { line ->
        def matcher = line =~ pattern
        if (matcher) {
            percent_aligned = matcher[0][1].toFloat()
        }
    }

    return percent_aligned
}

workflow ALIGN_STAR {
    take:
    ch_samples: Channel<Tuple> // (meta, reads, args)
    index: Value<Path>
    gtf: Value<Path?>
    star_ignore_sjdbgtf: Boolean // when using pre-built STAR indices do not re-extract and use splice junctions from the GTF file
    fasta: Value<Path?>
    fai: Value<Path?>
    use_sentieon_star: Boolean // whether star alignment is accelerated with Sentieon
    use_parabricks_star: Boolean // whether star alignment (and mark duplicates) is accelerated with Parabricks
    skip_markduplicates: Boolean // whether to skip marking duplicates
    tool_args: Record // the samtools index options, forwarded

    main:

    //
    // Map reads with STAR. The three modules share one output layout:
    // (meta, raw_bams, bam_sorted, bam_sorted_aligned, bam_unsorted, transcriptome_bam, unmapped, sam, junction, spl_junc_tab, read_per_gene_tab, wig, bedgraph, orig_bai, qc_metrics, duplicate_metrics, star(log_final, log_out, log_progress, tab))
    //
    def ch_star_out
    if (use_sentieon_star) {
        // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
        ch_star_out = SENTIEON_STAR_ALIGN(ch_samples.map { t -> tuple(t) }.combine(index).combine(gtf).combine(channel.value(star_ignore_sjdbgtf)))
    } else if (use_parabricks_star) {
        ch_star_out = PARABRICKS_RNA_FQ2BAM(
            ch_samples
                .map { t -> tuple(t) }
                .combine(fasta)
                .combine(index)
                .map { meta, reads, args, fasta_file, index_dir -> tuple(meta, reads, fasta_file, index_dir, true, !skip_markduplicates, args) }
        )
    } else {
        ch_star_out = STAR_ALIGN(ch_samples.map { t -> tuple(t) }.combine(index).combine(gtf).combine(channel.value(star_ignore_sjdbgtf)))
    }

    // A run that produced no BAM drops out of the downstream channels.
    ch_star = ch_star_out.filter { _meta, raw_bams, _bam_sorted, _bam_sorted_aligned, _bam_unsorted, _transcriptome_bam, _unmapped, _sam, _junction, _spl_junc_tab, _read_per_gene_tab, _wig, _bedgraph, _orig_bai, _qc_metrics, _duplicate_metrics, _star ->
        !raw_bams.isEmpty()
    }

    //
    // Sort, index BAM file and run samtools stats, flagstat and idxstats
    //
    // (meta, bam, bai, samtools(stats, flagstat, idxstats))
    def ch_sorted = BAM_SORT_STATS_SAMTOOLS(
        ch_star.map { meta, raw_bams, _bam_sorted, _bam_sorted_aligned, _bam_unsorted, _transcriptome_bam, _unmapped, _sam, _junction, _spl_junc_tab, _read_per_gene_tab, _wig, _bedgraph, _orig_bai, _qc_metrics, _duplicate_metrics, _star ->
            tuple(meta, raw_bams)
        },
        fasta,
        fai,
        tool_args
    )

    // A remainder join pads a missing side with a single null, so each side's payload is nested
    // to keep the tuple width fixed: (meta, star_out, sorted)
    ch_star_sorted = ch_star
        .map { meta, raw_bams, bam_sorted, bam_sorted_aligned, bam_unsorted, transcriptome_bam, unmapped, sam, junction, spl_junc_tab, read_per_gene_tab, wig, bedgraph, orig_bai, qc_metrics, duplicate_metrics, star ->
            tuple(meta, tuple(raw_bams, bam_sorted, bam_sorted_aligned, bam_unsorted, transcriptome_bam, unmapped, sam, junction, spl_junc_tab, read_per_gene_tab, wig, bedgraph, orig_bai, qc_metrics, duplicate_metrics, star))
        }
        .join(ch_sorted.map { meta, bam, bai, samtools -> tuple(meta, tuple(bam, bai, samtools)) }, by: 0, remainder: true)
    ch_star_sorted.subscribe { meta, star_out, sorted ->
        if( star_out == null || sorted == null ) {
            error "Sample '${meta.id}' is missing its STAR sorted BAM result"
        }
    }

    // (meta, raw_bams, bam_sorted, bam_sorted_aligned, bam_unsorted, transcriptome_bam, unmapped, sam, junction, spl_junc_tab, read_per_gene_tab, wig, bedgraph, orig_bai, qc_metrics, duplicate_metrics, star(log_final, log_out, log_progress, tab), bam, bai, samtools(stats, flagstat, idxstats), aligner, percent_mapped)
    ch_results = ch_star_sorted
        .filter { _meta, star_out, sorted -> star_out != null && sorted != null }
        .map { meta, star_out, sorted ->
            def (raw_bams, bam_sorted, bam_sorted_aligned, bam_unsorted, transcriptome_bam, unmapped, sam, junction, spl_junc_tab, read_per_gene_tab, wig, bedgraph, orig_bai, qc_metrics, duplicate_metrics, star) = star_out
            def (bam, bai, samtools) = sorted
            def log_final = star[0]
            tuple(
                meta, raw_bams, bam_sorted, bam_sorted_aligned, bam_unsorted, transcriptome_bam, unmapped, sam, junction, spl_junc_tab, read_per_gene_tab, wig, bedgraph, orig_bai, qc_metrics, duplicate_metrics, star,
                bam, bai, samtools,
                'star', getStarPercentMapped(params, log_final)
            )
        }

    emit:
    ch_results
}
