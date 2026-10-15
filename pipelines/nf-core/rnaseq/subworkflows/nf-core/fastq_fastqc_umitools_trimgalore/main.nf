nextflow.enable.types = true

//
// Read QC, UMI extraction and trimming
//

include { FASTQC           } from '../../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT } from '../../../modules/nf-core/umitools/extract/main'
include { TRIMGALORE       } from '../../../modules/nf-core/trimgalore/main'

//
// Function that parses TrimGalore log output file to get total number of reads after trimming
//
def getTrimGaloreReadsAfterFiltering(log_file) {
    def total_reads = 0
    def filtered_reads = 0
    log_file.eachLine { line ->
        def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
        def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s*([\d\.]+)/
        if (total_reads_matcher) {
            total_reads = total_reads_matcher[0][1].toFloat()
        }
        if (filtered_reads_matcher) {
            filtered_reads = filtered_reads_matcher[0][1].toFloat()
        }
    }
    return total_reads - filtered_reads
}

workflow FASTQ_FASTQC_UMITOOLS_TRIMGALORE {
    take:
    ch_reads: Channel // records with meta, reads; their other fields are carried through to the output
    skip_fastqc: Boolean // true/false
    with_umi: Boolean // true/false
    skip_umi_extract: Boolean // true/false
    skip_trimming: Boolean // true/false
    umi_discard_read: Integer // 0, 1 or 2
    min_trimmed_reads: Integer // > 0
    tool_args: Record // umi_extract, trimgalore: UMI-tools extract and TrimGalore options

    main:
    // Each stage that runs joins its outputs onto the incoming record, overwriting the null
    // placeholders of the fields it owns. Joins on meta happen before the UMI read discard
    // changes it, or after on both sides.
    ch_results = ch_reads.map { r ->
        r + record(
            fastqc_raw_html:   null,
            fastqc_raw_zip:    null,
            umi:               null,
            trim:              null,
            num_trimmed_reads: null
        )
    }

    if (!skip_fastqc) {
        ch_fastqc = FASTQC(ch_reads.map { r -> record(meta: r.meta, reads: r.reads) })
        ch_results = ch_results.join(
            ch_fastqc.map { r -> record(meta: r.meta, fastqc_raw_html: r.html, fastqc_raw_zip: r.zip) },
            by: 'meta'
        )
    }

    if (with_umi && !skip_umi_extract) {
        ch_umi = UMITOOLS_EXTRACT(ch_reads.map { r -> record(meta: r.meta, reads: r.reads, args: tool_args.umi_extract) })
        ch_results = ch_results
            .join(ch_umi.map { r -> record(meta: r.meta, umi: record(log: r.log, reads: r.reads)) }, by: 'meta')
            .map { r ->
                // Discard R1 / R2 if required
                def discard = umi_discard_read in [1, 2] && !r.meta.single_end
                r + record(
                    meta:  discard ? r.meta + record(single_end: true) : r.meta,
                    reads: discard ? [r.umi.reads[umi_discard_read % 2]] : r.umi.reads
                )
            }
    }

    if (!skip_trimming) {
        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        ch_trim = TRIMGALORE(ch_results.map { r -> record(meta: r.meta, reads: r.reads, args: tool_args.trimgalore) }).map { r ->
            def num_reads = r.log.isEmpty()
                ? (min_trimmed_reads as Float) + 1
                : getTrimGaloreReadsAfterFiltering(r.log[-1]) as Float
            record(
                meta:              r.meta,
                reads:             num_reads >= (min_trimmed_reads as Float) ? r.reads : null,
                // TrimGalore reports and unpaired reads are all optional outputs
                trim: record(
                    html:     r.html.isEmpty() ? null : r.html,
                    zip:      r.zip.isEmpty() ? null : r.zip,
                    log:      r.log.isEmpty() ? null : r.log,
                    json:     r.json.isEmpty() ? null : r.json,
                    unpaired: r.unpaired.isEmpty() ? null : r.unpaired
                ),
                num_trimmed_reads: num_reads
            )
        }
        ch_results = ch_results.join(ch_trim, by: 'meta')
    }

    emit:
    ch_results
}
