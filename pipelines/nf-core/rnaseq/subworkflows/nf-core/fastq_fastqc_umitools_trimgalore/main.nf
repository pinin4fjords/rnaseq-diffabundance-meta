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
    ch_reads: Channel<Tuple> // (meta, reads)
    skip_fastqc: Boolean // true/false
    with_umi: Boolean // true/false
    skip_umi_extract: Boolean // true/false
    skip_trimming: Boolean // true/false
    umi_discard_read: Integer // 0, 1 or 2
    min_trimmed_reads: Integer // > 0
    tool_args: Record // umi_extract, trimgalore: UMI-tools extract and TrimGalore options

    main:
    // Each stage that does not run fills its elements with null. Joins on meta happen before
    // the UMI read discard changes it, or after on both sides.

    // (meta, reads, fastqc_raw_html, fastqc_raw_zip)
    ch_fastqc_raw = ch_reads.map { meta, reads -> tuple(meta, reads, null, null) }
    if (!skip_fastqc) {
        ch_fastqc_raw = ch_reads.join(FASTQC(ch_reads), by: 0)
    }

    // (meta, reads, fastqc_raw_html, fastqc_raw_zip, umi(log, reads))
    ch_umi_extracted = ch_fastqc_raw.map { meta, reads, fastqc_raw_html, fastqc_raw_zip ->
        tuple(meta, reads, fastqc_raw_html, fastqc_raw_zip, null)
    }
    if (with_umi && !skip_umi_extract) {
        // (meta, reads, log)
        ch_umi = UMITOOLS_EXTRACT(ch_reads.map { meta, reads -> tuple(meta, reads, tool_args.umi_extract) })
        ch_umi_extracted = ch_fastqc_raw
            .join(ch_umi, by: 0)
            .map { meta, _reads, fastqc_raw_html, fastqc_raw_zip, umi_reads, umi_log ->
                // Discard R1 / R2 if required
                def discard = umi_discard_read in [1, 2] && !meta.single_end
                tuple(
                    discard ? meta + [single_end: true] : meta,
                    discard ? [umi_reads[umi_discard_read % 2]] : umi_reads,
                    fastqc_raw_html,
                    fastqc_raw_zip,
                    tuple(umi_log, umi_reads)
                )
            }
    }

    // (meta, reads, fastqc_raw_html, fastqc_raw_zip, umi(log, reads), trim(html, zip, log, json, unpaired), num_trimmed_reads)
    ch_results = ch_umi_extracted.map { meta, reads, fastqc_raw_html, fastqc_raw_zip, umi ->
        tuple(meta, reads, fastqc_raw_html, fastqc_raw_zip, umi, null, null)
    }
    if (!skip_trimming) {
        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        // (meta, reads, trim(html, zip, log, json, unpaired), num_trimmed_reads)
        ch_trim = TRIMGALORE(
            ch_umi_extracted.map { meta, reads, _fastqc_raw_html, _fastqc_raw_zip, _umi -> tuple(meta, reads, tool_args.trimgalore) }
        ).map { meta, reads, log, json, unpaired, html, zip ->
            def num_reads = log.isEmpty()
                ? (min_trimmed_reads as Float) + 1
                : getTrimGaloreReadsAfterFiltering(log[-1]) as Float
            tuple(
                meta,
                num_reads >= (min_trimmed_reads as Float) ? reads : null,
                // TrimGalore reports and unpaired reads are all optional outputs
                tuple(
                    html.isEmpty() ? null : html,
                    zip.isEmpty() ? null : zip,
                    log.isEmpty() ? null : log,
                    json.isEmpty() ? null : json,
                    unpaired.isEmpty() ? null : unpaired
                ),
                num_reads
            )
        }
        ch_results = ch_umi_extracted
            .join(ch_trim, by: 0)
            .map { meta, _reads, fastqc_raw_html, fastqc_raw_zip, umi, trimmed_reads, trim, num_trimmed_reads ->
                tuple(meta, trimmed_reads, fastqc_raw_html, fastqc_raw_zip, umi, trim, num_trimmed_reads)
            }
    }

    emit:
    // (meta, reads, fastqc_raw_html, fastqc_raw_zip, umi(log, reads), trim(html, zip, log, json, unpaired), num_trimmed_reads)
    ch_results
}
