nextflow.enable.types = true

//
// Read QC, UMI extraction and trimming
//
include { FASTQC as FASTQC_RAW  } from '../../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIM } from '../../../modules/nf-core/fastqc/main'
include { UMITOOLS_EXTRACT      } from '../../../modules/nf-core/umitools/extract/main'
include { FASTP                 } from '../../../modules/nf-core/fastp/main'

//
// Function that parses fastp json output file to get total number of reads after trimming
//

def getFastpReadsAfterFiltering(json_file, min_num_reads) {

    if (workflow.stubRun) {
        return min_num_reads
    }

    def json = (new groovy.json.JsonSlurper().parseText(json_file.text) as Map)['summary'] as Map
    return (json['after_filtering'] as Map)['total_reads'] as Long
}

def getFastpAdapterSequence(json_file) {
    // Handle stub runs
    if (workflow.stubRun) {
        return ""
    }

    def json = new groovy.json.JsonSlurper().parseText(json_file.text) as Map
    try {
        return (json['adapter_cutting'] as Map)['read1_adapter_sequence'] as String
    }
    catch (Exception _ex) {
        return ""
    }
}

workflow FASTQ_FASTQC_UMITOOLS_FASTP {
    take:
    ch_reads: Channel<Tuple> // (meta, reads, adapter_fasta)
    skip_fastqc: Boolean // true/false
    with_umi: Boolean // true/false
    skip_umi_extract: Boolean // true/false
    umi_discard_read: Integer // 0, 1 or 2
    skip_trimming: Boolean // true/false
    save_trimmed_fail: Boolean // true/false
    save_merged: Boolean // true/false
    min_trimmed_reads: Integer // > 0
    tool_args: Record // umi_extract, fastp: UMI-tools extract and fastp options

    main:
    // Each stage that does not run fills its elements with null. Joins on meta happen before
    // the UMI read discard changes it, or after on both sides.

    // (meta, reads)
    ch_raw_reads = ch_reads.map { meta, reads, _adapter_fasta -> tuple(meta, reads) }

    // (meta, reads, adapter_fasta, fastqc_raw_html, fastqc_raw_zip)
    ch_fastqc_raw = ch_reads.map { meta, reads, adapter_fasta -> tuple(meta, reads, adapter_fasta, null, null) }
    if (!skip_fastqc) {
        ch_fastqc_raw = ch_reads.join(FASTQC_RAW(ch_raw_reads), by: 0)
    }

    // (meta, reads, adapter_fasta, fastqc_raw_html, fastqc_raw_zip, umi(log, reads))
    ch_umi_extracted = ch_fastqc_raw.map { meta, reads, adapter_fasta, fastqc_raw_html, fastqc_raw_zip ->
        tuple(meta, reads, adapter_fasta, fastqc_raw_html, fastqc_raw_zip, null)
    }
    if (with_umi && !skip_umi_extract) {
        // (meta, reads, log)
        ch_umi = UMITOOLS_EXTRACT(ch_raw_reads.map { meta, reads -> tuple(meta, reads, tool_args.umi_extract) })
        ch_umi_extracted = ch_fastqc_raw
            .join(ch_umi, by: 0)
            .map { meta, _reads, adapter_fasta, fastqc_raw_html, fastqc_raw_zip, umi_reads, umi_log ->
                // Discard R1 / R2 if required
                def discard = umi_discard_read in [1, 2] && !meta.single_end
                tuple(
                    discard ? meta + [single_end: true] : meta,
                    discard ? [umi_reads[umi_discard_read % 2]] : umi_reads,
                    adapter_fasta,
                    fastqc_raw_html,
                    fastqc_raw_zip,
                    tuple(umi_log, umi_reads)
                )
            }
    }

    // (meta, reads, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi(log, reads),
    //  trim(html, json, log, reads_fail, reads_merged), adapter_seq, num_trimmed_reads)
    ch_results = ch_umi_extracted.map { meta, reads, _adapter_fasta, fastqc_raw_html, fastqc_raw_zip, umi ->
        tuple(meta, reads, fastqc_raw_html, fastqc_raw_zip, null, null, umi, null, null, null)
    }

    if (!skip_trimming) {
        //
        // Filter FastQ files based on minimum trimmed read count after adapter trimming
        //
        ch_fastp_input = ch_umi_extracted.map { meta, reads, adapter_fasta, _fastqc_raw_html, _fastqc_raw_zip, _umi ->
            tuple(meta, reads, adapter_fasta, tool_args.fastp, false, save_trimmed_fail, save_merged)
        }
        // (meta, reads, trim(html, json, log, reads_fail, reads_merged), adapter_seq, num_trimmed_reads)
        ch_trim = FASTP(ch_fastp_input).map { meta, reads, json, html, log, reads_fail, reads_merged ->
            // FASTP reads are optional, so a sample can lack a read count
            def num_reads = reads.isEmpty() ? null : getFastpReadsAfterFiltering(json, min_trimmed_reads as Long)
            tuple(
                meta,
                num_reads != null && num_reads >= (min_trimmed_reads as Long) ? reads : null,
                tuple(html, json, log, reads_fail.isEmpty() ? null : reads_fail, reads_merged),
                getFastpAdapterSequence(json),
                num_reads
            )
        }
        // (meta, reads, fastqc_raw_html, fastqc_raw_zip, umi(log, reads), trim(html, json, log, reads_fail, reads_merged), adapter_seq, num_trimmed_reads)
        ch_trimmed = ch_umi_extracted
            .join(ch_trim, by: 0)
            .map { meta, _reads, _adapter_fasta, fastqc_raw_html, fastqc_raw_zip, umi, trimmed_reads, trim, adapter_seq, num_trimmed_reads ->
                tuple(meta, trimmed_reads, fastqc_raw_html, fastqc_raw_zip, umi, trim, adapter_seq, num_trimmed_reads)
            }

        if (!skip_fastqc) {
            // Samples below min_trimmed_reads are not passed to FASTQC_TRIM
            // (meta, html, zip)
            ch_fastqc_trim = FASTQC_TRIM(
                ch_trimmed
                    .filter { _meta, reads, _fastqc_raw_html, _fastqc_raw_zip, _umi, _trim, _adapter_seq, _num_trimmed_reads -> reads != null }
                    .map { meta, reads, _fastqc_raw_html, _fastqc_raw_zip, _umi, _trim, _adapter_seq, _num_trimmed_reads -> tuple(meta, reads) }
            )

            // Every FASTQC_TRIM sample is in ch_trimmed, so only the FastQC side can be missing
            // from the remainder join, and its two outputs are nested to keep the tuple width fixed.
            ch_results = ch_trimmed
                .join(ch_fastqc_trim.map { meta, html, zip -> tuple(meta, tuple(html, zip)) }, by: 0, remainder: true)
                .map { meta, reads, fastqc_raw_html, fastqc_raw_zip, umi, trim, adapter_seq, num_trimmed_reads, fastqc_trim ->
                    tuple(
                        meta,
                        reads,
                        fastqc_raw_html,
                        fastqc_raw_zip,
                        fastqc_trim != null ? fastqc_trim[0] : null,
                        fastqc_trim != null ? fastqc_trim[1] : null,
                        umi,
                        trim,
                        adapter_seq,
                        num_trimmed_reads
                    )
                }
        }
        else {
            ch_results = ch_trimmed.map { meta, reads, fastqc_raw_html, fastqc_raw_zip, umi, trim, adapter_seq, num_trimmed_reads ->
                tuple(meta, reads, fastqc_raw_html, fastqc_raw_zip, null, null, umi, trim, adapter_seq, num_trimmed_reads)
            }
        }
    }

    emit:
    // (meta, reads, fastqc_raw_html, fastqc_raw_zip, fastqc_trim_html, fastqc_trim_zip, umi(log, reads),
    //  trim(html, json, log, reads_fail, reads_merged), adapter_seq, num_trimmed_reads)
    ch_results
}
