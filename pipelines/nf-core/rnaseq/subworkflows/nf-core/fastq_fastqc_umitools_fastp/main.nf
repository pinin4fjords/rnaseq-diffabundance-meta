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
    ch_reads: Channel // records with meta, reads and adapter_fasta; their other fields are carried through to the output
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
    // Each stage that runs joins its outputs onto the incoming record, overwriting the null
    // placeholders of the fields it owns. Joins on meta happen before the UMI read discard
    // changes it, or after on both sides.
    ch_results = ch_reads.map { r ->
        r + record(
            fastqc_raw_html:   null,
            fastqc_raw_zip:    null,
            fastqc_trim_html:  null,
            fastqc_trim_zip:   null,
            umi:               null,
            trim:              null,
            adapter_seq:       null,
            num_trimmed_reads: null
        )
    }

    if (!skip_fastqc) {
        ch_fastqc_raw = FASTQC_RAW(ch_reads.map { r -> record(meta: r.meta, reads: r.reads) })
        ch_results = ch_results.join(
            ch_fastqc_raw.map { r -> record(meta: r.meta, fastqc_raw_html: r.html, fastqc_raw_zip: r.zip) },
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
        ch_fastp_input = ch_results.map { r ->
            record(
                meta:                 r.meta,
                reads:                r.reads,
                adapter_fasta:        r.adapter_fasta,
                args:                 tool_args.fastp,
                discard_trimmed_pass: false,
                save_trimmed_fail:    save_trimmed_fail,
                save_merged:          save_merged
            )
        }
        ch_trim = FASTP(ch_fastp_input).map { r ->
            // FASTP reads are optional, so a sample can lack a read count
            def num_reads = r.reads.isEmpty() ? null : getFastpReadsAfterFiltering(r.json, min_trimmed_reads as Long)
            record(
                meta:              r.meta,
                reads:             num_reads != null && num_reads >= (min_trimmed_reads as Long) ? r.reads : null,
                trim:              record(
                    html:         r.html,
                    json:         r.json,
                    log:          r.log,
                    reads_fail:   r.reads_fail.isEmpty() ? null : r.reads_fail,
                    reads_merged: r.reads_merged
                ),
                adapter_seq:       getFastpAdapterSequence(r.json),
                num_trimmed_reads: num_reads
            )
        }
        ch_results = ch_results.join(ch_trim, by: 'meta')

        if (!skip_fastqc) {
            ch_fastqc_trim = FASTQC_TRIM(ch_results.filter { r -> r.reads != null }.map { r -> record(meta: r.meta, reads: r.reads) })

            // Samples below min_trimmed_reads are not passed to FASTQC_TRIM
            ch_results = ch_results.join(
                ch_fastqc_trim.map { r -> record(meta: r.meta, fastqc_trim_html: r.html, fastqc_trim_zip: r.zip) },
                by: 'meta',
                remainder: true
            )
        }
    }

    emit:
    ch_results
}
