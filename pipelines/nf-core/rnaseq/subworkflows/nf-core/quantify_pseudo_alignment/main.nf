nextflow.enable.types = true

//
// Pseudoalignment and quantification with Salmon or Kallisto
//

include { SALMON_QUANT     } from '../../../modules/nf-core/salmon/quant'
include { KALLISTO_QUANT   } from '../../../modules/nf-core/kallisto/quant'

include { QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT } from '../quant_tximport_summarizedexperiment'

//
// Salmon quant options: the library type is the given one if it is valid, else derived from the sample strandedness
//
def salmonQuantArgs(libtype: String?, extra_args: String?, meta: Map) -> String {
    def valid_libtypes = [ 'A', 'U', 'SF', 'SR', 'IS', 'IU', 'ISF', 'ISR', 'OS', 'OU', 'OSF', 'OSR', 'MS', 'MU', 'MSF', 'MSR' ]
    def lib_type = 'A'
    if (libtype) {
        if (valid_libtypes.contains(libtype)) {
            lib_type = libtype
        }
        else {
            log.info("[Salmon Quant] Invalid library type specified '--libType=${libtype}', defaulting to auto-detection with '--libType=A'.")
        }
    }
    else {
        def strandedness = meta.strandedness
        lib_type = meta.single_end
            ? (strandedness == 'forward' ? 'SF' : strandedness == 'reverse' ? 'SR' : 'U')
            : (strandedness == 'forward' ? 'ISF' : strandedness == 'reverse' ? 'ISR' : 'IU')
    }
    return [ "--libType=${lib_type}", extra_args ?: '' ].join(' ').trim()
}

//
// Kallisto quant options: the strandedness of the sample unless the extra options set one
//
def kallistoQuantArgs(extra_args: String?, meta: Map) -> String {
    def extra = extra_args ?: ''
    def strandedness = (extra.contains('--fr-stranded') || extra.contains('--rf-stranded')) ? '' :
        meta.strandedness == 'forward' ? '--fr-stranded' :
        meta.strandedness == 'reverse' ? '--rf-stranded' : ''
    return [ strandedness, extra ].join(' ').trim()
}

workflow QUANTIFY_PSEUDO_ALIGNMENT {
    take:
    samplesheet: Value<Path>
    ch_samples: Channel<Tuple> // (meta, reads)
    index: Value<Path?>
    transcript_fasta: Value<Path?>
    gtf: Value<Path>
    gtf_id_attribute: String // GTF gene ID attribute
    gtf_extra_attribute: String // GTF alternative gene attribute (e.g. gene_name)
    pseudo_aligner: String // kallisto or salmon
    kallisto_quant_fraglen: Integer? // Estimated fragment length required by Kallisto in single-end mode
    kallisto_quant_fraglen_sd: Integer? // Estimated standard error for fragment length required by Kallisto in single-end mode
    skip_merge: Boolean // skip cross-sample merging, run tximport per-sample
    salmon_libtype: String? // Salmon library type; derived from the sample strandedness when null
    tool_args: Record // extra Salmon and Kallisto quant options (salmon_quant, kallisto_quant) and the prefix of the merged SummarizedExperiment files (se_prefix, default: the sample id)

    main:

    //
    // Quantify and merge counts across samples
    //
    // ch_salmon: (meta, quant_dir, json_info, lib_format_counts), ch_kallisto: (meta, quant_dir, json_info, log)
    def ch_salmon: Channel<Tuple>   = channel.empty()
    def ch_kallisto: Channel<Tuple> = channel.empty()
    if (pseudo_aligner == 'salmon') {
        ch_salmon = SALMON_QUANT(
            ch_samples
                .map { meta, reads -> tuple(meta, reads, salmonQuantArgs(salmon_libtype, tool_args.salmon_quant, meta)) }
                .combine(index)
                .combine(gtf)
                .combine(transcript_fasta)
        )
    } else {
        // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
        ch_kallisto = KALLISTO_QUANT(
            ch_samples
                .map { t -> tuple(t) }
                .combine(index)
                .combine(gtf)
                .map { meta, reads, index_dir, gtf_file ->
                    tuple(meta, reads, index_dir, gtf_file, null, kallisto_quant_fraglen, kallisto_quant_fraglen_sd, kallistoQuantArgs(tool_args.kallisto_quant, meta))
                }
        )
    }
    // (meta, quants)
    def ch_quant_dirs = ch_salmon.mix(ch_kallisto).map { meta, quant_dir, _json_info, _log_or_counts -> tuple(meta, [ quant_dir ]) }

    //
    // Post-process quantifications with tximport and SummarizedExperiment
    //
    ch_quant_merged = QUANT_TXIMPORT_SUMMARIZEDEXPERIMENT(
        samplesheet,
        ch_quant_dirs,
        gtf,
        gtf_id_attribute,
        gtf_extra_attribute,
        pseudo_aligner,
        skip_merge,
        tool_args
    )

    emit:
    // per sample, when pseudo_aligner is salmon: (meta, quant_dir, json_info, lib_format_counts)
    salmon:   Channel<Tuple> = ch_salmon
    // per sample, when pseudo_aligner is kallisto: (meta, quant_dir, json_info, log)
    kallisto: Channel<Tuple> = ch_kallisto
    // one row per sample under skip_merge, a single 'all_samples' row otherwise:
    // (meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, tpm_transcript, counts_transcript, lengths_transcript, tx2gene_augmented, merged_gene_rds, merged_transcript_rds, tx2gene)
    merged:   Channel<Tuple> = ch_quant_merged
}
