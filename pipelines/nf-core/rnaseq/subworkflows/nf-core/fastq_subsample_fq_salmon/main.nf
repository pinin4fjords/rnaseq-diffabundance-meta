nextflow.enable.types = true

//
// Sub-sample FastQ files and pseudo-align with Salmon
//      can be used to infer strandedness of library
//

include { SALMON_INDEX } from '../../../modules/nf-core/salmon/index/main'
include { FQ_SUBSAMPLE } from '../../../modules/nf-core/fq/subsample/main'
include { SALMON_QUANT } from '../../../modules/nf-core/salmon/quant/main'

workflow FASTQ_SUBSAMPLE_FQ_SALMON {
    take:
    ch_samples: Channel<Tuple> // (meta, reads)
    ch_genome_fasta: Value<Path?> // decoys for the Salmon index, absent when the genome is not provided
    ch_transcript_fasta: Value<Path>
    ch_gtf: Value<Path>
    ch_index: Value<Path?> // prebuilt Salmon index, used when make_index is false
    make_index: Boolean // Whether to create salmon index before running salmon quant
    tool_args: Record // salmon_index_strandedness: Salmon index options for the index of the strandedness inference

    main:

    //
    // Create Salmon index if required
    //
    if (make_index) {
        ch_index_built = SALMON_INDEX(
            ch_transcript_fasta
                .combine(ch_genome_fasta)
                .map { transcript_fasta, genome_fasta ->
                    tuple([id: 'salmon_index'], transcript_fasta, genome_fasta, tool_args.salmon_index_strandedness)
                }
        )
        ch_index_ref = ch_index_built.map { _meta, index -> index }
        ch_index_file = ch_index_ref
    }
    else {
        ch_index_ref = ch_index
        ch_index_file = ch_index.map { _index -> null }
    }

    //
    // Sub-sample FastQ files with fq
    //
    // (meta, reads)
    ch_subsampled = FQ_SUBSAMPLE(ch_samples)

    //
    // Pseudo-alignment with Salmon
    //
    // (meta, quant_dir, json_info, lib_format_counts)
    ch_quant = SALMON_QUANT(
        ch_subsampled
            .map { meta, reads -> tuple(meta, reads, null) }
            .combine(ch_index_ref)
            .combine(ch_gtf)
            .combine(ch_transcript_fasta)
    )

    // (meta, reads, quant_dir, json_info, lib_format_counts)
    ch_results = ch_subsampled.join(ch_quant, by: 0)

    emit:
    samples:     Channel<Tuple> = ch_results // (meta, reads, quant_dir, json_info, lib_format_counts)
    index_built: Value<Path?>   = ch_index_file // null unless the index was built here
}
