nextflow.enable.types = true

include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'

workflow BAM_STRINGTIE_MERGE {
    take:
    ch_bams: Channel<Tuple> // (meta, bam, lrbam)
    mode: Value<List<String>>
    chrgtf: Value<Path?>

    main:

    // (meta, transcript_gtf, abundance, coverage_gtf, ballgown)
    ch_assemblies = STRINGTIE_STRINGTIE(
        ch_bams
            .combine(mode)
            .combine(chrgtf)
    )

    // (meta, gtf, assemblies), assemblies holding every STRINGTIE_STRINGTIE output tuple
    ch_to_merge = ch_assemblies
        .collect()
        .flatMap { assemblies ->
            assemblies.isEmpty() ? [] : [
                tuple(
                    [id: 'stringtie_merge'],
                    assemblies.collect { _meta, transcript_gtf, _abundance, _coverage_gtf, _ballgown -> transcript_gtf }.toSorted { f -> f.name },
                    assemblies
                )
            ]
        }

    // (meta, merged_gtf)
    ch_merged = STRINGTIE_MERGE(ch_to_merge.map { meta, gtf, _assemblies -> tuple(meta, gtf) }.combine(chrgtf))

    // (meta, gtf, assemblies, merged_gtf)
    ch_results = ch_to_merge.join(ch_merged, by: 0)

    emit:
    ch_results
}
