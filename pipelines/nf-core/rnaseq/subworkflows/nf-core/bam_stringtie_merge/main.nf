nextflow.enable.types = true

include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'

workflow BAM_STRINGTIE_MERGE {
    take:
    ch_bams: Channel<Record>
    mode: Value<List<String>>
    chrgtf: Value<Path?>

    main:

    ch_assemblies = STRINGTIE_STRINGTIE(
        ch_bams
            .map { r -> record(meta: r.meta, bam: r.bam, lrbam: r.lrbam) }
            .combine(mode: mode, annotation_gtf: chrgtf)
    )

    ch_to_merge = ch_assemblies
        .collect()
        .flatMap { assemblies ->
            assemblies.isEmpty() ? [] : [
                record(
                    meta:       record(id: 'stringtie_merge'),
                    gtf:        assemblies.collect { r -> r.transcript_gtf }.toSorted { f -> f.name },
                    assemblies: assemblies
                )
            ]
        }

    ch_merged = STRINGTIE_MERGE(ch_to_merge.map { r -> record(meta: r.meta, gtf: r.gtf) }.combine(annotation_gtf: chrgtf))

    ch_results = ch_to_merge.join(ch_merged.map { r -> record(meta: r.meta, merged_gtf: r.merged_gtf) }, by: 'meta')

    emit:
    ch_results
}
