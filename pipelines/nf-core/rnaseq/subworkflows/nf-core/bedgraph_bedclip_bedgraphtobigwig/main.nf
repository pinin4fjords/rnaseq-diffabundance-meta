nextflow.enable.types = true

//
// Run bedClip and bedGraphToBigWig
//

include { UCSC_BEDCLIP          } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG {
    take:
    ch_bedgraph: Channel<Record>
    sizes: Value<Path> // chrom.sizes

    main:

    //
    // Clip bedGraph file
    //
    ch_clipped = UCSC_BEDCLIP(ch_bedgraph.map { r -> record(meta: r.meta, bedgraph: r.bedgraph) }.combine(sizes: sizes))

    //
    // Convert bedGraph to bigWig
    //
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG(ch_clipped.map { r -> record(meta: r.meta, bedgraph: r.bedgraph) }.combine(sizes: sizes))

    ch_results = ch_clipped.join(ch_bigwig, by: 'meta')

    emit:
    ch_results
}
