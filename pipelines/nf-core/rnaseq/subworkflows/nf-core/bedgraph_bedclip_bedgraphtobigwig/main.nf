nextflow.enable.types = true

//
// Run bedClip and bedGraphToBigWig
//

include { UCSC_BEDCLIP          } from '../../../modules/nf-core/ucsc/bedclip/main'
include { UCSC_BEDGRAPHTOBIGWIG } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow BEDGRAPH_BEDCLIP_BEDGRAPHTOBIGWIG {
    take:
    ch_bedgraph: Channel<Tuple> // (meta, bedgraph)
    sizes: Value<Path> // chrom.sizes

    main:

    //
    // Clip bedGraph file
    //
    // (meta, bedgraph)
    // Process outputs and join results are lists, which combine nests rather than extends; tuple(t) turns them into tuples.
    ch_clipped = UCSC_BEDCLIP(ch_bedgraph.map { t -> tuple(t) }.combine(sizes))

    //
    // Convert bedGraph to bigWig
    //
    // (meta, bigwig)
    ch_bigwig = UCSC_BEDGRAPHTOBIGWIG(ch_clipped.map { t -> tuple(t) }.combine(sizes))

    // (meta, bedgraph, bigwig)
    ch_results = ch_clipped.join(ch_bigwig, by: 0)

    emit:
    ch_results
}
