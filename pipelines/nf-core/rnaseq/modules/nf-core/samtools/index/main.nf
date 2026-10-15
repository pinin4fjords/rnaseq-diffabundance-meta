nextflow.enable.types = true

process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e994bf4eb3731150511a14f5706b7bdfd64df1b6d40898fff334286c027e0859/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'}"

    input:
    tuple(
        meta: Map,
        bam:  Path,
        args: String?
    )

    output:
    // (meta, bai)
    tuple(meta, file('*.{bai,csi,crai}'))

    topic:
    tuple(task.process, 'samtools', eval("samtools version | sed '1!d;s/.* //'")) >> 'versions'

    script:
    def tool_args = task.ext.args ?: args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus} \\
        ${tool_args} \\
        ${bam}
    """

    stub:
    def tool_args = task.ext.args ?: ''
    def extension = bam.getExtension() == 'cram'
        ? "crai"
        : tool_args.contains("-c") ? "csi" : "bai"
    """
    touch ${bam}.${extension}
    """
}
