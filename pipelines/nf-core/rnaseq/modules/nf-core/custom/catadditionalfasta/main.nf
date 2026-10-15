nextflow.enable.types = true

process CUSTOM_CATADDITIONALFASTA {
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'quay.io/biocontainers/python:3.12' }"

    input:
    record(
        meta:      Record,
        fasta:     Path,
        gtf:       Path,
        prefix:    String?,
        add_fasta: Path,
        biotype:   String
    )

    output:
    record(
        meta:  meta,
        fasta: file("out/${out_prefix}.fasta"),
        gtf:   file("out/${out_prefix}.gtf")
    )

    topic:
    file('versions.yml') >> 'versions'

    script:
    out_prefix = task.ext.prefix ?: prefix ?: "${meta.id}"

    template 'fasta2gtf.py'

    stub:
    out_prefix = task.ext.prefix ?: prefix ?: "${meta.id}"
    """
    mkdir out
    touch out/${out_prefix}.fasta
    touch out/${out_prefix}.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
