nextflow.enable.types = true

process CUSTOM_MULTIQCCUSTOMBIOTYPE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12.12' :
        'quay.io/biocontainers/python:3.12.12' }"

    input:
    tuple(
        meta:   Map,
        counts: Path,
        header: Path
    )

    output:
    // (meta, tsv, rrna)
    tuple(
        meta,
        file('*biotype_counts_mqc.tsv'),
        file('*biotype_counts_rrna_mqc.tsv')
    )

    topic:
    file('versions.yml') >> 'versions'

    script:
    template 'mqc_features_stat.py'

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.biotype_counts_mqc.tsv
    touch ${prefix}.biotype_counts_rrna_mqc.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
