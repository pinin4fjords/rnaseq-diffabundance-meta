nextflow.enable.types = true

process TXIMETA_TXIMPORT {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bd/bdba33f8ad1b2df156f8f6775279fb217ce0f8233a3dc637337245de9ad2f29f/data' :
        'community.wave.seqera.io/library/bioconductor-tximeta_jq:78bccd386c46a07c' }"

    input:
    tuple(
        meta:       Map,
        quants:     List<Path>,
        tx2gene:    Path,
        quant_type: String
    )

    stage:
    stageAs quants, 'quants/*'

    output:
    // (meta, tpm_gene, counts_gene, lengths_gene, counts_gene_length_scaled, counts_gene_scaled, tpm_transcript, counts_transcript, lengths_transcript, tx2gene_augmented)
    tuple(
        meta,
        file("*gene_tpm.tsv"),
        file("*gene_counts.tsv"),
        file("*gene_lengths.tsv"),
        file("*gene_counts_length_scaled.tsv"),
        file("*gene_counts_scaled.tsv"),
        file("*transcript_tpm.tsv"),
        file("*transcript_counts.tsv"),
        file("*transcript_lengths.tsv"),
        file("*tx2gene_augmented.tsv")
    )

    topic:
    file('versions.yml') >> 'versions'

    script:
    template 'tximport.r'

    stub:
    """
    touch ${meta.id}.gene_tpm.tsv
    touch ${meta.id}.gene_counts.tsv
    touch ${meta.id}.gene_counts_length_scaled.tsv
    touch ${meta.id}.gene_counts_scaled.tsv
    touch ${meta.id}.gene_lengths.tsv
    touch ${meta.id}.transcript_tpm.tsv
    touch ${meta.id}.transcript_counts.tsv
    touch ${meta.id}.transcript_lengths.tsv
    touch ${meta.id}.tx2gene_augmented.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bioconductor-tximeta: \$(Rscript -e "library(tximeta); cat(as.character(packageVersion('tximeta')))")
    END_VERSIONS
    """
}
