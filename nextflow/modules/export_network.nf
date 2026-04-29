// =============================================================================
// Nextflow Module: EXPORT_NETWORK
// Runs scripts/03_export_network.R
// =============================================================================

process EXPORT_NETWORK {
    label 'high_mem'
    tag   "tom_cytoscape"

    publishDir "${params.outdir}/03_network_export", mode: 'copy'

    input:
    path network_objects
    val  skip_cytoscape
    path r_script

    output:
    path "tom_objects.rds",    emit: tom_objects
    path "cytoscape_files/",   emit: cytoscape_dir

    script:
    """
    Rscript $r_script \\
        ${network_objects} \\
        . \\
        ${skip_cytoscape}
    """
}
