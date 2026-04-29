// =============================================================================
// Nextflow Module: BUILD_NETWORK
// Runs scripts/02_build_network.R
// This is the most resource-intensive step (labeled 'high_mem').
// =============================================================================

process BUILD_NETWORK {
    label 'high_mem'
    tag   "wgcna_network"

    publishDir "${params.outdir}/02_network", mode: 'copy'

    input:
    path filtered_expr
    path gene_annotation_rds
    path r_script

    output:
    path "network_objects.rds",    emit: network_objects
    path "colormodule.csv",        emit: colormodule
    path "summary_colors.csv",     emit: summary_colors
    path "color_modules/",         emit: color_modules_dir
    path "plots/",                 emit: plots_dir
    path "TOM/",                   emit: tom_dir

    script:
    """
    Rscript $r_script \\
        ${filtered_expr} \\
        ${gene_annotation_rds} \\
        .
    """
}
