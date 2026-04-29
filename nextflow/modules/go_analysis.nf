// =============================================================================
// Nextflow Module: GO_ANALYSIS
// Runs scripts/05_go_analysis.R
// GO enrichment is computed via hypergeometric test (phyper) using
// org.Hs.eg.db and GO.db annotations. No external GO files needed.
// =============================================================================

process GO_ANALYSIS {
    label 'small'
    tag   "go_analysis"

    publishDir "${params.outdir}/05_go_analysis", mode: 'copy'

    input:
    path network_objects
    path r_script

    output:
    path "Goterms_22q13.csv",  emit: go_terms_22q13
    path "top5_go_all.csv",    emit: top5_go
    path "filtered_go/",       emit: filtered_go_dir
    path "plots/",             emit: plots_dir

    script:
    """
    Rscript $r_script \\
        ${network_objects} \\
        .
    """
}
