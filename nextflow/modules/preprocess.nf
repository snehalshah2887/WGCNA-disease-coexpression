// =============================================================================
// Nextflow Module: PREPROCESS
// Runs scripts/01_preprocess.R
// =============================================================================

process PREPROCESS {
    label 'small'
    tag   "preprocessing"

    publishDir "${params.outdir}/01_preprocessed", mode: 'copy'

    input:
    path expression
    path r_script

    output:
    path "gene_expression_filtered.csv", emit: filtered_expr
    path "gene_annotation.rds",          emit: gene_annotation_rds

    script:
    """
    Rscript $r_script \\
        ${expression} \\
        .
    """
}
