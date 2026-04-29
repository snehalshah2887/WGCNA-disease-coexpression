// =============================================================================
// Nextflow Module: SUPPLEMENT_TABLE
// Runs scripts/04_supplement_table.R
// =============================================================================

process SUPPLEMENT_TABLE {
    label 'small'
    tag   "supplement_table"

    publishDir "${params.outdir}/04_supplement_table", mode: 'copy'

    input:
    path tom_objects
    path network_objects
    path cytoscape_dir
    path r_script

    output:
    path "Table22q13.csv", emit: supplement_table

    script:
    """
    Rscript $r_script \\
        ${tom_objects} \\
        ${network_objects} \\
        ${cytoscape_dir} \\
        .
    """
}
