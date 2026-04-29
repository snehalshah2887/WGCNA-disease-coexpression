// =============================================================================
// Nextflow Module: GENE_RANKING
// Runs scripts/06_gene_ranking.R
//
// Produces one ranked CSV per disease pair (PMS vs ASD, ID, Seizures,
// Hypotonia, LangImp) and one combined wide-format summary.
// =============================================================================

process GENE_RANKING {
    label 'medium'
    tag   "gene_ranking"

    publishDir "${params.outdir}/06_gene_ranking", mode: 'copy'

    input:
    path network_objects
    path tom_objects
    path r_script

    output:
    path "ranked_*_PMS.csv",             emit: ranked_by_disease
    path "ranked_all_diseases_PMS.csv",  emit: ranked_summary
    path "top5_PMS_disease_genes.csv",   emit: top5_summary

    script:
    """
    Rscript $r_script \\
        ${network_objects} \\
        ${tom_objects} \\
        .
    """
}
