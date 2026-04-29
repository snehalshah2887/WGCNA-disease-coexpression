#!/usr/bin/env nextflow
// =============================================================================
// 22q13 / Phelan-McDermid Syndrome WGCNA Pipeline — Main Workflow
// Nextflow DSL2
//
// Usage:
//   nextflow run nextflow/main.nf \
//     --expression data/input/gene_expression_filtered.csv \
//     --outdir     results/ \
//     -profile docker
//
// To skip the Cytoscape live-session steps (default):
//   --skip_cytoscape true
// =============================================================================

// ---------------------------------------------------------------------------
// Print help message
// ---------------------------------------------------------------------------
if (params.help) {
    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║       22q13 WGCNA Pipeline — Help                               ║
    ╚══════════════════════════════════════════════════════════════════╝

    REQUIRED PARAMETERS
      --expression        Path to combined gene expression CSV
                          (columns 1-10: annotation, 11-534: BrainSpan samples)
                          (e.g. data/input/gene_expression_filtered.csv)

    OPTIONAL PARAMETERS
      --outdir            Output directory (default: results/)
      --skip_cytoscape    Skip live Cytoscape network steps (default: true)
      --help              Show this message

    PROFILES
      -profile docker        Run inside the 22q13-wgcna Docker container
      -profile singularity   Run with Singularity
      -profile local         Run locally (all R packages must be installed)
      -profile slurm         SLURM HPC execution with Singularity

    EXAMPLE
      nextflow run nextflow/main.nf \\
        --expression data/input/gene_expression_filtered.csv \\
        -profile docker
    """.stripIndent()
    exit 0
}

// ---------------------------------------------------------------------------
// Validate required inputs
// ---------------------------------------------------------------------------
if (!params.expression) {
    error "Please provide --expression (path to combined gene expression CSV)"
}

// ---------------------------------------------------------------------------
// Import modules
// ---------------------------------------------------------------------------
include { PREPROCESS       } from './modules/preprocess'
include { BUILD_NETWORK    } from './modules/build_network'
include { EXPORT_NETWORK   } from './modules/export_network'
include { SUPPLEMENT_TABLE } from './modules/supplement_table'
include { GO_ANALYSIS      } from './modules/go_analysis'
include { GENE_RANKING     } from './modules/gene_ranking'

// ---------------------------------------------------------------------------
// Main workflow
// ---------------------------------------------------------------------------
workflow {

    log.info """
    ╔══════════════════════════════════════════════════════════════════╗
    ║         22q13 WGCNA Pipeline                                    ║
    ╚══════════════════════════════════════════════════════════════════╝
    Expression data : ${params.expression}
    GO enrichment   : AnnotationDbi + phyper (org.Hs.eg.db / GO.db)
    Output dir      : ${params.outdir}
    Skip Cytoscape  : ${params.skip_cytoscape}
    """.stripIndent()

    // Input data channels
    ch_expression = Channel.fromPath(params.expression, checkIfExists: true)

    // R script channels — staged into each process work directory so they are
    // accessible both in Docker (mounted work dir) and in local execution
    ch_script_01 = Channel.fromPath("${projectDir}/../scripts/01_preprocess.R",    checkIfExists: true)
    ch_script_02 = Channel.fromPath("${projectDir}/../scripts/02_build_network.R", checkIfExists: true)
    ch_script_03 = Channel.fromPath("${projectDir}/../scripts/03_export_network.R",checkIfExists: true)
    ch_script_04 = Channel.fromPath("${projectDir}/../scripts/04_supplement_table.R", checkIfExists: true)
    ch_script_05 = Channel.fromPath("${projectDir}/../scripts/05_go_analysis.R",   checkIfExists: true)
    ch_script_06 = Channel.fromPath("${projectDir}/../scripts/06_gene_ranking.R",  checkIfExists: true)

    // -----------------------------------------------------------------------
    // Step 1: Preprocess
    // -----------------------------------------------------------------------
    PREPROCESS(ch_expression, ch_script_01)

    // -----------------------------------------------------------------------
    // Step 2: Build co-expression network, annotate modules, run enrichment
    // -----------------------------------------------------------------------
    BUILD_NETWORK(
        PREPROCESS.out.filtered_expr,
        PREPROCESS.out.gene_annotation_rds,
        ch_script_02
    )

    // -----------------------------------------------------------------------
    // Step 3: TOM calculation + Cytoscape network export
    // -----------------------------------------------------------------------
    EXPORT_NETWORK(
        BUILD_NETWORK.out.network_objects,
        Channel.value(params.skip_cytoscape),
        ch_script_03
    )

    // -----------------------------------------------------------------------
    // Step 4: Supplemental co-expression table
    // -----------------------------------------------------------------------
    SUPPLEMENT_TABLE(
        EXPORT_NETWORK.out.tom_objects,
        BUILD_NETWORK.out.network_objects,
        EXPORT_NETWORK.out.cytoscape_dir,
        ch_script_04
    )

    // -----------------------------------------------------------------------
    // Step 5: GO term analysis (AnnotationDbi + phyper — no external files needed)
    // -----------------------------------------------------------------------
    GO_ANALYSIS(
        BUILD_NETWORK.out.network_objects,
        ch_script_05
    )

    // -----------------------------------------------------------------------
    // Step 6: Gene ranking
    // -----------------------------------------------------------------------
    GENE_RANKING(
        BUILD_NETWORK.out.network_objects,
        EXPORT_NETWORK.out.tom_objects,
        ch_script_06
    )
}

workflow.onComplete {
    log.info """
    Pipeline completed!
    Status    : ${workflow.success ? 'SUCCESS' : 'FAILED'}
    Duration  : ${workflow.duration}
    Output dir: ${params.outdir}
    """.stripIndent()
}

workflow.onError {
    log.error """
    Pipeline FAILED.
    Error     : ${workflow.errorMessage}
    Work dir  : ${workflow.workDir}
    Tip       : cd to the failed process work dir and run: bash .command.run
    """.stripIndent()
}
