#' Annotate SCE
#'
#' Add additional metadata to a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment object
#' @param org Whether it is human or mouse data
#' @param add_anno Whether to add feature annotation?
#' @param host BioMart host to use
#' @param calc_qc Whether to calculate QC metrics?
#' @param calc_cpm Whether to calculate Counts Per Million?
#' @param cell_cycle Whether to assign cell cycle?
#' @param BPPARAM A BiocParallelParam object to use for parallel processing
#' @param verbose Whether to print progress messages?
#'
#' @return SingleCellExperiment object with additional annotation
annotateSCE <- function(sce, org = c("human", "mouse"), add_anno = FALSE,
                        host = "www.ensembl.org", calc_qc = FALSE,
                        calc_cpm = FALSE, cell_cycle = FALSE,
                        BPPARAM = BiocParallel::SerialParam(),
                        verbose = FALSE) {

    if (add_anno) {
        if (verbose) {message("Adding feature annotation...")}
        sce <- addFeatureAnnos(sce, org, host)
    }

    if (cell_cycle) {
        if (verbose) {message("Assigning cell cycle phases...")}
        sce <- addCellCycle(sce, org, BPPARAM, verbose)
    }

    if (calc_cpm) {
        if (verbose) {message("Calculating CPM...")}
        SingleCellExperiment::cpm(sce) <- scater::calculateCPM(sce,
                                                               use_size_factors = FALSE)
    }

    if (calc_qc) {
        if (verbose) {message("Calculating QC metrics...")}
        if ("chromosome_name" %in%
            colnames(SummarizedExperiment::rowData(sce))) {
            is_mt <- SummarizedExperiment::rowData(sce)$chromosome_name == "MT"
            sce <- scater::calculateQCMetrics(sce,
                                              feature_controls = list(MT = is_mt),
                                              BPPARAM = BPPARAM)
        } else {
            sce <- scater::calculateQCMetrics(sce, BPPARAM = BPPARAM)
        }
    }

    return(sce)
}


#' Add feature annotations
#'
#' Add feature annotations to an SCE object. Just calls the
#' \code{getBMFeatureAnnos} function in \code{scater} with a specific set of
#' parameters.
#'
#' @param sce SCE object to add annotations to.
#' @param org Whether it is human or mouse data
#' @param host BioMart host to use
#'
#' @return SCE with addtional feature annotations.
addFeatureAnnos <- function(sce, org, host) {

    if (org == "human") {
        dataset <- "hsapiens_gene_ensembl"
        symbol <- "hgnc_symbol"
    } else if (org == "mouse") {
        dataset <- "mmusculus_gene_ensembl"
        symbol <- "mgi_symbol"
    } else {
        stop("org should be 'human' or 'mouse'")
    }

    sce <- scater::getBMFeatureAnnos(sce,
                                     ids <- SummarizedExperiment::rowData(sce)$ID,
                                     filters = "ensembl_gene_id",
                                     attributes = c("ensembl_gene_id",
                                                    "entrezgene",
                                                    "external_gene_name",
                                                    symbol,
                                                    "chromosome_name",
                                                    "description",
                                                    "gene_biotype",
                                                    "percentage_gene_gc_content"),
                                     dataset = dataset,
                                     host = host)

    SummarizedExperiment::rowData(sce)$description <- gsub("\\s\\[.*\\]", "",
                                                           SummarizedExperiment::rowData(sce)$description)

    return(sce)
}


#' Add cell cycle
#'
#' Add cell cycle phase to cells in an SCE object. Calculated using the
#' \code{\link[scran]{cyclone}} function.
#'
#' @param sce SCE object.
#' @param org Whether it is human or mouse data
#' @param BPPARAM A BiocParallelParam object to use for parallel processing
#' @param verbose Whether to print progress messages?
#'
#' @return SCE object with assigned cell cycles
addCellCycle <- function(sce, org, BPPARAM, verbose) {

    if (org == "human") {
        cc_pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds",
                                        package = "scran"))
    } else {
        cc_pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds",
                                        package = "scran"))
    }

    cycles <- scran::cyclone(sce, pairs = cc_pairs,
                             gene.names = SummarizedExperiment::rowData(sce)$ID,
                             BPPARAM = BPPARAM,
                             verbose = verbose)
    SummarizedExperiment::colData(sce)$G1Score <- cycles$scores$G1
    SummarizedExperiment::colData(sce)$SScore <- cycles$scores$G1
    SummarizedExperiment::colData(sce)$G2MScore <- cycles$scores$G2M
    SummarizedExperiment::colData(sce)$CellCycle <- cycles$phases

    return(sce)
}
