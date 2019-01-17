#' Read 10x data
#'
#' Read 10x data into a SingleCellExperiment object
#'
#' @param path Directory containing Cell Ranger counts matrix
#' @param dataset Name of the dataset
#' @param verbose Whether to print progress messages?
#'
#' @details
#' Based on DropletUtils::read10xCounts and adapted for Cell Ranger 3 output
#'
#' @return SingleCellExperiment object containing 10x data
read10x <- function(path, dataset, verbose = FALSE) {

    if (verbose) {message("Reading gene info...")}
    gene_info <- readr::read_tsv(file.path(path, "features.tsv.gz"),
                                 col_names = c("ID", "Symbol", "Type"),
                                 col_types = readr::cols(
                                     .default = readr::col_character()
                                 ))
    gene_info <- S4Vectors::DataFrame(gene_info)
    gene_info$Name <- scater::uniquifyFeatureNames(gene_info$ID,
                                                   gene_info$Symbol)
    rownames(gene_info) <- gene_info$Name

    if (verbose) {message("Reading cell info...")}
    cell_names <- readr::read_lines(file.path(path, "barcodes.tsv.gz"))
    n_cells <- length(cell_names)
    cell_ids <- paste0("Cell", 1:n_cells)

    cell_info <- S4Vectors::DataFrame(Cell = cell_ids,
                                      Dataset = rep(dataset, n_cells),
                                      Barcode = stringr::str_sub(cell_names,
                                                                 end = -3),
                                      Sample = stringr::str_sub(cell_names,
                                                                start = -1),
                                      row.names = NULL)

    if (verbose) {message("Reading expression matrix...")}
    mat <- as(Matrix::readMM(file.path(path, "matrix.mtx.gz")), "dgCMatrix")
    colnames(mat) <- cell_info$Cell

    if (verbose) {message("Creating SingleCellExperiment...")}
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = mat),
        rowData = gene_info,
        colData = cell_info
    )

    if (verbose) {message("Done!")}
    return(sce)
}
