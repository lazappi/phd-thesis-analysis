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


#' Read Alevin matrix
#'
#' Read Alevin data into a matrix
#'
#' @param path Directory containing Alevin output
#' @param type Whether to read quant or tier matrix
#' @param verbose Whether to print progress messages?
#'
#' @details
#' Based on tximport::readAlevin
#'
#' @return Sparse matrix containing Alevin data
readAlevinMat <- function(path, type = c("quant", "tier"), verbose = FALSE) {

    type <- match.arg(type)

    if (verbose) {message("Locating files...")}
    barcode_file <- file.path(path, "quants_mat_rows.txt")
    gene_file <- file.path(path, "quants_mat_cols.txt")

    switch(type,
           quant = {
               matrix_file <- file.path(path, "quants_mat.gz")
           },
           tier = {
               matrix_file <- file.path(path, "quants_tier_mat.gz")
           })

    for (f in c(barcode_file, gene_file, matrix_file)) {
        if (!file.exists(f)) {
            stop("Missing file ", f, call. = FALSE)
        }
    }

    if (verbose) {message("Reading cell barcodes...")}
    barcodes <- readLines(barcode_file)

    if (verbose) {message("Reading gene names...")}
    genes <- readLines(gene_file)

    n_cells <- length(barcodes)
    n_genes <- length(genes)
    mat <- matrix(nrow = n_genes, ncol = n_cells,
                  dimnames = list(genes, barcodes))

    if (verbose) {message("Reading matrix...")}
    mat_con <- gzcon(file(matrix_file, "rb"))
    for (j in seq_len(n_cells)) {
        mat[, j] <- readBin(mat_con, double(), endian = "little", n = n_genes)
    }
    close(mat_con)

    if (verbose) {message("Coverting to sparse matrix...")}
    mat <- as(mat, "dgCMatrix")

    if (verbose) {message("Done!")}
    return(mat)
}


#' Read Alevin
#'
#' Read Alevin data into a SingleCellExperiment object
#'
#' @param paths Directories containing Alevin output
#' @param dataset Name of the dataset
#' @param verbose Whether to print progress messages?
#'
#' @return SingleCellExperiment object containing Alevin data
readAlevin <- function(paths, dataset, verbose = FALSE) {

    mats <- lapply(seq_along(paths), function(idx) {
        if (verbose) {message("Reading sample ", idx, "...")}
        mat <- readAlevinMat(paths[idx], verbose = verbose)
        colnames(mat) <- paste(colnames(mat), idx, sep = "-")

        return(mat)
    })
    if (verbose) {message("Combining matrices...")}
    mat <- do.call("cbind", mats)

    cell_names <- colnames(mat)
    n_cells <- length(cell_names)
    cell_ids <- paste0("Cell", 1:n_cells)

    cell_info <- S4Vectors::DataFrame(Cell = cell_ids,
                                      Dataset = rep(dataset, n_cells),
                                      Barcode = stringr::str_sub(cell_names,
                                                                 end = -3),
                                      Sample = stringr::str_sub(cell_names,
                                                                start = -1),
                                      row.names = NULL)
    colnames(mat) <- cell_info$Cell

    gene_info <- S4Vectors::DataFrame(ID = rownames(mat))

    if (verbose) {message("Creating SingleCellExperiment...")}
    sce <- SingleCellExperiment::SingleCellExperiment(
        list(counts = mat),
        rowData = gene_info,
        colData = cell_info
    )

    if (verbose) {message("Done!")}
    return(sce)
}
