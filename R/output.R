#' Get download markdown link
#'
#' Convert an output file name and location to a URL that can be used to
#' download the file.
#'
#' @param file name of the output file
#' @param folder name of the directory in the output directory containing the
#' output file
#'
#' @return Markdown URL link to the file
getDownloadLink <- function(file, folder = NULL) {
    remote <- workflowr::wflow_git_remote(verbose = FALSE)["origin"]

    url <- gsub(":", "/", remote)
    url <- gsub("git@", "http://", url)
    url <- gsub(".git", "", url, fixed = TRUE)
    url <- paste(url, "raw/master/output", sep = "/")

    if (is.null(folder)) {
        url <- paste(url, file, sep = "/")
    } else {
        url <- paste(url, folder, file, sep = "/")
    }

    link <- glue::glue("[{file}]({url})")

    return(link)
}


#' Write gene table
#'
#' Save a gene table to a file
#'
#' @param gene_table gene table to save
#' @param path file path to save location
#'
#' @details
#' data.frame objects will be saved as a (zipped) CSV file. List objects will
#' be saved in XLSX format.
writeGeneTable <- function(gene_table, path) {

    if (is.data.frame(gene_table)) {
        zip_path <- paste0(path, ".zip")
        if (file.exists(zip_path)) {file.remove(zip_path)}
        readr::write_csv(gene_table, path, na = "")
        zip(zip_path, path, flags = "-q -j")
        invisible(file.remove(path))
    } else {
        writexl::write_xlsx(gene_table, path)
    }
}
