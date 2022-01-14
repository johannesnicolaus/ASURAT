#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Add metadata of variables and samples.
#'
#' This function adds metadata of variables and samples.
#'
#' @param sce A SingleCellExperiment object.
#' @param mitochondria_symbol A string representing for mitochondrial genes.
#'   This function computes percents of reads that map to the mitochondrial
#'   genes. Examples are `^MT-`, `^mt-`, etc.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @export
#'
add_metadata_asrt <- function(sce, mitochondria_symbol){
  mat <- assay(sce, "counts")
  #--------------------------------------------------
  # Variable metadata
  #--------------------------------------------------
  rowData(sce)$nSamples <- apply(mat, 1, function(x) sum(x > 0))
  #--------------------------------------------------
  # Sample metadata
  #--------------------------------------------------
  sce$nReads <- as.integer(apply(mat, 2, function(x) sum(x)))
  sce$nGenes <- as.integer(apply(mat, 2, function(x) sum(x > 0)))
  sce$percMT <- apply(mat[grepl(mitochondria_symbol, rownames(mat)), ],
                      2, function(x) 100 * sum(x)) /
                  apply(mat, 2, function(x) sum(x))

  return(sce)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Remove variables based on expression profiles across samples.
#'
#' This function removes low expressed variable data.
#'
#' @param sce A SingleCellExperiment object.
#' @param min_nsamples An integer. This function removes variables for which
#'   the numbers of non-zero expressing samples are less than this value.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @export
#'
remove_variables_asrt <- function(sce, min_nsamples){
  mat <- assay(sce, "counts")
  inds <- which(apply(mat, 1, function(x) sum(x > 0)) >= min_nsamples)

  return(sce[inds, ])
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Remove samples based on expression profiles across variables.
#'
#' This function removes sample data by setting minimum and maximum threshold
#'   values for the metadata.
#'
#' @param sce A SingleCellExperiment object.
#' @param min_nReads A minimum threshold value of the number of reads.
#' @param max_nReads A maximum threshold value of the number of reads.
#' @param min_nGenes A minimum threshold value of the number of non-zero
#'   expressed genes.
#' @param max_nGenes A maximum threshold value of the number of non-zero
#'   expressed genes.
#' @param min_percMT A minimum threshold value of the percent of reads that map
#'   to mitochondrial genes.
#' @param max_percMT A maximum threshold value of the percent of reads that map
#'   to mitochondrial genes.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @export
#'
remove_samples_asrt <- function(
  sce, min_nReads = 0, max_nReads = 1e+10, min_nGenes = 0, max_nGenes = 1e+10,
  min_percMT = 0, max_percMT = 100
){
  mat <- assay(sce, "counts")
  inds_1 <- which((sce$nReads >= min_nReads) & (sce$nReads <= max_nReads))
  inds_2 <- which((sce$nGenes >= min_nGenes) & (sce$nGenes <= max_nGenes))
  inds_3 <- which((sce$percMT >= min_percMT) & (sce$percMT <= max_percMT))
  inds <- intersect(intersect(inds_1, inds_2), inds_3)

  return(sce[, inds])
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Remove variables based on the mean expression levels across samples.
#'
#' This function removes variable data such that the mean expression levels
#'   across samples are less than `min_meannReads`.
#'
#' @param sce A SingleCellExperiment object.
#' @param min_meannReads An integer. This function removes variables for which
#'   the mean read counts are less than this value.
#'
#' @return A SingleCellExperiment object.
#' @import SingleCellExperiment
#' @import SummarizedExperiment
#' @export
#'
remove_variables_second <- function(sce, min_meannReads){
  mat <- assay(sce, "counts")
  inds <- which(apply(mat, 1, mean) >= min_meannReads)

  return(sce[inds, ])
}
