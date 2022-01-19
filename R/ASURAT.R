#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Functional annotation-driven unsupervised clustering of single-cell data.
#'
#' ASURAT is a computational pipeline for single-cell or spatial
#'    transcriptome data analyses, which simultaneously performs unsupervised
#'    clustering and biological interpretation in terms of cell type, disease,
#'    biological process, and signaling pathway activity.
#'    Inputting single-cell RNA-seq (scRNA-seq) data and knowledge-based
#'    databases, some of which have already been prepared, ASURAT transforms
#'    a read count table into an original multivariate table,
#'    termed a sign-by-sample matrix (SSM). By analyzing SSMs, users can
#'    cluster cells (or spatial spots) to aid their interpretation.
#'
#' @docType package
#' @name ASURAT
NULL
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Print future work.
#'
#' This function prints future work.
#'
#' @return A sentences.
#' @export
#'
print_future_asrt <- function(){
  message("* Create biological sentences by combining multiple signs.")
  message("  - Explore a way to characterize samples by sentences.")
  message("* Find a feedback structure in the system of linguistics.")
  message("  - What is le langage in a precise mathematical sense?")
  message("  - Find an analogy between le langue and cell systems.")
  message("  - Investigate a catastrophe structure in linguistics.")
}
