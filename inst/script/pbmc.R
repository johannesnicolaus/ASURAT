#-----------------------------------------------------------------------------80
# 
#-----------------------------------------------------------------------------80
#' Produce a read count table.
#'
#' This function produces a read count table from 10x Cell Ranger outputs.
#'
#' @param Directory path.
#'
#' @return 
#'
read_matrix_10xdata_v2 <- function(path_dir){
  mat <- Seurat::Read10X(data.dir = path_dir, gene.column = 2,
                         unique.features = TRUE, strip.suffix = FALSE)
  return(mat)
}
#-----------------------------------------------------------------------------80
# 
#-----------------------------------------------------------------------------80
#' Create a small read count table by processing conservative quality control.
#'
path_dir <- "inst/extdata/pbmc_4000/filtered_gene_bc_matrices/GRCh38/"
mat <- read_matrix_10xdata_v2(path_dir = path_dir)

nr <- apply(mat, 2, sum) ; tmp <- mat[, which((4010 < nr) & (nr < 20000))]
nc <- apply(tmp, 1, sum) ; tmp <- tmp[which((260 < nc) & (nc < 200000)), ]

pbmc_counts <- as.matrix(tmp)
save(pbmc_counts, file = "data/pbmc_counts.rda")
