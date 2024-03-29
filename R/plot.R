#-----------------------------------------------------------------------------80
# 
#-----------------------------------------------------------------------------80
#' Visualize a three-dimensional data with labels and colors.
#'
#' This function visualizes a three-dimensional data with labels and colors.
#'
#' @param dataframe3D A dataframe with three columns.
#' @param labels NULL or a vector of labels of all the samples,
#'   corresponding to colors.
#' @param colors NULL or a vector of colors of all the samples,
#'   corresponding to labels.
#' @param theta Angle of the plot.
#' @param phi Angle of the plot.
#' @param title Title.
#' @param xlabel x-axis label.
#' @param ylabel y-axis label.
#' @param zlabel z-axis label.
#'
#' @return A scatter3D object in plot3D package.
#' @import plot3D
#' @export
#'
#' @examples
#' data(pbmcs_eg)
#' mat <- SingleCellExperiment::reducedDim(pbmcs_eg$CM, "UMAP")[, 1:3]
#' dataframe3D <- as.data.frame(mat)
#' labels <- SummarizedExperiment::colData(pbmcs_eg$CM)$seurat_clusters
#' plot_dataframe3D(dataframe3D = dataframe3D, labels = labels, colors = NULL,
#'                  theta = 45, phi = 20, title = "PBMC (CO & MSigDB)",
#'                  xlabel = "UMAP_1", ylabel = "UMAP_2", zlabel = "UMAP_3")
#'
plot_dataframe3D <- function(
  dataframe3D = NULL, labels = NULL, colors = NULL, theta = 30, phi = 30,
  title = "", xlabel = "", ylabel = "", zlabel = ""
){
  if(is.null(labels)){
    scatter3D(dataframe3D[, 1], dataframe3D[, 2], dataframe3D[, 3],
              main = title, xlab = xlabel, ylab = ylabel, zlab = zlabel,
              box = TRUE, bty = "b2", axes = TRUE, nticks = 5,
              theta = theta, phi = phi, pch = 16, cex = 0.5, alpha = 1.0,
              col = "black", colvar = NA, colkey = FALSE)
  }else{
    if(is.null(colors)){
      myggcolor <- function(n, l = 65){
        hues <- seq(15, 375, length = n + 1)
        grDevices::hcl(h = hues, l = l, c = 100)[seq_len(n)]
      }
      colors <- myggcolor(length(unique(labels))) ; colors <- colors[labels]
      dataframe3D$label <- labels ; dataframe3D$color <- colors
      dataframe3D <- dataframe3D[order(dataframe3D$label), ]
      scatter3D(dataframe3D[, 1], dataframe3D[, 2], dataframe3D[, 3],
                main = title, xlab = xlabel, ylab = ylabel, zlab = zlabel,
                box = TRUE, bty = "b2", axes = TRUE, nticks = 5,
                theta = theta, phi = phi, pch = 16, cex = 0.5, alpha = 1.0,
                col = dataframe3D$color, colvar = NA, colkey = FALSE)
      graphics::legend("bottomright", legend = unique(dataframe3D$label),
                       pch = 16, col = unique(dataframe3D$color), cex = 1.2,
                       inset = c(0.02))
    }else{
      dataframe3D$label <- labels ; dataframe3D$color <- colors
      dataframe3D <- dataframe3D[order(dataframe3D$label), ]
      scatter3D(dataframe3D[, 1], dataframe3D[, 2], dataframe3D[, 3],
                main = title, xlab = xlabel, ylab = ylabel, zlab = zlabel,
                box = TRUE, bty = "b2", axes = TRUE, nticks = 5,
                theta = theta, phi = phi, pch = 16, cex = 0.5, alpha = 1.0,
                col = dataframe3D$color, colvar = NA, colkey = FALSE)
      graphics::legend("bottomright", legend=unique(dataframe3D$label),
                       pch = 16, col = unique(dataframe3D$color), cex = 1.2,
                       inset = c(0.02))
    }
  }
}
#-----------------------------------------------------------------------------80
# 
#-----------------------------------------------------------------------------80
#' Visualize multivariate data by heatmaps.
#'
#' This function visualizes multivariate data by heatmaps.
#'
#' @param ssm_list A list of sign-by-sample matrices.
#' @param gem_list A list of gene-by-sample matrices.
#' @param ssmlabel_list NULL or a list of dataframes of sample (cell)
#'   labels and colors.
#'   The length of the list must be as same as that of ssm_list, and
#'   the order of labels in each list must be as same as those in ssm_list.
#' @param gemlabel_list NULL or a list of dataframes of sample (cell)
#'   annotations and colors.
#'   The length of the list must be as same as that of gem_list, and
#'   the order of labels in each list must be as same as those in gem_list.
#' @param nrand_samples Number of samples (cells) used for random sampling.
#' @param show_row_names TRUE or FALSE: if TRUE, row names are shown.
#' @param title Title.
#'
#' @return A ComplexHeatmap object.
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#' @export
#'
#' @examples
#' data(pbmcs_eg)
#' mat_CM <- SummarizedExperiment::assay(pbmcs_eg$CM, "counts")
#' mat_GO <- SummarizedExperiment::assay(pbmcs_eg$GO, "counts")
#' mat_KG <- SummarizedExperiment::assay(pbmcs_eg$KG, "counts")
#' ssm_list <- list(SSM_COMSig = mat_CM, SSM_GO = mat_GO, SSM_KEGG = mat_KG)
#' se <- SingleCellExperiment::altExp(pbmcs_eg$CM, "logcounts")
#' mat <- SummarizedExperiment::assay(se, "counts")
#' se <- SingleCellExperiment::altExp(pbmcs_eg$CM, "logcounts")
#' gem_list <- list(GeneExpr = SummarizedExperiment::assay(se, "counts"))
#' labels <- list() ; ssmlabel_list <- list()
#' for(i in seq_along(pbmcs_eg)){
#'   fa <- SummarizedExperiment::colData(pbmcs_eg[[i]])$seurat_clusters
#'   labels[[i]] <- data.frame(label = fa)
#'   colors <- rainbow(length(unique(labels[[i]]$label)))[labels[[i]]$label]
#'   labels[[i]]$color <- colors
#'   ssmlabel_list[[i]] <- labels[[i]]
#' }
#' names(ssmlabel_list) <- c("Label_COMSig", "Label_GO", "Label_KEGG")
#' phases <- SummarizedExperiment::colData(pbmcs_eg$CM)$Phase
#' label_CC <- data.frame(label = phases, color = NA)
#' gemlabel_list <- list(CellCycle = label_CC)
#' plot_multiheatmaps(ssm_list = ssm_list, gem_list = gem_list,
#'                    ssmlabel_list = ssmlabel_list,
#'                    gemlabel_list = gemlabel_list, nrand_samples = 50,
#'                    show_row_names = FALSE, title = "PBMC")
#'
plot_multiheatmaps <- function(
  ssm_list = NULL, gem_list = NULL, ssmlabel_list = NULL, gemlabel_list = NULL,
  nrand_samples = NULL, show_row_names = FALSE, title = NULL
){
  #--------------------------------------------------
  # Error handling
  #--------------------------------------------------
  if(dim(ssm_list[[1]])[2] > 2000){
    warning("Huge sample size> 2000. It is recommended to use nrand_samples.")
  }
  #--------------------------------------------------
  # Set names of dataframes.
  #--------------------------------------------------
  if(!is.null(ssmlabel_list)){
    for(i in seq_along(ssmlabel_list)){
      colnames(ssmlabel_list[[i]]) <- c("label", "color")
    }
  }
  if(!is.null(gemlabel_list)){
    for(i in seq_along(gemlabel_list)){
      colnames(gemlabel_list[[i]]) <- c("label", "color")
    }
  }
  #--------------------------------------------------
  # Fix ComplexHeatmap parameters.
  #--------------------------------------------------
  if(show_row_names){
    nrow_max <- 1
    for(i in seq_along(ssm_list)){
      if(nrow_max < dim(ssm_list[[i]])[1]){
        nrow_max <- dim(ssm_list[[i]])[1]
      }
    }
    for(i in seq_along(gem_list)){
      if(nrow_max < dim(gem_list[[i]])[1]){
        nrow_max <- dim(gem_list[[i]])[1]
      }
    }
    if(nrow_max > 200){
      message("Row names are removed because the number of rows> 200.")
      show_row_names <- FALSE
      row_names_side <- NULL
    }else{
      row_names_side <- "left"
    }
  }
  show_row_dend <- FALSE
  cluster_row_slices <- FALSE
  show_column_names <- FALSE
  column_dend_side <- "top"
  cluster_column_slices <- FALSE
  border <- FALSE
  #--------------------------------------------------
  # Random sampling
  #--------------------------------------------------
  if(!is.null(nrand_samples)){
    inds <- sample(ncol(ssm_list[[1]]), size = nrand_samples, replace = FALSE)
    for(i in seq_along(ssm_list)){
      ssm_list[[i]] <- ssm_list[[i]][, inds]
    }
    for(i in seq_along(ssmlabel_list)){
      ssmlabel_list[[i]] <- ssmlabel_list[[i]][inds, ]
    }
    for(i in seq_along(gem_list)){
      gem_list[[i]] <- gem_list[[i]][, inds]
    }
    for(i in seq_along(gemlabel_list)){
      gemlabel_list[[i]] <- gemlabel_list[[i]][inds, ]
    }
  }
  #--------------------------------------------------
  # Compute heatmaps of sign-by-sample matrices.
  #--------------------------------------------------
  p <- c()
  ha <- list()
  for(i in seq_along(ssm_list)){
    if((!is.element(NA, ssmlabel_list[[i]]$label)) &
       (!is.null(ssmlabel_list[[i]]$label))){
      inds <- order(ssmlabel_list[[i]]$label)
      if((!is.element(NA, ssmlabel_list[[i]]$color)) &
         (!is.null(ssmlabel_list[[i]]$color))){
        mycolor <- unique(ssmlabel_list[[i]][inds, ]$color)
      }else{
        myggcolor <- function(n, l = 65){
          hues <- seq(15, 375, length = n + 1)
          grDevices::hcl(h = hues, l = l, c = 100)[seq_len(n)]
        }
        mycolor <- myggcolor(length(unique(ssmlabel_list[[i]]$label)))
      }
      names(mycolor) <- unique(ssmlabel_list[[i]][inds, ]$label)
      mycolor <- list(mycolor)
      names(mycolor) <- names(ssmlabel_list)[[i]]
      mylabel <- list(ssmlabel_list[[i]]$label)
      names(mylabel) <- names(ssmlabel_list)[[i]]
      ha[[i]] <- HeatmapAnnotation(df = mylabel, col = mycolor,
                                   annotation_name_side = "left")

      if(show_row_names){
        q <- Heatmap(as.matrix(ssm_list[[i]]), top_annotation = ha[[i]],
                     # Row data
                     row_names_side = row_names_side, show_row_dend = FALSE,
                     cluster_row_slices = cluster_row_slices,
                     row_gap = unit(1.0, "mm"),
                     # Column data
                     show_column_names = show_column_names,
                     column_dend_side = column_dend_side,
                     column_split = ssmlabel_list[[i]]$label,
                     column_gap = unit(1.5, "mm"),
                     cluster_column_slices = cluster_column_slices,
                     column_title = title,
                     # Option
                     name = names(ssm_list)[[i]], border = border)
      }else{
        q <- Heatmap(as.matrix(ssm_list[[i]]), top_annotation = ha[[i]],
                     # Row data
                     show_row_names = show_row_names,
                     show_row_dend = show_row_dend,
                     cluster_row_slices = cluster_row_slices,
                     row_gap = unit(1.0, "mm"),
                     # Column data
                     show_column_names = show_column_names,
                     column_dend_side = column_dend_side,
                     column_split = ssmlabel_list[[i]]$label,
                     column_gap = unit(1.5, "mm"),
                     cluster_column_slices = cluster_column_slices,
                     column_title = title,
                     # Option
                     name = names(ssm_list)[[i]], border = border)
      }
    }else{
      if(show_row_names){
        q <- Heatmap(as.matrix(ssm_list[[i]]), top_annotation = NULL,
                     # Row data
                     row_names_side = row_names_side, show_row_dend = FALSE,
                     cluster_row_slices = cluster_row_slices,
                     row_gap = unit(1.0, "mm"),
                     # Column data
                     show_column_names = show_column_names,
                     column_dend_side = column_dend_side,
                     column_split = ssmlabel_list[[i]]$label,
                     column_gap = unit(1.5, "mm"),
                     cluster_column_slices = cluster_column_slices,
                     column_title = title,
                     # Option
                     name = names(ssm_list)[[i]], border = border)
      }else{
        q <- Heatmap(as.matrix(ssm_list[[i]]), top_annotation = NULL,
                     # Row data
                     show_row_names = show_row_names,
                     show_row_dend = show_row_dend,
                     cluster_row_slices = cluster_row_slices,
                     row_gap = unit(1.0, "mm"),
                     # Column data
                     show_column_names = show_column_names,
                     column_dend_side = column_dend_side,
                     column_split = ssmlabel_list[[i]]$label,
                     column_gap = unit(1.5, "mm"),
                     cluster_column_slices = cluster_column_slices,
                     column_title = title,
                     # Option
                     name = names(ssm_list)[[i]], border = border)
      }
    }
    p <- p %v% q
  }
  #--------------------------------------------------
  # Compute heatmaps of gene-by-sample matrices.
  #--------------------------------------------------
  ha2 <- list()
  for(i in seq_along(gem_list)){
    if((!is.element(NA, gemlabel_list[[i]]$label)) &
       (!is.null(gemlabel_list[[i]]$label))){
      inds <- order(gemlabel_list[[i]]$label)
      if((!is.element(NA, gemlabel_list[[i]]$color)) &
         (!is.null(gemlabel_list[[i]]$color))){
        mycolor <- unique(gemlabel_list[[i]][inds, ]$color)
      }else{
        myggcolor <- function(n, l = 65){
          hues <- seq(15, 375, length = n + 1)
          grDevices::hcl(h = hues, l = l, c = 100)[seq_len(n)]
        }
        mycolor <- myggcolor(length(unique(gemlabel_list[[i]]$label)))
      }
      names(mycolor) <- unique(gemlabel_list[[i]][inds, ]$label)
      mycolor <- list(mycolor)
      names(mycolor) <- names(gemlabel_list)[[i]]
      mylabel <- list(gemlabel_list[[i]]$label)
      names(mylabel) <- names(gemlabel_list)[[i]]
      ha2[[i]] <- HeatmapAnnotation(df = mylabel, col = mycolor,
                                    annotation_name_side = "left")

      if(show_row_names){
        q <- Heatmap(as.matrix(gem_list[[i]]), top_annotation = ha2[[i]],
                     # Row data
                     row_names_side = row_names_side, show_row_dend = FALSE,
                     cluster_row_slices = cluster_row_slices,
                     row_gap = unit(1.0, "mm"),
                     # Column data
                     show_column_names = show_column_names,
                     column_dend_side = column_dend_side,
                     column_split = ssmlabel_list[[i]]$label,
                     column_gap = unit(1.5, "mm"),
                     cluster_column_slices = cluster_column_slices,
                     column_title = title,
                     # Option
                     name = names(gem_list)[[i]], border = border)
      }else{
        q <- Heatmap(as.matrix(gem_list[[i]]), top_annotation = ha2[[i]],
                     # Row data
                     show_row_names = show_row_names,
                     show_row_dend = show_row_dend,
                     cluster_row_slices = cluster_row_slices,
                     row_gap = unit(1.0, "mm"),
                     # Column data
                     show_column_names = show_column_names,
                     column_dend_side = column_dend_side,
                     column_split = ssmlabel_list[[i]]$label,
                     column_gap = unit(1.5, "mm"),
                     cluster_column_slices = cluster_column_slices,
                     column_title = title,
                     # Option
                     name = names(gem_list)[[i]], border = border)
      }
    }else{
      if(show_row_names){
        q <- Heatmap(as.matrix(gem_list[[i]]), top_annotation = NULL,
                     # Row data
                     row_names_side = row_names_side, show_row_dend = FALSE,
                     cluster_row_slices = cluster_row_slices,
                     row_gap = unit(1.0, "mm"),
                     # Column data
                     show_column_names = show_column_names,
                     column_dend_side = column_dend_side,
                     column_split = ssmlabel_list[[i]]$label,
                     column_gap = unit(1.5, "mm"),
                     cluster_column_slices = cluster_column_slices,
                     column_title = title,
                     # Option
                     name = names(gem_list)[[i]], border = border)
      }else{
        q <- Heatmap(as.matrix(gem_list[[i]]), top_annotation = NULL,
                     # Row data
                     show_row_names = show_row_names,
                     show_row_dend = show_row_dend,
                     cluster_row_slices = cluster_row_slices,
                     row_gap = unit(1.0, "mm"),
                     # Column data
                     show_column_names = show_column_names,
                     column_dend_side = column_dend_side,
                     column_split = ssmlabel_list[[i]]$label,
                     column_gap = unit(1.5, "mm"),
                     cluster_column_slices = cluster_column_slices,
                     column_title = title,
                     # Option
                     name = names(gem_list)[[i]], border = border)
      }
    }
    p <- p %v% q
  }

  return(p)
}

