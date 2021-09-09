#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_metaData <- function(df, title, title_size, xlabel, ylabel){
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2]), color="black", alpha=1, size=1) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"))
  p <- ggMarginal(p, type="histogram", margins="both", size=5, col="black",
                  fill="gray")
  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
find_inflection <- function(log_appr){
  I <- length(log_appr$x) - 1
  tmp <- c()
  for(i in 1:I){
    d <- abs((log_appr$y[i+1] - log_appr$y[i]) / (log_appr$x[i+1] - log_appr$x[i]))
    tmp <- rbind(tmp, c(log_appr$x[i], log_appr$y[i], d))
  }
  tmp <- as.data.frame(tmp)
  colnames(tmp) <- c("log_x", "log_y", "derivative")
  infp <- as.numeric(tmp[which(tmp$derivative == max(tmp$derivative))[1],][1:2])
  infp <- 10^(infp)
  return(infp)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_metaData_wAnnotation <- function(
  vec, target_xRange, n, title, title_size, xlabel, ylabel
){
  #--------------------------------------------------
  # Prepare a data frame to be output
  #--------------------------------------------------
  df <- as.numeric(sort(vec, decreasing = TRUE))
  df <- data.frame(x = 1:length(df), y = df)
  #--------------------------------------------------
  # Error criteria
  #--------------------------------------------------
  if(target_xRange[2] > max(which(df$y > 0)))
    stop("Error: `target_xRange` is too wide. In addition, it cannot include zero points.")
  #--------------------------------------------------
  # Find an inflection point within the `target_xRange`.
  #--------------------------------------------------
  df_target_xRange <- df[target_xRange[1]:target_xRange[2],]
  log_appr <- approx(log10(df_target_xRange$x), log10(df_target_xRange$y), n = n)
  annoP <- find_inflection(log_appr)
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2]), color="black", alpha=1, size=2) +
    geom_line(aes(x=10^(log_appr$x), y=10^(log_appr$y)), color="red", alpha=1, size=1) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=16, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="none") +
    scale_x_log10(limits=c(-NA, NA)) + scale_y_log10(limits=c(-NA, NA)) +
    geom_point(aes(x=annoP[1], y=annoP[2]), color="red", alpha=1, size=4) +
    geom_text_repel(
      aes(x=annoP[1], y=annoP[2]),
      label=paste("(", trunc(annoP[1]), ", ", trunc(annoP[2]), ")", sep = ""),
      size=5,
      colour="red",
      box.padding=unit(0.5, "lines"),
      point.padding=unit(0.5, "lines"),
      force=10,
      segment.color="red"
    )
  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_meanReads_wAnnotation <- function(
  obj, threshold, digits, title, title_size, xlabel, ylabel
){
  #--------------------------------------------------
  # Prepare a data frame to be output
  #--------------------------------------------------
  mat <- as.matrix(obj[["data"]][["raw"]])
  df <- apply(mat, 1, mean)
  df <- as.numeric(sort(df, decreasing = TRUE))
  df <- data.frame(x = 1:length(df), y = df)
  #--------------------------------------------------
  # Set the annotation point
  #--------------------------------------------------
  ind <- max(which(abs(df$y - threshold) == min(abs(df$y - threshold))))
  annoP <- c(df[ind,]$x, df[ind,]$y)
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2]), color="black", alpha=1, size=2) +
    geom_point(aes(x=annoP[1], y=annoP[2]), color="blue", alpha=1, size=4) +
    geom_hline(aes(yintercept=threshold), color="blue", linetype="dashed", size=0.75) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=16, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(0.5,2,0.5,0.5), "lines"), legend.position="none") +
    scale_x_log10(limits=c(-NA, NA)) + scale_y_log10(limits=c(-NA, NA)) +
    geom_point(aes(x=annoP[1], y=annoP[2]), color="blue", alpha=1, size=4) +
    geom_text_repel(
      aes(x=annoP[1], y=annoP[2]),
      label=paste("(", trunc(annoP[1]), ", ", round(annoP[2], digits=digits), ")", sep = ""),
      size=5,
      colour="blue",
      box.padding=unit(0.5, "lines"),
      point.padding=unit(0.5, "lines"),
      force=10,
      segment.color="blue"
    )
  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_Heatmap_GenexSamp <- function(obj, genes, method, show_nReads, title, name){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- as.matrix(obj[["data"]][["normalized"]])
  if(length(genes) == 1){
    mat <- t(as.matrix(mat[which(rownames(mat) %in% genes),]))
    rownames(mat) <- genes
  }else{
    mat <- as.matrix(mat[which(rownames(mat) %in% genes),])
  }
  tmp <- c()
  goi <- c()
  for(i in 1:length(genes)){
    if(genes[i] %in% rownames(mat)){
      tmp <- rbind(tmp, mat[which(rownames(mat) == genes[i]),])
      goi <- rbind(goi, genes[i])
    }
  }
  #--------------------------------------------------
  # Exception handling
  #--------------------------------------------------
  if(length(goi) == 0){
    stop("None of the input genes are included in your data.")
  }
  #--------------------------------------------------
  rownames(tmp) <- goi
  mat <- as.matrix(tmp)
  set.seed(8)
  col_hc <- hclust(dist(t(mat)), method = method)
  #--------------------------------------------------
  # Heatmap
  # Note: the positions of dendrogram leaves are
  # slightly adjusted by the gaps between slices.
  #--------------------------------------------------
  ht_opt$message = FALSE
  p <- Heatmap(
    mat,
    column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
    name=name,
    cluster_columns=col_hc, cluster_rows = FALSE,
    show_row_names=TRUE, row_names_side="right", show_row_dend=FALSE,
    show_column_names=FALSE, column_dend_side="top",
    show_parent_dend_line=FALSE
  )
  if(show_nReads){
    mtx <- t(obj[["sample"]][["nReads"]])
    rownames(mtx) <- "nReads"
    p_add <- Heatmap(
      mtx,
      name="nReads",
      cluster_columns=col_hc,
      show_row_names=TRUE, row_names_side="right", show_row_dend=FALSE,
      show_column_names=FALSE, show_column_dend=FALSE,
      col=colorRamp2(c(min(mtx), max(mtx)), c("cyan", "magenta"))
    )
    p <- p %v% p_add
  }

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_Heatmap_SignxSamp <- function(
  obj, data_type, category, algo_name, method, show_nReads, title, name,
  show_rownames_sign, show_rownames_label, show_rownames_nReads,
  default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  sign <- "sign"
  slot_name <- paste(data_type, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[slot_name]][[category]])
  set.seed(8)
  row_hc <- hclust(dist(mat), method = method)
  col_hc <- hclust(dist(t(mat)), method = method)
  if(is.null(algo_name)){
    labels <- NULL
  }else{
    slot_name <- paste(algo_name, "_", data_type, "_", category, sep = "")
    labels <- obj[["sample"]][[slot_name]]
  }
  #--------------------------------------------------
  # Heatmap
  # Note: the positions of dendrogram leaves are
  # slightly adjusted by the gaps between slices.
  #--------------------------------------------------
  ht_opt$message = FALSE
  if(is.null(algo_name)){
    p <- Heatmap(
      mat,
      column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
      name=name,
      cluster_rows=row_hc,
      cluster_columns=col_hc,
      show_row_names=show_rownames_sign, row_names_side="right", row_dend_side="left",
      show_column_names=FALSE, column_dend_side="top",
      show_parent_dend_line=FALSE
    )
  }else{
    if(!is.null(labels)){
      if(default_color){
        #------------------------------
        # Option 1: ggplot's default
        #------------------------------
        tmp <- unique(sort(labels))
        ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
          hues <- seq(15, 375, length=n+1)
          hcl(h=hues, l=l, c=100)[1:n]
        }
        my_colors <- ggColorHue(length(tmp))
        names(my_colors) <- tmp
      }else{
        #------------------------------
        # Option 2: rainbow colors
        #------------------------------
        tmp <- unique(sort(labels))
        my_colors <- rainbow(length(tmp))
        names(my_colors) <- tmp
      }
      ha <- HeatmapAnnotation(
        Label = labels, col = list(Label = my_colors), annotation_name_side = "left"
      )
      column_split <- labels
      p <- Heatmap(
        mat,
        column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
        name=name,
        cluster_rows=row_hc,
        column_split=column_split, column_gap=unit(1.5, "mm"),
        border=FALSE,
        show_row_names=show_rownames_sign, row_names_side="right", row_dend_side="left",
        show_column_names=FALSE, column_dend_side="top",
        show_parent_dend_line=FALSE, top_annotation=ha
      )
    }else{
      column_split <- labels
      p <- Heatmap(
        mat,
        column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
        name=name,
        cluster_rows=row_hc,
        column_split=column_split, column_gap=unit(1.5, "mm"),
        border=FALSE,
        show_row_names=show_rownames_sign, row_names_side="right", row_dend_side="left",
        show_column_names=FALSE, column_dend_side="top",
        show_parent_dend_line=FALSE
      )
    }
  }
  if(show_nReads){
    mtx <- t(obj[["sample"]][["nReads"]])
    rownames(mtx) <- "nReads"
    q <- Heatmap(
      mtx,
      name="nReads",
      cluster_columns=col_hc,
      show_row_names=show_rownames_nReads, row_names_side="right", show_row_dend=FALSE,
      show_column_names=FALSE, show_column_dend=FALSE,
      col=colorRamp2(c(min(mtx), max(mtx)), c("cyan", "magenta"))
    )
    p <- p %v% q
  }

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_tsne_sign <- function(
  obj, data_type, category, algo_name, theta = 40, phi = 40, title, title_size,
  xlabel, ylabel, zlabel, default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["reduction"]][["tsne"]][[data_type]][[category]][["Y"]]
  df <- as.data.frame(mat)
  if(is.null(algo_name) == FALSE){
    slot_name <- paste(algo_name, "_", data_type, "_", category, sep = "")
    df$label <- obj[["sample"]][[slot_name]]
  }
  dim_tsne <- dim(mat)[2]
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(dim_tsne == 2){
    if(is.null(df$label)){
      p <- ggplot() +
        geom_point(aes(x=df[,1], y=df[,2]), color="black", size=0.5, alpha=1.0) +
        labs(title=title, x=xlabel, y=ylabel) +
        theme_classic(base_size=20, base_family="Helvetica") +
        theme(plot.title=element_text(hjust=0.5, size=title_size),
          plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    }else{
      if(default_color){
        #------------------------------
        # Option 1: ggplot's default
        #------------------------------
        n_groups <- length(unique(df$label))
        ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
          hues <- seq(15, 375, length=n+1)
          hcl(h=hues, l=l, c=100)[1:n]
        }
        my_colors <- ggColorHue(n_groups)
      }else{
        #------------------------------
        # Option 2: rainbow colors
        #------------------------------
        n_groups <- length(unique(df$label))
        my_colors <- rainbow(n_groups)
      }
      #------------------------------
      # Positions of labels
      #------------------------------
      cluster_num <- sort(unique(df$label))
      cog <- c()  # Center of gravity
      for(i in cluster_num){
        cog <- rbind(cog, cbind(
          mean(df[which(df$label == i),][,1]), mean(df[which(df$label == i),][,2])))
      }
      cog <- data.frame(x = cog[,1], y = cog[,2], label = cluster_num)
      #------------------------------
      # ggplot
      #------------------------------
      p <- ggplot() +
        geom_point(aes(x=df[,1], y=df[,2], color=as.factor(df$label)), size=0.5, alpha=1.0) +
        labs(title=title, x=xlabel, y=ylabel, colour="") +
        scale_colour_manual(values=my_colors) +
        theme_classic(base_size=20, base_family="Helvetica") +
        theme(plot.title=element_text(hjust=0.5, size=title_size),
          plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
        geom_text(aes(x=cog[,1], y=cog[,2], label=sprintf("%d", cluster_num)), size=6) +
        guides(colour = guide_legend(override.aes = list(size=4)))
    }
    return(p)
  }else if(dim_tsne == 3){
    if(is.null(df$label)){
      #------------------------------
      # scatter3D in plot3d package
      #------------------------------
      scatter3D(df[,1], df[,2], df[,3],
        main=title, xlab=xlabel, ylab=ylabel, zlab=zlabel, cex.main=title_size,
        box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
        theta=theta, phi=phi,
        pch=16, cex=0.5, alpha=1.0, col="black", colvar=NA, colkey=FALSE)
    }else{
      if(default_color){
        #------------------------------
        # Option 1: ggplot's default
        #------------------------------
        n_groups <- length(unique(df$label))
        ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
          hues <- seq(15, 375, length=n+1)
          hcl(h=hues, l=l, c=100)[1:n]
        }
        my_colors <- ggColorHue(n_groups)[df$label]
        label_colors <- ggColorHue(n_groups)
      }else{
        #------------------------------
        # Option 2: rainbow colors
        #------------------------------
        n_groups <- length(unique(df$label))
        my_colors <- rainbow(n_groups)[df$label]
        label_colors <- rainbow(n_groups)
      }
      #------------------------------
      # scatter3D in plot3d package
      #------------------------------
      par(mar = c(1.5,1.5,1.5,1.5))
      par(oma = c(1.5,1.5,1.5,1.5))
      scatter3D(df[,1], df[,2], df[,3],
        main=title, xlab=xlabel, ylab=ylabel, zlab=zlabel, cex.main=title_size,
        box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
        theta=theta, phi=phi,
        pch=16, cex=0.5, alpha=1.0, col=my_colors, colvar=NA, colkey=FALSE)
      legend("bottomright", legend=1:n_groups, pch=16, col=label_colors, cex=1, inset=c(0.02))
    }
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_umap_sign <- function(
  obj, data_type, category, algo_name, theta = 40, phi = 40, title, title_size,
  xlabel, ylabel, zlabel, default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["reduction"]][["umap"]][[data_type]][[category]][["layout"]]
  df <- as.data.frame(mat)
  if(is.null(algo_name) == FALSE){
    slot_name <- paste(algo_name, "_", data_type, "_", category, sep = "")
    df$label <- obj[["sample"]][[slot_name]]
  }
  dim_umap <- dim(mat)[2]
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(dim_umap == 2){
    if(is.null(df$label)){
      p <- ggplot() +
        geom_point(aes(x=df[,1], y=df[,2]), color="black", size=0.5, alpha=1.0) +
        labs(title=title, x=xlabel, y=ylabel) +
        theme_classic(base_size=20, base_family="Helvetica") +
        theme(plot.title=element_text(hjust=0.5, size=title_size),
          plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    }else{
      if(default_color){
        #------------------------------
        # Option 1: ggplot's default
        #------------------------------
        n_groups <- length(unique(df$label))
        ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
          hues <- seq(15, 375, length=n+1)
          hcl(h=hues, l=l, c=100)[1:n]
        }
        my_colors <- ggColorHue(n_groups)
      }else{
        #------------------------------
        # Option 2: rainbow colors
        #------------------------------
        n_groups <- length(unique(df$label))
        my_colors <- rainbow(n_groups)
      }
      #------------------------------
      # Positions of labels
      #------------------------------
      cluster_num <- sort(unique(df$label))
      cog <- c()  # Center of gravity
      for(i in cluster_num){
        cog <- rbind(cog, cbind(
          mean(df[which(df$label == i),][,1]), mean(df[which(df$label == i),][,2])))
      }
      cog <- data.frame(x = cog[,1], y = cog[,2], label = cluster_num)
      #------------------------------
      # ggplot
      #------------------------------
      p <- ggplot() +
        geom_point(aes(x=df[,1], y=df[,2], color=as.factor(df$label)), size=0.5, alpha=1.0) +
        labs(title=title, x=xlabel, y=ylabel, colour="") +
        scale_colour_manual(values=my_colors) +
        theme_classic(base_size=20, base_family="Helvetica") +
        theme(plot.title=element_text(hjust=0.5, size=title_size),
          plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
#        geom_text(aes(x=cog[,1], y=cog[,2], label=sprintf("%d", cluster_num)), size=6) +
        guides(colour = guide_legend(override.aes = list(size=4)))
    }
    return(p)
  }else if(dim_umap == 3){
    if(is.null(df$label)){
      #------------------------------
      # scatter3D in plot3d package
      #------------------------------
      scatter3D(df[,1], df[,2], df[,3],
        main=title, xlab=xlabel, ylab=ylabel, zlab=zlabel, cex.main=title_size,
        box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
        theta=theta, phi=phi,
        pch=16, cex=0.5, alpha=1.0, col="black", colvar=NA, colkey=FALSE)
    }else{
      if(default_color){
        #------------------------------
        # Option 1: ggplot's default
        #------------------------------
        n_groups <- length(unique(df$label))
        ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
          hues <- seq(15, 375, length=n+1)
          hcl(h=hues, l=l, c=100)[1:n]
        }
        my_colors <- ggColorHue(n_groups)[df$label]
        label_colors <- ggColorHue(n_groups)
      }else{
        #------------------------------
        # Option 2: rainbow colors
        #------------------------------
        n_groups <- length(unique(df$label))
        my_colors <- rainbow(n_groups)[df$label]
        label_colors <- rainbow(n_groups)
      }
      #------------------------------
      # scatter3D in plot3d package
      #------------------------------
      par(mar = c(1.5,1.5,1.5,1.5))
      par(oma = c(1.5,1.5,1.5,1.5))
      scatter3D(df[,1], df[,2], df[,3],
        main=title, xlab=xlabel, ylab=ylabel, zlab=zlabel, cex.main=title_size,
        box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
        theta=theta, phi=phi,
        pch=16, cex=0.5, alpha=1.0, col=my_colors, colvar=NA, colkey=FALSE)
      legend("bottomright", legend=1:n_groups, pch=16, col=label_colors, cex=1, inset=c(0.02))
    }
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_dmap_sign <- function(
  obj, data_type, category, algo_name, dims, theta = 40, phi = 40,
  title, title_size, xlabel, ylabel, zlabel, default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["reduction"]][["dmap"]][[data_type]][[category]]@eigenvectors
  df <- as.data.frame(mat)
  if(is.null(algo_name) == FALSE){
    slot_name <- paste(algo_name, "_", data_type, "_", category, sep = "")
    df$label <- obj[["sample"]][[slot_name]]
  }
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(length(dims) == 2){
    if(is.null(df$label)){
      p <- ggplot() +
        geom_point(aes(x=df[,dims[1]], y=df[,dims[2]]), color="black", size=0.5, alpha=1.0) +
        labs(title=title, x=xlabel, y=ylabel) +
        theme_classic(base_size=20, base_family="Helvetica") +
        theme(plot.title=element_text(hjust=0.5, size=title_size),
          plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    }else{
      if(default_color){
        #------------------------------
        # Option 1: ggplot's default
        #------------------------------
        n_groups <- length(unique(df$label))
        ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
          hues <- seq(15, 375, length=n+1)
          hcl(h=hues, l=l, c=100)[1:n]
        }
        my_colors <- ggColorHue(n_groups)
      }else{
        #------------------------------
        # Option 2: rainbow colors
        #------------------------------
        n_groups <- length(unique(df$label))
        my_colors <- rainbow(n_groups)
      }
      #------------------------------
      # Positions of labels
      #------------------------------
      cluster_num <- sort(unique(df$label))
      cog <- c()  # Center of gravity
      for(i in cluster_num){
        cog <- rbind(cog, cbind(
          mean(df[which(df$label == i),][,1]), mean(df[which(df$label == i),][,2])))
      }
      cog <- data.frame(x = cog[,1], y = cog[,2], label = cluster_num)
      #------------------------------
      # ggplot
      #------------------------------
      p <- ggplot() +
        geom_point(aes(x=df[,dims[1]], y=df[,dims[2]], color=as.factor(df$label)), size=0.5, alpha=1.0) +
        labs(title=title, x=xlabel, y=ylabel, colour="") +
        scale_colour_manual(values=my_colors) +
        theme_classic(base_size=20, base_family="Helvetica") +
        theme(plot.title=element_text(hjust=0.5, size=title_size),
          plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
        geom_text(aes(x=cog[,1], y=cog[,2], label=sprintf("%d", cluster_num)), size=6) +
        guides(colour = guide_legend(override.aes = list(size=4)))
    }
    return(p)
  }else if(length(dims) == 3){
    if(is.null(df$label)){
      #------------------------------
      # scatter3D in plot3d package
      #------------------------------
      scatter3D(df[,dims[1]], df[,dims[2]], df[,dims[3]],
        main=title, xlab=xlabel, ylab=ylabel, zlab=zlabel,
        cex.main=title_size,
        box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
        theta=theta, phi=phi,
        pch=16, cex=0.5, alpha=1.0, col="black", colvar=NA, colkey=FALSE)
    }else{
      if(default_color){
        #------------------------------
        # Option 1: ggplot's default
        #------------------------------
        n_groups <- length(unique(df$label))
        ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
          hues <- seq(15, 375, length=n+1)
          hcl(h=hues, l=l, c=100)[1:n]
        }
        my_colors <- ggColorHue(n_groups)[df$label]
        label_colors <- ggColorHue(n_groups)
      }else{
        #------------------------------
        # Option 2: rainbow colors
        #------------------------------
        n_groups <- length(unique(df$label))
        my_colors <- rainbow(n_groups)[df$label]
        label_colors <- rainbow(n_groups)
      }
      #------------------------------
      # scatter3D in plot3d package
      #------------------------------
      par(mar = c(2.0,2.0,2.0,2.0))
      par(oma = c(2.0,2.0,2.0,2.0))
      scatter3D(df[,1], df[,2], df[,3],
        main=title, xlab=xlabel, ylab=ylabel, zlab=zlabel, cex.main=title_size,
        box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
        theta=theta, phi=phi,
        pch=16, cex=0.5, alpha=1.0, col=my_colors, colvar=NA, colkey=FALSE)
      legend("bottomright", legend=1:n_groups, pch=16, col=label_colors, cex=1.0, inset=c(0.02))
    }
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_ElasticTree_dmap <- function(
  obj, data_type, category, theta = 40, phi = 40, title, title_size,
  xlabel, ylabel, zlabel, default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  tmp <- obj[["clustering"]][["merlot"]][[data_type]][[category]][["ElasticTree"]]
  #--------------------------------------------------
  # Prepare data frames to be output
  #--------------------------------------------------
  df <- as.data.frame(tmp[["CellCoords"]])
  df$branch <- tmp[["Cells2Branches"]]
  dg <- as.data.frame(tmp[["Nodes"]])
  dh <- as.data.frame(tmp[["EndpointsCoords"]])
  di <- tmp[["BranchpointsCoords"]]
  if(is.null(dim(di))){
    di <- as.data.frame(t(tmp[["BranchpointsCoords"]]))
  }else{
    di <- as.data.frame(tmp[["BranchpointsCoords"]])
  }
  dj <- list()
  for(j in unique(sort(tmp[["Cells2Branches"]]))){
    dj[[j]] <- as.data.frame(dg[tmp[["Branches"]][[j]],])
  }
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(default_color){
    #------------------------------
    # Option 1: ggplot's default
    #------------------------------
    n_groups <- length(unique(df$branch))
    ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
      hues <- seq(15, 375, length=n+1)
      hcl(h=hues, l=l, c=100)[1:n]
    }
    my_colors <- ggColorHue(n_groups)
  }else{
    #------------------------------
    # Option 2: rainbow colors
    #------------------------------
    n_groups <- length(unique(df$branch))
    my_colors <- rainbow(n_groups)
  }
  if(dim(tmp[["CellCoords"]])[2] == 2){
    p <- ggplot() +
      geom_point(aes(x=df[,1], y=df[,2], color=as.factor(df$branch)), size=3.0, alpha=0.3)
    for(j in unique(sort(tmp[["Cells2Branches"]]))){
      p <- p + geom_line(
        aes_string(x=dj[[j]][,1], y=dj[[j]][,2]), color="black", size=1.0, alpha=1.0)
    }
    p <- p +
      geom_point(aes(x=dh[,1], y=dh[,2]), color="black", size=7.0, alpha=0.5) +
      geom_point(aes(x=di[,1], y=di[,2]), color="black", size=7.0, alpha=0.5) +
      geom_point(aes(x=dg[,1], y=dg[,2]), color="black", size=2.0, alpha=1.0) +
      labs(title=title, x=xlabel, y=ylabel, colour="") +
      scale_fill_manual(values=my_colors) +
      theme_classic(base_size=20, base_family="Helvetica") +
      theme(plot.title=element_text(hjust=0.5, size=title_size),
        plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    return(p)
  }else if(dim(tmp[["CellCoords"]])[2] == 3){
    #------------------------------
    # scatter3D in plot3d package
    #------------------------------
    scatter3D(df[,1], df[,2], df[,3],
      main=title, xlab="DC_1", ylab="DC_2", zlab="DC_3", cex.main=title_size,
      box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
      theta=theta, phi=phi,
      pch=16, cex=0.75, alpha=0.1, col=my_colors[df$branch], colvar=NA, colkey=FALSE)
    for(j in unique(sort(tmp[["Cells2Branches"]]))){
      scatter3D(dj[[j]][,1], dj[[j]][,2], dj[[j]][,3],
        type="l", lwd=2, alpha=1, col="black", add=TRUE)
    }
    scatter3D(dh[,1], dh[,2], dh[,3], pch=16, cex=3, alpha=0.25, col="black", add=TRUE)
    scatter3D(di[,1], di[,2], di[,3], pch=16, cex=3, alpha=0.40, col="black", add=TRUE)
    scatter3D(dg[,1], dg[,2], dg[,3], pch=16, cex=1, alpha=1, col="black", add=TRUE)
    legend("bottomright", legend=unique(sort(df$branch)), pch=16,
      col=my_colors[1:n_groups], cex=1.3, inset=c(0.02))
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_Pseudotime_dmap <- function(
  obj, data_type, category, theta = 40, phi = 40, title, title_size,
  xlabel, ylabel, zlabel
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  tmp <- obj[["clustering"]][["merlot"]][[data_type]][[category]]
  #--------------------------------------------------
  # Prepare data frames to be output
  #--------------------------------------------------
  df <- as.data.frame(tmp[["ElasticTree"]][["CellCoords"]])
  df$pseudotime <- tmp[["Pseudotimes"]][["Times_cells"]]
  df$color <- (df$pseudotime - min(df$pseudotime)) /
    (max(df$pseudotime) - min(df$pseudotime))
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(dim(tmp[["ElasticTree"]][["CellCoords"]])[2] == 2){
    p <- ggplot() +
      geom_point(aes(x=df[,1], y=df[,2], color=df$color), size=2, alpha=1) +
      scale_colour_gradientn(colours = plasma(50)) +
      labs(title=title, x=xlabel, y=ylabel, color="") +
      theme_classic(base_size=20, base_family="Helvetica") +
      theme(plot.title=element_text(hjust=0.5, size=title_size),
        plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    return(p)
  }else if(dim(tmp[["ElasticTree"]][["CellCoords"]])[2] == 3){
    #------------------------------
    # scatter3D in plot3d package
    #------------------------------
    legend.col <- function(col, lev){
      opar <- par
      n <- length(col)
      bx <- par("usr")
      box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
        bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 40)
      box.cy <- c(bx[3], bx[3])
      box.sy <- (bx[4] - bx[3]) / n
      xx <- rep(box.cx, each = 2)
      par(xpd = TRUE)
      for(i in 1:n){
        yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * (i)),
          box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
      }
      par(new = TRUE)
      plot(0, 0, type = "n", ylim = c(min(lev), max(lev)),
        yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
      axis(side = 4, las = 2, tick = FALSE, line = .25)
      par <- opar
    }
    par(mar = c(2.0,2.0,2.0,2.0))
    par(oma = c(2.0,2.0,2.0,2.0))
    scatter3D(df[,1], df[,2], df[,3],
      main=title, xlab=xlabel, ylab=ylabel, zlab=zlabel, cex.main=title_size,
      box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
      theta=theta, phi=phi,
      pch=16, cex=1.0, alpha=0.8, col=plasma(50)[as.numeric(cut(df$color, breaks=50))],
      colvar=NA, colkey=FALSE)
    legend.col(col=plasma(50), lev = df$color)
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_bargraph_sign <- function(
  obj, data_type, category, algo_name, title, title_size, xlabel, ylabel, ymax
){
  #--------------------------------------------------
  # Set data frame to be output
  #--------------------------------------------------
  slot_name <- paste(algo_name, "_", data_type, "_", category, sep = "")
  tmp <- obj[["sample"]]
  colnames(tmp)[which(colnames(tmp) == slot_name)] <- "Count"
  df <- tmp %>% group_by(Count) %>% summarise (n=n())
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot(df) +
   geom_bar(aes(x=as.factor(Count), y=n),
      color="black", fill="black", stat="identity", width=0.7) +
    theme_classic(base_size=18, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size)) +
    geom_text(aes(x=as.factor(Count), y=n, label=sprintf("%d", n)), size=6, vjust=-0.5) +
    ylim(0,ymax) + labs(title=title, x=xlabel, y=ylabel)

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_violin_signScore <- function(
  obj, sign_name, data_type_for_label, category_for_label, algo_name_for_label,
  data_type_for_expr, category_for_expr, title, title_size, default_color = TRUE
){
  #--------------------------------------------------
  # Set data frame to be output
  #--------------------------------------------------
  sign <- "sign"
  slot_name <- paste(data_type_for_expr, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[slot_name]][[category_for_expr]])
  slot_name <- paste(algo_name_for_label, "_", data_type_for_label, "_",
                     category_for_label, sep = "")
  df <- data.frame(label = obj[["sample"]][[slot_name]],
                   expr = mat[which(rownames(mat) == sign_name),])
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(default_color){
    #------------------------------
    # Option 1: ggplot's default
    #------------------------------
    n_groups <- length(unique(df$label))
    ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
      hues <- seq(15, 375, length=n+1)
      hcl(h=hues, l=l, c=100)[1:n]
    }
    my_colors <- ggColorHue(n_groups)
  }else{
    #------------------------------
    # Option 2: rainbow colors
    #------------------------------
    n_groups <- length(unique(df$label))
    my_colors <- rainbow(n_groups)
  }
  p <- ggplot() +
    geom_violin(aes(x=as.factor(df[,1]), y=df[,2], fill=as.factor(df[,1])),
      trim=FALSE, size=0.5) +
    geom_boxplot(aes(x=as.factor(df[,1]), y=df[,2]), width=0.15, alpha=0.6) +
    labs(title=title,
      x=paste("Label (", data_type_for_label, ": ", category_for_label, ")", sep = ""),
      y="Sign score", fill="") +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="none")
  p <- p + scale_fill_manual(values=my_colors)

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_tsne_geneExpression <- function(
  obj, gene_names, zscore = FALSE, data_type_for_tsne, category_for_tsne,
  theta, phi, title, title_size, label_name, xlabel, ylabel, zlabel
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["reduction"]][["tsne"]][[data_type_for_tsne]][[category_for_tsne]][["Y"]]
  df <- as.data.frame(mat)
  dim_tsne <- dim(mat)[2]
  mat <- as.matrix(obj[["data"]][["normalized"]])
  if(length(gene_names) == 1){
    df$expr <- mat[which(rownames(mat) == gene_names),]
  }else{ 
    df$expr <- apply(mat[which(rownames(mat) %in% gene_names),], 2, mean)
  }
  if(zscore){
    stdev <- sd(df$expr)
    df$expr <- (df$expr - mean(df$expr)) / stdev
  }
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(dim_tsne == 2){
    p <- ggplot() +
      geom_point(aes(x=df[,1], y=df[,2], color=df$expr), size=0.5, alpha=1.0) +
      labs(title=title, x=xlabel, y=ylabel, color=label_name) +
      scale_colour_gradientn(colours=c("blue", "white", "red")) +
      theme_classic(base_size=20, base_family="Helvetica") +
      theme(plot.title=element_text(hjust=0.5, size=title_size),
        plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    return(p)
  }else if(dim_tsne == 3){
    #------------------------------
    # scatter3D in plot3d package
    #------------------------------
    legend.col <- function(col, lev){
      opar <- par
      n <- length(col)
      bx <- par("usr")
      box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
        bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 40)
      box.cy <- c(bx[3], bx[3])
      box.sy <- (bx[4] - bx[3]) / n
      xx <- rep(box.cx, each = 2)
      par(xpd = TRUE)
      for(i in 1:n){
        yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * (i)),
          box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
      }
      par(new = TRUE)
      plot(0, 0, type = "n", ylim = c(min(lev), max(lev)),
        yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
      axis(side = 4, las = 2, tick = FALSE, line = .25)
      par <- opar
    }
    par(mar = c(1.5,1.5,2.5,1.5))
    par(oma = c(1.5,1.5,2.5,1.5))
    scatter3D(df[,1], df[,2], df[,3],
      main=paste(sign_name, "\n", title, sep=""), xlab=xlabel, ylab=ylabel, zlab=zlabel,
      cex.main=title_size,
      box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
      theta=theta, phi=phi,
      pch=16, cex=0.5, alpha=1.0,
      col=colorRampPalette(c("blue", "white", "red"))(50)[as.numeric(cut(df$expr, breaks=50))],
      colvar=NA, colkey=FALSE)
    legend.col(col=colorRampPalette(c("blue", "white", "red"))(50), lev = df$expr)
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_tsne_signScore <- function(
  obj, sign_name, data_type_for_tsne, category_for_tsne, data_type_for_expr,
  category_for_expr, theta, phi, title, title_size, label_name, xlabel, ylabel,
  zlabel
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["reduction"]][["tsne"]][[data_type_for_tsne]][[category_for_tsne]][["Y"]]
  df <- as.data.frame(mat)
  dim_tsne <- dim(mat)[2]
  sign <- "sign"
  slot_name <- paste(data_type_for_expr, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[slot_name]][[category_for_expr]])
  df$expr <- mat[which(rownames(mat) == sign_name),]
#  stdev <- sd(df$expr)
#  df$expr <- (df$expr - mean(df$expr)) / stdev
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(dim_tsne == 2){
    p <- ggplot() +
      geom_point(aes(x=df[,1], y=df[,2], color=df$expr), size=0.5, alpha=1.0) +
      labs(title=title, x=xlabel, y=ylabel, color=label_name) +
      scale_colour_gradientn(colours=c("blue", "white", "red")) +
      theme_classic(base_size=20, base_family="Helvetica") +
      theme(plot.title=element_text(hjust=0.5, size=title_size),
        plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    return(p)
  }else if(dim_tsne == 3){
    #------------------------------
    # scatter3D in plot3d package
    #------------------------------
    legend.col <- function(col, lev){
      opar <- par
      n <- length(col)
      bx <- par("usr")
      box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
        bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 40)
      box.cy <- c(bx[3], bx[3])
      box.sy <- (bx[4] - bx[3]) / n
      xx <- rep(box.cx, each = 2)
      par(xpd = TRUE)
      for(i in 1:n){
        yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * (i)),
          box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
      }
      par(new = TRUE)
      plot(0, 0, type = "n", ylim = c(min(lev), max(lev)),
        yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
      axis(side = 4, las = 2, tick = FALSE, line = .25)
      par <- opar
    }
    par(mar = c(1.5,1.5,2.5,1.5))
    par(oma = c(1.5,1.5,2.5,1.5))
    scatter3D(df[,1], df[,2], df[,3],
      main=paste(sign_name, "\n", title, sep=""), xlab=xlabel, ylab=ylabel, zlab=zlabel,
      cex.main=title_size,
      box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
      theta=theta, phi=phi,
      pch=16, cex=0.5, alpha=1.0,
      col=colorRampPalette(c("blue", "white", "red"))(50)[as.numeric(cut(df$expr, breaks=50))],
      colvar=NA, colkey=FALSE)
    legend.col(col=colorRampPalette(c("blue", "white", "red"))(50), lev = df$expr)
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_umap_geneExpression <- function(
  obj, gene_names, zscore = FALSE, data_type_for_umap, category_for_umap,
  theta, phi, title, title_size, label_name, xlabel, ylabel, zlabel
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["reduction"]][["umap"]][[data_type_for_umap]][[category_for_umap]][["layout"]]
  df <- as.data.frame(mat)
  dim_umap <- dim(mat)[2]
  mat <- as.matrix(obj[["data"]][["normalized"]])
  if(length(gene_names) == 1){
    df$expr <- mat[which(rownames(mat) == gene_names),]
  }else{ 
    df$expr <- apply(mat[which(rownames(mat) %in% gene_names),], 2, mean)
  }
  if(zscore){
    stdev <- sd(df$expr)
    df$expr <- (df$expr - mean(df$expr)) / stdev
  }
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(dim_umap == 2){
    p <- ggplot() +
      geom_point(aes(x=df[,1], y=df[,2], color=df$expr), size=0.5, alpha=1.0) +
      labs(title=title, x=xlabel, y=ylabel, color=label_name) +
      scale_colour_gradientn(colours=c("blue", "white", "red")) +
      theme_classic(base_size=20, base_family="Helvetica") +
      theme(plot.title=element_text(hjust=0.5, size=title_size),
        plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    return(p)
  }else if(dim_umap == 3){
    #------------------------------
    # scatter3D in plot3d package
    #------------------------------
    legend.col <- function(col, lev){
      opar <- par
      n <- length(col)
      bx <- par("usr")
      box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
        bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 40)
      box.cy <- c(bx[3], bx[3])
      box.sy <- (bx[4] - bx[3]) / n
      xx <- rep(box.cx, each = 2)
      par(xpd = TRUE)
      for(i in 1:n){
        yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * (i)),
          box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
      }
      par(new = TRUE)
      plot(0, 0, type = "n", ylim = c(min(lev), max(lev)),
        yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
      axis(side = 4, las = 2, tick = FALSE, line = .25)
      par <- opar
    }
    par(mar = c(1.5,1.5,2.5,1.5))
    par(oma = c(1.5,1.5,2.5,1.5))
    scatter3D(df[,1], df[,2], df[,3],
      main=paste(sign_name, "\n", title, sep=""), xlab=xlabel, ylab=ylabel, zlab=zlabel,
      cex.main=title_size,
      box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
      theta=theta, phi=phi,
      pch=16, cex=0.5, alpha=1.0,
      col=colorRampPalette(c("blue", "white", "red"))(50)[as.numeric(cut(df$expr, breaks=50))],
      colvar=NA, colkey=FALSE)
    legend.col(col=colorRampPalette(c("blue", "white", "red"))(50), lev = df$expr)
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_umap_signScore <- function(
  obj, sign_name, data_type_for_umap, category_for_umap, data_type_for_expr,
  category_for_expr, theta, phi, title, title_size, label_name, xlabel, ylabel,
  zlabel
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["reduction"]][["umap"]][[data_type_for_umap]][[category_for_umap]][["layout"]]
  df <- as.data.frame(mat)
  dim_umap <- dim(mat)[2]
  sign <- "sign"
  slot_name <- paste(data_type_for_expr, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[slot_name]][[category_for_expr]])
  df$expr <- mat[which(rownames(mat) == sign_name),]
#  stdev <- sd(df$expr)
#  df$expr <- (df$expr - mean(df$expr)) / stdev
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(dim_umap == 2){
    p <- ggplot() +
      geom_point(aes(x=df[,1], y=df[,2], color=df$expr), size=0.5, alpha=1.0) +
      labs(title=title, x=xlabel, y=ylabel, color=label_name) +
      scale_colour_gradientn(colours=c("blue", "white", "red")) +
      theme_classic(base_size=20, base_family="Helvetica") +
      theme(plot.title=element_text(hjust=0.5, size=title_size),
        plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    return(p)
  }else if(dim_umap == 3){
    #------------------------------
    # scatter3D in plot3d package
    #------------------------------
    legend.col <- function(col, lev){
      opar <- par
      n <- length(col)
      bx <- par("usr")
      box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
        bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 40)
      box.cy <- c(bx[3], bx[3])
      box.sy <- (bx[4] - bx[3]) / n
      xx <- rep(box.cx, each = 2)
      par(xpd = TRUE)
      for(i in 1:n){
        yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * (i)),
          box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
      }
      par(new = TRUE)
      plot(0, 0, type = "n", ylim = c(min(lev), max(lev)),
        yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
      axis(side = 4, las = 2, tick = FALSE, line = .25)
      par <- opar
    }
    par(mar = c(1.5,1.5,2.5,1.5))
    par(oma = c(1.5,1.5,2.5,1.5))
    scatter3D(df[,1], df[,2], df[,3],
      main=paste(sign_name, "\n", title, sep=""), xlab=xlabel, ylab=ylabel, zlab=zlabel,
      cex.main=title_size,
      box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
      theta=theta, phi=phi,
      pch=16, cex=0.5, alpha=1.0,
      col=colorRampPalette(c("blue", "white", "red"))(50)[as.numeric(cut(df$expr, breaks=50))],
      colvar=NA, colkey=FALSE)
    legend.col(col=colorRampPalette(c("blue", "white", "red"))(50), lev = df$expr)
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_dmap_geneExpression <- function(
  obj, gene_names, zscore = FALSE, data_type_for_dmap, category_for_dmap,
  dim_dmap, theta, phi, title, title_size, label_name, xlabel, ylabel, zlabel
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["reduction"]][["dmap"]][[data_type_for_dmap]][[category_for_dmap]]@eigenvectors
  df <- as.data.frame(mat)
  mat <- as.matrix(obj[["data"]][["normalized"]])
  if(length(gene_names) == 1){
    df$expr <- mat[which(rownames(mat) == gene_names),]
  }else{ 
    df$expr <- apply(mat[which(rownames(mat) %in% gene_names),], 2, mean)
  }
  if(zscore){
    stdev <- sd(df$expr)
    df$expr <- (df$expr - mean(df$expr)) / stdev
  }
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(length(dim_dmap) == 2){
    p <- ggplot() +
      geom_point(aes(x=df[,1], y=df[,2], color=df$expr), size=0.5, alpha=1.0) +
      labs(title=title, x=xlabel, y=ylabel, color=label_name) +
      scale_colour_gradientn(colours=c("blue", "white", "red")) +
      theme_classic(base_size=20, base_family="Helvetica") +
      theme(plot.title=element_text(hjust=0.5, size=title_size),
        plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    return(p)
  }else if(length(dim_dmap) == 3){
    #------------------------------
    # scatter3D in plot3d package
    #------------------------------
    legend.col <- function(col, lev){
      opar <- par
      n <- length(col)
      bx <- par("usr")
      box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
        bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 40)
      box.cy <- c(bx[3], bx[3])
      box.sy <- (bx[4] - bx[3]) / n
      xx <- rep(box.cx, each = 2)
      par(xpd = TRUE)
      for(i in 1:n){
        yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * (i)),
          box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
      }
      par(new = TRUE)
      plot(0, 0, type = "n", ylim = c(min(lev), max(lev)),
        yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
      axis(side = 4, las = 2, tick = FALSE, line = .25)
      par <- opar
    }
    par(mar = c(2,2,3,2))
    par(oma = c(2,2,3,2))
    scatter3D(df[,1], df[,2], df[,3],
      main=title, xlab=xlabel, ylab=ylabel, zlab=zlabel,
      cex.main=title_size,
      box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
      theta=theta, phi=phi,
      pch=16, cex=0.5, alpha=1.0,
      col=colorRampPalette(c("blue", "white", "red"))(50)[as.numeric(cut(df$expr, breaks=50))],
      colvar=NA, colkey=FALSE)
    legend.col(col=colorRampPalette(c("blue", "white", "red"))(50), lev = df$expr)
  }else{
    stop("Error: length(dim_dmap) must be 2 or 3.")
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_dmap_signScore <- function(
  obj, sign_name, data_type_for_dmap, category_for_dmap, data_type_for_expr,
  category_for_expr, dim_dmap, theta, phi, title, title_size, label_name, xlabel,
  ylabel, zlabel
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["reduction"]][["dmap"]][[data_type_for_dmap]][[category_for_dmap]]@eigenvectors
  df <- as.data.frame(mat)
  sign <- "sign"
  slot_name <- paste(data_type_for_expr, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[slot_name]][[category_for_expr]])
  df$expr <- mat[which(rownames(mat) == sign_name),]
#  stdev <- sd(df$expr)
#  df$expr <- (df$expr - mean(df$expr)) / stdev
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(length(dim_dmap) == 2){
    p <- ggplot() +
      geom_point(aes(x=df[,1], y=df[,2], color=df$expr), size=0.5, alpha=1.0) +
      labs(title=title, x=xlabel, y=ylabel, color=label_name) +
      scale_colour_gradientn(colours=c("blue", "white", "red")) +
      theme_classic(base_size=20, base_family="Helvetica") +
      theme(plot.title=element_text(hjust=0.5, size=title_size),
        plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
    return(p)
  }else if(length(dim_dmap) == 3){
    #------------------------------
    # scatter3D in plot3d package
    #------------------------------
    legend.col <- function(col, lev){
      opar <- par
      n <- length(col)
      bx <- par("usr")
      box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000, 
        bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 40)
      box.cy <- c(bx[3], bx[3]) 
      box.sy <- (bx[4] - bx[3]) / n
      xx <- rep(box.cx, each = 2)
      par(xpd = TRUE)
      for(i in 1:n){
        yy <- c(box.cy[1] + (box.sy * (i - 1)), box.cy[1] + (box.sy * (i)),
          box.cy[1] + (box.sy * (i)), box.cy[1] + (box.sy * (i - 1)))
        polygon(xx, yy, col = col[i], border = col[i])
      }
      par(new = TRUE) 
      plot(0, 0, type = "n", ylim = c(min(lev), max(lev)),
        yaxt = "n", ylab = "", xaxt = "n", xlab = "", frame.plot = FALSE)
      axis(side = 4, las = 2, tick = FALSE, line = .25)
      par <- opar
    }
    par(mar = c(2,2,3,2))
    par(oma = c(2,2,3,2))
    scatter3D(df[,1], df[,2], df[,3],
      main=title, xlab=xlabel, ylab=ylabel, zlab=zlabel,
      cex.main=title_size,
      box=TRUE, bty="b2", axes=TRUE, nticks=5,# ticktype="detailed",
      theta=theta, phi=phi,
      pch=16, cex=0.5, alpha=1.0,
      col=colorRampPalette(c("blue", "white", "red"))(50)[as.numeric(cut(df$expr, breaks=50))],
      colvar=NA, colkey=FALSE)
    legend.col(col=colorRampPalette(c("blue", "white", "red"))(50), lev = df$expr)
  }else{
    stop("Error: length(dim_dmap) must be 2 or 3.")
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_pseudotime_vs_signscore <- function(
  obj, sign_name, data_type_for_tree, category_for_tree, data_type_for_expr,
  category_for_expr, title, title_size, label_name, xlabel, ylabel,
  default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  tmp <- obj[["clustering"]][["merlot"]][[data_type_for_tree]][[category_for_tree]]
  slot_name <- paste(data_type_for_expr, "xSample", sep = "")
  mat <- obj[["sign"]][[slot_name]][[category_for_expr]]
  n_groups <- length(tmp[["Pseudotimes"]][["Branches"]])
  #--------------------------------------------------
  # Prepare data frames to be output
  #--------------------------------------------------
  df <- numeric(0)
  df <- data.frame(
    pseudotime = (tmp[["Pseudotimes"]][["Times_cells"]] - min(tmp[["Pseudotimes"]][["Times_cells"]])) /
      (max(tmp[["Pseudotimes"]][["Times_cells"]]) - min(tmp[["Pseudotimes"]][["Times_cells"]])),
    branch = tmp[["Pseudotimes"]][["Cells2Branches"]],
    expression = as.numeric(mat[which(rownames(mat) == sign_name),])
  )
  rownames(df) <- rownames(tmp[["ElasticTree"]][["CellCoords"]])
  dg <- numeric(0)
  for(i in 1:length(tmp[["Pseudotimes"]][["Branches"]])){
    tmp2 <- df[which(df$branch == i),]
    pt <- unique(sort(tmp2$pseudotime))
    for(j in 1:length(pt)){
      dg <- rbind(dg, c(
        pt[j], i, mean(tmp2[which(tmp2$pseudotime == pt[j]),]$expression),
        ifelse(
          length(tmp2[which(tmp2$pseudotime == pt[j]),]$expression) >= 3,
          sd(tmp2[which(tmp2$pseudotime == pt[j]),]$expression),
          0
        )))
    }
  }
  dg <- as.data.frame(dg)
  colnames(dg) <- c("pseudotime", "branch", "expression", "sd")
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(default_color){
    #------------------------------
    # Option 1: ggplot's default
    #------------------------------
    ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
      hues <- seq(15, 375, length=n+1)
      hcl(h=hues, l=l, c=100)[1:n]
    }
    my_colors <- ggColorHue(n_groups)
  }else{
    #------------------------------
    # Option 2: rainbow colors
    #------------------------------
    my_colors <- rainbow(n_groups)
  }
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,3], color=as.factor(df$branch)), size=2, alpha=0.1) +
    geom_line(aes(x=dg[,1], y=dg[,3], color=as.factor(dg[,2])), size=2, alpha=1) +
    geom_ribbon(
      aes(
        x=dg[,1], y=dg[,3], ymin=dg[,3]-dg[,4], ymax=dg[,3]+dg[,4],
        color=as.factor(dg[,2]), fill=as.factor(dg[,2])
      ), size=0, alpha=0.3
    ) +
    scale_colour_manual(values=my_colors) +
    scale_fill_manual(values=my_colors) +
    labs(title=title, x=xlabel, y=ylabel, color=label_name, fill=label_name) +
    theme_classic(base_size=25, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_manifold_sign <- function(
  obj, batch, plot_type, data_type, category, algo_name, title, title_size,
  xlabel, ylabel, default_color
){
  #--------------------------------------------------
  # Prepare a data frame to be output
  #--------------------------------------------------
  if(plot_type == "tsne"){
    mat <- obj[["reduction"]][["tsne"]][[data_type]][[category]][["Y"]][,1:2]
  }else if(plot_type == "umap"){
    mat <- obj[["reduction"]][["umap"]][[data_type]][[category]][["layout"]][,1:2]
  }else if(plot_type == "dmap"){
    mat <- obj[["reduction"]][["dmap"]][[data_type]][[category]]@eigenvectors[,1:2]
  }else{
    stop("Check plot_type.")
  }
  if(batch){
    df <- data.frame(x = mat[,1], y = mat[,2],
                     label = obj[["sample"]][["orig_ident"]])
  }else{
    slot_name <- paste(algo_name, "_", data_type, "_", category, sep = "")
    df <- data.frame(x = mat[,1], y = mat[,2],
                     label = obj[["sample"]][[slot_name]])
  }
  if(default_color){
    #------------------------------
    # Option 1: ggplot's default
    #------------------------------
    n_groups <- length(unique(df$label))
    ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
      hues <- seq(15, 375, length=n+1)
      hcl(h=hues, l=l, c=100)[1:n]
    }
    my_colors <- ggColorHue(n_groups)
  }else{
    #------------------------------
    # Option 2: rainbow colors
    #------------------------------
    n_groups <- length(unique(df$label))
    my_colors <- rainbow(n_groups)
  }
  if(batch){
    #------------------------------
    # ggplot
    #------------------------------
    dg <- df[which(df$label == sort(unique(df$label))[2]),]
    p <- ggplot() +
      geom_point(aes(x=df[,1], y=df[,2], color=as.factor(df$label)), size=0.1, alpha=1.0) +
      geom_point(aes(x=dg[,1], y=dg[,2], color=as.factor(dg$label)), size=0.1, alpha=1.0) +
      labs(title=title, x=xlabel, y=ylabel, colour="") +
      scale_colour_manual(values=my_colors) +
      theme_classic(base_size=20, base_family="Helvetica") +
      theme(plot.title=element_text(hjust=0.5, size=title_size),
        plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
      guides(colour = guide_legend(override.aes = list(size=4)))
  }else{
    #------------------------------
    # Positions of labels
    #------------------------------
    cluster_num <- sort(unique(df$label))
    cog <- c()  # Center of gravity
    for(i in cluster_num){
      cog <- rbind(cog, cbind(
        mean(df[which(df$label == i),][,1]), mean(df[which(df$label == i),][,2])))
    }
    cog <- data.frame(x = cog[,1], y = cog[,2], label = cluster_num)
    #------------------------------
    # ggplot
    #------------------------------
    p <- ggplot() +
      geom_point(aes(x=df[,1], y=df[,2], color=as.factor(df$label)), size=0.1, alpha=1.0) +
      labs(title=title, x=xlabel, y=ylabel, colour="") +
      scale_colour_manual(values=my_colors) +
      theme_classic(base_size=20, base_family="Helvetica") +
      theme(plot.title=element_text(hjust=0.5, size=title_size),
        plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
      guides(colour = guide_legend(override.aes = list(size=4)))
  }

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_MultiHeatmaps_SignxSamp <- function(
  obj, data_types, categories, algo_names, show_labels, splitting, method,
  show_nReads, title, names, show_rownames_sign, show_rownames_label,
  show_rownames_nReads, default_colors
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  sign <- "sign"
  mat <- list()
  labels <- list()
  for(i in 1:length(data_types)){
    slot_name <- paste(data_types[i], "xSample", sep = "")
    mat[[i]] <- as.matrix(obj[[sign]][[slot_name]][[categories[i]]])
    slot_name <- paste(algo_names[i], "_", data_types[i], "_", categories[i], sep = "")
    labels[[i]] <- obj[["sample"]][[slot_name]]
  }
  set.seed(8)
  col_hc <- hclust(dist(t(mat[[1]])), method = method)
  #--------------------------------------------------
  # Heatmap
  # Note: the positions of dendrogram leaves are
  # slightly adjusted by the gaps between slices.
  #--------------------------------------------------
  ht_opt$message = FALSE
  i <- 1
  if(splitting){
    if(show_labels[i]){
      if(default_colors[i]){
        #------------------------------
        # Option 1: ggplot's default
        #------------------------------
        tmp <- unique(sort(labels[[i]]))
        ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
          hues <- seq(15, 375, length=n+1)
          hcl(h=hues, l=l, c=100)[1:n]
        }
        my_colors <- ggColorHue(length(tmp))
        names(my_colors) <- tmp
      }else{
        #------------------------------
        # Option 2: rainbow colors
        #------------------------------
        tmp <- unique(sort(labels[[i]]))
        my_colors <- rainbow(length(tmp))
        names(my_colors) <- tmp
      }
      ha <- HeatmapAnnotation(
        Label = labels[[i]], col = list(Label = my_colors),
        annotation_name_side = "left"
      )
      column_split <- labels[[i]]
      p <- Heatmap(
        mat[[i]],
        column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
        name=names[i],
        column_split=column_split, column_gap=unit(1.5, "mm"),
        border=FALSE,
        show_row_names=show_rownames_sign, row_names_side="right", show_row_dend=FALSE,
#        row_title=paste("Sign (", data_types[i], ")", sep = ""),
        show_column_names=FALSE, column_dend_side="top",
        show_parent_dend_line=FALSE, top_annotation=ha
      )
    }else{
      p <- Heatmap(
        mat[[i]],
        column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
        name=names[i],
        column_split=column_split, column_gap=unit(1.5, "mm"),
        border=FALSE,
        show_row_names=show_rownames_sign, row_names_side="right", show_row_dend=FALSE,
        row_title=paste("Sign (", data_types[i], ")", sep = ""),
        show_column_names=FALSE, column_dend_side="top",
        show_parent_dend_line=FALSE
      )
    }
  }else{
    if(show_labels[i]){
      tmp <- unique(sort(labels[[i]]))
      my_colors <- rainbow(length(tmp))
      names(my_colors) <- tmp
      ha <- HeatmapAnnotation(
        Label = labels[[i]], col = list(Label = my_colors),
        annotation_name_side = "left"
      )
      p <- Heatmap(
        mat[[i]],
        column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
        name=names[i],
        cluster_columns=col_hc,
        show_row_names=show_rownames_sign, row_names_side="right", row_dend_side="left",
#        row_title=paste("Sign (", data_types[i], ")", sep = ""),
        show_column_names=FALSE, column_dend_side="top",
        show_parent_dend_line=FALSE, top_annotation=ha
      )
    }else{
      p <- Heatmap(
        mat[[i]],
        column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
        name=names[i],
        cluster_columns=col_hc,
        show_row_names=show_rownames_sign, row_names_side="right", show_row_dend=FALSE,
#        row_title=paste("Sign (", data_types[i], ")", sep = ""),
        show_column_names=FALSE, column_dend_side="top",
        show_parent_dend_line=FALSE
      )
    }
  }
  for(i in 2:length(data_types)){
    if(show_labels[i]){
      if(default_colors[i]){
        #------------------------------
        # Option 1: ggplot's default
        #------------------------------
        tmp <- unique(sort(labels[[i]]))
        ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
          hues <- seq(15, 375, length=n+1)
          hcl(h=hues, l=l, c=100)[1:n]
        }
        my_colors <- ggColorHue(length(tmp))
        names(my_colors) <- tmp
      }else{
        #------------------------------
        # Option 2: rainbow colors
        #------------------------------
        tmp <- unique(sort(labels[[i]]))
        my_colors <- rainbow(length(tmp))
        names(my_colors) <- tmp
      }
      if(i == 2){
        ha <- HeatmapAnnotation(
          Label_2 = labels[[i]], col = list(Label_2 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 3){
        ha <- HeatmapAnnotation(
          Label_3 = labels[[i]], col = list(Label_3 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 4){
        ha <- HeatmapAnnotation(
          Label_4 = labels[[i]], col = list(Label_4 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 5){
        ha <- HeatmapAnnotation(
          Label_5 = labels[[i]], col = list(Label_5 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 6){
        ha <- HeatmapAnnotation(
          Label_6 = labels[[i]], col = list(Label_6 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 7){
        ha <- HeatmapAnnotation(
          Label_7 = labels[[i]], col = list(Label_7 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 8){
        ha <- HeatmapAnnotation(
          Label_8 = labels[[i]], col = list(Label_8 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 9){
        ha <- HeatmapAnnotation(
          Label_9 = labels[[i]], col = list(Label_9 = my_colors), annotation_name_side = "left"
        )
      }
      q <- Heatmap(
        mat[[i]],
        column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
        name=names[i],
        cluster_columns=col_hc,
        show_row_names=show_rownames_sign, row_names_side="right", show_row_dend=FALSE,
#        row_title=paste("Sign (", data_types[i], ")", sep = ""),
        show_column_names=FALSE, column_dend_side="top",
        show_parent_dend_line=FALSE, top_annotation=ha
      )
    }else{
       q <- Heatmap(
        mat[[i]],
        column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
        name=names[i],
        cluster_columns=col_hc,
        show_row_names=show_rownames_sign, row_names_side="right", show_row_dend=FALSE,
#        row_title=paste("Sign (", data_types[i], ")", sep = ""),
        show_column_names=FALSE, column_dend_side="top",
        show_parent_dend_line=FALSE
      )
    }
    p <- p %v% q
  }
  if(show_nReads){
    mtx <- t(obj[["sample"]][["nReads"]])
    rownames(mtx) <- "nReads"
    q <- Heatmap(
      mtx,
      name="nReads",
      cluster_columns=col_hc,
      show_row_names=show_rownames_nReads, row_names_side="left", show_row_dend=FALSE,
      show_column_names=FALSE, show_column_dend=FALSE,
      col=colorRamp2(c(min(mtx), max(mtx)), c("cyan", "magenta"))
    )
    p <- p %v% q
  }

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_MultiBargraphs_sign <- function(
  obj, data_type_1, category_1, algo_name_1, data_type_2, category_2, algo_name_2,
  title, title_size, cbar_title, xlabel, ylabel, ymax, default_color
){
  #--------------------------------------------------
  # Set data frame to be output
  #--------------------------------------------------
  slot_name_1 <- paste(algo_name_1, "_", data_type_1, "_", category_1, sep = "")
  slot_name_2 <- paste(algo_name_2, "_", data_type_2, "_", category_2, sep = "")
  tmp <- obj[["sample"]]
  colnames(tmp)[which(colnames(tmp) == slot_name_1)] <- "cnt_1"
  colnames(tmp)[which(colnames(tmp) == slot_name_2)] <- "cnt_2"
  df <- tmp %>% group_by(cnt_1, cnt_2) %>% summarise (n=n())
  #--------------------------------------------------
  # Error control
  #--------------------------------------------------
  y <- 0
  n_labels <- length(unique(df$cnt_1))
  for(i in 1:n_labels){
    y1 <- sum(df[which(df$cnt_1 == i),]$n)
    if(y1 >= y){
      y <- y1
    }
  }
  if(y > ymax){
    stop("ymax must be greater than or equal to the maximum value of y.")
  }
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(default_color){
    #------------------------------
    # Option 1: ggplot's default
    #------------------------------
    n_groups <- length(unique(df$cnt_2))
    ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
      hues <- seq(15, 375, length=n+1)
      hcl(h=hues, l=l, c=100)[1:n]
    }
    my_colors <- ggColorHue(n_groups)
  }else{
    #------------------------------
    # Option 2: rainbow colors
    #------------------------------
    n_groups <- length(unique(df$cnt_2))
    my_colors <- rainbow(n_groups)
  }
  p <- ggplot() +
    geom_bar(aes(x=as.factor(df$cnt_1), y=df$n, fill=as.factor(df$cnt_2)),
      stat="identity", width=0.7) +
    theme_classic(base_size=18, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      axis.text.x=element_text(vjust=0.5, angle=0)) +
    ylim(0,ymax) + labs(title=title, x=xlabel, y=ylabel, fill=cbar_title)
  p <- p + scale_fill_manual(values=my_colors)

  return(p)
} 
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_MultiHeatmaps_SignxSamp_GenexSamp <- function(
  obj, topn_genes, topn_signs, data_types, categories, algo_names, show_labels,
  method, show_nReads, title, names_signs, name_genes,
  show_rownames_sign, show_rownames_label, show_rownames_nReads, default_colors,
  show_cellcycle
){
  #--------------------------------------------------
  # Definition for sign-by-sample matrices
  #--------------------------------------------------
  sign <- "sign"
  mat <- list()
  col_labels <- list()
  row_labels <- list()
  for(i in 1:length(data_types)){
    slot_name <- paste(data_types[i], "xSample", sep = "")
    mat[[i]] <- as.matrix(obj[[sign]][[slot_name]][[categories[i]]])
    mat[[i]] <- mat[[i]][which(rownames(mat[[i]]) %in% topn_signs[[i]]$sign), ]
    slot_name <- paste(algo_names[i], "_", data_types[i], "_", categories[i], sep = "")
    col_labels[[i]] <- obj[["sample"]][[slot_name]]
    row_labels[[i]] <- topn_signs[[i]]$label
  }
  tmp <- list()
  soi <- list()
  for(i in 1:length(data_types)){
    tmp2 <- c()
    soi2 <- c()
    for(j in 1:nrow(topn_signs[[i]])){
      tmp2 <- rbind(tmp2, mat[[i]][which(rownames(mat[[i]]) == topn_signs[[i]]$sign[j]),])
      soi2 <- rbind(soi2, topn_signs[[i]]$sign[j])
    }
    rownames(tmp2) <- soi2
    tmp[[i]] <- tmp2
    soi[[i]] <- soi2
  }
  mat <- tmp
  
  set.seed(8)
  col_hc <- hclust(dist(t(mat[[1]])), method = method)
  #--------------------------------------------------
  # Heatmap (sign-by-sample matrices)
  # Note: the positions of dendrogram leaves are
  # slightly adjusted by the gaps between slices.
  #--------------------------------------------------
  ht_opt$message = FALSE
  i <- 1
  if(show_labels[i]){
    if(default_colors[i]){
      #------------------------------
      # Option 1: ggplot's default
      #------------------------------
      tmp <- unique(sort(col_labels[[i]]))
      ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
        hues <- seq(15, 375, length=n+1)
        hcl(h=hues, l=l, c=100)[1:n]
      }
      my_colors <- ggColorHue(length(tmp))
      names(my_colors) <- tmp
    }else{
      #------------------------------
      # Option 2: rainbow colors
      #------------------------------
      tmp <- unique(sort(col_labels[[i]]))
      my_colors <- rainbow(length(tmp))
      names(my_colors) <- tmp
    }
    ha <- HeatmapAnnotation(
      Label = col_labels[[i]], col = list(Label = my_colors),
      annotation_name_side = "left"
    )
    column_split <- col_labels[[i]]
    p <- Heatmap(
      mat[[i]],
      column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
      name=names_signs[i],
      column_split=column_split, column_gap=unit(1.5, "mm"),
      row_split=row_labels[[i]], row_gap=unit(1.0, "mm"),
      border=FALSE,
      show_row_names=show_rownames_sign, row_names_side="left", show_row_dend=FALSE,
#      row_title=paste("Sign (", data_types[i], ")", sep = ""),
      show_column_names=FALSE, column_dend_side="top",
      show_parent_dend_line=FALSE, top_annotation=ha
    )
  }
  for(i in 2:length(data_types)){
    if(show_labels[i]){
      if(default_colors[i]){
        #------------------------------
        # Option 1: ggplot's default
        #------------------------------
        tmp <- unique(sort(col_labels[[i]]))
        ggColorHue <- function(n, l=65){   # To produce default colors of ggplot2
          hues <- seq(15, 375, length=n+1)
          hcl(h=hues, l=l, c=100)[1:n]
        }
        my_colors <- ggColorHue(length(tmp))
        names(my_colors) <- tmp
      }else{
        #------------------------------
        # Option 2: rainbow colors
        #------------------------------
        tmp <- unique(sort(col_labels[[i]]))
        my_colors <- rainbow(length(tmp))
        names(my_colors) <- tmp
      }
      if(i == 2){
        ha <- HeatmapAnnotation(
          Label_2 = col_labels[[i]], col = list(Label_2 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 3){
        ha <- HeatmapAnnotation(
          Label_3 = col_labels[[i]], col = list(Label_3 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 4){
        ha <- HeatmapAnnotation(
          Label_4 = col_labels[[i]], col = list(Label_4 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 5){
        ha <- HeatmapAnnotation(
          Label_5 = col_labels[[i]], col = list(Label_5 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 6){
        ha <- HeatmapAnnotation(
          Label_6 = col_labels[[i]], col = list(Label_6 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 7){
        ha <- HeatmapAnnotation(
          Label_7 = col_labels[[i]], col = list(Label_7 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 8){
        ha <- HeatmapAnnotation(
          Label_8 = col_labels[[i]], col = list(Label_8 = my_colors), annotation_name_side = "left"
        )
      }
      if(i == 9){
        ha <- HeatmapAnnotation(
          Label_9 = col_labels[[i]], col = list(Label_9 = my_colors), annotation_name_side = "left"
        )
      }
      q <- Heatmap(
        mat[[i]],
        column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
        name=names_signs[i],
        cluster_columns=col_hc,
        show_row_names=show_rownames_sign, row_names_side="left", show_row_dend=FALSE,
#        row_title=paste("Sign (", data_types[i], ")", sep = ""),
        show_column_names=FALSE, column_dend_side="top",
        show_parent_dend_line=FALSE, top_annotation=ha
      )
    }else{
      q <- Heatmap(
        mat[[i]],
        column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
        name=names_signs[i],
        cluster_columns=col_hc,
        row_split=row_labels[[i]], row_gap=unit(1.0, "mm"),
        show_row_names=show_rownames_sign, row_names_side="left", show_row_dend=FALSE,
#        row_title=paste("Sign (", data_types[i], ")", sep = ""),
        show_column_names=FALSE, column_dend_side="top",
        show_parent_dend_line=FALSE
      )
    }
    p <- p %v% q
  }
  #--------------------------------------------------
  # Definition for gene-by-sample matrices
  #--------------------------------------------------
  res <- c()
  mat <- t(scale(t(obj[["data"]][["normalized"]])))
  if(nrow(topn_genes) == 1){
    mat <- t(as.matrix(mat[which(rownames(mat) %in% topn_genes$gene),]))
    rownames(mat) <- topn_genes$gene
  }else{
    mat <- as.matrix(mat[which(rownames(mat) %in% topn_genes$gene),])
  }
  tmp <- c()
  goi <- c()
  for(i in 1:length(topn_genes$gene)){
    if(topn_genes$gene[i] %in% rownames(mat)){
      tmp <- rbind(tmp, mat[which(rownames(mat) == topn_genes$gene[i]),])
      goi <- rbind(goi, topn_genes$gene[i])
    }
  }
  rownames(tmp) <- goi
  mat <- as.matrix(tmp)
  row_labels <- topn_genes$label
  #--------------------------------------------------
  # Heatmap (gene-by-sample matrices)
  # Note: the positions of dendrogram leaves are
  # slightly adjusted by the gaps between slices.
  #--------------------------------------------------  
  if(show_cellcycle){
    tmp <- unique(sort(obj[["sample"]][["phase"]]))
    my_colors <- rainbow(length(tmp)) 
#    my_colors <- c("#0000FF", "#00FF00", "#FF0000") 
    names(my_colors) <- tmp 
    ha <- HeatmapAnnotation(
      Label_cc = obj[["sample"]][["phase"]], col = list(Label_cc = my_colors), 
      annotation_name_side = "left" 
    ) 
    q <- Heatmap(
      mat,
      column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
      name=name_genes,
      cluster_columns=col_hc,
      row_split=row_labels, row_gap=unit(1.0, "mm"),
      border=FALSE,
      show_row_names=TRUE, row_names_side="left", show_row_dend=FALSE,
      show_column_names=FALSE, column_dend_side="top",
      show_parent_dend_line=FALSE, top_annotation=ha
    )
  }else{
    q <- Heatmap(
      mat,
      column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
      name=name_genes,
      cluster_columns=col_hc,
      row_split=row_labels, row_gap=unit(1.0, "mm"),
      border=FALSE,
      show_row_names=TRUE, row_names_side="left", show_row_dend=FALSE,
      show_column_names=FALSE, column_dend_side="top",
      show_parent_dend_line=FALSE
    )
  }
  p <- p %v% q
  if(show_nReads){
    mtx <- t(obj[["sample"]][["nReads"]])
    rownames(mtx) <- "nReads"
    q <- Heatmap(
      mtx,
      name="nReads",
      cluster_columns=col_hc,
      show_row_names=show_rownames_nReads, row_names_side="left", show_row_dend=FALSE,
      show_column_names=FALSE, show_column_dend=FALSE,
      col=colorRamp2(c(min(mtx), max(mtx)), c("cyan", "magenta"))
    )
    p <- p %v% q
  }

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80

