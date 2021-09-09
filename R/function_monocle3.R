#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_tsne_monocle3 <- function(
  obj, title, title_size, xlabel, ylabel, default_color
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj@int_colData@listData[["reducedDims"]]@listData[["tSNE"]]
  df <- as.data.frame(mat)
  df$label <- obj@clusters@listData[["tSNE"]][["clusters"]]
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
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
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_umap_monocle3 <- function(
  obj, title, title_size, xlabel, ylabel, default_color
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj@int_colData@listData[["reducedDims"]]@listData[["UMAP"]]
  df <- as.data.frame(mat)
  df$label <- obj@clusters@listData[["UMAP"]][["clusters"]]
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
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
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_dmap_monocle3 <- function(
  obj, sigma = "local", distance = "euclidean", pca_dim
){
  #--------------------------------------------------
  # DiffusionMap
  #--------------------------------------------------
  set.seed(8)
  if(is.null(pca_dim)){
    mat <- as.matrix(obj@assays@data@listData[["logcounts"]])
    res <- DiffusionMap(t(mat), sigma = sigma, distance = distance)
  }else{
    mat <- obj@int_colData@listData[["reducedDims"]]@listData[["PCA"]]
    res <- DiffusionMap(mat[,1:pca_dim], sigma = sigma, distance = distance)
  }
  obj@int_colData@listData[["reducedDims"]]@listData[["dmap"]] <- res@eigenvectors

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_dmap_monocle3 <- function(
  obj, title, title_size, xlabel, ylabel, default_color
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj@int_colData@listData[["reducedDims"]]@listData[["dmap"]]
  df <- as.data.frame(mat)
  df$label <- obj@clusters@listData[["UMAP"]][["clusters"]]
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
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
}

