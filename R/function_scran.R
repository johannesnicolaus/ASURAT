#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
process_001_scran <- function(obj_sce){
  #--------------------------------------------------
  # Create a list
  #--------------------------------------------------
  tmp <- obj_sce
  obj <- list()
  obj[["sce"]] <- tmp
  #--------------------------------------------------
  # Quick clustering of cells
  #--------------------------------------------------
  clusters <- quickCluster(obj[["sce"]])
  #--------------------------------------------------
  # Normalize
  #--------------------------------------------------
  obj[["sce"]] <- computeSumFactors(obj[["sce"]], clusters = clusters)
  obj[["sce"]] <- logNormCounts(obj[["sce"]])
  #--------------------------------------------------
  # Perform a variance modeling
  #--------------------------------------------------
  obj[["dec"]] <- modelGeneVar(obj[["sce"]])

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_variance_modeling_scran <- function(
  obj, title, title_size, xlabel, ylabel
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df1 <- data.frame(x = obj[["dec"]]$mean, y = obj[["dec"]]$total)
  df2 <- data.frame(x = obj[["dec"]]$mean, y = obj[["dec"]]$tech)
  #--------------------------------------------------
  # plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df1[,1], y=df1[,2]), color="black", alpha=1, size=2) +
    geom_line(aes(x=df2[,1], y=df2[,2]), color="red", alpha=1, size=2) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
          plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="none")

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_pcaEigen_scran <- function(obj, title, title_size, xlabel, ylabel){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  tmp <- obj[["sce"]]@int_colData@listData[["reducedDims"]]@listData[["PCA"]]
  tmp <- apply(tmp, 2, sd)
  df <- data.frame(num = 1:length(tmp), eigen = tmp)
  #--------------------------------------------------
  # plot
  #--------------------------------------------------
  p <- ggplot(df, aes(x=df[,1], y=df[,2])) +
    geom_point(shape=21, color="black", alpha=1, size=1, stroke=1) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=15, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="none")
  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
process_002_scran <- function(obj, auto_pca = FALSE, num_hvg, k){
  #--------------------------------------------------
  # Denoise PCA
  #--------------------------------------------------
  if(auto_pca){
    obj[["sce"]] <- denoisePCA(obj[["sce"]], obj[["dec"]], subset.row = num_hvg)
  }
  #--------------------------------------------------
  # KNN graph-based clustering
  #--------------------------------------------------
  set.seed(8)
  g <- buildSNNGraph(obj[["sce"]], use.dimred = "PCA", k = k, type = "rank")
  c <- igraph::cluster_louvain(g)$membership
  obj[["sce"]][["cluster"]] <- as.factor(c)

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_tsne_scran <- function(obj, pca_dim, tsne_dim = 2){
  set.seed(8)
  if(is.null(pca_dim)){
    mat <- obj[["sce"]]@assays@data@listData[["logcounts"]]
    res <- Rtsne(t(mat), dim = tsne_dim, pca = FALSE)
  }else{
    mat <- obj[["sce"]]@int_colData@listData[["reducedDims"]]@listData[["PCA"]]
    res <- Rtsne(mat[,1:pca_dim], dim = tsne_dim, pca = FALSE)
  }
  obj[["sce"]]@int_colData@listData[["reducedDims"]]@listData[["tsne"]] <- res

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_tsne_scran <- function(
  obj, title, title_size, xlabel, ylabel, default_color
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["sce"]]@int_colData@listData[["reducedDims"]]@listData[["tsne"]][["Y"]]
  df <- as.data.frame(mat)
  df$label <- obj[["sce"]][["cluster"]]
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
do_umap_scran <- function(obj, pca_dim, umap_dim = 2){
  set.seed(8)
  if(is.null(pca_dim)){
    mat <- obj[["sce"]]@assays@data@listData[["logcounts"]]
    res <- umap(t(mat), n_components = umap_dim)
  }else{
    mat <- obj[["sce"]]@int_colData@listData[["reducedDims"]]@listData[["PCA"]]
    res <- umap(mat[,1:pca_dim], n_components = umap_dim)
  }
  obj[["sce"]]@int_colData@listData[["reducedDims"]]@listData[["umap"]] <- res
  
  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_umap_scran <- function(
  obj, title, title_size, xlabel, ylabel, default_color
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["sce"]]@int_colData@listData[["reducedDims"]]@listData[["umap"]][["layout"]]
  df <- as.data.frame(mat)
  df$label <- obj[["sce"]][["cluster"]]
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
do_dmap_scran <- function(
  obj, sigma = "local", distance = "euclidean", pca_dim
){
  #--------------------------------------------------
  # DiffusionMap
  #--------------------------------------------------
  set.seed(8)
  if(is.null(pca_dim)){
    mat <- obj[["sce"]]@assays@data@listData[["logcounts"]]
    res <- DiffusionMap(t(mat), sigma = sigma, distance = distance)
  }else{
    mat <- obj[["sce"]]@int_colData@listData[["reducedDims"]]@listData[["PCA"]]
    res <- DiffusionMap(mat[,1:pca_dim], sigma = sigma, distance = distance)
  }
  obj[["sce"]]@int_colData@listData[["reducedDims"]]@listData[["dmap"]] <- res

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_dmap_scran <- function(
  obj, title, title_size, xlabel, ylabel, default_color
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- obj[["sce"]]@int_colData@listData[["reducedDims"]]@listData[["dmap"]]@eigenvectors
  df <- as.data.frame(mat)
  df$label <- obj[["sce"]][["cluster"]]
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
process_003_scran <- function(obj, direction = "up", pval.type = "all"){
  #------------------------------
  # Pairwise T test
  #------------------------------
  obj[["pwtt"]] <- pairwiseTTests(
    x = logcounts(obj[["sce"]]), groups = obj[["sce"]][["cluster"]],
    direction = "up"
  )
  obj[["cmb"]] <- combineMarkers(
    de.lists = obj[["pwtt"]][["statistics"]], pairs = obj[["pwtt"]][["pairs"]],
    pval.type = "all"
  )
  #------------------------------
  # Create a result table
  #------------------------------
  res <- list()
  for(i in 1:length(obj[["cmb"]]@listData)){
    tmp <- data.frame(
      label = as.integer(i),
      gene = obj[["cmb"]]@listData[[i]]@rownames,
      pval = obj[["cmb"]]@listData[[i]]@listData[["p.value"]],
      FDR = obj[["cmb"]]@listData[[i]]@listData[["FDR"]],
      logFC = obj[["cmb"]]@listData[["1"]]@listData[["summary.logFC"]]
    )
    res[[i]] <- tmp
  }
  obj[["stat"]] <- res
  #------------------------------
  # Concatenate all the results
  #------------------------------
  res <- c()
  for(i in 1:length(obj[["stat"]])){
    res <- rbind(res, obj[["stat"]][[i]])
  }
  obj[["stat"]][["all"]] <- res

  return(obj)
}

