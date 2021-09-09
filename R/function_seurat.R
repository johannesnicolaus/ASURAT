#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
process_001_seurat <- function(obj, nfeatures){
  #--------------------------------------------------
  # Normalize the data
  #--------------------------------------------------
  obj <- NormalizeData(obj, normalization.method = "LogNormalize")
  #--------------------------------------------------
  # Perform a variance stabilizing transform (VST).
  # As mentioned in Cruz and Wishart, Cancer Inform. 2, 59-77, 2006.
  # the sample-per-variable feature ratio is set as 5:1.
  # n <- round(0.2 * obj@assays[["RNA"]]@counts@Dim[2])
  #--------------------------------------------------
  n <- nfeatures
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n)
  #--------------------------------------------------
  # Scale the data
  #--------------------------------------------------
  obj <- ScaleData(obj)
  #--------------------------------------------------
  # Principal component analysis
  #--------------------------------------------------
  obj <- RunPCA(obj, features = VariableFeatures(obj))

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
process_002_seurat <- function(obj, pc, resolution){
  set.seed(8)
  #--------------------------------------------------
  # Sample clustering
  #--------------------------------------------------
  obj <- FindNeighbors(obj, reduction = "pca", dim = 1:pc)
  obj <- FindClusters(obj, resolution = resolution)
  #--------------------------------------------------
  # t-SNE
  #--------------------------------------------------
  obj <- RunTSNE(
    obj, dims.use = 1:2, reduction = "pca", dims = 1:pc,
    do.fast = FALSE, perplexity = 30
  )
  #--------------------------------------------------
  # UMAP
  #--------------------------------------------------
  obj <- RunUMAP(obj, dims = 1:pc)

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_IntegrateData_seurat <- function(
  obj_1, obj_2, nfeatures_1, nfeatures_2, nfeatures,
  normalization.method = LogNormalize
){
  # ----------------------------------------
  # Set up Seurat objects using normalized data
  # ----------------------------------------
  obj_1_seurat <- CreateSeuratObject(
    counts = obj_1[["data"]][["bayNorm"]][["Bay_out"]],
    project = obj_1[["history"]][["make_asurat_obj"]][["obj_name"]])
  obj_2_seurat <- CreateSeuratObject(
    counts = obj_2[["data"]][["bayNorm"]][["Bay_out"]],
    project = obj_2[["history"]][["make_asurat_obj"]][["obj_name"]])
  obj_1_seurat <- NormalizeData(obj_1_seurat,
                                normalization.method = normalization.method,
                                scale.factor = 10000)
  obj_2_seurat <- NormalizeData(obj_2_seurat,
                                normalization.method = normalization.method,
                                scale.factor = 10000)
  # ****************************************
  # Memo: log(x + 1) and NormalizeData()
  # produce totally different results
  # ****************************************
  obj_1_seurat <- FindVariableFeatures(obj_1_seurat, selection.method = "vst",
                                       nfeatures = nfeatures_1)
  obj_2_seurat <- FindVariableFeatures(obj_2_seurat, selection.method = "vst",
                                       nfeatures = nfeatures_2)
  # ----------------------------------------
  # Select features that are repeatedly variable
  # across datasets for integration
  # ----------------------------------------
  features <- SelectIntegrationFeatures(
    object.list = list(obj_1_seurat, obj_2_seurat), nfeatures = nfeatures)
  # ----------------------------------------
  # Identify anchors
  # ----------------------------------------
  pdac_a_anchors <- FindIntegrationAnchors(
    object.list = list(obj_1_seurat, obj_2_seurat), anchor.features = features)
  # ----------------------------------------
  # Create an integrated data assay
  # ----------------------------------------
  obj_seurat <- IntegrateData(anchorset = pdac_a_anchors,
                              normalization.method = normalization.method)
  DefaultAssay(obj_seurat) <- "integrated"

  return(obj_seurat)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_dmap_seurat <- function(
  obj, sigma = "local", distance = "euclidean", pca_dim
){
  #--------------------------------------------------
  # DiffusionMap
  #--------------------------------------------------
  set.seed(8)
  if(is.null(pca_dim)){
    mat <- as.matrix(obj@assays[["RNA"]]@scale.data)
    res <- DiffusionMap(t(mat), sigma = sigma, distance = distance)
  }else{
    mat <- as.matrix(obj@reductions[["pca"]]@cell.embeddings)
    res <- DiffusionMap(mat[,1:pca_dim], sigma = sigma, distance = distance)
  }
  obj@reductions[["dmap"]] <- res

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_enrichGO_seurat <- function(
  obj, padj_cutoff, qval_cutoff_enrichGO, orgdb, cutoff, ont
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  tmp <- obj@misc[["markers"]]
  category_names <- c("MF", "BP", "CC", "ALL")
  res <- list()
  cluster_names <- as.character(unique(tmp$cluster))
  #--------------------------------------------------
  # enrichGO()
  #--------------------------------------------------
  for(cluster in cluster_names){
    for(category in category_names){
      df <- tmp[which(tmp$cluster == cluster),]
      goi <- df[which(df$p_val_adj <= padj_cutoff), ]$gene
      goi <- bitr(goi, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
      ego <- enrichGO(
        gene = goi$ENTREZID, OrgDb = orgdb, ont = category,
        pAdjustMethod = "BH", qvalueCutoff = qval_cutoff_enrichGO, readable = TRUE
      )
      res[[cluster]][[category]] <- ego
    }
  }
  obj@misc[["enrichGO"]] <- res
  #--------------------------------------------------
  # simplify()
  #--------------------------------------------------
  couster_names <- names(res)
  for(cluster in cluster_names){
    df <- res[[cluster]]
    category_names_woALL <- setdiff(names(df), "ALL")
    for(category in category_names_woALL){
      res[[cluster]][[category]] <- simplify(
        x = df[[category]],
        cutoff = cutoff,
        by = "p.adjust",
        select_fun = min,
        measure = "Wang",
        semData = NULL
      )
    }
  }
  obj@misc[["enrichGO_simplified"]] <- res

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_variableFeatures_seurat <- function(obj, n, title, title_size, xlabel, ylabel){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  top_n <- head(VariableFeatures(obj), n = n)
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- LabelPoints(plot=VariableFeaturePlot(obj), points=top_n, repel=TRUE, xnudge=0, ynudge=0) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=16, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_pcaEigen_seurat <- function(obj, title, title_size, xlabel, ylabel){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  pca_eigen <- apply(obj@reductions[["pca"]]@cell.embeddings, 2, sd)
  I <- length(pca_eigen)
  df <- data.frame(dim = 1:I, eigen = pca_eigen)
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(
      aes(x=df$dim, y=df$eigen),
      shape=21, color="black", fill="white", alpha=1, size=1, stroke=1
    ) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=15, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="none")

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_bargraph_seurat <- function(obj, title, title_size, xlabel, ylabel, ymax){
  #--------------------------------------------------
  # Set data frame to be output
  #--------------------------------------------------
  tmp <- as.data.frame(obj@active.ident)
  names(tmp) <- "Count"
  df <- tmp %>% group_by(Count) %>% summarise (n=n())
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot(df) +
    geom_bar(aes(x=as.factor(Count), y=n),
      color="black", fill="black", stat="identity", width=0.7) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size)) +
    geom_text(aes(x=as.factor(Count), y=n, label=sprintf("%d", n)), size=6, vjust=-0.5) +
    ylim(0,ymax) + labs(title=title, x=xlabel, y=ylabel)

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_tsne_seurat <- function(
  obj, title, title_size, xlabel, ylabel, default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df <- data.frame(
    x = obj@reductions[["tsne"]]@cell.embeddings[,1],
    y = obj@reductions[["tsne"]]@cell.embeddings[,2],
    label = obj@active.ident
  )
  colnames(df) <- c("x", "y", "label")
  n_groups <- length(unique(obj@active.ident))
  #--------------------------------------------------
  # Positions of labels
  #--------------------------------------------------
  cluster_num <- as.integer(as.character(sort(unique(df$label))))
  cog <- c()  # Center of gravity
  for(i in cluster_num){
    cog <- rbind(cog, cbind(
      mean(df[which(df$label == i),]$x), mean(df[which(df$label == i),]$y)))
  }
  cog <- data.frame(x = cog[,1], y = cog[,2], label = cluster_num)
  #--------------------------------------------------
  # Color
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
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, colour="") +
    scale_colour_manual(values=my_colors) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
    geom_text(aes(x=cog[,1], y=cog[,2], label=sprintf("%d", cluster_num)), size=6) +
    guides(colour = guide_legend(override.aes = list(size=4)))

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_umap_seurat <- function(
  obj, title, title_size, xlabel, ylabel, default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df <- data.frame(
    x = obj@reductions[["umap"]]@cell.embeddings[,1],
    y = obj@reductions[["umap"]]@cell.embeddings[,2],
    label = obj@active.ident
  )
  colnames(df) <- c("x", "y", "label")
  n_groups <- length(unique(obj@active.ident))
  #--------------------------------------------------
  # Positions of labels
  #--------------------------------------------------
  cluster_num <- as.integer(as.character(sort(unique(df$label))))
  cog <- c()  # Center of gravity
  for(i in cluster_num){
    cog <- rbind(cog, cbind(
      mean(df[which(df$label == i),]$x), mean(df[which(df$label == i),]$y)))
  }
  cog <- data.frame(x = cog[,1], y = cog[,2], label = cluster_num)
  #--------------------------------------------------
  # Color
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
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, colour="") +
    scale_colour_manual(values=my_colors) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
    geom_text(aes(x=cog[,1], y=cog[,2], label=sprintf("%d", cluster_num)), size=6) +
    guides(colour = guide_legend(override.aes = list(size=4)))

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_dmap_seurat <- function(
  obj, title, title_size, xlabel, ylabel, default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df <- data.frame(
    x = obj@reductions[["dmap"]]@eigenvectors[,1],
    y = obj@reductions[["dmap"]]@eigenvectors[,2],
    label = obj@active.ident
  )
  colnames(df) <- c("x", "y", "label")
  n_groups <- length(unique(obj@active.ident))
  #--------------------------------------------------
  # Positions of labels
  #--------------------------------------------------
  cluster_num <- as.integer(as.character(sort(unique(df$label))))
  cog <- c()  # Center of gravity
  for(i in cluster_num){
    cog <- rbind(cog, cbind(
      mean(df[which(df$label == i),]$x), mean(df[which(df$label == i),]$y)))
  }
  cog <- data.frame(x = cog[,1], y = cog[,2], label = cluster_num)
  #--------------------------------------------------
  # Color
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
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, colour="") +
    scale_colour_manual(values=my_colors) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
    geom_text(aes(x=cog[,1], y=cog[,2], label=sprintf("%d", cluster_num)), size=6) +
    guides(colour = guide_legend(override.aes = list(size=4)))

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_tsne_nReads_seurat <- function(obj, title, title_size, xlabel, ylabel){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df <- data.frame(
    x = obj@reductions[["tsne"]]@cell.embeddings[,1],
    y = obj@reductions[["tsne"]]@cell.embeddings[,2],
    nReads = log(obj@meta.data[["nCount_RNA"]])
  )
  stdev <- sd(df$nReads)
  df$nReads <- (df$nReads - mean(df$nReads)) / stdev
  colnames(df) <- c("x", "y", "nReads")
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, color="z-log") +
    scale_colour_gradientn(colours=c("gray90","blue")) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_Heatmap_seurat <- function(
  obj, topn_genes, method, show_nReads, title, name, show_rownames_nReads,
  default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- as.matrix(obj@assays[["RNA"]]@data)
  mat <- t(scale(t(mat)))
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
  col_labels <- obj@active.ident
  row_labels <- topn_genes$label
  set.seed(8)
  col_hc <- hclust(dist(t(mat)), method = method)
  #--------------------------------------------------
  # Heatmap
  # Note: the positions of dendrogram leaves are
  # slightly adjusted by the gaps between slices.
  #--------------------------------------------------
  ht_opt$message = FALSE
  if(default_color){
    #------------------------------
    # Option 1: ggplot's default
    #------------------------------
    tmp <- unique(sort(col_labels))
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
    tmp <- unique(sort(col_labels))
    my_colors <- rainbow(length(tmp))
    names(my_colors) <- tmp
  }
  col_ha <- HeatmapAnnotation(
    Label = col_labels, col = list(Label = my_colors), annotation_name_side = "left"
  )
  p <- Heatmap(
    mat,
    column_title=title, column_title_gp=gpar(fontsize=16, fontface="bold"),
    name=name,
    column_split=col_labels, column_gap=unit(1.5, "mm"),
    row_split=row_labels, row_gap=unit(1.5, "mm"),
    border=FALSE,
    show_row_names=TRUE, row_names_side="left", show_row_dend=FALSE,
    show_column_names=FALSE, column_dend_side="top",
    show_parent_dend_line=TRUE, top_annotation=col_ha
  )
  if(show_nReads){
    mtx <- t(obj@meta.data[["nCount_RNA"]])
    rownames(mtx) <- "nReads"
    p_add <- Heatmap(
      mtx,
      name="nReads",
      cluster_columns=col_hc,
      show_row_names=TRUE, row_names_side="left", show_row_dend=FALSE,
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
plot_manifold2d_seurat <- function(
  obj, batch, plot_type, title, title_size, xlabel, ylabel, default_color
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df <- data.frame(x = obj@reductions[[plot_type]]@cell.embeddings[,1],
                   y = obj@reductions[[plot_type]]@cell.embeddings[,2])
  if(batch){
    df$label <- obj@meta.data[["orig.ident"]]
  }else{
    df$label <- obj@meta.data[["seurat_clusters"]]
  }
  colnames(df) <- c("x", "y", "label")
  n_groups <- length(unique(df$label))
  #--------------------------------------------------
  # Color
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
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, colour="") +
    scale_colour_manual(values=my_colors) +
    theme_classic(base_size=20, base_family="Helvetica") +  
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
    guides(colour = guide_legend(override.aes = list(size=4)))

  return(p)
 
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_manifold2d_seurat <- function(
  obj, batch, plot_type, title, title_size, xlabel, ylabel, default_color
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df <- data.frame(x = obj@reductions[[plot_type]]@cell.embeddings[,1],
                   y = obj@reductions[[plot_type]]@cell.embeddings[,2])
  if(batch){
    df$label <- obj@meta.data[["orig.ident"]]
  }else{
    df$label <- obj@meta.data[["seurat_clusters"]]
  }
  colnames(df) <- c("x", "y", "label")
  n_groups <- length(unique(df$label))
  #--------------------------------------------------
  # Color
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
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  if(batch){
    size <- 0.25
  }else{
    size <- 0.25
  }
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=size, alpha=1) +
    labs(title=title, x=xlabel, y=ylabel, colour="") +
    scale_colour_manual(values=my_colors) +
    theme_classic(base_size=20, base_family="Helvetica") +  
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
    guides(colour = guide_legend(override.aes = list(size=4)))

  return(p)
 
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_integrated_geneExpression_seurat <- function(
  obj, gene_name, plot_type, title, title_size, xlabel, ylabel, label_name
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df <- data.frame(x = obj@reductions[[plot_type]]@cell.embeddings[,1],
                   y = obj@reductions[[plot_type]]@cell.embeddings[,2])
  mat <- as.matrix(obj@assays[["integrated"]]@scale.data)
  mat <- mat[which(rownames(mat) == gene_name),]
  df$expr <- mat
  colnames(df) <- c("x", "y", "expr")
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, color=label_name) +
    scale_colour_gradientn(colours=c("blue", "white", "red")) +
    theme_classic(base_size=20, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_cc_umap_seurat <- function(
  obj, title, title_size, xlabel, ylabel, default_color = TRUE
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df <- data.frame(
    x = obj@reductions[["umap"]]@cell.embeddings[,1],
    y = obj@reductions[["umap"]]@cell.embeddings[,2],
    label = obj@meta.data[["Phase"]]
  )
  colnames(df) <- c("x", "y", "label")
  n_groups <- length(unique(df$label))
  #--------------------------------------------------
  # Color
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
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
#  my_colors <- c("#0000FF", "#00FF00", "#FF0000", "#8000FF") 
#  my_colors <- c("#0000FF", "#00FF00", "#FF0000") 
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,3]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, colour="") +
    scale_colour_manual(values=my_colors) +
    theme_classic(base_size=20, base_family="Helvetica") +  
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right") +
    guides(colour = guide_legend(override.aes = list(size=4)))

  return(p)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_ccscore_umap_seurat <- function(obj, title, title_size, xlabel, ylabel){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  df <- data.frame(
    x = obj@reductions[["umap"]]@cell.embeddings[,1],
    y = obj@reductions[["umap"]]@cell.embeddings[,2],
    s = obj@meta.data[["S.Score"]],
    t = obj@meta.data[["G2M.Score"]]
  )
  df$s <- (df$s - median(df$s)) / (IQR(df$s) / 1.3489)
  df$t <- (df$t - median(df$t)) / (IQR(df$t) / 1.3489)
  df$z <- df$s + df$t
  df$z <- (df$z - min(df$z)) / (max(df$z) - min(df$z))
  #--------------------------------------------------
  # Plot
  #--------------------------------------------------
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2], color=df[,5]), size=0.5, alpha=1.0) +
    labs(title=title, x=xlabel, y=ylabel, colour="") +
    scale_colour_gradientn(colours=c("blue", "green", "red")) +
    theme_classic(base_size=20, base_family="Helvetica") +  
    theme(plot.title=element_text(hjust=0.5, size=title_size),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")

  return(p)
}

