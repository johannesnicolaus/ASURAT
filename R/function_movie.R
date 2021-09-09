#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_quickQC_sign_movie <- function(
  obj, data_type, category, min_ngenes = 2, max_ngenes = 1000
){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  sign <- "sign"
  obj_genes <- obj[["variable"]][["symbol"]]
  obj_geneIDs <- obj[["variable"]][["entrez"]]
  #--------------------------------------------------
  # Select genes existing in user's data set
  #--------------------------------------------------
  tmp <- obj[[sign]][[data_type]][[category]]
  if(is.na(tmp$Gene[1]) & is.na(tmp$Gene[2]) & is.na(tmp$Gene[3])){
    for(i in 1:nrow(tmp)){
      #------------------------------
      # GeneID: NCBI gene (formerly Entrezgene) ID
      #------------------------------
      geneIDs <- unlist(strsplit(tmp$GeneID[i], "/"))
      geneIDs <- geneIDs[which(geneIDs %in% intersect(geneIDs, obj_geneIDs))]
      tmp$GeneID[i] <- paste(geneIDs, collapse = "/")
      tmp$Count[i] <- as.integer(length(geneIDs))
      #------------------------------
      # Gene: gene symbol
      #------------------------------
      genes <- c()
      for(g in geneIDs)
        genes <- c(genes, obj_genes[which(obj_geneIDs == g)])
      tmp$Gene[i] <- paste(genes, collapse = "/")
    }
  }else if(is.na(tmp$GeneID[1]) & is.na(tmp$GeneID[2]) & is.na(tmp$GeneID[3])){
    for(i in 1:nrow(tmp)){
      #------------------------------
      # Gene: gene symbol
      #------------------------------
      genes <- unlist(strsplit(tmp$Gene[i], "/"))
      genes <- genes[which(genes %in% intersect(genes, obj_genes))]
      tmp$Gene[i] <- paste(genes, collapse = "/")
      tmp$Count[i] <- as.integer(length(genes))
      #------------------------------
      # GeneID: NCBI gene (formerly Entrezgene) ID
      #------------------------------
      geneIDs <- c()
      for(g in genes){
        geneIDs <- c(geneIDs, obj_geneIDs[which(obj_genes == g)])
      }
      tmp$GeneID[i] <- paste(geneIDs, collapse = "/")
    }
  }else{
    for(i in 1:nrow(tmp)){
      #------------------------------
      # Count
      #------------------------------
      genes <- unlist(strsplit(tmp$Gene[i], "/"))
      genes <- genes[which(genes %in% intersect(genes, obj_genes))]
      tmp$Count[i] <- as.integer(length(genes))
    }
  }
  obj[[sign]][[data_type]][[category]] <- tmp
  #--------------------------------------------------
  # Filter out the IDs including too few or too many genes
  #--------------------------------------------------
  tmp <- obj[[sign]][[data_type]][[category]]
  res <- tmp[which((tmp$Count >= min_ngenes) & (tmp$Count <= max_ngenes)),]
  rownames(res) <- 1:nrow(res)
  obj[[sign]][[data_type]][[category]] <- as.data.frame(res)

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
separate_variables_sign_movie <- function(
  obj, obj_cor, data_type, category, method, th_posi, th_nega
){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  sign <- "sign"
  cormat <- obj_cor[[method]]
  #--------------------------------------------------
  # Separate genes
  #--------------------------------------------------
  df <- obj[[sign]][[data_type]][[category]]
  res <- data.frame(
    matrix(ncol = 12, nrow = 0, dimnames = list(NULL, c(
      "ID",
      "Description",
      "Count_strg",
      "Count_vari",
      "Count_weak",
      "Corr_strg",
      "Corr_vari",
      "Corr_weak",
      "Gene_strg",
      "Gene_vari",
      "Gene_weak",
      "GeneID"
    ))))
  for(i in 1:nrow(df)){
    #--------------------------------------------------
    # Definitions
    #--------------------------------------------------
    genes <- unlist(strsplit(df$Gene[i], "/"))
    if(length(genes) <= 1)
      next
    inds <- which(rownames(cormat) %in% genes)
    mat <- as.matrix(cormat[inds, inds])
    #--------------------------------------------------
    # Decomposition of correlation matrix
    #--------------------------------------------------
    diag(mat) <- -99
    inds_posi <- which(apply(mat, 2, function(x) sum(x >= th_posi)) > 0)
    diag(mat) <- 99
    inds_nega <- which(apply(mat, 2, function(x) sum(x <= th_nega)) > 0)
    #--------------------------------------------------
    #
    #--------------------------------------------------
    if(length(inds_posi) == 0){
      genes_strg <- NA
      genes_vari <- NA
      genes_weak <- genes
      mean_strg <- NA
      mean_vari <- NA
      inds <- which(rownames(mat) %in% genes_weak)
      mtx <- mat[inds, inds]
      mean_weak <- ifelse(length(inds) <= 1, NA, mean(mtx[mtx <= 1]))
    }else if(length(inds_nega) >= 2){
      #------------------------------
      # Definitions
      #------------------------------
      inds <- union(inds_posi, inds_nega)
      tmp <- as.matrix(mat[inds, inds])
      gset <- list()
      mean <- list()
      #------------------------------
      # Pam clustering
      #------------------------------
      diag(tmp) <- 1
      set.seed(8)
      clust <- pam(tmp, k = 2)
      diag(tmp) <- 99
      for(j in 1:2){
        gset[[j]] <- names(clust[["clustering"]])[which(clust[["clustering"]] == j)]
        inds <- which(rownames(tmp) %in% gset[[j]])
        mtx <- tmp[inds, inds]
        mean[[j]] <- ifelse(length(inds) <= 1, NA, mean(mtx[mtx <= 1]))
      }
      #------------------------------
      # Definitions
      #------------------------------
      means <- unlist(mean)
      inds <- order(means, decreasing = TRUE)
      genes_strg <- gset[[inds[1]]]
      genes_vari <- gset[[inds[2]]]
      mean_strg <- means[inds[1]]
      mean_vari <- means[inds[2]]
      genes_weak <- setdiff(genes, union(gset[[1]], gset[[2]]))
      inds <- which(rownames(mat) %in% genes_weak)
      mtx <- mat[inds, inds]
      mean_weak <- ifelse(length(inds) <= 1, NA, mean(mtx[mtx <= 1]))
    }else{
      mtx <- mat[inds_posi, inds_posi]
      #------------------------------
      # Definitions
      #------------------------------
      genes_strg <- rownames(mtx)
      genes_vari <- NA
      genes_weak <- setdiff(genes, genes_strg)
      mean_strg <- mean(mtx[mtx <= 1])
      mean_vari <- NA
      inds <- which(rownames(mat) %in% genes_weak)
      mtx <- mat[inds, inds]
      mean_weak <- ifelse(length(inds) <= 1, NA, mean(mtx[mtx <= 1]))
    }
    #--------------------------------------------------
    # Re-separation
    #--------------------------------------------------
    if((!(is.na(mean_strg))) && (mean_strg < th_posi)){
      genes_strg <- NA
      genes_vari <- NA
      genes_weak <- genes
      mean_strg <- NA
      mean_vari <- NA
      inds <- which(rownames(mat) %in% genes_weak)
      mtx <- mat[inds, inds]
      mean_weak <- ifelse(length(inds) <= 1, NA, mean(mtx[mtx <= 1]))
    }
    #--------------------------------------------------
    # res
    #--------------------------------------------------
    res <- rbind(res, data.frame(
      ID = df$ID[i],
      Description = df$Description[i],
      Count_strg = as.integer(ifelse(is.element(NA, genes_strg), 0, length(genes_strg))),
      Count_vari = as.integer(ifelse(is.element(NA, genes_vari), 0, length(genes_vari))),
      Count_weak = as.integer(ifelse(is.element(NA, genes_weak), 0, length(genes_weak))),
      Corr_strg = mean_strg,
      Corr_vari = mean_vari,
      Corr_weak = mean_weak,
      Gene_strg = ifelse(is.element(NA, genes_strg), NA, paste(genes_strg, collapse = "/")),
      Gene_vari = ifelse(is.element(NA, genes_vari), NA, paste(genes_vari, collapse = "/")),
      Gene_weak = ifelse(is.element(NA, genes_weak), NA, paste(genes_weak, collapse = "/")),
      GeneID = df$GeneID[i]
    ))
  }
  #--------------------------------------------------
   # Store the results
  #--------------------------------------------------
  obj[[sign]][[data_type]][[category]] <- res

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
select_sign_movie <- function(obj, data_type, category, min_cnt, min_cnt_weak = 2){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  sign <- "sign"
  slot_name <- paste("select_", data_type, sep = "")
  #--------------------------------------------------
  # 
  #--------------------------------------------------
  df <- obj[[sign]][[data_type]][[category]]
  #============================================================
  # User-defined criteria:
  #   `|` and `&` mean "or" and "and", respectively.
  #============================================================
  # (Start) ===================================================

  inds <- which(
    (df$Count_strg + df$Count_vari >= min_cnt) &
    df$Count_weak >= min_cnt_weak
  )

  # (End) =====================================================
  #============================================================
  df <- df[inds,]
  if(length(inds) != 0)
    rownames(df) <- 1:nrow(df)
  obj[[sign]][[data_type]][[category]] <- df

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
compute_SemSim_sign_movie <- function(
  obj, data_type, category, measure, orgdb, treeTable = NULL, IC
){
  #------------------------------
  # Definitions
  #------------------------------
  sign <- "sign"
  slot_name <- paste(data_type, "_Sim", sep = "")
  #------------------------------
  # Compute semantic similarities
  #------------------------------
  #------------------------------
  # Add IC
  #------------------------------
  df <- obj[[sign]][[data_type]][[category]]
  tmp <- IC[[category]]
  for(i in 1:nrow(df)){
    ind <- which(tmp$ID == df$ID[i])
    if(length(ind) == 0)
      df$IC[i] <- -log(0)
    else
      df$IC[i] <- tmp[ind,]$IC
  }
  obj[[sign]][[data_type]][[category]] <- df
  #------------------------------
  # Reorder: necessary for the next step
  #------------------------------
  df <- df[order(df$IC, decreasing = FALSE),]

  #==================================================
  # DO
  #==================================================
  if(data_type == "DO"){
    library(DOSE)
    #------------------------------
    # doSim
    #------------------------------
    simmat <- doSim(df$ID, df$ID, measure = measure)
    obj[[sign]][[slot_name]][[category]][["matrix"]] <- simmat  # Note: class(simmat) is "matrix"
  }
  #==================================================
  # CO
  #==================================================
  if(data_type == "CO"){
    simmat <- matrix(0, nrow = nrow(df), ncol = nrow(df))
    tree <- treeTable[[category]]
    for(i in 1:(nrow(df)-1)){
      for(j in (i+1):nrow(df)){
        #------------------------------
        # Information content (IC)
        #------------------------------
        IC_i <- tmp[which(tmp$ID == df$ID[i]),]$IC
        IC_j <- tmp[which(tmp$ID == df$ID[j]),]$IC
        #------------------------------
        # IC of most informative common ancestor (MICA)
        #------------------------------
        ancestors_i <- tree[which(tree$child == df$ID[i]),]$parent
        ancestors_j <- tree[which(tree$child == df$ID[j]),]$parent
        common_ancestors <- intersect(ancestors_i, ancestors_j)
        common_ancestors <- data.frame(ID = common_ancestors, IC = NA)
        for(n in 1:nrow(common_ancestors))
          common_ancestors$IC[n] <- tmp[which(tmp$ID == common_ancestors$ID[n]),]$IC
        IC_MICA <- common_ancestors[order(common_ancestors$IC, decreasing = TRUE),]$IC[1]
        #------------------------------
        # Resnik, Lin
        #------------------------------
        if(measure == "Resnik")
          simmat[i,j] <- IC_MICA
        else if(measure == "Lin")
          simmat[i,j] <- 2 * IC_MICA / (IC_i + IC_j)
        else
          stop("Error: check the mesure.")
      }
    }
    simmat[lower.tri(simmat)] <- simmat[upper.tri(simmat)]
    diag(simmat) <- 1
    rownames(simmat) <- df$ID
    colnames(simmat) <- df$ID
    obj[[sign]][[slot_name]][[category]][["matrix"]] <- simmat
  }
  #==================================================
  # GO
  #==================================================
  if(data_type == "GO"){
    library(GO.db)
    library(GOSemSim)
    #------------------------------
    # mgoSim
    #------------------------------
    simdata <- godata(OrgDb = orgdb, ont = category, computeIC = TRUE)
    simmat <- mgoSim(df$ID, df$ID, semData = simdata, measure = measure, combine = NULL)
    obj[[sign]][[slot_name]][[category]][["matrix"]] <- simmat  # Note: class(simmat) is "matrix"
  }

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
reduce_sign_movie <- function(
  obj, data_type, category, threshold, keep_rareID
){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  sign <- "sign"
  slot_name <- paste(data_type, "_Sim", sep = "")
  category_names <- names(obj[[sign]][[data_type]])
  #--------------------------------------------------
  # Remove IDs
  #--------------------------------------------------
  #--------------------------------------------------
  # Definition of report
  #--------------------------------------------------
  df <- obj[[sign]][[data_type]][[category]]
  simmat <- obj[[sign]][[slot_name]][[category]][["matrix"]]
  I <- dim(simmat)[1]
  if(I < 2)
     stop("Error: similarity matrix is too small")
  report <- data.frame(
    matrix(ncol = 11, nrow = 0, dimnames = list(NULL, c(
      "Similarity",
      "Kept_ID",
      "Removed_ID",
      "Kept_Description",
      "Removed_Description",
      "Kept_Count",
      "Removed_Count",
      "Kept_Gene",
      "Removed_Gene",
      "Kept_IC",
      "Removed_IC"
    ))))
  flag <- rep(0, dim(simmat)[1])
  #--------------------------------------------------
  #
  #--------------------------------------------------
  if(keep_rareID == TRUE){
    for(i in I:2){
      J <- dim(simmat)[2] - (I - i + 1)  # Must be set before the next command
      if(flag[i] == 1)                   # Must be set after the previous command
        next
      for(j in 1:J){
        if(is.na(simmat[i,j]))
          next
        if(abs(simmat[i,j]) >= threshold){
          ind_kept <- which(df$ID == rownames(simmat)[i])
          ind_removed <- which(df$ID == rownames(simmat)[j])
          report <- rbind(report, data.frame(
            Similarity = simmat[i,j],
            Kept_ID = df[ind_kept,]$ID,
            Removed_ID = df[ind_removed,]$ID,
            Kept_Description = df[ind_kept,]$Description,
            Removed_Description = df[ind_removed,]$Description,
            Kept_Count = df[ind_kept,]$Count_strg +
              df[ind_kept,]$Count_vari + df[ind_kept,]$Count_weak,
            Removed_Count = df[ind_removed,]$Count_strg +
              df[ind_removed,]$Count_vari + df[ind_removed,]$Count_weak,
            Kept_Gene = paste(df[ind_kept,]$Gene_strg, df[ind_kept,]$Gene_vari,
              df[ind_kept,]$Gene_weak, sep = "||"),
            Removed_Gene = paste(df[ind_removed,]$Gene_strg, df[ind_removed,]$Gene_vari,
              df[ind_removed,]$Gene_weak, sep = "||"),
            Kept_IC = df[ind_kept,]$IC,
              Removed_IC = df[ind_removed,]$IC
          ))
          flag[j] <- 1
        }
      }
    }
  }else{
    I <- dim(simmat)[1] - 1
    J <- dim(simmat)[2]
    for(i in 1:I){
      j0 <- i + 1       # Must be set before the next command
      if(flag[i] == 1)  # Must be set after the previous command
        next
      for(j in j0:J){
        if(is.na(simmat[i,j]))
          next
        if(abs(simmat[i,j]) >= threshold){
          ind_kept <- which(df$ID == rownames(simmat)[i])
          ind_removed <- which(df$ID == rownames(simmat)[j])
          report <- rbind(report, data.frame(
            Similarity = simmat[i,j],
            Kept_ID = df[ind_kept,]$ID,
            Removed_ID = df[ind_removed,]$ID,
            Kept_Description = df[ind_kept,]$Description,
            Removed_Description = df[ind_removed,]$Description,
            Kept_Count = df[ind_kept,]$Count_strg +
              df[ind_kept,]$Count_vari + df[ind_kept,]$Count_weak,
            Removed_Count = df[ind_removed,]$Count_strg +
              df[ind_removed,]$Count_vari + df[ind_removed,]$Count_weak,
            Kept_Gene = paste(df[ind_kept,]$Gene_strg, df[ind_kept,]$Gene_vari,
              df[ind_kept,]$Gene_weak, sep = "||"),
            Removed_Gene = paste(df[ind_removed,]$Gene_strg, df[ind_removed,]$Gene_vari,
              df[ind_removed,]$Gene_weak, sep = "||"),
            Kept_IC = df[ind_kept,]$IC,
            Removed_IC = df[ind_removed,]$IC
          ))
          flag[j] <- 1
        }
      }
    }
  }
  #--------------------------------------------------
  # Remove IDs
  #--------------------------------------------------
  if(length(report$Removed_ID) >= 1){
    kept_IDs <- setdiff(df$ID, report$Removed_ID)
    res <- df[which(df$ID %in% kept_IDs),]
    rownames(res) <- 1:nrow(res)
  }else{
    res <- df
  }
  obj[[sign]][[data_type]][[category]] <- res

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
make_signxsample_matrix_movie <- function(
  obj, data_type, category, weight_strg = 0.5, weight_vari = 0.5
){
  #--------------------------------------------------
  # Error
  #--------------------------------------------------
  if((weight_strg < 0) | (weight_strg > 1) | (weight_vari < 0) | (weight_vari > 1))
    stop("Error: weight_* must be from 0 to 1.")
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  sign <- "sign"
  category_names <- names(obj[[sign]][[data_type]])
  #--------------------------------------------------
  # Sign-by-sample matrices
  #--------------------------------------------------
  mat <- obj[["data"]][["centered"]]
  res <- c()
  df <- obj[[sign]][[data_type]][[category]]
  if(nrow(df) != 0){
    for(i in 1:nrow(df)){
      #------------------------------
      # Strongly correlated gene set
      #------------------------------
      if(df$Count_strg[i] != 0)
        genes_strg <- unlist(strsplit(df$Gene_strg[i], "/"))
      else
        genes_strg <- NA
      #------------------------------
      # Variably correlated gene set
      #------------------------------
      if(df$Count_vari[i] != 0)
        genes_vari <- unlist(strsplit(df$Gene_vari[i], "/"))
      else
        genes_vari <- NA
      #------------------------------
      # Weakly correlated gene set
      #------------------------------
      if(df$Count_weak[i] != 0)
        genes_weak <- unlist(strsplit(df$Gene_weak[i], "/"))
      else
        genes_weak <- NA
      #------------------------------
      # Sign-by-sample matrix
      #------------------------------
      if((df$Count_strg[i] >= 2) & (df$Count_vari[i] < 2) & (df$Count_weak[i] >= 2)){
        tmp_strg <- as.data.frame(mat[which(rownames(mat) %in% genes_strg),])
        tmp_weak <- as.data.frame(mat[which(rownames(mat) %in% genes_weak),])
        vec_strg <- weight_strg * apply(tmp_strg, 2, mean) +
          (1 - weight_strg) * apply(tmp_weak, 2, mean)
        tmp <- as.matrix(t(vec_strg))
        rownames(tmp) <- paste(df$ID[i], "_", "S", sep = "")
      }
      if((df$Count_strg[i] >= 2) & (df$Count_vari[i] >= 2) & (df$Count_weak[i] >= 2)){
        tmp_strg <- as.data.frame(mat[which(rownames(mat) %in% genes_strg),])
        tmp_vari <- as.data.frame(mat[which(rownames(mat) %in% genes_vari),])
        tmp_weak <- as.data.frame(mat[which(rownames(mat) %in% genes_weak),])
        vec_strg <- weight_strg * apply(tmp_strg, 2, mean) +
          (1 - weight_strg) * apply(tmp_weak, 2, mean)
        vec_vari <- weight_vari * apply(tmp_vari, 2, mean) +
          (1 - weight_vari) * apply(tmp_weak, 2, mean)
        tmp <- rbind(vec_strg, vec_vari)
        rownames(tmp) <- c(paste(df$ID[i], "_", "S", sep = ""), paste(df$ID[i], "_", "V", sep = ""))
      }
      res <- rbind(res, tmp)
    }
  }
  slot_name <- paste(data_type, "xSample", sep = "")
  obj[[sign]][[slot_name]][[category]] <- as.data.frame(res)

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_pca_sign_movie <- function(obj, data_type, category, pca_dim, tsne_dim){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  sign <- "sign"
  slot_name <- paste(data_type, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[slot_name]][[category]])
  #--------------------------------------------------
  # PCA
  #--------------------------------------------------
  set.seed(8)
  obj[["reduction"]][["pca"]][[data_type]][[category]] <- prcomp(t(mat), scale = FALSE)

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_tsne_sign_movie <- function(obj, data_type, category, pca_dim, tsne_dim){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  sign <- "sign"
  slot_name <- paste(data_type, "xSample", sep = "")
  #--------------------------------------------------
  # tSNE
  #--------------------------------------------------
  set.seed(8)
  if(is.null(pca_dim)){
    mat <- as.matrix(obj[[sign]][[slot_name]][[category]])
    obj[["reduction"]][["tsne"]][[data_type]][[category]] <-
      Rtsne(t(mat), dim = tsne_dim, pca = FALSE)
  }else{
    mat <- obj[["reduction"]][["pca"]][[data_type]][[category]][["x"]]
    obj[["reduction"]][["tsne"]][[data_type]][[category]] <-
      Rtsne(mat[,1:pca_dim], dim = tsne_dim, pca = FALSE)
  }

  return(obj)
} 
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_umap_sign_movie <- function(obj, data_type, category, pca_dim, umap_dim){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  sign <- "sign"
  slot_name <- paste(data_type, "xSample", sep = "")
  #--------------------------------------------------
  # UMAP
  #--------------------------------------------------
  set.seed(8)
  if(is.null(pca_dim)){
    mat <- as.matrix(obj[[sign]][[slot_name]][[category]])
    obj[["reduction"]][["umap"]][[data_type]][[category]] <-
      umap(t(mat), n_components = umap_dim)
  }else{
    mat <- obj[["reduction"]][["pca"]][[data_type]][[category]][["x"]]
    obj[["reduction"]][["umap"]][[data_type]][[category]] <-
      umap(mat[,1:pca_dim], n_components = umap_dim)
  }

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_DiffusionMap_movie <- function(
  obj, data_type, category, sigma = "local", distance = "euclidean"
){
  #--------------------------------------------------
  # Definitions
  #--------------------------------------------------
  algo_name <- "dmap"
  sign <- "sign"
  text <- paste(data_type, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[text]][[category]])
  #--------------------------------------------------
  # DiffusionMap
  #--------------------------------------------------
  set.seed(8)
  res <- DiffusionMap(t(mat), sigma = sigma, distance = distance)
  #--------------------------------------------------
  # Store the results
  #--------------------------------------------------
  obj[["reduction"]][[algo_name]][[data_type]][[category]] <- res

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_pca_sign_movie <- function(
  obj, data_type, category, p1, p2, p3, title, xlabel, ylabel, dirname, cnt
){
  nsigns <- nrow(obj[["sign"]][[data_type]][[category]])
  mat <- obj[["reduction"]][["pca"]][[data_type]][[category]][["x"]]
  filename <- paste(dirname, "/figure_", cnt, ".png", sep="")
  df <- data.frame(x = mat[,1], y = mat[,2])
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2]), color="black", size=0.5, alpha=0.4) +
#    xlim(-50, 50) + ylim(-50, 50) +
    annotate("text", x=-Inf, y=Inf, label=paste("th_p = ", format(p1, nsmall=2), sep=""),
      hjust=-0.16, vjust=1.5, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("th_n = ", format(p2, nsmall=2), sep=""),
      hjust=-0.16, vjust=3.2, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("min_cnt = ", p3, sep=""),
      hjust=-0.15, vjust=4.9, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("nsigns = ", nsigns, sep=""),
      hjust=-0.15, vjust=6.6, size=4) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=15, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
  ggsave(file=filename, plot=p, dpi=150, width=3.8, height=4)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_tsne_sign_movie <- function(
  obj, data_type, category, p1, p2, p3, title, xlabel, ylabel, dirname, cnt
){
  nsigns <- nrow(obj[["sign"]][[data_type]][[category]])
  mat <- obj[["reduction"]][["tsne"]][[data_type]][[category]][["Y"]]
  filename <- paste(dirname, "/figure_", cnt, ".png", sep="")
  df <- data.frame(x = mat[,1], y = mat[,2])
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2]), color="black", size=0.5, alpha=0.4) +
#    xlim(-50, 50) + ylim(-50, 50) +
    annotate("text", x=-Inf, y=Inf, label=paste("th_p = ", format(p1, nsmall=2), sep=""),
      hjust=-0.16, vjust=1.5, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("th_n = ", format(p2, nsmall=2), sep=""),
      hjust=-0.16, vjust=3.2, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("min_cnt = ", p3, sep=""),
      hjust=-0.15, vjust=4.9, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("nsigns = ", nsigns, sep=""),
      hjust=-0.15, vjust=6.6, size=4) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=15, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
  ggsave(file=filename, plot=p, dpi=150, width=3.8, height=4)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_umap_sign_movie <- function(
  obj, data_type, category, p1, p2, p3, title, xlabel, ylabel, dirname, cnt
){
  nsigns <- nrow(obj[["sign"]][[data_type]][[category]])
  mat <- obj[["reduction"]][["umap"]][[data_type]][[category]][["layout"]]
  filename <- paste(dirname, "/figure_", cnt, ".png", sep="")
  df <- data.frame(x = mat[,1], y = mat[,2])
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2]), color="black", size=0.5, alpha=0.4) +
#    xlim(-15, 15) + ylim(-15, 15) +
    annotate("text", x=-Inf, y=Inf, label=paste("th_p = ", format(p1, nsmall=2), sep=""),
      hjust=-0.16, vjust=1.5, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("th_n = ", format(p2, nsmall=2), sep=""),
      hjust=-0.16, vjust=3.2, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("min_cnt = ", p3, sep=""),
      hjust=-0.15, vjust=4.9, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("nsigns = ", nsigns, sep=""),
      hjust=-0.15, vjust=6.6, size=4) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=15, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
  ggsave(file=filename, plot=p, dpi=150, width=3.8, height=4)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
plot_DiffusionMap_movie <- function(
  obj, data_type, category, p1, p2, p3, title, xlabel, ylabel, dirname, cnt
){
  nsigns <- nrow(obj[["sign"]][[data_type]][[category]])
  mat <- obj[["reduction"]][["dmap"]][[data_type]][[category]]@eigenvectors
  filename <- paste(dirname, "/figure_", cnt, ".png", sep="")
  df <- data.frame(x = mat[,1], y = mat[,2])
  p <- ggplot() +
    geom_point(aes(x=df[,1], y=df[,2]), color="black", size=0.5, alpha=0.4) +
#    xlim(-15, 15) + ylim(-15, 15) +
    annotate("text", x=-Inf, y=Inf, label=paste("th_p = ", format(p1, nsmall=2), sep=""),
      hjust=-0.16, vjust=1.5, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("th_n = ", format(p2, nsmall=2), sep=""),
      hjust=-0.16, vjust=3.2, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("min_cnt = ", p3, sep=""),
      hjust=-0.15, vjust=4.9, size=4) +
    annotate("text", x=-Inf, y=Inf, label=paste("nsigns = ", nsigns, sep=""),
      hjust=-0.15, vjust=6.6, size=4) +
    labs(title=title, x=xlabel, y=ylabel) +
    theme_classic(base_size=15, base_family="Helvetica") +
    theme(plot.title=element_text(hjust=0.5),
      plot.margin=unit(c(1.5,1.5,1.5,1.5), "lines"), legend.position="right")
  ggsave(file=filename, plot=p, dpi=150, width=3.8, height=4)
}
