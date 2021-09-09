#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_quickQC_sign <- function(obj, data_type, min_ngenes = 2, max_ngenes = 1000){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  sign <- "sign"
  category_names <- names(obj[[sign]][[data_type]])
  obj_genes <- obj[["variable"]][["symbol"]]
  obj_geneIDs <- obj[["variable"]][["entrez"]]
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- paste("do_quickQC_", data_type, sep = "")
  obj[["history"]][[slot_name]][["min_ngenes"]] <- min_ngenes
  obj[["history"]][[slot_name]][["max_ngenes"]] <- max_ngenes
  #--------------------------------------------------
  # Select genes existing in user's data set
  #--------------------------------------------------
  for(category in category_names){
    tmp <- obj[[sign]][[data_type]][[category]]
    if(is.na(tmp$Gene[1]) & is.na(tmp$Gene[2]) & is.na(tmp$Gene[3]) &
       is.na(tmp$Gene[4]) & is.na(tmp$Gene[5]) & is.na(tmp$Gene[6])){
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
        for(g in geneIDs){
          genes <- c(genes, obj_genes[which(obj_geneIDs == g)])
        }
        tmp$Gene[i] <- paste(genes, collapse = "/")
      }
    }else if(
      is.na(tmp$GeneID[1]) & is.na(tmp$GeneID[2]) & is.na(tmp$GeneID[3]) &
      is.na(tmp$GeneID[4]) & is.na(tmp$GeneID[5]) & is.na(tmp$GeneID[6])
    ){
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
  }
  #--------------------------------------------------
  # Filter out the IDs including too few or too many genes
  #--------------------------------------------------
  for(category in category_names){
    tmp <- obj[[sign]][[data_type]][[category]]
    res <- tmp[which((tmp$Count >= min_ngenes) & (tmp$Count <= max_ngenes)),]
    rownames(res) <- 1:nrow(res)
    obj[[sign]][[data_type]][[category]] <- as.data.frame(res)
  }

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
separate_variables_sign <- function(
  obj, obj_cor, data_type, method, th_posi, th_nega
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  sign <- "sign"
  category_names <- names(obj[[sign]][[data_type]])
  cormat <- obj_cor[[method]]
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- paste("separate_variables_", data_type, sep = "")
  obj[["history"]][[slot_name]][["method"]] <- method
  obj[["history"]][[slot_name]][["th_posi"]] <- th_posi
  obj[["history"]][[slot_name]][["th_nega"]] <- th_nega
  #--------------------------------------------------
  # Separate genes
  #--------------------------------------------------
  for(category in category_names){
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
      # Definition
      #--------------------------------------------------
      genes <- unlist(strsplit(df$Gene[i], "/"))
      if(length(genes) <= 1){
        next
      }
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
        # Definition
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
        # Definition
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
        # Definition
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
        Count_strg = as.integer(
          ifelse(is.element(NA, genes_strg), 0, length(genes_strg))),
        Count_vari = as.integer(
          ifelse(is.element(NA, genes_vari), 0, length(genes_vari))),
        Count_weak = as.integer(
          ifelse(is.element(NA, genes_weak), 0, length(genes_weak))),
        Corr_strg = mean_strg,
        Corr_vari = mean_vari,
        Corr_weak = mean_weak,
        Gene_strg = ifelse(
          is.element(NA, genes_strg), NA, paste(genes_strg, collapse = "/")),
        Gene_vari = ifelse(
          is.element(NA, genes_vari), NA, paste(genes_vari, collapse = "/")),
        Gene_weak = ifelse(
          is.element(NA, genes_weak), NA, paste(genes_weak, collapse = "/")),
        GeneID = df$GeneID[i]
      ))
    }
    #--------------------------------------------------
    # Store the results
    #--------------------------------------------------
    obj[[sign]][[data_type]][[category]] <- res
  }

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
select_sign <- function(obj, data_type, min_cnt, min_cnt_weak = 2){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  sign <- "sign"
  category_names <- names(obj[[sign]][[data_type]])
  slot_name <- paste("select_", data_type, sep = "")
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  obj[["history"]][[slot_name]][["min_cnt"]] <- min_cnt
  obj[["history"]][[slot_name]][["min_cnt_weak"]] <- min_cnt_weak
  #--------------------------------------------------
  # 
  #--------------------------------------------------
  for(category in category_names){
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
  }

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
compute_SemSim_sign <- function(obj, data_type, measure, orgdb, treeTable = NULL, IC){
  #------------------------------
  # Definition
  #------------------------------
  sign <- "sign"
  slot_name <- paste(data_type, "_Sim", sep = "")
  category_names <- names(obj[[sign]][[data_type]])
  #------------------------------
  # History
  #------------------------------
  text <- paste("compute_SemSim_", data_type, sep = "")
  obj[["history"]][[text]][["measure"]] <- measure
  #------------------------------
  # Compute semantic similarities
  #------------------------------
  for(category in category_names){
    #------------------------------
    # Add IC
    #------------------------------
    df <- obj[[sign]][[data_type]][[category]]
    tmp <- IC[[category]]
    if(nrow(df) == 0){
      next
    }
    for(i in 1:nrow(df)){
      ind <- which(tmp$ID == df$ID[i])
      if(length(ind) == 0){
        df$IC[i] <- Inf
      }else{
        df$IC[i] <- tmp[ind,]$IC
      }
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
      #------------------------------
      # mgoSim
      #------------------------------
      simdata <- godata(OrgDb = orgdb, ont = category, computeIC = TRUE)
      simmat <- mgoSim(df$ID, df$ID, semData = simdata, measure = measure, combine = NULL)
      obj[[sign]][[slot_name]][[category]][["matrix"]] <- simmat  # Note: class(simmat) is "matrix"
    }
  }

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
reduce_sign <- function(obj, data_type, threshold, keep_rareID){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  sign <- "sign"
  slot_name <- paste(data_type, "_Sim", sep = "")
  category_names <- names(obj[[sign]][[data_type]])
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  text <- paste("reduce_", data_type, sep = "")
  obj[["history"]][[text]][["threshold"]] <- threshold
  obj[["history"]][[text]][["keep_rareID"]] <- keep_rareID
  #--------------------------------------------------
  # Remove IDs
  #--------------------------------------------------
  for(category in category_names){
    #--------------------------------------------------
    # Definition of report
    #--------------------------------------------------
    df <- obj[[sign]][[data_type]][[category]]
    if(nrow(df) == 0){
      next
    }
    simmat <- obj[[sign]][[slot_name]][[category]][["matrix"]]
    I <- dim(simmat)[1]
    if(I < 2){
      next
    }
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
    # Store the results
    #--------------------------------------------------
    obj[[sign]][[slot_name]][[category]][["report"]] <- report

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
  }

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
manual_exclusion_sign <- function(obj, data_type, keywords){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  sign <- "sign"
  category_names <- names(obj[[sign]][[data_type]])
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- paste("manual_exclusion_", data_type, sep = "")
  obj[["history"]][[slot_name]][["keywords"]] <- keywords
  #--------------------------------------------------
  # manual curation
  #--------------------------------------------------
  for(category in category_names){
    df <- obj[[sign]][[data_type]][[category]]
    if(nrow(df) != 0){
      df <- df[!((grepl(keywords, df$ID)) | (grepl(keywords, df$Description, ignore.case=TRUE))),]
      rownames(df) <- 1:nrow(df)
    }
    obj[[sign]][[data_type]][[category]] <- df
  }

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
manual_selection_sign <- function(obj, data_type, keywords){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  sign <- "sign"
  category_names <- names(obj[[sign]][[data_type]])
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- paste("manual_selection_", data_type, sep = "")
  obj[["history"]][[slot_name]][["keywords"]] <- keywords
  #--------------------------------------------------
  # manual curation
  #--------------------------------------------------
  for(category in category_names){
    df <- obj[[sign]][[data_type]][[category]]
    if(nrow(df) != 0){
      df <- df[(grepl(keywords, df$ID)) | (grepl(keywords, df$Description, ignore.case=TRUE)),]
      rownames(df) <- 1:nrow(df)
    }
    obj[[sign]][[data_type]][[category]] <- df
  }

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
make_signxsample_matrix <- function(
  obj, data_type, weight_strg = 0.5, weight_vari = 0.5
){
  #--------------------------------------------------
  # Error
  #--------------------------------------------------
  if((weight_strg < 0) | (weight_strg > 1) | (weight_vari < 0) |
     (weight_vari > 1)){
    stop("Error: weight_* must be from 0 to 1.")
  }
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  sign <- "sign"
  category_names <- names(obj[[sign]][[data_type]])
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  text <- paste("make_signxsample_matrix_", data_type, sep = "")
  obj[["history"]][[text]][["weight_strg"]] <- weight_strg
  obj[["history"]][[text]][["weight_vari"]] <- weight_vari
  #--------------------------------------------------
  # Sign-by-sample matrices
  #--------------------------------------------------
  mat <- obj[["data"]][["centered"]]
  for(category in category_names){
    res <- c()
    df <- obj[[sign]][[data_type]][[category]]
    if(nrow(df) != 0){
      for(i in 1:nrow(df)){
        #------------------------------
        # Strongly correlated gene set
        #------------------------------
        if(df$Count_strg[i] != 0){
          genes_strg <- unlist(strsplit(df$Gene_strg[i], "/"))
        }else{
          genes_strg <- NA
        }
        #------------------------------
        # Variably correlated gene set
        #------------------------------
        if(df$Count_vari[i] != 0){
          genes_vari <- unlist(strsplit(df$Gene_vari[i], "/"))
        }else{
          genes_vari <- NA
        }
        #------------------------------
        # Weakly correlated gene set
        #------------------------------
        if(df$Count_weak[i] != 0){
          genes_weak <- unlist(strsplit(df$Gene_weak[i], "/"))
        }else{
          genes_weak <- NA
        }
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
  }

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_pca_sign <- function(obj, data_type, category){
  #--------------------------------------------------
  # Definition
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
do_tsne_sign <- function(obj, data_type, category, pca_dim, tsne_dim){
  #--------------------------------------------------
  # Definition
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
do_umap_sign <- function(obj, data_type, category, pca_dim, umap_dim){
  #--------------------------------------------------
  # Definition
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
do_dmap_sign <- function(
  obj, data_type, category, sigma = "local", distance = "euclidean", pca_dim
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  algo_name <- "dmap"
  sign <- "sign"
  text <- paste(data_type, "xSample", sep = "")
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  obj[["history"]][["dmap"]][[data_type]][[category]][["sigma"]] <- sigma
  obj[["history"]][["dmap"]][[data_type]][[category]][["distance"]] <- distance
  #--------------------------------------------------
  # DiffusionMap
  #--------------------------------------------------
  set.seed(8)
  if(is.null(pca_dim)){
    mat <- as.matrix(obj[[sign]][[text]][[category]])
    res <- DiffusionMap(t(mat), sigma = sigma, distance = distance)
  }else{
    mat <- obj[["reduction"]][["pca"]][[data_type]][[category]][["x"]]
    res <- DiffusionMap(mat[,1:pca_dim], sigma = sigma, distance = distance)
  }
  #--------------------------------------------------
  # Store the results
  #--------------------------------------------------
  obj[["reduction"]][[algo_name]][[data_type]][[category]] <- res

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
cluster_samples_pam_sign <- function(obj, data_type, category, k){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  algo_name <- "pam"
  sign <- "sign"
  text <- paste(data_type, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[text]][[category]])
  slot_name <- paste(algo_name, "_", data_type, "_", category, sep = "")
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  text <- paste("cluster_samples_", algo_name, "_", data_type, "_", category, sep = "")
  obj[["history"]][[text]][["k"]] <- k
  #--------------------------------------------------
  # PAM
  #--------------------------------------------------
  set.seed(8)
  res <- pam(x = t(mat), k = k)
  #--------------------------------------------------
  # Store the results
  #--------------------------------------------------
  obj[["clustering"]][[algo_name]][[data_type]][[category]] <- res
  obj[["sample"]][[slot_name]] <- as.integer(res$cluster)

  return(obj)
} 
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
cluster_samples_hclustCutree_sign <- function(obj, data_type, category, method, k){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  algo_name <- "hclustCutree"
  sign <- "sign"
  text <- paste(data_type, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[text]][[category]])
  slot_name <- paste(algo_name, "_", data_type, "_", category, sep = "")
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  text <- paste("cluster_samples_", algo_name, "_", data_type, "_", category, sep = "")
  obj[["history"]][[text]][["method"]] <- method
  obj[["history"]][[text]][["k"]] <- k
  #--------------------------------------------------
  # hclust and cutree
  #--------------------------------------------------
  set.seed(8)
  hc <- hclust(dist(t(mat)), method = method)
  res <- cutree(hc, k = k)
  #--------------------------------------------------
  # Store the results
  #--------------------------------------------------
  obj[["clustering"]][[algo_name]][[data_type]][[category]][["hclust"]] <- hc
  obj[["clustering"]][[algo_name]][[data_type]][[category]][["cutree"]] <- res
  obj[["sample"]][[slot_name]] <- as.integer(res)

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
cluster_samples_seuratFindClusters_sign <- function(
  obj, data_type, category,
  reduction = "pca", dims, k.param, prune.SNN, resolution, algorithm
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  algo_name <- "seuratFindClusters"
  sign <- "sign"
  text <- paste(data_type, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[text]][[category]])
  slot_name <- paste(algo_name, "_", data_type, "_", category, sep = "")
  #--------------------------------------------------
  # Error control
  #--------------------------------------------------
  if(length(dims) > nrow(obj[[sign]][[text]][[category]])){
    stop("\"dims\" must be less than the number of rows of the SSM.")
  }
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  text <- paste("cluster_samples_", algo_name, "_", data_type, "_", category, sep = "")
  obj[["history"]][[text]][["reduction"]] <- reduction
  obj[["history"]][[text]][["dims"]] <- dims
  obj[["history"]][[text]][["k.param"]] <- k.param
  obj[["history"]][[text]][["prune.SNN"]] <- prune.SNN
  obj[["history"]][[text]][["resolution"]] <- resolution
  obj[["history"]][[text]][["algorithm"]] <- algorithm
  #--------------------------------------------------
  # Set Seurat object
  #--------------------------------------------------
  res <- CreateSeuratObject(
    counts = mat, project = obj[["history"]][["make_asurat_obj"]][["obj_name"]])
  res <- ScaleData(res, features = rownames(res))
  res <- RunPCA(res, features = rownames(res))
  #--------------------------------------------------
  # FindNeighbors
  #--------------------------------------------------
  set.seed(8)
  res <- FindNeighbors(
    res, reduction = reduction, dims = dims, k.param = k.param, prune.SNN = prune.SNN)
  #--------------------------------------------------
  # FindClusters
  #--------------------------------------------------
  res <- FindClusters(res, resolution = resolution, algorithm = algorithm)
  #--------------------------------------------------
  # Change the names of clusters by adding 1
  # e.g., change from `0` to `1`
  #--------------------------------------------------
  SetIdent(res, value = "seurat_clusters")
  num <- as.integer(as.character(Idents(res))) + 1
  label <- factor(num, levels = 1:max(num))
  #--------------------------------------------------
  # Store the results
  #--------------------------------------------------
  obj[["clustering"]][[algo_name]][[data_type]][[category]] <- label
  obj[["sample"]][[slot_name]] <- label

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_CalculateScaffoldTree_dmap <- function(
  obj, data_type, category, dims, NEndpoints, random_seed
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  algo_name <- "merlot"
  #--------------------------------------------------
  # CalculateScaffoldTree
  #--------------------------------------------------
  tmp <- obj[["reduction"]][["dmap"]][[data_type]][[category]]
  tmp <- CalculateScaffoldTree(
    CellCoordinates = tmp@eigenvectors[,dims],
    NEndpoints = NEndpoints, random_seed
  )
  obj[["clustering"]][[algo_name]][[data_type]][[category]][["ScaffoldTree"]] <- tmp

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_CalculateElasticTree_dmap <- function(obj, data_type, category, N_yk){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  algo_name <- "merlot"
  #--------------------------------------------------
  # CalculateElasticTree
  #--------------------------------------------------
  tmp <- obj[["clustering"]][[algo_name]][[data_type]][[category]]
  tmp <- CalculateElasticTree(
    ScaffoldTree = tmp[["ScaffoldTree"]],
    N_yk = N_yk
  )
  obj[["clustering"]][[algo_name]][[data_type]][[category]][["ElasticTree"]] <- tmp

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_SignSpaceEmbedding_dmap <- function(obj, data_type, category){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  algo_name <- "merlot"
  slot_name <- paste(data_type, "xSample", sep = "")
  mat <- t(obj[["sign"]][[slot_name]][[category]])
  ept <- obj[["clustering"]][[algo_name]][[data_type]][[category]][["ElasticTree"]]
  #--------------------------------------------------
  # GenesSpaceEmbedding
  #--------------------------------------------------
  tmp <- GenesSpaceEmbedding(ExpressionMatrix = mat, ElasticTree = ept)
  obj[["clustering"]][[algo_name]][[data_type]][[category]][["EmbeddedTree"]] <- tmp

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_CalculatePseudotimes_dmap <- function(obj, data_type, category, T0){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  algo_name <- "merlot"
  #--------------------------------------------------
  # CalculatePseudotimes
  #--------------------------------------------------
  tmp <- CalculatePseudotimes(
    InputTree = obj[["clustering"]][[algo_name]][[data_type]][[category]][["EmbeddedTree"]],
    T0 = T0
  )
  obj[["clustering"]][[algo_name]][[data_type]][[category]][["Pseudotimes"]] <- tmp

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
cluster_samples_merlot_sign <- function(obj, data_type, category){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  algo_name <- "merlot"
  sign <- "sign"
  text <- paste(data_type, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[text]][[category]])
  slot_name <- paste(algo_name, "_", data_type, "_", category, sep = "")
  #--------------------------------------------------
  # Store the results
  #--------------------------------------------------
  res <- obj[["clustering"]][[algo_name]][[data_type]][[category]][["Pseudotimes"]][["Cells2Branches"]]
  obj[["sample"]][[slot_name]] <- as.integer(res)

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
compute_nswap_001 <-
"int compute_nswap_001(NumericVector enco){
  int  swaps = 0;

  for(int i=0; i<(enco.length()-1); ++i){
    for(int j=i+1; j<enco.length(); ++j){
      if(enco[i] > enco[j])
        swaps++;
    }
  }
  return(swaps);
}"
Rcpp::cppFunction(compute_nswap_001)
#-----------------------------------------------------------------------------80
# This function computes the number of swaps of adjacent elements required for
# transforming vector 1 to vector 2, which have the same length and same amount
# of given numbers (e.g., vector1 = c(0,2,1), vector2 = c(1,2,0)).
# The following R code is inspired by Kendall tau distance algorithm.
#-----------------------------------------------------------------------------80
compute_nswaps <- function(vec1, vec2){
  I <- length(vec1)

  # (1)
  dict <- c()
  for(w in unique(vec1)){
    val <- which(vec1 == w)
    dict <- c(dict, list(val))
  }
  names(dict) <- unique(vec1)

  # (2)
  cnt <- c()
  for(c in names(dict)){
    cnt <- c(cnt, list(1))
  }
  names(cnt) <- names(dict)
  enco <- vec2
  for(i in 1:I){
    enco[i] <- dict[[as.character(vec2[i])]][cnt[[as.character(vec2[i])]]]
    cnt[[as.character(vec2[i])]] <- cnt[[as.character(vec2[i])]] + 1
  }

  # (3): since this step is time consuming, we utilize `002_codes/utilities.cpp`
  swaps <- compute_nswap_001(enco)

  return(swaps)
}
#-----------------------------------------------------------------------------80
#
# This is find_marker_sign, not find_marker_gene !!!!!
#
#-----------------------------------------------------------------------------80
find_marker_sign <- function(
  obj, data_type_for_label, category_for_label, algo_name_for_label,
  data_type_for_expr, category_for_expr, labels_1, labels_2
){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  sign <- "sign"
  #--------------------------------------------------
  # Divide samples
  #--------------------------------------------------
  slot_name <- paste(algo_name_for_label, "_", data_type_for_label, "_",
    category_for_label, sep = "")
  tmp <- obj[["sample"]][[slot_name]]
  inds_1 <- which(tmp %in% labels_1)
  inds_2 <- which(tmp %in% labels_2)
  pop_1 <- obj[["sample"]][["barcode"]][inds_1]
  pop_2 <- obj[["sample"]][["barcode"]][inds_2]
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  slot_name <- paste(data_type_for_expr, "xSample", sep = "")
  mat <- as.matrix(obj[[sign]][[slot_name]][[category_for_expr]])
  mat <- mat[,union(inds_1, inds_2)]
  tmp <- obj[[sign]][[data_type_for_expr]][[category_for_expr]]
  #--------------------------------------------------
  # Prepare a data frame to be output
  #--------------------------------------------------
  res <- data.frame(
    Label = as.factor(paste(labels_1, collapse = "/")),
    Rank = rep(NA, dim(mat)[1]),
    Sign = NA,
    ParentDescription = NA,
    Gene = NA,
    sep_I = NA
  )
  for(i in 1:dim(mat)[1]){
    #--------------------------------------------------
    # Sign
    #--------------------------------------------------
    res$Sign[i] <- rownames(mat)[i]
    #--------------------------------------------------
    # ParentDescription
    #--------------------------------------------------
    res$ParentDescription[i] <-
      tmp[which(tmp$ID == strsplit(res$Sign[i], split='_')[[1]][1]),]$Description
    #--------------------------------------------------
    # Gene
    #--------------------------------------------------
    if(strsplit(res$Sign[i], split='_')[[1]][2] == "S"){
      res$Gene[i] <- tmp[which(tmp$ID == strsplit(res$Sign[i], split='_')[[1]][1]),]$Gene_strg
    }else{
      res$Gene[i] <- tmp[which(tmp$ID == strsplit(res$Sign[i], split='_')[[1]][1]),]$Gene_vari
    }
    #--------------------------------------------------
    # Edit distance between two (0,1)-vectors
    #--------------------------------------------------
    data <- mat[i,]
    vec_1 <- sort(data, decreasing = FALSE)
    vec_1 <- ifelse(names(vec_1) %in% pop_1, 1, 0)
    vec_2 <- sort(vec_1, decreasing = FALSE) # vec_3 = (0, 0, ..., 1, 1)
    vec_3 <- sort(vec_1, decreasing = TRUE)  # vec_2 = (1, 1, ..., 0, 0)
    dist1 <- compute_nswaps(vec_1, vec_2)
    dist2 <- compute_nswaps(vec_1, vec_3)
    res$sep_I[i] <- round(1 - 2 * dist1 / (dist1 + dist2), digits = 6)
  }
  #--------------------------------------------------
  # Arrange the data frame in order of res$sep_I
  #--------------------------------------------------
  inds <- order(res$sep_I, decreasing = TRUE)
  res <- res[inds,]
  res$Rank <- 1:length(inds)
  #--------------------------------------------------
  # Output
  #--------------------------------------------------
  slot_name <- paste(
    "Label_", data_type_for_label, "_", category_for_label, "_",
    paste(labels_1, collapse = "/"), "_vs_", paste(labels_2, collapse = "/"),
    sep = ""
  )
  obj[["marker"]][[data_type_for_expr]][[category_for_expr]][[slot_name]] <- res
  return(obj)
}
#-----------------------------------------------------------------------------80
#
# This is find_marker_gene, not find_marker_sign !!!!!
#
#-----------------------------------------------------------------------------80
find_marker_gene <- function(
  obj, data_type_for_label, category_for_label, algo_name_for_label,
  labels_1, labels_2, parametric
){
  #--------------------------------------------------
  # Divide samples
  #--------------------------------------------------
  slot_name <- paste(algo_name_for_label, "_", data_type_for_label, "_",
    category_for_label, sep = "")
  tmp <- obj[["sample"]][[slot_name]]
  inds_1 <- which(tmp %in% labels_1)
  inds_2 <- which(tmp %in% labels_2)
  pop_1 <- obj[["sample"]][["barcode"]][inds_1]
  pop_2 <- obj[["sample"]][["barcode"]][inds_2]
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  mat <- as.matrix(obj[["data"]][["normalized"]])
  #--------------------------------------------------
  # Error control
  #--------------------------------------------------
  if((length(inds_1) <= 1) || (length(inds_2) <= 1)){
    stop("Error: one of the subpopulations is too small.")
  }
  #--------------------------------------------------
  # Prepare a data frame to be output
  #--------------------------------------------------
  res <- data.frame(
    Label = as.factor(paste(labels_1, collapse = "/")),
    Rank = rep(NA, dim(mat)[1]),
    Gene = NA,
    p_val = NA,
    p_adj = NA
  )
  if(parametric){
    for(i in 1:dim(mat)[1]){
      #--------------------------------------------------
      # Gene
      #--------------------------------------------------
      res$Gene[i] <- rownames(mat)[i]
      #--------------------------------------------------
      # p_val
      #--------------------------------------------------
      values_1 <- mat[, inds_1][i,]
      values_2 <- mat[, inds_2][i,]
      res$p_val[i] <- t.test(x = values_1, y = values_2,
        var.equal = FALSE, paired = FALSE)[["p.value"]]
    }
  }else{
    for(i in 1:dim(mat)[1]){
      #--------------------------------------------------
      # Gene
      #--------------------------------------------------
      res$Gene[i] <- rownames(mat)[i]
      #--------------------------------------------------
      # p_val
      #--------------------------------------------------
      values_1 <- mat[, inds_1][i,]
      values_2 <- mat[, inds_2][i,]
      res$p_val[i] <- wilcox.exact(x = values_1, y = values_2,
        paired = FALSE)[["p.value"]]
    }
  }
  #--------------------------------------------------
  # p_adj
  #--------------------------------------------------
  res$p_adj <- p.adjust(res$p_val, method = "BH")
  #--------------------------------------------------
  # Arrange the data frame in order of res$p_val
  #--------------------------------------------------
  inds <- order(res$p_val, decreasing = FALSE)
  res <- res[inds,]
  res$Rank <- 1:length(inds)
  #--------------------------------------------------
  # Output
  #--------------------------------------------------
  slot_name <- paste(
    "Label_", data_type_for_label, "_", category_for_label, "_",
    paste(labels_1, collapse = "/"), "_vs_", paste(labels_2, collapse = "/"),
    sep = ""
  )
  obj[["marker"]][["gene"]][[data_type_for_label]][[category_for_label]][[slot_name]] <- res
  return(obj)
}
#-----------------------------------------------------------------------------80
#
# This is auto_find_marker_sign, not auto_find_marker_gene !!!!!
#
#-----------------------------------------------------------------------------80
auto_find_marker_sign <- function(obj, data_type_for_label, category_for_label, 
  algo_name_for_label, data_type_for_expr, category_for_expr
){
  #--------------------------------------------------
  # Initialize
  #--------------------------------------------------
  obj[["marker"]][[data_type_for_expr]][[category_for_expr]] <- NULL
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  slot_name <- paste(algo_name_for_label, "_", data_type_for_label, "_",
    category_for_label, sep = "")
  labels <- unique(sort(obj[["sample"]][[slot_name]]))
  #--------------------------------------------------
  # Looping
  #--------------------------------------------------
  res <- c()
  for(i in labels){
    labels_1 <- i
    labels_2 <- setdiff(labels, i)
    tmp <- find_marker_sign(
      obj = obj, data_type_for_label = data_type_for_label,
      category_for_label = category_for_label,
      algo_name_for_label = algo_name_for_label,
      data_type_for_expr = data_type_for_expr,
      category_for_expr = category_for_expr,
      labels_1 = labels_1, labels_2 = labels_2)
    slot_name <- paste(
      "Label_", data_type_for_label, "_", category_for_label, "_",
      paste(labels_1, collapse = "/"), "_vs_", paste(labels_2, collapse = "/"), sep = ""
    )
    res[[slot_name]] <- tmp[["marker"]][[data_type_for_expr]][[category_for_expr]][[slot_name]]
  }
  res[["all"]] <- c()
  for(name in names(res)){
    res[["all"]] <- rbind(res[["all"]], res[[name]])
  }
  rownames(res[["all"]]) <- 1:nrow(res[["all"]])
  #--------------------------------------------------
  # Output
  #--------------------------------------------------
  obj[["marker"]][[data_type_for_expr]][[category_for_expr]] <- res

  return(obj)
}
#-----------------------------------------------------------------------------80
#
# This is auto_find_marker_gene, not auto_find_marker_sign !!!!!
#
#-----------------------------------------------------------------------------80
auto_find_marker_gene <- function(
  obj, data_type_for_label, category_for_label, algo_name_for_label, parametric
){
  #--------------------------------------------------
  # Initialize
  #--------------------------------------------------
  obj[["marker"]][["gene"]][[data_type_for_label]][[category_for_label]] <- NULL
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  slot_name <- paste(
    algo_name_for_label, "_", data_type_for_label, "_", category_for_label,
    sep = ""
  )
  labels <- unique(sort(obj[["sample"]][[slot_name]]))
  #--------------------------------------------------
  # Looping
  #--------------------------------------------------
  res <- c()
  for(i in labels){
    labels_1 <- i
    labels_2 <- setdiff(labels, i)
    tmp <- find_marker_gene(
      obj = obj, data_type_for_label = data_type_for_label,
      category_for_label = category_for_label,
      algo_name_for_label = algo_name_for_label,
      labels_1 = labels_1, labels_2 = labels_2,
      parametric = parametric
    )
    slot_name <- paste(
      "Label_", data_type_for_label, "_", category_for_label, "_",
      paste(labels_1, collapse = "/"), "_vs_", paste(labels_2, collapse = "/"), sep = ""
    )
    res[[slot_name]] <-
      tmp[["marker"]][["gene"]][[data_type_for_label]][[category_for_label]][[slot_name]]
  }
  res[["all"]] <- c()
  for(name in names(res)){
    res[["all"]] <- rbind(res[["all"]], res[[name]])
  }
  rownames(res[["all"]]) <- 1:nrow(res[["all"]])
  #--------------------------------------------------
  # Output
  #--------------------------------------------------
  obj[["marker"]][["gene"]][[data_type_for_label]][[category_for_label]] <- res

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
concatenate_obj_sign <- function(obj1, obj2){
  #--------------------------------------------------
  # history
  #--------------------------------------------------
  hist_names_1 <- names(obj1[["history"]])
  hist_names_2 <- names(obj2[["history"]])
  for(names in setdiff(hist_names_2, hist_names_1)){
    obj1[["history"]][[names]] <- obj2[["history"]][[names]]
  }
  #--------------------------------------------------
  # sample
  #--------------------------------------------------
  samp_names_1 <- names(obj1[["sample"]])
  samp_names_2 <- names(obj2[["sample"]])
  for(name in setdiff(samp_names_2, samp_names_1)){
    obj1[["sample"]][[name]] <- obj2[["sample"]][[name]]
  }
  #--------------------------------------------------
  # misc
  #--------------------------------------------------
  misc_names_1 <- names(obj1[["misc"]])
  misc_names_2 <- names(obj2[["misc"]])
  for(name in setdiff(misc_names_2, misc_names_1)){
    obj1[["misc"]][[name]] <- obj2[["misc"]][[name]]
  }
  #--------------------------------------------------
  # sign
  #--------------------------------------------------
  sign_names_1 <- names(obj1[["sign"]])
  sign_names_2 <- names(obj2[["sign"]])
  for(name in setdiff(sign_names_2, sign_names_1)){
    obj1[["sign"]][[name]] <- obj2[["sign"]][[name]]
  }
  #--------------------------------------------------
  # clustering
  #--------------------------------------------------
  clustering_names_1 <- names(obj1[["clustering"]])
  clustering_names_2 <- names(obj2[["clustering"]])
  clustering_names <- union(clustering_names_1, clustering_names_2)
  for(name in clustering_names){
    if(length(obj2[["clustering"]][[name]]) != 0){
      category_names <- names(obj2[["clustering"]][[name]])
      for(cname in category_names){
        obj1[["clustering"]][[name]][[cname]] <-
          obj2[["clustering"]][[name]][[cname]]
      }
    }
  }
  #--------------------------------------------------
  # reduction
  #--------------------------------------------------
  reduction_names_1 <- names(obj1[["reduction"]])
  reduction_names_2 <- names(obj2[["reduction"]])
  reduction_names <- union(reduction_names_1, reduction_names_2)
  for(name in reduction_names){
    if(length(obj2[["reduction"]][[name]]) != 0){
      category_names <- names(obj2[["reduction"]][[name]])
      for(cname in category_names){
        obj1[["reduction"]][[name]][[cname]] <- obj2[["reduction"]][[name]][[cname]]
      }
    }
  }
  #--------------------------------------------------
  # marker
  #--------------------------------------------------
  marker_names_1 <- names(obj1[["marker"]])
  marker_names_2 <- names(obj2[["marker"]])
  marker_names <- union(marker_names_1, marker_names_2)
  for(name in marker_names){
    if(length(obj2[["marker"]][[name]]) != 0){
      category_names <- names(obj2[["marker"]][[name]])
      for(cname in category_names){
        obj1[["marker"]][[name]][[cname]] <- obj2[["marker"]][[name]][[cname]]
      }
    }
  }

  return(obj1)
}
