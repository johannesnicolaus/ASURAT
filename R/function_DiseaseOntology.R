#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
collect_DO <- function(){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  data(DO2EG)
  #--------------------------------------------------
  # enrichDO
  #--------------------------------------------------
  res <- enrichDO(
    unlist(DO2EG),
    ont = "DO",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 0,
    maxGSSize = 10000000,
    qvalueCutoff = 1,
    readable = FALSE
  )

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
format_DO <- function(dict){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  dict <- list(disease = dict)
  category_names <- names(dict)
  #--------------------------------------------------
  # Reformat
  #--------------------------------------------------
  res <- list()
  for(category in category_names){
    tmp <- dict[[category]]@result
    map <- unique(data.frame(
      ID = tmp[["ID"]],
      Description = tmp[["Description"]],
      Count = tmp[["Count"]],
      Gene = NA,
      GeneID = tmp[["geneID"]]
    ))
    res[[category]] <- map
  }

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
compute_IC_DO <- function(tidy){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  category_names <- names(tidy)
  #--------------------------------------------------
  # Compute information contents (ICs)
  #--------------------------------------------------
  IC_DO <- list()
  for(category in category_names){
    #------------------------------
    # Definition
    #------------------------------
    df <- tidy[[category]]
    all_geneIDs <- c()
    for(i in 1:nrow(df))
      all_geneIDs <- c(all_geneIDs, unlist(strsplit(df$GeneID[i], "/")))
    all_geneIDs <- unique(all_geneIDs)
    #--------------------------------------------------
    # Prepare (count, child, parent)-table
    #--------------------------------------------------
    tmp <- toTable(DOOFFSPRING) ; names(tmp) <- c("child", "parent")
    tmp <- data.frame(num_annot_child = NA, child = tmp$child, parent = tmp$parent)
    tmp <- tmp[which((tmp$child %in% df$ID) & (tmp$parent %in% df$ID)),]
    rownames(tmp) <- 1:nrow(tmp)
    #--------------------------------------------------
    # Count the number of times a gene is annotated with the children
    #--------------------------------------------------
    for(i in 1:nrow(tmp)){
      cnt <- df[which(df$ID == tmp$child[i]),]$Count
      if(length(cnt) != 0)
        tmp$num_annot_child[i] <- cnt
      else
        tmp$num_annot_child[i] <- 0
    }
    #--------------------------------------------------
    # Add parents without child
    #--------------------------------------------------
    tmp <- rbind(tmp, data.frame(
      num_annot_child = 0,
      child = NA,
      parent = setdiff(df$ID, unique(tmp$parent))
    ))
    #--------------------------------------------------
    # Compute the sum of `tmp$num_annot_child` for each parent
    #--------------------------------------------------
    ids <- unique(tmp$parent)
    res <- data.frame(parent = ids, num_annot_parent = NA, sum_annot_child = NA)
    for(i in 1:length(ids)){
      cnt <- df[which(df$ID == ids[i]),]$Count
      if(length(cnt) != 0)
        res$num_annot_parent[i] <- cnt
      else
        res$num_annot_parent[i] <- 0
      res$sum_annot_child[i] <- sum(tmp[which(tmp$parent == ids[i]),]$num_annot_child)
    }
    #--------------------------------------------------
    # Compute IC for each parent
    #--------------------------------------------------
    res$Freq_parent <- res$num_annot_parent + res$sum_annot_child
    if(category == "disease"){
      ind <- which(res$parent == "DOID:4")
      if(length(ind) == 0)
        stop("Error: tidy_DO must include root ontology term.")
      res$Freq_root <- res[ind,]$Freq_parent
    }
    res$Prob <- res$Freq_parent / res$Freq_root
    res$IC <- -log(res$Prob)
    IC_DO[[category]] <- data.frame(ID = res$parent, IC = res$IC)
  }

  return(IC_DO)
}
