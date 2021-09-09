#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_groupGO <- function(genes, orgdb, ont, level){
  res <- groupGO(
    gene = genes,
    OrgDb = orgdb,
    keyType = "ENTREZID",
    ont = ont,
    level = level,
    readable = FALSE
  )
  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
collect_GO <- function(orgdb, all_geneIDs){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  category_names <- c("MF", "BP", "CC")

  #--------------------------------------------------
  # groupGO
  #--------------------------------------------------
  res <- list()
  for(category in category_names){
    tmp <- c()
    level <- 1
    while(1){
      ggo <- try(
        do_groupGO(gene = all_geneIDs, orgdb = orgdb, ont = category, level = level),
        silent = TRUE
      )
      if(class(ggo) == "try-error"){
        break
      }else{
        tmp <- c(tmp, list(ggo))
        level <- level + 1
      }
    }
    res[[category]] <- tmp
  }

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
format_GO <- function(dict){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  category_names <- names(dict)

  #--------------------------------------------------
  # Reformat
  #--------------------------------------------------
  res <- list()
  for(category in category_names){
    tmp <- c()
    level_names <- length(dict[[category]])
    for(level in 1:level_names){
      tmp <- rbind(tmp, as.data.frame(dict[[category]][[level]]@result))
    }
    tmp <- unique(tmp)  # Notice: the right hand side must be data.frame
    df <- data.frame(
      ID = tmp$ID,
      Description = tmp$Description,
      Count = tmp$Count,
      Gene = NA,
      GeneID = tmp$geneID
    )
    df$Count <- as.numeric(df$Count)
    res[[category]] <- df
  }

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
compute_IC_GO <- function(tidy){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  category_names <- names(tidy)
  #--------------------------------------------------
  # Compute information contents (ICs)
  #--------------------------------------------------
  IC_GO <- list()
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
    if(category == "MF") tmp <- toTable(GOMFOFFSPRING)
    if(category == "BP") tmp <- toTable(GOBPOFFSPRING)
    if(category == "CC") tmp <- toTable(GOCCOFFSPRING)
    names(tmp) <- c("child", "parent")
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
    if(category == "MF"){
      ind <- which(res$parent == "GO:0003674")
      if(length(ind) == 0)
        stop("Error: tidy_GO must include root ontology term.")
      res$Freq_root <- res[ind,]$Freq_parent
    }else if(category == "BP"){
      ind <- which(res$parent == "GO:0008150")
      if(length(ind) == 0)
        stop("Error: tidy_GO must include root ontology term.")
      res$Freq_root <- res[ind,]$Freq_parent
    }else if(category == "CC"){
      ind <- which(res$parent == "GO:0005575")
      if(length(ind) == 0)
        stop("Error: tidy_GO must include root ontology term.")
      res$Freq_root <- res[ind,]$Freq_parent
    }
    res$Prob <- res$Freq_parent / res$Freq_root
    res$IC <- -log(res$Prob)
    IC_GO[[category]] <- data.frame(ID = res$parent, IC = res$IC)
  }

  return(IC_GO)
}
