#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_cellTypeToGenes <- function(data, orgdb){
  res <- suppressMessages(cellTypeToGenes(data, orgDb = orgdb, gotab = allGOterms))
  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
collect_CO <- function(orgdb){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  co <- ontoProc::getCellOnto()
  co <- data.frame(ID = co[["id"]], Description = co[["name"]])
  co <- co[which(!is.na(co$Description)),]
  co <- co[(grepl("CL:", co$ID)),]
  #--------------------------------------------------
  # Collect COs
  #--------------------------------------------------
  res <- data.frame(
    matrix(ncol = 5, nrow = 0, dimnames = list(NULL, c(
      "ID",
      "Description",
      "Symbol",
      "GO",
      "Evidence"
    ))))
  for(i in 1:nrow(co)){
    tmp <- try(do_cellTypeToGenes(co$Description[i], orgdb = orgdb), silent = TRUE)
    if(class(tmp) == "try-error")
      next
    else if(nrow(tmp) != 0){
      for(j in 1:nrow(tmp)){
        res <- rbind(res, data.frame(
          ID = co$ID[i],
          Description = co$Description[i],
          Symbol = tmp$SYMBOL[j],
          GO = tmp$GO[j],
          Evidence = tmp$EVIDENCE[j]
        ))
      }
    }
  }

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
format_CO <- function(dict){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  dict <- list(cell = dict)
  category_names <- names(dict)
  #--------------------------------------------------
  # Reformat
  #--------------------------------------------------
  res <- list()
  for(category in category_names){
    tmp <- dict[[category]]
    map <- unique(data.frame(
      ID = tmp[["ID"]],
      Description = tmp[["Description"]],
      Count = NA,
      Gene = NA,
      GeneID = NA
    ))
    for(i in 1:nrow(map)){
      #------------------------------
      # Gene and Count
      #------------------------------
      genes <- unique(tmp[which(tmp[["ID"]] == map[["ID"]][i]),]$Symbol)
      map$Gene[i] <- paste(genes, collapse = "/")
      map$Count[i] <- length(genes)
    }
    rownames(map) <- 1:nrow(map)
    res[[category]] <- map
  }

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
find_descendants <- function(id, co, map){
  children <- co[["children"]][[id]]
  if(length(children) == 0){
    return(NA)
  }else{
    for(child in children)
      map <- c(map, c(child, find_descendants(child, co, map)))
    return(setdiff(map, NA))
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
make_treeTable_CO <- function(tidy){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  category_names <- names(tidy)
  category_names_woALL <- setdiff(category_names, "ALL")
  co <- ontoProc::getCellOnto()
  #--------------------------------------------------
  # Compute information contents (ICs)
  #--------------------------------------------------
  res <- list() 
  for(category in category_names_woALL){
    #------------------------------
    # Definition
    #------------------------------
    df <- tidy[[category]]

    #--------------------------------------------------
    # Prepare (count, child, parent)-table
    #--------------------------------------------------
    map <- c()
    tmp <- c()
    for(i in 1:nrow(df))
      tmp <- rbind(
        tmp,
        data.frame(child = find_descendants(df$ID[i], co, map), parent = df$ID[i])
      ) 
    tmp <- tmp[which((tmp$child %in% df$ID) & (tmp$parent %in% df$ID)),]
    rownames(tmp) <- 1:nrow(tmp)
    res[[category]] <- tmp
  }

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
compute_IC_CO <- function(tidy, treeTable){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  category_names <- names(treeTable)

  #--------------------------------------------------
  # Compute information contents (ICs)
  #--------------------------------------------------
  IC_CO <- list()
  for(category in category_names){
    #------------------------------
    # Definition
    #------------------------------
    df <- tidy[[category]]
    all_genes <- c()
    for(i in 1:nrow(df))
      all_genes <- c(all_genes, unlist(strsplit(df$Gene[i], "/")))
    all_genes <- unique(all_genes)
    tmp <- treeTable[[category]]
    tmp <- data.frame(num_annot_child = NA, child = tmp$child, parent = tmp$parent)
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
    if(category == "cell"){
      ind <- which(res$parent == "CL:0000000")
      if(length(ind) == 0)
        stop("Error: tidy_CO must include root ontology term.")
      res$Freq_root <- res[ind,]$Freq_parent
    }
    res$Prob <- res$Freq_parent / res$Freq_root
    res$IC <- -log(res$Prob)
    IC_CO[[category]] <- data.frame(ID = res$parent, IC = res$IC)
  }

  return(IC_CO)
}
