#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_keggGet <- function(data){
  return(keggGet(data))
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
collect_KEGG <- function(organism, category_names){
  #--------------------------------------------------
  # collect KEGG
  #--------------------------------------------------
  res <- list()
  for(category in category_names){
    ids <- keggLink(category, organism)  # Note: names(ids) are gene IDs
    map <- data.frame(
      ID = ids,
      Description = NA,
      KEGG_geneID = names(ids)
    )
    map$NCBI_geneID <- NA
    flags <- c()
    I <- length(map$ID)
    for(i in 1:I){
      Sys.sleep(0.1)
      #--------------------------------------------------
      # Print the current process
      #--------------------------------------------------
      if((i == 1) || (i %% floor(0.25 * I) == 0 && i < 0.95 * I) || (i == I)){
        text <- paste("Now processing ", i, "/", I, " for ", category, "...\n", sep = "")
        cat(text)
      }
      #--------------------------------------------------
      # (i) Description
      #--------------------------------------------------
      category_id <- map$ID[i]
      if(category == "module")
        category_id <- str_extract(category_id, "(?<=_)(.*)")
      tmp_description <- try(do_keggGet(category_id), silent = TRUE)
      if(class(tmp_description) == "try-error"){
        flags <- c(flags, i)
        next
      }
      #--------------------------------------------------
      # (ii) NCBI_geneID
      #--------------------------------------------------
      tmp_gene <- try(do_keggGet(map$KEGG_geneID[i]), silent = TRUE)
      if(class(tmp_gene) == "try-error"){
        flags <- c(flags, i)
        next
      }
      #------------------------------
      # (i)
      #------------------------------
      name <- tmp_description[[1]][["NAME"]]
      map$Description[i] <- ifelse(is.null(name), "No_record", name)
      #------------------------------
      # (ii)
      #------------------------------
      dblinks <- tmp_gene[[1]][["DBLINKS"]]
      id <- dblinks[grep("NCBI-GeneID", dblinks)]
      id <- gsub("NCBI-GeneID: ", "", id)
      map$NCBI_geneID[i] <- id
    }
    res[[category]][["success"]] <- map
    #--------------------------------------------------
    # Rescue the failures
    #--------------------------------------------------
    if(length(flags) > 0){
      failure <- c()
      for(i in flags){
        #------------------------------
        # (i) Description
        #------------------------------
        category_id <- map$ID[i]
        if(category == "module")
          category_id <- str_extract(category_id, "(?<=_)(.*)")
        tmp_description <- try(do_keggGet(category_id), silent = TRUE)
        cnt = 0
        while(class(tmp_description) == "try-error"){
          cnt <- cnt + 1
          if(cnt >= 10)
            break
          Sys.sleep(1)
          tmp_description <- try(do_keggGet(category_id), silent = TRUE)
        }
        if(cnt >= 10){
          failure <- rbind(failure, c("ID", map$ID[i]))
          next
        }
        #------------------------------
        # (ii) NCBI_geneID
        #------------------------------
        tmp_gene <- try(do_keggGet(map$KEGG_geneID[i]), silent = TRUE)
        cnt = 0
        while(class(tmp_gene) == "try-error"){
          cnt <- cnt + 1
          if(cnt >= 10)
            break
          Sys.sleep(1)
          tmp_gene <- try(do_keggGet(map$KEGG_geneID[i]), silent = TRUE)
        }
        if(cnt >= 10){
          failure <- rbind(failure, c("KEGG_geneID", map$KEGG_geneID[i]))
          next
        }
        #----------
        # (i)
        #----------
        name <- tmp_description[[1]][["NAME"]]
        map$Description[i] <- ifelse(is.null(name), "No_record", name)
        #----------
        # (ii)
        #----------
        dblinks <- tmp_gene[[1]][["DBLINKS"]]
        id <- dblinks[grep("NCBI-GeneID", dblinks)]
        id <- gsub("NCBI-GeneID: ", "", id)
        map$NCBI_geneID[i] <- id
      }
      res[[category]][["success"]] <- map
      tmp <- data.frame(
        matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c(
          "Input_category_for_keggGet", "Input_for_keggGet"
        ))))
      tmp <- rbind(tmp, data.frame(
        Input_category_for_keggGet = unique(failure)[,1],
        Input_for_keggGet = unique(failure)[,2]
      ))
      res[[category]][["failure"]] <- tmp
    }
  }

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
format_KEGG <- function(dict){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
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
      genes <- unique(tmp[which(tmp[["ID"]] == map[["ID"]][i]),]$NCBI_geneID)
      map$GeneID[i] <- paste(genes, collapse = "/")
      map$Count[i] <- length(genes)
    }
    rownames(map) <- 1:nrow(map)
    res[[category]] <- map
  }

  return(res)
}
