#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Reformat a Disease Ontology database.
#'
#' This function reformats a Disease Ontology database.
#'
#' @param dict A res slot of enrichDO().
#' @param all_geneIDs All genes in the form of ENTREZ ID.
#' @param orgdb A genome annotation package.
#'
#' @return A formatted database.
#' @export
#'
format_DO_asrt <- function(dict, all_geneIDs, orgdb){
  #--------------------------------------------------
  # Compute information contents.
  # See computeIC() in DOSE package.
  #--------------------------------------------------
  doids <- BiocGenerics::toTable(DO.db::DOTERM)
  doterms <- doids$do_id
  docount <- table(doterms)
  doids <- names(docount)  #unique(doterms)
  xx <- as.list(DO.db::DOOFFSPRING)
  cnt <- vapply(doids, function(x){
    n <- docount[xx[[x]]]
    docount[x] + sum(n[!is.na(n)])
  }, integer(1))
  names(cnt) <- doids
  p <- cnt / sum(docount)
  IC <- -log(p)
  dict$IC <- NA
  for(i in seq_len(nrow(dict))){
    if(is.element(dict$ID[i], names(IC))){
      dict$IC[i] <- IC[which(names(IC) == dict$ID[i])]
    }else{
      dict$IC[i] <- 99
    }
  }
  df <- data.frame(ID = dict$ID, Description = dict$Description, IC = dict$IC,
                   Count = dict$Count, Gene = NA, GeneID = dict$geneID)
  res <- list(disease = unique(df))
  #--------------------------------------------------
  # Fix the slots of gene symbols and ENTREZ Gene IDs.
  #--------------------------------------------------
  dictionary <- AnnotationDbi::select(orgdb, key = unique(all_geneIDs),
                                      columns = "SYMBOL",
                                      keytype = "ENTREZID")
  for(k in seq_len(length(res))){
    for(i in seq_len(nrow(res[[k]]))){
      genes <- c()
      geneIDs <- c()
      g <- unlist(strsplit(res[[k]]$GeneID[i], "/"))
      if(length(g) == 0){
        next
      }
      for(j in seq_len(length(g))){
        ind <- which(dictionary$ENTREZID == g[j])
        if(!is.na(dictionary[ind, ]$SYMBOL[1])){
          genes <- c(genes, dictionary[ind, ]$SYMBOL[1])
          geneIDs <- c(geneIDs, g[j])
        }
      }
      res[[k]]$Gene[i] <- paste(genes, collapse = "/")
      res[[k]]$GeneID[i] <- paste(geneIDs, collapse = "/")
      res[[k]]$Count[i] <- as.integer(length(geneIDs))
    }
  }
  tidy <- list()
  categories <- names(res)
  for(k in seq_len(length(categories))){
    tidy[[categories[k]]] <- res[[k]]
  }
  #--------------------------------------------------
  # Compute a similarity matrix.
  #--------------------------------------------------
  categories <- names(res)
  for(k in seq_len(length(categories))){
    df <- res[[k]]
    simmat <- DOSE::doSim(df$ID, df$ID, measure = "Jiang")
    tidy[["similarity_matrix"]][[categories[k]]] <- simmat
  }

  return(tidy)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Do cellTypeToGenes().
#'
#' This function performs cellTypeToGenes().
#'
#' @param data Description.
#' @param orgdb A genome annotation package.
#'
#' @return Results of cellTypeToGenes().
#'
do_cellTypeToGenes <- function(data, orgdb){
  res <- suppressMessages(
    ontoProc::cellTypeToGenes(data, orgDb = orgdb,
                              gotab = ontoProc::allGOterms))
  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Collect a Cell Ontology database.
#'
#' This function performs cellTypeToGenes() for collecting a Cell Ontology
#'   database.
#'
#' @param orgdb A genome annotation package.
#'
#' @return Results from getCellOnto().
#' @export
#'
collect_CO_asrt <- function(orgdb){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  co <- ontoProc::getCellOnto()
  co <- data.frame(ID = co[["id"]], Description = co[["name"]])
  co <- co[which(!is.na(co$Description)), ]
  co <- co[(grepl("CL:", co$ID)), ]
  #--------------------------------------------------
  # Collect CO terms.
  #--------------------------------------------------
  res <- c("ID", "Description", "Symbol", "GO", "Evidence")
  res <- data.frame(matrix(ncol = 5, nrow = 0, dimnames = list(NULL, res)))
  for(i in seq_len(nrow(co))){
    tmp <- try(do_cellTypeToGenes(co$Description[i], orgdb = orgdb),
               silent = TRUE)
    if(class(tmp) == "try-error"){
      next
    }else if(nrow(tmp) != 0){
      for(j in seq_len(nrow(tmp))){
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
#' Find all the descendants for parent terms.
#'
#' This function finds all the descendants for parent terms.
#'
#' @param id ID.
#' @param co A result of getCellOnto().
#' @param map A map data.
#'
#' @return A map data.
#'
find_descendants_asrt <- function(id, co, map){
  children <- co[["children"]][[id]]
  if(length(children) == 0){
    return(NA)
  }else{
    for(i in seq_len(length(children))){
      map <- c(map, c(children[i], find_descendants_asrt(children[i], co, map)))
    }
    return(setdiff(map, NA))
  }
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Make a parent-child relation table for CO terms.
#'
#' This function makes a parent-child relation table for CO terms.
#'
#' @param tidy A list of result of format_CO().
#'
#' @return Parent-child relation table.
#'
make_treeTable_CO_asrt <- function(tidy){
  categories_woALL <- setdiff(names(tidy), "ALL")
  co <- ontoProc::getCellOnto()
  res <- list() 
  for(k in seq_len(length(categories_woALL))){
    df <- tidy[[categories_woALL[k]]]
    map <- c() ; tmp <- c()
    for(i in seq_len(nrow(df))){
      dg <- data.frame(child = find_descendants_asrt(df$ID[i], co, map),
                       parent = df$ID[i])
      tmp <- rbind(tmp, dg) 
    }
    res[[categories_woALL[k]]] <- tmp
  }

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Compute information contents.
#'
#' This function computes information contents for Cell Ontology terms.
#'
#' @param dict A result of format_CO().
#' @param tidy A list of result of format_CO().
#' @param treeTable A result of make_treeTable_CO_asrt().
#'
#' @return A list of result of format_CO().
#' @seealso Mistry and Pavlidis, BMC Bioinformatics, 2008.
#'
compute_IC_CO_asrt <- function(dict, tidy, treeTable){
  categories_woALL <- setdiff(names(tidy), "ALL")
  for(k in seq_len(length(categories_woALL))){
    #------------------------------
    # Definition
    #------------------------------
    all_genes <- unique(dict$Symbol)
    res <- tidy[[categories_woALL[k]]]
    tmp <- data.frame(num_annot_child = NA,
                      child = treeTable[[categories_woALL[k]]]$child,
                      parent = treeTable[[categories_woALL[k]]]$parent)
    #--------------------------------------------------
    # Count the number of times that a gene is annotated with children.
    #--------------------------------------------------
    for(i in seq_len(nrow(tmp))){
      cnt <- res[which(res$ID == tmp$child[i]), ]$Count
      if(length(cnt) != 0){
        tmp$num_annot_child[i] <- cnt
      }else{
        tmp$num_annot_child[i] <- 0
      }
    }
    #--------------------------------------------------
    # Compute the sum of `tmp$num_annot_child` for each parent.
    #--------------------------------------------------
    ids <- unique(tmp$parent)
    tbl <- data.frame(parent = ids, num_annot_parent = NA, sum_annot_child = NA)
    for(i in seq_len(length(ids))){
      cnt <- res[which(res$ID == ids[i]), ]$Count
      if(length(cnt) != 0){
        tbl$num_annot_parent[i] <- cnt
      }else{
        tbl$num_annot_parent[i] <- 0
      }
      tbl$sum_annot_child[i] <-
        sum(tmp[which(tmp$parent == ids[i]), ]$num_annot_child)
    }
    #--------------------------------------------------
    # Compute IC for each parent.
    #--------------------------------------------------
    tbl$Freq_parent <- tbl$num_annot_parent + tbl$sum_annot_child
    if(categories_woALL[k] == "cell"){
      ind <- which(tbl$parent == "CL:0000000")
      if(length(ind) == 0){
        stop("tidy_CO must include root ontology term.")
      } 
      tbl$Freq_root <- tbl[ind, ]$Freq_parent
    }
    tbl$Prob <- tbl$Freq_parent / tbl$Freq_root
    tbl$IC <- -log(tbl$Prob)
    if(identical(tbl$parent, res$ID)){
      res$IC <- tbl$IC
    }else{
      stop("IDs are inconsistent. Check the code of compute_IC_CO_asrt().")
    }
    tidy[[categories_woALL[k]]] <- res
  }
  return(tidy)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Reformat the result of collect_CO().
#'
#' This function reformats the result of collect_DO().
#'
#' @param dict A result of collect_CO().
#' @param orgdb A genome annotation package.
#'
#' @return A formatted database.
#' @export
#'
format_CO_asrt <- function(dict, orgdb){
  #--------------------------------------------------
  # Reformat dict.
  #--------------------------------------------------
  df <- data.frame(ID = dict$ID, Description = dict$Description, IC = NA,
                   Count = NA, Gene = NA, GeneID = NA)
  df <- unique(df)
  for(i in seq_len(nrow(df))){
    genes <- unique(dict[which(dict$ID == df$ID[i]), ]$Symbol)
    df$Gene[i] <- paste(genes, collapse = "/")
    df$Count[i] <- length(genes)
  }
  rownames(df) <- seq_len(nrow(df))
  res <- list(cell = df)
  #--------------------------------------------------
  # Compute information contents, defined in
  # Mistry and Pavlidis, BMC Bioinformatics, 2008.
  #--------------------------------------------------
  treeTable <- make_treeTable_CO_asrt(tidy = res)
  res <- compute_IC_CO_asrt(dict = dict, tidy = res, treeTable = treeTable)
  #--------------------------------------------------
  # Fix the slots of gene symbols and ENTREZ Gene IDs.
  #--------------------------------------------------
  dictionary <- AnnotationDbi::select(orgdb, key = unique(dict$Symbol),
                                      columns = "ENTREZID", keytype = "SYMBOL")
  categories <- names(res)
  for(k in seq_len(length(categories))){
    for(i in seq_len(nrow(res[[categories[k]]]))){
      genes <- c()
      geneIDs <- c()
      g <- unlist(strsplit(res[[categories[k]]]$Gene[i], "/"))
      if(length(g) == 0){
        next
      }
      for(j in seq_len(length(g))){
        ind <- which(dictionary$SYMBOL == g[j])
        if(!is.na(dictionary[ind, ]$ENTREZID[1])){
          geneIDs <- c(geneIDs, dictionary[ind, ]$ENTREZID[1])
          genes <- c(genes, g[j])
        }
      }
      res[[categories[k]]]$Gene[i] <- paste(genes, collapse = "/")
      res[[categories[k]]]$GeneID[i] <- paste(geneIDs, collapse = "/")
      res[[categories[k]]]$Count[i] <- as.integer(length(geneIDs))
    }
  }
  tidy <- list()
  categories <- names(res)
  for(k in seq_len(length(categories))){
    tidy[[categories[k]]] <- res[[categories[k]]]
  }
  #--------------------------------------------------
  # Compute a similarity matrix.
  #--------------------------------------------------
  categories <- names(res)
  for(k in seq_len(length(categories))){
    df <- res[[categories[k]]]
    simmat <- matrix(0, nrow = nrow(df), ncol = nrow(df))
    tree <- treeTable[[categories[k]]]
    for(i in seq_len(nrow(df) - 1)){
      for(j in seq(i + 1, nrow(df))){
        #------------------------------
        # IC of most informative common ancestor (MICA)
        #------------------------------
        ancestors_i <- tree[which(tree$child == df$ID[i]), ]$parent
        ancestors_j <- tree[which(tree$child == df$ID[j]), ]$parent
        common_ancestors <- intersect(ancestors_i, ancestors_j)
        if(length(common_ancestors) == 0){
          simmat[i, j] <- 0
        }else{
          common_ancestors <- data.frame(ID = common_ancestors, IC = NA)
          for(n in seq_len(nrow(common_ancestors))){
            common_ancestors$IC[n] <-
              df[which(df$ID == common_ancestors$ID[n]), ]$IC
          }
          inds <- order(common_ancestors$IC, decreasing = TRUE)
          IC_MICA <- common_ancestors[inds, ]$IC[1]
          #------------------------------
          # Jiang's method
          #------------------------------
          simmat[i, j] <- 2 * IC_MICA / (df$IC[i] + df$IC[j])
        }
      }
    }
    simmat[lower.tri(simmat)] <- simmat[upper.tri(simmat)]
    diag(simmat) <- 1
    rownames(simmat) <- df$ID
    colnames(simmat) <- df$ID
    tidy[["similarity_matrix"]][[categories[k]]] <- simmat
  }

  return(tidy)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Do groupGO().
#'
#' This function performs groupGO().
#'
#' @param genes A gene set.
#' @param orgdb A genome annotation package.
#' @param ont Category of Gene Ontology.
#' @param level Level of Gene Ontology.
#'
#' @return Results of groupGO().
#'
do_groupGO <- function(genes, orgdb, ont, level){
  res <- clusterProfiler::groupGO(
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
#' Collect Gene Ontology database.
#'
#' This function collects Gene Ontology database.
#'
#' @param orgdb A genome annotation package.
#'
#' @return Results from groupGO().
#' @export
#'
collect_GO_asrt <- function(orgdb){
  #--------------------------------------------------
  # Preparation
  #--------------------------------------------------
  if(identical(orgdb, org.Hs.eg.db::org.Hs.eg.db)){
    all_geneIDs <- as.list(org.Hs.eg.db::org.Hs.egGO)
  }else if(identical(orgdb, org.Mm.eg.db::org.Mm.eg.db)){
    all_geneIDs <- as.list(org.Mm.eg.db::org.Mm.egGO)
  }else{
    stop("Currently, only org.Hs.eg.db and org.Mm.eg.db are acceptable.")
  }
  all_geneIDs <- names(all_geneIDs[!is.na(all_geneIDs)])
  categories <- c("MF", "BP", "CC")
  #--------------------------------------------------
  # groupGO
  #--------------------------------------------------
  res <- list()
  for(k in seq_len(length(categories))){
    tmp <- c()
    level <- 1
    while(1){
      ggo <- try(
        do_groupGO(genes = all_geneIDs, orgdb = orgdb, ont = categories[k],
                   level = level),
        silent = TRUE
      )
      if(class(ggo) == "try-error"){
        break
      }else{
        tmp <- c(tmp, list(ggo))
        level <- level + 1
      }
    }
    res[[categories[k]]] <- tmp
  }

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Reformat the result of collect_GO().
#'
#' This function reformats the result of collect_GO().
#'
#' @param dict A result of collect_GO().
#' @param orgdb A genome annotation package.
#'
#' @return A formatted database.
#' @export
#'
format_GO_asrt <- function(dict, orgdb){
  #--------------------------------------------------
  # Reformat dict.
  #--------------------------------------------------
  categories <- names(dict)
  res <- list()
  for(k in seq_len(length(categories))){
    tmp <- c()
    levels <- length(dict[[categories[k]]])
    for(lv in seq_len(levels)){
      tmp <- rbind(tmp, as.data.frame(dict[[categories[k]]][[lv]]@result))
    }
    tmp <- unique(tmp)  # Notice: the right hand side must be data.frame
    df <- data.frame(ID = tmp$ID, Description = tmp$Description, IC = NA,
                     Count = tmp$Count, Gene = NA, GeneID = tmp$geneID)
    df$Count <- as.integer(df$Count)
    res[[categories[k]]] <- df
  }
  #--------------------------------------------------
  # Compute information contents.
  #--------------------------------------------------
  res_godata <- list()
  for(k in seq_len(length(categories))){
    res_godata[[categories[k]]] <- GOSemSim::godata(OrgDb = orgdb,
                                                    ont = categories[k],
                                                    computeIC = TRUE)
  }
  for(k in seq_len(length(categories))){
    simdata <- res_godata[[categories[k]]]
    df <- res[[categories[k]]]
    for(i in seq_len(nrow(df))){
      if(is.element(df$ID[i], names(simdata@IC))){
        df$IC[i] <- simdata@IC[which(names(simdata@IC) == df$ID[i])]
        df$IC[i] <- ifelse(is.infinite(df$IC[i]), 99, df$IC[i])
      }else{
        df$IC[i] <- 99
      }
    }
    res[[categories[k]]] <- df
  }
  #--------------------------------------------------
  # Fix the slots of gene symbols and ENTREZ Gene IDs.
  #--------------------------------------------------
  geneIDs_MF <- res[["MF"]][which(res[["MF"]]$ID == "GO:0003674"), ]$GeneID
  geneIDs_MF <- unlist(strsplit(geneIDs_MF, "/"))
  geneIDs_BP <- res[["BP"]][which(res[["BP"]]$ID == "GO:0008150"), ]$GeneID
  geneIDs_BP <- unlist(strsplit(geneIDs_BP, "/"))
  geneIDs_CC <- res[["CC"]][which(res[["CC"]]$ID == "GO:0005575"), ]$GeneID
  geneIDs_CC <- unlist(strsplit(geneIDs_CC, "/"))
  geneIDs <- unique(c(geneIDs_MF, geneIDs_BP, geneIDs_CC))
  dictionary <- AnnotationDbi::select(orgdb, key = geneIDs, columns = "SYMBOL",
                                      keytype = "ENTREZID")
  categories <- names(res)
  for(k in seq_len(length(categories))){
    for(i in seq_len(nrow(res[[categories[k]]]))){
      genes <- c()
      geneIDs <- c()
      g <- unlist(strsplit(res[[categories[k]]]$GeneID[i], "/"))
      if(length(g) == 0){
        next
      }
      for(j in seq_len(length(g))){
        ind <- which(dictionary$ENTREZID == g[j])
        if(!is.na(dictionary[ind, ]$SYMBOL[1])){
          genes <- c(genes, dictionary[ind, ]$SYMBOL[1])
          geneIDs <- c(geneIDs, g[j])
        }
      }
      res[[categories[k]]]$Gene[i] <- paste(genes, collapse = "/")
      res[[categories[k]]]$GeneID[i] <- paste(geneIDs, collapse = "/")
      res[[categories[k]]]$Count[i] <- as.integer(length(geneIDs))
    }
  }
  tidy <- list()
  for(k in seq_len(length(categories))){
    tidy[[categories[k]]] <- res[[categories[k]]]
  }
  #--------------------------------------------------
  # Compute a similarity matrix.
  #--------------------------------------------------
  for(k in seq_len(length(categories))){
    df <- res[[categories[k]]]
    simmat <- GOSemSim::mgoSim(df$ID, df$ID,
                               semData = res_godata[[categories[k]]],
                               measure = "Jiang", combine = NULL)
    tidy[["similarity_matrix"]][[categories[k]]] <- simmat
  }
  
  return(tidy)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Do keggGet().
#'
#' This function performs keggGet().
#'
#' @param data KEGG identifiers.
#'
#' @return Results of keggGet().
#'
do_keggGet <- function(data){
  return(KEGGREST::keggGet(data))
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Collect KEGG database.
#'
#' This function collects KEGG database.
#'
#' @param organism An identifier of organism.
#' @param categories Category name.
#'
#' @return Results from keggGet().
#' @export
#'
collect_KEGG_asrt <- function(organism, categories){
  #--------------------------------------------------
  # Collect KEGG terms.
  #--------------------------------------------------
  res <- list()
  for(k in seq_len(length(categories))){
    ids <- KEGGREST::keggLink(categories[k], organism) #names(ids) are gene IDs
    map <- data.frame(
      ID = ids,
      Description = NA,
      KEGG_geneID = names(ids)
    )
    map$NCBI_geneID <- NA
    flags <- c()
    I <- length(map$ID)
    for(i in seq_len(I)){
      Sys.sleep(0.1)
      #--------------------------------------------------
      # Print the current process.
      #--------------------------------------------------
      if((i == 1) || (i %% floor(0.25 * I) == 0 && i < 0.95 * I) || (i == I)){
        text <- paste("Now processing ", i, "/", I, " for ",
                      categories[k], "...\n", sep = "")
        cat(text)
      }
      #--------------------------------------------------
      # (i) Description
      #--------------------------------------------------
      category_id <- map$ID[i]
      if(categories[k] == "module")
        category_id <- stringr::str_extract(category_id, "(?<=_)(.*)")
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
    res[[categories[k]]][["success"]] <- map
    #--------------------------------------------------
    # Rescue the failures.
    #--------------------------------------------------
    if(length(flags) > 0){
      failure <- c()
      for(i in flags){
        #------------------------------
        # (i) Description
        #------------------------------
        category_id <- map$ID[i]
        if(categories[k] == "module")
          category_id <- stringr::str_extract(category_id, "(?<=_)(.*)")
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
      res[[categories[k]]][["success"]] <- map
      tmp <- data.frame(
        matrix(ncol = 2, nrow = 0, dimnames = list(NULL, c(
          "Input_category_for_keggGet", "Input_for_keggGet"
        ))))
      tmp <- rbind(tmp, data.frame(
        Input_category_for_keggGet = unique(failure)[, 1],
        Input_for_keggGet = unique(failure)[, 2]
      ))
      res[[categories[k]]][["failure"]] <- tmp
    }
  }

  return(res)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
#' Reformat collect_KEGG().
#'
#' This function reformats collect_KEGG().
#'
#' @param dict A result of collect_KEGG().
#' @param orgdb A genome annotation package.
#'
#' @return A formatted database.
#' @export
#'
format_KEGG_asrt <- function(dict, orgdb){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  categories <- names(dict)

  #--------------------------------------------------
  # Reformat
  #--------------------------------------------------
  res <- list()
  for(k in seq_len(length(categories))){
    tmp <- dict[[categories[k]]]
    map <- unique(data.frame(
      ID = tmp[["ID"]],
      Description = tmp[["Description"]],
      IC = NA,
      Count = NA,
      Gene = NA,
      GeneID = NA
    ))
    for(i in seq_len(nrow(map))){
      #------------------------------
      # Gene and Count
      #------------------------------
      genes <- unique(tmp[which(tmp[["ID"]] == map[["ID"]][i]), ]$NCBI_geneID)
      map$GeneID[i] <- paste(genes, collapse = "/")
      map$Count[i] <- length(genes)
    }
    rownames(map) <- seq_len(nrow(map))
    res[[categories[k]]] <- map
  }
  #--------------------------------------------------
  # Fix the slots of gene symbols and ENTREZ Gene IDs.
  #--------------------------------------------------
  for(k in seq_len(length(categories))){
    geneIDs <- unique(dict[[categories[k]]][["NCBI_geneID"]])
    dictionary <- AnnotationDbi::select(orgdb, key = geneIDs,
                                        columns = "SYMBOL",
                                        keytype = "ENTREZID")
    for(i in seq_len(nrow(res[[categories[k]]]))){
      genes <- c()
      geneIDs <- c()
      g <- unlist(strsplit(res[[categories[k]]]$GeneID[i], "/"))
      if(length(g) == 0){
        next
      }
      for(j in seq_len(length(g))){
        ind <- which(dictionary$ENTREZID == g[j])
        if(!is.na(dictionary[ind, ]$SYMBOL[1])){
          genes <- c(genes, dictionary[ind, ]$SYMBOL[1])
          geneIDs <- c(geneIDs, g[j])
        }
      }
      res[[categories[k]]]$Gene[i] <- paste(genes, collapse = "/")
      res[[categories[k]]]$GeneID[i] <- paste(geneIDs, collapse = "/")
      res[[categories[k]]]$Count[i] <- as.integer(length(geneIDs))
    }
  }
  tidy <- list()
  categories <- names(res)
  for(k in seq_len(length(categories))){
    tidy[[categories[k]]] <- res[[categories[k]]]
  }

  return(tidy)
}

