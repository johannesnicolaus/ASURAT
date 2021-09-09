#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
make_asurat_obj <- function(mat, obj_name){
  #--------------------------------------------------
  # Definition
  #--------------------------------------------------
  obj <- list(
    history = c(),
    variable = data.frame(
      symbol = rownames(mat),
      entrez = NA
    ),
    sample = data.frame(
      barcode = colnames(mat),
      orig_ident = obj_name
    ),
    data = c(),
    misc = c()
  )
  obj[["data"]][["raw"]] <- as.data.frame(mat)
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- "make_asurat_obj"
  obj[["history"]][[slot_name]][["obj_name"]] <- obj_name

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
do_quickQC <- function(obj, min_nsamples, mitochondria_symbol){
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- "do_quickQC"
  obj[["history"]][[slot_name]][["min_nsamples"]] <- min_nsamples
  #--------------------------------------------------
  # Reduce variables
  #--------------------------------------------------
  tmp <- as.matrix(obj[["data"]][["raw"]])
  inds <- which(apply(tmp, 1, function(x) sum(x>0)) >= min_nsamples)
  mat <- tmp[inds,]
  obj[["data"]][["raw"]] <- as.data.frame(mat)
  #--------------------------------------------------
  # Other information
  #--------------------------------------------------
  obj[["variable"]] <- obj[["variable"]][inds,]
  obj[["sample"]][["nReads"]] <- as.integer(apply(mat, 2, function(x) sum(x)))
  obj[["sample"]][["nGenes"]] <- as.integer(apply(mat, 2, function(x) sum(x > 0)))
  obj[["sample"]][["percent_MT"]] <- apply(
    mat[grepl(mitochondria_symbol, rownames(mat)),], 2, function(x) 100 * sum(x)
  ) / apply(mat, 2, function(x) sum(x))
  tmp <- obj[["sample"]][which(is.nan(obj[["sample"]]$percent_MT)),]
  if(nrow(tmp) > 0){
    tmp$percent_MT <- 0.0
  }
  obj[["sample"]][which(is.nan(obj[["sample"]]$percent_MT)),] <- tmp

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
trim_samples <- function(obj,
  min_nReads = 0, max_nReads = 1e+10,
  min_nGenes = 0, max_nGenes = 1e+10,
  min_percent_MT = 0, max_percent_MT = 100
){
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- "trim_samples"
  obj[["history"]][[slot_name]][["min_nReads"]] <- min_nReads
  obj[["history"]][[slot_name]][["max_nReads"]] <- max_nReads
  obj[["history"]][[slot_name]][["min_nGenes"]] <- min_nGenes
  obj[["history"]][[slot_name]][["max_nGenes"]] <- max_nGenes
  obj[["history"]][[slot_name]][["min_percent_MT"]] <- min_percent_MT
  obj[["history"]][[slot_name]][["max_percent_MT"]] <- max_percent_MT
  #--------------------------------------------------
  # Reduce count matrices
  #--------------------------------------------------
  inds_1 <- which(
    (obj[["sample"]][["nReads"]] >= min_nReads) &
    (obj[["sample"]][["nReads"]] <= max_nReads)
  )
  inds_2 <- which(
    (obj[["sample"]][["nGenes"]] >= min_nGenes) &
    (obj[["sample"]][["nGenes"]] <= max_nGenes)
  )
  inds_3 <- which(
    (obj[["sample"]][["percent_MT"]] >= min_percent_MT) &
    (obj[["sample"]][["percent_MT"]] <= max_percent_MT)
  )
  inds <- intersect(intersect(inds_1, inds_2), inds_3)
  obj[["data"]][["raw"]] <- as.data.frame(obj[["data"]][["raw"]][,inds])
  #--------------------------------------------------
  # The others
  #--------------------------------------------------
  obj[["sample"]] <- obj[["sample"]][inds,]

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
trim_variables <- function(obj, min_meanReads){
  #--------------------------------------------------
  # History
  #--------------------------------------------------
  slot_name <- "trim_variables"
  obj[["history"]][[slot_name]][["min_meanReads"]] <- min_meanReads
  #--------------------------------------------------
  # Reduce count matrix
  #--------------------------------------------------
  tmp <- as.matrix(obj[["data"]][["raw"]])
  inds <- which(apply(tmp, 1, mean) >= min_meanReads)
  mat <- tmp[inds, ]
  obj[["data"]][["raw"]] <- as.data.frame(mat)
  #--------------------------------------------------
  # The others
  #--------------------------------------------------
  obj[["variable"]] <- obj[["variable"]][inds,]

  return(obj)
}
#-----------------------------------------------------------------------------80
#
#-----------------------------------------------------------------------------80
convert_seurat2asurat <- function(obj_seurat, orgdb, obj_name){
  # ----------------------------------------
  # Set up ASURAT objects using normalized data
  # ----------------------------------------
  mat <- as.matrix(obj_seurat@assays[["integrated"]]@data)
  obj <- make_asurat_obj(mat = mat, obj_name = obj_name)
  obj[["sample"]][["orig_ident"]] <- obj_seurat@meta.data[["orig.ident"]]
  # ----------------------------------------
  # Prepare normalized and centered data slot
  # ----------------------------------------
  obj[["data"]][["raw"]] <- NULL # Remove this because the data is normalized
  obj[["data"]][["normalized"]] <- as.data.frame(mat)
  cmat <- sweep(mat, 1, apply(mat, 1, mean), FUN = "-")
  obj[["data"]][["centered"]] <- as.data.frame(cmat)
  # ----------------------------------------
  # Convert gene symbols into entrez IDs
  # ----------------------------------------
  dictionary <- AnnotationDbi::select(orgdb,
                                      key = obj[["variable"]][["symbol"]],
                                      columns = "ENTREZID", keytype = "SYMBOL")
  dictionary <- dictionary[!duplicated(dictionary$SYMBOL), ]
  names(dictionary) <- c("symbol", "entrez")
  obj[["variable"]] <- dictionary

  return(obj)
}
