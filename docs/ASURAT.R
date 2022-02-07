## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = FALSE------------------------------------------------------------
#  devtools::install_github("keita-iida/ASURAT", upgrade = "never")

## ---- message = FALSE, warning = FALSE----------------------------------------
library(ASURAT)
library(SingleCellExperiment)
library(SummarizedExperiment)

## ---- eval = FALSE------------------------------------------------------------
#  mat <- as.matrix(assay(sce, "logcounts"))
#  assay(sce, "centered") <- sweep(mat, 1, apply(mat, 1, mean), FUN = "-")

## -----------------------------------------------------------------------------
urlpath <- "https://github.com/keita-iida/ASURAT_DB/blob/main/transcriptome/"
load(url(paste0(urlpath, "pbmc_counts.rda?raw=true")))

## -----------------------------------------------------------------------------
pbmc_counts[1:5, 1:3]

## -----------------------------------------------------------------------------
pbmc <- SingleCellExperiment(assays = list(counts = pbmc_counts),
                             rowData = data.frame(gene = rownames(pbmc_counts)),
                             colData = data.frame(cell = colnames(pbmc_counts)))

## -----------------------------------------------------------------------------
dim(pbmc)

## -----------------------------------------------------------------------------
pbmc <- add_metadata(sce = pbmc, mitochondria_symbol = "^MT-")

## -----------------------------------------------------------------------------
pbmc <- remove_variables(sce = pbmc, min_nsamples = 10)

## -----------------------------------------------------------------------------
dataframe2D <- data.frame(x = colData(pbmc)$nReads, y = colData(pbmc)$nGenes)
plot_dataframe2D(dataframe2D = dataframe2D) +
  ggplot2::labs(x = "Number of reads", y = "Number of genes") +
  ggplot2::theme_classic(base_size = 20)

## -----------------------------------------------------------------------------
dataframe2D <- data.frame(x = colData(pbmc)$nReads, y = colData(pbmc)$percMT)
plot_dataframe2D(dataframe2D = dataframe2D) +
  ggplot2::labs(x = "Number of reads", y = "Perc of MT reads") +
  ggplot2::theme_classic(base_size = 20)

## -----------------------------------------------------------------------------
pbmc <- remove_samples(sce = pbmc, min_nReads = 0, max_nReads = 1e+10,
                       min_nGenes = 0, max_nGenes = 1e+10,
                       min_percMT = NULL, max_percMT = NULL)

## -----------------------------------------------------------------------------
dataframe2D <- data.frame(x = 1:nrow(rowData(pbmc)),
                          y = sort(rowData(pbmc)$nSamples, decreasing = TRUE))
plot_dataframe2D(dataframe2D = dataframe2D) +
  ggplot2::labs(x = "Rank of genes", y = "Mean read counts") +
  ggplot2::theme_classic(base_size = 20)

## -----------------------------------------------------------------------------
pbmc <- remove_variables_second(sce = pbmc, min_meannReads = 0.01)

## ---- results = "hide", message = FALSE---------------------------------------
bayout <- bayNorm::bayNorm(Data = counts(pbmc), mode_version = TRUE)
assay(pbmc, "normalized") <- bayout$Bay_out

## -----------------------------------------------------------------------------
assay(pbmc, "logcounts") <- log(assay(pbmc, "normalized") + 1)

## -----------------------------------------------------------------------------
mat <- assay(pbmc, "logcounts")
assay(pbmc, "centered") <- sweep(mat, 1, apply(mat, 1, mean), FUN = "-")

## -----------------------------------------------------------------------------
urlpath <- "https://github.com/keita-iida/ASURAT_DB/blob/main/genes2bioterm/"
load(url(paste0(urlpath, "20220108_human_COMSig.rda?raw=TRUE"))) # CO & MSigDB
load(url(paste0(urlpath, "20201213_human_GO_red.rda?raw=TRUE"))) # GO
load(url(paste0(urlpath, "20201213_human_KEGG.rda?raw=TRUE")))   # KEGG

## -----------------------------------------------------------------------------
pbmc_cormat <- cor(t(assay(pbmc, "centered")), method = "spearman")

## -----------------------------------------------------------------------------
sname <- "logcounts"
altExp(pbmc, sname) <- SummarizedExperiment(list(counts = assay(pbmc, sname)))

## ---- message = FALSE, warning = FALSE----------------------------------------
dictionary <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                                    key = rownames(pbmc),
                                    columns = "ENTREZID", keytype = "SYMBOL")
dictionary <- dictionary[!duplicated(dictionary$SYMBOL), ]
rowData(pbmc)$geneID <- dictionary$ENTREZID

## -----------------------------------------------------------------------------
pbmcs <- list(CM = pbmc, GO = pbmc, KG = pbmc)
metadata(pbmcs$CM) <- list(sign = human_COMSig[["cell"]])
metadata(pbmcs$GO) <- list(sign = human_GO[["BP"]])
metadata(pbmcs$KG) <- list(sign = human_KEGG[["pathway"]])

## -----------------------------------------------------------------------------
pbmcs$CM <- remove_signs(sce = pbmcs$CM, min_ngenes = 2, max_ngenes = 1000)
pbmcs$GO <- remove_signs(sce = pbmcs$GO, min_ngenes = 2, max_ngenes = 1000)
pbmcs$KG <- remove_signs(sce = pbmcs$KG, min_ngenes = 2, max_ngenes = 1000)

## -----------------------------------------------------------------------------
set.seed(8)
pbmcs$CM <- cluster_genesets(sce = pbmcs$CM, cormat = pbmc_cormat,
                             th_posi = 0.20, th_nega = -0.20)
pbmcs$GO <- cluster_genesets(sce = pbmcs$GO, cormat = pbmc_cormat,
                             th_posi = 0.24, th_nega = -0.20)
pbmcs$KG <- cluster_genesets(sce = pbmcs$KG, cormat = pbmc_cormat,
                             th_posi = 0.20, th_nega = -0.20)

## -----------------------------------------------------------------------------
pbmcs$CM <- create_signs(sce = pbmcs$CM, min_cnt_strg = 2, min_cnt_vari = 2)
pbmcs$GO <- create_signs(sce = pbmcs$GO, min_cnt_strg = 2, min_cnt_vari = 2)
pbmcs$KG <- create_signs(sce = pbmcs$KG, min_cnt_strg = 2, min_cnt_vari = 2)

## -----------------------------------------------------------------------------
pbmcs$GO <- remove_signs_redundant(
  sce = pbmcs$GO, similarity_matrix = human_GO$similarity_matrix$BP,
  threshold = 0.80, keep_rareID = TRUE)

## -----------------------------------------------------------------------------
keywords <- "Covid19|foofoo|hogehoge"
pbmcs$KG <- remove_signs_manually(sce = pbmcs$KG, keywords = keywords)

## -----------------------------------------------------------------------------
keywords <- "lung|pleural|respiratory"
test <- select_signs_manually(sce = pbmcs$CM, keywords = keywords)

## -----------------------------------------------------------------------------
pbmcs$CM <- makeSignMatrix(sce = pbmcs$CM, weight_strg = 0.5, weight_vari = 0.5)
pbmcs$GO <- makeSignMatrix(sce = pbmcs$GO, weight_strg = 0.5, weight_vari = 0.5)
pbmcs$KG <- makeSignMatrix(sce = pbmcs$KG, weight_strg = 0.5, weight_vari = 0.5)

## -----------------------------------------------------------------------------
rbind(head(assay(pbmcs$CM, "counts")[, 1:3], n = 4),
      tail(assay(pbmcs$CM, "counts")[, 1:3], n = 4))

## -----------------------------------------------------------------------------
set.seed(1)
for(i in seq_len(length(pbmcs))){
  res <- Rtsne::Rtsne(t(assay(pbmcs[[i]], "counts")), dim = 2, pca = FALSE)
  reducedDim(pbmcs[[i]], "TSNE") <- res[["Y"]]
}

## -----------------------------------------------------------------------------
set.seed(1)
res <- umap::umap(t(assay(pbmcs$CM, "counts")), n_components = 3)
reducedDim(pbmcs$CM, "UMAP") <- res[["layout"]]

## -----------------------------------------------------------------------------
dataframe2D <- as.data.frame(reducedDim(pbmcs$CM, "TSNE"))
plot_dataframe2D(dataframe2D = dataframe2D) +
  ggplot2::labs(title = "PBMC (CO & MSigDB)", x = "tSNE_1", y = "tSNE_2") +
  ggplot2::theme_classic(base_size = 15)

## -----------------------------------------------------------------------------
dataframe2D <- as.data.frame(reducedDim(pbmcs$GO, "TSNE"))
plot_dataframe2D(dataframe2D = dataframe2D) +
  ggplot2::labs(title = "PBMC (GO)", x = "tSNE_1", y = "tSNE_2") +
  ggplot2::theme_classic(base_size = 15)

## -----------------------------------------------------------------------------
dataframe2D <- as.data.frame(reducedDim(pbmcs$KG, "TSNE"))
plot_dataframe2D(dataframe2D = dataframe2D) +
  ggplot2::labs(title = "PBMC (KEGG)", x = "tSNE_1", y = "tSNE_2") +
  ggplot2::theme_classic(base_size = 15)

## ---- message = FALSE, warning = FALSE, results = "hide"----------------------
resolutions <- c(0.15, 0.15, 0.25)
for(i in seq_len(length(pbmcs))){
  pbmc_surt <- Seurat::as.Seurat(pbmcs[[i]], counts = "counts", data = "counts")
  pbmc_surt[["SSM"]] <- Seurat::CreateAssayObject(counts = as.matrix(assay(pbmcs[[i]], "counts")))
  Seurat::DefaultAssay(pbmc_surt) <- "SSM"
  pbmc_surt <- Seurat::ScaleData(pbmc_surt, features = rownames(pbmc_surt))
  pbmc_surt <- Seurat::RunPCA(pbmc_surt, features = rownames(pbmc_surt))
  pbmc_surt <- Seurat::FindNeighbors(pbmc_surt, reduction = "pca", dims = 1:20)
  pbmc_surt <- Seurat::FindClusters(pbmc_surt, resolution = resolutions[i])
  pbmc_temp <- Seurat::as.SingleCellExperiment(pbmc_surt)
  colData(pbmcs[[i]])$seurat_clusters <- colData(pbmc_temp)$seurat_clusters
}

## -----------------------------------------------------------------------------
labels <- colData(pbmcs$CM)$seurat_clusters
dataframe2D <- as.data.frame(reducedDim(pbmcs$CM, "TSNE"))
plot_dataframe2D(dataframe2D = dataframe2D, labels = labels, colors = NULL) +
  ggplot2::labs(title = "PBMC (CO & MSigDB)", x = "tSNE_1", y = "tSNE_2", color = "") +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4)))

## -----------------------------------------------------------------------------
labels <- colData(pbmcs$GO)$seurat_clusters
dataframe2D <- as.data.frame(reducedDim(pbmcs$GO, "TSNE"))
plot_dataframe2D(dataframe2D = dataframe2D, labels = labels, colors = NULL) +
  ggplot2::labs(title = "PBMC (GO)", x = "tSNE_1", y = "tSNE_2", color = "") +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4)))

## -----------------------------------------------------------------------------
labels <- colData(pbmcs$KG)$seurat_clusters
dataframe2D <- as.data.frame(reducedDim(pbmcs$KG, "TSNE"))
plot_dataframe2D(dataframe2D = dataframe2D, labels = labels, colors = NULL) +
  ggplot2::labs(title = "PBMC (KEGG)", x = "tSNE_1", y = "tSNE_2", color = "") +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4)))

## ---- message = FALSE, warning = FALSE----------------------------------------
pbmc_surt <- Seurat::as.Seurat(pbmcs$CM, counts = "counts", data = "counts")
pbmc_surt[["GEM"]] <- Seurat::CreateAssayObject(counts = as.matrix(assay(altExp(pbmcs$CM), "counts")))
Seurat::DefaultAssay(pbmc_surt) <- "GEM"
pbmc_surt <- Seurat::ScaleData(pbmc_surt, features = rownames(pbmc_surt))
pbmc_surt <- Seurat::RunPCA(pbmc_surt, features = rownames(pbmc_surt))
pbmc_surt <- Seurat::CellCycleScoring(pbmc_surt,
                                      s.features = Seurat::cc.genes$s.genes,
                                      g2m.features = Seurat::cc.genes$g2m.genes)
pbmc_temp <- Seurat::as.SingleCellExperiment(pbmc_surt)
colData(pbmcs$CM)$Phase <- colData(pbmc_temp)$Phase

## ---- message = FALSE, results = "hide"---------------------------------------
for(i in seq_len(length(pbmcs))){
  set.seed(1)
  labels <- colData(pbmcs[[i]])$seurat_clusters
  pbmcs[[i]] <- compute_sepI_all(sce = pbmcs[[i]], labels = labels,
                                 nrand_samples = 1000)
}

## -----------------------------------------------------------------------------
vname <- "MSigID:92-S"
pbmc_sub <- pbmcs$CM[rownames(pbmcs$CM) %in% vname, ]
labels <- colData(pbmc_sub)$seurat_clusters
dataframe1D <- as.data.frame(t(assay(pbmc_sub, "counts")))
plot_violin(dataframe1D = dataframe1D, labels = labels, colors = NULL) +
  ggplot2::labs(title = paste0(vname, "\n", "B cell (CD79A, CD22, ...)"),
                x = "Cluster (CO & MSigDB)", y = "Sign score", fill = "Cluster") +
  ggplot2::theme_classic(base_size = 20)

## ---- message = FALSE, warning = FALSE----------------------------------------
pbmc_surt <- Seurat::as.Seurat(pbmcs$CM, counts = "counts", data = "counts")
#If there is "expressions" data in `altExp(sce)`, set Seurat default assay.
Seurat::DefaultAssay(pbmc_surt) <- "logcounts"
pbmc_surt <- Seurat::SetIdent(pbmc_surt, value = "seurat_clusters")
res <- Seurat::FindAllMarkers(pbmc_surt, only.pos = TRUE,
                              min.pct = 0.25, logfc.threshold = 0.25)
metadata(pbmcs$CM)$marker_genes$all <- res

## -----------------------------------------------------------------------------
vname <- "BANK1"
expr_sub <- altExp(pbmcs$CM, "logcounts")
expr_sub <- expr_sub[rownames(expr_sub) == vname]
labels <- colData(pbmcs$CM)$seurat_clusters
dataframe1D <- as.data.frame(t(assay(expr_sub, "counts")))
plot_violin(dataframe1D = dataframe1D, labels = labels, colors = NULL) +
  ggplot2::labs(title = "BANK1", x = "Cluster (CO & MSigDB)",
                y = "Normalized expression", fill = "Cluster") +
  ggplot2::theme_classic(base_size = 20)

## -----------------------------------------------------------------------------
# Significant signs
marker_signs <- list()
for(i in seq_len(length(pbmcs))){
  marker_signs[[i]] <- metadata(pbmcs[[i]])$marker_signs$all
  marker_signs[[i]] <- dplyr::group_by(marker_signs[[i]], Ident_1)
  marker_signs[[i]] <- dplyr::slice_max(marker_signs[[i]], sepI, n = 1)
}
# Significant genes
marker_genes_CM <- metadata(pbmcs$CM)$marker_genes$all
marker_genes_CM <- dplyr::group_by(marker_genes_CM, cluster)
marker_genes_CM <- dplyr::slice_min(marker_genes_CM, p_val_adj, n = 1)

## -----------------------------------------------------------------------------
# ssm_list
pbmcs_sub <- list() ; ssm_list <- list()
for(i in seq_len(length(pbmcs))){
  pbmcs_sub[[i]] <- pbmcs[[i]][rownames(pbmcs[[i]]) %in% marker_signs[[i]]$SignID, ]
  ssm_list[[i]] <- assay(pbmcs_sub[[i]], "counts")
}
names(ssm_list) <- c("SSM_COMSig", "SSM_GO", "SSM_KEGG")
# gem_list
expr_sub <- altExp(pbmcs$CM, "logcounts")
expr_sub <- expr_sub[rownames(expr_sub) %in% marker_genes_CM$gene]
gem_list <- list(GeneExpr = assay(expr_sub, "counts"))
# ssmlabel_list
labels <- list() ; ssmlabel_list <- list()
for(i in seq_len(length(pbmcs))){
  labels[[i]] <- data.frame(label = colData(pbmcs_sub[[i]])$seurat_clusters)
  labels[[i]]$color <- NA
  ssmlabel_list[[i]] <- labels[[i]]
}
names(ssmlabel_list) <- c("Label_COMSig", "Label_GO", "Label_KEGG")
# gemlabel_list
label_CC <- data.frame(label = colData(pbmcs$CM)$Phase, color = NA)
gemlabel_list <- list(CellCycle = label_CC)

## ---- message = FALSE, warning = FALSE----------------------------------------
set.seed(1)
plot_multiheatmaps(ssm_list = ssm_list, gem_list = gem_list,
                   ssmlabel_list = ssmlabel_list, gemlabel_list = gemlabel_list,
                   nSamples = 100, show_row_names = TRUE, title = "PBMC")

## -----------------------------------------------------------------------------
vname <- "MSigID:162-S"
pbmc_sub <- pbmcs$CM[rownames(pbmcs$CM) %in% vname, ]
labels <- colData(pbmc_sub)$seurat_clusters
colors <- rainbow(length(unique(labels)))[labels]
dataframe1D <- as.data.frame(t(assay(pbmc_sub, "counts")))
plot_violin(dataframe1D = dataframe1D, labels = labels, colors = colors) +
  ggplot2::labs(title = paste0(vname, "\n", "Naive T cell (CD3E, TRAC, ...)"),
                x = "Cluster (CO & MSigDB)", y = "Sign score", fill = "Cluster") +
  ggplot2::theme_classic(base_size = 20)

## -----------------------------------------------------------------------------
labels <- colData(pbmcs$CM)$seurat_clusters
colors <- rainbow(length(unique(labels)))[labels]
dataframe2D <- as.data.frame(reducedDim(pbmcs$CM, "TSNE"))
plot_dataframe2D(dataframe2D = dataframe2D, labels = labels, colors = colors) +
  ggplot2::labs(title = "PBMC (CO & MSigDB)", x = "tSNE_1", y = "tSNE_2", color = "") +
  ggplot2::theme_classic(base_size = 15) +
  ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 4)))

## -----------------------------------------------------------------------------
dataframe3D <- as.data.frame(reducedDim(pbmcs$CM, "UMAP")[, 1:3])
labels <- colData(pbmcs$CM)$seurat_clusters
plot_dataframe3D(dataframe3D = dataframe3D, labels = labels, colors = NULL,
                 theta = 45, phi = 20, title = "PBMC (CO & MSigDB)",
                 xlabel = "UMAP_1", ylabel = "UMAP_2", zlabel = "UMAP_3")

## -----------------------------------------------------------------------------
# ssm_list
ssm_list <- list(SSM_COMSig = assay(pbmcs$CM, "counts"),
                 SSM_GO     = assay(pbmcs$GO, "counts"),
                 SSM_KEGG   = assay(pbmcs$KG, "counts"))
# gem_list
gem_list <- list(GeneExpr = assay(altExp(pbmcs$CM, "logcounts"), "counts"))
# ssmlabel_list
labels <- list() ; ssmlabel_list <- list()
for(i in seq_len(length(pbmcs))){
  labels[[i]] <- data.frame(label = colData(pbmcs[[i]])$seurat_clusters)
  labels[[i]]$color <- rainbow(length(unique(labels[[i]]$label)))[labels[[i]]$label]
  ssmlabel_list[[i]] <- labels[[i]]
}
names(ssmlabel_list) <- c("Label_COMSig", "Label_GO", "Label_KEGG")
# gemlabel_list
label_CC <- data.frame(label = colData(pbmcs$CM)$Phase, color = NA)
gemlabel_list <- list(CellCycle = label_CC)

## -----------------------------------------------------------------------------
set.seed(1)
plot_multiheatmaps(ssm_list = ssm_list, gem_list = gem_list,
                   ssmlabel_list = ssmlabel_list, gemlabel_list = gemlabel_list,
                   nSamples = 100, show_row_names = FALSE, title = "PBMC")

## -----------------------------------------------------------------------------
sessionInfo()

