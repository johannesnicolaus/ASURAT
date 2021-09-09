# ASURAT

ASURAT is a single-cell RNA sequencing (scRNA-seq) data analysis pipeline, developed for simultaneously clustering cells and biological interpretation.

Introduction, documentation, and tutorial can be found at

https://keita-iida.github.io/ASURAT/

Our preprint was released on bioRxiv (12th June, 2021).

https://www.biorxiv.org/content/10.1101/2021.06.09.447731v1

<br>

## ASURAT's workflow

1. Preprocessing (data quality control, normalization, etc.)
2. Collecting databases (DBs), but the following DBs were downloaded in December 2020: Disease Ontology (DO), Cell Ontology (CO), and Gene Ontology (GO) databases, Kyoto Encyclopedia of Genes and Genomes (KEGG), and Reactome.
3. Creating signs (sign is a general biological term representing cell type and various biological functions)
4. Creating sign-by-sample matrices (SSMs)
5. Unsupervised clustering of cells
6. Finding significant signs
7. Multiple sing analysis

<br>

## Quick start from a Seurat object

Although the above [URL](https://keita-iida.github.io/ASURAT/) does not assume Seurat-based analyses, it might be beneficial to begin with a Seurat object `obj` including `obj@assays[["RNA"]]@counts`.

Load a Seurat object for human scRNA-seq data (below is an example).

```R
test <- readRDS(file = "backup/test.rds")  # Seurat object
```

Below are stopgap installations (see [Chapter 1](https://keita-iida.github.io/ASURAT/) for all the requirements).
Note that users need to replace `org.Hs.eg.db` when analyzing other animal's data.

```R
library(tidyverse)                # For efficient handling of data.frame
library(org.Hs.eg.db)             # For using human genome annotation package
library(Seurat)                   # For using Seurat
source("R/function_general.R")		# ASURAT's function
```

Create an ASURAT object.

```R
test <- make_asurat_obj(mat = test@assays[["RNA"]]@counts, obj_name = "test")
```

Convert gene symbols into Entrez IDs by using `org.Hs.eg.db` package.

```R
dictionary <- AnnotationDbi::select(org.Hs.eg.db,
                                    key = test[["variable"]][["symbol"]],
                                    columns = c("ENTREZID"), keytype = "SYMBOL")
dictionary <- dictionary[!duplicated(dictionary$SYMBOL), ]
names(dictionary) <- c("symbol", "entrez")
test[["variable"]] <- dictionary
```

The following function `log1p_data()` performs log transform of the input data with a pseudo count `eps`.

```R
log1p_data <- function(obj, eps){
  obj[["history"]][["log1p_data"]][["eps"]] <- eps
  mat <- as.matrix(obj[["data"]][["raw"]])
  lmat <- log(mat + eps)
  obj[["data"]][["normalized"]] <- as.data.frame(lmat)
  return(obj)
}

test <- log1p_data(obj = test, eps = 1)
```

The following function `center_data()` centralizes the input data on a gene-by-gene basis.

```R
center_data <- function(obj){
  mat <- as.matrix(obj[["data"]][["log1p"]])
  cmat <- sweep(mat, 1, apply(mat, 1, mean), FUN = "-")
  obj[["data"]][["centered"]] <- as.data.frame(cmat)
  return(obj)
}

test <- center_data(obj = test)
```

The following function `do_cor_variables()` computes a correlation matrix from the input data.
Users can choose a measure of correlation coefficient by setting `method` (vector form is also accepted but not recommended due to the file size) such as `pearson`, `spearman`, and `kendall`.

```R
do_cor_variables <- function(obj, method){
  res <- list()
  tmat <- t(obj[["data"]][["normalized"]])
  for(m in method){
    res <- c(res, list(cor(tmat, method = m)))
  }
  names(res) <- method
  return(res)
}

test_cor <- do_cor_variables(obj = test, method = c("spearman"))
```

Save the objects. Please note that the suffixes of the following filenames, such as `09`, `005`, and `006`, are only for identifying the computational steps (there is no special significance).

```R
saveRDS(test, file = "backup/00_005_test_normalized.rds")
saveRDS(test_cor, file = "backup/00_006_test_correlation.rds")
```

### ASURAT using several databases

Go to [Chapter 9](https://keita-iida.github.io/ASURAT/asurat-using-disease-ontology-database.html) for inferring cell types or diseases based on Disease Ontology (DO) or Cell Ontology (CO).

Go to [Chapter 10](https://keita-iida.github.io/ASURAT/asurat-using-cell-ontology-database-optional.html) for inferring biological functions based on Gene Ontology (GO).

Go to [Chapter 11](https://keita-iida.github.io/ASURAT/asurat-using-gene-ontology-database-optional.html) for inferring pathway activities based on Encyclopedia of Genes and Genomes (KEGG).

Go to [Chapter 12](https://keita-iida.github.io/ASURAT/asurat-using-kegg-optional.html) for multi-layered analyses across all the signs.
