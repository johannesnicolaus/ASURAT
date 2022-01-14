# ASURAT

ASURAT is a computational pipeline, implemented in the R programming language, for single-cell RNA sequencing (scRNA-seq) data analysis.

<img src="figures/figure_00_0000.png" width="500px">

Our preprint was released on bioRxiv (12th October 2021).

https://www.biorxiv.org/content/10.1101/2021.06.09.447731v3



<br>

## The latest version of ASURAT (version 0.1.0)

ASURAT (version 0.1.0) is currently in the beta phase of its development.
This means that we are almost ready to submit it to Bioconductor.

This version is maintained in the following URL:

https://github.com/keita-iida/ASURAT_0.1.0



<br>

## What can ASURAT do?

Using ASURAT, users can transform single-cell RNA-seq data into novel sign-by-sample matrices (SSMs), in which rows and columns stand for biological terms (i.e., signs) and samples (cells), respectively. Here, biological terms are defined by functional gene sets, in terms of cell type, disease, functions, and signaling pathway activity. Thus, users can cluster single-cell transcriptomes from multiple biological aspects and obtain the most interpretable results.

https://user-images.githubusercontent.com/50622599/132992121-95c344d7-4262-45bf-83b8-25febd41365b.mov

Concatenating multiple SSMs, gene expression matrices, and other information such as cell cycle, users can characterize individual cells in columns from multiple biological terms in rows, facilitating the interpretation.
