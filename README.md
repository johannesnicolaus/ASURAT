# ASURAT
## Brief introduction
ASURAT is a computational pipeline, implemented in the R programming language, for single-cell RNA sequencing (scRNA-seq) data analysis.
Using ASURAT, one can simultaneously perform unsupervised clustering and biological interpretation in terms of cell type, disease, biological process, and signaling pathway activity.

## Graphical abstract
<img src="figures/figure_00_0001.png" width="600px">

## Vignettes
Well-documented vignette and tutorial are available from the following URL:

* https://keita-iida.github.io/ASURAT/,

in which we analyze public scRNA-seq data of human small cell lung cancer (SCLC) and peripheral blood mononuclear cells (PBMCs).

## Installation
One can install ASURAT by the following code:

```{r}
devtools::install_github("keita-iida/ASURAT", upgrade = "never")
```

## Requirements
ASURAT requires to install the following R packages:

* SingleCellExperiment
* SummerizedExperiment
* cluster

## Preprint
https://www.biorxiv.org/content/10.1101/2021.06.09.447731v3.full

## Author
Keita Iida

## License
[GPL-3](https://github.com/keita-iida/ASURAT/blob/main/LICENSE)
