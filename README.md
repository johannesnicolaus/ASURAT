# ASURAT
## Brief introduction
ASURAT is a computational pipeline, implemented in the R programming language, for single-cell RNA sequencing (scRNA-seq) data analysis.
Using ASURAT, one can simultaneously perform unsupervised clustering and biological interpretation in terms of cell type, disease, biological process, and signaling pathway activity.

## Graphical abstract
<img src="figures/figure_00_0001.png" width="500px">

## Vignettes
Well-documented vignette and tutorial is available from the following URL:

* https://keita-iida.github.io/ASURAT/,

in which we analyze public scRNA-seq data of human small cell lung cancer (SCLC) and peripheral blood mononuclear cells (PBMCs).

## Installation
One can install ASURAT by the following code:

```{r}
devtools::install_github("keita-iida/ASURAT", upgrade = "never")
```

## Author
Keita Iida

## License
[GPL-3](https://github.com/keita-iida/ASURAT/blob/main/LICENSE)
