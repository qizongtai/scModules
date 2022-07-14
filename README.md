
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scModules

<img align="right" width="108" height="125" src="man/figures/scModules1.png">

<!-- badges: start -->
<!-- badges: end -->

**scModules** is an R package for advanced processing and visualization
of single-cell RNA-seq and TCR-seq data. It integrates multiple genomic
information (e.g., copy number aberration (CNA)) and statistical
algorithms (e.g., non-negative matrix factorization (NMF)) to provide
semi-automated end-to-end analyses that are NOT available in the popular
single-cell package **Seurat**. Currently, four modules are developed.
The methodologies have been tested in several high-impact publications
including Cell and Science. scModules enhances our ability to extract
accurate genomic and cellular information from single-cell data and
discover new cancer biology.

-   [An Integrative Model of Cellular States, Plasticity, and Genetics
    for Glioblastoma (Neftel, Laffy et al., 2019,
    *Cell*)](https://doi.org/10.1016/j.cell.2019.06.024)
-   [Developmental and oncogenic programs in H3K27M gliomas dissected by
    single-cell RNA-seq (Filbin, Tirosh et al., 2018,
    *Science*)](https://doi.org/10.1126/science.aao4750)
-   [Decoupling genetics, lineages, and microenvironment in IDH-mutant
    gliomas by single-cell RNA-seq (Venteicher, Tirosh et al., 2017,
    *Science*)](https://doi.org/10.1126/science.aai8478)

## Installation

You can install the development version of scModules from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("qizongtai/scModules")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
if (!require(devtools)) {
  install.packages("devtools")
}

devtools::install_github("qizongtai/scModules", ref = "development")
library(scModules)
```

## Modules Summary

-   `Module1: Cell type annotation (Marker genes based)`
-   `Module2: Malignant cell identification by copy number aberration`
-   `Module3: Metaprogram (cell state) identification by NMF`
-   `Module4: Immune (TCR/BCR) repertoire analysis`

### Module1

![1_scRNA_celltype_workflow](https://user-images.githubusercontent.com/33009124/177924123-b77d89d4-fc91-4673-8ca3-1823942e7d36.PNG)

### Module2

![2_scRNA pipeline infercna](https://user-images.githubusercontent.com/33009124/177924157-cda90bf3-4953-4c3f-9ab5-7ba03c6222c9.PNG)



infercna aims to provide functions for inferring CNA values from
scRNA-seq data and related queries. 
- `infercna()` to infer copy-number alterations from single-cell RNA-seq data 
- `refCorrect()` to convert relative CNA values to absolute values + computed in `infercna()` if reference cells are provided 
- `cnaPlot()` to plot a heatmap of CNA values 
- `cnaScatterPlot()` to visualise malignant and non-malignant cell subsets 
- `cnaCor()` a parameter to identify cells with high CNAs + computed in `cnaScatterPlot()` 
- `cnaSignal()` a second parameter to identify cells with high CNAs + computed in `cnaScatterPlot()` 
- `findMalignant()` to find malignant subsets of cells 
- `findClones()` to identify genetic subclones 
- `fitBimodal()` to fit a bimodal gaussian distribution + used in `findMalignant()` + used in `findClones()` 
- `filterGenes()` to filter genes by their genome features 
- `splitGenes()` to split genes by their genome features 
- `orderGenes()` to order genes by their genomic position - `useGenome()` to change the default genome configured with infercna 
- `addGenome()` to configure infercna with a new genome specified by the user

### Module3

Functions are ready and need to be intergared into the R package. **to
be continued**

### Module4

Functions are ready and need to be intergared into the R package. **to
be continued**

------------------------------------------------------------------------

## Version History

-   July 08, 2022 Version 0.1.0: Initial release; essential functions
    for streamlined scRNA-seq analysis and visualization
