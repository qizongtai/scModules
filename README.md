# scModules

<img align="right" width="108" height="125" src="Scillus.png">

---

[![Build Status](https://travis-ci.com/xmc811/Scillus.svg?branch=master)](https://travis-ci.com/xmc811/Scillus)
[![Build status](https://ci.appveyor.com/api/projects/status/dkq1xn6574kqgs0s/branch/master?svg=true)](https://ci.appveyor.com/project/xmc811/scillus/branch/master)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)


scModules is an R package for advanced processing and visualization of single-cell RNA-seq and TCR-seq data.  
It provides advanced analyses that NOT available in the popular single-cell package Seurat and currently has four functional modules as below. The methodology behind modules has been tried and tested in several high-impact publications, including Cell and Science. scModules enhances our ability to extract accurate and meaningful information from single-cell data and discover interesting cancer biology

- `Module1: Cell type annotation (Marker genes based)`
- `Module2: Malignant cell identification by copy number aberration`
- `Module3: Metaprogram (cell state) identification`
- `Module4: Immune (TCR/BCR) repertoire analysis`

Please use the following code to install and load the package:

```R
if (!require(devtools)) {
  install.packages("devtools")
}

devtools::install_github("qizongtai/scModules", ref = "development")
library(scModules)
```

---

### Version History

* August 25, 2021

  Version 0.1.0: Initial release; essential functions for streamlined scRNA-seq analysis and visualization
