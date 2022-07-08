# scModules

<img align="right" width="108" height="125" src="Scillus.png">

---

[![Build Status](https://travis-ci.com/xmc811/Scillus.svg?branch=master)](https://travis-ci.com/xmc811/Scillus)
[![Build status](https://ci.appveyor.com/api/projects/status/dkq1xn6574kqgs0s/branch/master?svg=true)](https://ci.appveyor.com/project/xmc811/scillus/branch/master)


Scillus is a Seurat and ggplot2 wrapper for enhanced processing and visualization of scRNA-seq data. For detailed usage, please visit the Scillus documentation website: [scillus.netlify.com](http://scillus.netlify.com) [![Netlify Status](https://api.netlify.com/api/v1/badges/eadbcb9a-16d1-4a9a-9e50-c0e8d4104ddc/deploy-status)](https://app.netlify.com/sites/scillus/deploys)


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
