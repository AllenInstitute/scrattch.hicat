# hicat
Iterative clustering of transcriptomic data in R
## Installation

hicat has several dependencies, including one from BioConductor:
```
source("https://bioconductor.org/biocLite.R")
biocLite("WGCNA")
biocLite("limma")
library(devtools)
devtools::install_github("JinmiaoChenLab/Rphenograph")
devtools::install_github("AllenInstitute/hicat",auth_token = "96cb6605b9e7d395b5b3e3e3c04f0eb001cf4674")
```
