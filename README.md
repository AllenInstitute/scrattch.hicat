# scrattch.hicat: Hierarchical, Iterative Clustering for Analysis of Transcriptomics 

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/The_cat_-_an_introduction_to_the_study_of_backboned_animals%2C_especially_mammals_%281881%29_%2820577850762%29.jpg/640px-The_cat_-_an_introduction_to_the_study_of_backboned_animals%2C_especially_mammals_%281881%29_%2820577850762%29.jpg" alt="A hicat" width="300px"/>

Master: [![Travis build status](https://travis-ci.org/AllenInstitute/scrattch.hicat.svg?branch=master)](https://travis-ci.org/AllenInstitute/scrattch.hicat)
[![Coverage status](https://codecov.io/gh/AllenInstitute/scrattch.hicat/branch/master/graph/badge.svg)](https://codecov.io/github/AllenInstitute/scrattch.hicat?branch=master)

Dev: [![Travis build status](https://travis-ci.org/AllenInstitute/scrattch.hicat.svg?branch=dev)](https://travis-ci.org/AllenInstitute/scrattch.hicat)
[![Coverage status](https://codecov.io/gh/AllenInstitute/scrattch.hicat/branch/dev/graph/badge.svg)](https://codecov.io/github/AllenInstitute/scrattch.hicat?branch=dev)

## Installation

`scrattch.hicat` has several dependencies, including two from BioConductor and one from Github:
```
source("https://bioconductor.org/biocLite.R")
biocLite("limma")

devtools::install_github("JinmiaoChenLab/Rphenograph")
```

Once these dependencies are installed, `scrattch.hicat` can be installed with:
```
devtools::install_github("AllenInstitute/scrattch.hicat")
```

## Vignettes

[An overview of the main functions in `scrattch.hicat`](http://htmlpreview.github.io/?https://github.com/AllenInstitute/scrattch.hicat/blob/master/vignettes/scrattch.hicat_release.html)  

## Tutorials

[An interactive walkthrough of the major steps in clustering for `scrattch.hicat`.](https://taxonomy.shinyapps.io/scrattch_tutorial/)  

## Roadmap

The next few updates to `scrattch.hicat` will be aimed at getting code testing in place for major clustering functions:  
0.0.22: Current version; Tests in place for de.genes.R functions.  
0.0.23: Tests in place for cluster.R functions.  
0.1.0: Vignette re-integrated; Adding pkgdown page; Update to Master branch.

Previous updates:
0.0.21: Added TravisCI and `covr` integration.

## The `scrattch` suite

`scrattch.hicat` is one component of the [scrattch](https://github.com/AllenInstitute/scrattch/) suite of packages for Single Cell RNA-seq Analysis for Transcriptomic Type CHaracterization from the Allen Institute.

## License

The license for this package is available on Github at: https://github.com/AllenInstitute/scrattch.hicat/blob/master/LICENSE

## Level of Support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.

## Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/scrattch.hicat/blob/master/CONTRIBUTION

#### Image attribution:
By Internet Archive Book Images [No restrictions], via Wikimedia Commons
