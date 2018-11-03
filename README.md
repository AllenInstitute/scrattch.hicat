# scrattch.hicat: Hierarchical, Iterative Clustering for Analysis of Transcriptomics 

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/c/cb/The_cat_-_an_introduction_to_the_study_of_backboned_animals%2C_especially_mammals_%281881%29_%2820577850762%29.jpg/640px-The_cat_-_an_introduction_to_the_study_of_backboned_animals%2C_especially_mammals_%281881%29_%2820577850762%29.jpg" alt="A hicat" width="300px"/>

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

[An overview of the main functions in `scrattch.hicat`](http://htmlpreview.github.io/?https://github.com/AllenInstitute/scrattch.hicat/blob/master/vignettes/scrattch.hicat.html)  

## Tutorials

[An interactive walkthrough of the major steps in clustering for `scrattch.hicat`.](https://taxonomy.shinyapps.io/scrattch_tutorial/)  

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
