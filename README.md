[![Travis-CI Build Status](https://travis-ci.org/ppapasaikas/griph.svg?branch=master)](https://travis-ci.org/ppapasaikas/griph) 
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ppapasaikas/griph?branch=master&svg=true)](https://ci.appveyor.com/project/ppapasaikas/griph?branch=master)

# griph: Graph Inference of Population Heterogeneity
Panagiotis Papasaikas, Michael Stadler, Atul Sethi  
April 2017  

*griph* (Graph Inference of Population Heterogeneity) is an R package for the analysis of
single cell RNA-sequencing data. It can be used to automatically identify different cell
types or states, even in the presence of confounding sources of variance such as
cell cycle stages or batch effects.

## System Requirements
In order to use *griph*, you need:  

- R (www.r-project.org, tested with version >= 3.4)  
- R packages: igraph, QUIC, coop, corpcor, foreach, doParallel, parallel, bigmemory, RColorBrewer, Rcpp (>= 0.12.9), RcppProgress (>= 0.2.1), RcppArmadillo (>= 0.7.200.2.0)
- A C++ compiler that supports C++11 (e.g. gcc 4.8 or clang 3.3)

To benefit from OpenMP-based parallelization, your compiler must support OpenMP (optional, e.g. gcc 4.8 or clang 4.0)

In order to run *griph*, especially when using large single cell datasets (more than
5,000 cells), we recommend available memory (RAM) of at least 8 GB.

For a dataset of about 10,000 cells and 20,000 genes, *griph* takes about three minutes
on a standard desktop computer with 4 CPU cores and a total memory of 8 GB.

## Installation
*griph* is free (GPL3) and currently available through https://github.com. It can be installed from R using:  

```r
library(devtools)
install_git("git://github.com/ppapasaikas/griph.git", subdir = "griph")
#Or on a windows machine substitute the previous line with:
install_git("https://github.com/ppapasaikas/griph.git", subdir = "griph")
```

## Documentation and Examples
... including step-by-step analysis examples of multiple datasets provided by *griph*
are available from the vignette: (https://ppapasaikas.github.io/griph/)  
and from within R:  

```r
library(griph)
vignette(griph)
```

![Buettner](https://raw.githubusercontent.com/ppapasaikas/griph/master/griph_example.png)
