[![Travis-CI Build Status](https://travis-ci.org/ppapasaikas/griph.svg?branch=master)](https://travis-ci.org/ppapasaikas/griph) 
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ppapasaikas/griph?branch=master&svg=true)](https://ci.appveyor.com/project/ppapasaikas/griph?branch=master)

# griph: Graph Inference of Population Heterogeneity
Panagiotis Papasaikas, Michael Stadler, Atul Sethi  
April 2017  

*griph* (Graph Inference of Population Heterogeneity) is an R package for the analysis of
single cell RNA-sequencing data. It can be used to automatically identify different cell
types or states, even in the presence of confounding sources of variance such as
cell cycle stages or batch effects.

## Installation
*griph* is free and currently available through https://github.com and can be installed using:  

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
