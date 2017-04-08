griph: Graph Inference of Population Heterogeneity
================
Panagiotis Papasaikas, Michael Stadler
8 April 2017

Introduction
============

Graph Inference of Population Heterogeneity

Sample datasets included in this package
========================================

Buettner et al., Nature Biotechnology 2015
------------------------------------------

Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells. Buettner F, Natarajan KN, Casale FP, Proserpio V, Scialdone A, Theis FJ, Teichmann SA, Marioni JC, Stegle O. Nat Biotechnol. 2015 Feb;33(2):155-60. doi: 10.1038/nbt.3102.

-   Description: 288 mouse ES cells (Fluidigm C1), classified into three cell cycle stages
-   PubMed: <http://www.ncbi.nlm.nih.gov/pubmed/25599176>
-   DOI: <https://doi.org/10.1038/nbt.3102>
-   Data: <http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2805/>
-   Filtering: all cells, top 10,000 genes

``` r
fname <- system.file("extdata", "buettner_top10k.rds", package = "griph")
M <- readRDS(fname)
dim(M)
```

    ## [1] 10000   288

``` r
label <- attr(M, "label")
table(label)
```

    ## label
    ##  G1 G2M   S 
    ##  96  96  96

Kolodziejck et al., Cell Stem Cell 2015
---------------------------------------

NA

-   Description: 704 mouse ES cells (Fluidigm C1), three culture conditions: serum + LIF (lif), 2i + LIF (2i) and alternative 2i + LIF (a2i)
-   PubMed: <https://www.ncbi.nlm.nih.gov/pubmed/26431182>
-   DOI: <https://doi.org/10.1016/j.stem.2015.09.011>
-   Data: <http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-2600/>, <http://www.ebi.ac.uk/teichmann-srv/espresso>

``` r
fname <- system.file("extdata", "kolodziejck_top10k.rds", package = "griph")
M <- readRDS(fname)
dim(M)
```

    ## [1] 10000   704

``` r
label <- attr(M, "label")
table(label)
```

    ## label
    ## lif  2i a2i 
    ## 250 295 159

Usoskin et al., Nat Neurosci. 2015
----------------------------------

NA

-   Description: 799 mouse cells (622 neurons) from lumbar DRGs (custom protocol), different sensory neuron subtypes
-   PubMed: <https://www.ncbi.nlm.nih.gov/pubmed/25420068>
-   DOI: <https://doi.org/10.1038/nn.3881>
-   Data: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59739>

``` r
fname <- system.file("extdata", "usoskin_top10k.rds", package = "griph")
M <- readRDS(fname)
dim(M)
```

    ## [1] 10000   799

``` r
label <- attr(M, "label")
table(label)
```

    ## label
    ## Central, unsolved                NF        NF outlier               NoN 
    ##                39               139                 9               109 
    ##       NoN outlier                NP               PEP                TH 
    ##                 2               169                81               233 
    ##        TH outlier 
    ##                18

``` r
label2 <- attr(M, "label2")
table(label2)
```

    ## label2
    ## Central, unsolved               NF1             NF2/3             NF4/5 
    ##                39                31                60                48 
    ##        NF outlier               NoN       NoN outlier               NP1 
    ##                 9               109                 2               125 
    ##               NP2              PEP1              PEP2                TH 
    ##                44                64                17               233 
    ##        TH outlier 
    ##                18

``` r
label3 <- attr(M, "label3")
table(label3)
```

    ## label3
    ## Central, unsolved               NF1               NF2               NF3 
    ##                39                31                48                12 
    ##               NF4               NF5        NF outlier               NoN 
    ##                22                26                 9               109 
    ##       NoN outlier               NP1               NP2               NP3 
    ##                 2               125                32                12 
    ##              PEP1              PEP2                TH        TH outlier 
    ##                64                17               233                18

Zeisel et al., Science 2015
---------------------------

NA

-   Description: 3005 mouse brain cells (Fluidigm C1), different neuron subtypes
-   PubMed: <https://www.ncbi.nlm.nih.gov/pubmed/25700174>
-   DOI: <https://doi.org/10.1126/science.aaa1934>
-   Data: <http://linnarssonlab.org/cortex/>, <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60361>

``` r
fname <- system.file("extdata", "zeisel_top10k.rds", package = "griph")
M <- readRDS(fname)
dim(M)
```

    ## [1] 10000  3005

``` r
label <- attr(M, "label")
table(label)
```

    ## label
    ## astrocytes_ependymal    endothelial-mural         interneurons 
    ##                  224                  235                  290 
    ##            microglia     oligodendrocytes        pyramidal CA1 
    ##                   98                  820                  939 
    ##         pyramidal SS 
    ##                  399

``` r
label2 <- attr(M, "label2")
table(label2)
```

    ## label2
    ##    Astro1    Astro2   CA1Pyr1   CA1Pyr2 CA1PyrInt   CA2Pyr2   Choroid 
    ##        68        61       380       447        49        41        10 
    ##   ClauPyr     Epend      Int1     Int10     Int11     Int12     Int13 
    ##         5        20        12        21        10        21        15 
    ##     Int14     Int15     Int16      Int2      Int3      Int4      Int5 
    ##        22        18        20        24        10        15        20 
    ##      Int6      Int7      Int8      Int9      Mgl1      Mgl2    (none) 
    ##        22        23        26        11        17        16       189 
    ##    Oligo1    Oligo2    Oligo3    Oligo4    Oligo5    Oligo6     Peric 
    ##        45        98        87       106       125       359        21 
    ##      Pvm1      Pvm2   S1PyrDL  S1PyrL23   S1PyrL4   S1PyrL5  S1PyrL5a 
    ##        32        33        81        74        26        16        28 
    ##   S1PyrL6  S1PyrL6b    SubPyr     Vend1     Vend2      Vsmc 
    ##        39        21        22        32       105        62

Quickstart: A sample analysis
=============================

We will use the mouse ES cell data from Buettner et al. for this example. The dataset is included in the package, including know cell cycle labels:

``` r
M <- readRDS(system.file("extdata", "buettner_top10k.rds", package = "griph"))
dim(M) # genes by cells
```

    ## [1] 10000   288

``` r
label.true <- attr(M, "label")
table(label.true)
```

    ## label.true
    ##  G1 G2M   S 
    ##  96  96  96

Cell types can be identified using :

``` r
library(griph)
res <- SC_cluster(M, ClassAssignment = label.true, plotG = FALSE)
```

    ## Preprocessing...

    ## Calculating Pairwise and Diffused Similarities... 1 / 1

    ## Estimating Graph Structure...

    ## Calculating edge weights and knn-based pruning...

    ## Pruning based on global node similarity...1/1

    ## Keeping...32.7% of edges

    ## Detecting Graph Communities...

    ## Done...

    ## Elapsed Time:  2.911

``` r
table(res$MEMB)
```

    ## 
    ##  1  2  3  4  5  6  7  8  9 
    ## 94 75 32 31 24  8 12  6  6

The automatic generation of plots ( and arguments) is switched off here. identified 9 cell types in the data, and the confusion matrix returned in the result summarizes how they relate to the known cell types or states:

``` r
res$ConfMatrix # confusion matrix
```

    ##          true
    ## predicted G1 G2M  S
    ##         1 13   5 76
    ##         2  2  73  0
    ##         3 28   2  2
    ##         4 29   0  2
    ##         5 24   0  0
    ##         6  0   0  8
    ##         7  0  10  2
    ##         8  0   0  6
    ##         9  0   6  0

``` r
res$miscl # misclassification error
```

    ## [1] 0.09722222

When comparing different classifications with each other, they may differ in granularity, yet be consistent with each other. This can be seen for example in the above confusion matrix for predicted classes *3*, *4* and *5* that all contain mostly *G1* cells. takes this into account when calculating classification error by assigning each identified class to the major known class (*G1* in this example) and only counting the 6 cells in *3*, *4* and *5* as wrongly classified that are not *G1*:

The obtained results can be visualized for example by plotting the graph obtained from :

``` r
plot(res$GRAO)
```

![](https://github.com/ppapasaikas/griph/griph/README_files/figure-markdown_github/buettner_graph-1.png)

or also by applying a dimensionality reduction/projection technique such as tSNE to the affinity matrix returned by :

``` r
library(Rtsne)
out <- Rtsne(as.dist(1-res$DISTM), pca = FALSE, perplexity = 5)
plot(out$Y, col=res$MEMB, pch=as.numeric(label.true))
```

![](https://github.com/ppapasaikas/griph/griph/README_files/figure-markdown_github/buettner_tsne-1.png)

Session info
============

Here is the output of sessionInfo() on the system on which this document was compiled:

    ## R version 3.3.2 (2016-10-31)
    ## Platform: x86_64-apple-darwin13.4.0 (64-bit)
    ## Running under: macOS Sierra 10.12.4
    ## 
    ## locale:
    ## [1] C
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] Rtsne_0.11      BiocStyle_2.2.1 rmarkdown_1.4   griph_0.1.0    
    ## [5] igraph_1.0.1   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.10        bigmemory.sri_0.1.3 knitr_1.15.1       
    ##  [4] magrittr_1.5        nnls_1.4            doParallel_1.0.10  
    ##  [7] coop_0.6-0          foreach_1.4.3       bigmemory_4.5.19   
    ## [10] stringr_1.2.0       tools_3.3.2         QUIC_1.1           
    ## [13] parallel_3.3.2      rNMF_0.5.0          corpcor_1.6.9      
    ## [16] htmltools_0.3.5     iterators_1.0.8     gtools_3.5.0       
    ## [19] yaml_2.1.14         digest_0.6.12       rprojroot_1.2      
    ## [22] codetools_0.2-15    evaluate_0.10       stringi_1.1.5      
    ## [25] backports_1.0.5

