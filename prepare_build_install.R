library(Rcpp)
library(devtools)

# compile attributes for Rcpp code
# will re-create griph/src/RcppExports.cpp and griph/R/RcppExports.R
Rcpp::compileAttributes(pkgdir = "./griph")

# create to for native routine registration (should be put into R_init_griph.c)
# (might eventually be obsolete, onece Rcpp::compileAttributes does this, see: https://github.com/RcppCore/Rcpp/issues/636)
tools::package_native_routine_registration_skeleton("./griph")

# re-create documentation using roxygen2::roxygenize
# see: browseVignettes("roxygen2")
devtools::document(pkg = "./griph") # will not overwrite existing NAMESPACE
#roxygen2::roxygenize(package.dir = "./griph") # crashes

# build package tar-ball
devtools::clean_dll(pkg = "./griph") # clean compiled objects from /src/ directory
devtools::build(pkg = "./griph")
devtools::build(pkg = "./griph", vignettes = FALSE)

# install package
devtools::install(pkg = "./griph")
devtools::install(pkg = "./griph", build_vignettes = TRUE)
devtools::install(pkg = "./griph", reload = TRUE, quick = TRUE)

# ... directly from github
devtools::install_git("git://github.com/ppapasaikas/griph.git", subdir = "griph")
devtools::install_git("git://github.com/ppapasaikas/griph.git", subdir = "griph", branch = "Testing")

# setup unit tests with testthat
# This will create ‘tests/testthat.R’, ‘tests/testthat/’ and add testthat to the suggested packages.
if (!file.exists("./griph/tests"))
  devtools::use_testthat(pkg = "./griph")

# check package
# library(BiocCheck)
# BiocCheck(package = "./griph") # currently does not work
devtools::test(pkg = "./griph") # automatically reloads code/tests
devtools::check(pkg = "./griph")

# run examples
devtools::run_examples(pkg = "./griph")

# play with griph
library(griph)
vignette(package="griph")
help(package = "griph")
package?griph
?HellingerMat
?sparsify
example(HellingerMat)

x <- matrix(0:8, nrow = 3)
format(sum(HellingerMat(x)), digits = 10)
sum(is.finite(HellingerMat(x - 3)))

# ... GriphResult methods
M <- readRDS(system.file("extdata", "buettner_top10k.rds", package = "griph"))
label <- attr(M, "label")
res <- SC_cluster(M, ClassAssignment = label)
getClass("GriphResult")
names(res)
obj <- new(Class = "GriphResult", DM = M, is.cor = FALSE, ClassAssignment = res$MEMB.true,
           BatchAssignment = factor(),
           MEMB = res$MEMB, DISTM = res$DISTM, ConfMatrix = res$ConfMatrix,
           miscl = res$miscl, GRAO = res$GRAO, plotGRAO = res$plotGRAO)
obj
show(obj)
slotNames(obj)
slot(obj, "ConfMatrix")
slot(obj, "ClassAssignment")

# ... test data from pacakge
list.files(system.file("extdata", package = "griph"))

# ... ... Buettner
M <- readRDS(system.file("extdata", "buettner_top10k.rds", package = "griph"))
label <- attr(M, "label")
res <- SC_cluster(M, ClassAssignment = label)
res$miscl # 0.09722222
g <- plotGraph(res)
g <- plotGraph(res, mark.type = "true")
g <- plotGraph(res, mark.type = "predicted")
g <- plotGraph(res, fill.type = "predicted", line.type = "none")
g <- plotGraph(res, collapse.type = "true")
g <- plotGraph(res, collapse.type = "predicted")

# nice graph visualization examples: https://rpubs.com/kateto/netviz
g2 <- igraph::simplify( igraph::contract(res$GRAO, res$MEMB) ) # just plot one vertex per cell types
plot(g2)

par(mfrow = c(1,2)); g.true <- plotGraph(res, fill.type = "true", line.type = "none"); g.pred <- plotGraph(res, fill.type = "pred", line.type = "none")
g <- plotGraph(res)
g <- plotGraph(res, forceRecalculation = TRUE, maxG = 50)
# g <- plotGraph(res, image.format = NA)
# g <- plotGraph(res, image.format = "png")
# igraph::plot.igraph(g, asp=0, vertex.label=NA, edge.lty=0)

# ... ... Kolodziejck
M <- readRDS(system.file("extdata", "kolodziejck_top10k.rds", package = "griph"))
label <- attr(M, "label")
res <- SC_cluster(M, ClassAssignment = label)
res$miscl # 0

par(mfrow=c(1,2))
plotGraph(res, fill.type = "true", line.type = "none")
plotGraph(res, fill.type = "pred", line.type = "none")

wc <- cluster_infomap(res$GRAO)
wc$membership
plot(wc, res$GRAO)

# ... ... Usoskin
M <- readRDS(system.file("extdata", "usoskin_top10k.rds", package = "griph"))
label <- attr(M, "label")
#label <- attr(M, "label2")
#label <- attr(M, "label3")
res <- SC_cluster(M, ClassAssignment = label)
res$miscl # 0.1451815
griph:::classError(res$MEMB, attr(M, "label"))  # 0.1451815
griph:::classError(res$MEMB, attr(M, "label2")) # 0.2941176
griph:::classError(res$MEMB, attr(M, "label3")) # 0.3128911

# ... ... Zeisel
M <- readRDS(system.file("extdata", "zeisel_top10k.rds", package = "griph"))
label <- attr(M, "label")
#label <- attr(M, "label2")
#res <- SC_cluster(M, ClassAssignment = label) # not run yet - takes a long time
#res$miscl # 
griph:::classError(res$MEMB, attr(M, "label"))  #
griph:::classError(res$MEMB, attr(M, "label2")) #

res$ConfMatrix




# ... test data
load(file="/work/gbioinfo/papapana/FMI_groups/SingleCell_Analysis/Data/Test_1_mECS.RData",verbose=TRUE)    # Buettneret al Nat. Biotechn. 2015: Embryonic stem cell under different cell cycle stages, http://www.ncbi.nlm.nih.gov/pubmed/25599176
load(file="/work/gbioinfo/papapana/FMI_groups/SingleCell_Analysis/Data/Test_2_Kolod.RData",verbose=TRUE)   # Kolodziejck et. al Cell Stem Cell 2015: pluripotent cells under different environment conditions, http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4595712/
load(file="/work/gbioinfo/papapana/FMI_groups/SingleCell_Analysis/Data/Test_3_Pollen.RData",verbose=TRUE)  # Pollen et. al Botechn, 2014. 11 classes. Distinct cell types, http://www.ncbi.nlm.nih.gov/pubmed/25086649
load(file="/work/gbioinfo/papapana/FMI_groups/SingleCell_Analysis/Data/Test_4_Usoskin.RData",verbose=TRUE) # Usoskin et. al Nat NeuroSc, 2014. Neuronal cells with sensory subtypes, http://www.ncbi.nlm.nih.gov/pubmed/25420068
load(file="/work/gbioinfo/papapana/FMI_groups/SingleCell_Analysis/Data/Zelsel.RData",verbose=TRUE)         # Zeisel et al., Science, 2015. Brain structure. Cell types in the mouse cortex and hippocampus revealed by single-cell RNA-seq, https://www.ncbi.nlm.nih.gov/pubmed/25700174

Data=Test_1_mECS    # Fluidigm C1, 182 cells (<5s)
Data=Test_2_Kolod   # Smart-seq2,  704 cells (~1min)
Data=Test_3_Pollen  # 249 cells (~8s)
Data=Test_4_Usoskin # 622 cells (~1min)
Data=Zelsel         # 3005 cells (~82min)

M=Data$in_X         # Load Count Data. These are NOT raw counts.
c=Data$n_clust      # Number of CLusters
ClassAssignment=as.vector(unlist(Data$true_labs)) # True Assignment
ncol(M)
summary(as.vector(M))

## res <- SC_cluster(M, ClassAssignment = ClassAssignment)
res <- SC_cluster(2^M, ClassAssignment = ClassAssignment) # M is in log space (assume base 2)

names(res)
res$ConfMatrix
res$miscl

# benchmark
#library(microbenchmark)
#x <- matrix(rpois(5000*300, 5), ncol = 300) # 5000 genes by 300 cells
#microbenchmark(
#    { res1 <- griph::HellingerMat(t(x)); res1 <- res1 + t(res1) },
#    { res2 <- griph::HellingerMatColumns(x); res2 <- res2 + t(res2) },
#    { res3 <- griph::HellingerMatColumns2(x) },
#    times = 10
#)

#library(microbenchmark)
#x <- matrix(rpois(2000*200, 5), ncol = 200) # 2000 genes by 200 cells
#res1 <- griph:::JSDmat(t(x)); res1 <- res1 + t(res1)
#res2 <- griph:::JSDmatColumns(x)
#identical(res1, res2) # TRUE
#microbenchmark(
#    { res1 <- griph:::JSDmat(t(x)); res1 <- res1 + t(res1) },
#    { res2 <- griph:::JSDmatColumns(x) },
#    times = 3
#)
