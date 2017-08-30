library(Rcpp)
library(devtools)

# compile attributes for Rcpp code
# will re-create griph/src/RcppExports.cpp and griph/R/RcppExports.R
Rcpp::compileAttributes(pkgdir = "./griph")

# create to for native routine registration (should be put into R_init_griph.c)
# (might eventually be obsolete, onece Rcpp::compileAttributes does this, see: https://github.com/RcppCore/Rcpp/issues/636)
tools::package_native_routine_registration_skeleton("./griph", character_only = FALSE)

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

# build vignette only
#devtools::build_vignettes(pkg = "./griph/") # will build vignette, clean up the ./vignettes directory and copy the output files to ./inst/doc

# ... directly from github
devtools::install_git("git://github.com/ppapasaikas/griph.git", subdir = "griph")
#devtools::install_git("git://github.com/ppapasaikas/griph.git", subdir = "griph", branch = "Testing")

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
vignette(package = "griph")
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

# ... test data from package
list.files(system.file("extdata", package = "griph"))

##### structureScore #####
library(griph)
M <- readRDS(system.file("extdata", "buettner_top10k.rds", package = "griph"))
label <- attr(M, "label")
res  <- SC_cluster(M, ClassAssignment = label, plot = FALSE)

ss <- clusteringScore(M, label, score.type = "sdLogFC") # published labels
ss2 <- clusteringScore(M, res$MEMB, score.type = "sdLogFC") # griph labels
ss3 <- clusteringScore(M, seq.int(ncol(M)), score.type = "sdLogFC") # each cell it's own label
ss4 <- clusteringScore(M, sample(5, ncol(M), replace = TRUE), score.type = "sdLogFC") # random labels (k=5)
ss$score.norm
ss2$score.norm
ss3$score.norm
ss4$score.norm
#plot(density(ss$score.rand), lwd = 2, xlim = c(0.22, 0.34))
#abline(v = ss$score.obs, col = "red")
#lines(density(ss2$score.rand), col="gray")

##### ... ... Buettner #####
M <- readRDS(system.file("extdata", "buettner_top10k.rds", package = "griph"))
label <- attr(M, "label")
res <- griph_cluster(M, ClassAssignment = label, plot_ = FALSE)
plotGraph(res, fill.type = "predicted", line.type = "none")
plotTsne(res, fill.type = "none", line.type = "true")
plotLVis(res, fill.type = label, line.type = "predicted")
plotLVis(res, fill.type = "true", line.type = "none", mark.type = "predicted")
plotLVis(res, fill.type = colSums(M > 0), line.type = "none", mark.type = "predicted")
plotLVis(res, fill.type = colMeans(M > 0), line.type = "true")

M2 <- M
library(org.Mm.eg.db)
eg <- select(x = org.Mm.eg.db, keys = rownames(M2), keytype = "ENSEMBL", columns = "ENTREZID")
eg <- eg[!duplicated(eg$ENSEMBL), ]
M2 <- M2[eg$ENSEMBL,]
rownames(M2) <- as.character(eg$ENTREZID)

cc   <- predictCellCycle(M2, org = "mouse.Whitfield", cor_thr = 0.2, granularity = "low")
ccr  <- predictCellCycle(M2, org = "mouse.Whitfield", cor_thr = 0.2, granularity = "low", refine_iter = 50)
cc2  <- predictCellCycle(M2, org = "mouse.Ishida",    cor_thr = 0.2, granularity = "low")
cc2r <- predictCellCycle(M2, org = "mouse.Ishida",    cor_thr = 0.2, granularity = "low", refine_iter = 50)
table(known = label, predicted = cc)
table(known = label, predicted = ccr)
table(known = label, predicted = cc2)
table(known = label, predicted = cc2r)
griph:::classError(label, cc)
griph:::classError(label, ccr)
griph:::classError(label, cc2)
griph:::classError(label, cc2r)
griph:::classError(label, c("G1.S"="G1", "S"="S", "G2"="G2M", "G2.M"="G2M", "M.G1"="G1")[as.character(cc)])
griph:::classError(label, c("G1.S"="G1", "S"="S", "G2"="G2M", "G2.M"="G2M", "M.G1"="G1")[as.character(ccr)])
chisq.test(table(label, cc))
chisq.test(table(label, ccr))
chisq.test(table(label, cc2))
chisq.test(table(label, cc2r))

clusteringScore(M, label)$score.norm
clusteringScore(M, res$MEMB)$score.norm
clusteringScore(M, cc)$score.norm
clusteringScore(M, ccr)$score.norm
clusteringScore(M, cc2)$score.norm
clusteringScore(M, cc2r)$score.norm

res  <- SC_cluster(M, ClassAssignment = label, plot = FALSE)
resP <- SC_cluster(M, ClassAssignment = label, plot = FALSE, use.par = TRUE, ncores = 8)
for(nm in names(res))
    cat(nm,":",identical(res[[nm]], resP[[nm]]), "\n")
all.equal(res$DISTM, resP$DISTM)
isomorphic(res$GRAO, resP$GRAO)
isomorphic(res$plotGRAO, resP$plotGRAO)

res$miscl # 0.1736111
res$ConfMatrix

griph:::classError(res$MEMB, label)
griph:::classError(res$MEMB, cc)
griph:::classError(res$MEMB, c("G1.S" = "G1", "S" = "S", "G2" = "G2M",
                               "G2.M" = "G2M", "M.G1" = "G1")[as.character(cc)])

res <- griph_cluster(M, ClassAssignment = label, plot = FALSE)
res$miscl


g <- plotGraph(res)
g <- plotGraph(res, fill.type = "custom", custom.class = cc, fill.col = "Set1")
g <- plotGraph(res, mark.type = "true")
g <- plotGraph(res, mark.type = "predicted")
g <- plotGraph(res, mark.type = "custom", custom.class = cc)
g <- plotGraph(res, fill.type = "predicted", line.type = "none")
g <- plotGraph(res, collapse.type = "true")
g <- plotGraph(res, collapse.type = "predicted")

res2 <- griph_cluster(M, ref.iter = 1, ClassAssignment = label, plot = FALSE)
g2 <- plotGraph(res2)
g2 <- plotGraph(res2, line.type = "true", fill.type = "custom", fill.col = "Set1",
                custom.class = cc, draw.edges = TRUE)

library(largeVis)
res3 <- largeVis(res$DISTM, dim = 2, K = 20, save_edges = FALSE)
plot(t(res3$coords), col = RColorBrewer::brewer.pal(length(unique(res$MEMB)), "Set1")[res$MEMB])
str(res3$coords)

res4 <- projectKNNs(Matrix(res$DISTM), dim = 2, verbose = TRUE, gamma = 10, alpha = 0.02, rho = 0.5)
plot(t(res4), col = RColorBrewer::brewer.pal(length(unique(res$MEMB)), "Set1")[res$MEMB])

# nice graph visualization examples: https://rpubs.com/kateto/netviz
g2 <- igraph::simplify( igraph::contract(res$GRAO, res$MEMB) ) # just plot one vertex per cell types
plot(g2)

par(mfrow = c(1,2)); g.true <- plotGraph(res, fill.type = "true", line.type = "none"); g.pred <- plotGraph(res, fill.type = "pred", line.type = "none")
g <- plotGraph(res)
g <- plotGraph(res, forceRecalculation = TRUE, maxG = 50)
# g <- plotGraph(res, image.format = NA)
# g <- plotGraph(res, image.format = "png")
# igraph::plot.igraph(g, asp=0, vertex.label=NA, edge.lty=0)



##### ... ... Kolodziejck #####
library(griph)
M <- readRDS(system.file("extdata", "kolodziejck_top10k.rds", package = "griph"))
label <- attr(M, "label")

M <- M[grep("^ENSM",rownames(M)),] # remove diagnostic/ERCC rows

M2 <- M
library(org.Mm.eg.db)
eg <- select(x=org.Mm.eg.db, keys=rownames(M2), keytype = "ENSEMBL", columns = "ENTREZID")
eg <- eg[!duplicated(eg$ENSEMBL), ]
M2 <- M2[eg$ENSEMBL,]
rownames(M2) <- as.character(eg$ENTREZID)

cc1 <- predictCellCycle(M2, org="mouse.Whitfield", cor_thr = 0.3, granularity = "low")
cc2 <- predictCellCycle(M2, org="mouse.Ishida", cor_thr = 0.3, granularity = "low")
table(cc1, cc2)
griph:::classError(cc1, cc2)
cc <- cc2

res <- SC_cluster(M2, ClassAssignment = label, plot = FALSE, use.par = TRUE, ncores = 8)
res2 <- SC_cluster(M2, ClassAssignment = label, BatchAssignment = cc, plot = FALSE, use.par = TRUE, ncores = 8)

par(mfrow=c(1,2))
g <- plotGraph(res)
g2 <- plotGraph(res2)

par(mfrow=c(1,2))
g <- plotGraph(res, fill.type = "predicted", line.type = "none", mark.type = "true")
g <- plotGraph(res, fill.type = "custom",    line.type = "none", fill.col = "Set1", custom.class = cc)

par(mfrow=c(1,2))
g3 <- plotGraph(res, collapse.type = "predicted")
g4 <- plotGraph(res2, collapse.type = "predicted")

par(mfrow=c(1,2))
t1 <- plotTsne(res, fill.type = "predicted", line.type = "none", mark.type = "true")
t2 <- plotTsne(res, fill.type = "custom",    line.type = "none", fill.col = "Set1", custom.class = cc)

table(res$MEMB, cc)
table(res2$MEMB, cc)
table(res$MEMB, res2$MEMB)

res$miscl # 0
res$ConfMatrix

par(mfrow=c(1,2))
plotGraph(res, fill.type = "true", line.type = "none")
plotGraph(res, fill.type = "pred", line.type = "none")

wc <- cluster_infomap(res$GRAO)
wc$membership
plot(wc, res$GRAO)




##### ... ... Usoskin #####
M <- readRDS(system.file("extdata", "usoskin_top10k.rds", package = "griph"))
label <- attr(M, "label")
#label <- attr(M, "label2")
#label <- attr(M, "label3")
res <- SC_cluster(M, ClassAssignment = label)
res$miscl # 0.1451815
griph:::classError(res$MEMB, attr(M, "label"))  # 0.1451815
griph:::classError(res$MEMB, attr(M, "label2")) # 0.2941176
griph:::classError(res$MEMB, attr(M, "label3")) # 0.3128911




##### ... ... Zeisel #####
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
# library(microbenchmark)
# 
# set.seed(0)
# M1 <- matrix(rpois(5000*300, 8), ncol=300)
# M2 <- matrix(rpois(5000*300, 8), ncol=300)
# bres <- microbenchmark(
#     { res1 <- griph:::PHellingerMat(M1, M2) },
#     { res2 <- griph:::PHellingerMatOMP(M1, M2) },
#     { res3 <- griph:::PHellingerMatOMP(M1, M2, 1) },
#     { res4 <- griph:::PHellingerMatOMP(M1, M2, 2) },
#     { res5 <- griph:::PHellingerMatOMP(M1, M2, 4) },
#     { res5 <- griph:::PHellingerMatOMP(M1, M2, 6) },
#     { res6 <- griph:::PHellingerMatOMP(M1, M2, 8) },
#     times = 5
# )
# bres <- microbenchmark(
#     { res1 <- griph:::PCanberraMat(M1, M2) },
#     { res2 <- griph:::PCanberraMatOMP(M1, M2) },
#     { res3 <- griph:::PCanberraMatOMP(M1, M2, 1) },
#     { res4 <- griph:::PCanberraMatOMP(M1, M2, 2) },
#     { res5 <- griph:::PCanberraMatOMP(M1, M2, 4) },
#     { res5 <- griph:::PCanberraMatOMP(M1, M2, 6) },
#     { res6 <- griph:::PCanberraMatOMP(M1, M2, 8) },
#     times = 5
# )
# bres <- microbenchmark(
#     { res1 <- cor(M1, M2) },
#     { res2 <- griph:::PPearsonMatOMP(M1, M2) },
#     { res3 <- griph:::PPearsonMatOMP(M1, M2, 1) },
#     { res4 <- griph:::PPearsonMatOMP(M1, M2, 2) },
#     { res5 <- griph:::PPearsonMatOMP(M1, M2, 4) },
#     { res5 <- griph:::PPearsonMatOMP(M1, M2, 6) },
#     { res6 <- griph:::PPearsonMatOMP(M1, M2, 8) },
#     times = 3
# )
# cl <- makeCluster(4)
# registerDoParallel(cl)
# bres <- microbenchmark(
#     #{ registerDoParallel(cl) }, # negletable, takes ~0.1 millisec
#     { registerDoParallel(cl[1]); res1 <- griph:::FlashPHellinger(M1, M2) },
#     { res2 <- griph:::FlashPHellingerOMP(M1, M2, 1) },
#     { registerDoParallel(cl[1:2]); res3 <- griph:::FlashPHellinger(M1, M2) },
#     { res4 <- griph:::FlashPHellingerOMP(M1, M2, 2) },
#     { registerDoParallel(cl[1:4]); res5 <- griph:::FlashPHellinger(M1, M2) },
#     { res6 <- griph:::FlashPHellingerOMP(M1, M2, 4) },
#     times = 5
# )
# bres <- microbenchmark(
#     { registerDoParallel(cl[1]); res1 <- griph:::FlashPCanberra(M1, M2) },
#     { res2 <- griph:::FlashPCanberraOMP(M1, M2, 1) },
#     { registerDoParallel(cl[1:2]); res3 <- griph:::FlashPCanberra(M1, M2) },
#     { res4 <- griph:::FlashPCanberraOMP(M1, M2, 2) },
#     { registerDoParallel(cl[1:4]); res5 <- griph:::FlashPCanberra(M1, M2) },
#     { res6 <- griph:::FlashPCanberraOMP(M1, M2, 4) },
#     times = 3
# )
# bres <- microbenchmark(
#     { registerDoParallel(cl[1]); res1 <- griph:::FlashPPearsonCor(M1, M2) },
#     { res2 <- griph:::FlashPPearsonCorOMP(M1, M2, 1) },
#     { registerDoParallel(cl[1:2]); res3 <- griph:::FlashPPearsonCor(M1, M2) },
#     { res4 <- griph:::FlashPPearsonCorOMP(M1, M2, 2) },
#     { registerDoParallel(cl[1:4]); res5 <- griph:::FlashPPearsonCor(M1, M2) },
#     { res6 <- griph:::FlashPPearsonCorOMP(M1, M2, 4) },
#     times = 3
# )
# bres <- microbenchmark(
#     { registerDoParallel(cl[1]); res1 <- griph:::FlashPSpearmanCor(M1, M2) },
#     { res2 <- griph:::FlashPSpearmanCorOMP(M1, M2, 1) },
#     { registerDoParallel(cl[1:2]); res3 <- griph:::FlashPSpearmanCor(M1, M2) },
#     { res4 <- griph:::FlashPSpearmanCorOMP(M1, M2, 2) },
#     { registerDoParallel(cl[1:4]); res5 <- griph:::FlashPSpearmanCor(M1, M2) },
#     { res6 <- griph:::FlashPSpearmanCorOMP(M1, M2, 4) },
#     times = 3
# )
# identical(res1, res2) # summary(as.vector(res1 - res2))
# identical(res1, res3)
# identical(res1, res4)
# identical(res1, res5)
# identical(res1, res6)
# bres
# sbres <- summary(bres)
# sbres[1,"min"] / sbres[4,"min"] # speed-up
# plot(c(NA,NA,1,2,4,6,8), sbres[1,"min"] / sbres[,"min"], xlab="No. threads", ylab="Speedup"); abline(a=0, b=1)

library(microbenchmark)
set.seed(0)
M1 <- as(matrix(sample(0:10, 5000*300, TRUE, c(.5, rep(.05,10))), ncol=300), "sparseMatrix")

microbenchmark(
    { res1 <- cov(as.matrix(M1)) },
    { res2 <- SPearsonMatOMP(M1, 1) },
    { res3 <- SPearsonMatOMP(M1, 2) },
    times = 3
)
identical(res1, res2)
identical(res1, res3)


library(microbenchmark)
set.seed(0)
M1 <- matrix(sample(0:10, 5000*300, TRUE, c(.5, rep(.05,10))), ncol=300)
M2 <- matrix(sample(0:10, 5000*300, TRUE, c(.5, rep(.05,10))), ncol=300)
S1 <- as(M1, "sparseMatrix")
S2 <- as(M2, "sparseMatrix")

microbenchmark(
    { res1 <- PCanberraMatOMP(M1, M2, 2) },
    { res2 <- griph:::ssPCanberraMatOMP(S1, S2, 2) }, # ~20% slower than dense
    { res3 <- griph:::sdPCanberraMatOMP(S1, M2, 2) }, # ~10% faster than dense
    { res4 <- t(griph:::sdPCanberraMatOMP(S2, M1, 2)) },
    times = 5
)
identical(res1, res2)
identical(res1, res3)
identical(res1, res4)
identical(res3, res4)

# griph:::PCanberraMatOMP  (M1[1:4,1,drop=FALSE], M2[1:4,1:2], 1)
# griph:::ssPCanberraMatOMP(S1[1:4,1,drop=FALSE], S2[1:4,1:2], 1)
# griph:::sdPCanberraMatOMP(S1[1:4,1,drop=FALSE], M2[1:4,1:2], 1)
# t(griph:::sdPCanberraMatOMP(S2[1:4,1:2], M1[1:4,1,drop=FALSE], 1))
# S1[1:4,1,drop=FALSE]
# S2[1:4,1:2]
# j=1; k=1; abs(M1[1:4,j]-M2[1:4,k])/(M1[1:4,j]+M2[1:4,k])
# j=1; k=2; abs(M1[1:4,j]-M2[1:4,k])/(M1[1:4,j]+M2[1:4,k])

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
