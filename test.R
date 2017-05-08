M <- readRDS(system.file("extdata", "buettner_top10k.rds", package = "griph"))
dim(M) # genes by cells
trueLabel <- attr(M, "label")
table(trueLabel)


for(i in 1:4){
res=griph_cluster(M,ref.iter=0,ClassAssignment = trueLabel, plotG = TRUE, fsuffix='buettner0',use.par=FALSE)
cat("\n\n\n","Miscl:",res$miscl,"\n", "Nclust:", length(unique(res$MEMB)), "\n\n\n")
}