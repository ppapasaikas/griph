if(getRversion() >= "2.15.1")  utils::globalVariables(c("i","M")) # let codetools know about "i" and "M" in foreach()


# Fast, memory efficient computation of Shrinked Correlation Estimates for big matrices using BLAS implementations, parallelization and bigmemory objects 
# More optimized BLAS versions (ATLAS, openBLAS) should offer further improvements.
# Input is a matrix, the desired number of cores and the lambda parameter controlling shrinkage intensity
FlashShrinkCor <- function (DM, ncores=NULL, lambda=0, w=rep(1,nrow(DM)), verbose=FALSE) { 
  DM   <- as.big.matrix(DM)
  desc <- bigmemory::describe(DM)
  NCOL <- desc@description["totalCols"][[1]]
  
  if (!is.null(ncores))
      warning("Ignoring 'ncores' argument - using the value from getDoParWorkers().")
  ncores <- getDoParWorkers()
  nblocks <- ncores 

  corMAT <- matrix(data = 0, nrow = NCOL, ncol = NCOL)
  corMAT <- as.big.matrix(corMAT, type="double")
  bigCOR <- bigmemory::describe(corMAT)
  ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
  SPLIT <- split(1:NCOL, ceiling(seq_along(1:NCOL) / ceiling(NCOL / nblocks)))
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  if (nrow(COMBS)>1) { COMBS <- COMBS[-which(COMBS[,1] == COMBS[,2]),] }
  COMBS <- unique(COMBS)
  ## iterate through each block combination, calculate correlation matrix between blocks and store them:
  results <- foreach(i = 1:nrow(COMBS), .export=c('cor.shrink'), .packages=('bigmemory')) %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    x <- bigmemory::attach.big.matrix(desc)
    corMAT <- bigmemory::attach.big.matrix(bigCOR)
    v <- unique(c(G1, G2))
    COR <- cor.shrink(x[,v], verbose = verbose, lambda = lambda, w = w)
    corMAT[v, v] <- COR
    COR <- NULL
    return(0L)
  }
  return(as.matrix(corMAT))
}




# Parallel regularized estimates of the inverse covarinace matrix using the QUIC implementation.
FlashGlasso <- function (DM, ncores=NULL, rho=0.3, tol=1e-3, maxIter=100, msg=0) {
  DM   <- as.big.matrix(DM)
  desc <- describe(DM)
  NCOL <- desc@description["totalCols"][[1]]
  NROW <- desc@description["totalRows"][[1]]
  
  if (!is.matrix(rho) && length(rho) != 1 && length(rho) != NCOL)  { stop("Wrong number of elements in rho") }
  if (is.vector(rho)) { rho <- matrix(sqrt(rho)) %*% sqrt(rho) }
  if (length(rho) == 1) {rho <- matrix(rho, ncol = NCOL, nrow = NCOL)}

  if (!is.null(ncores))
      warning("Ignoring 'ncores' argument - using the value from getDoParWorkers().")
  ncores <- getDoParWorkers()
  nblocks <- 4*ncores 
  
  PrecMAT <- matrix(data = 0, nrow = NCOL, ncol = NCOL)
  PrecMAT <- as.big.matrix(PrecMAT, type = "double")
  bigPREC <- describe(PrecMAT)
  ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
  SPLIT <- split(1:NCOL, ceiling(seq_along(1:NCOL) / ceiling(NCOL / nblocks)))
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  if (nrow(COMBS)>1) { COMBS <- COMBS[-which(COMBS[,1] == COMBS[,2]),] }
  COMBS <- unique(COMBS)

  ## iterate through each block combination, calculate sparse inv. cov. matrix between blocks and store them:
  results <- foreach(i = 1:nrow(COMBS), .export=c('QUIC'), .packages=('bigmemory'), .combine='c') %dopar% {
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    x <- attach.big.matrix(desc)
    PrecMAT <- attach.big.matrix(bigPREC)
    v <- unique(c(G1, G2))
    r <- rho[v, v]
    PREC <- QUIC(x[v, v], rho = r, maxIter = maxIter, tol = tol, msg = msg)$X
    PrecMAT[v, v] <- PREC
    PREC <- NULL
    return(0L)
  }
  return(as.matrix(PrecMAT))
}




#Parallel computation of personalized-pagerank based distances
FlashPPR <- function (G, ncores=NULL, df=0.75) { 
  if (!isTRUE(is.igraph(G))) {  
    if (!isSymmetric(G) )  { stop("!! G should be either a graph object or a symmetric matrix !!") }
    G <- igraph::graph.adjacency(G[1:nrow(G), 1:nrow(G)], mode = "max", weighted = TRUE, diag = FALSE)
  }
  
  NCOL <- length(V(G))
  L <- length(V(G))

  if (!is.null(ncores))
      warning("Ignoring 'ncores' argument - using the value from getDoParWorkers().")
  ncores <- getDoParWorkers()
  nblocks <- ncores 

  PR <- diag(nrow = L)
  ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
  SPLIT <- split(1:L, ceiling(seq_along(1:L) / ceiling(L / nblocks)))
  ## iterate through each block combination, calculate kernel distance between blocks and store them:
  results <- foreach(i = 1:length(SPLIT), .export=c('page_rank'), .combine='c') %dopar% {
    vids <- SPLIT[[i]]
    vals <- sapply(vids, function(x) page_rank(G, vids = c(1:L), personalized = PR[,x], damping=df)$vector)
    vals <- log2(unlist(vals) + (0.01 / L))
  }
  PR[1:length(PR)] <- results
  PR[lower.tri(PR)] <- (PR[lower.tri(PR)] + t(PR)[lower.tri(t(PR))]) / 2
  PR[upper.tri(PR)] <- t(PR)[upper.tri(t(PR))]
  PR <- 2^PR
  return(PR)
}




# Fastcomputation of Pearson's correlation between two big matrices 
# Inputs are the two matrices and the desired number of cores
# The function assumes that matrix with the largest dimension is DM2
FlashPPearsonCor <- function (DM1, DM2, ncores=NULL) { 
    if (!is.null(ncores))
        warning("Ignoring 'ncores' argument - using the value from getDoParWorkers().")
    ncores <- getDoParWorkers()
    nblocks <- 2*ncores 
    ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
    SPLIT <- split(1:ncol(DM2), ceiling(seq_along(1:ncol(DM2)) / ceiling(ncol(DM2) / nblocks)))
    DM2.list <- lapply(1:length(SPLIT), function(x) DM2[, SPLIT[[x]], drop = FALSE])
    ## iterate through each block combination, calculate correlation matrix between blocks and store them:
    resMAT <- foreach(M = DM2.list, .combine='cbind') %dopar% {
        vals <- cor(DM1, M)
    }
    return(as.matrix(resMAT))
}



# Fastcomputation of Pearson's correlation between two big matrices (OpenMP version)
# Inputs are the two matrices and the desired number of cores
FlashPPearsonCorOMP <- function (DM1, DM2, ncores=NULL) { 
    ncores <- if (is.null(ncores)) getDoParWorkers() else ncores
    resMAT <- PPearsonMatOMP(DM1, DM2, ncores)
    return(resMAT)
}



# Fastcomputation of Spearman's correlation between two big matrices 
# Inputs are the two matrices and the desired number of cores
# The function assumes that matrix with the largest dimension is
FlashPSpearmanCor <- function (DM1, DM2, ncores=NULL) { 
    if (!is.null(ncores))
        warning("Ignoring 'ncores' argument - using the value from getDoParWorkers().")
    ncores <- getDoParWorkers()
    nblocks <- 2*ncores 
    ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
    SPLIT <- split(1:ncol(DM2), ceiling(seq_along(1:ncol(DM2)) / ceiling(ncol(DM2) / nblocks)))
    DM2.list <- lapply(1:length(SPLIT), function(x) DM2[, SPLIT[[x]], drop = FALSE])
    ## iterate through each block combination, calculate correlation matrix between blocks and store them:
    resMAT <- foreach(M = DM2.list, .combine='cbind') %dopar% {
        #vals <- cor(DM1, M,method="spearman")
        vals <- PSpcor(DM1, M)
    }
    return(as.matrix(resMAT))
}




# Fastcomputation of Spearman's correlation between two big matrices (OpenMP version)
# Inputs are the two matrices and the desired number of cores
FlashPSpearmanCorOMP <- function (DM1, DM2, ncores=NULL) { 
    ncores <- if (is.null(ncores)) getDoParWorkers() else ncores
    R1 <- vapply(c(1:ncol(DM1)), function (x) rank(DM1[,x]), FUN.VALUE=double(length=nrow(DM1)))
    R2 <- vapply(c(1:ncol(DM2)), function (x) rank(DM2[,x]), FUN.VALUE=double(length=nrow(DM2)))
    resMAT <- PPearsonMatOMP(R1, R2, ncores)
    return(resMAT)
}




# Fastcomputation of Hellinger Distance between two big matrices 
# Inputs are the two matrices and the desired number of cores
# The function assumes that matrix with the largest dimension is
FlashPHellinger <- function (DM1, DM2, ncores=NULL) { 
    if (!is.null(ncores))
        warning("Ignoring 'ncores' argument - using the value from getDoParWorkers().")
    ncores <- getDoParWorkers()
    nblocks <- 2*ncores 
    ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
    SPLIT <- split(1:ncol(DM2), ceiling(seq_along(1:ncol(DM2)) / ceiling(ncol(DM2) / nblocks)))
    DM2.list <- lapply(1:length(SPLIT), function(x) DM2[, SPLIT[[x]], drop = FALSE])
    ## iterate through each block combination, calculate correlation matrix between blocks and store them:
    resMAT <- foreach (M = DM2.list, .combine='cbind') %dopar% {
        vals <- PHellingerMat(DM1, M)
    }
    return(as.matrix(resMAT))
}




# Fastcomputation of Hellinger Distance between two big matrices (OpenMP version)
# Inputs are the two matrices and the desired number of cores
FlashPHellingerOMP <- function (DM1, DM2, ncores=NULL) { 
    ncores <- if (is.null(ncores)) getDoParWorkers() else ncores
    resMAT <- PHellingerMatOMP(DM1, DM2, ncores)
    return(resMAT)
}




# Fastcomputation of Canberra distances between two big matrices 
# Inputs are the two matrices and the desired number of cores
# The function assumes that matrix with the largest dimension is
FlashPCanberra <- function (DM1, DM2, ncores=NULL) {
    if (!is.null(ncores))
        warning("Ignoring 'ncores' argument - using the value from getDoParWorkers().")
    ncores <- getDoParWorkers()
    nblocks <- 2*ncores 
    ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
    SPLIT <- split(1:ncol(DM2), ceiling(seq_along(1:ncol(DM2)) / ceiling(ncol(DM2) / nblocks)))
    DM2.list <- lapply(1:length(SPLIT), function(x) DM2[, SPLIT[[x]], drop = FALSE])
    ## iterate through each block combination, calculate correlation matrix between blocks and store them:
    resMAT <- foreach(M = DM2.list, .combine='cbind' ) %dopar% {
        vals <- PCanberraMat(DM1, M)
    }
    return(as.matrix(resMAT))
}




# Fastcomputation of Canberra distances between two big matrices (OpenMP version)
# Inputs are the two matrices and the desired number of cores
FlashPCanberraOMP <- function (DM1, DM2, ncores=NULL) {
    ncores <- if (is.null(ncores)) getDoParWorkers() else ncores
    resMAT <- PCanberraMatOMP(DM1, DM2, ncores)
    return(resMAT)
}




# Fast computation of Pearson's correlation for sparse matrices
# Input is a sparse matrix (see package "Matrix") and the desired number of cores
FlashSPearsonCor <- function (DM, ncores=NULL) { 
    if (!is.null(ncores))
        warning("Ignoring 'ncores' argument - using the value from getDoParWorkers().")
    ncores <- getDoParWorkers()
    nblocks <- 2 * ncores 
    SPLIT <- split(1:ncol(DM), ceiling(seq_along(1:ncol(DM)) / ceiling(ncol(DM) / nblocks)))
    ## iterate through each block combination, calculate correlation matrix between blocks and store them:
    corMAT <- foreach(i = 1:length(SPLIT), .export=c('sparse.cor','psparse.cor'),
                       .packages=('Matrix'), .combine='cbind' ) %dopar% {
        v <- SPLIT[[i]]
        return(psparse.cor(DM, DM[, v, drop = FALSE]))
    }
    return(as.matrix(corMAT))
}




