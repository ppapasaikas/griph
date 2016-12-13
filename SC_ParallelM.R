# Fast, memory efficient computation of Pearson's correlation for big matrices using BLAS implementations, parallelization and bigmemory objects 
# In a 16core machine with R' BLAS implementeation a 10Kx10K input matrix is processed in ~5' (~10x speed-up over the built-in impementation). 
# More optimized BLAS versions (ATLAS, openBLAS) should offer further improvements.
# Input is a matrix and the desired number of cores
FlashPearsonCor <- function (DM,ncores="all") { 
  require(foreach)
  require(doParallel)
  require (bigmemory)
  require (coop)      #(pcor) Fastest
  
  DM=as.big.matrix(DM)
  desc<-describe(DM)
  NCOL=desc@description["totalCols"][[1]]
  
  if(ncores=="all"){
    ncores = detectCores()
    ncores=min(48,ncores,ceiling(NCOL/200))
  } else{
    ncores=min(48,ncores,detectCores(),ceiling(NCOL/200))
  }
  
  nblocks = ncores 
  cl<-makeCluster(ncores)
  registerDoParallel(cl)
  corMAT <- matrix(data=0,nrow=NCOL,ncol=NCOL, dimnames=list(desc@description["colNames"][[1]],desc@description["colNames"][[1]]) )
  corMAT <- as.big.matrix(corMAT,type="double")
  bigCOR <- describe(corMAT)
  ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
  SPLIT=split(1:NCOL, ceiling(seq_along(1:NCOL)/ceiling(NCOL/nblocks)  ))
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate correlation matrix between blocks and store them:
  results<-foreach(i = 1:nrow(COMBS),.export=c('pcor'), .packages=('bigmemory') ) %dopar%{
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    DM<-attach.big.matrix(desc)
    corMAT<-attach.big.matrix(bigCOR)
    
    if(COMB[1]==COMB[2] || nrow(COMBS)< 20) {
      v=unique(c(G1,G2))
      COR <- pcor(as.matrix(DM[,v]))
      corMAT[v, v] <- COR
    }
    else{
      COR <- cor(as.matrix(DM[,G1]), as.matrix(DM[,G2]))
      corMAT[G1,G2]=COR
      corMAT[G2,G1]=t(COR)
    }
    
    COR <- NULL
    return(0L)
  }
  stopCluster(cl)
  return(as.matrix(corMAT))
}



# Fast, memory efficient computation of Shrinked Correlation Estimates for big matrices using BLAS implementations, parallelization and bigmemory objects 
# More optimized BLAS versions (ATLAS, openBLAS) should offer further improvements.
# Input is a matrix, the desired number of cores and the lambda parameter controlling shrinkage intensity
FlashShrinkCor <- function (DM,ncores="all",l=0.5) { 
  require(foreach)
  require(doParallel)
  require (bigmemory)
  #require (coop)      #(pcor) Fastest
  #require (HiClimR)  #(fastCor) Memory Efficient
  require (corpcor)  #(cor.shrink)  Incorporates weights and Shrinkage
  
  DM=as.big.matrix(DM)
  desc<-describe(DM)
  NCOL=desc@description["totalCols"][[1]]
  
  if(ncores=="all"){
    ncores = detectCores()
    ncores=min(48,ncores,ceiling(NCOL/150))
  } else{
    ncores=min(48,ncores,detectCores(),ceiling(NCOL/150))
  }
  
  nblocks = ncores 
  cl<-makeCluster(ncores)
  registerDoParallel(cl)
  corMAT <- matrix(data=0,nrow=NCOL,ncol=NCOL, dimnames=list(desc@description["colNames"][[1]],desc@description["colNames"][[1]]) )
  corMAT <- as.big.matrix(corMAT,type="double")
  bigCOR <- describe(corMAT)
  ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
  SPLIT=split(1:NCOL, ceiling(seq_along(1:NCOL)/ceiling(NCOL/nblocks)  ))
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  if (nrow(COMBS)>1) {COMBS <- COMBS[-which(COMBS[,1]==COMBS[,2]),]}
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate correlation matrix between blocks and store them:
  results<-foreach(i = 1:nrow(COMBS),.export=c('cor.shrink'), .packages=('bigmemory') ) %dopar%{
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    x<-attach.big.matrix(desc)
    corMAT<-attach.big.matrix(bigCOR)
    v=unique(c(G1,G2))
    COR <- cor.shrink(x[,v],verbose=FALSE,lambda=l)
    corMAT[v, v] <- COR
    COR <- NULL
    return(0L)
  }
  stopCluster(cl)
  return(as.matrix(corMAT))
}



# Fast, memory efficient computation of Spearman's correlation for big matrices using BLAS implementations, parallelization and bigmemory objects 
# In a 16core machine with R' BLAS implementeation a 10Kx10K input matrix is processed in ~5' (~10x speed-up over the built-in impementation). 
# More optimized BLAS versions (ATLAS, openBLAS) should offer further improvements.
# Input is a matrix and the desired number of cores
FlashSpearmanCor <- function (DM,ncores="all") { 
  require(foreach)
  require(doParallel)
  require (bigmemory)
  require (coop)      #(pcor) Fastest
  #require (HiClimR)  #(fastCor) Memory Efficient
  #require (corpcor)  #(cor.shrink)  Incorporates weights and Shrinkage
  
  DM=as.big.matrix(DM)
  desc<-describe(DM)
  NCOL=desc@description["totalCols"][[1]]
  NROW=desc@description["totalRows"][[1]]
  if(ncores=="all"){
    ncores = detectCores()
    ncores=min(48,ncores,ceiling(NCOL/200))
  } else{
    ncores=min(48,ncores,detectCores(),ceiling(NCOL/200))
  }
  
  nblocks = ncores 
  cl<-makeCluster(ncores)
  registerDoParallel(cl)
  corMAT <- matrix(data=0,nrow=NCOL,ncol=NCOL, dimnames=list(desc@description["colNames"][[1]],desc@description["colNames"][[1]]) )
  corMAT <- as.big.matrix(corMAT,type="double")
  bigCOR <-  describe(corMAT)
  RankMAT <- matrix(data=0,nrow=NROW,ncol=NCOL, dimnames=list(desc@description["rowNames"][[1]],desc@description["colNames"][[1]]) )
  RankMAT <- as.big.matrix(RankMAT,type="double")
  bigRank <- describe(RankMAT)
  
  ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
  SPLIT=split(1:NCOL, ceiling(seq_along(1:NCOL)/ceiling(NCOL/nblocks)  ))
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  COMBS <- unique(COMBS)
  
  ## iterate through each block and fill the rank matrix:
  results<-foreach(i = 1:length(SPLIT), .packages=('bigmemory') ) %dopar%{
    M<-attach.big.matrix(desc)
    RankMAT<-attach.big.matrix(bigRank)
    RANK=vapply( SPLIT[[i]],function (x) rank(M[,x]),FUN.VALUE =( double(length=nrow(M) ))  )
    RankMAT[,SPLIT[[i]]]=RANK
    RANK <- NULL
    return(0L)
  }
  
  ## iterate through each block combination, calculate correlation matrix between blocks and store them:
  results<-foreach(i = 1:nrow(COMBS),.export=c('pcor'), .packages=('bigmemory') ) %dopar%{
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    RankMAT<-attach.big.matrix(bigRank)
    corMAT<-attach.big.matrix(bigCOR)
    
    if(COMB[1]==COMB[2] || nrow(COMBS)< 20) {
      v=unique(c(G1,G2))
      COR <- pcor(as.matrix(RankMAT[,v]))
      corMAT[v, v] <- COR
    }
    else{
      COR <- cor(as.matrix(RankMAT[,G1]), as.matrix(RankMAT[,G2]))
      corMAT[G1,G2]=COR
      corMAT[G2,G1]=t(COR)
    }
    COR <- NULL
    return(0L)
  }
  stopCluster(cl)
  return(as.matrix(corMAT))
}



# Parallel regularized estimates of the inverse covarinace matrix using the QUIC implementation.
FlashGlasso <- function (DM,ncores="all",rho=0.3) { 
  require(foreach)
  require(doParallel)
  require (bigmemory)
  require (QUIC)  # Fast implementation of glasso for precision (inverse covariance) matrix estimation
  
  DM=as.big.matrix(DM)
  desc<-describe(DM)
  NCOL=desc@description["totalCols"][[1]]
  NROW=desc@description["totalRows"][[1]]
  
  if(ncores=="all"){
    ncores = detectCores()
    ncores=min(48,ncores,ceiling(NCOL/150))
  } else{
    ncores=min(48,ncores,detectCores(),ceiling(NCOL/150))
  }
  
  
  if (!is.matrix(rho) && length(rho) != 1 && length(rho) !=NCOL)  {  stop("Wrong number of elements in rho") }
  if (is.vector(rho)) { rho <- matrix(sqrt(rho)) %*% sqrt(rho) }
  if (length(rho) == 1) {rho <- matrix(rho, ncol = NCOL, nrow = NCOL)}
  
  
  nblocks = ncores 
  cl<-makeCluster(ncores)
  registerDoParallel(cl)
  PrecMAT <- matrix(data=0,nrow=NCOL,ncol=NCOL, dimnames=list(desc@description["colNames"][[1]],desc@description["colNames"][[1]]) )
  PrecMAT <- as.big.matrix(PrecMAT,type="double")
  bigPREC <- describe(PrecMAT)
  ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
  SPLIT=split(1:NCOL, ceiling(seq_along(1:NCOL)/ceiling(NCOL/nblocks)  ))
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  if (nrow(COMBS)>1) {COMBS <- COMBS[-which(COMBS[,1]==COMBS[,2]),]}
  COMBS <- unique(COMBS)
  
  ## iterate through each block combination, calculate sparse inv. cov. matrix between blocks and store them:
  results<-foreach(i = 1:nrow(COMBS),.export=c('QUIC'), .packages=('bigmemory') ) %dopar%{
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    x<-attach.big.matrix(desc)
    PrecMAT<-attach.big.matrix(bigPREC)
    v=unique(c(G1,G2))
    r=rho[v,v]
    PREC <- QUIC(x[v,v],rho=r,maxIter=25,tol=5e-02)
    PrecMAT[v, v] <- PREC$X
    PREC <- NULL
    return(0L)
  }
  stopCluster(cl)
  return(as.matrix(PrecMAT))
}



#Parallel regularized non-negative matrix factorization
FlashRNMF <- function (DM,ncores="all") { 
  require(foreach)
  require(doParallel)
  require (bigmemory)
  require (rNMF)  # Fast implementation of glasso for precision (inverse covariance) matrix estimation
  
  DM=as.big.matrix(DM)
  desc<-describe(DM)
  NCOL=desc@description["totalCols"][[1]]
  NROW=desc@description["totalRows"][[1]]
  
  if(ncores=="all"){
    ncores = detectCores()
    ncores=min(48,ncores,ceiling(NCOL/150))
  } else{
    ncores=min(48,ncores,detectCores(),ceiling(NCOL/150))
  }
  nblocks = ncores 
  cl<-makeCluster(ncores)
  registerDoParallel(cl)
  DMimp <- matrix(data=0,nrow=NROW,ncol=NCOL, dimnames=list(desc@description["rowNames"][[1]],desc@description["colNames"][[1]]) )
  DMimp <- as.big.matrix(DMimp,type="double")
  bigDMimp <- describe(DMimp)
  ## split column numbers into 'nblocks' groups
  SPLIT=split(1:NCOL, ceiling(seq_along(1:NCOL)/ceiling(NCOL/nblocks)  ))

  ## iterate through each block and apply rNMF:
  results<-foreach(i = 1:length(SPLIT),.export=c('rnmf'), .packages=('bigmemory') ) %dopar%{
    G <- SPLIT[[i]]
    x<-attach.big.matrix(desc)
    DMimp<-attach.big.matrix(bigDMimp)
    IMP <- rnmf(x[,G],k=6,alpha = 0.15,tol=1e-2,maxit=10,showprogress=FALSE,quiet=TRUE)$fit
    DMimp[, G] <- IMP
    IMP <- NULL
    return(0L)
  }
  stopCluster(cl)
  return(as.matrix(DMimp))
}



#Parallel computation of Gaussian Kernel Distance. For a NxD mstrix (N variables, D features) returns an NxN distance matrix
FlashGKRNL <- function (DM,ncores="all",sigma=1) { 
  require(foreach)
  require(doParallel)
  require (bigmemory)
  require (KRLS)  # Gaussian Kernel Distance
  
  DM=as.big.matrix(DM)
  desc<-describe(DM)
  NROW=desc@description["totalRows"][[1]]
  
  if(ncores=="all"){
    ncores = detectCores()
    ncores=min(48,ncores,ceiling(NROW/150))
  } else{
    ncores=min(48,ncores,detectCores(),ceiling(NROW/150))
  }
  
  nblocks = ncores 
  cl<-makeCluster(ncores)
  registerDoParallel(cl)
  DMAT <- matrix(data=0,nrow=NROW,ncol=NROW, dimnames=list(desc@description["rowNames"][[1]],desc@description["rowNames"][[1]]) )
  DMAT <- as.big.matrix(DMAT,type="double")
  bigD <- describe(DMAT)
  ## split column numbers into 'nblocks' groups and create all unique combinations of blocks
  SPLIT=split(1:NROW, ceiling(seq_along(1:NROW)/ceiling(NROW/nblocks)  ))
  COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
  COMBS <- t(apply(COMBS, 1, sort))
  if (nrow(COMBS)>1) {COMBS <- COMBS[-which(COMBS[,1]==COMBS[,2]),]}
  COMBS <- unique(COMBS)
  ## iterate through each block combination, calculate kernel distance between blocks and store them:
  results<-foreach(i = 1:nrow(COMBS),.export=c('gausskernel'), .packages=('bigmemory') ) %dopar%{
    COMB <- COMBS[i, ]
    G1 <- SPLIT[[COMB[1]]]
    G2 <- SPLIT[[COMB[2]]]
    x<-attach.big.matrix(desc)
    DMAT<-attach.big.matrix(bigD)
    v=unique(c(G1,G2))
    DIST <- gausskernel(x[v,],sigma=sigma)
    DMAT[v, v] <- DIST
    DIST <- NULL
    return(0L)
  }
  stopCluster(cl)
  return(as.matrix(DMAT))
}




