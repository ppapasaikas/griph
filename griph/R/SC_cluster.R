#' Combined distance calculations.
#' 
#' @description Calculate and combine several between-cell distance measures to
#'     obtain robust cell-to-cell distances.
#' 
#' @param M gene-by-cell count matrix.
#' @param PPearsonCor Function to calculate Pearson correlation between the columns of two matrices.
#' @param PSpearmanCor Function to calculate Spearman correlation between the columns of two matrices.
#' @param PHellinger Function to calculate Hellinger Distance between the columns of two matrices.
#' @param PCanberra Function to calculate Canberra Distance between the columns of two matrices.
#' @param ShrinkCor Function to calculate shrinkage correlation.
#' 
#' @return cell-by-cell distance matrix.
WScor <- function (M, PPearsonCor, PSpearmanCor, PHellinger, PCanberra, ShrinkCor=ShrinkCor   ) {
        
    
    nBulks=min(1500, ceiling(5*ncol(M)) )
    FBsize=3
    Gcounts=colSums(M)
    HighQual=c(1:ncol(M))[-which(Gcounts <  quantile(Gcounts,0.05)  )]
    rep.ind=rep(HighQual, ceiling(nBulks*FBsize/length(HighQual))         )
    SMPL=sample(rep.ind,(nBulks*FBsize),replace=FALSE    )
    FB=matrix(0,nrow(M),nBulks)
    if (FBsize>1){
        index=1
        for (c in 1:nBulks) {
        FB[,c]=rowSums(M[,SMPL[ index:(index+FBsize-1)]   ])
        index=index+FBsize
        }
    }
    else {FB=M[,SMPL[1:min(nBulks,ncol(M))]]}


    D=PPearsonCor(log2(FB+1),log2(M+1))
    R=vapply(c(1:ncol(D)),function (x) rank(D[,x]),FUN.VALUE=double(length=nrow(D) ) )  #pearson's cor    

    ######## Counts per Million:
    CellCounts=colSums(FB)
    FB=sweep(FB,2,CellCounts,FUN="/")
    FB=FB*1e6
    CellCounts=colSums(M)
    M=sweep(M,2,CellCounts,FUN="/")
    M=M*1e6 
    
    Dt=PCanberra( log2(FB+1),log2(M+1) )   #
    Dt=1-( (Dt-min(Dt))/ diff(range(Dt)) )
    D=D+Dt
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #canberra

    Dt=PSpearmanCor(FB,M)
    D=D+Dt
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #spearman's cor 
    
    Dt=PHellinger(FB,M)
    Dt=1-( (Dt-min(Dt))/ diff(range(Dt)) )
    D=D+Dt
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #Hellinger distance
    
    R=R/4
    WM=D/4
    
    ##### Kernelize average distance
    ave=mean(WM[which(WM>0)])
    GK=exp(- ( ((1-WM)^2) / ((1-ave)^2) ) ) 
    
    ##### Compute weights of Fake bulks based on local density
    W=apply(GK ,1,sum )^(3) #
    minW=quantile(W,probs=seq(0,1,0.1))[[2]]
    W=((minW)/(W+1*minW))  #
    W=W/max(W)
    R=R^2

    R=ShrinkCor( R ,verbose=FALSE,lambda=0, w=W ) 
    return(as(R,"matrix"))
}





Spcor <- function (M) { #Spearman's correlation using coop
    R=vapply(c(1:ncol(M)),function (x) rank(M[,x]),FUN.VALUE=double(length=nrow(M) ) )
    R=coop::pcor(R)
    return(R)
}


PSpcor <- function (M1,M2) { #Spearman's correlation using coop
    R=stats::cor(M1,M2,method="spearman")
    return(R)
}


canberra <- function (M) {
    D=as.matrix(dist(t(M),method="canberra" )) #
    return(D)
}


Qglasso <- function (C,rho=0.5,tol=1e-3,maxIter=100,msg=0) {
    X<-QUIC::QUIC(C,rho=rho,tol=tol,maxIter=maxIter,msg=msg)$X 
    return(X)
}


PPR <- function (G,df=0.75){
    if (!isTRUE(igraph::is.igraph(G))) {  
        if (!isSymmetric(G) )  {stop ("!! G should be either a graph object or a symmetric matrix !!")}   
        G=igraph::graph.adjacency(G [1:nrow(G),1:nrow(G)],mode=c("max"),weighted=TRUE,diag=FALSE)
    }
    L=length(V(G))
    PR=diag(nrow=L)
    vals=sapply(1:L, function(x) igraph::page_rank(G, vids=c(1:L),personalized=PR[,x],damping=df )$vector   )
    PR[1:length(PR)]=log2(unlist(vals)+(0.01/L))
    PR[lower.tri(PR)]=(PR[lower.tri(PR)]+t(PR)[lower.tri(t(PR))])/2
    PR[upper.tri(PR)]=t(PR)[upper.tri(t(PR))]
    PR=2^PR
    return(PR)
}







#' @title Unsupervised graph-based clustering of single cell RNA-seq data.
#' 
#' @description \code{SC_cluster} takes a gene-by-cell matrix with read counts as
#'     input and groups individual cells into clusters of similar cell (cell types).
#'     The number of cell types is automatically determined.
#' 
#' @details 
#' TODO
#' 
#' @param DM Count data (n genes-by-k cells) or directly a correlation k-by-k matrix.
#'     Required argument.
#' @param use.par If \code{TRUE}, use parallel versions of distance calculation
#'     functions based on \code{\link[foreach]{foreach}} (see details).
#' @param ncores a numeric(1) or character(1), either specifying the number of
#'     parallel jobs to use, or \code{"all"} (default) to use up to 90 percent of
#'     the available cores. Ignored if \code{use.par==FALSE}.
#' @param is.cor If \code{TRUE}, \code{DM} is assumed to be a correlation matrix.
#'     Otherwise (the default), \code{DM} is assumed to be a genes-by-cells count
#'     matrix, and a correlation matrix will be computed.
#' @param filter T/F Filter genes according to cv=f( mean ) fitting. Default: \code{TRUE}.
#' @param do.glasso T/F  Sparsify cell correlation matrix using graphical lasso. Default: \code{TRUE}.
#' @param rho Inverse covariance matrix regularization (graph sparsification) parameter -> [0,1].
#'     Default=0.25.  The parameter is then automatically scaled to account for
#'     number of variables or converted into a matrix and adjusted according to
#'     the batch.penalty factor to account for BatchAssignment (if given).
#' @param pr.iter A numberic scalar defining the number of iterations to use for
#'     re-weighting the graph edges using personalized PageRank node similarity score.
#'     Set to zero to deactivate.
#' @param batch.penalty [0,1] rho scaling factor for enforcing topological constraints
#'     variables according to \code{BatchAssignment}. For penalty p  -> rho_same_batch=rho^(1-p),
#'     rho_diff_batch=rho^(1+p). It is ignored if \code{BatchAssignment==NULL}. 
#' @param comm.method  Community detection algorithm. See igraph "communities" 
#' @param ncom Forces the community detection algorithm to a fixed number of communities.
#'     If \code{NULL} (default) then the optimal determined number of clusters of
#'     the used  community detection algorithm.
#' @param ClassAssignment If available a numeric vector of length \code{k} with numeric
#'     class labels (e.g-> c(1,2,1,1,1,2,3,3,3,1,2))
#' @param BatchAssignment If available a numeric vector of length \code{k} with numeric
#'     batch labels (e.g-> c(1,2,1,1,1,2,3,3,3,1,2))
#' @param plotG if \code{TRUE} plots the resulting graph
#' @param maxG Approximate maximal number of vertices to include when plotting the graph.
#' @param fsuffix A suffix added to the file names of output plots. If not given
#'     it will use a random 5 character string.
#' @param image.format Specifies the format of the created images. Currently only pdf and png filetypes are supported.
#' @param ... Additional arguments passed from griph_cluster
#' 
#' @return Currently a list with the clustering results.

SC_cluster <- function(DM, use.par=FALSE,ncores="all",is.cor = FALSE,
                       filter = FALSE, do.glasso=TRUE, rho = 0.25, pr.iter = 1, batch.penalty = 0.5,
                       comm.method=igraph::cluster_infomap,ncom=NULL,
                       ClassAssignment = rep(1,ncol(DM)), BatchAssignment = NULL,
                       plotG = TRUE, maxG = 2500, fsuffix = RandString(), image.format='png', ...){
    
    #######Internal parameters for testing puproses only:  
    qnt=8 #Max gene expression decile for imputation (e.g 8->bottom 70% of the genes are imputed) 
    rho=rho+( (ncol(DM)/1e9)^0.2) #Scale rho for number of cells. MAKE SURE rho is <=1
    rho=min(0.9,rho)
    #######Switch to parallelized functions if use.par=TRUE
    
    PPearsonCor=stats::cor
    PSpearmanCor=PSpcor
    PHellinger=PHellingerMat
    PCanberra=PCanberraMat
    
    ShrinkCor=corpcor::cor.shrink
    Glasso=Qglasso
    PPRank=PPR
    
    if (isTRUE(use.par)) {  
        PPearsonCor=FlashPPearsonCor
        PSpearmanCor=FlashPSpearmanCor
        PHellinger=FlashPHellinger
        PCanberra=FlashPCanberra 
        
        ShrinkCor=FlashShrinkCor
        Glasso=FlashGlasso
        PPRank=FlashPPR 
    }
    
    
    ##### Strip dimnames:
    CellIds=colnames(DM)
    dimnames(DM)=list(NULL,NULL)
    
    if (!isTRUE(is.cor)) {  
        
        message("Preprocessing...", appendLF = FALSE)
        
        AllZeroRows=which  ( rowSums(DM)<1e-9 )
        if(length(AllZeroRows)>0){
            DM=DM[-AllZeroRows , ] 
        }
        
        
        ##########  Remove invariant genes:
        meanDM=mean(DM)
        nSD=apply(DM,1,function(x) sd(x)/meanDM)
        ConstRows=which   ( nSD < 0.25 )
        if(length(ConstRows)>0){
            DM=DM[-ConstRows , ]
        }
        
        #############CV=f(mean) -based filtering:
        CellCounts=colSums(DM)
        nDM=sweep(DM,2,CellCounts,FUN="/")
        nDM=DM*10000 #Counts per 10K
        if(is.null(filter)){
            medianComplexity=median(apply(DM,2,function(x) sum(x>0))) 
            filter=ifelse( medianComplexity > 2500,TRUE,FALSE)
            message("Median Library Complexity: ",medianComplexity," --> Gene Filtering: ", filter ,"\r")
            
        }
        if (filter){
            X1=log2(rowMeans(nDM+1/ncol(nDM)))
            Y1=apply(nDM,1,function(x) log2(sd(x)/mean(x+1/length(x))+1/length(x) )  )
            m=nls(Y1 ~ a*X1+b, start=list(a=-5,b=-10)  )
            Yhat=predict(m)
            DM=DM[which(Y1 > Yhat),]
        }
        
        NoData <- which(colSums(DM) == 0)
        if (length(NoData > 0)){
            DM <- DM[,-NoData]
            ClassAssignment <- ClassAssignment[-NoData]
            BatchAssignment <- BatchAssignment[-NoData]
            CellCounts <- CellCounts[-NoData]
        }
        message("done")
        
        C=list()

        message("Calculating Pairwise and Diffused Similarities...", appendLF = FALSE)
        
        C[[2]]=WScor(DM, PPearsonCor=PPearsonCor, PSpearmanCor=PSpearmanCor, PHellinger=PHellinger, PCanberra=PCanberra, ShrinkCor=ShrinkCor )
        nDM<-NULL
        
        message("done")
    }
    
    
    else {
        C=list()
        C[[2]]=DM
        
    }
    
    Cuse=C[[2]]

    if (do.glasso==TRUE){
    
        ############### glasso-based graph structure estimation: #####################
        message("Estimating Graph Structure...", appendLF = FALSE)
    
        RHO=matrix(rho,nrow=nrow(Cuse),ncol=ncol(Cuse) )

        ###### Mutual k-nn based pruning as per Harel and Koren 2001 SIGKDD:
        k=min(max(5*sqrt(ncol(Cuse)),50) ,floor(ncol( Cuse))/1.5)  
        kN=get.knn(Cuse,k=k )
        for(i in 1:ncol(  RHO  )){
            RHO[ -kN[,i],i]=min(1,1.5*rho)
            RHO[i,-kN[,i] ]=min(1,1.5*rho)
        }
        kN=NULL
    
        if(!is.null(BatchAssignment)){
            rL=min(1,rho^(1-batch.penalty))
            rS=rho^(1+batch.penalty)
            RHO=mapply (  function(r,c) {if (BatchAssignment[r]==BatchAssignment[c]) {rL} else { rS } },row(Cuse),col(Cuse) )
            RHO=matrix(RHO,nrow=nrow(Cuse),ncol=ncol(Cuse) )
            ###### Mutual k-nn based pruning as per Harel and Koren 2001 SIGKDD:
            k=min(max( floor(5*sqrt(ncol(Cuse))),50) ,floor(ncol( Cuse))/1.5)  
            kN=get.knn(Cuse,k=k )
            for(i in 1:ncol(  RHO  )){
                RHO[ -kN[,i],i]=pmin(rep(1,length(RHO[ -kN[,i],i])  ),1.5*RHO[ -kN[,i],i])
                RHO[i,-kN[,i] ]=pmin(rep(1,length(RHO[ -kN[,i],i])  ),1.5*RHO[i,-kN[,i] ])
            }
        }

        
        C=NULL
        tol=5e-02
        maxIter=40
        if (ncol(Cuse)<200) {tol=1e-02;maxIter=80}
        if (ncol(Cuse)>800) {tol=1e-01;maxIter=20}
        X<-Glasso(Cuse,rho=RHO,tol=tol,maxIter=maxIter,msg=0)  #0.1 for C4/ 0.3 for C3 / 0.5 for C2
        ADJ= -X
        ADJ=Matrix::Matrix(data=ADJ,sparse=TRUE,dimnames=NULL)
        RHO<-NULL
        message("done")
    }
    
    else{
    #####Assumes Cuse comes sparse...    
    ADJ=Cuse    
    }
    
    ######## Graph weights:
    message("Calculating edge weights and knn-based pruning...", appendLF = FALSE)
    ave=mean(Cuse[which(ADJ>0)])

    ADJ[which(ADJ>0)]=exp(- ( ((1-Cuse[which(ADJ>0)])^2) / ((1-ave)^2) ) )   #Kernelize distance according to Haren and Koren 2001 section 3
    message("\n0\n",is(ADJ,'sparseMatrix'),"\n")
    
    ADJ[which(ADJ < 0)]=0
    diag(ADJ)=1
    message("\n1\n",is(ADJ,'sparseMatrix'),"\n")
    
    ###### Mutual k-nn based pruning as per Harel and Koren 2001 SIGKDD:
    k=min(max( floor(5*sqrt(ncol(ADJ))) ,100) ,floor(ncol( ADJ))/1.5)  
    kN=get.knn(ADJ,k=k )
    for(i in 1:ncol(  ADJ  )){
        ADJ[ -kN[,i],i]=0
        ADJ[i,-kN[,i] ]=0
    }
    apply(ADJ,2,function(x) x[] )
    kN=NULL
    message("\n2\n",is(ADJ,'sparseMatrix'),"\n")
    
    
    
    GRAO<-igraph::graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
    message("done")
    
    message("\n3\n",is(ADJ,'sparseMatrix'),"\n")
    
    
    if(pr.iter > 0) {
        for (i in 1:pr.iter){
            message("Pruning based on global node similarity: ",i," / ",pr.iter, "\r", appendLF = FALSE)
            flush.console()
            df=0.75
            PR=PPRank(GRAO,df=df)
            ADJ=PR
            ADJ[which(PR < (0.01/(ncol(ADJ))  ) ) ]=0
            ###### Mutual k-nn based pruning as per Harel and Koren 2001 SIGKDD:
            k=min(max(5*sqrt(ncol(Cuse)),100) ,floor(ncol( Cuse))/1.5)   
            kN=get.knn(PR,k=k )
            for(i in 1:ncol(ADJ)){
                ADJ[ -kN[,i],i]=0
                ADJ[i,-kN[,i] ]=0
            }
            ave=mean(Cuse[which(ADJ>0)])
            ADJ[which(ADJ>0)]=exp(- ( ((1-Cuse[which(ADJ>0)])^2) / ((1-ave)^2) ) )   #According to Harel and Koren 2001 SIGKDD section 3
            GRAO<-igraph::graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
        }
    }


    if (do.glasso==TRUE){
    ADJ[which(X>=0)]=0  #Reinstate 0 state from glasso
    X<-NULL
    }
    
    message("\n4\n",is(ADJ,'sparseMatrix'),"\n")
    
    pct <- 1
    if (median(igraph::degree(GRAO)) > 10) {
        pct <- min(1,1/(median(igraph::degree(GRAO))^0.1 )  )
        message("\tkeeping ", round(100*pct,1), "% of edges")
        ADJtemp <- apply(ADJ,2,function(x) sparsify(x,pct) )
        GRAO <- igraph::graph.adjacency(ADJtemp, mode=c("max"), weighted=TRUE, diag=FALSE)
        ADJtemp <- NULL
    }
    
    message("\n5\n",is(ADJ,'sparseMatrix'),"\n")
    
    
    GRAO<-igraph::set.vertex.attribute(GRAO, "class", value=as.character(ClassAssignment))

    ######## COMMUNITY DETECTION #########
    message("Detecting Graph Communities...")
    memb <- comm.method(GRAO)
    
    #if called with ncom:
    if (!is.null(ncom)) {
        min.csize <-max(4, ceiling(0.2 * sqrt(length(memb$membership))))
        if (!is.null(memb$memberships) &
            max(apply(memb$memberships, 1, function(x) sum(table(x) > min.csize))) >= ncom) {    #check if multi level louvain was used and a level has enough substantial clusters
            nclust <- apply(memb$memberships, 1, function(x) sum(table(x) > min.csize))
            memb$membership <- memb$memberships[which.min(abs(nclust - ncom)),]

        } else {     # otherwise run simple kmeans on a 6D projection generated by largeVis
            message(memb$algorithm, " returned too few clusters, switching to k-NN...")
            x <- projectKNNs(as_adj(GRAO, names=FALSE, sparse=TRUE), dim=6,
                             sgd_batches=min(20000 / (sum(ADJ!=0) / 2), 0.99),
                             M=3, gamma=10, alpha=0.1, useDegree=TRUE)

            ncenters <- ncom
            for (i in 0:4) { # try up to five times to get ncom clusters with at least min.csize members
                memb$membership <- kmeans(t(x), centers=ncenters)$cluster
                nclust <- sum(table(memb$membership) > min.csize)

                if (nclust >= ncom) {
                    break
                } else {
                    ncenters <- ncenters + (ncom - nclust)
                }
            }
        }
    }
    csize <- table(memb$membership)
    ConfMatrix <- table(predicted=memb$membership, true=ClassAssignment)
    message("done")
    
    V(GRAO)$membership <- memb$membership
    V(GRAO)$community.size <- csize[memb$membership]
    E(GRAO)$weight <- edge.weights(memb, GRAO, weight.within=2, weight.between=0.5)
    
    ######### Optimal Mapping between true and estimated class assignments: ##########
    mapping <- mapLabelsGreedy(memb$membership, ClassAssignment)
    misclErr <- classError(memb$membership, ClassAssignment, mapping)
    
    #### Add back Cell Ids to igraph object, ADJ, MEMB and prepare return value
    
    dimnames(ADJ) <- list(CellIds,CellIds)
    dimnames(Cuse) <- list(CellIds,CellIds)
    names(memb$membership) <- CellIds
    
    V(GRAO)$labels=CellIds

    message("\n6\n",is(ADJ,'sparseMatrix'),"\n")
    
    
    ret <- list(MEMB=memb$membership, MEMB.true=ClassAssignment,
                DISTM=ADJ, CORM=Cuse, ConfMatrix=ConfMatrix,
                miscl=misclErr, GRAO=GRAO, plotGRAO=NULL)

    ######### graph visualization
    if (plotG)
        ret[["plotGRAO"]] <- plotGraph(ret, maxG = maxG, fsuffix = fsuffix,
                                       image.format = image.format, quiet = FALSE)
    
    
    
    return(ret)
}
