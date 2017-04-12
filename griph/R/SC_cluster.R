#' Combined distance calculations.
#' 
#' @description Calculate and combine several between-cell distance measures to
#'     obtain robust cell-to-cell distances.
#' 
#' @param M gene-by-cell count matrix.
#' @param C1 cell-by-cell (correlation?) matrix.
#' @param CanberraDist Function to calculate Canberra distance.
#' @param SpearmanCor Function to calculate Spearman rank correlation.
#' @param HellingerDist Function to calculate Hellinger distance.
#' @param ShrinkCor Function to calculate shrinkage correlation.
#' 
#' @return cell-by-cell distance matrix.
WScor <- function (M, C1=matrix(1,ncol(M),ncol(M) ),
                   CanberraDist=CanberraDist, SpearmanCor=SpearmanCor,
                   HellingerDist=HellingerDist, ShrinkCor=ShrinkCor  ) { #Spearman's weighted correlation using cor.shrink
    D=C1
    R=vapply(c(1:ncol(C1)),function (x) rank(C1[,x]),FUN.VALUE=double(length=nrow(C1) ) )  #pearson's cor
    
    Dt=CanberraDist(log2(M+1)) #
    #Dt=as.matrix(dist(t(log2(M+1)),method="canberra" )) #
    Dt=1-( (Dt-min(Dt))/ diff(range(Dt)) )
    D=D+Dt
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #canberra
    
    Dt=SpearmanCor((M))
    D=D+Dt
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #spearman's cor 
    
    Dt=HellingerDist(M)
    Dt=1-( (Dt-min(Dt))/ diff(range(Dt)) )
    D=D+Dt
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #Hellinger distance
    
    R=R/4
    WM=D/4
    
    ##### Kernelize average distance
    ave=mean(WM[which(WM>0)])
    GK=exp(- ( ((1-WM)^2) / ((1-ave)^2) ) ) 
    
    ##### Compute weights based on local density
    W=apply(GK ,1,sum )^(3) #
    minW=quantile(W,probs=seq(0,1,0.1))[[2]]
    W=((minW)/(W+1*minW))  #
    W=W/max(W)
    R=ShrinkCor( R^2 ,verbose=FALSE,lambda=0,w=W ) 
    return(as(R,"matrix"))
}


Spcor <- function (M) { #Spearman's correlation using coop
    R=vapply(c(1:ncol(M)),function (x) rank(M[,x]),FUN.VALUE=double(length=nrow(M) ) )
    R=coop::pcor(R)
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
#' @param impute Perform imputation prior to clustering. Default: \code{FALSE}.
#' @param filter T/F Filter genes according to cv=f( mean ) fitting. Default: \code{TRUE}.
#' @param diffuse.iter Number of similarity difusion iterations -> [2,5]/ Default=2.
#' @param rho Inverse covariance matrix regularization (graph sparsification) parameter -> [0,1].
#'     Default=0.25.  The parameter is then automatically scaled to account for
#'     number of variables or converted into a matrix and adjusted according to
#'     the batch.penalty factor to account for BatchAssignment (if given).
#' @param batch.penalty [0,1] rho scaling factor for enforcing topological constraints
#'     variables according to \code{BatchAssignment}. For penalty p  -> rho_same_batch=rho^(1-p),
#'     rho_diff_batch=rho^(1+p). It is ignored if \code{BatchAssignment==NULL}. 
#' @param ClassAssignment If available a numeric vector of length \code{k} with numeric
#'     class labels (e.g-> c(1,2,1,1,1,2,3,3,3,1,2))
#' @param BatchAssignment If available a numeric vector of length \code{k} with numeric
#'     batch labels (e.g-> c(1,2,1,1,1,2,3,3,3,1,2))
#' @param plotG if \code{TRUE} plots the resulting graph
#' @param maxG Approximate maximal number of vertices to include when plotting the graph.
#' @param plotSP if \code{TRUE} plots the speactral projection of the graph in 2D
#' @param fsuffix A suffix added to the file names of output plots. If not given
#'     it will use a random 5 character string.
#' @param image.format Specifies the format of the created images. Currently only pdf and png filetypes are supported.

#' 
#' @return Currently a list with the clustering results.

SC_cluster <- function(DM, use.par=FALSE,ncores="all",is.cor = FALSE, impute = FALSE, filter = FALSE, rho = 0.25,
                       diffuse.iter = 2, batch.penalty = 0.5,
                       ClassAssignment = rep(1,ncol(DM)), BatchAssignment = NULL,
                       plotG = TRUE, maxG = 2500, plotSP = FALSE, fsuffix = RandString(), image.format='png' ){
    
    #######Internal parameters for testing puproses only:  
    comm.method=igraph::cluster_infomap  # Community detection algorithm. See igraph "communities" 
    ncom=NULL  #Forces the community detection algorithm to a fixed number of communities. Only possible with hierarchical methods (like fast_greedy),
    qnt=8 #Max gene expression decile for imputation (e.g 8->bottom 70% of the genes are imputed) 
    rho=rho+( (ncol(DM)/1e9)^0.2) #Scale rho for number of cells. MAKE SURE rho is <=1
    rho=min(0.9,rho)
    

    #######Switch to parallelized functions if use.par=TRUE
    PearsonCor=coop::pcor
    SpearmanCor=Spcor
    CanberraDist=canberra
    ShrinkCor=corpcor::cor.shrink
    HellingerDist=HellingerMat
    Glasso=Qglasso
    RNMF=rNMF::rnmf
    PPRank=PPR
    
    if (isTRUE(use.par)) {  
        PearsonCor=FlashPearsonCor
        SpearmanCor=FlashSpearmanCor
        CanberraDist=FlashCanberra
        ShrinkCor=FlashShrinkCor
        HellingerDist=FlashHellinger
        Glasso=FlashGlasso
        RNMF=FlashRNMF
        PPRank=FlashPPR 
        #####Register cluster here, then stop cluster at the end of processing
        if(ncores=="all"){
            ncores = parallel::detectCores()
            ncores=min(48,floor(0.9*ncores),ceiling(ncol(DM)/200))
        } else{
            ncores=min(48,ncores,floor(0.9*parallel::detectCores()),ceiling(ncol(DM)/200))
        }
        cl<-parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
    }
    
    
    
    if(!is.null(BatchAssignment)){
        if (length(BatchAssignment) != ncol(DM) ) {stop ("!Length of BatchAssignment vector should have length == ncol(DM)")}
    }
    
    ptm=proc.time()
    
    
    if (!isTRUE(is.cor)) {  
        
        if (impute==TRUE){
            message("Imputing...", appendLF = FALSE)
            DMimp=RNMF(DM,k =6,alpha = 0.15,tol=1e-2,maxit=10,showprogress=FALSE,quiet=TRUE)$fit
            GF=rowSums(DM)/sum(rowSums(DM))
            QNT=quantile(GF,probs=seq(0,1,0.1))
            Low=which(GF <= QNT[qnt]  )  
            DM[Low,]=DMimp[Low,]
            DMimp=NULL
            message("done")
        }
        
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
        
        NoData=which(colSums(DM)==0)
        if (length(NoData>0)){
            DM=DM[,-NoData]
            ClassAssignment=ClassAssignment[,-NoData]
            BatchAssignment=BatchAssignment[,-NoData]
            CellCounts=CellCounts[-NoData]
        }
        #
        
        ##### Strip dimnames:
        CellIds=colnames(DM)
        dimnames(DM)=NULL
        
        C=list()
        C[[1]]=PearsonCor(log2(DM+1))
        
        message("done")
    }
    
    
    else {
        C=list()
        C[[1]]=DM
    }
    
    averageDist=NULL
    if(is.null(diffuse.iter)) {
        diffuse.iter=2  
        averageDist=TRUE
    }
    
    if (diffuse.iter > 1) {
        for(order in 2:diffuse.iter){
            message("Calculating Pairwise and Diffused Similarities: ", order-1, " / ", diffuse.iter-1)
            C[[order]]=WScor( nDM,C1=C[[1]], CanberraDist=CanberraDist,
                              SpearmanCor=SpearmanCor, HellingerDist=HellingerDist, ShrinkCor=ShrinkCor  )
            
            if (order >2) {C[[order-1]]=NA}
        }
    }
    nDM<-NULL
    
    ####Normalize Initial Correlation matrix:
    W=pmax(1e-1,colMeans(C[[1]]))/mean(colMeans(C[[1]]))
    W=sqrt(W) %o% sqrt(W)
    C[[1]]=( C[[1]] / W )
    
    #####Metric by averaging:
    if(!is.null(averageDist)){
        Cuse=(C[[1]]+C[[2]])/2
    }
    
    else{  
        Cuse=as.matrix(C[[diffuse.iter]])
    }
    
    ClassAssignment.numeric <- as.numeric(factor(ClassAssignment, levels=unique(ClassAssignment)))

    ############### glasso-based graph structure estimation: #####################
    message("Estimating Graph Structure...", appendLF = FALSE)

    RHO=matrix(rho,nrow=nrow(Cuse),ncol=ncol(Cuse) )
    
    if(!is.null(BatchAssignment)){
        rL=min(1,rho^(1-batch.penalty))
        rS=rho^(1+batch.penalty)
        RHO=mapply (  function(r,c) {if (BatchAssignment[r]==BatchAssignment[c]) {rL} else { rS } },row(Cuse),col(Cuse) )
        RHO=matrix(RHO,nrow=nrow(Cuse),ncol=ncol(Cuse) )
    }
    
    
    ###### Mutual k-nn based pruning as per Harel and Koren 2001 SIGKDD:
    kvect=rep(  min(max(5*sqrt(ncol(Cuse)),50) ,floor(ncol( Cuse))/1.5)   ,    ncol(Cuse)  )
    kN=get.knn(Cuse,k=kvect )
    for(i in 1:ncol(  RHO  )){
        RHO[ -kN[,i],i]=min(1,1.5*rho)
        RHO[i,-kN[,i] ]=min(1,1.5*rho)
    }
    
    
    C=NULL
    tol=5e-02
    maxIter=40
    if (ncol(Cuse)<200) {tol=1e-02;maxIter=80}
    if (ncol(Cuse)>800) {tol=1e-01;maxIter=20}
    X<-Glasso(Cuse,rho=RHO,tol=tol,maxIter=maxIter,msg=0)  #0.1 for C4/ 0.3 for C3 / 0.5 for C2
    ADJ= -X
    X<-NULL
    
    message("done")
    
    ######## Graph weights:
    message("Calculating edge weights and knn-based pruning...", appendLF = FALSE)
    ave=mean(Cuse[which(ADJ>0)])
    ADJ[which(ADJ>0)]=exp(- ( ((1-Cuse[which(ADJ>0)])^2) / ((1-ave)^2) ) )   #Kernelize distance according to Haren and Koren 2001 section 3
    ADJ[which(ADJ < 0)]=0
    diag(ADJ)=1
    
    
    ###### Mutual k-nn based pruning as per Harel and Koren 2001 SIGKDD:
    kvect=rep(  min(max(5*sqrt(ncol(ADJ)),100) ,floor(ncol( ADJ))/1.5)   ,    ncol(ADJ)  )
    kN=get.knn(ADJ,k=kvect )
    for(i in 1:ncol(  ADJ  )){
        ADJ[ -kN[,i],i]=0
        ADJ[i,-kN[,i] ]=0
    }

    GRAO<-igraph::graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
    message("done")
    
    niter=1
    for (i in 1:niter){
        message("Pruning based on global node similarity: ",i," / ",niter, "\r", appendLF = FALSE)
        flush.console()
        df=0.75
        PR=PPRank(GRAO,df=df)
        ADJ=PR
        ADJ[which(PR < (0.01/(ncol(ADJ))  ) ) ]=0
        ###### Mutual k-nn based pruning as per Harel and Koren 2001 SIGKDD:
        kvect=rep( min(max(5*sqrt(ncol(Cuse)),100) ,floor(ncol( Cuse))/1.5)    ,ncol(ADJ)  )
        kN=get.knn(PR,k=kvect )
        for(i in 1:ncol(ADJ)){
            ADJ[ -kN[,i],i]=0
            ADJ[i,-kN[,i] ]=0
        }
        
        #PR=PR/ ( max(PR[upper.tri(PR)])+0.01/ncol(ADJ)  )
        #diag(PR)=1
        ave=mean(Cuse[which(ADJ>0)])
        ADJ[which(ADJ>0)]=exp(- ( ((1-Cuse[which(ADJ>0)])^2) / ((1-ave)^2) ) )   #According to Harel and Koren 2001 SIGKDD section 3
        #ADJ=(ADJ*(PR))
        GRAO<-igraph::graph.adjacency(ADJ,mode=c("max"),weighted=TRUE,diag=FALSE)
    }
    message("")
    
    
    pct <- 1
    if (median(igraph::degree(GRAO)) > 8) {
        pct <- min(1,1/(median(igraph::degree(GRAO))^0.25 )  )
        message("\tkeeping ", round(100*pct,1), "% of edges")
        ADJtemp <- apply(ADJ,1,function(x) sparsify(x,pct) )
        GRAO <- igraph::graph.adjacency(ADJtemp, mode=c("max"), weighted=TRUE, diag=FALSE)
        ADJtemp <- NULL
    }
    
    
    #GRAO<-igraph::set.vertex.attribute(GRAO, "class", value=ClassAssignment.numeric)
    GRAO<-igraph::set.vertex.attribute(GRAO, "class", value=as.character(ClassAssignment))
    Cuse<-NULL
    
    
    
    
    ######## COMMUNITY DETECTION #########
    message("Detecting Graph Communities...", appendLF = FALSE)
    memb=comm.method(GRAO)
    if (!is.null(ncom)) {
        memb$membership=cut_at(memb,no=ncom)  
    }
    csize=table(memb$membership)
    ConfMatrix <- table(predicted=memb$membership, true=ClassAssignment)
    message("done")

    V(GRAO)$membership <- memb$membership
    V(GRAO)$community.size <- csize[memb$membership]
    E(GRAO)$weight <- edge.weights(memb, GRAO, weight.within=2, weight.between=0.5)
    
    GRAO<-igraph::set.vertex.attribute(GRAO, "labels", value=CellIds)
    GRAOp<-NULL
    
    
    ######### Optimal Mapping between true and estimated class assignments: ##########
    mapping <- mapLabelsGreedy(memb$membership, ClassAssignment)
    misclErr <- classError(memb$membership, ClassAssignment, mapping)
    
    
    #### Add back Cell Ids to igraph object, ADJ, MEMB and prepare return value
    dimnames(ADJ) <- list(CellIds,CellIds)
    names(memb$membership) <- CellIds
    ret <- list(MEMB=memb$membership, MEMB.true=ClassAssignment,
                DISTM=ADJ, specp=NULL, ConfMatrix=ConfMatrix,
                miscl=misclErr, GRAO=GRAO, plotGRAO=NULL)

        
    ######### graph visualization
    if (plotG)
        ret[["plotGRAO"]] <- plotGraph(ret, maxG = maxG, fsuffix = fsuffix,
                                       image.format = image.format, quiet = FALSE)

    
    
    ########## SPECTRAL PROJECTION ########
    if (plotSP) {
        message("Computing Spectral Projection and Rendering...")
        colbar  <- gg_color_hue(  length(table(memb$membership) ) )
        comps <-memb$membership
        symbar <- c(21,24,22,25,23,c(0:14))
        class=ClassAssignment.numeric
        A=ADJ
        diag(A)=max(A)
        
        
        "%^%" <- function(M, power)
            with(eigen(M), vectors %*% (values^power * solve(vectors)))
        
        D <- diag(apply(A, 1, sum)) # sum rows
        U <- D-A
        #L <- diag(nrow(A)) - solve(D) %*% A #simple Laplacian
        #L <- (D %^% (-1/2)) %*% U %*% (D %^% (-1/2))  # normalized Laplacian
        L <-solve(D) %*% U #Generalized Laplacian
        
        k   <- 2
        evL <- eigen(U, symmetric=TRUE)
        Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
        Z=apply(Z,2,function(x) (x-min(x))/ (max(x)-min(x)) )
        #Z=apply(Z,2,function(x) jitter(x,amount=max(IQR(x),0.25) /10 ))
        D<-NULL; U<-NULL; L<-NULL; evL<-NULL
        
        mZ1=mean(Z[,1])
        mZ2=mean(Z[,2])
        outl=which(  abs(Z[,1]-mZ1) > 0.9 | abs(Z[,2]-mZ2)>0.9    )
        #message(length(outl),"\r")
        flush.console() 
        
        
        if (length(outl)>0 && length(outl) < 8 ) {
            A=A[-outl,-outl]
            comps=comps[-outl]
            class=class[-outl]
            D <- diag(apply(A, 1, sum)) # sum rows
            U <- D-A
            L <-solve(D) %*% U #Generalized Laplacian
            evL <- eigen(U, symmetric=TRUE)
            Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
            Z=apply(Z,2,function(x) (x-min(x))/ (max(x)-min(x)) )
            Z=apply(Z,2,function(x) jitter(x,amount=max(IQR(x),0.25) /10 ))
            D<-NULL; U<-NULL; L<-NULL; evL<-NULL
        }
        A<-NULL;
        
        #PCA=prcomp(x=Z,scale=T) 
        
        if(image.format=='pdf'){
            fname=paste('specProj_',fsuffix,'.pdf',sep="")
            pdf(fname,12,10)
        }
        else {
            fname=paste('specProj_',fsuffix,'.png',sep="")
            png(fname,12,10,units="in",res=300)   
        }
        
        par(mar=c(5.1, 4.1, 4.1, 14.1), xpd=TRUE)
        
        #plot(PCA$x[,1:2], col=colbar[comps], pch=symbar[class],cex=10/(nrow(Z)^0.35 ),bty='L',xlab="SD1",ylab="SD2",
        #     main="Spectral Projection")
        
        plot(Z, col=colbar[comps], pch=symbar[class],cex=10/(nrow(Z)^0.35 ),bty='L',xlab="SD1",ylab="SD2",
             main="Spectral Projection")
        
        
        legend("topright",inset=c(-0.35,0),title=" Pred.      True        ",
               legend= c(sort(unique(comps)),"" ,sort(unique(class))   ), 
               col=c(colbar[sort(unique(comps))],"white", rep("black",length(unique(class)) ) ), 
               pch=c(rep(20 ,length(unique(comps)) ),1, symbar[sort(unique(class))]),
               bty="n", border=F, ncol=2, text.width=0.02)
        dev.off()
        
        ret[["specp"]] <- Z
    }
    
    #########################################
    Te=(proc.time() - ptm)[3]
    Te=signif(Te,digits=6)
    message("Finished (Elapsed Time: ", Te, ")")

    
    ##### Stop registered cluster:
    if (isTRUE(use.par) & getDoParRegistered()) {  
        stopCluster(cl)
    }
    
    return(ret)
}


#' @title Visualize griph result as a graph.
#' 
#' @description Plot a graph obtained from \code{\link{SC_cluster}}, allowing to
#'     apply graph subsampling and control coloring.
#' 
#' @param gr A \code{griph} result, as returned by \code{\link{SC_cluster}}.
#' @param maxG Approximate maximal number of vertices to include when plotting the graph.
#' @param fill.type Type of fill color, one of \code{predicted} (predicted class labels, default),
#'     \code{true} (true class labels, if available) or \code{none} (no fill color).
#' @param line.type Type of line color, one of \code{true} (true class labels, if available, default),
#'     \code{predicted} (predicted class labels) or \code{none} no fill color.
#' @param fill.col Color vector defining the palette to use for vertex fill coloring.
#' @param line.col Color vector defining the palette to use for vertex outline coloring.
#' @param seed Random number seed to make graph layout deterministic.
#' @param fsuffix A suffix added to the file names of output plots. If not given
#'     it will use a random 5 character string. Ignored if \code{image.format} is \code{NULL}.
#' @param image.format Specifies the format of the created image. Currently supported are
#'     \code{\link{pdf}}, \code{\link{png}} or \code{NA}. If \code{NA} (the default), the plot
#'     is rendered on the currently opened plotting device.
#' @param forceRecalculation If \code{TRUE}, recalculate plotting-optimized graph
#'     even if it is already contained in \code{gr}.
#' @param quiet If \code{TRUE}, do not report on progress.
#' 
#' @return The plot-optimized version of the graph as an \code{igraph} object.
plotGraph <- function(gr, maxG=2500,
                      fill.type=c("predicted","true","none"), line.type=c("true","predicted","none"),
                      fill.col=c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#FFFFBF","#E6F598","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"),
                      #line.col=c("gold","maroon","green","blue","red","black","purple","darkorange","darkslategray","brown"),
                      line.col=c("#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F"),
                      seed=91919,
                      fsuffix=RandString(), image.format=NA,
                      forceRecalculation=FALSE, quiet=FALSE) {
    message("Computing Graph Layout and Rendering...")
    
    # get varaibles from gr
    GRAO <- gr$GRAO
    csize <- table(gr$MEMB)
    fill.type <- match.arg(fill.type)
    line.type <- match.arg(line.type)
    pct <- 1
    
    # get plot-optimized graph
    if(is.null(gr$plotGRAO) || forceRecalculation) {
        if (length(V(GRAO) ) > 1.25*maxG ) {
            if(!quiet) {
                message("\tRemark: Graph too large (>",maxG, " vertices). A sampled subgraph of ~", maxG, " vertices will be plotted")
                #message("\t$GRAO element of the results list will contain the complete graph object")
            }
    
            ###### Sample well-connected seeds from the members of each community 
            DEG <- igraph::degree(GRAO)
            snowball_seeds <- c()
            for (i in 1:length(csize)){
                if (csize[i]<5) {next}
                members <- which(gr$MEMB==i)
                minDegree <- quantile(DEG[members])[2]-1
                maxDegree <- quantile(DEG[members])[4]+1
                #seedn <- ceiling(csize[i]/max(csize))
                seedn <- ceiling (5 *(    sqrt(csize[i]-4)   /sqrt( max(csize)-4)   )  )
                seedn <- min(seedn,floor(csize[i]/4) )
                #message("minD:",minDegree, " maxD:",maxDegree," csize:", csize[i]  ,"\r")
                if (seedn > 1){
                    module_seeds <- sample(which(gr$MEMB==i & DEG >= floor(minDegree) &
                                                 DEG <= ceiling(maxDegree) ), seedn)
                } else {
                    module_seeds <- sample(which(gr$MEMB==i & DEG==max(DEG[members])), 1)    
                }
        
                snowball_seeds <- unique(c(snowball_seeds,module_seeds))
            }
    
            snowball <- c()
            seed.ego_size <- 0
            while (length(snowball) < maxG/2){
                seed.ego_size <- seed.ego_size+1  
                snowball <- unique(unlist(igraph::ego(GRAO,seed.ego_size,snowball_seeds)))
            }
            if (length(snowball) > 1.25*maxG && seed.ego_size > 1 ) {
                seed.ego_size <- seed.ego_size-1
                snowball <- unique(unlist(igraph::ego(GRAO,seed.ego_size,snowball_seeds))) 
            }
    
            GRAOp <- igraph::induced.subgraph(GRAO,sort(snowball) )
    
            if(!quiet)
                message("\tUsed vertices: ", length(V(GRAOp)),"  seed_size: ",seed.ego_size)
    
        } else {
            GRAOp <- GRAO
        }

        ######## Prune graph for better plot output
        if (median(igraph::degree(GRAOp)) > 4 ) {
            pct <- min(1,1/sqrt(0.1*median(igraph::degree(GRAOp)) ) )
            ADJp <- as.matrix(igraph::get.adjacency(GRAOp,attr='weight'))
            ADJp <- apply(ADJp,1,function(x) sparsify(x,pct) )
            ADJp[which(abs(ADJp) >0)] <- 1
            GRAOtemp <- igraph::graph.adjacency(ADJp,mode=c("max"),weighted=NULL,diag=FALSE)
            GRAOp <- igraph::intersection(GRAOp,GRAOtemp,byname=FALSE)
            GRAOtemp <- NULL
            ADJp <- NULL
        }

        GRAOp <- igraph::delete_vertices(GRAOp,which(igraph::ego_size(GRAOp,3) < 6))
        ###Delete Vertices from communites with few members:
        min.csize <- ceiling(0.25*sqrt(length(V(GRAO))))
        GRAOp <- igraph::delete_vertices(GRAOp,which( V(GRAOp)$community.size < min.csize ))  

        if(!quiet)
            message("\tRemark: Nodes from communities with <",min.csize, " members will not be displayed.")
        
    } else {
        if(!quiet)
            message("using existing plot-optimized graph")
        GRAOp <- gr$plotGRAO
    }
    if(!quiet)
        message("\tdisplaying ",round(100*pct,1), "% of edges")

    # get colors
    class.pred <- factor(V(GRAOp)$membership, levels=sort(unique(V(GRAOp)$membership)))
    #class.pred.numeric <- as.numeric(class.pred)
    class.true <- factor(V(GRAOp)$class, levels=unique(V(GRAOp)$class))

    fillColorPalette <- switch(fill.type,
                               predicted=grDevices::colorRampPalette(fill.col)(nlevels(class.pred)),
                               true=grDevices::colorRampPalette(fill.col)(nlevels(class.true)),
                               none=NA)
    fillColor <- switch(fill.type,
                        predicted=fillColorPalette[as.numeric(class.pred)],
                        true=fillColorPalette[as.numeric(class.true)],
                        none=rep(NA, length(V(GRAOp))))
    lineColorPalette <- switch(line.type,
                               predicted=grDevices::colorRampPalette(line.col)(nlevels(class.pred)),
                               true=grDevices::colorRampPalette(line.col)(nlevels(class.true)),
                               none=NA)
    lineColor <- switch(line.type,
                        predicted=lineColorPalette[as.numeric(class.pred)],
                        true=lineColorPalette[as.numeric(class.true)],
                        none=rep("black", length(V(GRAOp))))
    
    # set some more graph attributes
    #V(GRAOp)$classcolor <- line.col[V(GRAOp)$class]
    V(GRAOp)$classcolor <- lineColor
    V(GRAOp)$size <- 10 /(length(V(GRAOp))/60 )^0.3
    V(GRAOp)$cex <- V(GRAOp)$size / 3
    V(GRAOp)$frame.width <- 2 /(length(V(GRAOp))/60 )^0.3
    E(GRAOp)$width <- E(GRAOp)$weight / sqrt((length(V(GRAOp))/60 ))
    #colbar  <- gg_color_hue(  length(table(V(GRAOp)$membership) ) )
    #colbar <- grDevices::colorRampPalette(fill.col)(length(table(V(GRAOp)$membership)))
    #names(colbar) <- paste("c",sort(unique(V(GRAOp)$membership)),sep=""  )
    #V(GRAOp)$color <- colbar[   paste("c",V(GRAOp)$membership,sep=""  )   ]
    V(GRAOp)$color <- fillColor
    
    set.seed(seed = seed)
    l <- igraph::layout_with_fr(graph = GRAOp)
    # igraph::add.vertex.shape("fcircle", clip=igraph.shape.noclip, plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))

    if(!is.na(image.format)) {
        if(image.format=='pdf') {
            fname <- paste('graph_',fsuffix,'.pdf',sep="")
            pdf(file = fname, width = 12, height = 10)
        } else if (image.format=="png") {
            fname <- paste('graph_',fsuffix,'.png',sep="")
            png(filename = fname, width = 12, height = 10, units = "in", res = 300)   
        }
        if(!quiet)
            message("\tsaving graph to ",fname)
    }

    par(mar=c(5.1, 4.1, 4.1, 14.1), xpd=TRUE)
    #igraph::plot.igraph(GRAOp, layout=l, asp=0, vertex.label=NA, edge.lty=0,
    #                    vertex.frame.color=igraph::V(GRAOp)$classcolor, vertex.shape="fcircle",
    #                    vertex.frame.width=igraph::V(GRAOp)$frame.width, edge.curved=TRUE )
    #legend("topright",inset=c(-0.35,0),title=" Pred.      True        ",
    #       legend= c(sort(unique(V(GRAOp)$membership)),"" ,sort(unique(class.pred.numeric))   ), 
    #       col=c(colbar,"white", line.col[1:length(unique(class.pred.numeric))] ), 
    #       pch=c(rep(20 ,length(unique( V(GRAOp)$membership  )) ),1, rep(21 ,length(unique(class.pred.numeric))  )   ),
    #       bty="n", border=F, ncol=2, text.width=0.02)
    plot(l[,1], l[,2], col=lineColor, bg=fillColor, pch=21, cex=2.5, lwd=2, axes=FALSE, xlab="", ylab="")
    if(fill.type != "none") {
        lgd <- legend(x = par("usr")[2]+12*par("cxy")[1], y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
                      pch = 21, cex = 1, pt.cex = 2.5, col = if(line.type=="none") "black" else NA, pt.bg = fillColorPalette,
                      title = fill.type, legend = switch(fill.type,
                                                         predicted=levels(class.pred),
                                                         true=levels(class.true)))
    } else {
        lgd <- list(rect=list(left=par("usr")[2]+12*par("cxy")[1]))
    }
    if(line.type != "none") {
        legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
               pch = 21, cex = 1, pt.cex = 2.5, col = lineColorPalette, pt.bg = "white",
               title = line.type, legend = switch(line.type,
                                                  predicted=levels(class.pred),
                                                  true=levels(class.true)))
    }

    if(!is.na(image.format))
        dev.off()
    
    return(GRAOp)
}
