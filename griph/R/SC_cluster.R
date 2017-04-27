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
#' @param ClassAssignment If available a numeric vector of length \code{k} with numeric
#'     class labels (e.g-> c(1,2,1,1,1,2,3,3,3,1,2))
#' @param BatchAssignment If available a numeric vector of length \code{k} with numeric
#'     batch labels (e.g-> c(1,2,1,1,1,2,3,3,3,1,2))
#' @param plotG if \code{TRUE} plots the resulting graph
#' @param maxG Approximate maximal number of vertices to include when plotting the graph.
#' @param fsuffix A suffix added to the file names of output plots. If not given
#'     it will use a random 5 character string.
#' @param image.format Specifies the format of the created images. Currently only pdf and png filetypes are supported.

#' 
#' @return Currently a list with the clustering results.

SC_cluster <- function(DM, use.par=FALSE,ncores="all",is.cor = FALSE,
                       impute = FALSE, filter = FALSE, rho = 0.25, pr.iter = 1, batch.penalty = 0.5,
                       ClassAssignment = rep(1,ncol(DM)), BatchAssignment = NULL,
                       plotG = TRUE, maxG = 2500, fsuffix = RandString(), image.format='png' ){
    
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
    
    
    if (length(ClassAssignment) != ncol(DM))
        stop ("length(ClassAssignment) must be equal to ncol(DM)")
    if(!is.null(BatchAssignment) && length(BatchAssignment) != ncol(DM))
        stop ("length(BatchAssignment) must be equal to ncol(DM)")

    ptm=proc.time()
    

    ### wrap code in tryCatch block, ensuring that stopCluster(cl) is called even when a condition is raised
    tryCatch({

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

            message("done")

            ##### Strip dimnames:
            CellIds=colnames(DM)
            dimnames(DM)=NULL
            
            C=list()
            C[[1]]=PearsonCor(log2(DM+1))
            
            message("Calculating Pairwise and Diffused Similarities...", appendLF = FALSE)
            C[[2]]=WScor(nDM,C1=C[[1]], CanberraDist=CanberraDist,
                         SpearmanCor=SpearmanCor, HellingerDist=HellingerDist, ShrinkCor=ShrinkCor)
            nDM<-NULL
            
            message("done")
        }
        
        
        else {
            C=list()
            ####Normalize Initial Correlation matrix:
            cmDM <- colMeans(DM)
            W <- pmax(1e-1, cmDM) / mean(cmDM)
            W <- sqrt(W) %o% sqrt(W)
            C[[2]]=( DM / W )
        }
        
        Cuse=as.matrix(C[[2]])

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
        
        if(pr.iter > 0) {
            for (i in 1:pr.iter){
                message("Pruning based on global node similarity: ",i," / ",pr.iter, "\r", appendLF = FALSE)
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
        
    }, # end of tryCatch expression, cluster object cl not needed anymore
    
    # exception handlers could be defined here
    
    finally = {
        ##### Stop registered cluster:
        if (isTRUE(use.par) & getDoParRegistered())
            stopCluster(cl)
    })

    
    ######## COMMUNITY DETECTION #########
    message("Detecting Graph Communities...", appendLF = FALSE)
    memb=comm.method(GRAO)
    if (!is.null(ncom)) {
        if(is_hierarchical(memb)) {
            memb$membership=cut_at(memb,no=ncom)  
        } else {
            warning("ignoring 'ncom' when using non-hierachical community detection algorithm")
        }
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
                DISTM=ADJ, ConfMatrix=ConfMatrix,
                miscl=misclErr, GRAO=GRAO, plotGRAO=NULL)

        
    ######### graph visualization
    if (plotG)
        ret[["plotGRAO"]] <- plotGraph(ret, maxG = maxG, fsuffix = fsuffix,
                                       image.format = image.format, quiet = FALSE)


    #########################################
    Te=(proc.time() - ptm)[3]
    Te=signif(Te,digits=6)
    message("Finished (Elapsed Time: ", Te, ")")

    
    return(ret)
}


#' @title Visualize griph result as a graph.
#' 
#' @description Plot a graph obtained from \code{\link{SC_cluster}}, allowing to
#'     apply graph subsampling and control coloring.
#' 
#' @param gr A \code{griph} result, as returned by \code{\link{SC_cluster}}.
#' @param maxG Approximate maximal number of vertices to include when plotting
#'     the graph.
#' @param fill.type Type of fill color, one of \code{predicted} (predicted class
#'     labels, default), \code{true} (true class labels, if available), \code{none}
#'     (no fill color) or \code{custom} (use \code{custom.col}).
#' @param line.type Type of line color, one of \code{true} (true class labels, if
#'     available, default), \code{predicted} (predicted class labels), \code{none}
#'     for no fill color or \code{custom} (use \code{custom.col}).
#' @param mark.type Type of cell class defnition to mark using polygons,
#'     one of \code{none} (no polygons, the default), \code{predicted} (draw
#'     polygons around cells with the same predicted class label), \code{true}
#'     (polygons around cells with the same true class label, if available) or
#'     \code{custom} (polygons around cells with the same \code{custom.class}).
#' @param collapse.type Type of cell class to use for graph simplification,
#'     by combining cells of the same class into a single vertex. If set to a value
#'     other than \code{"none"}, the same value will also be used for \code{fill.type}
#'     and \code{line.type} is ignored.
#' @param fill.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for vertex fill coloring.
#' @param line.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for vertex outline coloring.
#' @param mark.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for cell class polygon marking.
#' @param custom.class Factor, character or numberic vector of the same length or
#'     with names corresponding to names(gr$MEMB) to use for custom cell classification
#'     (used if \code{fill.type} and/or \code{line.type} is set to "custom").
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
                      fill.type=c("predicted","true","none","custom"),
                      line.type=c("true","predicted","none","custom"),
                      mark.type=c("none","predicted","true","custom"),
                      collapse.type=c("none","predicted","true"),
                      fill.col="Spectral",
                      line.col="Dark2",
                      mark.col="Pastel1",
                      custom.class=factor(rep(1, length(gr$MEMB))),
                      seed=91919,
                      fsuffix=RandString(), image.format=NA,
                      forceRecalculation=FALSE, quiet=FALSE) {
    if(!quiet)
        message("Computing Graph Layout and Rendering...")
    
    # get varaibles from gr
    GRAO <- gr$GRAO
    MEMB <- gr$MEMB
    csize <- table(MEMB)
    
    # digest arguments
    fill.type <- match.arg(fill.type)
    line.type <- match.arg(line.type)
    mark.type <- match.arg(mark.type)
    collapse.type <- match.arg(collapse.type)
    if(is.null(names(custom.class))) {
        stopifnot(length(custom.class)==length(MEMB))
        names(custom.class) <- names(MEMB)
    }
    if(length(fill.col)==1) {
        stopifnot(fill.col %in% rownames(RColorBrewer::brewer.pal.info))
        fill.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[fill.col, "maxcolors"], fill.col)
    }
    if(length(line.col)==1) {
        stopifnot(line.col %in% rownames(RColorBrewer::brewer.pal.info))
        line.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[line.col, "maxcolors"], line.col)
    }
    if(length(mark.col)==1) {
        stopifnot(mark.col %in% rownames(RColorBrewer::brewer.pal.info))
        mark.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[mark.col, "maxcolors"], mark.col)
    }

    # global plotting paramterers
    pct <- 1
    my.pch <- 21L # should be in 21:25
    my.pt.cex <- 2.5
    my.pt.lwd <- if(line.type == "none" || collapse.type != "none") 1.0 else 2.5
    edge.lwd.max <- 12.0
    edge.col <- "#33333355"
    
    # get plot-optimized graph
    if(collapse.type != "none") {
        if(!quiet)
            message("\tcollapsing cells by ",collapse.type," class")
        if(fill.type != collapse.type && !quiet)
            message("\tSetting fill.type to collapse.type for collapse.type != 'none'.")
        if((line.type != "none" || mark.type != "none") && !quiet)
            message("\tSetting line.type and mark.type to 'none' for collapse.type != 'none'.")
        fill.type <- collapse.type
        line.type <- "none"
        mark.type <- "none"
        class.collapse <- switch(collapse.type,
                                 predicted=factor(V(GRAO)$membership, levels=sort(as.numeric(unique(V(GRAO)$membership)))),
                                 true=factor(V(GRAO)$class, levels=unique(V(GRAO)$class)))
        GRAOp <- igraph::simplify( igraph::contract(GRAO, class.collapse) )
        V(GRAOp)$membership <- V(GRAOp)$class <- levels(class.collapse)
        V(GRAOp)$size <- 10 /(length(V(GRAOp))/60 )^0.3 *(as.numeric(csize/median(csize)))^0.5
        
    } else if(is.null(gr$plotGRAO) || forceRecalculation) {
        if (length(V(GRAO) ) > 1.25*maxG ) {
            if(!quiet)
                message("\tRemark: Graph too large (>",maxG, " vertices). A sampled subgraph of ~", maxG, " vertices will be plotted")

            ###### Sample well-connected seeds from the members of each community 
            DEG <- igraph::degree(GRAO)
            snowball_seeds <- c()
            for (i in 1:length(csize)){
                if (csize[i]<5) {next}
                members <- which(MEMB==i)
                minDegree <- quantile(DEG[members])[2]-1
                maxDegree <- quantile(DEG[members])[4]+1
                #seedn <- ceiling(csize[i]/max(csize))
                seedn <- ceiling (5 *(    sqrt(csize[i]-4)   /sqrt( max(csize)-4)   )  )
                seedn <- min(seedn,floor(csize[i]/4) )
                #message("minD:",minDegree, " maxD:",maxDegree," csize:", csize[i]  ,"\r")
                if (seedn > 1){
                    module_seeds <- sample(which(MEMB==i & DEG >= floor(minDegree) &
                                                 DEG <= ceiling(maxDegree) ), seedn)
                } else {
                    module_seeds <- sample(which(MEMB==i & DEG==max(DEG[members])), 1)    
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
    class.pred <- factor(V(GRAOp)$membership, levels=sort(as.numeric(unique(V(GRAO)$membership))))
    class.true <- factor(V(GRAOp)$class, levels=unique(V(GRAO)$class))
    class.none <- factor(rep(NA, length(V(GRAOp))))
    class.custom <- factor(custom.class[V(GRAOp)$labels], levels=unique(custom.class))

    class.fill <- switch(fill.type,
                         predicted=class.pred,
                         true=class.true,
                         none=class.none,
                         custom=class.custom)
    fillColorPalette <- if(nlevels(class.fill) > length(fill.col)) grDevices::colorRampPalette(fill.col)(nlevels(class.fill)) else fill.col
    fillColor <- fillColorPalette[as.numeric(class.fill)]

    class.line <- switch(line.type,
                         predicted=class.pred,
                         true=class.true,
                         none=class.none,
                         custom=class.custom)
    lineColorPalette <- if(nlevels(class.line) > length(line.col)) grDevices::colorRampPalette(line.col)(nlevels(class.line)) else line.col
    lineColor <- lineColorPalette[as.numeric(class.line)]

    class.mark <- switch(mark.type,
                         predicted=class.pred,
                         true=class.true,
                         none=class.none,
                         custom=class.custom)
    markElements <- split(seq_along(class.mark), class.mark)
    markColor <- if(nlevels(class.mark) > length(mark.col)) grDevices::colorRampPalette(mark.col)(nlevels(class.mark)) else mark.col[1:nlevels(class.mark)]

    # set some more graph attributes
    V(GRAOp)$classcolor <- lineColor
    V(GRAOp)$size <- if(collapse.type == "none") 10 /(length(V(GRAOp))/60 )^0.3 else V(GRAOp)$size
    V(GRAOp)$cex <- V(GRAOp)$size / 3
    V(GRAOp)$frame.width <- 2 /(length(V(GRAOp))/60 )^0.3
    E(GRAOp)$width <- E(GRAOp)$weight / sqrt((length(V(GRAOp))/60 ))
    V(GRAOp)$color <- fillColor
    
    # compute graph layout
    set.seed(seed = seed)
    l <- igraph::layout_with_fr(graph = GRAOp)
    # igraph::add.vertex.shape("fcircle", clip=igraph.shape.noclip, plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))

    # open output file
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

    # setup plot coordinate system
    par(mar=c(5.1, 4.1, 4.1, 14.1), xpd=TRUE)
    plot(l[,1], l[,2], type = "n", axes=FALSE, xlab="", ylab="")
    
    # add mark polygons
    if(mark.type != "none") {
        for(j in which(lengths(markElements) > 0)) {
            xy <- l[markElements[[j]], , drop = FALSE]
            off <- par("cxy")[2]*1
            #avg <- matrix(colMeans(xy), ncol = 2, nrow = nrow(xy), byrow = TRUE)
            #pp <- xy + sign(xy - avg) * off
            pp <- rbind(xy,
                        cbind(xy[, 1] - off, xy[, 2]),
                        cbind(xy[, 1] + off, xy[, 2]),
                        cbind(xy[, 1], xy[, 2] - off), 
                        cbind(xy[, 1], xy[, 2] + off))
            cl <- igraph::convex_hull(pp)
            graphics::xspline(cl$rescoords, shape = 0.5, open = FALSE,
                              col = paste0(markColor[j], "66"),
                              border = adjust.color(markColor[j], 0.5))
            
        }
    }
    
    # add edges (currently only for collapase.type != "none")
    if(collapse.type != "none") {
        el <- igraph::as_edgelist(GRAOp, names = FALSE)
        graphics::segments(x0 = l[,1][el[,1]], y0 = l[,2][el[,1]],
                           x1 = l[,1][el[,2]], y1 = l[,2][el[,2]],
                           col = edge.col, lwd = E(GRAOp)$weight /max(E(GRAOp)$weight) * edge.lwd.max)
    }
    
    # add vertices
    points(l[,1], l[,2], col = if(line.type=="none") "black" else lineColor,
           bg = fillColor, pch = my.pch, lwd = my.pt.lwd,
           cex=my.pt.cex * if(collapse.type == "none") 1.0 else (as.numeric(csize/median(csize)))^0.5)

    # add legend(s)
    if(mark.type != "none") {
        lgd <- legend(x = par("usr")[2]+12*par("cxy")[1], y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
                      cex = 1, fill = paste0(markColor, "66"), title = mark.type, legend = levels(class.mark))
    } else {
        lgd <- list(rect=list(left=par("usr")[2]+12*par("cxy")[1]))
    }
    if(fill.type != "none" && nlevels(class.fill) > 0) {
        lgd <- legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
                      pch = my.pch, pt.lwd = my.pt.lwd, cex = 1, pt.cex = my.pt.cex,
                      col = if(line.type=="none") "black" else "white", pt.bg = fillColorPalette,
                      title = fill.type, legend = levels(class.fill))
    }
    if(line.type != "none" && nlevels(class.line) > 0) {
        legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
               pch = my.pch, pt.lwd = my.pt.lwd, cex = 1, pt.cex = my.pt.cex,
               col = lineColorPalette, pt.bg = "white",
               title = line.type, legend = levels(class.line))
    }

    # close output file
    if(!is.na(image.format))
        dev.off()
    
    return(GRAOp)
}
