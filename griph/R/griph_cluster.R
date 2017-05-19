#' Combined distance calculations.
#' 
#' @description Calculate and combine several between-cell distance measures to
#'     obtain robust cell-to-cell distances.
#' 
#' @param M gene-by-cell count matrix.
#' @param FB gene-by-FakeBulks gene sums count matrix.
#' @param ShrinkCor Function to calculate shrinkage correlation.
#' @return cell-by-cell distance matrix.
WScorFB <- function (M,FB, ShrinkCor=ShrinkCor   ) {
    message(ncol(M),nrow(M),"\r")
    message(ncol(FB),nrow(FB),"\r")
    
    CellIds=colnames(M)
    dimnames(M)=NULL
    message("1","\n")
    D=cor(log2(FB+1),log2(M+1))
    R=vapply(c(1:ncol(D)),function (x) rank(D[,x]),FUN.VALUE=double(length=nrow(D) ) )  #pearson's cor    
    
    ######## Counts per 10K:
    CellCounts=colSums(FB)
    FB=sweep(FB,2,CellCounts,FUN="/")
    FB=FB*10000
    CellCounts=colSums(M)
    M=sweep(M,2,CellCounts,FUN="/")
    M=M*10000 
    
    Dt=PCanberraMat( log2(FB+1),log2(M+1) )   #
    Dt=1-( (Dt-min(Dt))/ diff(range(Dt)) )
    D=D+Dt
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #canberra
    Dt=cor(FB,M,method="spearman")
    D=D+Dt
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #spearman's cor 
    Dt=PHellingerMat(FB,M)
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
        R=sweep(R,2,colMeans(R),"-")
        #R=R*(W^0.4)
        R=buildEdgeMatrix( R ,distance_method="Cosine",K=ceiling(1e2+sqrt(ncol(R))  )  )
        R=sparseMatrix(i=R$i,j=R$j,x=1-(R$x/2),dims=attr(R,"dims"),dimnames=list(CellIds,CellIds))

    return(R)    
}



















#' @title Wrapper function for SC_cluster
#' 
#' @description \code{griph_cluster} takes a gene-by-cell matrix with read counts as
#'     input and groups individual cells into clusters of similar cell (cell types).
#'     The number of cell types is automatically determined.
#' 
#' @details 
#' TODO
#' 
#' @param DM Count data (n genes-by-k cells) or directly a correlation k-by-k matrix.
#'     Required argument.
#' @param SamplingSize Number of sampled cells in initialization step. 
#' @param ref.iter Number of clustering refinement iterations.  If set to 0 only the clustering initialization 
#'     step is performed to \code{min(SamplingSize,ncol(DM))} cells. 
#' @param use.par If \code{TRUE}, use parallel versions of distance calculation
#'     functions based on \code{\link[foreach]{foreach}} (see details).
#' @param ncores a numeric(1) or character(1), either specifying the number of
#'     parallel jobs to use, or \code{"all"} (default) to use up to 90 percent of
#'     the available cores. Ignored if \code{use.par==FALSE}.
#' @param filter T/F Filter genes according to cv=f( mean ) fitting. Default: \code{TRUE}.
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
#' @param ncom Forces detection of a fixed number of clusters (communities).
#' @param plotG if \code{TRUE} plots the resulting graph
#' @param maxG Approximate maximal number of vertices to include when plotting the graph.
#' @param fsuffix A suffix added to the file names of output plots. If not given
#'     it will use a random 5 character string.
#' @param image.format Specifies the format of the created images. Currently only pdf and png filetypes are supported.
#' 
#' @return Currently a list with the clustering results.
#' 



griph_cluster <- function(DM, SamplingSize=750,ref.iter=1,use.par=FALSE,ncores="all",
                          filter = FALSE, rho = 0.25, batch.penalty = 0.5,
                          ClassAssignment = rep(1,ncol(DM)), BatchAssignment = NULL, ncom=NULL,
                          plotG = TRUE, maxG = 2500, fsuffix = RandString(), image.format='png'){
    
    ptm=proc.time() #Start clock
    
    #######Define functions if use.par=FALSE
    SpearmanCor=Spcor
    ShrinkCor=corpcor::cor.shrink
    
    if (length(ClassAssignment) != ncol(DM))
        stop ("length(ClassAssignment) must be equal to ncol(DM)")
    if(!is.null(BatchAssignment) && length(BatchAssignment) != ncol(DM))
        stop ("length(BatchAssignment) must be equal to ncol(DM)")
    
    # Register cluster here, remove regiastation block from SC_cluster
    if (isTRUE(use.par)) {
        #######Switch to parallelized functions if use.par=TRUE
        SpearmanCor=FlashSpearmanCor
        ShrinkCor=FlashShrinkCor
        
        if(ncores=="all"){
            ncores = parallel::detectCores()
            ncores=min(48,floor(0.9*ncores),ceiling(ncol(DM)/200))
        } else{
            ncores=min(48,ncores,floor(0.9*parallel::detectCores()),ceiling(ncol(DM)/200))
        }
        cl<-parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
    }
    
    params <-as.list(environment())
    params$plotG=FALSE  
    ### wrap code in tryCatch block, ensuring that stopCluster(cl) is called even when a condition is raised  
    tryCatch({    
        for (i in 0:ref.iter) { 
            if (i==0) {

                #Set the number of initialization clusters to smth reasonable given the number of cells:
                params$ncom=min(0.5*ceiling(sqrt(ncol(DM))),16)
                params$ncom=max(params$ncom,5)

                if (ref.iter==0){
                params$ncom=ncom    
                }
                
                if (ncol(DM)>SamplingSize){
                SMPL=sample(1:ncol(DM),SamplingSize)
                params$DM=DM[,SMPL]
                params$ClassAssignment=ClassAssignment[SMPL]
                    if (!is.null(BatchAssignment)){
                    params$BatchAssignment=BatchAssignment[SMPL]   
                    }
                }
                
                else{
                params$DM=DM
                }
                
                cluster.res <- do.call(SC_cluster, c( params,list(comm.method=igraph::cluster_louvain,pr.iter=1 ) )     )
            }
            else {
                message("\n\nRefining Cluster Structure...\n", appendLF = FALSE)
                params$ncom=ncom
                params$is.cor=TRUE
                params$ClassAssignment=ClassAssignment
                params$BatchAssignment=BatchAssignment  
                
                message("MISCL","\n",cluster.res$miscl,"\n")
                
                ####### construct cell2cell correlation matrix using the current cluster.res: ########
                memb=cluster.res$MEMB
                min.csize <-max(4, ceiling(0.2*sqrt(length(memb)) ) )
                nclust=length(unique(memb) )
                good.clust=as.vector(which(table(memb)>=min.csize) )

                if (length(good.clust)<3){
                    message("\nToo few (<3) substantial clusters. Using fake bulks to refine clusters not possible\n Reverting to previous iteration...\n", appendLF = FALSE)
                    break  
                }
                
                message("\nUsing ", length(good.clust) ," fake bulks to refine clusters...\n", appendLF = FALSE)
                

                FakeBulk=matrix(0,nrow(DM),length(good.clust))
                for (c in 1:length(good.clust)) {
                    clust=good.clust[c]
                    FakeBulk[,c]=rowSums(DM[,names(memb)][,memb==clust])
                }
                
                ###### Calculate distances of all the cells to the FakeBulks:
                message("Calculating Cell Distances to Cluster Centroids (Bulks)...", appendLF = FALSE)
                params$DM=WScorFB(DM[,names(memb)],FakeBulk,ShrinkCor=ShrinkCor)
                message("done")
                
                cluster.res <- do.call(SC_cluster, c(params,list(comm.method=igraph::cluster_louvain,do.glasso=FALSE,pr.iter=0) ) )
            }
            gc() #Call garbage collector
        }
        
        if (plotG==TRUE){    
            plotGraph(cluster.res, maxG = maxG, fsuffix = fsuffix,image.format = image.format, quiet = FALSE)
        }
        
    }, # end of tryCatch expression, cluster object cl not needed anymore    
    finally = { 
        ##### Stop registered cluster:
        if (isTRUE(use.par) & foreach::getDoParRegistered())
            parallel::stopCluster(cl)
    })
    
    
    ################## Stop the clock #######################
    Te=(proc.time() - ptm)[3]
    Te=signif(Te,digits=6)
    message("Finished (Elapsed Time: ", Te, ")")
    cluster.res$CORM=NULL
    return(cluster.res)    
}

