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
#' @param is.cor If \code{TRUE}, \code{DM} is assumed to be a correlation matrix.
#'     Otherwise (the default), \code{DM} is assumed to be a genes-by-cells count
#'     matrix, and a correlation matrix will be computed.
#' @param ref.iter Number of clustering refinement iterations.  
#' @param use.par If \code{TRUE}, use parallel versions of distance calculation
#'     functions based on \code{\link[foreach]{foreach}} (see details).
#' @param ncores a numeric(1) or character(1), either specifying the number of
#'     parallel jobs to use, or \code{"all"} (default) to use up to 90 percent of
#'     the available cores. Ignored if \code{use.par==FALSE}.
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
#' 



griph_cluster <- function(DM, is.cor = FALSE,ref.iter=0,use.par=FALSE,ncores="all",
                          impute = FALSE, filter = FALSE, rho = 0.25, pr.iter = 1, batch.penalty = 0.5,
                          ClassAssignment = rep(1,ncol(DM)), BatchAssignment = NULL,
                          plotG = TRUE, maxG = 2500, fsuffix = RandString(), image.format='png'){
    
    ptm=proc.time() #Start clock
    
    SpearmanCor=Spcor
    
    if (length(ClassAssignment) != ncol(DM))
        stop ("length(ClassAssignment) must be equal to ncol(DM)")
    if(!is.null(BatchAssignment) && length(BatchAssignment) != ncol(DM))
        stop ("length(BatchAssignment) must be equal to ncol(DM)")
    
    # Register cluster here, remove regiastation block from SC_cluster
    if (isTRUE(use.par)) {
        #######Switch to parallelized functions if use.par=TRUE
        SpearmanCor=FlashSpearmanCor
        
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
    
    ### wrap code in tryCatch block, ensuring that stopCluster(cl) is called even when a condition is raised  
    tryCatch({    
        for (i in 0:ref.iter) { 
            if (i<ref.iter){
                params$plotG=FALSE   
            }
            else {params$plotG=plotG}
            
            if (i==0) {
                params$DM=DM
                params$pr.iter=1
                cluster.res <- do.call("SC_cluster2",params)
            }
            else {
                message("\n\nRefining Cluster Structure...\n", appendLF = FALSE)
                params$is.cor=TRUE
                params$pr.iter=0
                ####### construct cell2cell correlation matrix using the current cluster.res: ########
                memb=cluster.res$MEMB
                min.csize <-max(4, ceiling(0.25*sqrt(length(memb)) ) )
                nclust=length(unique(memb) )
                good.clust=as.vector(which(table(memb)>=min.csize) )
                if (length(good.clust)<3){
                message("\nToo few substantial clusters (<3). Using fake bulks to refine clusters not possible\n Reverting to previous iteration...\n", appendLF = FALSE)
                break  
                }

                message("\nUsing ", length(good.clust) ," fake bulks to refine clusters...\n", appendLF = FALSE)
 
                FakeBulk=matrix(0,nrow(DM),length(good.clust))
                for (c in 1:length(good.clust)) {
                    clust=good.clust[c]
                    FakeBulk[,c]=rowSums(DM[,memb==clust])
                }
                learnStep=1-(1/sqrt(length(good.clust)))
                
                params$DM=((1-learnStep)*cluster.res$CORM+learnStep*SpearmanCor(cor(log2(FakeBulk+1),log2(DM+1)   )))
                cluster.res <- do.call("SC_cluster2",c(params,list(comm.method=igraph::cluster_louvain) ) )
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

