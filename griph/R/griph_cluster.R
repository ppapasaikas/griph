#' Combined distance calculations.
#' 
#' @description Calculate and combine several between-cell distance measures to
#'     obtain robust cell-to-cell distances.
#' 
#' @param M gene-by-cell count matrix.
#' @param FB gene-by-FakeBulks gene sums count matrix.
#' @param PPearsonCor Function to calculate Pearson correlation between the columns of two matrices.
#' @param PSpearmanCor Function to calculate Spearman correlation between the columns of two matrices.
#' @param PHellinger Function to calculate Hellinger Distance between the columns of two matrices.
#' @param PCanberra Function to calculate Canberra Distance between the columns of two matrices.
#' @param ShrinkCor Function to calculate shrinkage correlation.
#' @return cell-by-cell distance matrix.
WScorFB <- function (M,FB, PSpearmanCor, PPearsonCor, PHellinger, PCanberra, ShrinkCor=ShrinkCor   ) {
    
    CellIds=colnames(M)
    dimnames(M)=NULL
    
    ptm1=proc.time() #Start clock
    
    D=PPearsonCor(log2(FB+1),log2(M+1))
    R=vapply(c(1:ncol(D)),function (x) rank(D[,x]),FUN.VALUE=double(length=nrow(D) ) )  #pearson's cor    
    
    Te1=signif((proc.time() - ptm1)[3],digits=6)
    message("\nP (Elapsed Time: ", Te1, ")")
    
    ######## Counts per 100K:
    CellCounts=colSums(FB)
    FB=sweep(FB,2,CellCounts,FUN="/")
    FB=FB*100000
    CellCounts=colSums(M)
    M=sweep(M,2,CellCounts,FUN="/")
    M=M*100000 
    
    ptm1=proc.time() #Start clock
    Dt=PCanberra( log2(FB+1),log2(M+1) )   
    Dt=1-( (Dt-min(Dt))/ diff(range(Dt)) )
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #canberra
    Te1=signif((proc.time() - ptm1)[3],digits=6)
    message("C (Elapsed Time: ", Te1, ")")
    
    ptm1=proc.time() #Start clock
    Dt=PSpearmanCor(FB,M)
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #spearman's cor 
    Te1=signif((proc.time() - ptm1)[3],digits=6)
    message("S (Elapsed Time: ", Te1, ")")
    
    ptm1=proc.time() #Start clock
    Dt=PHellinger(FB,M)
    Dt=1-( (Dt-min(Dt))/ diff(range(Dt)) )
    R=R+vapply(c(1:ncol(Dt)),function (x) rank(Dt[,x]),FUN.VALUE=double(length=nrow(Dt) ) )  #Hellinger distance
    Te1=signif((proc.time() - ptm1)[3],digits=6)
    message("H (Elapsed Time: ", Te1, ")")
    
    ptm1=proc.time() #Start clock
    R=(R/4)^2
    R=sweep(R,2,colMeans(R),"-")
    R=buildEdgeMatrix( R ,distance_method="Cosine",K=min( 250, ceiling(1e2+sqrt(ncol(R)) ) )   )
    R=Matrix::sparseMatrix(i=R$i,j=R$j,x=1-(R$x/2),dims=attr(R,"dims"),dimnames=list(CellIds,CellIds))
    Te1=signif((proc.time() - ptm1)[3],digits=6)
    message("bEM (Elapsed Time: ", Te1, ")")
    
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
#' @param SamplingSize Number of sampled cells in initialization step. If NULL it is set to \code{max(750,ncores*250)}
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



griph_cluster <- function(DM, SamplingSize= NULL,ref.iter=1,use.par=TRUE,ncores="all",
                          filter = TRUE, rho = 0.25, batch.penalty = 0.5,
                          ClassAssignment = rep(1,ncol(DM)), BatchAssignment = NULL, ncom=NULL,
                          plotG = TRUE, maxG = 2500, fsuffix = RandString(), image.format='png'){
    if (ref.iter==0 && !is.null(SamplingSize) && ncol(DM) > SamplingSize)
        warning("only ",SamplingSize," of ",ncol(DM)," cells selected for clustering")
    
    ptm=proc.time() #Start clock
    
    params <-as.list(environment())
    params$plotG=FALSE 
    
    #Make sure DM comes with rownames/colnames:
    if (is.null(rownames(DM))) {rownames(DM)=c(1:nrow(DM)); params$DM=DM} 
    if (is.null(colnames(DM))) {colnames(DM)=c(1:ncol(DM)); params$DM=DM} 
    
    
    if (ncol(DM)<300) {use.par=FALSE; params$use.par=FALSE} 
    
    
    if (is.null(SamplingSize)){
    params$SamplingSize=2000    
    }
    
    #######Define functions if use.par=FALSE
    SPearsonCor <- sparse.cor
    PPearsonCor <- stats::cor
    PSpearmanCor <- PSpcor
    PHellinger <- PHellingerMat
    PCanberra <- PCanberraMat
    ShrinkCor <- corpcor::cor.shrink
    
    if (length(ClassAssignment) != ncol(DM))
        stop ("length(ClassAssignment) must be equal to ncol(DM)")
    if(!is.null(BatchAssignment) && length(BatchAssignment) != ncol(DM))
        stop ("length(BatchAssignment) must be equal to ncol(DM)")
    
    # Register cluster here, remove registration block from SC_cluster
    if (isTRUE(use.par)) {
        #######Switch to parallelized functions if use.par=TRUE
        SPearsonCor <- FlashSPearsonCor
        PPearsonCor <- if (checkOpenMP()) FlashPPearsonCorOMP else FlashPPearsonCor
        PSpearmanCor <- if (checkOpenMP()) FlashPSpearmanCorOMP else FlashPSpearmanCor
        PHellinger <- if (checkOpenMP()) FlashPHellingerOMP else FlashPHellinger
        PCanberra <- if (checkOpenMP()) FlashPCanberraOMP else FlashPCanberra 
        ShrinkCor <- FlashShrinkCor
        
        if(ncores=="all"){
            ncores = parallel::detectCores()
            ncores=min(48,floor(0.9*ncores),ceiling(ncol(DM)/200))
        } else{
            ncores=min(ncores,parallel::detectCores(),ceiling(ncol(DM)/200))
        }
        cl<-parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
        
        if (is.null(SamplingSize)){
        params$SamplingSize=max(params$SamplingSize, 250 * ncores )    
        }
    }
    
    ### wrap code in tryCatch block, ensuring that stopCluster(cl) is called even when a condition is raised  
    tryCatch({    
        for (i in 0:ref.iter) { 
            if (i==0) {

                if (ref.iter==0){
                params$ncom=ncom    
                }
                
                Gcounts=colSums(DM>0)
                LowQual=which(Gcounts <  quantile(Gcounts,0.01)  )
                
                if ( (ncol(DM) -length(LowQual) )  > params$SamplingSize ){
                SMPL=sample(1:ncol(DM)[-LowQual],params$SamplingSize)
                }
                
                else{
                SMPL=c(1:ncol(DM))[-LowQual] 
                }
                
                
                message("Preprocessing...", appendLF = FALSE)
                
                ########## Remove ghost cells (cells without detected genes):
                NoData <- which(colSums( DM[,SMPL] ) == 0)
                if (length(NoData > 0)){
                SMPL=SMPL[-NoData]
                }

                params$DM=DM[,SMPL]
                
                ########## Remove no-show genes:
                AllZeroRows=which  ( rowSums( params$DM )<1e-9 )
                if(length(AllZeroRows)>0){
                params$DM=params$DM[-AllZeroRows , ] 
                }
                ##########  Remove invariant genes:
                meanDM=mean( params$DM )
                nSD=apply(params$DM,1,function(x) sd(x)/meanDM)
                ConstRows=which   ( nSD < 0.25 )
                if(length(ConstRows)>0){
                params$DM=params$DM[-ConstRows , ]
                }
                message("\nRemoved ", length(c(ConstRows,AllZeroRows)), " uninformative (invariant/no-show) gene(s)...\n", appendLF = FALSE)

                
                ##########  Remove promiscuous cells (this only affects the sampling iteration):
                DMS=as(params$DM,"dgCMatrix") #filteredGenes x SMPL
                DMS@x=log2(DMS@x+1)
                cM=SPearsonCor(DMS)
                sum.cM=(colSums(cM)-1)/2
                Y1=sum.cM
                X1=log2(Gcounts[SMPL])
                m=nls(Y1 ~ a*X1+b, start=list(a=-5,b=-10)  )
                Yhat=predict(m)
                exclude=which(Y1/Yhat >  quantile(Y1/Yhat,0.5)  & sum.cM >   quantile(sum.cM,0.25 )     )
                fraction=min( (ncol(DM)^2)/1e07,0.9)
                exclude=sample(exclude, ceiling(length(exclude)*fraction) )
                SMPL=SMPL[-c(exclude)]
                cat("::\n",length(SMPL), "\n")
                params$DM=params$DM[,-c(exclude)]

                
                
                #############CV=f(mean) -based filtering:
                CellCounts=colSums(  params$DM )
                nDM=sweep( params$DM ,2,CellCounts,FUN="/")
                #nDM=params$DM*10000 #Counts per 10K
                nDM=nDM*1e9 #Counts per Billion
                
                if(is.null(filter)){
                    medianComplexity=median(apply( params$DM ,2,function(x) sum(x>0))) 
                    filter=ifelse( medianComplexity > 2500,TRUE,FALSE)
                    message("Median Library Complexity: ",medianComplexity," --> Gene Filtering: ", filter ,"\r")
                    
                }
                if (filter){
                    message("\nFiltering Genes...", appendLF = FALSE)
                    X1=log2(rowMeans(nDM+1/ncol(nDM)))
                    Y1=apply(nDM,1,function(x) log2(sd(x)/mean(x+1/length(x))+1/length(x) )  )
                    m=nls(Y1 ~ a*X1+b, start=list(a=-5,b=-10)  )
                    Yhat=predict(m)
                    params$DM=params$DM[which(Y1 > 0.9*Yhat),]
                    message("Removed an additional ",nrow(nDM)-nrow(params$DM), " gene(s) with low variance\n", appendLF = FALSE)
                }
                
                genelist<-rownames(params$DM)
                
                params$ClassAssignment=ClassAssignment[SMPL]
                if (!is.null(BatchAssignment)){
                    params$BatchAssignment=BatchAssignment[SMPL]   
                }
                
                message("done")
                
                
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
                min.csize <-max(4, ceiling(0.25*sqrt(length(memb)) ) )
                nclust=length(unique(memb) )
                good.clust=as.vector(which(table(memb)>=min.csize) )

                if (length(good.clust)<2){
                    message("\nNotice: No substantial clusters found. This might indicate unstructured data...\n", appendLF = FALSE)
                    #break  
                }
                else{
                    message("\n", length(good.clust)," substantial clusters found...\n", appendLF = FALSE)
                }
                
                message("\nBootstrapping to refine clusters...\n", appendLF = FALSE)
                
                #Number of boostrapping samples:
                Nboot.Smpls=min(ceiling(150/(length(good.clust)^2) ), 80)
                Nboot.Smpls=max(Nboot.Smpls,2)
                bootS.size=Nboot.Smpls^(-0.5)
                
                FakeBulk=matrix(0,length(genelist),length(good.clust)*Nboot.Smpls)
                r=0
                for (c in 1:length(good.clust)) {
                    clust=good.clust[c]
                    for (b in 1:Nboot.Smpls){
                        r=r+1
                        cluster.sample=sample(which(memb==clust),ceiling(sum(memb==clust)*bootS.size)+1 ,replace=TRUE    )
                        FakeBulk[,r]=rowSums(DM[genelist,names(memb)][,cluster.sample])
                    }
                }
                

                ###### Calculate distances of all the cells to the FakeBulks:
                message("Calculating Cell Distances to Cluster Centroids (Bulks)...", appendLF = FALSE)
                params$DM <- WScorFB(DM[genelist,], FakeBulk,PSpearmanCor=PSpearmanCor, PPearsonCor=PPearsonCor, PHellinger=PHellinger, PCanberra=PCanberra, ShrinkCor=ShrinkCor)
                message("done")
                cluster.res <- do.call(SC_cluster, c(params,list(comm.method=igraph::cluster_louvain,do.glasso=FALSE,pr.iter=0) ) )
                cluster.res$GeneList<-genelist        
            }
            gc() #Call garbage collector
        }
        
        ######Top FeatureGenes:
        
        
        ######Mark Doublets:
        #if (markDoublets==TRUE){
        #memb=cluster.res$MEMB
        #n=length(memb)
        #nsq=n^2
        #TBL=table(memb)
        #nclust=length(TBL)
        #between=sapply(1:nclust, function(x)  rowMeans(cluster.res$DISTM[,memb==x])    )
        #expected=matrix(0,nclust,nclust)
        #    for (c1 in 1: nclust){
        #        for (c2 in 1: nclust){
        #        expected[c1,c2]=(1+sum(between[,c1] * between[,c2] > 0)   ) 
        #        }
        #    }
        
        #    get.mixr <- function (x) {
        #    o=order(x,decreasing=TRUE)
        #    mixr=x[o[2]]/(x[o[1]]+1e-3)
        #    mixr=1/ ((expected[ o[1],o[2] ]) )
        #    }
        #mixing.ratio=apply(between,1,get.mixr )
        #}
        
        
        
        if (plotG==TRUE){    
            cluster.res[["plotGRAO"]] <- plotGraph(cluster.res, maxG = maxG, fsuffix = fsuffix,image.format = image.format, quiet = FALSE)
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
    return(cluster.res)    
}

