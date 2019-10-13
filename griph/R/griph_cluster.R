#' Combined distance calculations.
#' 
#' @description Calculate and combine several between-cell distance measures to
#'     obtain robust cell-to-cell distances.
#' 
#' @param M gene-by-cell count matrix.
#' @param K Number of nearest neighbors to retain
#' @param FB gene-by-FakeBulks gene sums count matrix.
#' @param PPearsonCor Function to calculate Pearson correlation between the columns of two matrices.
#' @param PSpearmanCor Function to calculate Spearman correlation between the columns of two matrices.
#' @param PHellinger Function to calculate Hellinger Distance between the columns of two matrices.
#' @param PCanberra Function to calculate Canberra Distance between the columns of two matrices.
#' @param ShrinkCor Function to calculate shrinkage correlation.
#' @return cell-by-cell distance matrix.
WScorFB <- function(M, FB, K=50, PSpearmanCor, PPearsonCor, PHellinger, PCanberra, ShrinkCor=ShrinkCor) {
    
    CellIds <- colnames(M)
    dimnames(M) <- NULL
    
    ptm1 <- proc.time() #Start clock
    
    D <- PPearsonCor(log2(FB + 1), log2(M + 1))
    R <- vapply(c(1:ncol(D)), function(x) rank(D[,x]), FUN.VALUE = double(length = nrow(D)))  # Pearson's cor    
    
    Te1 <- signif((proc.time() - ptm1)[3], digits = 6)
    message("\n..PearsonC: (Elapsed Time: ", Te1, ")")
    
    ######## Counts per 100K:
    CellCounts <- colSums(FB)
    FB <- sweep(FB, 2, CellCounts, FUN = "/")
    FB <- FB * 100000
    CellCounts <- colSums(M)
    M <- sweep(M, 2, CellCounts, FUN = "/")
    M <- M * 100000 
    
    ptm1 <- proc.time() #Start clock
    Dt <- PCanberra(log2(FB + 1), log2(M + 1))   
    Dt <- 1 - ((Dt - min(Dt)) / diff(range(Dt)))
    R <- R + vapply(c(1:ncol(Dt)), function(x) rank(Dt[,x]), FUN.VALUE = double(length = nrow(Dt)))  # Canberra
    Te1 <- signif((proc.time() - ptm1)[3], digits = 6)
    message("..Canberra: (Elapsed Time: ", Te1, ")")
    
    ptm1 <- proc.time() #Start clock
    Dt <- PSpearmanCor(FB, M)
    R <- R + vapply(c(1:ncol(Dt)), function(x) rank(Dt[,x]), FUN.VALUE = double(length = nrow(Dt)))  # Spearman's cor
    Te1 <- signif((proc.time() - ptm1)[3], digits = 6)
    message("..SpearmanC: (Elapsed Time: ", Te1, ")")
    
    ptm1 <- proc.time() #Start clock
    Dt <- PHellinger(FB, M)
    Dt <- 1 - ((Dt - min(Dt)) / diff(range(Dt)))
    R <- R + vapply(c(1:ncol(Dt)), function(x) rank(Dt[,x]), FUN.VALUE = double(length = nrow(Dt)))  # Hellinger distance
    Te1 <- signif((proc.time() - ptm1)[3], digits = 6)
    message("..Hellinger: (Elapsed Time: ", Te1, ")")
    
    ptm1 <- proc.time() #Start clock
    R <- (R / 4)^2
    #R <- sweep(R, 2, colMeans(R), "-")
    #K=min(max(floor(0.25 * ( sqrt(ncol(R)) + nrow(R) )) , 10), floor(ncol(R)) / 1.5) 
    R <- buildEdgeMatrix(R, distance_method = "Cosine", K = K, threads=foreach::getDoParWorkers()    )    #
    R <- Matrix::sparseMatrix(i = R$i, j = R$j, x = 1 - (R$x / 2), dims = attr(R,"dims"), dimnames = list(CellIds,CellIds))
    Te1 <- signif((proc.time() - ptm1)[3], digits = 6)
    message("Calculating KNNGraph (Elapsed Time: ", Te1, ")")
    
    return(R)    #This is NOT a symmetric matrix
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
#' @param K Number of nearest neighbors to consider. If \code{NULL} (default) it is set to to \code{0.25*sqrt(k)+4*D}, D being
#'     the dataset dimensionality determined in the initialization (sampling) step. Smaller K values result in more disconnected
#'     cell clusters and vice versa. Typically should be in the range [20,200]. 
#' @param SamplingSize Number of sampled cells in initialization step. If NULL it is set to \code{min(2000,k/2)}
#' @param ref.iter Number of clustering refinement iterations.  If set to 0 only the clustering initialization 
#'     step is performed to \code{min(SamplingSize,ncol(DM))} cells. 
#' @param use.par If \code{TRUE}, use parallel versions of distance calculation
#'     functions based on \code{\link[foreach]{foreach}} (see details).
#' @param ncores a numeric(1) or character(1), either specifying the number of
#'     parallel jobs to use, or \code{"all"} (default) to use up to 90 percent of
#'     the available cores. Ignored if \code{use.par==FALSE}.
#' @param filter TRUE/FALSE or a positive numeric values. Determines  whether filtering for retaining only overdispersed genes should be performed.
#'     If FALSE no gene filtering is performed. If TRUE (default) only the top 25% overdispersed genes are retained.
#'     If set to a number between 0 and 1 than this fraction of genes is retained.
#'     If set to an integer > 1 than this number of genes is retained
#### Filter genes according to cv=f( mean ) fitting. Default: \code{TRUE}.
#' @param rho Inverse covariance matrix regularization (graph sparsification) parameter -> [0,1].
#'     Default=0.25.  The parameter is then automatically scaled to account for
#'     number of variables or converted into a matrix and adjusted according to
#'     the batch.penalty factor to account for BatchAssignment (if given).
#' @param batch.penalty [0,1] rho scaling factor for enforcing topological constraints
#'     variables according to \code{BatchAssignment}. For penalty p  -> rho_same_batch=rho^(1-p),
#'     rho_diff_batch=rho^(1+p). It is ignored if \code{BatchAssignment==NULL}. 
#' @param seed Set seed for reproducible results.
#' @param ClassAssignment If available a numeric vector of length \code{k} with numeric
#'     class labels (e.g-> c(1,2,1,1,1,2,3,3,3,1,2))
#' @param BatchAssignment If available a numeric vector of length \code{k} with numeric
#'     batch labels (e.g-> c(1,2,1,1,1,2,3,3,3,1,2))
#' @param ncom Forces detection of a fixed number of clusters (communities).
#' @param plot_ if \code{TRUE} plots the obtained graph using \code{plotLVis}. The resulting
#'     2D-embedding is stored in the \code{plotLVis} component of the result.
#' @param maxG Approximate maximal number of vertices to include when plotting the graph.
#' @param fsuffix A suffix added to the file names of output plots. If \code{NULL} (default),
#'     it will use a random 5 character string.
#' @param image.format Specifies the format of the created images. Currently only pdf and png filetypes are supported.
#' 
## @param comm.method  Community detection algorithm. See igraph "communities". By default multilevel louvain is used. 
#' 
#' @return Currently a list with the clustering results:
#' @return MEMB: A vector of membership assignment for each cell
#' @return MEMB.true: A vector with the true cell assignment if one was specified.
#' @return DISTM: A sparse (dgCMatrix) k x k cell distance matrix 
#' @return ConfMatrix: Confusion matrix based on the MEMB and MEMB.true vectors.
#' @return miscl: Misclassification error based on the MEMB and MEMB.true vectors.
#' @return GRAO: An igraph graph object modelling the cell population
#' @return plotLVis: The plotLVis projection of the graph if plot_ was set to TRUE. NULL otherwise
#' @return complete_MEMB: The igraph communities object returned by the community detection algorithm 
#' @return GeneList: A character vector with the genes that survived the filtering.

griph_cluster <- function(DM, K=NULL, SamplingSize= NULL, ref.iter = 3, use.par = TRUE, ncores = "all",
                          filter = TRUE, rho = 0.25, batch.penalty = 0.5, seed = 127350,
                          ClassAssignment = rep(1,ncol(DM)), BatchAssignment = NULL, ncom = NULL,
                          plot_ = TRUE, maxG = 2500, fsuffix = NULL, image.format='png')
                          {
    if (ref.iter == 0 && !is.null(SamplingSize) && ncol(DM) > SamplingSize)
        warning("only ", SamplingSize," of ", ncol(DM)," cells selected for clustering")
    
    ptm <- proc.time() #Start clock
    set.seed(seed = seed) #Set seed for reproducible results
    params <- as.list(environment())
    params$plot_ <- FALSE 
    
    #Make sure DM comes with rownames/colnames:
    if (is.null(rownames(DM))) {
        rownames(DM) <- c(1:nrow(DM))
        params$DM <- DM
    } 
    if (is.null(colnames(DM))) {
        colnames(DM) <- c(1:ncol(DM))
        params$DM <- DM
    } 
    if (ncol(DM) < 300) {
        use.par <- FALSE
        params$use.par <- FALSE
    } 
    if (is.null(SamplingSize)) {
        params$SamplingSize <- max(500,min(2000, ceiling(ncol(DM)/2)) )
    }
    
    #######Define functions if use.par=FALSE
    SPearsonCor <- sparse.cor
    PPearsonCor <- stats::cor
    PSpearmanCor <- PSpcor
    PHellinger <- PHellingerMat
    PCanberra <- PCanberraMat
    ShrinkCor <- corpcor::cor.shrink
    
    if (length(ClassAssignment) != ncol(DM))
        stop("length(ClassAssignment) must be equal to ncol(DM)")
    if (!is.null(BatchAssignment) && length(BatchAssignment) != ncol(DM))
        stop("length(BatchAssignment) must be equal to ncol(DM)")
    
    # Register cluster here, remove registration block from SC_cluster
    if (isTRUE(use.par)) {
        #######Switch to parallelized functions if use.par=TRUE
        SPearsonCor <- FlashSPearsonCor
        PPearsonCor <- if (checkOpenMP()) FlashPPearsonCorOMP else FlashPPearsonCor
        PSpearmanCor <- if (checkOpenMP()) FlashPSpearmanCorOMP else FlashPSpearmanCor
        PHellinger <- if (checkOpenMP()) FlashPHellingerOMP else FlashPHellinger
        PCanberra <- if (checkOpenMP()) FlashPCanberraOMP else FlashPCanberra 
        ShrinkCor <- FlashShrinkCor
        
        if (ncores == "all") {
            ncores <- parallel::detectCores()
            ncores <- min( ceiling(0.9 * ncores), ceiling(ncol(DM) / 200))
        } else {
            ncores <- min(ncores, parallel::detectCores(), ceiling(ncol(DM) / 200))
        }
        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)
        
        #if (is.null(SamplingSize)) {
        #    params$SamplingSize <- max(params$SamplingSize, 250 * ncores )    
        #}
    }
    
    ### wrap code in tryCatch block, ensuring that stopCluster(cl) is called even when a condition is raised  
    tryCatch({   
        current.iter=0
        iter.number=0
        continue=TRUE
        while (continue==TRUE) { 
            if (current.iter == 0) {
                params$ncom <- ncom

                Gcounts <- colSums(DM > 0)
                LowQual <- which(Gcounts <= quantile(Gcounts, 0.01))
                
                if ((ncol(DM) - length(LowQual))  > params$SamplingSize) {
                    SMPL <- sample(1:ncol(DM)[-LowQual], params$SamplingSize)
                } else {
                    SMPL <- c(1:ncol(DM))[-LowQual] 
                }
                
                
                message("Preprocessing...", appendLF = FALSE)
                
                ########## Remove ghost cells (cells without detected genes):
                NoData <- which(colSums(DM[,SMPL]) == 0)
                if (length(NoData > 0)) {
                    SMPL <- SMPL[-NoData]
                }
                params$DM <- DM[, SMPL]
                
                ########## Remove no-show genes:
                AllZeroRows <- which(rowSums(params$DM) < 1e-9)
                if (length(AllZeroRows) > 0) {
                    params$DM <- params$DM[-AllZeroRows, ] 
                }
                ##########  Remove invariant  (completely flat) genes:
                meanDM <- rowMeans(params$DM)
                sdDM <-  apply(params$DM, 1, function(x) sd(x) )
                nSD <- sdDM/mean(meanDM)
                ConstRows <- which(nSD < 1e-2 | sdDM < 0.5)
                if (length(ConstRows) > 0) {
                    params$DM <- params$DM[-ConstRows , ]
                }

                message("\nRemoved ", length(c(ConstRows,AllZeroRows)), " uninformative (completely flat) gene(s)...\n", appendLF = FALSE)

                
                ##########  Remove promiscuous cells (this only affects the sampling iteration):
                DMS <- as(params$DM, "dgCMatrix") #filteredGenes x SMPL
                DMS@x <- log2(DMS@x + 1)
                cM <- SPearsonCor(DMS)
                sum.cM <- (colSums(cM) - 1) / 2
                Y1 <- sum.cM
                X1 <- log2(Gcounts[SMPL])
                #m <- nls(Y1 ~ a * X1 + b, start = list(a = -5, b = -10))
                m <- lm(Y1 ~ X1)
                Yhat <- predict(m)
                exclude <- which(Y1 / Yhat > quantile(Y1 / Yhat, 0.25) & sum.cM > quantile(sum.cM, 0.25))
                fraction <- min( ((ncol(DM)^2) / 1e06), 0.9)
                exclude <- sample(exclude, ceiling(length(exclude) * fraction))
                SMPL <- SMPL[-c(exclude)]
                params$DM <- params$DM[, -c(exclude)]
                
                
                #############Filtering to retain only overdispersed genes: 
                if (filter) {
                    if (!is.numeric(filter)) { keep <- ceiling(0.5*nrow(params$DM)) }
                    else if (filter <  0)    { stop("filter cannot be a negative number") } 
                    else if (filter <= 1 )   { keep <- ceiling(nrow(params$DM) * filter) }
                    else if (filter > 1 )    { keep <- ceiling(filter) }
                    if (keep < nrow(params$DM)) {   
                        message("\nFiltering Genes...", appendLF = FALSE)
                        norm_GeneDispersion <- select_variable_genes(params$DM)
                        disp_cut_off <- sort(norm_GeneDispersion,decreasing = TRUE)[keep]
                        use <- norm_GeneDispersion >= disp_cut_off
                        fraction <- signif((100 * sum(use)) / nrow(params$DM), digits = 3)
                        params$DM <- params$DM[use,]
                        message("\nRetained the top ", fraction, "% overdispersed gene(s)\n", appendLF = FALSE)
                    }
                }
                
                genelist <- rownames(params$DM)
                
                params$ClassAssignment <- ClassAssignment[SMPL]
                if (!is.null(BatchAssignment)) {
                    params$BatchAssignment <- BatchAssignment[SMPL]   
                }
                
                message("...done")
                
                cluster.res <- do.call(SC_cluster, c(params, list(pr.iter = 1, iter.number=iter.number)))
                if (current.iter >= ref.iter) {continue=FALSE}
                current.iter <- current.iter +1
                iter.number <- iter.number +1

            } else {
                
                message("\n\nRefining Cluster Structure...\n", appendLF = FALSE)
                
                params$ncom <- ncom
                
                params$is.cor <- TRUE
                params$ClassAssignment <- ClassAssignment
                params$BatchAssignment <- BatchAssignment  


                ####### construct cell2cell correlation matrix using the current cluster.res: ########
                memb <- cluster.res$MEMB
                min.csize <- max(4, ceiling(0.25*sqrt(length(memb))))
                nclust <- length(unique(memb))
                good.clust <- as.vector(which(table(memb) >= min.csize))

                if (length(good.clust) < 2) {
                    message("\nNotice: No substantial clusters found. This might indicate unstructured data...\n", appendLF = FALSE)
                    #break  
                } else {
                    message("\n", length(good.clust)," substantial clusters found...\n", appendLF = FALSE)
                }
                
                message("\nBootstrapping to refine clusters...\n", appendLF = FALSE)
                
                #Number of boostrapping samples:
                Nboot.Smpls <- min(ceiling(250 / (length(good.clust)^2)), 80)
                Nboot.Smpls <- max(Nboot.Smpls, 16)
                bootS.size <- Nboot.Smpls^(-0.4)
                
                FakeBulk <- matrix(0, length(genelist), length(good.clust) * Nboot.Smpls)
                r <- 0
                for (c in 1:length(good.clust)) {
                    clust <- good.clust[c]
                    ssize <- ceiling(sum(memb == clust) * bootS.size) + 1
                    ssize <- min(100,ssize)
                    for (b in 1:Nboot.Smpls) {
                        r <- r + 1
                        
                        cluster.sample <- sample(which(memb == clust), ssize ,replace = TRUE) 
                        FakeBulk[,r] <- rowMeans(DM[genelist, names(memb)][, cluster.sample])
                    }
                }
                

                #### pass K for mutual NN to  WScorFB, SC_cluster. Should be the same as K passed in buildEdgeMatrix of WScorFB
                if (is.null(K)) {
                Kmnn=min(max(floor(0.25 * ( sqrt(ncol(DM)) + ncol(FakeBulk) )) , 10), floor(ncol(DM)) / 1.5) 
                }
                else{
                Kmnn=K
                }
                

                ###### Calculate distances of all the cells to the FakeBulks:
                message("Calculating Cell Distances to Cluster Centroids (Bulks)...", appendLF = FALSE)
                params$DM <- WScorFB(DM[genelist,], FakeBulk,PSpearmanCor = PSpearmanCor,
                                     PPearsonCor = PPearsonCor, PHellinger = PHellinger,
                                     PCanberra = PCanberra, ShrinkCor = ShrinkCor, K=Kmnn)
                message("done")

                cluster.res <- do.call(SC_cluster, c(params, list(do.glasso = FALSE, pr.iter = 0, Kmnn=Kmnn, iter.number=iter.number) ) )
                cluster.res$GeneList <- genelist   
                
                if (current.iter >= ref.iter) {continue=FALSE}
                current.iter <- current.iter + 1
                iter.number <- iter.number+1
                
            }
            gc() #Call garbage collector
        }
        
        
        ######Top FeatureGenes:
        # PLACEHOLDER
        #######################
        

        ######### graph visualization
        if (plot_ == TRUE) {
            if (is.null(fsuffix))
                fsuffix <- RandString()
            # cluster.res[["plotGRAO"]] <- plotGraph(cluster.res, maxG = maxG, fsuffix = fsuffix,
            #                                        image.format = image.format, quiet = FALSE)
            cluster.res[["plotLVis"]] <- plotLVis(cluster.res, fsuffix = fsuffix, image.format = image.format, quiet = FALSE)
        }
    
        
            
    }, # end of tryCatch expression, cluster object cl not needed anymore    
    
    
    finally = { 
        ##### Stop registered cluster:
        if (isTRUE(use.par) & foreach::getDoParRegistered())
            parallel::stopCluster(cl)
    })
    
    
    ################## Stop the clock #######################
    Te <- (proc.time() - ptm)[3]
    Te <- signif(Te, digits <- 6)
    message("Finished (Elapsed Time: ", Te, ")")
    return(cluster.res)    
}

