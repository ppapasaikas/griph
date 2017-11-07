#' @title Pseudo time ordering of cell groups
#' 
#' @description Calculates a pseudotime ordering of groups of cells using the affinity matrix  
#'     and membership assignemnt returned from \code{\link{griph_cluster}} or the affinity matrix 
#'     and a custom membership assignment.
#' 
#' @param gr A \code{griph} result, as returned by \code{\link{griph_cluster}}.
#' @param memb A membership vector of length equal to the number of cells present in gr.


griph_order <- function(gr, memb=NULL )
{
    
    if (is.null(memb) ){
    memb <- gr$MEMB
    
    }else {
        if (length(memb) != ncol(gr$DISTM) ) {
        stop("length(memb) must be equal to ncol(gr$DISTM)")   
        }
    }
    
    min.csize <- max(4, ceiling(0.25*sqrt(length(memb))))
    good.clust <- as.vector(which(table(memb) >= min.csize))
    
    if (length(good.clust) < 2) {
        stop("\nToo few (<2) substantial clusters present. No ordering done...\n", appendLF = FALSE)
        #break  
    } else {
        message("\n", length(good.clust)," substantial clusters present...\n", appendLF = FALSE)
    }
    
    use.nodes=memb %in% names(table(memb))[good.clust]
    new.memb=memb[use.nodes]
    GRAOi=igraph::induced_subgraph(gr$GRAO,vids = use.nodes)
    
    class.collapse <- factor(new.memb[use.nodes], levels = sort(unique( new.memb[use.nodes] ))  )
    CollapsedDISTM=matrix(0,length(good.clust),length(good.clust),dimnames = list( levels( class.collapse ), levels( class.collapse)  )  )

    DISTM=1-Matrix::as.matrix( gr$DISTM[use.nodes,use.nodes] )
    diag(DISTM)=0
    
    for (i in 1:length(good.clust)){
        for (j in 1:length(good.clust)){ 
            if (i>=j) {next}
            cl1=levels( class.collapse )  [i]
            cl2=levels( class.collapse )  [j]
            CollapsedDISTM[i,j]= mean(DISTM[ new.memb==cl1,new.memb==cl2 ], trim=0.005 ) 
            CollapsedDISTM[j,i]=   CollapsedDISTM[i,j]
        }
    }
    
    CollapsedDISTM=CollapsedDISTM-min(CollapsedDISTM[upper.tri(CollapsedDISTM)])
    diag(CollapsedDISTM)=0
    CollapsedDISTM=as.matrix(dist(CollapsedDISTM))
    GRAOp <- igraph::graph.adjacency(CollapsedDISTM, mode = c("max"), weighted = TRUE, diag = FALSE)
    MST=igraph::mst(GRAOp)
    
gr$MST=MST
gr$MST.memb=memb
return(gr)
}



