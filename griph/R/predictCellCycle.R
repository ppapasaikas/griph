#' Predict (infer) cell cycle stage.
#' 
#' @description Predict (infert) cell cycle stage for individual cells using
#'     sets of known cell cycle-regulated genes.
#' 
#' @param cnt gene-by-cell count matrix (raw or normalized counts). Rownames must
#'     be Entrez Gene identifiers.
#' @param org Character scalar specifying the species; one of "human"
#'     (default) or "mouse".
#' @param cor_thr Numeric scalar value in [0,1] used to filter known cell cycle
#'     genes (see Details).
#' @param granularity Character scaler selecting the level of detail for cell
#'     cycle classification. \code{"high"} (default) will call 10 different cell
#'     cycle stages, \code{"low"} will call only 5 different stages. 
#' 
#' @details \code{predictCellCycle} uses genes known to peak in transcription at
#'     given cell cycle stages (for human obtained from Whitfield et al., 2002, dowloaded
#'     from broadinstitute.org, see references; for mouse 1-to-1 orthologs of the
#'     human genes, as defined in Homologene release 68 were used).
#'     For each cell, five normalized cell cycle stage scores are calculated, and
#'     the cell is assigned to the one that most closely resembles the expected
#'     profile based on correlation.
#' 
#' @references Whitfield et al., "Identification of genes periodically expressed
#'     in the human cell cycle and their expression in tumors.", Mol Bio Cell, 2002,
#'     \url{http://www.ncbi.nlm.nih.gov/pubmed/12058064},
#'     G1_S genes: \url{http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=WHITFIELD_CELL_CYCLE_G1_S},
#'     G2 genes: \url{http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=WHITFIELD_CELL_CYCLE_G2},
#'     G2_M genes: \url{http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=WHITFIELD_CELL_CYCLE_G2_M},
#'     M_G1 genes: \url{http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=WHITFIELD_CELL_CYCLE_M_G1},
#'     S genes: \url{http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=WHITFIELD_CELL_CYCLE_S}
#'
#' @return A factor of length \code{ncol(cnt)} (the number of cells) with
#'     predicted cell cycle stages.
predictCellCycle <- function(cnt, org = c("human", "mouse"), cor_thr = 0.3,
                             granularity = c("high", "low")) {
    
    # load gene sets
    org <- match.arg(org)
    gset <- switch(org,
                   human = readRDS(system.file("extdata", "CellCycleGenesets_Human.rds", package = "griph")),
                   mouse = readRDS(system.file("extdata", "CellCycleGenesets_Mouse.rds", package = "griph")))
    
    # Clean and Nomalize count matrix
    cnt <- t(t(cnt) /colSums(cnt) *1e6)
    cnt <- log2(cnt+1)
    cnt <- cnt[rowSums(cnt>0) > 0,] # Remove genes not measured in any cell

    # check if rownames(cnt) match to gene sets
    fractionMatched <- mean(unlist(gset) %in% rownames(cnt))
    if(fractionMatched < 0.1)
        stop(if(fractionMatched == 0) "None " else paste0("Only ", round(fractionMatched*100,1), "% "),
             "of the known cell cycle genes are found in the count matrix.",
             "Did you use ",org," Entrez Gene identifiers as row names?")
    
    # Clean Genesets
    gset <- lapply(gset, function(x) x[x %in% rownames(cnt)])
  
    # Calculate cell cycle phase-specific scores:
    phaseSpecificScore <- matrix(NA, nrow=length(gset), ncol=ncol(cnt), dimnames=list(names(gset),colnames(cnt)))
    for(i in 1:nrow(phaseSpecificScore))
      phaseSpecificScore[i,] <- colMeans(cnt[gset[[i]],],na.rm=FALSE)

    # Remove all genes whose correation (expression) with their parent gene set is <= cor_thr
    for(i in 1:nrow(phaseSpecificScore)) {
        x <- cor(phaseSpecificScore[i,], t(cnt[gset[[i]],]))
        gset[[i]] <- gset[[i]][x >= cor_thr]
    }

    # RE-Calculate cell cycle phase-specific scores:
    for(i in 1:nrow(phaseSpecificScore))
        phaseSpecificScore[i,] <- colMeans(cnt[gset[[i]],],na.rm=FALSE)
  
    phaseSpecificScoreNorm1 <- t(scale(t(phaseSpecificScore))) # 1. scale phases
    phaseSpecificScoreNorm2 <- scale(phaseSpecificScoreNorm1)  # 2. scale cells

    # Define cell cycle patterns and classify the cells based on the maximal correlation of the 
    # phase-specific scores with these patterns.
    granularity <- match.arg(granularity)
    if (granularity == "high") {
        phasePatterns <- cbind(pat01=c(G1.S=1, S=0, G2=0, G2.M=0, M.G1=0), # only G1/S
                               pat02=c(G1.S=1, S=1, G2=0, G2.M=0, M.G1=0), # both G1/S and S
                               pat03=c(G1.S=0, S=1, G2=0, G2.M=0, M.G1=0), # only S
                               pat04=c(G1.S=0, S=1, G2=1, G2.M=0, M.G1=0), # both S and G2
                               pat05=c(G1.S=0, S=0, G2=1, G2.M=0, M.G1=0), # only G2
                               pat06=c(G1.S=0, S=0, G2=1, G2.M=1, M.G1=0), # G2 and G2/M
                               pat07=c(G1.S=0, S=0, G2=0, G2.M=1, M.G1=0), # only G2/M
                               pat08=c(G1.S=0, S=0, G2=0, G2.M=1, M.G1=1), # G2/M and M/G1
                               pat09=c(G1.S=0, S=0, G2=0, G2.M=0, M.G1=1), # only M/G1
                               pat10=c(G1.S=1, S=0, G2=0, G2.M=0, M.G1=1)) # both M/G1 and G1/S
      pat2cc <- c(pat01="G1.S", pat02="G1.S/S", pat03="S", pat04="S/G2", pat05="G2", 
                  pat06="G2/G2.M", pat07="G2.M", pat08="G2.M/M.G1", pat09="M.G1", pat10="M.G1/G1.S")
    } else {
        phasePatterns <- cbind(pat01=c(G1.S=1, S=0, G2=0, G2.M=0, M.G1=0), # only G1/S
                               pat02=c(G1.S=0, S=1, G2=0, G2.M=0, M.G1=0), # only S
                               pat03=c(G1.S=0, S=0, G2=1, G2.M=0, M.G1=0), # only G2
                               pat04=c(G1.S=0, S=0, G2=0, G2.M=1, M.G1=0), # only G2/M
                               pat05=c(G1.S=0, S=0, G2=0, G2.M=0, M.G1=1)) # only M/G1
        pat2cc <- c(pat01="G1.S", pat02="S", pat03="G2", pat04="G2.M", pat05="M.G1")
    }

    # Assigning cells to different cell cycle patterns
    patCor <- cor(phaseSpecificScoreNorm2, phasePatterns)
    patCells <- colnames(phasePatterns)[unlist(apply(patCor,1,which.max))]
  
    return(factor(unname(pat2cc[patCells]), levels = unname(pat2cc)))
}

