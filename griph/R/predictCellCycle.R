#' Predict (infer) cell cycle stage.
#' 
#' @description Predict (infert) cell cycle stage for individual cells using
#'     sets of known cell cycle-regulated genes.
#' 
#' @param cnt gene-by-cell count matrix (raw or normalized counts). Rownames must
#'     be Entrez Gene identifiers.
#' @param org Character scalar specifying the species, which selects a corresponding
#'     set of known cell cycle associated genes; one of "human.Whitfield" (default, 
#'     use genes from Whitfield et al., 2002) or "mouse".
#' @param cor_thr Numeric scalar value in [0,1] used to filter known cell cycle
#'     genes (see Details). Partial matching is accepted if sufficient for unique
#'     match (e.g. \code{"human"} selects "human.Whitfield").
#' @param granularity Character scalar selecting the level of detail for cell
#'     cycle classification. \code{"high"} (default) will call as many more different
#'     cell cycle stages, depending on the available training genes; \code{"low"}
#'     will call fewer different stages (prefered for data with many non-observed genes). 
#' @param refine_iter Numberic scalar giving the maximum number of refinement iterations.
#' 
#' @details \code{predictCellCycle} uses genes known to peak in transcription at
#'     given cell cycle stages (for \code{org="human.Whitfield"} obtained from
#'     Whitfield et al., 2002; for mouse, \code{"mouse.Whitfield"}
#'     uses 1-to-1 orthologs of the \code{human.Whitfield} genes, as defined in
#'     Homologene release 68, \code{"mouse.Ishida"} uses the genes from Ishida et
#'     al., 2001; see references).
#'     For each cell, five normalized cell cycle stage scores are calculated, and
#'     the cell is assigned to the one that most closely resembles the expected
#'     profile based on correlation. Optionally, labels are further refined
#'     (\code{refine_iter} parameter) by iteratively estimating new cell cycle
#'     score profiles based on estimated cell labels, and re-assigning cells to
#'     the profile with the highest correlation.
#' 
#' @references Whitfield et al., "Identification of genes periodically expressed
#'     in the human cell cycle and their expression in tumors.", Mol Bio Cell, 2002,
#'     \url{http://www.ncbi.nlm.nih.gov/pubmed/12058064}, gene sets WHITFIELD_CELL_CYCLE_G1_S,
#'     WHITFIELD_CELL_CYCLE_G2, WHITFIELD_CELL_CYCLE_G2_M, WHITFIELD_CELL_CYCLE_M_G1 and
#'     WHITFIELD_CELL_CYCLE_S downloaded e.g. from \url{http://software.broadinstitute.org/gsea/msigdb/geneset_page.jsp?geneSetName=WHITFIELD_CELL_CYCLE_G1_S}
#'     
#'     Ishida et al., "Role for E2F in control of both DNA replication and mitotic
#'     functions as revealed from DNA microarray analysis.", Mol Cell Biol, 2001,
#'     \url{https://www.ncbi.nlm.nih.gov/pubmed/11416145}
#'
#' @return A factor of length \code{ncol(cnt)} (the number of cells) with
#'     predicted cell cycle stages.
predictCellCycle <- function(cnt, org = c("human.Whitfield", "mouse.Whitfield", "mouse.Ishida"),
                             cor_thr = 0.2, granularity = c("high", "low"),
                             refine_iter = 0) {
    
    # load gene sets
    org <- match.arg(org)
    gset <- switch(org,
                   human.Whitfield = readRDS(system.file("extdata", "CellCycleGenesets_Human_Whitfield2002.rds", package = "griph")),
                   mouse.Whitfield = readRDS(system.file("extdata", "CellCycleGenesets_Mouse_Whitfield2002.rds", package = "griph")),
                   mouse.Ishida    = readRDS(system.file("extdata", "CellCycleGenesets_Mouse_Ishida2001.rds",    package = "griph")))
    
    # Clean and Nomalize count matrix
    cnt <- t(t(cnt) /colSums(cnt) *1e6)
    cnt <- log2(cnt+1)
    cnt <- cnt[rowSums(cnt>0) > 0,] # Remove genes not measured in any cell

    # check if rownames(cnt) match to gene sets
    fractionMatched <- mean(unlist(gset) %in% rownames(cnt))
    if(fractionMatched < 0.1)
        stop(if(fractionMatched == 0) "None " else paste0("Only ", round(fractionMatched*100,1), "% "),
             "of the known cell cycle genes are found in the count matrix. ",
             "Did you use ",sub("[.].+$","",org)," Entrez Gene identifiers as row names?")
    
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
    if (org %in% c("human.Whitfield","mouse.Whitfield") && granularity == "high") {
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
    } else if (org %in% c("human.Whitfield","mouse.Whitfield") && granularity == "low") {
        phasePatterns <- cbind(pat01=c(G1.S=1, S=0, G2=0, G2.M=0, M.G1=0), # only G1/S
                               pat02=c(G1.S=0, S=1, G2=0, G2.M=0, M.G1=0), # only S
                               pat03=c(G1.S=0, S=0, G2=1, G2.M=0, M.G1=0), # only G2
                               pat04=c(G1.S=0, S=0, G2=0, G2.M=1, M.G1=0), # only G2/M
                               pat05=c(G1.S=0, S=0, G2=0, G2.M=0, M.G1=1)) # only M/G1
        pat2cc <- c(pat01="G1.S", pat02="S", pat03="G2", pat04="G2.M", pat05="M.G1")
    } else if (org == "mouse.Ishida" && granularity == "high") {
        phasePatterns <- cbind(pat01=c(G0=1, G1_Early=0, G1_Cycle=0, G1_Growth=0, G1.S_Cycle=0, G1.S_Growth=0, G2_Cycle=0), # only G0
                               pat02=c(G0=0, G1_Early=1, G1_Cycle=0, G1_Growth=0, G1.S_Cycle=0, G1.S_Growth=0, G2_Cycle=0), # only G1_Early
                               pat03=c(G0=0, G1_Early=0, G1_Cycle=1, G1_Growth=0, G1.S_Cycle=0, G1.S_Growth=0, G2_Cycle=0), # only G1_Cycle
                               pat04=c(G0=0, G1_Early=0, G1_Cycle=0, G1_Growth=1, G1.S_Cycle=0, G1.S_Growth=0, G2_Cycle=0), # only G1_Growth
                               pat05=c(G0=0, G1_Early=0, G1_Cycle=0, G1_Growth=0, G1.S_Cycle=1, G1.S_Growth=0, G2_Cycle=0), # only G1.S_Cycle
                               pat06=c(G0=0, G1_Early=0, G1_Cycle=0, G1_Growth=0, G1.S_Cycle=0, G1.S_Growth=1, G2_Cycle=0), # only G1.S_Growth
                               pat07=c(G0=0, G1_Early=0, G1_Cycle=0, G1_Growth=0, G1.S_Cycle=0, G1.S_Growth=0, G2_Cycle=1)) # only G2_Cycle
        pat2cc <- c(pat01="G0", pat02="G1_Early", pat03="G1_Cycle", pat04="G1_Growth",
                    pat05="G1.S_Cycle", pat06="G1.S_Growth", pat07="G2_Cycle")
        
    } else if (org == "mouse.Ishida" && granularity == "low") {
        phasePatterns <- cbind(pat01=c(G0=1, G1_Early=0, G1_Cycle=0, G1_Growth=0, G1.S_Cycle=0, G1.S_Growth=0, G2_Cycle=0), # only G0
                               pat02=c(G0=0, G1_Early=1, G1_Cycle=1, G1_Growth=1, G1.S_Cycle=0, G1.S_Growth=0, G2_Cycle=0), # G1_Early, G1_Cycle, G1_Growth
                               pat03=c(G0=0, G1_Early=0, G1_Cycle=0, G1_Growth=0, G1.S_Cycle=1, G1.S_Growth=1, G2_Cycle=0), # G1.S_Cycle, G1.S_Growth
                               pat04=c(G0=0, G1_Early=0, G1_Cycle=0, G1_Growth=0, G1.S_Cycle=0, G1.S_Growth=0, G2_Cycle=1)) # only G2_Cycle
        pat2cc <- c(pat01="G0", pat02="G1", pat03="G1.S", pat04="G2")
    } else {
        stop("don't know how to get expected phase patterns")
    }
        

    # Assigning cells to different cell cycle patterns
    patCor <- cor(phaseSpecificScoreNorm2, phasePatterns)
    patCells <- colnames(phasePatterns)[unlist(apply(patCor,1,which.max))]
    
    # optional iterative refinement of cell labels
    refine_rate <- 0.5
    while (refine_iter > 0) {
        #message("refining ",refine_iter)
        #message("\t: ",paste(table(patCells),collapse = ", "))
        patnms <- colnames(phasePatterns)
        phasePatterns <- (1 - refine_rate) * phasePatterns +
            refine_rate * do.call(cbind,
                                  lapply(patnms,
                                         function(patnm) {
                                             rowMeans(phaseSpecificScoreNorm2[, patCells == patnm, drop = FALSE])
                                         })
                                  )
        patCor <- cor(phaseSpecificScoreNorm2, phasePatterns)
        patCellsNew <- patnms[unlist(apply(patCor,1,which.max))]
        #message("\t: ",paste(table(patCellsNew),collapse = ", "))
        if (all(patCells == patCellsNew)) {
            break
        } else {
            refine_iter <- refine_iter - 1
            patCells <- patCellsNew
        }
    }
  
    return(factor(unname(pat2cc[patCells]), levels = unname(pat2cc)))
}

