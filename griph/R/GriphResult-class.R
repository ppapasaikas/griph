setOldClass("igraph") # make S3 class igraph visible to S4

#' An S4 class to represent the result of identified cell types
#' obtained by running \code{griph}.
#'
#' @slot DM A g-by-c matrix with raw counts for g genes and c cells, or a c-by-c
#'     matrix with correlation coefficients between pairs of cell
#' @slot is.cor A length-one logical. \code{TRUE} indicates that \code{DM} is a correlation matrix.
#' @slot ClassAssignment A numeric vector of length c with known cell types for
#'     the provided cells (columns of \code{M} matrix)
#' @slot BatchAssignment A numeric vector of length c with known batch assignments
#'      for the provided cells (columns of \code{M} matrix)
#' @slot MEMB A numeric vector of length c giving the cell type memberships
#'     identified by \code{griph} for the provided cells (columns of \code{M} matrix)
#' @slot DISTM A c-by-c matrix with similarities for each pair of cells
#' @slot ConfMatrix A n-by-m confusion matrix with the numbers of cells of known
#'     class 1..n (\code{ClassAssignment}, rows) classified by \code{griph} as cell type
#'     1..m (columns).
#' @slot miscl A length-one numerical value giving the misclassification error.
#' @slot GRAO An object of class \code{igraph} with c nodes representing cells
#'     and edges between similar cells.
#' @slot plotGRAO An object of class \code{igraph} with c nodes representing cells
#'     and edges between similar cells. This is a rendering optimized version of GRAO.
#'     Edges are pruned, vertex attributes for plotting are added and vertices are sampled
#'     if they exceed the maxG argument.  
GriphResult <- setClass("GriphResult",
                        slots = representation(DM = "matrix",
                                               is.cor = "logical",
                                               ClassAssignment = "factor",
                                               BatchAssignment = "factor",
                                               MEMB = "numeric",
                                               DISTM = "matrix",
                                               ConfMatrix = "table",
                                               miscl = "numeric",
                                               GRAO = "igraph",
                                               plotGRAO="igraph")
)

#' @param object Instance of class \code{GriphResult}.
#' @describeIn GriphResult Print a summary of a GriphResult object.
setMethod("show", "GriphResult", function(object) {
    DMtype <- if(object@is.cor) "correlation coefficients" else "gene counts"
    nClassKnown <- nlevels(object@ClassAssignment)
    nClassPred <- nlevels(object@MEMB)
    nBatch <- nlevels(object@BatchAssignment)
    objectname <- deparse(substitute(object)) # does not work for automatic printing...
    cat("GriphResult\n", sep="")
    cat("Options   : DM                   : ", nrow(object@DM), "-by-", ncol(object@DM), " matrix of ", DMtype,
      "\n            known cell types     : ", if (nClassKnown > 1) sprintf("yes (%d)",nClassKnown) else "no",
      "\n            known batches        : ", if (!is.null(object@BatchAssignment)) sprintf("yes (%d)",nBatch) else "no",
      "\nResults   : identified cell types: ", nClassPred,
      "\n            missclassification   : ", if (nClassKnown > 1) object@miscl else "not available", "\n\n", sep="")
    cat("Slots     : use slot(",objectname,", 'nm') to get the content of slot 'nm'",
      "\n            use slotNames(",objectname,") to get all available slot names", "\n\n", sep="")
})
