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
#' @slot specp A c-by-k matrix with the spectral projection of the c cells to k spectral components
#' @slot ConfMatrix A n-by-m confusion matrix with the numbers of cells of known
#'     class 1..n (\code{ClassAssignment}, rows) classified by \code{griph} as cell type
#'     1..m (columns).
#' @slot miscl A length-one numerical value giving the misclassification error.
#' @slot GRAO An object of class \code{igraph} with c nodes representing cells
#'     and edges between similar cells.
GriphResult <- setClass("GriphResult",
                        slots = representation(DM = "matrix",
                                               is.cor = "logical",
                                               ClassAssignment = "numeric",
                                               BatchAssignment = "numeric",
                                               MEMB = "numeric",
                                               DISTM = "matrix",
                                               specp = "matrix",
                                               ConfMatrix = "matrix",
                                               miscl = "numeric",
                                               GRAO = "igraph")
)

#' @param object Instance of class \code{GriphResult}.
#' @describeIn GriphResult Print a summary of a GriphResult object.
setMethod("show", "GriphResult", function(object) {
    DMtype <- if(object@is.cor) "correlation coefficients" else "gene counts"
    cat("GriphResult\n", sep="")
    cat("Options   : DM                   : ", nrow(object@DM), "-by-", ncol(object@DM), " matrix of ", DMtype,
      "\n            known cell types     : ", if (length(unique(object@ClassAssignment)) > 1) "yes" else "no",
      "\n            known batches        : ", if (!is.null(object@BatchAssignment)) "yes" else "no",
      "\nResults   : identified cell types: ", length(unique(object@MEMB)),
      "\n            missclassification   : ", object@miscl, "\n\n", sep="")
})
