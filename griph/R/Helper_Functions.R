#' Get k nearest neighbors of variables given a feature matrix. 
#' 
#' @param S similarity matrix (here: always square)
#' @param k integer giving the number of nearest
#'     neighbors to be returned. Defaults to the
#'     rounded square root of the number of rows (columns).
#' 
#' @return A \code{k} by \code{ncol(S)} matrix.
get.knn <- function(S, k = round(sqrt(nrow(S)) ) ) {
    diag(S) <- 0
    kN <- sapply(1:nrow(S), function(x) utils::tail(order(S[,x]), k))
    return(kN)
}




#' Keep k mutual nearest neighbors of variables given a sparse similarity matrix. 
#' 
#' @param S  a sparse similarity matrix (square, symmetric, non-negative) of class CsparseMatrix
#' @param k integer giving the number of mutual nearest meighbors
#'      Defaults to the rounded square root of the number of variables.
#' 
#' @return a sparse Matrix of class CsparseMatrix with at most \code{length(S@x)} non-zero elements.
keep.mknn <- function(S, k = round(sqrt(nrow(S) )  )  ) {
  n <- diff(S@p)
  lst <- split(S@x, rep.int(1:ncol(S), n))
  o <- lapply(lst, function(x) rank(-x))
  o <- unlist(o)
  S@x[o > k] <- 0
  
  S@x <- pmin(S@x, t(S)@x)
  #S=sqrt( S * t(S) )  #Slower than line above, only works for strictly positive S
  
  S <- Matrix::drop0((S))
  return( S ) 
 }



#' Keep k mutual nearest neighbors of variables given a sparse similarity matrix.  
#' By default 10 nearest neighbors of each node are retained to avoid node isolation
#' 
#' @param S  a sparse similarity matrix (square, symmetric, non-negative) of class CsparseMatrix
#' @param k integer giving the number of mutual nearest meighbors
#'      Defaults to the rounded square root of the number of variables.
#' @param quant  Edges with weights > quantile(edge_weights, quant)  are retained in order to avoid large breaks
#' 
#' @return a sparse Matrix of class CsparseMatrix with at most \code{length(S@x)} non-zero elements.

keep.mknn2 <- function(S, k = round(sqrt(nrow(S) ) ), quant=0.8   ) {
    n <- diff(S@p)
    thr=quantile(S@x,quant)
    lst <- split(S@x, rep.int(1:ncol(S), n))
    o <- lapply(lst, function(x) rank(-x))
    o <- unlist(o)
    xtemp=S@x
    S@x[o > k] <- 0
    S@x[xtemp > thr] <- xtemp[xtemp > thr] #Keep edges with weights > quantile(edge_weights, quant)   => Avoid large breaks
    S@x <- pmin(S@x, t(S)@x)
    S@x[o < 11] <- xtemp[o < 11]             #Keep the top 10 edges of every node                        => Avoid isolated nodes
    S@x <- pmax(S@x, t(S)@x)
    S <- Matrix::drop0((S))
    return( S ) 
}




#' Prune an affinity (weighted adjacency) matrix to keep only top \code{pct} fraction neighbors  
#' 
#' @param S a sparse similarity matrix (square, symmetric, non-negative) of class CsparseMatrix.
#' @param pct fraction of edges (per node) to keep.
#' 
#' @return a sparse Matrix of class CsparseMatrix with a fraction of \code{pct} non-zero elements per column.
sparsify <- function(S, pct = 0.1) {
    n <- diff(S@p)
    lst <- split(S@x, rep.int(1:ncol(S), n))
    o <- lapply(lst, function(x) {
        y <- rank(-x)
        if (length(y) > 2) {
            y[ y > max(3,min(100,length(y)*pct))] <- 0
        }
        return(y)
    })
    o <- unlist(o)
    S@x[o == 0] <- 0
    #S@x=pmin(S@x, t(S)@x)
    S <- Matrix::drop0((S))
    return( as(S, "dgCMatrix")  )
}



#' Prune an affinity (weighted adjacency) matrix to keep only top \code{pct} fraction neighbors  
#' 
#' @param S a sparse similarity matrix (square, symmetric, non-negative) of class CsparseMatrix.
#' @param quant fraction of edges to keep.
#' 
#' @return a sparse Matrix of class CsparseMatrix with a fraction of \code{pct} non-zero elements per column.


sparsify2 <- function(S, quant = 0.1) {
    n <- diff(S@p)
    lst <- split(S@x, rep.int(1:ncol(S), n))
    thr=quantile(S@x,quant)
    
    o <- lapply(lst, function(x) {
        y <- x
        if (length(y) > 9) {
            y[ y < thr] <- 0
        }
        return(y)
    })
    o <- unlist(o)
    S@x[o == 0] <- 0
    S@x=pmax(S@x, t(S)@x)
    S <- Matrix::drop0((S))
    return( as(S, "dgCMatrix")  )
}











#' Correlation for sparse matrices
#' 
#' @param x a sparse matrix with n columns
#' 
#' @return an n-by-n symmetric matrix with Pearson's correlation coefficients between
#'     the columns of the input matrix \code{x}.
sparse.cor <- function(x){
    n <- nrow(x)
    cMeans <- colMeans(x)
    covmat <- (as.matrix(crossprod(x)) - n*tcrossprod(cMeans))/(n - 1)
    sdvec <- sqrt(diag(covmat)) 
    cormat <- covmat/tcrossprod(sdvec)
}



#' Pairwise correlation for 2 sparse matrices (with the same number of rows)
#' 
#' @param x a sparse matrix.
#' @param y a sparse matrix (same number of rows as \code{x}).
#' 
#' @return a matrix with correlations between columns of the inputs \code{x} and \code{y}.
psparse.cor <- function(x, y){
    n <- nrow(x)
    v <- t(scale(y))
    sdvec <- sqrt(colSums(x^2)/(n - 1) - colSums(x)^2/(n^2 - n)) 
    cormat <- tcrossprod(t(x)/sdvec, v)/(n - 1)
    return( as(cormat,"matrix") ) 
}





#' This is the function for selection of overdispersed genes adapted from:
#' https://github.com/10XGenomics/single-cell-3prime-paper/
#'
#' @param m  a (protentially sparse) gene x cells count matrix
#' @return a vector of normalized (robust Z-scores) dispersion values, one per gene.
select_variable_genes <- function(m) {
    df <- data.frame(mean = rowMeans(m + 1/ncol(m)), cv = apply(m,1,sd) / rowMeans(m + 1/ncol(m)), var = apply(m,1,var))
    df$dispersion <- with(df, var/mean)
    df$mean_bin <- with(df, cut(mean, breaks = c(-Inf, unique(quantile(mean, seq(0.1,1,0.05), na.rm = TRUE)), Inf)))
    var_by_bin <- data.frame(mean_bin = factor(levels(df$mean_bin), levels = levels(df$mean_bin)),
                             bin_median = as.numeric(tapply(df$dispersion, df$mean_bin, stats::median)),
                             bin_mad = as.numeric(tapply(df$dispersion, df$mean_bin, stats::mad)))[table(df$mean_bin) > 0,]
    df$bin_disp_median <- var_by_bin$bin_median[match(df$mean_bin, var_by_bin$mean_bin)]
    df$bin_disp_mad <- var_by_bin$bin_mad[match(df$mean_bin, var_by_bin$mean_bin)]
    df$dispersion_norm <- with(df, (dispersion - bin_disp_median)/(bin_disp_mad + 0.01) )
    return(df$dispersion_norm)
}












#' Emulate ggplot2 color palette.
#' 
#' @param n numeric(1) specifying the number of colors to return.
#' 
#' @return A character verctor with \code{n} colors spreading the whole rainbow.
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Darken or lighten colors.
#' 
#' @param col Character vector with valid R colors that \code{col2rgb()} accepts.
#' @param f Factor controlling strength of color change. A value of 1.0 does not
#'     change the colors, values < 1.0 make them darker, values > 1.0 ligther.
#'
#' @return A character vector with altered colors of the same length as \code{col}.
adjust.color <- function(col, f) {
    col <- grDevices::col2rgb(col)
    if (f < 1.0) {        # shade colors (make them darker)
        col <- round(col*f)
    } else if (f > 1.0) { # tint colors (make them lighter)
        col <- round(pmin(col + (255 - col) / f, 255))
    }
    apply(col, 2, function(x) grDevices::rgb(x[1], x[2], x[3], maxColorValue = 255))
}

#' Generate a random alpha-numeric string
#' 
#' @details
#' Generate a character vector of \code{n} random strings of length \code{len}
#' consisting of alpha-numeric characters (0-9, a-z or A-Z), for example to use
#' as output files prefixes.
#' 
#' @param n numeric(1) specifying the number of strings to generate (default: 1)
#' @param len the number of characters per string (default: 5)
#' 
#' @return a character vector with \code{n} random strings with \code{len} characters each.
RandString <- function(n=1, len=5){
  randomString <- c(1:n)  # initialize vector
  for (i in 1:n)  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    len, replace = TRUE), collapse = "")
  }
  return(randomString)
}


#' Generate custom vertices for igraph plotting
#' 
#' @description 
#' A plotting function for graph nodes used in graph visualization, passed to
#' \code{igraph::add.vertex.shape(..., plot=mycircle)}.
#' 
#' @param coords The coordinates of the vertices, a matrix with two columns.
#' @param v The ids of the vertices to plot. It should match the number of rows in the coords argument.
#' @param params This is a function object that can be called to query graphical parameters.
#' 
#' @return This function is used for its side effect of plotting graph vertices.
#' 
#' @seealso \code{\link[igraph]{add.vertex.shape}}.
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN = function(x, y, bg, fg, size, lwd) {
           symbols(x = x, y = y, bg = bg, fg = fg, lwd = lwd,
                   circles = size, add = TRUE, inches = FALSE)
         })
}


#' @title Map class labels using a greedy heuristic.
#' 
#' @description 
#' Map class labels (levels) of a given classification \code{a} in a greedy
#' fashion to its corresponding labels (levels) of a second classifciation \code{b}.
#' 
#' @details 
#' Each label in \code{a} is assigned to the label from \code{b} that it most
#' frequently co-occurs with. If there are multiple equally frequent pairs,
#' an arbitrary one of them is selected.
#' 
#' @param a Either a numeric, charactor of factor with class labels, or if \code{b}
#'     is \code{NULL}, a matrix or two-dimensional array (e.g. returned by \code{table})
#'     corresponding to the confusion matrix.
#' @param b Either a numeric, charactor of factor with class labels of the same
#'     length as \code{a}, or \code{NULL} if \code{a} is the confusion matrix.
#'
#' @return A character vector with the mapping of labels from \code{a} (names) to
#'     labels in \code{b} (values, coerced to \code{character}).
#' 
#' @seealso \code{\link{mapLabelsExhaustive}} for an exhaustive version that is
#'     guaranteed to return the optimal mapping that minimizes the classification
#'     error. \code{\link{classError}} to calculate classification error. 
mapLabelsGreedy <- function(a, b=NULL) {
    if (is.null(b)) {
        cm <- a # assume a to be a confusion matrix
    } else {
        stopifnot(length(a) == length(b))
        a <- as.character(a)
        b <- as.character(b)
        cm <- table(a, b)
    }
    map <- rownames(cm)
    names(map) <- map
    ub <- colnames(cm)
    for (i in seq_along(map))
        map[i] <- ub[which.max(cm[i,])]
    return(map)
}

#' @title Map class labels using an exhaustive algorithm.
#' 
#' @description 
#' Map class labels (levels) of a given classification \code{a} in an exhaustive
#' fashion to its corresponding labels (levels) of a second classifciation \code{b}.
#' 
#' @details 
#' Each label in \code{a} is assigned to the label from \code{b} that minimizes
#' the classification error, exhaustively searching over all possible permutations
#' of label assignments.
#' 
#' @param a Either a numeric, charactor of factor with class labels, or if \code{b}
#'     is \code{NULL}, a matrix or two-dimensional array (e.g. returned by \code{table})
#'     corresponding to the confusion matrix.
#' @param b Either a numeric, charactor of factor with class labels of the same
#'     length as \code{a}, or \code{NULL} if \code{a} is the confusion matrix.
#'
#' @return A character vector with the mapping of labels from \code{a} (names) to
#'     labels in \code{b} (values, coerced to \code{character}).
#' 
#' @seealso \code{\link{mapLabelsGreedy}} for a greedy version that is
#'     faster especially for many classes. \code{\link{classError}} to calculate classification error. 
mapLabelsExhaustive <- function(a, b=NULL) {
    #.permutations <- function( x, prefix = c() ) {
    #    if(length(x) == 0 )
    #        return(prefix)
    #    do.call(rbind, sapply(1:length(x), FUN = function(idx) .permutations( x[-idx], c( prefix, x[idx])), simplify = FALSE))
    #}

    if (!is.element("gtools", utils::installed.packages()[,1]))
        stop('"mapLabelsExhaustive" requires the "gtools" package. Please install it with:\n\t',
             'install.packages("gtools")')
    
    if (is.null(b)) {
        cm <- a # assume a to be a confusion matrix
    } else {
        stopifnot(length(a) == length(b))
        a <- as.character(a)
        b <- as.character(b)
        cm <- table(a, b)
    }
    
    # permC <- .permutations(rownames(cm)) # does not work if length(unique(a)) != length(unique(b))
    permC <- gtools::permutations(n = ncol(cm), r = nrow(cm), v = colnames(cm), set = FALSE, repeats.allowed = TRUE)
    colnames(permC) <- rownames(cm)
    permC[which.min(sapply(seq.int(nrow(permC)), function(i) classError(a, b, permC[i,]))), ]
}

#' @title Classification error.
#' 
#' @description 
#' Calculate the error of a given classification relative to a second one.
#' 
#' @details
#' The classification error is calculated as the fraction of elements, for which
#' the label in \code{a} is not equal to the label given in \code{b}. Optionally,
#' a mapping between the two label sets can be specificed using \code{map}.
#' 
#' @param a A numeric, charactor of factor with class labels.
#' @param b A numeric, charactor of factor with class labels (same length as \code{a}).
#' @param map A mapping from labels in \code{a} to labels in \code{b}, for example
#'     returned by \code{mapLabelsGreedy}. If \code{map} is \code{NULL},
#'     \code{mapLabelsExhaustive} (if \code{exhaustive} is \code{TRUE}) or
#'     \code{mapLabelsGreedy} is called to create such a mapping.
#' @param exhaustive If \code{TRUE} and \code{map} is \code{NULL}, use \code{mapLabelsExhaustive}
#'     to get \code{map}. If \code{FALSE} and \code{map} is \code{NULL}, use \code{mapLabelsGreedy}.
#'     Otherwise, \code{exhaustive} is ignored.
#'     
#' @return A numerical value corresponding to the classification error.
#'     
#' @seealso \code{\link{mapLabelsExhaustive}} and \code{\link{mapLabelsGreedy}} for
#'     calculation of a mapping between \code{a} and \code{b} labels.
classError <- function(a, b, map=NULL, exhaustive=FALSE) {
    if (is.null(map))
        map <- if (exhaustive) mapLabelsExhaustive(a, b) else mapLabelsGreedy(a, b)
    am <- map[as.character(a)]
    mean(am != as.character(b))
}


#' @title Structuredness of a classification.
#' 
#' @description Calculate a score that measures the structuredness of a classification.
#' 
#' @details The score is the value of a statistic (e.g. the standard deviation of log2
#'     fold-change of average gene expression levels within a class over the mean
#'     of all classes) using the raw data and a given classification, relative to
#'     the same statistic obtained from the raw data and permuted classifications.
#' 
#' @param DM Gene-by-cell expression matrix (raw counts).
#' @param classification Factor, numerical or character vector with class labels.
#' @param score.type Statistical measure used to calculate the score. One of sdLogFC".
#' @param R Integer scalar defining the number of permutations on the classification
#'     to perform for normalization of the statistical measure.
#' 
#' @return A list with components "score.type" (the selected measurement statistic),
#'     "score.obs" (value of the measurment statistic), "score.rand" (values for
#'     randomized classifications) and "score.norm" (the
#'     normalized value of the measurment statistic).
clusteringScore <- function(DM, classification, score.type = c("sdLogFC"), R = 20) {
    # digest arguments
    stopifnot(is.matrix(DM))
    classification <- as.numeric(factor(classification))
    stopifnot(length(classification) == ncol(DM))
    score.type <- match.arg(score.type)
    
    # calculate measurment for real classification
    calcStat <- switch(score.type,
                       sdLogFC = function(x, cl) {
                           # create "fake" bulk profiles
                           jByClass <- split(seq.int(ncol(x)), cl)
                           xc <- do.call(cbind, lapply(jByClass, function(j) rowSums(x[, j, drop = FALSE])))
                           # normalize
                           xl <- log2(t(t(xc) / colSums(xc)) * 1e6 + 1.0) # log2(cpm+pseudocount)
                           xr <- xl - rowMeans(xl) # relative to the mean
                           # score
                           sd(xr)
                       })
    message("calculating score for real classification...", appendLF = FALSE)
    val.obs <- calcStat(DM, classification)
    message("done")
    
    # calculate measurment for randomized classifications
    message("calculating scores for ",R," randomized classifications...", appendLF = FALSE)
    val.rand <- unlist(lapply(seq.int(R), function(i) calcStat(DM, sample(classification))))
    message("done")
    
    # return results
    list(score.type = score.type, score.obs = val.obs, score.rand = val.rand,
         score.norm = val.obs / median(val.rand))
}




#' @title Edge reweighting to group nodes by community during plotting
#' 
#' @description 
#' Edges crossing communites are downweighted while edges within 
#' a community are unaffected (or up-weighted). As a result during plotting
#' with an \code{igraph} layout algorithm that takes into account edge
#' weights (such as \code{layout_with_fr}) communities are brought 
#' together and stand out better in the plot. This operation distorts
#' the "true node distances" to improve community visibility
#' 
#' @param community A community object returned by any community detection
#' function of \code{igraph}
#' @param network A graph object of class \code{igraph}
#' @param weight.within A scalar specifying the weight for edges connecting nodes
#' within the same community
#' @param weight.between A scalar specifying the weight for edges connecting nodes
#' between different communities (communities crossing edges)
#' @return A numeric vector of length \code{length(E(network))} with the new
#' edge weights
#' 
edge.weights <- function(community, network, weight.within = 10, weight.between = 1) {
    bridges <- igraph::crossing(communities = community, graph = network)
    weights <- ifelse(test = bridges, yes = weight.between, no = weight.within)
    weights <- weights * igraph::E(network)$weight
    return(weights) 
}


