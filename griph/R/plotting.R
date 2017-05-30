#' @title Visualize griph result as a graph.
#' 
#' @description Plot a graph obtained from \code{\link{SC_cluster}}, allowing to
#'     apply graph subsampling and control coloring.
#' 
#' @param gr A \code{griph} result, as returned by \code{\link{SC_cluster}}.
#' @param maxG Approximate maximal number of vertices to include when plotting
#'     the graph.
#' @param fill.type Type of fill color, one of \code{predicted} (predicted class
#'     labels, default), \code{true} (true class labels, if available), \code{none}
#'     (no fill color) or \code{custom} (use \code{custom.col}).
#' @param line.type Type of line color, one of \code{true} (true class labels, if
#'     available, default), \code{predicted} (predicted class labels), \code{none}
#'     for no fill color or \code{custom} (use \code{custom.col}).
#' @param mark.type Type of cell class defnition to mark using polygons,
#'     one of \code{none} (no polygons, the default), \code{predicted} (draw
#'     polygons around cells with the same predicted class label), \code{true}
#'     (polygons around cells with the same true class label, if available) or
#'     \code{custom} (polygons around cells with the same \code{custom.class}).
#' @param collapse.type Type of cell class to use for graph simplification,
#'     by combining cells of the same class into a single vertex. If set to a value
#'     other than \code{"none"}, the same value will also be used for \code{fill.type}
#'     and \code{line.type} is ignored.
#' @param fill.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for vertex fill coloring.
#' @param line.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for vertex outline coloring.
#' @param mark.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for cell class polygon marking.
#' @param custom.class Factor, character or numeric vector of the same length or
#'     with names corresponding to names(gr$MEMB) to use for custom cell classification
#'     (used if \code{fill.type} and/or \code{line.type} is set to "custom").
#' @param draw.edges If \code{NULL} (default), draw edges if \code{collapse.type != "none"}.
#'     \code{TRUE} or \code{FALSE} can be used to override the default.
#' @param seed Random number seed to make graph layout deterministic.
#' @param fsuffix A suffix added to the file names of output plots. If not given
#'     it will use a random 5 character string. Ignored if \code{image.format} is \code{NULL}.
#' @param image.format Specifies the format of the created image. Currently supported are
#'     \code{\link{pdf}}, \code{\link{png}} or \code{NA}. If \code{NA} (the default), the plot
#'     is rendered on the currently opened plotting device.
#' @param forceRecalculation If \code{TRUE}, recalculate plotting-optimized graph
#'     even if it is already contained in \code{gr}.
#' @param quiet If \code{TRUE}, do not report on progress.
#' 
#' @return Invisibly the plot-optimized version of the graph as an \code{igraph} object.
plotGraph <- function(gr, maxG=2500,
                      fill.type=c("predicted","true","none","custom"),
                      line.type=c("true","predicted","none","custom"),
                      mark.type=c("none","predicted","true","custom"),
                      collapse.type=c("none","predicted","true"),
                      fill.col="Spectral",
                      line.col="Dark2",
                      mark.col="Pastel1",
                      custom.class=factor(rep(1, length(gr$MEMB))),
                      draw.edges=NULL,
                      seed=91919,
                      fsuffix=RandString(), image.format=NA,
                      forceRecalculation=FALSE, quiet=FALSE) {
    if (!quiet)
        message("Computing Graph Layout and Rendering...")
    
    # get varaibles from gr
    GRAO <- gr$GRAO
    MEMB <- gr$MEMB
    csize <- table(MEMB)
    
    # digest arguments
    fill.type <- match.arg(fill.type)
    line.type <- match.arg(line.type)
    mark.type <- match.arg(mark.type)
    collapse.type <- match.arg(collapse.type)
    if (is.null(names(custom.class))) {
        stopifnot(length(custom.class) == length(MEMB))
        names(custom.class) <- names(MEMB)
    }
    if (length(fill.col) == 1) {
        stopifnot(fill.col %in% rownames(RColorBrewer::brewer.pal.info))
        fill.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[fill.col, "maxcolors"], fill.col)
    }
    if (length(line.col) == 1) {
        stopifnot(line.col %in% rownames(RColorBrewer::brewer.pal.info))
        line.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[line.col, "maxcolors"], line.col)
    }
    if (length(mark.col) == 1) {
        stopifnot(mark.col %in% rownames(RColorBrewer::brewer.pal.info))
        mark.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[mark.col, "maxcolors"], mark.col)
    }
    
    # global plotting paramterers
    pct <- 1
    my.pch <- 21L # should be in 21:25
    my.pt.cex <- 2.5
    my.pt.lwd <- if (line.type == "none" || collapse.type != "none") 1.0 else 2.5
    edge.lwd.max <- 12.0
    edge.col <- "#33333355"
    
    # get plot-optimized graph
    if (collapse.type != "none") {
        if (!quiet)
            message("\tcollapsing cells by ",collapse.type," class")
        if (fill.type != collapse.type && !quiet)
            message("\tSetting fill.type to collapse.type for collapse.type != 'none'.")
        if ((line.type != "none" || mark.type != "none") && !quiet)
            message("\tSetting line.type and mark.type to 'none' for collapse.type != 'none'.")
        fill.type <- collapse.type
        line.type <- "none"
        mark.type <- "none"
        class.collapse <- switch(collapse.type,
                                 predicted = factor(V(GRAO)$membership,
                                                    levels = sort(as.numeric(unique(V(GRAO)$membership)))),
                                 true = factor(V(GRAO)$class,
                                               levels = unique(V(GRAO)$class)))
        GRAOp <- igraph::simplify( igraph::contract(GRAO, class.collapse) )
        V(GRAOp)$membership <- V(GRAOp)$class <- levels(class.collapse)
        V(GRAOp)$size <- 10 / (length(V(GRAOp))/60 )^0.3 * (as.numeric(csize/median(csize)))^0.5
        
    } else if (is.null(gr$plotGRAO) || forceRecalculation) {
        if (length(V(GRAO) ) > 1.25*maxG ) {
            if (!quiet)
                message("\tGraph too large (>",maxG, " vertices). A sampled subgraph of ~", maxG, " vertices will be plotted", appendLF = FALSE)
            
            ###### Sample well-connected seeds from the members of each community 
            DEG <- igraph::degree(GRAO)
            snowball_seeds <- c()
            for (i in 1:length(csize)) {
                if (csize[i] < 5) {next}
                members <- which(MEMB == i)
                minDegree <- quantile(DEG[members])[2] - 1
                maxDegree <- quantile(DEG[members])[4] + 1
                #seedn <- ceiling(csize[i]/max(csize))
                seedn <- ceiling(5 * (sqrt(csize[i] - 4) / sqrt(max(csize) - 4)))
                seedn <- min(seedn,floor(csize[i]/4) )
                #message("minD:",minDegree, " maxD:",maxDegree," csize:", csize[i]  ,"\r")
                if (seedn > 1) {
                    module_seeds <- sample(which(MEMB == i & DEG >= floor(minDegree) &
                                                 DEG <= ceiling(maxDegree)), seedn)
                } else {
                    module_seeds <- sample(which(MEMB == i & DEG == max(DEG[members])), 1)    
                }
                
                snowball_seeds <- unique(c(snowball_seeds,module_seeds))
            }
            
            snowball <- c()
            seed.ego_size <- 0
            while (length(snowball) < maxG/2) {
                seed.ego_size <- seed.ego_size + 1  
                snowball <- unique(unlist(igraph::ego(GRAO,seed.ego_size,snowball_seeds)))
            }
            if (length(snowball) > 1.25*maxG && seed.ego_size > 1 ) {
                seed.ego_size <- seed.ego_size - 1
                snowball <- unique(unlist(igraph::ego(GRAO,seed.ego_size,snowball_seeds))) 
            }
            
            GRAOp <- igraph::induced.subgraph(GRAO,sort(snowball) )
            
            if (!quiet)
                message(" (used vertices: ", length(V(GRAOp)),"  seed_size: ",seed.ego_size,")")
            
        } else {
            if (!quiet)
                message("\tusing full graph for plotting")
            GRAOp <- GRAO
        }
        
        ######## Prune graph for better plot output
        if (median(igraph::degree(GRAOp)) > 4 ) {
            pct <- min(1, 1 / sqrt(0.1 * median(igraph::degree(GRAOp))))
            ADJp <- as.matrix(igraph::get.adjacency(GRAOp, attr = 'weight'))
            ADJp <- apply(ADJp, 1, function(x) sparsify(x,pct))
            ADJp[which(abs(ADJp) > 0)] <- 1
            GRAOtemp <- igraph::graph.adjacency(ADJp, mode = "max", weighted = NULL,
                                                diag = FALSE)
            GRAOp <- igraph::intersection(GRAOp, GRAOtemp, byname = FALSE)
            GRAOtemp <- NULL
            ADJp <- NULL
        }
        
        GRAOp <- igraph::delete_vertices(GRAOp, which(igraph::ego_size(GRAOp, 3) < 6))
        ###Delete Vertices from communites with few members:
        min.csize <- ceiling(0.25 * sqrt(length(V(GRAO))))
        GRAOp <- igraph::delete_vertices(GRAOp, which( V(GRAOp)$community.size < min.csize ))  
        
        if (!quiet)
            message("\tnodes from communities with <",min.csize, " members will not be displayed.")
        
    } else {
        if (!quiet)
            message("\tusing existing plot-optimized graph")
        GRAOp <- gr$plotGRAO
    }
    if (!quiet)
        message("\tdisplaying graph with ",length(V(GRAOp))," (",
                round(100 * length(V(GRAOp)) / length(V(GRAO)), 1), "%) vertices and ",
                length(E(GRAOp)), " (", round(100 * length(E(GRAOp)) / length(E(GRAO)), 1),
                "%) edges")
    
    # get colors
    class.pred <- factor(V(GRAOp)$membership, levels = sort(as.numeric(unique(V(GRAO)$membership))))
    class.true <- factor(V(GRAOp)$class, levels = unique(V(GRAO)$class))
    class.none <- factor(rep(NA, length(V(GRAOp))))
    class.custom <- factor(custom.class[V(GRAOp)$labels], levels = unique(custom.class))
    
    class.fill <- switch(fill.type,
                         predicted = class.pred,
                         true = class.true,
                         none = class.none,
                         custom = class.custom)
    fillColorPalette <- if (nlevels(class.fill) > length(fill.col)) grDevices::colorRampPalette(fill.col)(nlevels(class.fill)) else fill.col
    fillColor <- fillColorPalette[as.numeric(class.fill)]
    
    class.line <- switch(line.type,
                         predicted = class.pred,
                         true = class.true,
                         none = class.none,
                         custom = class.custom)
    lineColorPalette <- if (nlevels(class.line) > length(line.col)) grDevices::colorRampPalette(line.col)(nlevels(class.line)) else line.col
    lineColor <- lineColorPalette[as.numeric(class.line)]
    
    class.mark <- switch(mark.type,
                         predicted = class.pred,
                         true = class.true,
                         none = class.none,
                         custom = class.custom)
    markElements <- split(seq_along(class.mark), class.mark)
    markColor <- if (nlevels(class.mark) > length(mark.col)) grDevices::colorRampPalette(mark.col)(nlevels(class.mark)) else mark.col[1:nlevels(class.mark)]
    
    # set some more graph attributes
    V(GRAOp)$classcolor <- lineColor
    V(GRAOp)$size <- if (collapse.type == "none") 10 / (length(V(GRAOp))/60 )^0.3 else V(GRAOp)$size
    V(GRAOp)$cex <- V(GRAOp)$size / 3
    V(GRAOp)$frame.width <- 2 / (length(V(GRAOp))/60 )^0.3
    E(GRAOp)$width <- E(GRAOp)$weight / sqrt((length(V(GRAOp))/60 ))
    V(GRAOp)$color <- fillColor
    
    # compute graph layout
    set.seed(seed = seed)
    l <- igraph::layout_with_fr(graph = GRAOp)
    # igraph::add.vertex.shape("fcircle", clip=igraph.shape.noclip, plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))
    
    # open output file
    if (!is.na(image.format)) {
        if (image.format == 'pdf') {
            fname <- paste('graph_', fsuffix, '.pdf', sep = "")
            pdf(file = fname, width = 12, height = 10)
        } else if (image.format == "png") {
            fname <- paste('graph_', fsuffix, '.png', sep = "")
            png(filename = fname, width = 12, height = 10, units = "in", res = 300)   
        }
        if (!quiet)
            message("\tsaving graph to ", fname)
    }
    
    # setup plot coordinate system
    par(mar = c(5.1, 4.1, 4.1, 14.1), xpd = TRUE)
    plot(l[,1], l[,2], type = "n", axes = FALSE, xlab = "", ylab = "")
    
    # add mark polygons
    if (mark.type != "none") {
        for (j in which(lengths(markElements) > 0)) {
            xy <- l[markElements[[j]], , drop = FALSE]
            off <- par("cxy")[2] * 1
            #avg <- matrix(colMeans(xy), ncol = 2, nrow = nrow(xy), byrow = TRUE)
            #pp <- xy + sign(xy - avg) * off
            pp <- rbind(xy,
                        cbind(xy[, 1] - off, xy[, 2]),
                        cbind(xy[, 1] + off, xy[, 2]),
                        cbind(xy[, 1], xy[, 2] - off), 
                        cbind(xy[, 1], xy[, 2] + off))
            cl <- igraph::convex_hull(pp)
            graphics::xspline(cl$rescoords, shape = 0.5, open = FALSE,
                              col = paste0(markColor[j], "66"),
                              border = adjust.color(markColor[j], 0.5))
            
        }
    }
    
    # add edges (by default only for collapase.type != "none")
    if (isTRUE(draw.edges) || (is.null(draw.edges) && collapse.type != "none")) {
        el <- igraph::as_edgelist(GRAOp, names = FALSE)
        graphics::segments(x0 = l[,1][el[,1]], y0 = l[,2][el[,1]],
                           x1 = l[,1][el[,2]], y1 = l[,2][el[,2]],
                           col = edge.col, lwd = E(GRAOp)$weight / max(E(GRAOp)$weight) * edge.lwd.max)
    }
    
    # add vertices
    points(l[,1], l[,2], col = if (line.type == "none") "black" else lineColor,
           bg = fillColor, pch = my.pch, lwd = my.pt.lwd,
           cex = my.pt.cex * if (collapse.type == "none") 1.0 else (as.numeric(csize / median(csize)))^0.5)
    
    # add legend(s)
    if (mark.type != "none") {
        lgd <- legend(x = par("usr")[2] + 12*par("cxy")[1], y = par("usr")[4],
                      xjust = 1, yjust = 1, bty = "n", cex = 1,
                      fill = paste0(markColor, "66"), title = mark.type,
                      legend = levels(class.mark))
    } else {
        lgd <- list(rect = list(left = par("usr")[2] + 12*par("cxy")[1]))
    }
    if (fill.type != "none" && nlevels(class.fill) > 0) {
        lgd <- legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
                      pch = my.pch, pt.lwd = my.pt.lwd, cex = 1, pt.cex = my.pt.cex,
                      col = if (line.type == "none") "black" else "white", pt.bg = fillColorPalette,
                      title = fill.type, legend = levels(class.fill))
    }
    if (line.type != "none" && nlevels(class.line) > 0) {
        legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
               pch = my.pch, pt.lwd = my.pt.lwd, cex = 1, pt.cex = my.pt.cex,
               col = lineColorPalette, pt.bg = "white",
               title = line.type, legend = levels(class.line))
    }
    
    # close output file
    if (!is.na(image.format))
        dev.off()
    
    return(invisible(GRAOp))
}

#' @title Visualize griph result as a tSNE projection.
#' 
#' @description Plot a t-SNE projection of the affinity matrix obtained from
#'     \code{\link{SC_cluster}}, allowing to control coloring.
#' 
#' @param gr A \code{griph} result, as returned by \code{\link{SC_cluster}}.
#' @param fill.type Type of fill color, one of \code{predicted} (predicted class
#'     labels, default), \code{true} (true class labels, if available), \code{none}
#'     (no fill color) or \code{custom} (use \code{custom.col}).
#' @param line.type Type of line color, one of \code{true} (true class labels, if
#'     available, default), \code{predicted} (predicted class labels), \code{none}
#'     for no fill color or \code{custom} (use \code{custom.col}).
#' @param mark.type Type of cell class defnition to mark using polygons,
#'     one of \code{none} (no polygons, the default), \code{predicted} (draw
#'     polygons around cells with the same predicted class label), \code{true}
#'     (polygons around cells with the same true class label, if available) or
#'     \code{custom} (polygons around cells with the same \code{custom.class}).
#' @param fill.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for vertex fill coloring.
#' @param line.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for vertex outline coloring.
#' @param mark.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for cell class polygon marking.
#' @param custom.class Factor, character or numberic vector of the same length or
#'     with names corresponding to names(gr$MEMB) to use for custom cell classification
#'     (used if \code{fill.type} and/or \code{line.type} is set to "custom").
#' @param seed Random number seed to make t-SNE projection deterministic.
#' @param fsuffix A suffix added to the file names of output plots. If not given
#'     it will use a random 5 character string. Ignored if \code{image.format} is \code{NULL}.
#' @param image.format Specifies the format of the created image. Currently supported are
#'     \code{\link{pdf}}, \code{\link{png}} or \code{NA}. If \code{NA} (the default), the plot
#'     is rendered on the currently opened plotting device.
#' @param quiet If \code{TRUE}, do not report on progress.
#' @param ... additional arguments passed to \code{Rtsne}
#' 
#' @return Invisible the results of the t-SNE projection.
#' 
#' @seealso \code{Rtsne} used to calculate the t-SNE projection.
plotTsne <- function(gr,
                     fill.type=c("predicted","true","none","custom"),
                     line.type=c("true","predicted","none","custom"),
                     mark.type=c("none","predicted","true","custom"),
                     fill.col="Spectral",
                     line.col="Dark2",
                     mark.col="Pastel1",
                     custom.class=factor(rep(1, length(gr$MEMB))),
                     seed=91919,
                     fsuffix=RandString(), image.format=NA,
                     quiet=FALSE, ...) {
    if (!is.element("Rtsne", utils::installed.packages()[,1])) {
        stop('"plotTsne" requires the "Rtsne" package. Please install it with:\n\t',
             'install.packages("Rtsne")')
    }
    
    args <- list(...)
    # get varaibles from gr
    MEMB <- gr$MEMB
    MEMB.true <- gr$MEMB.true
    #csize <- table(MEMB)
    
    # digest arguments
    fill.type <- match.arg(fill.type)
    line.type <- match.arg(line.type)
    mark.type <- match.arg(mark.type)
    if (is.null(names(custom.class))) {
        stopifnot(length(custom.class) == length(MEMB))
        names(custom.class) <- names(MEMB)
    }
    if (length(fill.col) == 1) {
        stopifnot(fill.col %in% rownames(RColorBrewer::brewer.pal.info))
        fill.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[fill.col, "maxcolors"], fill.col)
    }
    if (length(line.col) == 1) {
        stopifnot(line.col %in% rownames(RColorBrewer::brewer.pal.info))
        line.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[line.col, "maxcolors"], line.col)
    }
    if (length(mark.col) == 1) {
        stopifnot(mark.col %in% rownames(RColorBrewer::brewer.pal.info))
        mark.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[mark.col, "maxcolors"], mark.col)
    }
    
    # global plotting paramterers
    my.pch <- 21L # should be in 21:25
    my.pt.cex <- 2.5
    my.pt.lwd <- if (line.type == "none") 1.0 else 2.5
    
    # get t-SNE projection
    if (!quiet)
        message("Computing t-SNE projection...")
    set.seed(seed = seed)
    if ( !("perplexity" %in% names(args))  ) {
        args$perplexity <- min(30,round(sqrt(nrow(gr$DISTM)) - 1))
    }
    res <- do.call(Rtsne::Rtsne, c(list(X = stats::as.dist(1 - gr$DISTM), pca = FALSE, is_distance = TRUE), args))
    
    # get colors
    class.pred <- factor(MEMB, levels = sort(as.numeric(unique(MEMB))))
    class.true <- factor(MEMB.true, levels = unique(MEMB.true))
    class.none <- factor(rep(NA, length(MEMB)))
    class.custom <- factor(custom.class, levels = unique(custom.class))
    
    class.fill <- switch(fill.type,
                         predicted = class.pred,
                         true = class.true,
                         none = class.none,
                         custom = class.custom)
    fillColorPalette <- if (nlevels(class.fill) > length(fill.col)) grDevices::colorRampPalette(fill.col)(nlevels(class.fill)) else fill.col
    fillColor <- fillColorPalette[as.numeric(class.fill)]
    
    class.line <- switch(line.type,
                         predicted = class.pred,
                         true = class.true,
                         none = class.none,
                         custom = class.custom)
    lineColorPalette <- if (nlevels(class.line) > length(line.col)) grDevices::colorRampPalette(line.col)(nlevels(class.line)) else line.col
    lineColor <- lineColorPalette[as.numeric(class.line)]
    
    class.mark <- switch(mark.type,
                         predicted = class.pred,
                         true = class.true,
                         none = class.none,
                         custom = class.custom)
    markElements <- split(seq_along(class.mark), class.mark)
    markColor <- if (nlevels(class.mark) > length(mark.col)) grDevices::colorRampPalette(mark.col)(nlevels(class.mark)) else mark.col[1:nlevels(class.mark)]
    
    # open output file
    if (!is.na(image.format)) {
        if (image.format == 'pdf') {
            fname <- paste('tSNE_',fsuffix,'.pdf',sep = "")
            pdf(file = fname, width = 12, height = 10)
        } else if (image.format == "png") {
            fname <- paste('tSNE_', fsuffix, '.png', sep = "")
            png(filename = fname, width = 12, height = 10, units = "in", res = 300)   
        }
        if (!quiet)
            message("\tsaving t-SNE plot to ",fname)
    }
    
    # setup plot coordinate system
    par(mar = c(5.1, 4.1, 4.1, 14.1), xpd = TRUE)
    plot(res$Y[,1], res$Y[,2], type = "n", axes = FALSE, xlab = "", ylab = "")
    
    # add mark polygons
    if (mark.type != "none") {
        for (j in which(lengths(markElements) > 0)) {
            xy <- res$Y[markElements[[j]], , drop = FALSE]
            off <- par("cxy")[2]*1
            pp <- rbind(xy,
                        cbind(xy[, 1] - off, xy[, 2]),
                        cbind(xy[, 1] + off, xy[, 2]),
                        cbind(xy[, 1], xy[, 2] - off), 
                        cbind(xy[, 1], xy[, 2] + off))
            cl <- igraph::convex_hull(pp)
            graphics::xspline(cl$rescoords, shape = 0.5, open = FALSE,
                              col = paste0(markColor[j], "66"),
                              border = adjust.color(markColor[j], 0.5))
            
        }
    }
    
    # add cells
    points(res$Y, col = if (line.type == "none") "black" else lineColor,
           bg = fillColor, pch = my.pch, lwd = my.pt.lwd, cex = my.pt.cex)
    
    # add legend(s)
    if (mark.type != "none") {
        lgd <- legend(x = par("usr")[2] + 12 * par("cxy")[1], y = par("usr")[4],
                      xjust = 1, yjust = 1, bty = "n", cex = 1,
                      fill = paste0(markColor, "66"), title = mark.type,
                      legend = levels(class.mark))
    } else {
        lgd <- list(rect = list(left = par("usr")[2] + 12 * par("cxy")[1]))
    }
    if (fill.type != "none" && nlevels(class.fill) > 0) {
        lgd <- legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
                      pch = my.pch, pt.lwd = my.pt.lwd, cex = 1, pt.cex = my.pt.cex,
                      col = if (line.type == "none") "black" else "white",
                      pt.bg = fillColorPalette,
                      title = fill.type, legend = levels(class.fill))
    }
    if (line.type != "none" && nlevels(class.line) > 0) {
        legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
               pch = my.pch, pt.lwd = my.pt.lwd, cex = 1, pt.cex = my.pt.cex,
               col = lineColorPalette, pt.bg = "white",
               title = line.type, legend = levels(class.line))
    }
    
    # close output file
    if (!is.na(image.format))
        dev.off()
    
    return(invisible(res))
}


#' @title Visualize griph result using a LargeVis algorithm (Tang et al., 2016).
#' 
#' @description Plot a largeVis projection of the affinity matrix obtained from
#'     \code{\link{SC_cluster}} or \code{\link{griph_cluster}}, allowing to control coloring.
#' 
#' @param gr A \code{griph} result, as returned by \code{\link{SC_cluster}}
#'     or \code{\link{griph_cluster}}.
#' @param fill.type Type of fill color, one of \code{predicted} (predicted class
#'     labels, default), \code{true} (true class labels, if available), \code{none}
#'     (no fill color) or \code{custom} (use \code{custom.col}).
#' @param line.type Type of line color, one of \code{true} (true class labels, if
#'     available, default), \code{predicted} (predicted class labels), \code{none}
#'     for no fill color or \code{custom} (use \code{custom.col}).
#' @param mark.type Type of cell class defnition to mark using polygons,
#'     one of \code{none} (no polygons, the default), \code{predicted} (draw
#'     polygons around cells with the same predicted class label), \code{true}
#'     (polygons around cells with the same true class label, if available) or
#'     \code{custom} (polygons around cells with the same \code{custom.class}).
#' @param fill.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for vertex fill coloring.
#' @param line.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for vertex outline coloring.
#' @param mark.col Character scalar with a \code{\link{RColorBrewer}} color palette name
#'     or color vector defining the palette to use for cell class polygon marking.
#' @param custom.class Factor, character or numberic vector of the same length or
#'     with names corresponding to names(gr$MEMB) to use for custom cell classification
#'     (used if \code{fill.type} and/or \code{line.type} is set to "custom").
#' @param seed Random number seed to make largeVis projection deterministic.
#' @param fsuffix A suffix added to the file names of output plots. If not given
#'     it will use a random 5 character string. Ignored if \code{image.format} is \code{NULL}.
#' @param image.format Specifies the format of the created image. Currently supported are
#'     \code{\link{pdf}}, \code{\link{png}} or \code{NA}. If \code{NA} (the default), the plot
#'     is rendered on the currently opened plotting device.
#' @param quiet If \code{TRUE}, do not report on progress.
#' @param ... additional arguments passed to \code{largeVis::projectKNNs}
#' 
#' @return Invisible the results of the largeVis projection.
#' 
#' @seealso \code{largeVis} used to calculate the largeVis projection.
plotLVis <- function(gr,
                     fill.type=c("predicted","true","none","custom"),
                     line.type=c("true","predicted","none","custom"),
                     mark.type=c("none","predicted","true","custom"),
                     fill.col="Spectral",
                     line.col="Dark2",
                     mark.col="Pastel1",
                     custom.class=factor(rep(1, length(gr$MEMB))),
                     seed=91919,
                     fsuffix=RandString(), image.format=NA,
                     quiet=FALSE, ...) {
    
    
    add.args=list(...)
    # get varaibles from gr
    MEMB <- gr$MEMB
    MEMB.true <- gr$MEMB.true
    csize <- table(MEMB)
    
    # digest arguments
    fill.type <- match.arg(fill.type)
    line.type <- match.arg(line.type)
    mark.type <- match.arg(mark.type)
    if(is.null(names(custom.class))) {
        stopifnot(length(custom.class)==length(MEMB))
        names(custom.class) <- names(MEMB)
    }
    if(length(fill.col)==1) {
        stopifnot(fill.col %in% rownames(RColorBrewer::brewer.pal.info))
        fill.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[fill.col, "maxcolors"], fill.col)
    }
    if(length(line.col)==1) {
        stopifnot(line.col %in% rownames(RColorBrewer::brewer.pal.info))
        line.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[line.col, "maxcolors"], line.col)
    }
    if(length(mark.col)==1) {
        stopifnot(mark.col %in% rownames(RColorBrewer::brewer.pal.info))
        mark.col <- RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[mark.col, "maxcolors"], mark.col)
    }
    
    # global plotting paramterers
    my.pch <- 21L # should be in 21:25
    my.pt.cex <- 2.5
    my.pt.lwd <- if(line.type == "none") 1.0 else 2.5
    
    # get largeVis projection
    if(!quiet)
        message("Computing largeVis projection...")
    
    set.seed(seed = seed)
    
    if (!is.element('sgd_batches',names(add.args))){
        add.args$sgd_batches=min(25000/(sum(gr$DISTM!=0)/2),0.99)  #!!!!!!Use GRAO instead!!! 
    }
    if (!is.element('M',names(add.args))){
        add.args$M=4    
    }
    if (!is.element('gamma',names(add.args))){
        add.args$gamma=1    
    }
    if (!is.element('alpha',names(add.args))){
        add.args$alpha=0.2
    }
    if (!is.element('useDegree',names(add.args))){
        add.args$useDegree=TRUE    
    }
    
    res <- do.call(projectKNNs, c( list(wij=as_adj(gr$GRAO,names=FALSE, sparse=TRUE),seed=seed),add.args )    )
    
    # get colors
    class.pred <- factor(MEMB, levels=sort(as.numeric(unique(MEMB))))
    class.true <- factor(MEMB.true, levels=unique(MEMB.true))
    class.none <- factor(rep(NA, length(MEMB)))
    class.custom <- factor(custom.class, levels=unique(custom.class))
    
    class.fill <- switch(fill.type,
                         predicted=class.pred,
                         true=class.true,
                         none=class.none,
                         custom=class.custom)
    fillColorPalette <- if(nlevels(class.fill) > length(fill.col)) grDevices::colorRampPalette(fill.col)(nlevels(class.fill)) else fill.col
    fillColor <- fillColorPalette[as.numeric(class.fill)]
    
    class.line <- switch(line.type,
                         predicted=class.pred,
                         true=class.true,
                         none=class.none,
                         custom=class.custom)
    lineColorPalette <- if(nlevels(class.line) > length(line.col)) grDevices::colorRampPalette(line.col)(nlevels(class.line)) else line.col
    lineColor <- lineColorPalette[as.numeric(class.line)]
    
    class.mark <- switch(mark.type,
                         predicted=class.pred,
                         true=class.true,
                         none=class.none,
                         custom=class.custom)
    markElements <- split(seq_along(class.mark), class.mark)
    markColor <- if(nlevels(class.mark) > length(mark.col)) grDevices::colorRampPalette(mark.col)(nlevels(class.mark)) else mark.col[1:nlevels(class.mark)]
    
    # open output file
    if(!is.na(image.format)) {
        if(image.format=='pdf') {
            fname <- paste('Lvis_',fsuffix,'.pdf',sep="")
            pdf(file = fname, width = 12, height = 10)
        } else if (image.format=="png") {
            fname <- paste('Lvis_',fsuffix,'.png',sep="")
            png(filename = fname, width = 12, height = 10, units = "in", res = 300)   
        }
        if(!quiet)
            message("\tsaving LargeVis plot to ",fname)
    }
    
    # setup plot coordinate system
    par(mar=c(5.1, 4.1, 4.1, 14.1), xpd=TRUE)
    plot(t(res)[,1], t(res)[,2], type = "n", axes=FALSE, xlab="", ylab="")
    
    # add mark polygons
    if(mark.type != "none") {
        for(j in which(lengths(markElements) > 0)) {
            xy <- t(res)[markElements[[j]], , drop = FALSE]
            off <- par("cxy")[2]*1
            pp <- rbind(xy,
                        cbind(xy[, 1] - off, xy[, 2]),
                        cbind(xy[, 1] + off, xy[, 2]),
                        cbind(xy[, 1], xy[, 2] - off), 
                        cbind(xy[, 1], xy[, 2] + off))
            cl <- igraph::convex_hull(pp)
            graphics::xspline(cl$rescoords, shape = 0.5, open = FALSE,
                              col = paste0(markColor[j], "66"),
                              border = adjust.color(markColor[j], 0.5))
            
        }
    }
    
    # add cells
    points(t(res), col = if(line.type=="none") "black" else lineColor,
           bg = fillColor, pch = my.pch, lwd = my.pt.lwd, cex=my.pt.cex)
    
    # add legend(s)
    if(mark.type != "none") {
        lgd <- legend(x = par("usr")[2]+12*par("cxy")[1], y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
                      cex = 1, fill = paste0(markColor, "66"), title = mark.type, legend = levels(class.mark))
    } else {
        lgd <- list(rect=list(left=par("usr")[2]+12*par("cxy")[1]))
    }
    if(fill.type != "none" && nlevels(class.fill) > 0) {
        lgd <- legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
                      pch = my.pch, pt.lwd = my.pt.lwd, cex = 1, pt.cex = my.pt.cex,
                      col = if(line.type=="none") "black" else "white", pt.bg = fillColorPalette,
                      title = fill.type, legend = levels(class.fill))
    }
    if(line.type != "none" && nlevels(class.line) > 0) {
        legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
               pch = my.pch, pt.lwd = my.pt.lwd, cex = 1, pt.cex = my.pt.cex,
               col = lineColorPalette, pt.bg = "white",
               title = line.type, legend = levels(class.line))
    }
    
    # close output file
    if(!is.na(image.format))
        dev.off()
    
    return(invisible(res))
}
