# checkColorType
# check if fill.type, line.type, mark.type argument is valid
#     val: type to check
checkColorType <- function(val, ncells) {
    if (!((is.character(val) && length(val) == 1 && val %in% c("none","predicted","true")) ||
          (is.numeric(val) && length(val) == ncells) ||
          (is.factor(val) && length(val) == ncells))) {
            stop(deparse(substitute(val)),
                 ' must be either a numeric(', ncells, '), ',
                 'a factor(', ncells, ') or one of "none", "predicted" or "true"')
    }
    return(TRUE)
}

# getPlotColors
# depending on val (fill.type, line.type, mark.type) and col (fill.col, line.col, mark.col),
# generate a vector of colors (palette and type as attributes)
#     gr: griph result
#     val: what to map to colors, one of:
#           "none": no color
#           "predicted": gr$MEMB (or V(graph)$membership)
#           "true": gr$MEMB.true (or V(GRAOp)$class)
#           numeric vector: linearly map to colors
#           factor: map each level to a different color
#     col: defines the color palette, either the name of a RColorBrewer palette, or a vector of R colors
#     graph: if non NULL, take 
#     NA.color: use for is.numeric(val) & is.na(val)
getPlotColors <- function(gr, val, col = "Spectral", graph = NULL, NA.color = "lightgray") {
    # digest arguments
    stopifnot(is.list(gr),
              all(c("MEMB","MEMB.true") %in% names(gr)),
              (is.character(val) && length(val) == 1 && val %in% c("none","predicted","true")) || is.numeric(val) || is.factor(val),
              is.character(col) && (length(col) > 1 || col %in% rownames(RColorBrewer::brewer.pal.info)),
              is.null(graph) || inherits(graph, "igraph"))
    ncells <- if (is.null(graph)) length(gr$MEMB) else length(igraph::V(graph))
    # map to colors
    coltype <- ""
    if ((is.character(val) && length(val) == 1 && val %in% c("none","predicted","true")) || is.factor(val)) {
        if (is.factor(val)) {
            if (is.null(graph)) {
            class.col <- val
            }
            else {class.col <- val[igraph:::V(gr$GRAO)$labels %in% igraph:::V(graph)$labels ]}
        } else {
            coltype <- val
            if (is.null(graph)) { # use gr
                class.col <- switch(val,
                                    predicted = factor(gr$MEMB, levels = sort(as.numeric(unique(gr$MEMB)))),
                                    true = factor(gr$MEMB.true, levels = unique(gr$MEMB.true)),
                                    none = factor(rep(NA, ncells)))
            } else {#               use graph
                class.col <- switch(val,
                                    predicted = factor(igraph::V(graph)$membership, levels = sort(as.numeric(unique(igraph::V(graph)$membership)))),
                                    true = factor(igraph::V(graph)$class, levels = unique(igraph::V(graph)$class)),
                                    none = factor(rep(NA, ncells)))
            }
        }
        nclass <- nlevels(class.col)
        if (length(col) == 1)
            col <- RColorBrewer::brewer.pal(min(RColorBrewer::brewer.pal.info[col, "maxcolors"], max(3, nclass)), col)[1:min(RColorBrewer::brewer.pal.info[col, "maxcolors"], nclass)]
        colpalette <- if (nclass > length(col)) grDevices::colorRampPalette(col)(nclass) else col
        names(colpalette) <- levels(class.col)
        cols <- colpalette[as.numeric(class.col)]
        names(cols) <- as.character(class.col)

    } else if (is.numeric(val)) {
        i <- !is.na(val)
        rng <- range(val[i])
        colpalette <- if (length(col) == 1) RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[col, "maxcolors"], col) else col
        names(colpalette) <- as.character(signif(seq(rng[1], rng[2], length.out = length(colpalette)), 3))
        cols <- rep(NA.color, length(val))
        cr <- colorRamp(colpalette)((val[i] - rng[1]) / (rng[2] - rng[1]))
        cols[i] <- grDevices::rgb(cr[,1], cr[,2], cr[,3], maxColorValue = 255)

    }
    
    attr(x = cols, which = "type") <- coltype
    attr(x = cols, which = "palette") <- colpalette
    return(cols)    
}

# drawLegends
# draw one or several legends (internal use only), depending on fill.X, line.X, mark.X
# and corresponding colors, consistently for all plotXXX functions
drawLegends <- function(mark.type, markColor, markColorType, markColorPalette,
                        fill.type, fillColor, fillColorType, fillColorPalette,
                        line.type, lineColor, lineColorType, lineColorPalette,
                        my.pch, my.pt.lwd, my.pt.cex) {
    
    # draw the legends
    # ... for mark
    if (mark.type != "none") {
        lgd <- legend(x = par("usr")[2] + 12 * par("cxy")[1], y = par("usr")[4],
                      xjust = 1, yjust = 1, bty = "n", cex = 1,
                      fill = paste0(markColorPalette, "66"), title = markColorType,
                      legend = names(markColorPalette))

    } else {
        lgd <- list(rect = list(left = par("usr")[2] + 12 * par("cxy")[1]))
    }
    # ... for fill
    if ((length(fill.type)==1 && fill.type != "none") || length(fill.type) > 1) {
        if (is.numeric(fill.type)) {

            # draw a continuous color legend
            llevs <- as.numeric(names(fillColorPalette))
            lcols <- fillColorPalette
            zlabs <- as.character(signif(seq(llevs[1], llevs[length(llevs)], length.out = 3), 2))
            pusr <- par('usr')
            pcxy <- par("cxy")
            xl <- lgd$rect$left - 2.0 * pcxy[1] - max(strwidth(zlabs, cex = 0.75))
            xr <- xl + 1.0 * pcxy[1]
            yb <- pusr[4] - 4.5 * pcxy[2]
            yt <- yb + 3.0 * pcxy[2]
            rect(xleft  = xl, ybottom = seq(yb, yt, length = 65)[-65],
                 xright = xr, ytop    = seq(yb, yt, length = 65)[-1],
                 col = grDevices::colorRampPalette(lcols)(64), border = NA)
            segments(x0 = xr, y0 = seq(yb, yt, length.out = 3),
                     x1 = xr + 0.3 * pcxy[1], y1 = seq(yb, yt, length.out = 3), col = "black")
            text(x = xr + 0.5 * pcxy[1], y = seq(yb, yt, length.out = 3),
                 labels = zlabs, adj = c(0, 0.5), srt = 0, cex = 0.75)
            lgd <- list(rect = list(left = xl - 0.5 * pcxy[1]))
        } else {

            # draw a discrete color legend
            lgd <- legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
                          pch = my.pch, pt.lwd = my.pt.lwd, cex = 1, pt.cex = my.pt.cex,
                          col = if (line.type[1] == "none") "black" else "white",
                          pt.bg = fillColorPalette,
                          title = fillColorType, legend = names(fillColorPalette))
        }
    }
    # ... for line
    if ((length(line.type)==1 && line.type != "none") || length(line.type) > 1) {

        legend(x = lgd$rect$left, y = par("usr")[4], xjust = 1, yjust = 1, bty = "n",
               pch = my.pch, pt.lwd = my.pt.lwd, cex = 1, pt.cex = my.pt.cex,
               col = lineColorPalette, pt.bg = "white",
               title = lineColorType, legend = names(lineColorPalette))
    }
}

# drawPolygons
# draw polygons (mark.type != "none")
#    e: ncell-by-2 matrix with x and y coordinates of 2D-embedding
#    markElements: list of indices for cells in each group
#    markColorPalette: colors for each group
drawPolygons <- function(e, markElements, markColorPalette) {
    for (j in which(lengths(markElements) > 0)) {
        xy <- e[markElements[[j]], , drop = FALSE]
        off <- par("cxy")[2] * 1
        pp <- rbind(xy,
                    cbind(xy[, 1] - off, xy[, 2]),
                    cbind(xy[, 1] + off, xy[, 2]),
                    cbind(xy[, 1], xy[, 2] - off), 
                    cbind(xy[, 1], xy[, 2] + off))
        cl <- igraph::convex_hull(pp)
        graphics::xspline(cl$rescoords, shape = 0.5, open = FALSE,
                          col = paste0(markColorPalette[j], "66"),
                          border = adjust.color(markColorPalette[j], 0.5))
    }
}

#' @title Visualize griph result as a graph.
#' 
#' @description Plot a graph obtained from \code{\link{SC_cluster}} or
#'     \code{\link{griph_cluster}}, allowing to apply graph subsampling and
#'     control coloring.
#' 
#' @param gr A \code{griph} result, as returned by \code{\link{SC_cluster}} or
#'     \code{\link{griph_cluster}}.
#' @param maxG Approximate maximal number of vertices to include when plotting
#'     the graph.
#' @param fill.type Type of fill color, one of \code{"predicted"} (predicted class
#'     labels, default), \code{"true"} (true class labels, if available), \code{"none"}
#'     (no fill color), a numeric vectors  or a factor whose values should be
#'     mapped to fill colors, using the palette given in \code{fill.col}.
#' @param line.type Type of line color, as in \code{fill.type}.
#' @param mark.type Type of cell class defnition to mark using polygons, as in \code{fill.type}.
#' @param collapse.type Type of cell class to use for graph simplification,
#'     by combining cells of the same class into a single vertex. If set to a value
#'     other than \code{"none"}, the same value will also be used for \code{fill.type}
#'     and \code{line.type} is ignored.
#' @param fill.col Fill color palette, either a character scalar with a \code{\link{RColorBrewer}}
#'     color palette name or color vector defining the palette to use for vertex fill coloring.
#' @param line.col Line color palette, as in \code{fill.col}.
#' @param mark.col Mark color palette, as in \code{fill.col}.
#' @param draw.edges If \code{NULL} (default), draw edges if \code{collapse.type != "none"}.
#'     \code{TRUE} or \code{FALSE} can be used to override the default.
#' @param seed Random number seed to make graph layout deterministic.
#' @param fsuffix A suffix added to the file names of output plots. If \code{NULL} (default),
#'     it will use a random 5 character string. Ignored if \code{image.format} is \code{NULL}.
#' @param image.format Specifies the format of the created image. Currently supported are
#'     \code{\link{pdf}}, \code{\link{png}} or \code{NA}. If \code{NA} (the default), the plot
#'     is rendered on the currently opened plotting device.
#' @param forceRecalculation If \code{TRUE}, recalculate plotting-optimized graph
#'     even if it is already contained in \code{gr}.
#' @param quiet If \code{TRUE}, do not report on progress.
#' @param plot.args a named list of  arguments for controlling specific plot points attributes (cex, lwd, pch -should between 21 and 25- )
#' 
#' @return Invisibly the plot-optimized version of the graph as an \code{igraph} object.
plotGraph <- function(gr, maxG=2500,
                      fill.type="predicted",
                      line.type="none",
                      mark.type="none",
                      collapse.type=c("none","predicted","true"),
                      fill.col="Spectral",
                      line.col="Dark2",
                      mark.col="Pastel1",
                      draw.edges=NULL,
                      seed=91919,
                      fsuffix=NULL, image.format=NA,
                      forceRecalculation=FALSE, quiet=FALSE,
                      plot.args=list(pch=21L, cex=1.0,lwd=2.5)
                      ) {
    if (!quiet)
        message("Computing Graph Layout and Rendering...")
    
    # digest arguments
    GRAO <- gr$GRAO
    MEMB <- gr$MEMB
    csize <- table(MEMB)
    checkColorType(fill.type, length(MEMB))
    checkColorType(line.type, length(MEMB))
    checkColorType(mark.type, length(MEMB))
    collapse.type <- match.arg(collapse.type)

    
    # global plotting paramterers (can come as arguments)
    plot.args$lwd <- if (line.type[1] == "none" || collapse.type[1] != "none") 0.1 else plot.args$lwd
    
    pct <- 1
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
                                 predicted = factor(igraph::V(GRAO)$membership,
                                                    levels = sort(as.numeric(unique(igraph::V(GRAO)$membership)))),
                                 true = factor(igraph::V(GRAO)$class,
                                               levels = unique(igraph::V(GRAO)$class)))
        GRAOp <- igraph::simplify( igraph::contract(GRAO, class.collapse) )
        igraph::V(GRAOp)$membership <- igraph::V(GRAOp)$class <- levels(class.collapse)
        igraph::V(GRAOp)$size <- 10 / (length(igraph::V(GRAOp)) / 60)^0.3 * (as.numeric(csize / median(csize)))^0.5
        
    } else if (is.null(gr$plotGRAO) || forceRecalculation) {
        if (length(igraph::V(GRAO)) > 1.25 * maxG ) {
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
        if (median(igraph::degree(GRAOp)) > 4) {
            pct <- min(1, 1 / sqrt(0.1 * median(igraph::degree(GRAOp))))
            ADJp <- igraph::as_adj(GRAOp, attr = 'weight', sparse = TRUE)
            ADJp <- sparsify(ADJp, pct)
            ADJp[which(abs(ADJp) > 0)] <- 1
            GRAOtemp <- igraph::graph.adjacency(ADJp, mode = "max", weighted = NULL,
                                                diag = FALSE)
            GRAOp <- igraph::intersection(GRAOp, GRAOtemp, byname = FALSE)
            GRAOtemp <- NULL
            ADJp <- NULL
        }
        
        GRAOp <- igraph::delete_vertices(GRAOp, which(igraph::ego_size(GRAOp, 3) < 6))
        ###Delete Vertices from communites with few members:
        min.csize <- ceiling(0.25 * sqrt(length(igraph::V(GRAO))))
        GRAOp <- igraph::delete_vertices(GRAOp, which(igraph::V(GRAOp)$community.size < min.csize))  
        
        if (!quiet)
            message("\tnodes from communities with <",min.csize, " members will not be displayed.")
        
    } else {
        if (!quiet)
            message("\tusing existing plot-optimized graph")
        GRAOp <- gr$plotGRAO$GRAO
    }
    if (!quiet)
        message("\tdisplaying graph with ",length(igraph::V(GRAOp))," (",
                round(100 * length(igraph::V(GRAOp)) / length(igraph::V(GRAO)), 1), "%) vertices and ",
                length(igraph::E(GRAOp)), " (", round(100 * length(igraph::E(GRAOp)) / length(igraph::E(GRAO)), 1),
                "%) edges")
    
    # get colors
    fillColor <- getPlotColors(gr = gr, val = fill.type, col = fill.col, graph = GRAOp)
    fillColorPalette <- attr(x = fillColor, which = "palette")
    fillColorType <- attr(x = fillColor, which = "type")
    
    lineColor <- getPlotColors(gr = gr, val = line.type, col = line.col, graph = GRAOp)
    lineColorPalette <- attr(x = lineColor, which = "palette")
    lineColorType <- attr(x = lineColor, which = "type")
    
    markColor <- getPlotColors(gr = gr, val = mark.type, col = mark.col, graph = GRAOp)
    markColorPalette <- attr(x = markColor, which = "palette")
    markColorType <- attr(x = markColor, which = "type")
    markElements <- split(seq_along(markColor), names(markColor))[names(markColorPalette)]

    # set some more graph attributes
    igraph::V(GRAOp)$classcolor <- lineColor
    igraph::V(GRAOp)$size <- if (collapse.type == "none") 10 / (length(igraph::V(GRAOp)) / 60)^0.3 else igraph::V(GRAOp)$size
    igraph::V(GRAOp)$cex <- igraph::V(GRAOp)$size / 3
    igraph::V(GRAOp)$frame.width <- 2 / (length(igraph::V(GRAOp)) / 6 )^0.3
    igraph::E(GRAOp)$width <- igraph::E(GRAOp)$weight / sqrt((length(igraph::V(GRAOp))/60 ))
    igraph::V(GRAOp)$color <- fillColor
    # compute graph layout
    set.seed(seed = seed)

    l <- igraph::layout_with_fr(graph = GRAOp)
    # igraph::add.vertex.shape("fcircle", clip=igraph.shape.noclip, plot=mycircle, parameters=list(vertex.frame.color=1, vertex.frame.width=1))

    # open output file
    if (!is.na(image.format)) {
        if (is.null(fsuffix))
            fsuffix <- RandString()
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
    xmin=max(  quantile(l[,1],0.1)-0.5*IQR(l[,1]), min(l[,1])  )
    xmax=min(  quantile(l[,1],0.9)+IQR(l[,1]), max(l[,1])  )
    ymin=max(  quantile(l[,2],0.1)-0.5*IQR(l[,2]), min(l[,2])  )
    ymax=min(  quantile(l[,2],0.9)+IQR(l[,2]), max(l[,2])  )
    plot(l[,1], l[,2], type = "n", axes = FALSE, xlab = "", ylab = "", xlim=c(xmin,xmax),ylim=c(ymin,ymax)    )
    
    
    
    # add mark polygons
    if (mark.type != "none")
        drawPolygons(e = l, markElements = markElements, markColorPalette = markColorPalette)
    
    
    # add edges (by default only for collapase.type != "none")
    if (isTRUE(draw.edges) || (is.null(draw.edges) && collapse.type != "none")) {
        el <- igraph::as_edgelist(GRAOp, names = FALSE)
        graphics::segments(x0 = l[,1][el[,1]], y0 = l[,2][el[,1]],
                           x1 = l[,1][el[,2]], y1 = l[,2][el[,2]],
                           col = edge.col, lwd = igraph::E(GRAOp)$weight / max(igraph::E(GRAOp)$weight) * edge.lwd.max)
    }
    
    # add vertices
    points(l[,1], l[,2], col = if (line.type[1] == "none") "black" else lineColor,
           bg = fillColor, pch = plot.args$pch, lwd = plot.args$lwd, 
           cex =  plot.args$cex * if (collapse.type == "none") 1.0 else (as.numeric(csize / median(csize)))^0.5)
    
    # add MST segments if MST is present
    if ( !is.null(gr$MST)) {
        for (e in 1:length(E(gr$MST ))   ) {
            source=names(igraph::tail_of ( gr$MST, e) )
            target=names(igraph::head_of ( gr$MST, e))
            if (exists("snowball")){
            source.nodes=which(sort(snowball) %in% which(gr$MST.memb==source) )
            target.nodes=which(sort(snowball) %in% which(gr$MST.memb==target) )
            }else{
            source.nodes=which(gr$MST.memb==source)
            target.nodes=which(gr$MST.memb==target)
            }
            X0=mean(l[source.nodes,1], trim=0.1)
            Y0=mean(l[source.nodes,2], trim=0.1)
            X1=mean(l[target.nodes,1], trim=0.1)
            Y1=mean(l[target.nodes,2], trim=0.1)
            segments (X0,Y0,X1,Y1,col="red",lwd=2  )
            points(c(X0,X1),c(Y0,Y1),col="red",bg="#555555",pch=21,cex=1.0)
        }
    }
    
    
    # add legend(s)
    drawLegends(mark.type, markColor, markColorType, markColorPalette,
                fill.type, fillColor, fillColorType, fillColorPalette,
                line.type, lineColor, lineColorType, lineColorPalette,
                plot.args$pch, plot.args$lwd, plot.args$cex)

    # close output file
    if (!is.na(image.format))
        dev.off()

    return(list(y = l[,1:2], GRAO = invisible(GRAOp)))

}





#' @title Visualize griph result as a t-SNE projection.
#' 
#' @description Plot a t-SNE projection of the affinity matrix obtained from
#'     \code{\link{SC_cluster}} or \code{\link{griph_cluster}}, allowing to
#'     control coloring.
#' 
#' @param gr A \code{griph} result, as returned by \code{\link{SC_cluster}} or
#'     \code{\link{griph_cluster}}.
#' @param fill.type Type of fill color, one of \code{"predicted"} (predicted class
#'     labels, default), \code{"true"} (true class labels, if available), \code{"none"}
#'     (no fill color), a numeric vectors  or a factor whose values should be
#'     mapped to fill colors, using the palette given in \code{fill.col}.
#' @param line.type Type of line color, as in \code{fill.type}.
#' @param mark.type Type of cell class defnition to mark using polygons, as in \code{fill.type}.
#' @param fill.col Fill color palette, either a character scalar with a \code{\link{RColorBrewer}}
#'     color palette name or color vector defining the palette to use for vertex fill coloring.
#' @param line.col Line color palette, as in \code{fill.col}.
#' @param mark.col Mark color palette, as in \code{fill.col}.
#' @param seed Random number seed to make t-SNE projection deterministic.
#' @param fsuffix A suffix added to the file names of output plots. If \code{NULL} (default),
#'     it will use a random 5 character string. Ignored if \code{image.format} is \code{NULL}.
#' @param image.format Specifies the format of the created image. Currently supported are
#'     \code{\link{pdf}}, \code{\link{png}} or \code{NA}. If \code{NA} (the default), the plot
#'     is rendered on the currently opened plotting device.
#' @param forceRecalculation If \code{TRUE}, recalculate t-SNE 2D-embedding
#'     even if it is already contained in \code{gr}.
#' @param quiet If \code{TRUE}, do not report on progress.
#' @param plot.args a named list of  arguments for controlling specific plot points attributes (cex, lwd, pch -should between 21 and 25- )
#' @param ... additional arguments passed to \code{Rtsne} 
#' 
#' @return Invisible the results of the t-SNE projection.
#' 
#' @seealso \code{Rtsne} used to calculate the t-SNE projection.
plotTsne <- function(gr,
                     fill.type="predicted",
                     line.type="none",
                     mark.type="none",
                     fill.col="Spectral",
                     line.col="Dark2",
                     mark.col="Pastel1",
                     seed=91919,
                     fsuffix=NULL, image.format=NA,
                     forceRecalculation=FALSE,
                     plot.args=list(pch=21L, cex=1.0,lwd=2.5),
                     quiet=FALSE, ...) {
    if (!is.element("Rtsne", utils::installed.packages()[,1])) {
        stop('"plotTsne" requires the "Rtsne" package. Please install it with:\n\t',
             'install.packages("Rtsne")')
    }
    
    # digest arguments
    add.args <- list(...)
    MEMB <- gr$MEMB
    checkColorType(fill.type, length(MEMB))
    checkColorType(line.type, length(MEMB))
    checkColorType(mark.type, length(MEMB))
    
    # global plotting paramterers (can come as arguments)
    plot.args$lwd <- if (line.type[1] == "none") 0.1 else plot.args$lwd
    

    # get t-SNE projection
    if (!forceRecalculation && "plotTsne" %in% names(gr) && !is.null(gr$plotTsne)) {
        if (!quiet)
            message("Using existing t-SNE embedding")
        res <- gr$plotTsne
    } else {
        if (!quiet)
            message("Computing t-SNE embedding...")
        
        set.seed(seed = seed)
        if (!("perplexity" %in% names(add.args)))
            add.args$perplexity <- min(30, round(sqrt(nrow(gr$DISTM)) - 1))
        res <- do.call(Rtsne::Rtsne, c(list(X = stats::as.dist(1 - gr$DISTM), pca = FALSE, is_distance = TRUE), add.args))
    }
    
    # get colors
    fillColor <- getPlotColors(gr = gr, val = fill.type, col = fill.col)
    fillColorPalette <- attr(x = fillColor, which = "palette")
    fillColorType <- attr(x = fillColor, which = "type")
    
    lineColor <- getPlotColors(gr = gr, val = line.type, col = line.col)
    lineColorPalette <- attr(x = lineColor, which = "palette")
    lineColorType <- attr(x = lineColor, which = "type")
    
    markColor <- getPlotColors(gr = gr, val = mark.type, col = mark.col)
    markColorPalette <- attr(x = markColor, which = "palette")
    markColorType <- attr(x = markColor, which = "type")
    markElements <- split(seq_along(markColor), names(markColor))[names(markColorPalette)]

    # open output file
    if (!is.na(image.format)) {
        if (is.null(fsuffix))
            fsuffix <- RandString()
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
    if (mark.type != "none")
        drawPolygons(e = res$Y, markElements = markElements, markColorPalette = markColorPalette)
    
    # add cells
    points(res$Y, col = if (line.type[1] == "none") "black" else lineColor,
           bg = fillColor, pch = plot.args$pch, lwd = plot.args$lwd, cex = plot.args$cex)
    
    # add MST segments if MST is present
    if ( !is.null(gr$MST)) {
        for (e in 1:length(E(gr$MST ))   ) {
            source=names(igraph::tail_of ( gr$MST, e) )
            target=names(igraph::head_of ( gr$MST, e))
            source.nodes=which(gr$MST.memb==source)
            target.nodes=which(gr$MST.memb==target)
            X0=mean(res$Y[source.nodes,1], trim=0.1)
            Y0=mean(res$Y[source.nodes,2], trim=0.1)
            X1=mean(res$Y[target.nodes,1], trim=0.1)
            Y1=mean(res$Y[target.nodes,2], trim=0.1)
            segments (X0,Y0,X1,Y1,col="red",lwd=2  )
            points(c(X0,X1),c(Y0,Y1),col="red",bg="#555555",pch=21,cex=1.0)
        }
    }
    
    
    
    
    # add legend(s)
    drawLegends(mark.type, markColor, markColorType, markColorPalette,
                fill.type, fillColor, fillColorType, fillColorPalette,
                line.type, lineColor, lineColorType, lineColorPalette,
                plot.args$pch, plot.args$lwd, plot.args$cex)
    
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
#' @param fill.type Type of fill color, one of \code{"predicted"} (predicted class
#'     labels, default), \code{"true"} (true class labels, if available), \code{"none"}
#'     (no fill color), a numeric vectors  or a factor whose values should be
#'     mapped to fill colors, using the palette given in \code{fill.col}.
#' @param line.type Type of line color, as in \code{fill.type}.
#' @param mark.type Type of cell class defnition to mark using polygons, as in \code{fill.type}.
#'     \code{custom} (polygons around cells with the same \code{custom.class}).
#' @param fill.col Fill color palette, either a character scalar with a \code{\link{RColorBrewer}}
#'     color palette name or color vector defining the palette to use for vertex fill coloring.
#' @param line.col Line color palette, as in \code{fill.col}.
#' @param mark.col Mark color palette, as in \code{fill.col}.
#' @param seed Random number seed to make largeVis projection deterministic.
#' @param fsuffix A suffix added to the file names of output plots. If \code{NULL} (default),
#'     it will use a random 5 character string. Ignored if \code{image.format} is \code{NULL}.
#' @param image.format Specifies the format of the created image. Currently supported are
#'     \code{\link{pdf}}, \code{\link{png}} or \code{NA}. If \code{NA} (the default), the plot
#'     is rendered on the currently opened plotting device.
#' @param forceRecalculation If \code{TRUE}, recalculate LargeVis 2D-embedding
#'     even if it is already contained in \code{gr}.
#' @param quiet If \code{TRUE}, do not report on progress.
#' @param plot.args a named list of  arguments for controlling specific plot points attributes (cex, lwd, pch -should between 21 and 25- )
#' @param ... additional arguments passed to \code{largeVis::projectKNNs} 
#' 
#' @return Invisible the results of the largeVis projection as a two-column matrix.
#' 
#' @seealso \code{largeVis} used to calculate the largeVis projection.
plotLVis <- function(gr,
                     fill.type="predicted",
                     line.type="none",
                     mark.type="none",
                     fill.col="Spectral",
                     line.col="Dark2",
                     mark.col="Pastel1",
                     seed=91919,
                     fsuffix=NULL, image.format=NA,
                     forceRecalculation=FALSE,
                     quiet=FALSE, plot.args=list(pch=21L, cex=1.0,lwd=2.5),
                     use.par = TRUE, ncores = "all",
                     ...) {
    # digest arguments
    add.args <- list(...)
    MEMB <- gr$MEMB
    checkColorType(fill.type, length(MEMB))
    checkColorType(line.type, length(MEMB))
    checkColorType(mark.type, length(MEMB))
    
    # global plotting paramterers (can come as arguments)
    plot.args$lwd <- if (length(line.type)==1 && line.type == "none") 0.1 else plot.args$lwd

    
    
    
    # get largeVis projection
    if (!forceRecalculation && "plotLVis" %in% names(gr) && !is.null(gr$plotLVis)) {
        if (!quiet)
            message("Using existing LargeVis embedding")
        res <- t(gr$plotLVis)
    } else {
        
        
        
        
        # Register cluster here if one not already registered.
        do.register=FALSE
        if (isTRUE(use.par) & !foreach::getDoParRegistered() ) {
            do.register=TRUE
            if (ncores == "all") {
                ncores <- parallel::detectCores()
            } else {
                ncores <- min(ncores, parallel::detectCores() )
            }
            cl <- parallel::makeCluster(ncores)
            doParallel::registerDoParallel(cl)
        }
        ### wrap code in tryCatch block, ensuring that stopCluster(cl) is called even when a condition is raised  
        tryCatch({    
        
        
            
            
        if (!quiet)
            message("Computing LargeVis embedding...")
    
        set.seed(seed = seed)
    
       if (!is.element('sgd_batches', names(add.args)))
            add.args$sgd_batches <- max(0.1 * sgdBatches(ncol(gr$DISTM), Matrix::nnzero(gr$DISTM)), 1e7)  #!!!!!!Use GRAO instead!!! 
        if (!is.element('M',names(add.args)))
            add.args$M <- 2    
        if (!is.element('gamma', names(add.args)))
            add.args$gamma <- 20    
        if (!is.element('alpha', names(add.args)))
            add.args$alpha <- 0.2
        if (!is.element('useDegree', names(add.args)))
            add.args$useDegree <- TRUE    

        #res <- do.call(projectKNNs, c( list(wij=igraph::as_adj(gr$GRAO,names=FALSE, attr = 'weight', sparse=TRUE),seed=seed),add.args )    )
        res <- do.call(projectKNNs, c(list(wij = gr$DISTM, seed = seed, threads=foreach::getDoParWorkers()   ), add.args))
        
        
        }, # end of tryCatch expression, cluster object cl not needed anymore    
       finally = { 
           ##### Stop registered cluster and unregister if initiated infunction (tracked by do.register):
           if (isTRUE(use.par) & foreach::getDoParRegistered() & do.register==TRUE)
               parallel::stopCluster(cl)
               env <- foreach:::.foreachGlobals
               rm(list=ls(name=env), pos=env)
       })    
        
    }
    

    # get colors
    fillColor <- getPlotColors(gr = gr, val = fill.type, col = fill.col)
    fillColorPalette <- attr(x = fillColor, which = "palette")
    fillColorType <- attr(x = fillColor, which = "type")
    
    lineColor <- getPlotColors(gr = gr, val = line.type, col = line.col)
    lineColorPalette <- attr(x = lineColor, which = "palette")
    lineColorType <- attr(x = lineColor, which = "type")
    
    markColor <- getPlotColors(gr = gr, val = mark.type, col = mark.col)
    markColorPalette <- attr(x = markColor, which = "palette")
    markColorType <- attr(x = markColor, which = "type")
    markElements <- split(seq_along(markColor), names(markColor))[names(markColorPalette)]

    # open output file
    if (!is.na(image.format)) {
        if (is.null(fsuffix))
            fsuffix <- RandString()
        if (image.format == 'pdf') {
            fname <- paste('Lvis_', fsuffix, '.pdf', sep = "")
            pdf(file = fname, width = 12, height = 10)
        } else if (image.format == "png") {
            fname <- paste('Lvis_', fsuffix, '.png', sep = "")
            png(filename = fname, width = 12, height = 10, units = "in", res = 300)   
        }
        if (!quiet)
            message("\tsaving LargeVis plot to ", fname)
    }
    
    # setup plot coordinate system
    par(mar = c(5.1, 4.1, 4.1, 14.1), xpd = TRUE)
    plot(t(res)[,1], t(res)[,2], type = "n", axes = FALSE, xlab = "", ylab = "")
    
    # add mark polygons
    if (mark.type != "none")
        drawPolygons(e = t(res), markElements = markElements, markColorPalette = markColorPalette)
    
    # add cells
    points(t(res), col = if (length(line.type)==1 && line.type == "none") "black" else lineColor,
           bg = fillColor, pch = plot.args$pch, lwd = plot.args$lwd, cex = plot.args$cex)
    
    
    # add MST segments if MST is present
    if ( !is.null(gr$MST)) {
    for (e in 1:length(E(gr$MST ))   ) {
        source=names(igraph::tail_of ( gr$MST, e) )
        target=names(igraph::head_of ( gr$MST, e))
        source.nodes=which(gr$MST.memb==source)
        target.nodes=which(gr$MST.memb==target)
        X0=mean(t(res)[source.nodes,1], trim=0.1)
        Y0=mean(t(res)[source.nodes,2], trim=0.1)
        X1=mean(t(res)[target.nodes,1], trim=0.1)
        Y1=mean(t(res)[target.nodes,2], trim=0.1)
        segments (X0,Y0,X1,Y1,col="red",lwd=2  )
        points(c(X0,X1),c(Y0,Y1),col="red",bg="#555555",pch=21,cex=1.0)
    }
    }
    
    # add legend(s)
    drawLegends(mark.type, markColor, markColorType, markColorPalette,
                fill.type, fillColor, fillColorType, fillColorPalette,
                line.type, lineColor, lineColorType, lineColorPalette,
                plot.args$pch, plot.args$lwd, plot.args$cex)
    
    # close output file
    if (!is.na(image.format))
        dev.off()
    
    rres <- t(res)
    dimnames(rres) <- list(names(gr$MEMB), c("x","y"))
    
    
    

    
    
    
    return(invisible(rres))
}
