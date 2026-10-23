.pm_network_plot_types <- function(x) {
  types <- c("network", if(!is.null(x$results$uva))"uva", "loadings", "loadings_matrix", "centrality", "tefi")
  if (!is.null(x$results$wording)) types <- c(types, "wording")
  if (!is.null(x$results$hierarchy)) types <- c(types, "hierarchy")
  if (!is.null(x$results$bootstrap)) types <- c(types, "ega_stability", "dimensions", "item_stability", "dimension_stability", "edge_accuracy", "centrality_stability")
  if (!is.null(x$results$invariance)) types <- c(types, "group_networks", "edge_differences", "loading_differences")
  types
}

#' Plot cached psychometric network results
#' @param x A `psychomatic_network` result.
#' @param type Figure type. See `network_plots()` for available figures.
#' @param language Figure language, defaults to the result language.
#' @param label_cex Node label size multiplier, default 1.2. Labels are not
#'   automatically shrunk to fit the nodes; use concise item identifiers.
#' @param ... Reserved for future graphical arguments.
#' @return Invisibly, a qgraph object for `"network"` and `"loadings"`, or the
#'   underlying results for other figures. No analysis is recomputed.
#' @export
plot.psychomatic_network <- function(x, type = "network", language = NULL, label_cex = 1.2, ...) {
  .pm_validate_scalar_number(label_cex,"label_cex",.3,4)
  language <- .pm_normalize_language(language, x$language)
  es <- language == "esp"
  tr <- function(esp, eng) if (es) esp else eng
  type <- match.arg(type, .pm_network_plot_types(x))
  r <- x$results
  old <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(old), add = TRUE)
  draw_network <- function(a, wc, labels, title, layout = "spring", maximum = NULL,
      edge_labels = FALSE, show_legend = TRUE, node_status = NULL) {
    .pm_network_require("qgraph")
    if (!length(wc) || all(is.na(wc))) {
      graphics::plot.new(); graphics::title(main=title)
      graphics::text(.5,.5,tr("Sin comunidades estimables\nTEFI y estabilidad no evaluables",
        "No estimable communities\nTEFI and stability not evaluable"),cex=.9)
      return(invisible(NULL))
    }
    groups <- split(seq_along(wc), wc)
    names(groups) <- if(prod(graphics::par('mfrow'))>1L) paste0('C',names(groups)) else paste(tr('Comunidad','Community'),names(groups))
    args <- list(input = a, layout = layout, groups = groups,
      labels = labels, title = title, theme = "colorblind", legend = show_legend,
      legend.mode='groups',legend.cex=.7,GLratio=4,
      color = unname(.pm_network_color_fill(.pm_network_community_colors(wc))),
      border.color=unname(.pm_network_community_colors(wc)[as.character(wc)]),
      posCol='#596773',negCol='#B44E35',negDashed=TRUE,label.color='#000000',
      vsize = 7, label.cex = label_cex, label.scale = FALSE,
      border.width = 1.5, mar = c(5,2,4,2),
      edge.labels = edge_labels, edge.label.cex = .75, edge.label.bg = 'white')
    if (!is.null(maximum)) args$maximum <- maximum
    if(!is.null(node_status)) {
      args$border.color <- ifelse(node_status=='difference','#B34D00','#666666')
      args$border.width <- ifelse(node_status=='difference',4,1.5)
      args$shape <- ifelse(node_status=='unavailable','square','circle')
    }
    graph <- do.call(qgraph::qgraph, args)
    if(prod(graphics::par('mfrow'))==1L) graphics::mtext(tr('Aristas: positivas continuas; negativas discontinuas. Grosor = |peso|.',
      'Edges: positive solid; negative dashed. Width = |weight|.'),side=1,line=3,cex=.65)
    invisible(graph)
  }
  bars <- function(values, labels, title, xlab, limits = NULL, colors = '#4477AA') {
    graphics::par(mar = c(5, max(5, min(12, max(nchar(labels))*.45)),4,2))
    args <- list(height = values, names.arg = labels, horiz = TRUE, las = 1,
      col = colors, border = NA, main = title, xlab = xlab, cex.names = .85)
    if (!is.null(limits)) args$xlim <- limits
    do.call(graphics::barplot, args)
  }
  if (type == "network") {
    return(invisible(draw_network(r$ega$network, r$ega$wc, rownames(r$ega$network),
      tr("Red EGA", "EGA network"),
      edge_labels=sum(r$ega$network[upper.tri(r$ega$network)]!=0)<=20L)))
  } else if(type == 'ega_stability') {
    graphics::layout(matrix(c(1,2,3,3),2,byrow=TRUE),heights=c(1,.16))
    draw_network(r$ega$network,r$ega$wc,rownames(r$ega$network),
      tr('Dimensionalidad | EGA','Dimensionality | EGA'),show_legend=FALSE)
    t <- r$bootstrap$structural$items
    .pm_network_stability_points(t$stability,t$item,.pm_network_item_membership(x,t$item),
      tr('Estabilidad de \u00edtems','Item stability'),es,legend=FALSE)
    graphics::par(mar=c(0,0,0,0));graphics::plot.new()
    palette <- .pm_network_community_colors(r$ega$wc)
    graphics::legend(.5,.8,xjust=.5,yjust=.5,bty='n',horiz=TRUE,
      legend=paste(tr('Comunidad','Community'),names(palette)),pch=21,
      pt.bg=.pm_network_color_fill(palette),col=palette,pt.cex=1.8,cex=.85)
    graphics::text(.5,.15,paste0(tr('Aristas = ','Edges = '),.pm_network_edge_kind(x,es),tr(
      '; replicaci\u00f3n = proporci\u00f3n de asignaciones reproducidas.',
      '; replication = proportion of reproduced assignments.')),cex=.7)
  } else if (type == "uva") {
    m <- r$uva$wto$matrix
    graphics::image(seq_len(nrow(m)), seq_len(ncol(m)), m, axes = FALSE,
      col = grDevices::hcl.colors(40, "YlOrRd"), xlab = "", ylab = "",
      main = tr("Redundancia: solapamiento topol\u00f3gico ponderado", "Redundancy: weighted topological overlap"))
    graphics::axis(1,seq_len(nrow(m)),rownames(m),las=2,cex.axis=.6)
    graphics::axis(2,seq_len(ncol(m)),colnames(m),las=2,cex.axis=.6)
  } else if (type == "loadings") {
    .pm_network_require('qgraph')
    m <- r$loadings$matrix; primary <- r$loadings$primary
    p <- nrow(m); d <- ncol(m)
    # Separate positions distinguish observed items from dimensions. Only
    # assigned primary loadings are drawn; all cross-loadings remain in the map.
    membership <- match(primary$dimension[match(rownames(m),primary$item)],colnames(m))
    if(anyNA(membership))stop('Primary loading dimensions are not assigned.',call.=FALSE)
    a <- matrix(0,p+d,p+d)
    values <- m[cbind(seq_len(p),membership)]
    values[!is.finite(values)] <- 0
    a[cbind(seq_len(p),p+membership)] <- values
    a <- a+t(a)
    order <- order(membership,seq_len(p))
    y <- numeric(p);y[order] <- seq(1,-1,length.out=p)
    dim_y <- vapply(seq_len(d),function(j)mean(y[membership==j]),double(1))
    layout <- rbind(cbind(rep(-1,p),y),cbind(rep(1,d),dim_y))
    return(invisible(qgraph::qgraph(a,layout=layout,
      labels=c(rownames(m),paste(tr('Dim','Dim'),colnames(m))),
      groups=split(seq_len(p+d),c(membership,seq_len(d))),theme='colorblind',
      color=unname(.pm_network_color_fill(.pm_network_community_colors(seq_len(d)))),
      border.color=unname(.pm_network_community_colors(seq_len(d))[as.character(c(membership,seq_len(d)))]),
      label.color='#000000',
      posCol='#4D4D4D',negCol='#D55E00',negDashed=TRUE,
      shape=c(rep('circle',p),rep('square',d)),vsize=c(rep(5,p),rep(7,d)),esize=7,
      label.cex=label_cex,label.scale=FALSE,legend=FALSE,edge.labels=TRUE,edge.label.cex=1,
      edge.label.position=.55,edge.label.bg='white',
      title=tr('Cargas primarias: \u00edtem\u2013dimensi\u00f3n','Primary loadings: item\u2013dimension'),
      mar=c(3,3,5,3))))
  } else if (type == "loadings_matrix") {
    m <- r$loadings$matrix
    limit <- max(abs(m[is.finite(m)]),.01)
    palette <- grDevices::colorRampPalette(c('#70224E','#F8FAFC','#244C73'))(101)
    breaks <- seq(-limit,limit,length.out=102)
    cells <- rep('#FFFFFF',length(m))
    valid <- is.finite(m)
    cells[valid] <- palette[findInterval(m[valid],breaks,all.inside=TRUE)]
    text_color <- matrix(.pm_network_contrast_text(cells),nrow=nrow(m))
    graphics::layout(matrix(c(1,2),nrow=1),widths=c(1,.10))
    graphics::par(mar = c(7,6,4,1),mgp=c(2.5,1,0))
    graphics::image(seq_len(ncol(m)),seq_len(nrow(m)),t(m),axes=FALSE,
      col=palette,breaks=breaks,
      xlab=tr("Dimensi\u00f3n","Dimension"),ylab="",main=tr("Cargas de red revisadas","Revised network loadings"),font.main=1)
    graphics::axis(1,seq_len(ncol(m)),colnames(m)); graphics::axis(2,seq_len(nrow(m)),rownames(m),las=1,cex.axis=.65)
    for (j in seq_len(ncol(m))) for(i in seq_len(nrow(m)))
      graphics::text(j,i,if(is.finite(m[i,j]))sprintf("%.2f",m[i,j])else "-",cex=.7,col=text_color[i,j])
    graphics::par(mar=c(7,0,4,4))
    ticks <- seq(-limit,limit,length.out=101)
    graphics::image(1,ticks,matrix(ticks,nrow=1),col=palette,breaks=breaks,
      axes=FALSE,xlab='',ylab='')
    ticks <- pretty(c(-limit,limit))
    graphics::axis(4,at=ticks[ticks>=-limit & ticks<=limit],las=1,cex.axis=.8)
    graphics::mtext(tr('Carga','Loading'),side=4,line=2.5,cex=.85)
  } else if (type == "centrality") {
    graphics::par(mfrow=c(1,3))
    for (metric in c("strength","betweenness","closeness")) {
      values <- r$centrality[[metric]]
      if (all(!is.finite(values))) { graphics::plot.new(); graphics::title(main=metric); graphics::text(.5,.5,tr("No evaluable","Not evaluable")) }
      else bars(values,r$centrality$item,metric,metric)
    }
  } else if (type == "tefi") {
    t <- r$tefi
    t <- t[is.finite(t$value),,drop=FALSE]
    if(!nrow(t)) {
      graphics::plot.new();graphics::title(main="TEFI / genTEFI")
      graphics::text(.5,.5,tr("No evaluable: no hay \u00edndices finitos",
        "Not evaluable: no finite indices"))
    } else {
      sets <- split(t,paste(t$comparison_set,t$index,sep=' / '))
      graphics::par(mfrow=c(ceiling(length(sets)/2),min(2,length(sets))))
      for (name in names(sets)) {
        panel <- sets[[name]]
        bars(panel$value,panel$model,name,panel$index[1])
      }
      graphics::mtext(tr('No comparar entre paneles: matrices o \u00edndices distintos.',
        'Do not compare across panels: different matrices or indices.'),side=1,line=3,cex=.7)
    }
  } else if (type == "wording") {
    graphics::par(mfrow=c(1,2))
    w <- r$wording
    draw_network(w$original$network,w$original$wc,rownames(w$original$network),"Raw EGA")
    draw_network(w$adjusted$EGA$network,w$adjusted$EGA$wc,rownames(w$adjusted$EGA$network),"riEGA")
  } else if (type == "hierarchy") {
    h <- r$hierarchy$native
    graphics::par(mfrow=c(1,2))
    draw_network(h$lower_order$network,h$lower_order$wc,rownames(h$lower_order$network),tr("Red de orden inferior","Lower-order network"))
    draw_network(h$higher_order$network,h$higher_order$wc,paste0('C',rownames(h$higher_order$network)),tr("Red de orden superior","Higher-order network"))
  } else if (type == "dimensions") {
    f <- r$bootstrap$structural$frequency
    bars(f[[2]],as.character(f[[1]]),tr("Dimensiones bootstrap","Bootstrap dimensionality"),tr("Frecuencia relativa","Relative frequency"),c(0,1))
  } else if (type == "item_stability") {
    t <- r$bootstrap$structural$items
    .pm_network_stability_points(t$stability,t$item,.pm_network_item_membership(x,t$item),
      tr('Estabilidad de \u00edtems','Item stability'),es)
  } else if (type == "dimension_stability") {
    t <- r$bootstrap$structural$dimensions
    .pm_network_stability_points(t$structural_consistency,as.character(t$dimension),t$dimension,
      tr('Consistencia estructural','Structural consistency'),es,axis_label=tr('Comunidades','Communities'))
  } else if (type == "edge_accuracy") {
    .pm_network_require("bootnet")
    p <- plot(r$bootstrap$accuracy, labels=FALSE, order="sample")
    if (inherits(p, "ggplot")) print(p)
  } else if (type == "centrality_stability") {
    .pm_network_require("bootnet")
    p <- plot(r$bootstrap$case, statistics=c("strength","betweenness","closeness"))
    if (inherits(p, "ggplot")) print(p)
  } else if (type == "group_networks") {
    graphs <- r$invariance$networks
    panels <- matrix(0L,ceiling(length(graphs)/2),min(2,length(graphs)))
    panels[seq_along(graphs)] <- seq_along(graphs)
    comparisons <- unique(r$invariance$loadings$comparison)
    annotate <- length(graphs)==2L && length(comparisons)==1L
    graphics::layout(rbind(panels,rep(length(graphs)+1L,ncol(panels))),heights=c(rep(1,nrow(panels)),if(annotate).25 else .16))
    .pm_network_require("qgraph")
    layout <- qgraph::qgraph(r$ega$network,layout="spring",DoNotPlot=TRUE)$layout
    maximum <- max(vapply(graphs,function(g)max(abs(g$network)),numeric(1)))
    states <- if(annotate).pm_network_invariance_status(x,comparisons[1],rownames(r$ega$network)) else NULL
    for(name in names(graphs)) draw_network(graphs[[name]]$network,r$ega$wc,rownames(r$ega$network),name,layout,maximum,show_legend=FALSE,node_status=states)
    graphics::par(mar=c(0,0,0,0));graphics::plot.new()
    palette <- .pm_network_community_colors(r$ega$wc)
    graphics::legend(x=.5,y=.8,xjust=.5,yjust=.5,bty='n',horiz=TRUE,cex=.85,
      legend=paste(tr('Comunidad','Community'),names(palette)),pch=21,pt.bg=.pm_network_color_fill(palette),col=palette,pt.cex=1.6)
    if(annotate) graphics::legend(.5,.46,xjust=.5,yjust=.5,bty='n',horiz=TRUE,cex=.7,
      legend=c(tr('Diferencia: p BH < .05','Difference: p BH < .05'),
        tr('Sin diferencia detectada: p BH >= .05','No difference detected: p BH >= .05'),
        tr('No evaluable','Not evaluable')),
      pch=c(21,21,22),pt.bg='white',col=c('#B34D00','#666666','#666666'),pt.lwd=c(3,1,1),pt.cex=1.5)
    graphics::text(.5,.12,tr('Colores = comunidades conjuntas; aristas negativas discontinuas; grosor = |peso|.',
      'Colors = pooled communities; negative edges dashed; width = |weight|.'),cex=.75)
  } else if (type == "edge_differences") {
    inv <- r$invariance
    graphics::par(mfrow=c(ceiling(length(inv$nct)/2),min(2,length(inv$nct))))
    for(name in names(inv$nct)) {
      n <- inv$nct[[name]]
      diff <- n$nw1-n$nw2
      tab <- inv$edges[inv$edges$comparison==name,]
      sig <- matrix(FALSE,nrow(diff),ncol(diff))
      for(i in seq_len(nrow(tab))) {
        a <- match(tab$item1[i],rownames(r$ega$network)); b <- match(tab$item2[i],rownames(r$ega$network))
        if (!is.na(a) && !is.na(b)) sig[a,b] <- sig[b,a] <- !is.na(tab$p_BH[i]) && tab$p_BH[i] < .05
      }
      diff[!sig] <- 0
      draw_network(diff,r$ega$wc,rownames(r$ega$network),paste(name,"BH < .05"))
    }
  } else if (type == "loading_differences") {
    .pm_network_loading_difference_plot(x,es)
  }
  invisible(r)
}

#' Export publication figures from cached network results
#' @param x A `psychomatic_network` object.
#' @param directory Existing output directory (created if absent).
#' @param types Figure types; NULL selects relevant publication figures,
#'   omitting unavailable wording diagnostics and descriptive TEFI/centrality
#'   figures. Supply `types=network_plots(x)` explicitly for all diagnostic types.
#'   With bootstrap results, `"ega_stability"` combines the EGA network and
#'   item replication plot. Standalone `"network"` and `"item_stability"`
#'   remain available explicitly. Two-group networks flag complete BH loading tests.
#' @param formats File formats, one or more of `"pdf"`, `"svg"`, `"png"`.
#' @param language `"eng"` or `"esp"`.
#' @param width,height Figure dimensions in inches.
#' @param dpi PNG resolution, default 300.
#' @param overwrite Allow existing files to be overwritten.
#' @param label_cex Node label size multiplier, default 1.2. Increase for
#'   figures that will be reduced to a small publication width.
#' @return Invisibly, paths of exported figures; without `directory`, available types.
#' @export
network_plots <- function(x, directory = NULL, types = NULL, formats = c("pdf","png"),
    language = NULL, width = 10, height = 7, dpi = 300, overwrite = FALSE,
    label_cex = 1.2) {
  if (!inherits(x,"psychomatic_network")) stop("`x` must be a psychomatic_network result.",call.=FALSE)
  available <- .pm_network_plot_types(x)
  if (is.null(directory)) return(invisible(available))
  if (is.null(types)) types <- .pm_network_publication_types(x)
  if (!is.character(types) || !length(types) || anyNA(types) || any(!types %in% available) || anyDuplicated(types)) stop("Unknown or duplicated figure types.",call.=FALSE)
  formats <- match.arg(formats,c("pdf","svg","png"),several.ok=TRUE)
  if (anyDuplicated(formats)) stop("Duplicate figure formats.",call.=FALSE)
  .pm_validate_scalar_number(width,"width",1,100); .pm_validate_scalar_number(height,"height",1,100)
  .pm_validate_scalar_number(dpi,"dpi",72,1200,integer=TRUE)
  .pm_validate_scalar_number(label_cex,"label_cex",.3,4)
  if (!is.character(directory) || length(directory)!=1L || is.na(directory) || !nzchar(directory)) stop("Invalid figure directory.",call.=FALSE)
  if (!is.logical(overwrite) || length(overwrite)!=1L || is.na(overwrite)) stop("Invalid overwrite flag.",call.=FALSE)
  paths <- as.vector(outer(types,formats,function(t,f)file.path(directory,paste0(t,".",f))))
  if (!overwrite && any(file.exists(paths))) stop("Figure files exist; use overwrite=TRUE explicitly.",call.=FALSE)
  if (!dir.exists(directory) && !dir.create(directory,recursive=TRUE)) stop("Cannot create figure directory.",call.=FALSE)
  for (format in formats) for(type in types) {
    path <- file.path(directory,paste0(type,".",format))
    .pm_network_draw_file(x,type,path,format,language,width,height,dpi,label_cex)
  }
  invisible(paths)
}

.pm_network_draw_file <- function(x,type,path,format="png",language=NULL,width=10,height=7,dpi=300,label_cex=1.2) {
  switch(format,png=grDevices::png(path,width=width,height=height,units="in",res=dpi),
    pdf=grDevices::pdf(path,width=width,height=height),svg=grDevices::svg(path,width=width,height=height))
  on.exit(grDevices::dev.off(),add=TRUE)
  plot(x,type=type,language=language,label_cex=label_cex)
  invisible(path)
}
