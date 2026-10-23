.pm_network_community_colors <- function(wc) {
  ids <- sort(unique(as.character(wc[!is.na(wc)])))
  # Paul Tol's bright scheme, reordered to separate the first five communities.
  palette <- c('#4477AA','#EE6677','#228833','#AA3377','#CCBB44','#66CCEE','#BBBBBB')
  colors <- if(length(ids)<=length(palette)) palette[seq_along(ids)] else grDevices::hcl.colors(length(ids),'Dark 3')
  stats::setNames(colors,ids)
}

.pm_network_color_fill <- function(colors) {
  rgb <- grDevices::col2rgb(colors)*.75 + 255*.25
  stats::setNames(grDevices::rgb(rgb[1,],rgb[2,],rgb[3,],maxColorValue=255),names(colors))
}

.pm_network_contrast_text <- function(colors) {
  rgb <- grDevices::col2rgb(colors)/255
  linear <- ifelse(rgb<=.04045,rgb/12.92,((rgb+.055)/1.055)^2.4)
  luminance <- colSums(linear*c(.2126,.7152,.0722))
  ifelse((1.05/(luminance+.05))>((luminance+.05)/.05),'#FFFFFF','#000000')
}

.pm_network_publication_types <- function(x) {
  available <- .pm_network_plot_types(x)
  types <- c(if(!is.null(x$results$bootstrap)) 'ega_stability' else 'network','loadings','loadings_matrix',
    if(!is.null(x$results$uva)) 'uva',
    if(identical(x$results$wording$status,'completed')) 'wording',
    'hierarchy','dimensions','dimension_stability',
    'edge_accuracy','centrality_stability','group_networks','loading_differences')
  if(is.data.frame(x$results$invariance$edges) && any(is.finite(x$results$invariance$edges$p_BH) & x$results$invariance$edges$p_BH<.05))
    types <- c(types,'edge_differences')
  intersect(types,available)
}

.pm_network_loading_plot_data <- function(table, n_perm) {
  table$available <- is.finite(table$p_BH) & !is.na(table$evaluable) & table$evaluable &
    !is.na(table$valid_permutations) & table$valid_permutations == n_perm
  table$significant <- table$available & table$p_BH < .05
  table
}

.pm_network_loading_difference_plot <- function(x, es=FALSE) {
  tr <- function(esp,eng)if(es)esp else eng
  t <- .pm_network_loading_plot_data(x$results$invariance$loadings,x$effective_settings$n_perm)
  comparisons <- unique(t$comparison)
  graphics::par(mfrow=c(ceiling(length(comparisons)/2),min(2,length(comparisons))))
  maxdiff <- max(abs(t$difference[is.finite(t$difference)]),.05)
  for (comparison in comparisons) {
    tab <- t[t$comparison==comparison,,drop=FALSE]
    y <- rev(seq_len(nrow(tab)))
    graphics::par(mar=c(5,6,4,2))
    graphics::plot(NA,xlim=c(-1.1*maxdiff,1.8*maxdiff),ylim=c(.5,nrow(tab)+1),
      axes=FALSE,xlab=tr('Diferencia de cargas (primer grupo menos segundo)','Loading difference (first group minus second)'),
      ylab='',main=comparison)
    ticks <- pretty(c(-maxdiff,maxdiff))
    graphics::axis(1,at=ticks[abs(ticks)<=1.1*maxdiff])
    graphics::axis(2,at=y,labels=tab$item,las=1,cex.axis=.8)
    graphics::abline(v=0,lty=2,col='#777777')
    finite <- is.finite(tab$difference)
    graphics::points(tab$difference[finite],y[finite],
      pch=ifelse(!tab$available[finite],4,ifelse(tab$significant[finite],19,1)),
      col=ifelse(tab$significant[finite],'#D55E00','#444444'),cex=1.2)
    if(any(!finite)) graphics::text(0,y[!finite],tr('No evaluable','Not evaluable'),cex=.7)
    p <- ifelse(tab$available,ifelse(tab$p_BH<.001,'<.001',sprintf('%.3f',tab$p_BH)), 'NA')
    graphics::text(1.4*maxdiff,y,labels=p,cex=.75)
    graphics::text(1.4*maxdiff,nrow(tab)+.65,labels='p BH',cex=.8,font=2)
    graphics::legend('topleft',bty='n',horiz=TRUE,cex=.65,
      legend=c('BH < .05','BH >= .05',tr('No evaluable','Not evaluable')),
      pch=c(19,1,4),col=c('#D55E00','#444444','#444444'))
  }
  invisible(t)
}
