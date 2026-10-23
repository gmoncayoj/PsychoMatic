.pm_network_cross_loadings <- function(loadings) {
  m <- loadings$matrix
  if(is.null(m) || !is.matrix(m) || ncol(m)<2L) return(NULL)
  primary <- loadings$primary
  assigned <- primary$dimension[match(rownames(m),primary$item)]
  column <- match(as.character(assigned),colnames(m))
  rows <- which(!is.na(column))
  m[cbind(rows,column[rows])] <- NA_real_
  colnames(m) <- paste0('D',colnames(m))
  data.frame(item=rownames(m),dimension=assigned,m,check.names=FALSE,row.names=NULL)
}

.pm_network_item_membership <- function(x,items) {
  wc <- x$results$ega$wc
  if(is.null(names(wc))) names(wc) <- rownames(x$results$ega$network)
  unname(wc[match(items,names(wc))])
}

.pm_network_invariance_status <- function(x,comparison,items) {
  table <- x$results$invariance$loadings
  status <- rep('unavailable',length(items))
  if(!is.data.frame(table) || !nrow(table)) return(status)
  table <- .pm_network_loading_plot_data(table[table$comparison==comparison,,drop=FALSE],
    x$effective_settings$n_perm)
  row <- match(items,table$item)
  valid <- which(!is.na(row) & table$available[row])
  status[valid] <- ifelse(table$significant[row[valid]],'difference','no_difference')
  status
}

.pm_network_stability_points <- function(values,labels,membership,title,es=FALSE,legend=TRUE,axis_label=NULL) {
  tr <- function(esp,eng)if(es)esp else eng
  if(!length(values)) {
    graphics::plot.new();graphics::title(main=title)
    graphics::text(.5,.5,tr('No evaluable','Not evaluable'))
    return(invisible(data.frame(item=character(),stability=numeric())))
  }
  palette <- .pm_network_community_colors(membership)
  color <- unname(palette[as.character(membership)])
  color[is.na(color)] <- '#777777'
  order <- order(match(as.character(membership),names(palette)),seq_along(labels),na.last=TRUE)
  values <- values[order];labels <- labels[order];color <- color[order]
  y <- rev(seq_along(values));n <- length(values)
  gap <- if(legend)max(.6,.5*ceiling(length(palette)/4)) else .4
  if(is.null(axis_label)) axis_label <- tr('\u00cdtems','Items')
  graphics::par(mar=c(4.5,max(5,min(12,max(nchar(labels))*.5)),4,1.5))
  graphics::plot(NA,xlim=c(-.025,1.09),ylim=c(.4,n+gap),axes=FALSE,
    xlab=tr('Replicaci\u00f3n','Replication'),ylab='',main=title,xaxs='i',font.main=1)
  graphics::axis(1,at=c(0,.25,.5,.75,1),labels=c('0.00','0.25','0.50','0.75','1.00'))
  graphics::axis(2,at=y,labels=labels,las=1,cex.axis=1)
  graphics::mtext(axis_label,side=2,line=max(3,min(10,max(nchar(labels))*.5)),cex=.85)
  valid <- is.finite(values)
  if(any(valid)) {
    graphics::segments(0,y[valid],values[valid],y[valid],col='#D5D5D5',lwd=1)
    graphics::points(values[valid],y[valid],pch=21,cex=4.5,
      bg=.pm_network_color_fill(color[valid]),col=color[valid],lwd=1.2)
    text <- ifelse(values[valid]==1,'1',ifelse(values[valid]==0,'0',sub('^0','',sprintf('%.3f',values[valid]))))
    graphics::text(values[valid],y[valid],text,cex=.9,col='#000000')
  }
  if(any(!valid)) {
    graphics::points(rep(0,sum(!valid)),y[!valid],pch=4,col='#777777')
    graphics::text(rep(.12,sum(!valid)),y[!valid],tr('No evaluable','Not evaluable'),cex=.65)
  }
  if(legend && length(palette)) graphics::legend(x=.5,y=n+gap-.05,xjust=.5,yjust=1,bty='n',
    legend=paste(tr('Comunidad','Community'),names(palette)),ncol=min(4,length(palette)),
    pch=21,pt.bg=.pm_network_color_fill(palette),col=palette,pt.cex=1.5,cex=.7)
  invisible(data.frame(item=labels,stability=values))
}
