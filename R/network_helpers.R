# The same estimator is called on raw rows by EGA, bootnet and NCT.
.pm_network_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE))
    stop(sprintf("Install optional package '%s' to run the requested network stage.", package), call. = FALSE)
}

.pm_network_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit(if (is.null(old)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  } else assign(".Random.seed", old, envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  force(code)
}

.pm_network_partition_evaluable <- function(ega) {
  length(ega$wc) > 0L && !anyNA(ega$wc) &&
    length(ega$n.dim) == 1L && is.finite(ega$n.dim) && ega$n.dim > 0
}

.pm_network_input <- function(data, items, group, reverse, ranges, missing) {
  spec <- if (inherits(data, "psychomatic_spec")) data else
    psychomatic_spec(data, items = items, group = group, reverse = reverse, ranges = ranges)
  if (inherits(data, "psychomatic_spec") &&
      any(vapply(list(items, group, reverse, ranges), Negate(is.null), logical(1))))
    stop("Supply items, group, reverse and ranges in the psychomatic_spec, or use a data frame.", call. = FALSE)
  raw <- spec$data[spec$items]
  if (anyDuplicated(names(raw)) || any(!nzchar(names(raw))))
    stop("Network items must have unique, non-empty names.", call. = FALSE)
  raw[] <- lapply(raw, function(x) if (is.ordered(x)) as.numeric(x) else x)
  if (length(raw) < 3L) stop("At least three network items are required.", call. = FALSE)
  if (any(!vapply(raw, is.numeric, logical(1))) ||
      any(vapply(raw, function(x) any(!is.finite(x) & !is.na(x)), logical(1))))
    stop("Items must be numeric or ordered, with finite observed values.", call. = FALSE)
  keyed <- raw
  for (item in spec$reverse) {
    bounds <- spec$ranges[[item]]
    if (is.null(bounds)) stop("Explicit theoretical `ranges` are required for reverse items.", call. = FALSE)
    observed <- raw[[item]][!is.na(raw[[item]])]
    if (any(observed < bounds[1] | observed > bounds[2]))
      stop(sprintf("Observed values of '%s' are outside its declared range.", item), call. = FALSE)
    keyed[[item]] <- sum(bounds) - raw[[item]]
  }
  groups <- if (is.null(spec$group)) NULL else spec$data[[spec$group]]
  keep <- if (missing == "listwise") stats::complete.cases(keyed) else
    rowSums(!is.na(keyed)) > 0L
  raw <- raw[keep, , drop = FALSE]; keyed <- keyed[keep, , drop = FALSE]
  if (nrow(raw) < 20L) stop("At least 20 usable observations are required.", call. = FALSE)
  if (any(vapply(raw, function(x) length(unique(x[!is.na(x)])) < 2L, logical(1))))
    stop("Items must vary in the analysis sample.", call. = FALSE)
  list(raw = raw, keyed = keyed, groups = if (is.null(groups)) NULL else groups[keep],
       spec = spec, omitted = sum(!keep), row_index = which(keep))
}

.pm_network_estimator <- function(data, gamma = .5, corr = "spearman", missing = "listwise") {
  args <- c(list(data=as.matrix(data),na.data=missing,gamma=gamma,
    model.selection="EBIC",verbose=FALSE),.pm_network_corr_args(corr,data))
  as.matrix(do.call(EGAnet::EBICglasso.qgraph,args))
}

.pm_network_estimator_function <- function(settings) {
  # PSOCK workers may resolve the package namespace from a different library.
  # Serialize only the native arguments, not a private PsychoMatic function or
  # the enclosing analysis (which also contains large cached bootstrap data).
  model <- .pm_null_default(settings$model,"glasso")
  args <- c(list(na.data=settings$missing,verbose=FALSE),
    if(model=="glasso")list(gamma=settings$gamma,model.selection="EBIC") else
      list(model=model,network.only=TRUE),if(settings$corr=="polychoric")
      list(corr="auto",ordinal.categories=.pm_null_default(settings$ordinal_categories,7L)) else
      list(corr=settings$corr))
  native <- if(model=="glasso")"EBICglasso.qgraph" else "network.estimation"
  fun <- function(data) as.matrix(do.call(getExportedValue("EGAnet",native),
    c(list(data=as.matrix(data)),args)))
  environment(fun) <- list2env(list(args=args,native=native),parent=baseenv())
  fun
}

.pm_network_corr_args <- function(corr, data, ordinal_categories = NULL) {
  if(corr != "polychoric")return(list(corr=corr))
  if(is.null(ordinal_categories))ordinal_categories <- max(7L,
    vapply(as.data.frame(data),function(x)length(unique(x[!is.na(x)])),integer(1)))
  if(ordinal_categories>11L)
    stop('The EGAnet polychoric estimator supports at most 11 ordinal categories per item.',call.=FALSE)
  list(corr="auto",ordinal.categories=ordinal_categories)
}

.pm_network_progress <- function(settings, text) {
  if(isTRUE(settings$verbose)) {
    cat(sprintf("[%s] %s\n",format(Sys.time(),"%H:%M:%S"),text))
    utils::flush.console()
  }
  invisible(NULL)
}

.pm_network_interrupted <- function(condition) {
  inherits(condition,"interrupt") || grepl("^Interrupted($|\\n)",conditionMessage(condition))
}

.pm_network_ega <- function(data, settings, variant = "EGA") {
  args <- c(list(data = as.matrix(data), na.data = settings$missing,
    seed = settings$seed, plot.EGA = FALSE),.pm_network_native_options(settings),
    .pm_network_corr_args(settings$corr,data,settings$ordinal_categories))
  do.call(getExportedValue("EGAnet", variant), args)
}

# Defaults allow cached objects and older private helper calls to remain usable.
.pm_network_native_options <- function(settings, hierarchy = FALSE) {
  model <- .pm_null_default(settings$model,"glasso")
  algorithm <- .pm_null_default(settings$algorithm,"louvain")
  args <- list(model=model)
  if(hierarchy) {
    args$lower.algorithm <- algorithm; args$higher.algorithm <- algorithm
  } else args$algorithm <- algorithm
  if(model=="glasso") args$gamma <- settings$gamma
  if(algorithm=="louvain") args$consensus.iter <- settings$consensus_iter else
    args$uni.method <- "LE"
  if(algorithm=="leiden") args$objective_function <- "modularity"
  args
}

.pm_network_edge_kind <- function(x, es = FALSE) {
  if(identical(x$effective_settings$model,"TMFG"))
    if(es)"correlaciones filtradas por TMFG" else "TMFG-filtered correlations" else
    if(es)"correlaciones parciales" else "partial correlations"
}

.pm_network_method <- function(settings, es = FALSE) {
  model <- .pm_null_default(settings$model,"glasso")
  algorithm <- .pm_null_default(settings$algorithm,"louvain")
  estimation <- if(model=="glasso")sprintf("glasso (EBICglasso, gamma = %.2f)",settings$gamma) else model
  community <- if(algorithm=="louvain")sprintf(if(es)"Louvain por consenso (%d iteraciones)" else
    "consensus Louvain (%d iterations)",settings$consensus_iter) else algorithm
  paste(estimation, settings$corr, community, sep="; ")
}

.pm_network_loading_magnitude <- function(x) {
  out <- rep(NA_character_, length(x)); ok <- is.finite(x)
  out[ok] <- c("below_small", "small", "moderate", "large")[
    findInterval(abs(x[ok]), c(.20, .35, .50)) + 1L]
  out
}

.pm_network_loadings <- function(ega) {
  native <- EGAnet::net.loads(ega, loading.method = "revised", ordered = "variable")
  m <- native$std
  membership <- as.character(ega$wc[rownames(m)])
  primary <- m[cbind(seq_len(nrow(m)), match(membership, colnames(m)))]
  table <- data.frame(item = rownames(m), dimension = membership, loading = primary,
    magnitude = .pm_network_loading_magnitude(primary), stringsAsFactors = FALSE)
  list(matrix = m, primary = table, native = native)
}

.pm_network_centrality <- function(network) {
  .pm_network_require("igraph")
  a <- abs(network); diag(a) <- 0
  graph <- igraph::graph_from_adjacency_matrix(a, mode = "undirected", weighted = TRUE, diag = FALSE)
  connected <- igraph::is_connected(graph)
  distance <- 1 / igraph::E(graph)$weight
  closeness <- if (connected) igraph::closeness(graph, weights = distance, normalized = FALSE)
    else rep(NA_real_, nrow(a))
  data.frame(item = rownames(a), strength = rowSums(a),
    betweenness = as.numeric(igraph::betweenness(graph, directed = FALSE, weights = distance, normalized = FALSE)),
    closeness = as.numeric(closeness), expected_influence = rowSums(network),
    connected = connected, stringsAsFactors = FALSE)
}

.pm_network_bh <- function(table) {
  table$p_BH <- NA_real_
  for (family in unique(table$family)) {
    idx <- which(table$family == family)
    table$p_BH[idx] <- stats::p.adjust(table$p_raw[idx], method = "BH")
  }
  table
}
