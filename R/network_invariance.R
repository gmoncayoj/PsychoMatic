.pm_network_invariance <- function(data, groups, ega, settings) {
  complete <- !is.na(groups)
  data <- data[complete, , drop = FALSE]; groups <- as.character(groups[complete])
  levels <- unique(groups)
  if (length(levels) < 2L) stop("Invariance requires at least two observed groups.", call. = FALSE)
  if (any(table(groups) < 20L)) stop("Every invariance group requires at least 20 usable observations.", call. = FALSE)
  if (any(vapply(levels, function(g) any(vapply(data[groups == g, , drop = FALSE],
      function(x) length(unique(x[!is.na(x)])) < 2L, logical(1))), logical(1))))
    stop("Every item must vary in every invariance group.", call. = FALSE)
  pairs <- utils::combn(levels, 2, simplify = FALSE)
  networks <- setNames(lapply(levels, function(g) .pm_network_ega(data[groups == g, , drop = FALSE], settings)), levels)
  # Structural membership is aligned by co-membership, never by arbitrary
  # community label numbers. Stability is measured against the pooled structure.
  config <- lapply(levels, function(g) {
    .pm_network_progress(settings,sprintf("Configural bootstrap: %s (%d resamples)",g,settings$n_boot))
    b <- .pm_network_boot_summary(.pm_network_boot_ega(data[groups == g, , drop = FALSE], settings), ega$wc)
    data.frame(group = g, b$items, same_empirical_partition =
      identical(unname(outer(networks[[g]]$wc, networks[[g]]$wc, "==")),
        unname(outer(ega$wc, ega$wc, "=="))),
      valid_bootstraps = b$valid)
  })
  config <- do.call(rbind, config)
  configural_supported <- all(is.finite(config$stability) & config$stability >= .70) &&
    all(config$same_empirical_partition)
  nct <- list(); global <- edge <- centrality <- data.frame()
  estimator <- .pm_network_estimator_function(settings)
  for (pair in pairs) {
    label <- paste(pair, collapse = " vs ")
    .pm_network_progress(settings,sprintf("NCT: %s (%d permutations; sequential)",label,settings$n_perm))
    native <- NetworkComparisonTest::NCT(as.matrix(data[groups == pair[1], , drop = FALSE]),
      as.matrix(data[groups == pair[2], , drop = FALSE]), it = settings$n_perm,
      test.edges = TRUE, test.centrality = TRUE,
      centrality = c("strength", "betweenness", "closeness"), estimator = estimator,
      p.adjust.methods = "none", progressbar = FALSE, verbose = FALSE)
    nct[[label]] <- native
    global <- rbind(global, data.frame(comparison = label,
      test = c("network_structure", "global_strength"),
      statistic = c(native$nwinv.real, native$glstrinv.real),
      p_raw = c(native$nwinv.pval, native$glstrinv.pval), family = "global"))
    edge <- rbind(edge, .pm_network_nct_edges(native, label, names(data)))
    m <- native$diffcen.pval
    if (!is.null(m)) centrality <- rbind(centrality, data.frame(comparison = label,
      item = rep(rownames(m), ncol(m)), metric = rep(colnames(m), each = nrow(m)),
      difference = as.vector(native$diffcen.real), p_raw = as.vector(m), family = "centrality"))
    connected <- all(vapply(networks[pair],function(g)all(.pm_network_centrality(g$network)$connected),logical(1)))
    if (!connected && nrow(centrality)) {
      undefined <- centrality$comparison == label & centrality$metric == "closeness"
      centrality$p_raw[undefined] <- NA_real_
      centrality$difference[undefined] <- NA_real_
    }
  }
  # Supplying a structure to EGAnet skips its configural reduction. We explicitly
  # test the original shared item set above and never silently drop unstable items.
  .pm_network_progress(settings,sprintf("EGA loadings: %d permutations; %d workers; fixed partition",
    settings$n_perm,settings$n_cores))
  native_metric <- .pm_network_metric_invariance(data,groups,ega$wc,settings)
  loading <- data.frame()
  for (i in seq_along(native_metric$groups$loadingsDifference)) {
    actual <- native_metric$groups$loadingsDifference[[i]]
    permutations <- native_metric$permutation$loadingsDifference[[i]]
    valid <- rowSums(is.finite(permutations))
    # Recompute unrounded Monte Carlo p from EVERY native permutation. EGAnet
    # 2.4.1 rounds its public p table and excludes the first draw in that table.
    p <- .pm_network_perm_p(actual, permutations, settings$n_perm)
    loading <- rbind(loading, data.frame(comparison = paste(pairs[[i]], collapse = " vs "),
      item = names(actual), dimension = as.character(ega$wc[names(actual)]),
      difference = as.numeric(actual), p_raw = p, valid_permutations = valid,
      evaluable = is.finite(p), family = "loadings"))
  }
  list(groups = data.frame(group = levels, n = as.integer(table(factor(groups, levels = levels)))),
    omitted_missing_group = sum(!complete), networks = networks, nct = nct,
    global = .pm_network_bh(global), edges = .pm_network_bh(edge),
    centrality = .pm_network_bh(centrality), loadings = .pm_network_bh(loading),
    configural = config, configural_supported = configural_supported,
    metric_status = if (configural_supported) "conditional_on_supported_configural_structure" else "exploratory_configural_not_supported",
    native_metric = native_metric, original_item_set = names(data),
    p_resolution = 1 / (settings$n_perm + 1),
    families = data.frame(family = c("global", "edges", "centrality", "loadings"),
      size = c(nrow(global), nrow(edge), nrow(centrality), nrow(loading)),
      scope = "All declared group pairs; one BH adjustment within each family"))
}

.pm_network_metric_invariance <- function(data,groups,structure,settings) {
  # The native test computes loads from each re-estimated network and the
  # supplied structure, ignoring newly detected communities. Walktrap/LE avoids
  # thousands of unused consensus Louvain runs; GLASSO, rows, fixed structure,
  # revised loadings and native permutation seeds remain the same. LE only
  # checks membership: unlike `expand`, it does not change the network matrix.
  args <- c(list(data=as.matrix(data),groups=groups,structure=structure,
    na.data=settings$missing,model=.pm_null_default(settings$model,"glasso"),
    algorithm="walktrap",uni.method="LE",
    loading.method="revised",iter=settings$n_perm,
    ncores=settings$n_cores,seed=settings$seed,verbose=FALSE),
    if(is.null(settings$model) || settings$model=="glasso")list(gamma=settings$gamma),
    .pm_network_corr_args(settings$corr,data,settings$ordinal_categories))
  do.call(EGAnet::invariance,args)
}

.pm_network_nct_edges <- function(native, label, items) {
  e <- native$einv.pvals
  labels <- function(x) {
    x <- as.character(x)
    if (all(x %in% items)) return(x)
    index <- suppressWarnings(as.integer(x))
    if (anyNA(index) || any(index < 1L | index > length(items)))
      stop("Unrecognized NCT edge labels; cannot align group networks.", call. = FALSE)
    items[index]
  }
  first <- labels(e[[1]]); second <- labels(e[[2]])
  difference <- native$nw1 - native$nw2
  data.frame(comparison = label, item1 = first, item2 = second, p_raw = e[[3]],
    difference = as.numeric(difference[cbind(match(first,items),match(second,items))]), family = "edges")
}

.pm_network_perm_p <- function(actual, permutations, requested) {
  valid <- rowSums(is.finite(permutations))
  p <- (1 + rowSums(abs(permutations) >= abs(actual), na.rm = TRUE)) / (valid + 1)
  p[!is.finite(actual) | valid != requested] <- NA_real_
  p
}
