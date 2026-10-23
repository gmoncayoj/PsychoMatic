# Reader summaries operate on cached estimates; filtering never changes BH families.
.pm_network_sections <- function(x, language = NULL, destination = "console",
                                 detail = c("report", "audit")) {
  detail <- match.arg(detail)
  out <- .pm_network_audit_sections(x, language, destination)
  if (detail == "audit") return(out)
  es <- out$language == "esp"
  tr <- function(esp, eng) if (es) esp else eng
  r <- x$results; sections <- list()
  add <- function(esp, eng, value) {
    if (is.null(value) || !length(value)) return(invisible(NULL))
    if (is.character(value)) value <- data.frame(Text = value)
    if (is.data.frame(value) && nrow(value))
      sections[[tr(esp, eng)]] <<- .pm_network_translate_table(value, es)
  }
  keep <- function(esp, eng) {
    value <- out$sections[[tr(esp, eng)]]
    if (!is.null(value)) sections[[tr(esp, eng)]] <<- value
  }
  add('Resumen', 'Summary', tr(
    sprintf('Se analizaron %d observaciones y %d de %d \u00edtems. EGA detect\u00f3 %d comunidades.',
      x$data_summary$n, x$data_summary$p_retained, x$data_summary$p_original, r$ega$n.dim),
    sprintf('The analysis used %d observations and %d of %d items. EGA detected %d communities.',
      x$data_summary$n, x$data_summary$p_retained, x$data_summary$p_original, r$ega$n.dim)))
  keep('M\u00e9todo', 'Method')
  primary <- r$loadings$primary
  if (is.data.frame(primary) && nrow(primary)) {
    communities <- lapply(unique(primary$dimension[!is.na(primary$dimension)]), function(d)
      data.frame(dimension = d, items = paste(primary$item[which(primary$dimension == d)], collapse = ', ')))
    if (length(communities)) add('Comunidades detectadas', 'Detected communities', do.call(rbind, communities))
    keep('Cargas primarias revisadas', 'Revised primary loadings')
    cross <- .pm_network_cross_loadings(r$loadings)
    if(!is.null(cross)) {
      add('Cargas cruzadas revisadas', 'Revised cross-loadings',cross)
      add('Lectura de cargas cruzadas','Cross-loading guide',tr(
        'D indica comunidad. Las celdas de la comunidad asignada se dejan vac\u00edas porque sus cargas primarias se reportan arriba; se conservan todas las cargas cruzadas y sus signos.',
        'D denotes community. Assigned-community cells are blank because primary loadings appear above; all cross-loadings and their signs are retained.'))
    }
    add('Interpretaci\u00f3n de cargas', 'Loading interpretation', tr(
      'Magnitudes orientativas de |carga|: .20 peque\u00f1a, .35 moderada, .50 grande (Christensen et al., 2025); no son criterios autom\u00e1ticos de eliminaci\u00f3n.',
      'Approximate |loading| magnitudes: .20 small, .35 moderate, .50 large (Christensen et al., 2025); these are not automatic deletion rules.'))
    add('Precisi\u00f3n de cargas','Loading precision',tr(
      'La magnitud se clasifica con valores sin redondear; una carga impresa como .500 puede estar ligeramente por debajo de .50.',
      'Magnitudes use unrounded values; a loading printed as .500 may be slightly below .50.'))
  }
  if (!is.null(r$uva)) {
    keep('Decisi\u00f3n UVA', 'UVA decision')
    keep('Pares redundantes', 'Redundant pairs')
    keep('\u00cdtems eliminados por UVA', 'Items removed by UVA')
    keep('Propuesta UVA y motivos de selecci\u00f3n', 'UVA proposal and selection reasons')
    keep('Interpretaci\u00f3n de UVA', 'UVA interpretation')
  }
  if (!is.null(r$bootstrap)) {
    b <- r$bootstrap
    keep('M\u00e9todo bootstrap', 'Bootstrap method')
    keep('R\u00e9plicas bootstrap', 'Bootstrap replications')
    keep('Consistencia estructural', 'Structural consistency')
    items <- b$structural$items
    if(is.data.frame(items) && nrow(items)) {
      assigned <- .pm_network_item_membership(x,items$item)
      add('Estabilidad de \u00edtems','Item stability',
        data.frame(item=items$item,dimension=as.character(assigned),stability=items$stability))
    }
    f <- b$structural$frequency
    if (is.data.frame(f) && all(c('dimensions','frequency') %in% names(f))) {
      frequency <- f$frequency[which(f$dimensions == r$ega$n.dim)]
      if (length(frequency) && is.finite(frequency[1])) add('Interpretaci\u00f3n bootstrap', 'Bootstrap interpretation', tr(
        sprintf('El n\u00famero de comunidades detectado por EGA (%d) se recuper\u00f3 en %.1f%% de las r\u00e9plicas v\u00e1lidas. La consistencia estructural indica la proporci\u00f3n que reprodujo exactamente cada comunidad; no equivale a la estabilidad media de sus \u00edtems (Christensen & Golino, 2021).', r$ega$n.dim, 100 * frequency[1]),
        sprintf('The number of communities detected by EGA (%d) was recovered in %.1f%% of valid replications. Structural consistency is the proportion reproducing each exact community; it differs from average item stability (Christensen & Golino, 2021).', r$ega$n.dim, 100 * frequency[1])))
    }
    items <- b$structural$items
    if (is.data.frame(items) && nrow(items) && all(c('item','stability') %in% names(items))) {
      weak <- items[!is.finite(items$stability) | items$stability < .70, , drop=FALSE]
      add('\u00cdtems que requieren revisi\u00f3n de estabilidad', 'Items requiring stability review',
        if (nrow(weak)) weak else tr(
          sprintf('Todos los \u00edtems alcanzaron estabilidad \u2265 .70; rango observado %.3f\u2013%.3f.',min(items$stability),max(items$stability)),
          sprintf('All items reached stability >= .70; observed range %.3f\u2013%.3f.',min(items$stability),max(items$stability))))
      add('Criterio de estabilidad', 'Stability criterion', tr(
        '.70 es una gu\u00eda operativa para revisar estabilidad de \u00edtems, no un umbral universal de validez. Las proporciones usan r\u00e9plicas v\u00e1lidas y el mismo conjunto de \u00edtems.',
        '.70 is an operational guide for reviewing item stability, not a universal validity cutoff. Proportions use valid replications and the fixed analysis item set.'))
    }
    keep('Coeficiente CS', 'CS coefficient')
    add('Interpretaci\u00f3n de centralidad', 'Centrality interpretation', tr(
      'CS < .25 indica estabilidad insuficiente; \u2265 .50 es preferible. No interpretar diferencias de intermediaci\u00f3n o cercan\u00eda cuando sus CS sean insuficientes (Epskamp et al., 2018).',
      'CS < .25 indicates inadequate stability; >= .50 is preferable. Avoid interpreting betweenness or closeness differences when their CS is inadequate (Epskamp et al., 2018).'))
  }
  if (is.data.frame(r$tefi)) {
    pooled <- r$tefi[which(r$tefi$model == 'EGA'), , drop=FALSE]
    add('Entrop\u00eda EGA', 'EGA entropy', pooled[intersect(c('model','index','value'),names(pooled))])
    add('Interpretaci\u00f3n TEFI', 'TEFI interpretation',tr(
      'TEFI expresa ajuste relativo: comparar solo modelos sobre la misma matriz y el mismo conjunto de \u00edtems. Un valor menor indica mejor ajuste relativo (Golino et al., 2021).',
      'TEFI measures relative fit: compare models only on the same matrix and item set. Lower values indicate better relative fit (Golino et al., 2021).'))
  }
  if (!is.null(r$wording) && identical(r$wording$status,'completed')) {
    keep('Diagn\u00f3stico de redacci\u00f3n', 'Wording diagnostic')
    add('Interpretaci\u00f3n riEGA', 'riEGA interpretation',tr(
      'riEGA usa respuestas originales sin recodificar. Los signos de asociaciones o un intercepto aleatorio no demuestran efectos de redacci\u00f3n ni identifican claves inversas.',
      'riEGA uses original unrecoded responses. Association signs or a random intercept do not establish wording effects or identify reverse keys.'))
  }
  if (!is.null(r$hierarchy)) {
    structure <- r$hierarchy$network_structure
    if(is.data.frame(structure) && 'dimensions' %in% names(structure)) structure$dimensions <- as.integer(structure$dimensions)
    add('Jerarqu\u00eda de redes (hierEGA)', 'Network hierarchy (hierEGA)',structure)
    keep('Interpretaci\u00f3n nativa exploratoria hierEGA', 'Native exploratory hierEGA interpretation')
    add('Interpretaci\u00f3n jer\u00e1rquica', 'Hierarchy interpretation', tr(
      'hierEGA eval\u00faa comunidades inferiores y superiores desde redes. Su implementaci\u00f3n es experimental y la interpretaci\u00f3n TEFI es exploratoria. El AFC suplementario, estimado en la misma muestra, no constituye confirmaci\u00f3n independiente (Jim\u00e9nez et al., 2025).',
      'hierEGA evaluates lower- and higher-order network communities. Its implementation is experimental and its TEFI interpretation is exploratory. Supplementary CFA on the same sample is not independent confirmation (Jim\u00e9nez et al., 2025).'))
    candidates <- r$hierarchy$candidates
    if (is.data.frame(candidates) && any(grepl('negative_variance', candidates$reason, fixed=TRUE))) {
      bad <- candidates$model[grepl('negative_variance', candidates$reason, fixed=TRUE)]
      add('Limitaci\u00f3n del AFC suplementario', 'Supplementary CFA limitation', tr(
        paste('Candidatos inadmisibles por varianza negativa (caso Heywood):', paste(bad,collapse=', '), '. Un CFI mayor no permite elegirlos.'),
        paste('Candidates inadmissible because of negative variance (Heywood case):', paste(bad,collapse=', '), '. Higher CFI does not justify selecting them.')))
    }
    add('Conclusi\u00f3n del AFC suplementario', 'Supplementary CFA conclusion',tr(
      paste('Selecci\u00f3n suplementaria:',r$hierarchy$recommendation,'.',.pm_network_reason_text(r$hierarchy$reason,TRUE)),
      paste('Supplementary selection:',r$hierarchy$recommendation,'.',r$hierarchy$reason)))
  }
  if (!is.null(r$invariance)) {
    i <- r$invariance
    keep('Grupos', 'Groups')
    add('M\u00e9todo de invarianza', 'Invariance method', tr(
      sprintf('NCT y permutaciones de cargas EGA (%d permutaciones; resoluci\u00f3n p = %.6f). BH se aplica por familia (global, aristas, centralidad, cargas), incluyendo todos los pares (Benjamini & Hochberg, 1995; van Borkulo et al., 2023; Jamison et al., 2024).', x$effective_settings$n_perm,i$p_resolution),
      sprintf('NCT and EGA loading permutations (%d permutations; p resolution = %.6f). BH is applied per family (global, edges, centrality, loadings), across all pairs (Benjamini & Hochberg, 1995; van Borkulo et al., 2023; Jamison et al., 2024).', x$effective_settings$n_perm,i$p_resolution)))
    global <- i$global
    if (is.data.frame(global) && all(c('test','p_BH') %in% names(global))) {
      add('Pruebas globales NCT', 'NCT global tests', global[intersect(c('comparison','test','statistic','p_BH'),names(global))])
      text <- vapply(seq_len(nrow(global)), function(j) {
        metric <- if (global$test[j] == 'global_strength') tr('fuerza global','global strength') else tr('estructura de red','network structure')
        comparison <- global$comparison[j]
        p <- global$p_BH[j]
        if (!is.finite(p)) return(tr(paste(comparison,':',metric,'no evaluable.'),paste(comparison,':',metric,'not evaluable.')))
        if (p < .05) tr(paste(comparison,': se detect\u00f3 diferencia de',metric,'(p BH < .05).'),paste(comparison,':',metric,'difference detected (p BH < .05).')) else
          tr(paste(comparison,': no se detect\u00f3 diferencia de',metric,'(p BH \u2265 .05).'),paste(comparison,': no',metric,'difference detected (p BH >= .05).'))
      }, character(1))
      add('Interpretaci\u00f3n NCT', 'NCT interpretation', text)
    }
    loads <- i$loadings
    supported <- isTRUE(i$configural_supported)
    add('Invarianza configural', 'Configural invariance', if (supported) tr(
      'La partici\u00f3n original recibi\u00f3 apoyo en todos los grupos seg\u00fan la comparaci\u00f3n de comunidades y estabilidad de \u00edtems (gu\u00eda operativa .70).',
      'The original partition was supported in all groups by community agreement and item stability (operational guide .70).') else tr(
      'No se apoy\u00f3 la estructura configural original; las diferencias de cargas se interpretan como exploratorias, no como evidencia confirmatoria de invarianza m\u00e9trica.',
      'The original configural structure was not supported; loading differences are exploratory, not confirmatory evidence of metric invariance.'))
    if (is.data.frame(loads) && nrow(loads) && all(c('p_BH','evaluable','valid_permutations') %in% names(loads))) {
      valid <- !is.na(loads$evaluable) & loads$evaluable & is.finite(loads$p_BH) &
        !is.na(loads$valid_permutations) & loads$valid_permutations == x$effective_settings$n_perm
      flagged <- loads[which(valid & loads$p_BH < .05), , drop=FALSE]
      add('\u00cdtems con diferencias de cargas (BH)', 'Items with loading differences (BH)',
        if(nrow(flagged)) flagged[intersect(c('comparison','item','dimension','difference','p_BH'),names(flagged))] else tr(
          'Ning\u00fan \u00edtem mostr\u00f3 diferencias de cargas con p BH < .05 entre las pruebas evaluables.',
          'No item showed loading differences with p BH < .05 among evaluable tests.'))
      add('Interpretaci\u00f3n de invarianza m\u00e9trica', 'Metric invariance interpretation',
        if(any(!valid)) tr('Evaluaci\u00f3n incompleta: hay pruebas de cargas no evaluables. No concluir invarianza para esos \u00edtems.',
          'Evaluation incomplete: some loading tests are not evaluable. Do not conclude invariance for those items.') else if(nrow(flagged)) tr(
          '\u00cdtems con posible falta de invarianza m\u00e9trica de red: revisar las diferencias y el apoyo configural.',
          'Items with possible metric network non-invariance: review the differences and configural support.') else tr(
          'No se detect\u00f3 falta de invarianza en las cargas de red. Esto no demuestra equivalencia ni invarianza escalar AFC.',
          'No metric network loading non-invariance was detected. This does not establish equivalence or CFA scalar invariance.'))
    } else add('Interpretaci\u00f3n de invarianza m\u00e9trica', 'Metric invariance interpretation',tr(
      'Las pruebas de cargas no est\u00e1n disponibles; no concluir invarianza m\u00e9trica.',
      'Loading tests are unavailable; do not conclude metric invariance.'))
  }
  stages <- r$stages
  unavailable <- stages[which(stages$status %in% c('failed','interrupted','not_evaluable')), , drop=FALSE]
  if(nrow(unavailable)) {
    # Describe unavailable scientific conclusions once; retain original conditions in audit.
    unavailable <- unavailable[!unavailable$stage %in% c('riEGA stability'),,drop=FALSE]
    notes <- vapply(seq_len(nrow(unavailable)), function(j) {
      if (unavailable$stage[j]=='riEGA') return(tr(
        'riEGA: partici\u00f3n residual no evaluable; no se concluye presencia ni ausencia de efectos de redacci\u00f3n.',
        'riEGA: residual partition not evaluable; neither the presence nor absence of wording effects can be concluded.'))
      paste0(unavailable$stage[j], ': ',
        if(es) c(failed='fallida',interrupted='interrumpida',not_evaluable='no evaluable')[[unavailable$status[j]]] else gsub('_',' ',unavailable$status[j]),
        tr('; no permite una conclusi\u00f3n para esa etapa.', '; no conclusion is available for this stage.'))
    },character(1))
    add('Evaluaciones no disponibles', 'Unavailable assessments',notes)
  }
  if ((!is.null(r$bootstrap) && x$effective_settings$n_boot < 1000L) ||
      (!is.null(r$invariance) && x$effective_settings$n_perm < 1000L))
    add('Precisi\u00f3n de remuestreo', 'Resampling precision',tr(
      'Se usaron recuentos reducidos de remuestreo; interpretar como diagn\u00f3stico y considerar la incertidumbre Monte Carlo antes del reporte definitivo.',
      'Reduced resampling counts were used; interpret as diagnostic and consider Monte Carlo uncertainty before final reporting.'))
  add('Detalle t\u00e9cnico', 'Technical detail',tr(
    'Los diagn\u00f3sticos completos se conservan en el objeto. Use print(fit, detail="audit") o psychomatic_report(fit, detail="audit") para consultarlos.',
    'Complete diagnostics remain in the object. Use print(fit, detail="audit") or psychomatic_report(fit, detail="audit") to inspect them.'))
  add('Referencias (APA 7)', 'References (APA 7)',.pm_network_report_references(x))
  out$sections <- sections
  if (destination != 'console') {
    types <- .pm_network_publication_types(x)
    out$figures <- Filter(function(figure)figure$type %in% types,out$figures)
  }
  out
}

.pm_network_report_references <- function(x) {
  refs <- .pm_reference_text(x$references)
  # Include the sources used by the displayed methods; the full bibliography is
  # retained on x and in the audit view.
  selected <- grepl('^Christensen,.*\\(2025\\)|^Golino,.*\\(202[01]\\)|^Christensen,.*\\(2024\\)|^Isvoranu,',refs)
  if (!is.null(x$results$uva)) selected <- selected | grepl('^Christensen,.*\\(2023\\)',refs)
  if (!is.null(x$results$bootstrap)) selected <- selected | grepl('^Christensen,.*\\(2021\\)|^Epskamp,',refs)
  if (!is.null(x$results$hierarchy)) selected <- selected | grepl('^Jim',refs)
  if (!is.null(x$results$invariance)) selected <- selected | grepl('^Benjamini,|^Jamison,|^van Borkulo,',refs)
  if (identical(x$results$wording$status,'completed')) selected <- selected | grepl('^Bulut,|^Garcia-Pardina,|^Faraci,',refs)
  refs[selected]
}
