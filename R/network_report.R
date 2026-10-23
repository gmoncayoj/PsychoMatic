.pm_network_audit_sections <- function(x, language = NULL, destination = "console") {
  language <- .pm_normalize_language(language, x$language)
  es <- language == "esp"; tr <- function(esp,eng) if(es)esp else eng
  r <- x$results; sections <- list()
  add <- function(name, table) {
    if (is.null(table) || !length(table)) return(invisible(NULL))
    if (is.character(table)) table <- data.frame(Text=table)
    if (is.matrix(table)) table <- .pm_matrix_table(table)
    if (!is.data.frame(table) || !nrow(table)) return(invisible(NULL))
    sections[[name]] <<- .pm_network_translate_table(table,es)
  }
  summary <- tr(
    sprintf("Se analizaron %d observaciones y %d de %d \u00edtems. EGA identific\u00f3 %d dimensiones. Estado del flujo: %s.",x$data_summary$n,x$data_summary$p_retained,x$data_summary$p_original,r$ega$n.dim,if(x$conclusion$status=="completed")"completado"else"parcial"),
    sprintf("The analysis used %d observations and %d of %d items. EGA identified %d dimensions. Workflow status: %s.",x$data_summary$n,x$data_summary$p_retained,x$data_summary$p_original,r$ega$n.dim,x$conclusion$status))
  add(tr("Resumen","Summary"),summary)
  seed_note <- if (identical(x$effective_settings$seed_source,"session_generated"))
    sprintf(tr("semilla de an\u00e1lisis %.0f generada desde el RNG de la sesi\u00f3n; sin valor fijo predeterminado",
      "analysis seed %.0f drawn from session RNG; no fixed default"),x$effective_settings$seed) else
    if(is.null(x$effective_settings$seed)) tr("sin semilla fija", "no fixed seed") else
      sprintf(tr("semilla %.0f", "seed %.0f"), x$effective_settings$seed)
  uva_method <- if(is.null(r$uva))tr("UVA no se ejecut\u00f3 (uva=NULL); no se evalu\u00f3 redundancia.",
    "UVA was not performed (uva=NULL); redundancy was not assessed.") else tr(
      sprintf("UVA (Christensen et al., 2023; wTO \u2265 %.2f) precedi\u00f3 al EGA; modo %s.",x$effective_settings$uva_cutoff,x$effective_settings$uva),
      sprintf("UVA (Christensen et al., 2023; wTO >= %.2f) preceded EGA; mode %s.",x$effective_settings$uva_cutoff,x$effective_settings$uva))
  add(tr("M\u00e9todo","Method"),paste(uva_method,tr(
    sprintf("Se estim\u00f3 una red con %s; %s. Las cargas revisadas siguen a Christensen et al. (2025); TEFI a Golino et al. (2021).",.pm_network_method(x$effective_settings,TRUE),seed_note),
    sprintf("A network used %s; %s. Revised loadings follow Christensen et al. (2025); TEFI follows Golino et al. (2021).",.pm_network_method(x$effective_settings),seed_note))))
  add(tr("Etapas y condiciones","Stages and conditions"),r$stages[intersect(c("stage","status","elapsed_seconds"),names(r$stages))])
  explained <- !is.na(r$stages$reason) & nzchar(r$stages$reason)
  add(tr("Motivos de etapas no evaluables o fallidas","Stage explanations"),r$stages[explained,c("stage","reason"),drop=FALSE])
  condition_columns <- function(t) intersect(c("stage","message"),names(t))
  if(nrow(x$errors)) add(tr("Errores de etapas","Stage errors"),unique(x$errors[condition_columns(x$errors)]))
  if(nrow(x$warnings)) add(tr("Advertencias","Warnings"),unique(x$warnings[condition_columns(x$warnings)]))
  if(!is.null(r$uva_selection))add(tr("Decisi\u00f3n UVA","UVA decision"),
    .pm_network_uva_decision_text(r$uva_selection$decision,es))
  if(!is.null(r$uva)) {
  add(tr("Selecci\u00f3n UVA","UVA selection"),data.frame(item=x$diagnostics$input$spec$items,
    retained=x$diagnostics$input$spec$items %in% r$retained_items))
  if(!is.null(r$uva_selection)) {
    s <- r$uva_selection
    add(tr("\u00cdtems eliminados por UVA","Items removed by UVA"),if(length(s$removed))
      data.frame(item=s$removed) else tr("No se elimin\u00f3 ning\u00fan \u00edtem.","No items were removed."))
    audit <- s$audit
    if(nrow(audit)) {
      audit$reason <- .pm_network_uva_reasons(audit,es)
      add(tr("Propuesta UVA y motivos de selecci\u00f3n","UVA proposal and selection reasons"),audit)
    }
  }
  pairs <- r$uva$wto$pairwise
  pairs <- pairs[is.finite(pairs$wto) & pairs$wto >= x$effective_settings$uva_cutoff,,drop=FALSE]
  add(tr("Pares redundantes","Redundant pairs"),if(nrow(pairs)) pairs else tr(
    sprintf("Ning\u00fan par alcanza el corte wTO \u2265 %.2f.",x$effective_settings$uva_cutoff),
    sprintf("No pairs meet the wTO cutoff >= %.2f.",x$effective_settings$uva_cutoff)))
  add(tr("Interpretaci\u00f3n de UVA","UVA interpretation"),tr(
    "UVA identifica posible dependencia local, no demuestra efectos de redacci\u00f3n ni menor calidad del \u00edtem eliminado. La decisi\u00f3n requiere revisar contenido y cobertura del constructo (Christensen et al., 2023).",
    "UVA identifies possible local dependence; it does not establish wording effects or poorer quality of the removed item. The decision requires review of content and construct coverage (Christensen et al., 2023)."))
  }
  add(tr("Cargas primarias revisadas","Revised primary loadings"),r$loadings$primary)
  add(tr("Cargas y cargas cruzadas","Loadings and cross-loadings"),r$loadings$matrix)
  add(tr("Centralidad","Centrality"),r$centrality)
  add(tr("Entrop\u00eda y conjuntos comparables","Entropy and comparison sets"),r$tefi)
  if(!is.null(r$wording)) {
    add(tr("Diagn\u00f3stico de redacci\u00f3n","Wording diagnostic"),r$wording$comparison)
    add(tr("Signos y claves declaradas","Signs and declared keys"),r$wording$candidates)
    if(!is.null(r$wording$stability)) {
      add(tr("Estabilidad riEGA","riEGA stability"),r$wording$stability$dimensions)
      add(tr("Estabilidad de \u00edtems riEGA","riEGA item stability"),r$wording$stability$items)
    }
  }
  if(!is.null(r$hierarchy)) {
    add(tr("M\u00e9todo jer\u00e1rquico","Hierarchy method"),tr(
      "hierEGA estima comunidades de orden inferior, puntuaciones de red y comunidades de orden superior (Jim\u00e9nez et al., 2025). Estas son evidencias exploratorias de redes. El AFC suplementario usa lavaan y la partici\u00f3n del EGA inicial de la muestra conjunta; no usa la partici\u00f3n inferior de hierEGA, que puede detectar otras comunidades. Su tabla no representa resultados de hierEGA ni confirmaci\u00f3n independiente. Con tres factores AFC de primer orden, el componente de segundo orden tiene cero grados de libertad: el mismo ajuste que el modelo correlacionado no demuestra superioridad. WLSMV no proporciona AIC/BIC. Una varianza negativa invalida un candidato, aunque su CFI sea mayor.",
      "hierEGA estimates lower-order communities, network scores and higher-order communities (Jim\u00e9nez et al., 2025). These are exploratory network results. Supplementary CFA uses lavaan and the initial pooled EGA partition; it does not use the hierEGA lower-order partition, which may have different communities. Its table is not output from hierEGA or independent confirmation. With three first-order CFA factors, the second-order component has zero degrees of freedom: equal fit with the correlated model does not establish superiority. WLSMV does not provide AIC/BIC. Negative variance makes a candidate inadmissible even with a higher CFI."))
    add(tr("Jerarqu\u00eda de redes (hierEGA)","Network hierarchy (hierEGA)"),r$hierarchy$network_structure)
    add(tr("Asignaciones jer\u00e1rquicas de red","Network hierarchy assignments"),r$hierarchy$network_assignments)
    interpretation <- r$hierarchy$network_interpretation
    if(es && !is.null(interpretation)) {
      pairs <- c("Based on lower ("="Seg\u00fan TEFI inferior (", ") and higher ("=") y superior (",
        ") TEFI, there is better fit for a "=") hay mejor ajuste relativo para una ",
        "correlated lower order factor structure than bifactor structure"="estructura inferior correlacionada que para una estructura bifactor",
        "bifactor structure than correlated lower order factor structure"="estructura bifactor que para una estructura inferior correlacionada")
      for(from in names(pairs))interpretation <- gsub(from,pairs[[from]],interpretation,fixed=TRUE)
    }
    add(tr("Interpretaci\u00f3n nativa exploratoria hierEGA","Native exploratory hierEGA interpretation"),interpretation)
    candidates <- r$hierarchy$candidates
    if(nrow(candidates))candidates$reason <- vapply(candidates$reason,function(reason) {
      codes <- trimws(strsplit(reason,";",fixed=TRUE)[[1]])
      translated <- c(post_check_failed=tr("Fall\u00f3 la comprobaci\u00f3n de admisibilidad","Admissibility post-check failed"),
        negative_variance=tr("Varianza negativa (caso Heywood)","Negative variance (Heywood case)"))
      mapped <- match(codes,names(translated));codes[!is.na(mapped)] <- translated[mapped[!is.na(mapped)]]
      paste(codes,collapse="; ")
    },character(1))
    add(tr("Candidatos AFC suplementarios","Supplementary CFA candidates"),candidates)
    if(!is.null(r$hierarchy$stability)) {
      add(tr("Estabilidad jer\u00e1rquica","Hierarchical stability"),r$hierarchy$stability$dimensions)
      add(tr("R\u00e9plicas jer\u00e1rquicas","Hierarchical replications"),r$hierarchy$stability$counts)
    }
    add(tr("Recomendaci\u00f3n exploratoria","Exploratory recommendation"),tr(
      paste("Candidato:",r$hierarchy$recommendation,".",.pm_network_reason_text(r$hierarchy$reason,TRUE),"Criterios: convergencia, admisibilidad, CFI \u2265 .95, RMSEA \u2264 .06, SRMR \u2264 .08; ventaja BIC \u226510 para preferir un modelo m\u00e1s complejo cuando es comparable. Los casos no identificables o sin ajuste aceptable se reportan como indeterminados."),
      paste("Candidate:",r$hierarchy$recommendation,".",r$hierarchy$reason)))
  }
  if(!is.null(r$bootstrap)) {
    b<-r$bootstrap
    bootstrap_method <- .pm_null_default(x$effective_settings$bootstrap_type,"resampling")
    add(tr("M\u00e9todo bootstrap","Bootstrap method"),paste(
      sprintf("bootEGA type = %s.",bootstrap_method),tr(
      "bootEGA estim\u00f3 la consistencia estructural; bootnet evalu\u00f3 precisi\u00f3n de aristas mediante remuestreo no param\u00e9trico y estabilidad de centralidad mediante eliminaci\u00f3n de casos (Christensen & Golino, 2021; Epskamp et al., 2018). El conjunto de \u00edtems del an\u00e1lisis se mantuvo fijo en las r\u00e9plicas.",
      "bootEGA assessed structural consistency; bootnet used nonparametric edge precision and case-drop centrality stability (Christensen & Golino, 2021; Epskamp et al., 2018). The analysis item set was held fixed across replications."),
      if(bootstrap_method=="parametric")tr(
        "Las r\u00e9plicas param\u00e9tricas son normales multivariadas continuas; no reproducen las categor\u00edas ordinales observadas. Con corr=polychoric, el estimador auto nativo trata las r\u00e9plicas continuas como continuas.",
        "Parametric replicates are continuous multivariate-normal data, not observed ordinal categories. With corr=polychoric, native auto correlation treats continuous replicates as continuous.") else tr(
        "Se remuestrearon participantes con reemplazo.","Participants were resampled with replacement.")))
    add(tr("R\u00e9plicas bootstrap","Bootstrap replications"),b$counts)
    add(tr("N\u00famero de dimensiones","Number of dimensions"),b$structural$frequency)
    add(tr("Estabilidad de \u00edtems","Item stability"),b$structural$items)
    add(tr("Consistencia estructural","Structural consistency"),b$structural$dimensions)
    add(tr("Precisi\u00f3n de aristas: IC percentil 95%","Edge precision: 95% percentile intervals"),b$edge_intervals)
    add(tr("Coeficiente CS","CS coefficient"),b$cs)
  }
  if(!is.null(r$invariance)) {
    i<-r$invariance
    add(tr("M\u00e9todo de comparaci\u00f3n","Comparison method"),tr(
      sprintf("Se aplic\u00f3 NCT (van Borkulo et al., 2023) y una prueba de cargas de EGA (Jamison et al., 2024), con %d permutaciones; resoluci\u00f3n m\u00ednima p = %.6f. BH se aplic\u00f3 una vez por familia, incluyendo todos los pares. Apoyo a la estructura configural original: %s. Las pruebas de cargas son exploratorias si falta este apoyo.",x$effective_settings$n_perm,i$p_resolution,if(i$configural_supported)"s\u00ed"else"no"),
      sprintf("NCT (van Borkulo et al., 2023) and EGA loading permutation tests (Jamison et al., 2024) used %d permutations; minimum p resolution = %.6f. BH was applied once per family across all pairs. Original configural structure supported: %s. Loading tests are exploratory when this support is absent.",x$effective_settings$n_perm,i$p_resolution,i$configural_supported)))
    add(tr("Grupos","Groups"),i$groups)
    add(tr("C\u00e1lculo de permutaciones de cargas","Loading permutation computation"),tr(
      "Se reestim\u00f3 la red con el estimador seleccionado en cada permutaci\u00f3n. Las cargas revisadas usan la partici\u00f3n original fija; Walktrap/LE evita consenso no utilizado y no selecciona una nueva estructura para esta prueba. NCT permanece secuencial; bootEGA, bootnet y las permutaciones de cargas usan n_cores.",
      "Networks were re-estimated with the selected estimator in every permutation. Revised loadings use the fixed original partition; Walktrap/LE avoids unused consensus and does not select a new structure for this test. NCT remains sequential; bootEGA, bootnet and loading permutations use n_cores."))
    add(tr("Familias BH","BH families"),i$families)
    add(tr("Estructura configural original","Original configural structure"),i$configural)
    add(tr("Pruebas globales NCT","NCT global tests"),i$global)
    add(tr("Diferencias de aristas","Edge differences"),i$edges)
    add(tr("Diferencias de centralidad","Centrality differences"),i$centrality)
    add(tr("Invarianza de cargas de red","Network loading invariance"),i$loadings)
  }
  add(tr("Interpretaci\u00f3n y l\u00edmites","Interpretation and limitations"),if(es) .pm_network_limitations_es() else x$diagnostics$limitations)
  add(tr("Referencias (APA 7)","References (APA 7)"),unique(x$references))
  figures <- if(destination=="console")list() else lapply(.pm_network_plot_types(x),function(type)
    list(type=type,title=.pm_network_figure_title(type,es),
      caption=.pm_network_figure_caption(type,es,x),
      draw=function(path) .pm_network_draw_file(x,type,path,language=language)))
  list(title=tr("An\u00e1lisis de redes psicom\u00e9tricas","Psychometric network analysis"),language=language,sections=sections,figures=figures)
}

.pm_network_translate_table <- function(t,es) {
  if(!es)return(t)
  semantic_columns <- names(t) %in% c("status", "magnitude", "metric", "interpretation", "model", "family", "test", "criterion","order")
  dict<-c(item="\u00cdtem",items="\u00cdtems",dimension="Dimensi\u00f3n",loading="Carga",magnitude="Magnitud",strength="Fuerza",
    elapsed_seconds="Segundos",order="Orden",separately_estimated="Estimado por separado",
    lower_dimension="Dimensi\u00f3n inferior",higher_dimension="Dimensi\u00f3n superior",
    df="Grados de libertad",higher_order_df="GL del componente de segundo orden",
    betweenness="Intermediaci\u00f3n",closeness="Cercan\u00eda",expected_influence="Influencia esperada",connected="Conectada",
    model="Modelo",index="\u00cdndice",value="Valor",comparison_set="Conjunto comparable",comparable="Comparable",
    stage="Etapa",status="Estado",reason="Motivo",retained="Retenido",dimensions="Dimensiones",
    proposed_remove="Eliminaci\u00f3n propuesta",representative="Representante conservado si se acepta",
    redundant_set="Conjunto redundante",criterion="Criterio",removed_score="Puntuaci\u00f3n del candidato",
    retained_score="Puntuaci\u00f3n del representante",tied="Empate",removed="Eliminado",
    same_partition="Misma partici\u00f3n",negative_links="Asociaciones negativas",declared_reverse="Clave inversa declarada",
    estimator="Estimador",converged="Convergencia",admissible="Admisible",requested="Solicitadas",valid="V\u00e1lidas",
    failed="Fallidas",stability="Estabilidad",structural_consistency="Consistencia estructural",
    average_item_stability="Estabilidad media de \u00edtems",edge="Arista",mean="Media",lower="L\u00edmite inferior",upper="L\u00edmite superior",
    metric="M\u00e9trica",interpretation="Interpretaci\u00f3n",group="Grupo",comparison="Comparaci\u00f3n",test="Prueba",
    statistic="Estad\u00edstico",family="Familia",difference="Diferencia",scope="Alcance",size="Tama\u00f1o",
    valid_permutations="Permutaciones v\u00e1lidas",package="Paquete",version="Versi\u00f3n",message="Mensaje",code="C\u00f3digo")
  idx<-match(names(t),names(dict));names(t)[!is.na(idx)]<-unname(dict[idx[!is.na(idx)]])
  values<-c(completed="completada",skipped="omitida",interrupted="interrumpida",failed="fallida",below_small="inferior a peque\u00f1a",small="peque\u00f1a",moderate="moderada",large="grande",
    lower_order="Orden inferior",higher_order="Orden superior",
    unstable="inestable",caution="cautela",preferable="preferible",not_evaluable="no evaluable",
    correlated="correlacionado",second_order="segundo orden",bifactor="bifactor",indeterminate="indeterminado",
    global="globales",edges="aristas",centrality="centralidad",loadings="cargas",network_structure="estructura de red",global_strength="fuerza global",
    min_external_max_wto="Menor wTO m\u00e1ximo externo",max_within_mean_wto="Mayor wTO medio interno")
  for(j in which(semantic_columns))if(is.character(t[[j]])){
    idx<-match(t[[j]],names(values));t[[j]][!is.na(idx)]<-unname(values[idx[!is.na(idx)]])
  }
  if("Motivo" %in% names(t)) t[["Motivo"]] <- .pm_network_reason_text(t[["Motivo"]],TRUE)
  t
}

.pm_network_reason_text <- function(text,es=FALSE) {
  if(!es)return(text)
  reasons <- c(
    "Skipped after user interruption."="Omitida tras la interrupci\u00f3n del usuario.",
    "The installed riEGA only passes its default seven-category threshold to the initial correlation estimator; riEGA was skipped to avoid mixing correlation methods for items with more than seven categories."=
      "La versi\u00f3n instalada de riEGA usa su umbral predeterminado de siete categor\u00edas al estimar las correlaciones iniciales; se omiti\u00f3 riEGA para evitar mezclar m\u00e9todos con \u00edtems de m\u00e1s de siete categor\u00edas.",
    "UVA was not requested (uva=NULL); redundancy was not assessed." =
      "UVA no se solicit\u00f3 (uva=NULL); no se evalu\u00f3 redundancia.",
    "The riEGA residual network has no estimable communities or unassigned items; its partition TEFI and structural stability are not evaluable. This does not prove wording effects." =
      "La red residual riEGA no tiene comunidades estimables o presenta \u00edtems sin asignar; no se pueden evaluar su TEFI de partici\u00f3n ni su estabilidad estructural. Esto no demuestra efectos de redacci\u00f3n.",
    "riEGA stability was not computed because the empirical residual partition is not evaluable." =
      "No se calcul\u00f3 la estabilidad riEGA porque la partici\u00f3n residual emp\u00edrica no es evaluable.",
    "Requires at least three dimensions with three indicators each." =
      "La comparaci\u00f3n AFC requiere al menos tres dimensiones con tres indicadores cada una despu\u00e9s de UVA. No se selecciona un modelo jer\u00e1rquico o bifactor cuando faltan indicadores.",
    "No candidate met all declared convergence, admissibility and fit heuristics." =
      "Ning\u00fan candidato cumpli\u00f3 todos los criterios declarados de convergencia, admisibilidad y ajuste.")
  index <- match(text,names(reasons))
  text[!is.na(index)] <- unname(reasons[index[!is.na(index)]])
  text
}

.pm_network_figure_title <- function(type,es) {
  eng<-c(network="EGA network and edge weights",uva="UVA redundancy",loadings="Primary item-dimension loadings",
    loadings_matrix="Revised loadings and cross-loadings",centrality="Centrality",
    tefi="Entropy comparison sets",wording="Wording diagnostic",hierarchy="Hierarchical structure",
    ega_stability="EGA dimensionality and item stability",dimensions="Bootstrap dimensionality",item_stability="Item stability",dimension_stability="Structural consistency",
    edge_accuracy="Edge precision",centrality_stability="Centrality stability",group_networks="Group networks",
    edge_differences="Edge differences (BH)",loading_differences="Loading differences")
  esp<-c(network="Red EGA y pesos de aristas",uva="Redundancia UVA",loadings="Cargas primarias de \u00edtem-dimensi\u00f3n",
    loadings_matrix="Cargas revisadas y cargas cruzadas",centrality="Centralidad",
    tefi="Conjuntos de entrop\u00eda",wording="Diagn\u00f3stico de redacci\u00f3n",hierarchy="Estructura jer\u00e1rquica",
    ega_stability="Dimensionalidad EGA y estabilidad de \u00edtems",dimensions="Dimensiones bootstrap",item_stability="Estabilidad de \u00edtems",dimension_stability="Consistencia estructural",
    edge_accuracy="Precisi\u00f3n de aristas",centrality_stability="Estabilidad de centralidad",group_networks="Redes por grupo",
    edge_differences="Diferencias de aristas (BH)",loading_differences="Diferencias de cargas")
  unname(if(es)esp[type]else eng[type])
}

.pm_network_figure_caption <- function(type,es,x=NULL) {
  tr<-function(esp,eng)if(es)esp else eng
  if(!is.null(x) && identical(x$effective_settings$model,"TMFG") && type %in% c('network','ega_stability'))
    return(tr('Las aristas muestran correlaciones filtradas por TMFG, no correlaciones parciales ni cargas. Los colores indican comunidades; las aristas negativas son discontinuas. Si hay un panel de estabilidad, muestra asignaciones reproducidas en r\u00e9plicas v\u00e1lidas.',
      'Edges show TMFG-filtered correlations, not partial correlations or loadings. Colors indicate communities; negative edges are dashed. A stability panel, when present, shows reproduced assignments in valid replicates.'))
  if(type=='ega_stability') return(tr(
    'Panel izquierdo: red EGA, con correlaciones parciales como aristas. Panel derecho: estabilidad de cada \u00edtem en r\u00e9plicas v\u00e1lidas, manteniendo fijo el conjunto de \u00edtems. Los colores representan comunidades conjuntas. No se eliminan \u00edtems por su estabilidad.',
    'Left panel: EGA network with partial-correlation edges. Right panel: item assignment stability over valid replicates with a fixed item set. Colors represent pooled communities. Stability does not automatically remove items.'))
  if(type=="network")return(tr("Los nodos identifican \u00edtems. Las aristas representan correlaciones parciales regularizadas, no cargas; las etiquetas num\u00e9ricas se muestran en redes dispersas. Los colores indican comunidades; las aristas negativas son discontinuas.",
    "Nodes identify items. Edges encode regularized partial correlations, not loadings; numeric labels appear in sparse graphs. Colors indicate communities. Negative edges are dashed."))
  if(type=="uva")return(tr("El mapa representa wTO; los pares por encima del umbral declarado requieren revisi\u00f3n de redundancia y cobertura de contenido.",
    "The map shows wTO; pairs above the declared threshold require redundancy and content-coverage review."))
  if(type=="loadings")return(tr("Las aristas \u00edtem-dimensi\u00f3n muestran cargas primarias revisadas no nulas, con su signo. Los cuadrados representan dimensiones, no \u00edtems ni factores AFC. Consulte las cargas cruzadas en el mapa y las tablas.",
    "Item-dimension edges show nonzero signed revised primary loadings. Squares represent dimensions, not items or CFA factors. Consult cross-loadings in the map and tables."))
  if(type=="loadings_matrix")return(tr("Se conservan signos y cargas cruzadas. Magnitudes orientativas: .20 peque\u00f1a, .35 moderada, .50 grande.",
    "Signs and cross-loadings are retained. Approximate magnitudes: .20 small, .35 moderate, .50 large."))
  if(type=="tefi")return(tr("Comparar \u00fanicamente valores del mismo conjunto y la misma matriz; un valor menor indica mejor ajuste relativo. genTEFI y TEFI residual tienen conjuntos distintos.",
    "Compare only within the same set and matrix; lower values indicate better relative fit. genTEFI and residual TEFI use separate sets."))
  if(type=="hierarchy")return(tr('Panel de orden inferior: nodos de \u00edtems. Panel de orden superior: nodos C que representan puntuaciones de comunidades inferiores. C en las leyendas identifica comunidades propias de cada orden. No es un diagrama AFC; no comparar pesos entre \u00f3rdenes.',
    'Lower-order panel: item nodes. Higher-order panel: C nodes representing lower-community network scores. C in legends identifies communities within each order. This is not a CFA diagram; do not compare weights across orders.'))
  if(type=="edge_accuracy")return(tr("Intervalos bootstrap percentiles del 95%; no constituyen pruebas de significaci\u00f3n de aristas regularizadas.",
    "95% bootstrap percentile intervals; these are not significance tests of regularized edges."))
  if(type=="centrality_stability")return(tr("Correlaciones con la centralidad original al eliminar casos; interpretar junto con CS (.25 m\u00ednimo, .50 preferible).",
    "Correlations with original centrality as cases are dropped; interpret with CS (.25 minimum, .50 preferable)."))
  if(type %in% c("item_stability","dimension_stability","dimensions"))return(tr(
    "Frecuencias y proporciones calculadas sobre r\u00e9plicas v\u00e1lidas, manteniendo fijo el conjunto de \u00edtems reportado.",
    "Frequencies and proportions are computed from valid replicates, holding the reported item set fixed."))
  if(type=="group_networks")return(tr("Los grupos comparten \u00edtems, posiciones y escala de pesos. Los colores representan comunidades del EGA conjunto. En una comparaci\u00f3n de dos grupos, bordes naranjas gruesos indican diferencias de cargas con p BH < .05; bordes grises, sin diferencia detectada; cuadrados, pruebas no evaluables. Sin apoyo configural, estas marcas son exploratorias. La no significaci\u00f3n no demuestra equivalencia. Aristas negativas discontinuas.",
    "Groups share items, positions and edge-weight scale. Colors identify pooled EGA communities. In a two-group comparison, thick orange borders flag loading differences with p BH < .05; gray borders indicate no difference detected; squares indicate unavailable tests. Without configural support, these flags are exploratory. Non-significance does not establish equivalence. Negative edges are dashed."))
  if(type=="edge_differences")return(tr("Diferencias firmadas de aristas del primer grupo menos el segundo; solo se muestran las que cumplen p BH < .05. Una figura vac\u00eda no demuestra equivalencia.",
    "Signed edge differences, first group minus second; only p BH < .05 edges are shown. An empty figure does not establish equivalence."))
  if(type=="loading_differences")return(tr("Diferencias firmadas del primer grupo menos el segundo. Puntos rellenos: p BH < .05; vac\u00edos: p BH \u2265 .05; cruces: pruebas no evaluables. No son intervalos de confianza. Interpretar como exploratorio si no hay apoyo configural.",
    "Signed differences, first group minus second. Filled points: p BH < .05; open points: p BH >= .05; crosses: unavailable tests. These are not confidence intervals. Interpret as exploratory without configural support."))
  tr("Diagn\u00f3stico exploratorio; interpretar junto con las tablas de ajuste, estabilidad y limitaciones.",
    "Exploratory diagnostic; interpret with fit, stability and limitation tables.")
}

.pm_network_limitations_es <- function() c(
  "EGA y los modelos jer\u00e1rquicos son exploratorios. Ajustar en la misma muestra no constituye confirmaci\u00f3n independiente.",
  "UVA selecciona \u00edtems por redundancia estad\u00edstica; la cobertura de contenido y los representantes retenidos requieren revisi\u00f3n.",
  "riEGA usa respuestas originales sin recodificar. Un intercepto aleatorio ajustado no prueba efectos de redacci\u00f3n ni identifica claves inversas.",
  "Cargas revisadas: |carga| <.20 inferior a peque\u00f1a; .20 peque\u00f1a; .35 moderada; .50 grande. Son gu\u00edas de magnitud, no criterios de eliminaci\u00f3n.",
  "TEFI expresa ajuste relativo sobre la misma matriz e \u00edtems. genTEFI y TEFI residual pertenecen a conjuntos de comparaci\u00f3n distintos.",
  "Fuerza suma pesos absolutos. Las distancias usan el inverso de pesos absolutos; cercan\u00eda no se define para redes desconectadas.",
  "CS <.25 es insuficiente; \u2265.50 es preferible. Intermediaci\u00f3n y cercan\u00eda suelen ser inestables.",
  "NCT compara estructura y fuerza global de redes. Las pruebas de cargas EGA eval\u00faan invarianza m\u00e9trica de red; ninguna demuestra invarianza escalar de AFC.",
  "BH controla la tasa de falsos descubrimientos dentro de familias declaradas que incluyen todos los pares. No significaci\u00f3n no equivale a equivalencia.",
  "Pocas r\u00e9plicas solo permiten diagn\u00f3sticos. Predeterminados para reporte: 1000 bootstrap y 5000 permutaciones; documentar incertidumbre Monte Carlo.",
  "La eliminaci\u00f3n completa es expl\u00edcita; las correlaciones por pares pueden producir matrices no definidas positivas. Evaluar el mecanismo de datos faltantes.")
