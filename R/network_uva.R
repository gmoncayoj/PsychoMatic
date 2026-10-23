# Reconstruct the documented native selection scores without making a second
# item-selection decision. Native keep/remove is authoritative; disagreement
# with the documented rule stops reduction rather than inventing a rationale.
.pm_network_uva_audit <- function(native) {
  empty <- data.frame(proposed_remove=character(), representative=character(),
    redundant_set=character(), criterion=character(), removed_score=double(),
    retained_score=double(), tied=logical(), removed=logical(), reason=character())
  if (!length(native$redundant)) return(empty)
  wto <- native$wto$matrix
  rows <- lapply(seq_along(native$redundant), function(i) {
    members <- unique(c(names(native$redundant)[i], native$redundant[[i]]))
    if (length(members)<2L || any(!members %in% rownames(wto)))
      stop('UVA returned an unsupported redundancy set; no reduction was applied.',call.=FALSE)
    if (length(members)==2L) {
      others <- setdiff(colnames(wto),members)
      if (!length(others)) stop('UVA selection requires items outside a redundant pair.',call.=FALSE)
      scores <- apply(wto[members,others,drop=FALSE],1L,max,na.rm=TRUE)
      selected <- which.min(scores)
      criterion <- 'min_external_max_wto'
    } else {
      # EGAnet includes the zero diagonal in the within-set mean.
      scores <- colMeans(wto[members,members,drop=FALSE],na.rm=TRUE)
      selected <- which.max(scores)
      criterion <- 'max_within_mean_wto'
    }
    if (any(!is.finite(scores)) || length(selected)!=1L)
      stop('UVA selection scores are not evaluable; no reduction was applied.',call.=FALSE)
    representative <- members[selected]
    removed <- members[-selected]
    if (!representative %in% native$keep_remove$keep ||
        any(!removed %in% native$keep_remove$remove))
      stop('Native UVA selection differs from its documented rule; inspect before reducing.',call.=FALSE)
    data.frame(proposed_remove=removed, representative=representative,
      redundant_set=paste(members,collapse=', '),criterion=criterion,
      removed_score=unname(scores[-selected]),retained_score=unname(scores[selected]),
      tied=sum(scores==scores[selected])>1L,removed=FALSE,reason='',stringsAsFactors=FALSE)
  })
  audit <- do.call(rbind,rows)
  if (!setequal(audit$proposed_remove,native$keep_remove$remove) ||
      any(audit$representative %in% audit$proposed_remove))
    stop('UVA returned conflicting selection sets; no reduction was applied.',call.=FALSE)
  audit$reason <- .pm_network_uva_reasons(audit,FALSE)
  audit
}

.pm_network_uva_reasons <- function(audit, es=FALSE) {
  vapply(seq_len(nrow(audit)),function(i) {
    a <- audit[i,]
    rule <- if(a$criterion=='min_external_max_wto') {
      if(es)'menor wTO m\u00e1ximo con los dem\u00e1s \u00edtems (excluido el par)' else
        'lowest maximum wTO to other items (excluding the pair)'
    } else {
      if(es)'mayor wTO medio dentro del conjunto redundante (incluye diagonal cero)' else
        'highest mean within-set wTO (including the zero diagonal)'
    }
    pattern <- if(es)if(isTRUE(a$removed))'Se elimin\u00f3 %s y se conserv\u00f3 %s por el %s: %.6f frente a %.6f, respectivamente.'else
      'Se propone eliminar %s y conservar %s por el %s: %.6f frente a %.6f, respectivamente.' else
      if(isTRUE(a$removed))'Removed %s and retained %s using the %s: %.6f versus %.6f, respectively.'else
        'Propose removing %s and retaining %s using the %s: %.6f versus %.6f, respectively.'
    text <- sprintf(pattern,
      a$proposed_remove,a$representative,rule,a$removed_score,a$retained_score)
    if(a$tied)text <- paste(text,if(es)
      'Empate: EGAnet conserva el primer candidato en su orden nativo; no implica mejor calidad.' else
      'Exact tie: EGAnet retains the first candidate in native order; this does not imply better quality.')
    text
  },character(1))
}

.pm_network_uva_select <- function(native,items,mode,cutoff) {
  if(is.null(mode))return(list(requested_mode=NULL,effective_mode=NULL,
    decision='not_requested',retained=items,removed=character(),pairs=NULL,
    audit=NULL,native_keep_remove=NULL))
  pairs <- native$wto$pairwise
  pairs <- pairs[is.finite(pairs$wto) & pairs$wto>=cutoff,,drop=FALSE]
  audit <- .pm_network_uva_audit(native)
  candidate_remove <- unique(audit$proposed_remove)
  decision <- if(!length(candidate_remove))'no_candidates' else
    if(mode=='removed')'explicit_removed'else'explicit_retain'
  remove <- if(mode=='removed')candidate_remove else character()
  audit$removed <- audit$proposed_remove %in% remove
  audit$reason <- .pm_network_uva_reasons(audit,FALSE)
  list(requested_mode=mode,effective_mode=mode,
    decision=decision,retained=items[!items %in% remove],removed=items[items %in% remove],
    pairs=pairs,audit=audit,native_keep_remove=native$keep_remove)
}

.pm_network_uva_decision_text <- function(decision,es=FALSE) {
  eng <- c(not_requested='UVA was not performed (uva=NULL); redundancy was not assessed and no items were removed.',
    explicit_removed='Explicit uva="removed": the justified native removal proposal was applied.',
    accepted='The user entered 1: the proposed reduction was accepted.',
    declined='The user entered 0: all items were retained.',
    empty_retained='Empty input/end of input: all items were retained.',
    explicit_reduce='Explicit uva="reduce": the proposed reduction was applied.',
    explicit_retain='Explicit uva="diagnostic": all items were retained.',
    noninteractive_retained='No interactive console: all items were retained. Use uva="reduce" for authorized reduction in scripts.',
    no_candidates='No redundant items were proposed for removal; all items were retained.')
  esp <- c(not_requested='UVA no se ejecut\u00f3 (uva=NULL); no se evalu\u00f3 redundancia ni se eliminaron \u00edtems.',
    explicit_removed='uva="removed" expl\u00edcito: se aplic\u00f3 la eliminaci\u00f3n nativa con motivos documentados.',
    accepted='El usuario ingres\u00f3 1: acept\u00f3 la reducci\u00f3n propuesta.',
    declined='El usuario ingres\u00f3 0: se conservaron todos los \u00edtems.',
    empty_retained='Respuesta vac\u00eda/fin de entrada: se conservaron todos los \u00edtems.',
    explicit_reduce='uva="reduce" expl\u00edcito: se aplic\u00f3 la reducci\u00f3n propuesta.',
    explicit_retain='uva="diagnostic" expl\u00edcito: se conservaron todos los \u00edtems.',
    noninteractive_retained='Sin consola interactiva: se conservaron todos los \u00edtems. Use uva="reduce" para autorizar la reducci\u00f3n en scripts.',
    no_candidates='No se propusieron \u00edtems redundantes para eliminar; se conservaron todos.')
  unname(if(es)esp[decision]else eng[decision])
}
