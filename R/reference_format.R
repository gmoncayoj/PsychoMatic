# Presentation-only APA runs. Raw bibliography strings are never altered.
.pm_reference_runs <- function(reference) {
  spans <- list()
  # The journal/volume tail is unambiguous in the package's APA article entries.
  pattern <- '^(.+\\([0-9]{4}[a-z]?\\)\\. .+\\. )([^\\.]+), ([0-9]+)(?=\\(|,|\\.)'
  match <- regexec(pattern, reference, perl=TRUE)
  parts <- regmatches(reference, match)[[1]]
  if (length(parts)) {
    start <- nchar(parts[2]) + 1L
    spans <- list(c(start, start + nchar(parts[3]) - 1L),
      c(start + nchar(parts[3]) + 2L, start + nchar(parts[3]) + 1L + nchar(parts[4])))
  } else {
    # Explicit titles avoid guessing whether a non-journal citation is a book,
    # chapter, software package or website. Includes books used in CFA/reliability.
    titles <- c('Confirmatory factor analysis for applied research',
      'Principles and practice of structural equation modeling',
      'Test theory: A unified treatment',
      'Latent variables analysis: Applications for developmental research',
      'Structural equation modeling: Present and future')
    for (title in titles) {
      start <- regexpr(title, reference, fixed=TRUE)[1]
      if (start > 0L) spans <- list(c(start, start + nchar(title) - 1L))
    }
  }
  chars <- strsplit(reference, '', fixed=TRUE)[[1]]
  italic <- rep(FALSE, length(chars))
  for (span in spans) italic[seq.int(span[1],span[2])] <- TRUE
  if (!length(chars)) return(data.frame(text=character(),italic=logical()))
  groups <- cumsum(c(TRUE, diff(italic) != 0))
  data.frame(text=vapply(split(chars,groups),paste,collapse='',FUN.VALUE=character(1)),
    italic=italic[!duplicated(groups)],row.names=NULL)
}

.pm_reference_ansi <- function() {
  option <- getOption('PsychoMatic.ansi')
  if (!is.null(option)) return(isTRUE(option))
  interactive() && !nzchar(Sys.getenv('NO_COLOR')) && Sys.getenv('TERM') != 'dumb' &&
    (nzchar(Sys.getenv('RSTUDIO')) || isatty(stdout()))
}

.pm_reference_console <- function(reference, ansi=.pm_reference_ansi(), width=NULL) {
  runs <- .pm_reference_runs(reference)
  render <- function(text, italic) paste0(ifelse(italic, paste0('\033[3m',text,'\033[23m'),text),collapse='')
  if (!ansi) return(if(is.null(width)) reference else paste(strwrap(reference,width=width),collapse='\n'))
  if (is.null(width)) return(render(runs$text,runs$italic))
  mask <- rep(runs$italic,nchar(runs$text))
  lines <- strwrap(reference,width=width)
  offset <- 0L
  styled <- vapply(lines,function(line) {
    start <- regexpr(line,substring(reference,offset+1L),fixed=TRUE)[1]
    if (start < 1L) return(line)
    start <- start + offset
    offset <<- start + nchar(line) - 1L
    chars <- strsplit(line,'',fixed=TRUE)[[1]]
    italics <- mask[seq.int(start,offset)]
    groups <- cumsum(c(TRUE,diff(italics)!=0))
    render(vapply(split(chars,groups),paste,collapse='',FUN.VALUE=character(1)),italics[!duplicated(groups)])
  },character(1))
  paste(styled,collapse='\n')
}

.pm_reference_section <- function(name) grepl('references|referencias',name,ignore.case=TRUE)

.pm_reference_columns <- function(table) intersect(c('reference','Reference','Referencia'),names(table))

.pm_reference_flextable <- function(ft, table) {
  for (column in .pm_reference_columns(table)) for (i in seq_len(nrow(table))) {
    if (is.na(table[[column]][i])) next
    runs <- .pm_reference_runs(as.character(table[[column]][i]))
    if (!any(runs$italic)) next
    chunks <- lapply(seq_len(nrow(runs)),function(j)
      flextable::as_chunk(runs$text[j],props=officer::fp_text_lite(italic=runs$italic[j])))
    ft <- flextable::compose(ft,i=i,j=column,value=do.call(flextable::as_paragraph,chunks))
  }
  ft
}

.pm_reference_word <- function(doc, reference, style='PsychoMatic Body') {
  runs <- .pm_reference_runs(reference)
  chunks <- lapply(seq_len(nrow(runs)),function(i)
    officer::ftext(runs$text[i],officer::fp_text_lite(italic=runs$italic[i])))
  officer::body_add_fpar(doc,do.call(officer::fpar,chunks),style=style)
}

.pm_reference_xml_escape <- function(text) {
  text <- gsub('&','&amp;',text,fixed=TRUE)
  text <- gsub('<','&lt;',text,fixed=TRUE)
  gsub('>','&gt;',text,fixed=TRUE)
}

.pm_reference_excel <- function(wb, references) {
  # openxlsx has no public rich-text writer. Replace only exact shared-string
  # entries after all cells have been written, preserving indices and counts.
  for (reference in unique(references)) {
    runs <- .pm_reference_runs(reference)
    if (!any(runs$italic)) next
    plain <- paste0('<si><t xml:space="preserve">',.pm_reference_xml_escape(reference),'</t></si>')
    index <- which(vapply(wb$sharedStrings,identical,logical(1),plain))
    if (!length(index)) next
    rich <- paste0('<si>',paste0('<r><rPr>',ifelse(runs$italic,'<i/>',''),
      '</rPr><t xml:space="preserve">',.pm_reference_xml_escape(runs$text),'</t></r>',collapse=''),'</si>')
    for (i in index) wb$sharedStrings[[i]] <- rich
  }
  invisible(wb)
}
