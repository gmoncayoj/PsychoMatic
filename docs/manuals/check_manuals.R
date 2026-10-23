# Structural parity gate for the complete translated reference, independent of PDF.
check_manuals <- function(root = ".", write_manifest = TRUE) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  canonical <- file.path(root, "man")
  translated <- file.path(root, "docs", "manuals", "es", "man")
  files <- sort(list.files(canonical, pattern = "[.]Rd$"))
  if (!identical(files, sort(list.files(translated, pattern = "[.]Rd$"))))
    stop("Translated Rd topic set differs: missing ",
         paste(setdiff(files, list.files(translated)), collapse = ", "))
  text <- function(x) paste(unlist(x, use.names = FALSE), collapse = "")
  block <- function(rd, tag) {
    idx <- which(vapply(rd, function(x) identical(attr(x, "Rd_tag"), tag), logical(1)))
    lapply(rd[idx], text)
  }
  arg_names <- function(rd) {
    args <- Filter(function(x) identical(attr(x, "Rd_tag"), "\\arguments"), rd)
    unlist(lapply(args, function(x) vapply(Filter(function(y)
      identical(attr(y, "Rd_tag"), "\\item"), x), function(y) text(y[[1L]]), "")),
      use.names = FALSE)
  }
  aliases <- character()
  rows <- lapply(files, function(f) {
    en_path <- file.path(canonical, f)
    es_path <- file.path(translated, f)
    en <- tools::parse_Rd(en_path, encoding = "UTF-8")
    es <- tools::parse_Rd(es_path, encoding = "UTF-8")
    for (tag in c("\\name", "\\alias", "\\usage", "\\keyword", "\\docType"))
      if (!identical(block(en, tag), block(es, tag)))
        stop("Translation changed ", tag, " in ", f)
    if (!identical(arg_names(en), arg_names(es))) stop("Argument names differ in ", f)
    # Parse all code, including dontrun templates; comments may be translated.
    en_ex <- tools::Rd2ex(en, out = tempfile(), commentDontrun = FALSE)
    es_ex <- tools::Rd2ex(es, out = tempfile(), commentDontrun = FALSE)
    en_code <- if (file.exists(en_ex)) parse(file = en_ex, keep.source = FALSE) else expression()
    es_code <- if (file.exists(es_ex)) parse(file = es_ex, keep.source = FALSE) else expression()
    unlink(c(en_ex, es_ex))
    if (!identical(en_code, es_code)) stop("Executable examples differ in ", f)
    aliases <<- c(aliases, unlist(block(en, "\\alias"), use.names = FALSE))
    data.frame(topic = text(block(en, "\\name")), file = f,
               canonical_md5 = unname(tools::md5sum(en_path)),
               spanish_md5 = unname(tools::md5sum(es_path)),
               stringsAsFactors = FALSE)
  })
  exports <- sub('^export\\((.*)\\)$', '\\1',
                 grep('^export\\(', readLines(file.path(root, "NAMESPACE")), value = TRUE))
  missing <- setdiff(exports, aliases)
  if (length(missing)) stop("Undocumented exports: ", paste(missing, collapse = ", "))
  manifest <- do.call(rbind, rows)
  if (write_manifest) utils::write.csv(manifest,
      file.path(root, "docs", "manuals", "coverage.csv"), row.names = FALSE,
      fileEncoding = "UTF-8")
  message(nrow(manifest), " topics; ", length(exports),
          " public functions; identical API and executable examples in EN/ES.")
  invisible(manifest)
}
if (sys.nframe() == 0L) check_manuals()
