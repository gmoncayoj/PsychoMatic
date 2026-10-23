# Standard R reference manuals; Spanish is an editorial translation, not an
# additional installed help tree. Run from the repository root with R + LaTeX.
args <- commandArgs(trailingOnly = TRUE)
english_only <- "--english-only" %in% args
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
stopifnot(file.exists(file.path(root, "DESCRIPTION")), dir.exists(file.path(root, "man")))
source(file.path(root, "docs", "manuals", "check_manuals.R"))
if (!english_only) check_manuals(root)
outdir <- file.path(root, "docs", "manuals")
work <- file.path(outdir, "work")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
version <- unname(read.dcf(file.path(root, "DESCRIPTION"))[1L, "Version"])
edition <- Sys.getenv("PSYCHOMATIC_MANUAL_DATE", "2026-10-20")
stopifnot(!is.na(as.Date(edition)))
Sys.setenv(R_PAPERSIZE = "a4")

if (!english_only) {
  desc <- read.dcf(file.path(root, "DESCRIPTION"))
  desc[1L, "Title"] <- "Flujos de trabajo psicométricos automatizados e informes"
  desc[1L, "Description"] <- paste(
    "Automatiza flujos psicométricos para investigadores aplicados: descriptivos",
    "de ítems, correlaciones, análisis factorial exploratorio y confirmatorio,",
    "modelos exploratorios de medición con ecuaciones estructurales, fiabilidad,",
    "redes psicométricas, invarianza multigrupo y optimización por alineación.",
    "Incluye procedimientos de análisis paralelo (Horn, 1965,",
    "<doi:10.1007/BF02289447>), diagnóstico de normalidad multivariante (Mardia,",
    "1970, <doi:10.1093/biomet/57.3.519>), cambios de ajuste en invarianza (Chen,",
    "2007, <doi:10.1080/10705510701301834>) y alineación (Asparouhov y Muthen,",
    "2014, <doi:10.1080/10705511.2014.919210>), entre otros. Devuelve objetos",
    "estructurados de R y exporta informes bilingües para documentar el análisis.")
  write.dcf(desc, file.path(outdir, "es", "DESCRIPTION"))
}

build_one <- function(language, input) {
  build <- file.path(work, tolower(language))
  # Rd2pdf recreates this explicitly bounded scratch directory on every run.
  stopifnot(startsWith(build, paste0(work, "/")))
  pdf <- file.path(outdir, paste0("PsychoMatic_", version, "_", language, ".pdf"))
  logfile <- file.path(work, paste0(tolower(language), "-rd2pdf.log"))
  status <- system2(file.path(R.home("bin"), "R"), c("CMD", "Rd2pdf",
      "--internals", "--no-preview", "--no-clean", "--force",
      shQuote(paste0("--build-dir=", build)), shQuote(paste0("--output=", pdf)),
      shQuote(input)), stdout = logfile, stderr = logfile)
  if (status != 0L) stop("Rd2pdf failed for ", language, "; see ", logfile)
  tex <- file.path(build, "Rd2.tex")
  lines <- readLines(tex, encoding = "UTF-8", warn = FALSE)
  # An explicit edition date avoids changing the displayed date across rebuilds.
  lines <- gsub("\\today", edition, lines, fixed = TRUE)
  # The default tight paragraph setting can clip long lists of API identifiers.
  # Permit wider interword spacing while retaining native Rd section/layout rules.
  at <- match("\\begin{document}", lines)
  lines <- append(lines, "\\sloppy", after = at)
  if (language == "ES") {
    labels <- c(Arguments = "Argumentos", Author = "Autoría", Description = "Descripción",
      Details = "Detalles", Examples = "Ejemplos", Note = "Nota", References = "Referencias",
      SeeAlso = "Véase también", Format = "Formato", Source = "Fuente",
      Usage = "Uso", Value = "Valor")
    preamble <- c("% Localized headings; installed R/Rd.sty is never modified.",
      vapply(names(labels), function(env) paste0("\\renewenvironment{", env,
        "}{\\begin{Rdsection}{", labels[[env]], "}}{\\end{Rdsection}}"), ""),
      "\\renewcommand{\\indexname}{Índice}")
    at <- match("\\begin{document}", lines)
    lines <- append(lines, preamble, after = at - 1L)
    lines <- gsub("\\Rdcontents{Contents}", "\\Rdcontents{Contenido}", lines, fixed = TRUE)
    lines <- gsub("Package `PsychoMatic'", "Paquete `PsychoMatic'", lines, fixed = TRUE)
    # DESCRIPTION is emitted as a native description list by Rd2pdf.
    labels_desc <- c(Title = "Título", Version = "Versión", Description = "Descripción",
      License = "Licencia", Encoding = "Codificación", Depends = "Depende de",
      Imports = "Importa", Suggests = "Sugiere", Author = "Autor", Maintainer = "Mantenedor",
      BugReports = "Informe de errores", Type = "Tipo", Date = "Fecha")
    for (label in names(labels_desc)) lines <- gsub(paste0("\\item[", label, "]"),
      paste0("\\item[", labels_desc[[label]], "]"), lines, fixed = TRUE)
  }
  writeLines(enc2utf8(lines), tex, useBytes = TRUE)
  old <- setwd(build)
  on.exit(setwd(old), add = TRUE)
  tools::texi2pdf("Rd2.tex", quiet = FALSE, index = TRUE,
    texinputs = file.path(R.home("share"), "texmf", "tex", "latex"))
  if (!file.copy("Rd2.pdf", pdf, overwrite = TRUE)) stop("Cannot save ", pdf)
  message("Saved ", pdf)
  invisible(pdf)
}
build_one("EN", root)
if (!english_only) build_one("ES", file.path(outdir, "es"))
