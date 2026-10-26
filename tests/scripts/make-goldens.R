# Regenerate the golden digest table used by testthat/test-10-goldens.R.
#
#   Rscript tests/scripts/make-goldens.R [library-path]
#
# With no library path the currently installed ProFound is used. To (re)create
# the reference table from a specific version, point it at a library holding
# that build, e.g. a clean checkout installed with:
#
#   git worktree add /tmp/ProFound-ref <ref>
#   R CMD INSTALL -l /tmp/pf_lib_ref /tmp/ProFound-ref
#   Rscript tests/scripts/make-goldens.R /tmp/pf_lib_ref
#
# Generating the goldens from the *older* build is what turns the suite into a
# real parity check for the speed-focused edits in src/.

# Locate the package root from the path this script was run from.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE),
                                      value = TRUE)[1])
root <- if (length(this_file) == 1L && nzchar(this_file) && file.exists(this_file)) {
  normalizePath(file.path(dirname(this_file), "..", ".."))
} else {
  normalizePath(getwd())
}
while (!file.exists(file.path(root, "DESCRIPTION")) && dirname(root) != root) {
  root <- dirname(root)
}
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("could not locate the package root; run from within the ProFound package")
}

pf_git_head <- function(dir) {
  if (!nzchar(Sys.which("git"))) {
    return("git unavailable")
  }
  out <- suppressWarnings(system2("git",
    c("-C", shQuote(dir), "rev-parse", "--short", "HEAD"),
    stdout = TRUE, stderr = TRUE))
  if (length(out) != 1L || grepl("fatal", out)) {
    return("not a git checkout")
  }
  out
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && nzchar(args[1])) {
  .libPaths(c(normalizePath(args[1], mustWork = TRUE), .libPaths()))
}
suppressPackageStartupMessages(library(ProFound))

tdir <- file.path(root, "tests", "testthat")
source(file.path(tdir, "helper-profound.R"))
source(file.path(tdir, "helper-cases.R"))

cat("ProFound used for goldens:", find.package("ProFound"), "\n")
cat("version:", as.character(packageVersion("ProFound")), "\n")

cases <- pf_cases()
stopifnot(!anyDuplicated(names(cases)))

rows <- lapply(names(cases), function(nm) {
  val <- tryCatch(cases[[nm]](), error = function(e) pf_error_val(e))
  note <- if (is.character(val) && length(val) == 1L && grepl("^ERROR:", val)) val else ""
  data.frame(case = nm, digest = pf_digest(val), note = note,
             size = as.numeric(object.size(pf_canon(val))), stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)

gdir <- file.path(tdir, "golden")
dir.create(gdir, showWarnings = FALSE, recursive = TRUE)
utils::write.csv(tab, file.path(gdir, "goldens.csv"), row.names = FALSE)

writeLines(c(
  "Golden digests for the ProFound regression case set (../helper-cases.R).",
  "",
  "Each digest is md5(serialize(value, version = 2)) of the value returned by",
  "that case after pf_canon(), which rounds numerics to 12 significant digits",
  "and relabels integer matrices into 1..k by first appearance. Watershed",
  "segment ids are an artefact of pixel processing order rather than physics,",
  "so canonicalising removes that freedom while still detecting any genuine",
  "change to segmentation geometry or photometry.",
  "",
  "Regenerate with:",
  "  Rscript tests/scripts/make-goldens.R [library-holding-the-reference-build]",
  "",
  "Provenance:",
  paste0("  built_with:  ProFound ", as.character(packageVersion("ProFound"))),
  paste0("  r_version:   ", paste0(R.version$major, ".", R.version$minor)),
  paste0("  generated:   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("  n_cases:     ", nrow(tab)),
  paste0("  source_head: ", pf_git_head(root)),
  "",
  "The checked-in reference table was produced from the committed v1.34.5",
  "source, i.e. BEFORE the speed-focused edits in src/ were committed, so the",
  "suite checks those edits against the numerical behaviour of the release that",
  "preceded them. Re-running this script against the current source is only",
  "meaningful once any intended behaviour change has been accepted.",
  "",
  "Deliberate exclusion: a `that` vector containing both negative and positive",
  "values makes this_in_that() index its lookup table out of bounds, producing",
  "heap-dependent results that differ between runs of the same build. Such",
  "input is not part of any case; see test-05-this-in-that.R.",
  "",
  "Suite validation and known pre-existing differences are documented inline in",
  "this file's current checked-in copy; regenerate it and re-add them if the",
  "reference build is ever rebased."
), file.path(gdir, "README.txt"))

cat("wrote", nrow(tab), "case digests to", file.path(gdir, "goldens.csv"), "\n")
