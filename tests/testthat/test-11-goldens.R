# Golden regression layer over the full case set in helper-cases.R.
#
# The digest table in golden/goldens.csv was produced from the committed
# ProFound v1.34.5 source, so this file checks the speed-focused edits in src/
# against the numerical behaviour of the release that preceded them, across
# ~138 scenarios spanning watershed, dilation, pixel coverage, weighting,
# deblended flux, interpolation, membership lookup, the sky layer, segmentation
# statistics, aperture photometry and profoundProFound itself (including on the
# packaged real FITS data).
#
# Digests are of pf_canon(value): numerics rounded to 12 significant digits and
# integer matrices relabelled by first appearance. That tolerates last-bit
# platform noise but not any meaningful change. test-10 adds the unrounded,
# bit-exact counterpart for the routines with rewritten C++.

skip_if_not_installed("imager")   # sky/RMS normalisation uses imager (Suggests)

test_that("the golden table is present and complete", {
  tab <- pf_golden_table()
  if (is.null(tab)) {
            skip("golden/goldens.csv not found; run tests/scripts/make-goldens.R")
  }
  cases <- pf_cases()
  expect_setequal(tab$case, names(cases))
  expect_identical(anyDuplicated(tab$case), 0L)
})

test_that("every case reproduces its v1.34.5 golden digest", {
  tab <- pf_golden_table()
  if (is.null(tab)) {
    skip("golden/goldens.csv not found; run tests/scripts/make-goldens.R")
  }
  cases <- pf_cases()
  expect_setequal(names(cases), tab$case)

  goldens <- stats::setNames(tab$digest, tab$case)
  fails <- character(0)
  errs <- character(0)

  for (nm in names(cases)) {
    val <- tryCatch(cases[[nm]](), error = function(e) pf_error_val(e))
    if (length(val) == 1L && is.character(val) && grepl("^ERROR:", val)) {
      errs <- c(errs, nm)
    }
    dig <- tryCatch(pf_digest(val), error = function(e) NA_character_)
    if (!isTRUE(identical(dig, goldens[[nm]]))) {
      fails <- c(fails, nm)
    }
  }

  if (length(fails)) {
    cat("\nCases differing from the v1.34.5 golden values:\n")
    for (nm in fails) {
      cat(sprintf("  %-28s golden=%s\n", nm, goldens[[nm]]))
    }
  }
  # A single informative failure listing every differing case.
  expect_identical(fails, character(0))
})

test_that("the golden set covers the rewritten C++ routines substantially", {
  tab <- pf_golden_table()
  if (is.null(tab)) {
    skip("golden/goldens.csv not found")
  }
  grouped <- function(pattern) sum(grepl(pattern, tab$case))
  # Each of these families maps to a function whose implementation changed.
  expect_gt(grouped("^water"), 20L)     # src/water.h union-find + counting pass
  expect_gt(grouped("^dilate"), 5L)     # src/dilate.cpp scan/scatter rewrite
  expect_gt(grouped("^ellipcover|^weight|^radial|^apercover|^polycover"), 8L)
  expect_gt(grouped("^akima"), 2L)      # src/IntpAkimaUniform2.h O(1) lookup
  expect_gt(grouped("^tit"), 5L)        # this_in_that (new ProFound caller)
  expect_gt(grouped("^ppf_"), 10L)      # profoundProFound end to end
  expect_gte(sum(grepl("^ERROR:", tab$note, useBytes = TRUE)), 0L)
})
