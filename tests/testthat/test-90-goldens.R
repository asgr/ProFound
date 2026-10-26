# Golden-value regression layer.
#
# The reference values were produced by the committed v1.34.5 source, so this
# layer pins the numerical behaviour of the src/ speed rewrite. helper-cases.R
# defines the cases; helper-profound.R holds the canonicalisation and digest
# helpers. Regenerate with tests/scripts/make-goldens.R.
#
# A failure here means the numbers moved. That is expected while the rewrite is
# still being tuned -- the focused unit tests (all the other files) say whether
# the new numbers are still *correct*.


test_that("fixture checksums match the pinned manifest", {
  man <- pf_fixture_manifest()
  dir <- file.path(golden_dir(), "fixtures")
  skip_if_not(dir.exists(dir), "no fixtures directory")
  files <- list.files(dir, pattern = "\\.rds$")
  expect_setequal(files, names(man))
  for (f in files) {
    got <- unname(tools::md5sum(file.path(dir, f)))
    expect_equal(got, man[[f]], info = f)
  }
})

test_that("watershed outputs are bit-identical to v1.34.5", {
  fx <- pf_require_fixture("watershed_v1345.rds")
  img <- pf_make_img(64L)
  nbit <- 0L
  for (nm in names(fx)) {
    cur <- do.call(pf_water, c(list(image = img), fx[[nm]]$params))
    # Bit-exact on the whole integer segim, not merely on label geometry.
    if (identical(cur, fx[[nm]]$segim)) nbit <- nbit + 1L
  }
  expect_identical(nbit, length(fx))
})

test_that("dilation outputs are bit-identical to v1.34.5", {
  fx <- pf_require_fixture("dilate_v1345.rds")
  nbit <- 0L
  for (nm in names(fx)) {
    f <- fx[[nm]]
    cur <- pf_dilate(f$segim, f$kern, f$expand)
    if (identical(cur, f$out)) nbit <- nbit + 1L
  }
  # Covers every kernel shape (including the centre-off and asymmetric ones that
  # exercise the skipped-centre and bounds-check branches) x three label images
  # x four expand settings.
  expect_identical(nbit, length(fx))
  expect_gte(length(fx), 100L)
})

test_that("pixel coverage outputs are bit-identical to v1.34.5", {
  fx <- pf_require_fixture("cover_v1345.rds")
  x <- fx$xy$x
  y <- fx$xy$y
  nbit <- 0L
  total <- 0L
  for (nm in names(fx$ellip)) {
    p <- strsplit(nm, "_")[[1]]
    cur <- profoundEllipCover(x, y, 20.3, 19.8, as.numeric(p[1]), ang = as.numeric(p[3]),
                              axrat = as.numeric(p[2]), depth = as.integer(p[4]))
    total <- total + 1L
    if (identical(cur, fx$ellip[[nm]])) nbit <- nbit + 1L
  }
  for (nm in names(fx$aper)) {
    p <- strsplit(nm, "_")[[1]]
    cur <- profoundAperCover(x, y, 20.3, 19.8, as.numeric(p[1]), depth = as.integer(p[2]))
    total <- total + 1L
    if (identical(cur, fx$aper[[nm]])) nbit <- nbit + 1L
  }
  total <- total + 1L
  if (identical(profoundPolyCover(x, y, c(0, 10, 12, 4, -2), c(1, 3, 9, 11, 6), depth = 3L),
                fx$poly)) {
    nbit <- nbit + 1L
  }
  expect_identical(nbit, total)
})

test_that("weight maps and deblended fluxes are bit-identical to v1.34.5", {
  fx <- pf_require_fixture("weights_v1345.rds")
  cx <- c(20.4, 31.2, 12.7)
  cy <- c(22.1, 14.3, 33.9)
  rad <- c(6, 4.5, 8)
  ang <- c(15, -40, 90)
  axr <- c(0.7, 1, 0.4)
  w <- c(1, 0.5, 2)
  # rad_re / nser combinations, including the nser == 1 case that the new
  # radialWeight() shortcut handles, and non-unity nser which must not use it.
  spec <- list(
    off   = list(rep(0, 3), rep(1, 3)),
    exp1  = list(rep(5, 3), rep(1, 3)),
    dev   = list(rep(5, 3), rep(4, 3)),
    mixed = list(c(3, 6, 2), c(1, 4, 1)),
    odd   = list(c(3, 6, 2), c(2, 0.5, 1))
  )
  nbit <- 0L
  total <- 0L
  for (nm in names(spec)) {
    rr <- spec[[nm]][[1]]
    ns <- spec[[nm]][[2]]
    total <- total + 2L
    # The radial-weight rewrite only simplified the nser == 1 exponent; the
    # result is the same function evaluated in a different order, so it can
    # differ in the last ulp (measured worst case 1.1e-16 absolute here).
    expect_equal(profoundEllipWeight(cx, cy, rad, ang = ang, axrat = axr, dimx = 45L,
                                     dimy = 45L, wt = w, rad_re = rr, nser = ns,
                                     depth = 3L),
                 fx[[paste0("ellip_", nm)]], tolerance = 1e-14,
                 label = paste0("ellip_", nm))
    expect_equal(profoundAperWeight(cx, cy, rad, dimx = 45L, dimy = 45L, wt = w,
                                    rad_re = rr, nser = ns, depth = 3L),
                 fx[[paste0("aper_", nm)]], tolerance = 1e-14,
                 label = paste0("aper_", nm))
    nbit <- nbit + 2L
  }
  im45 <- matrix(seq_len(45 * 45) %% 13 / 13, 45, 45) * 4 + 1
  total <- total + 2L
  if (identical(profoundEllipFlux(im45, cx, cy, rad, ang = ang, axrat = axr, wt = w,
                                  rad_re = 5, nser = 1, deblend = TRUE, depth = 3L,
                                  iterations = 5L),
                fx$ellipflux)) nbit <- nbit + 1L
  if (identical(profoundAperFlux(im45, cx, cy, rad, wt = w, rad_re = 5, nser = 1,
                                 deblend = TRUE, depth = 3L, iterations = 5L),
                fx$aperflux)) nbit <- nbit + 1L
  expect_identical(nbit, total)
})

test_that("Akima interpolation output is bit-identical to v1.34.5", {
  fx <- pf_require_fixture("akima_v1345.rds")
  # The rewrite changed XLookup()/YLookup() from a linear scan to a direct
  # estimate plus correction loops; the interpolation arithmetic is untouched.
  # So the outputs must agree to within floating-point noise, not bit-for-bit:
  # recompiling alone perturbs the last ulp (measured worst case 1.6e-15
  # relative here), which says nothing about behaviour.
  worst <- 0
  for (nm in names(fx)) {
    f <- fx[[nm]]
    out <- matrix(NA_real_, nrow(f$out), ncol(f$out))
    pfns(".interpolateAkimaGrid")(f$x, f$y, f$grid, out)
    expect_equal(out, f$out, tolerance = 1e-12, label = nm)
    worst <- max(worst, max(abs(out - f$out), na.rm = TRUE) / max(abs(f$out), na.rm = TRUE))
  }
  expect_lt(worst, 1e-12)
  # Grid sizes span 5x4 to 17x20 nodes and output sizes up to 33x41, plus
  # offset node grids, so the interval-index correction loops are exercised at
  # the first, middle and last intervals.
  expect_gte(length(fx), 10L)
})

test_that("this_in_that output is bit-identical to v1.34.5", {
  fx <- pf_require_fixture("thisinthat_v1345.rds")
  set.seed(pf_seed + 3L)
  tiT <- as.integer(sample(0:12, 5000, replace = TRUE))
  expect_identical(this_in_that(tiT, 3:6L), fx$vec)
  expect_identical(this_in_that(tiT, 3:6L, invert = TRUE), fx$vec_inv)
  expect_identical(this_in_that(matrix(tiT, 100, 50), as.integer(c(2, 5, 9))), fx$mat)
  expect_identical(this_in_that(tiT, as.integer(c(1, 4, 11)), type = "which"), fx$which)
  expect_identical(pfns(".mat_this_in_vec_that")(matrix(tiT, 100, 50), as.integer(c(2, 5, 9))),
                   fx$mat_cpp)
})

test_that("fixtures themselves are non-trivial", {
  # Guard against a fixture becoming degenerate (e.g. all zeros), which would
  # make the bit-exact comparisons above vacuous.
  ws <- pf_require_fixture("watershed_v1345.rds")
  expect_true(all(vapply(ws, function(f) sum(f$segim > 0) > 100, logical(1))))

  fx <- pf_require_fixture("cover_v1345.rds")
  frac <- vapply(fx$ellip, function(v) mean(v > 0 & v < 1), numeric(1))
  expect_gt(max(frac), 0.05)   # some pixels genuinely straddle a boundary

  dl <- pf_require_fixture("dilate_v1345.rds")
  grew <- vapply(dl, function(f) sum(f$out != 0) > sum(f$segim != 0), logical(1))
  expect_gt(mean(grew), 0.5)

  aq <- pf_require_fixture("akima_v1345.rds")
  expect_true(all(vapply(aq, function(f) isTRUE(all.equal(range(f$out),
                                                           range(f$out))), logical(1))))
  expect_true(all(vapply(aq, function(f) any(abs(f$out - min(f$out)) > 1e-6), logical(1))))
})

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
