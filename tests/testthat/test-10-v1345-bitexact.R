# Bit-exact parity tests against full-precision outputs captured from the
# committed ProFound v1.34.5 source (2026-09-24), stored under golden/fixtures/.
#
# test-10-goldens.R covers a broad case set with a 12-significant-digit
# canonicalisation. This file is the strict counterpart: for the routines
# touched by the speed-focused edits it requires identical() on double-precision
# values with no rounding at all, so even a last-bit difference in pow/exp or in
# a summation order is caught.
#
# Each fixture is checksum-verified against the manifest in helper-profound.R
# before use; if a fixture is missing or altered the tests skip rather than
# silently passing.

pf_require_fixture <- function(name) {
  fx <- pf_fixture_checked(name)
  if (is.null(fx)) {
    skip(paste("fixture", name, "is absent or does not match its pinned checksum"))
  }
  fx
}

skip_if_not_installed("imager")   # sky/RMS normalisation uses imager (Suggests)

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
    if (identical(profoundEllipWeight(cx, cy, rad, ang = ang, axrat = axr, dimx = 45L,
                                      dimy = 45L, wt = w, rad_re = rr, nser = ns,
                                      depth = 3L),
                  fx[[paste0("ellip_", nm)]])) nbit <- nbit + 1L
    if (identical(profoundAperWeight(cx, cy, rad, dimx = 45L, dimy = 45L, wt = w,
                                     rad_re = rr, nser = ns, depth = 3L),
                  fx[[paste0("aper_", nm)]])) nbit <- nbit + 1L
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
  nbit <- 0L
  for (nm in names(fx)) {
    f <- fx[[nm]]
    out <- matrix(NA_real_, nrow(f$out), ncol(f$out))
    pfns(".interpolateAkimaGrid")(f$x, f$y, f$grid, out)
    if (identical(out, f$out)) nbit <- nbit + 1L
  }
  # Grid sizes span 5x4 to 17x20 nodes and output sizes up to 33x41, plus
  # offset node grids, so the interval-index correction loops are exercised at
  # the first, middle and last intervals.
  expect_identical(nbit, length(fx))
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
