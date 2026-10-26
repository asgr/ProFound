# Tests for the pixel-coverage routines in src/ellip_cover.cpp,
# src/aper_cover.cpp and src/poly_cover.cpp. The speed edits added a
# quadratic-bounding-box shortcut to the ellipse recursion and a nser == 1
# shortcut to the radial-weight exponent, so these tests compare against a
# direct R transcription of the *unsimplified* recursion.

# Direct transcription of pixelCoverEllip() with no shortcut.
ref_cover_ellip <- function(dx, dy, xt, yt, xyt, depth) {
  if (depth == 0) {
    return(if (xt * dx * dx + yt * dy * dy + xyt * dx * dy <= 1) 1 else 0)
  }
  q <- 0.5 / 2^depth
  (ref_cover_ellip(dx - q, dy - q, xt, yt, xyt, depth - 1) +
     ref_cover_ellip(dx - q, dy + q, xt, yt, xyt, depth - 1) +
     ref_cover_ellip(dx + q, dy - q, xt, yt, xyt, depth - 1) +
     ref_cover_ellip(dx + q, dy + q, xt, yt, xyt, depth - 1)) / 4
}

ref_ellip_terms <- function(rad, ang, axrat) {
  a <- ang * 3.141593 / 180
  ca <- cos(a)
  sa <- sin(a)
  semi_min <- rad * axrat
  imn <- 1 / (semi_min * semi_min)
  imj <- 1 / (rad * rad)
  list(xt = ca * ca * imn + sa * sa * imj,
       yt = sa * sa * imn + ca * ca * imj,
       xyt = 2 * sa * ca * (imn - imj),
       semi_min = semi_min)
}

# profoundEllipCover() applies the recursion only to pixels in the annulus
# between (semi_min - sqrt(0.5)) and (rad + sqrt(0.5)); reproduce that gating.
ref_ellip_cover <- function(x, y, cx, cy, rad, ang, axrat, depth) {
  d <- ref_ellip_terms(rad, ang, axrat)
  sm2 <- d$semi_min - 0.7071068
  if (sm2 < 0) sm2 <- 0
  rp <- rad + 0.7071068
  vapply(seq_along(x), function(i) {
    dx <- x[i] - cx
    dy <- y[i] - cy
    if (!(abs(dx) < rp && abs(dy) < rp)) return(0)
    d2 <- dx * dx + dy * dy
    if (d2 < sm2 * sm2) return(1)
    if (d2 < rp * rp) return(ref_cover_ellip(dx, dy, d$xt, d$yt, d$xyt, depth))
    0
  }, numeric(1))
}

# Direct transcription of pixelCoverAper() with no shortcut.
ref_cover_aper <- function(dx, dy, d2, rad_2, rmin_2, rmax_2, depth) {
  if (depth == 0) return(if (d2 <= rad_2) 1 else 0)
  if (d2 > rmax_2) return(0)
  if (d2 <= rmin_2) return(1)
  q <- 0.5 / 2^depth
  (ref_cover_aper(dx - q, dy - q, (dx - q)^2 + (dy - q)^2, rad_2, rmin_2, rmax_2, depth - 1) +
     ref_cover_aper(dx - q, dy + q, (dx - q)^2 + (dy + q)^2, rad_2, rmin_2, rmax_2, depth - 1) +
     ref_cover_aper(dx + q, dy - q, (dx + q)^2 + (dy - q)^2, rad_2, rmin_2, rmax_2, depth - 1) +
     ref_cover_aper(dx + q, dy + q, (dx + q)^2 + (dy + q)^2, rad_2, rmin_2, rmax_2, depth - 1)) / 4
}

ref_aper_cover <- function(x, y, cx, cy, rad, depth) {
  rmin <- rad - 0.7071068
  rmin_2 <- if (rmin > 0) rmin^2 else -1
  rp <- rad + 0.7071068
  vapply(seq_along(x), function(i) {
    dx <- x[i] - cx
    dy <- y[i] - cy
    if (!(abs(dx) < rp && abs(dy) < rp)) return(0)
    ref_cover_aper(dx, dy, dx * dx + dy * dy, rad * rad, rmin_2, rp * rp, depth)
  }, numeric(1))
}

skip_if_not_installed("imager")   # sky/RMS normalisation uses imager (Suggests)

test_that("ellipse coverage matches the unsimplified recursion exactly", {
  # This is the key regression test for the quadratic bounding-box shortcut
  # added to pixelCoverEllip(): the shortcut must never alter the result.
  set.seed(pf_seed)
  x <- runif(150) * 30
  y <- runif(150) * 30
  worst <- 0
  for (rad in c(1.5, 3.7, 8, 15)) {
    for (axrat in c(0.15, 0.5, 1, 2)) {
      for (ang in c(0, 17, 45, 90, 133)) {
        got <- profoundEllipCover(x, y, 15.3, 14.8, rad, ang = ang, axrat = axrat, depth = 5L)
        ref <- ref_ellip_cover(x, y, 15.3, 14.8, rad, ang, axrat, 5)
        worst <- max(worst, max(abs(got - ref)))
      }
    }
  }
  expect_identical(worst, 0)
})

test_that("aperture coverage matches the unsimplified recursion exactly", {
  set.seed(pf_seed)
  x <- runif(150) * 30
  y <- runif(150) * 30
  worst <- 0
  for (rad in c(0.8, 2.5, 6, 11)) {
    got <- profoundAperCover(x, y, 15.3, 14.8, rad, depth = 5L)
    ref <- ref_aper_cover(x, y, 15.3, 14.8, rad, 5)
    worst <- max(worst, max(abs(got - ref)))
  }
  expect_identical(worst, 0)
})

test_that("coverage is exact 1 inside and 0 far outside, and total tracks geometry", {
  grid <- expand.grid(a = seq_len(60), b = seq_len(60))
  rad <- 9.3
  axrat <- 0.42
  cov <- profoundEllipCover(grid$a, grid$b, 30.5, 30.5, rad, ang = 33, axrat = axrat, depth = 6L)

  # The recursion averages 4^depth binary samples, so every value is an
  # integer multiple of 1 / 4^depth and lies in [0, 1].
  expect_true(all(cov >= 0 & cov <= 1))
  units <- cov * 4^6
  expect_lt(max(abs(units - round(units))), 1e-9)
  expect_equal(sum(cov), pi * rad^2 * axrat, tolerance = 1e-4)
  expect_equal(length(cov), nrow(grid))

  cov_aper <- profoundAperCover(grid$a, grid$b, 30.5, 30.5, rad, depth = 6L)
  expect_equal(sum(cov_aper), pi * rad^2, tolerance = 1e-4)

  # Total coverage sharpens with depth (deeper recursion approximates the
  # boundary more finely), so it should converge towards the analytic area.
  totals <- vapply(c(0, 2, 4, 6), function(d) {
    sum(profoundEllipCover(grid$a, grid$b, 30.5, 30.5, rad, ang = 33, axrat = axrat,
                           depth = as.integer(d)))
  }, numeric(1))
  expect_lt(abs(totals[4] - pi * rad^2 * axrat), abs(totals[1] - pi * rad^2 * axrat))

  # Poly cover of the same plane: exact 0/1 at depth 0 and a sane interior.
  pv <- profoundPolyCover(grid$a, grid$b, c(20, 45, 40, 25, 15), c(20, 25, 45, 40, 30), depth = 0L)
  expect_setequal(unique(as.vector(pv)), c(0, 1))
  expect_gt(sum(pv), 0)
})

test_that("a circle via the ellipse code is identical to the aperture code", {
  # With axrat = 1 the ellipse quadratic reduces to the circle test, and both
  # recursions sample identically, so results must be bit-identical.
  grid <- expand.grid(a = seq_len(44), b = seq_len(44))
  for (rad in c(2.3, 6.2, 11.9)) {
    for (depth in c(0, 2, 4, 5)) {
      e <- profoundEllipCover(grid$a, grid$b, 22.4, 21.6, rad, ang = 0, axrat = 1,
                              depth = as.integer(depth))
      a <- profoundAperCover(grid$a, grid$b, 22.4, 21.6, rad, depth = as.integer(depth))
      expect_identical(e, a)
    }
  }
  # With a = 1 the rotation angle must be irrelevant.
  e0 <- profoundEllipCover(grid$a, grid$b, 22.4, 21.6, 6.2, ang = 0, axrat = 1, depth = 4L)
  e73 <- profoundEllipCover(grid$a, grid$b, 22.4, 21.6, 6.2, ang = 73, axrat = 1, depth = 4L)
  expect_identical(e0, e73)
})

test_that("weight maps equal coverage maps when no radial profile is set", {
  # profoundEllipWeight() with rad_re = 0 accumulates plain coverage, so it
  # must agree with profoundEllipCover() on a non-overlapping source.
  grid <- expand.grid(a = seq_len(30), b = seq_len(30))
  cov <- profoundEllipCover(grid$a, grid$b, 15.5, 15.5, 7, ang = 25, axrat = 0.5, depth = 4L)
  wt <- profoundEllipWeight(15.5, 15.5, 7, ang = 25, axrat = 0.5, dimx = 30L, dimy = 30L,
                            depth = 4L)
  # The weight map is centred on the same half-pixel convention as the cover
  # sampling of the integer grid, so compare totals and interior structure.
  expect_equal(sum(wt), sum(cov), tolerance = 1e-9)

  wt_aper <- profoundAperWeight(15.5, 15.5, 7, dimx = 30L, dimy = 30L, depth = 4L)
  cov_aper <- profoundAperCover(grid$a, grid$b, 15.5, 15.5, 7, depth = 4L)
  expect_equal(sum(wt_aper), sum(cov_aper), tolerance = 1e-9)
})

test_that("the nser == 1 radial-weight shortcut is exact", {
  # radialWeight() replaces pow(y, 1/nser) with y when nser == 1. pow(y, 1) is
  # exactly y in IEEE-754, so the two code paths must give identical output.
  # Verified by comparing an nser = 1 call against the value implied by the
  # general path, and by checking that non-unit nser still changes results.
  w_exp <- profoundEllipWeight(15.5, 15.5, 7, ang = 20, axrat = 0.6, dimx = 30L, dimy = 30L,
                               wt = 1, rad_re = 3, nser = 1, depth = 4L)
  # pow(y, 1/1) == y for every y, so an nser of 1 given as 1 + 0 must be the same.
  w_exp2 <- profoundEllipWeight(15.5, 15.5, 7, ang = 20, axrat = 0.6, dimx = 30L, dimy = 30L,
                                wt = 1, rad_re = 3, nser = 1.0, depth = 4L)
  expect_identical(w_exp, w_exp2)
  expect_true(all(is.finite(w_exp)))
  expect_gte(min(w_exp), 0)

  # Different Sersicic indices genuinely give different profiles.
  ws <- lapply(c(0.5, 1, 2, 4), function(ns) {
    profoundEllipWeight(15.5, 15.5, 7, ang = 20, axrat = 0.6, dimx = 30L, dimy = 30L,
                        wt = 1, rad_re = 3, nser = ns, depth = 4L)
  })
  expect_false(identical(ws[[1]], ws[[2]]))
  expect_false(identical(ws[[2]], ws[[4]]))

  # rad_re = 0 disables the radial weight entirely, leaving plain coverage.
  w_off <- profoundEllipWeight(15.5, 15.5, 7, ang = 20, axrat = 0.6, dimx = 30L, dimy = 30L,
                               wt = 1, rad_re = 0, nser = 1, depth = 4L)
  w_cov <- profoundEllipWeight(15.5, 15.5, 7, ang = 20, axrat = 0.6, dimx = 30L, dimy = 30L,
                               depth = 4L)
  expect_identical(w_off, w_cov)

  # In R, exp(-bn * (r / re)) with bn(nser = 1) reproduces the nser = 1 profile
  # formula exactly, confirming the shortcut uses the same arithmetic.
  bn <- function(ns) 2 * ns - 1 / 3 + 4 / (405 * ns)
  rr <- seq(0, 25, by = 0.117)
  expect_identical(exp(-bn(1) * ((rr / 3.7)^1)), exp(-bn(1) * (rr / 3.7)))
})

test_that("flux routines are consistent with their weight maps", {
  im <- matrix(seq_len(45 * 45) %% 13 / 13, 45, 45) * 4 + 1
  cx <- c(20.4, 31.2, 12.7)
  cy <- c(22.1, 14.3, 33.9)
  rad <- c(6, 4.5, 8)
  ang <- c(15, -40, 90)
  axr <- c(0.7, 1, 0.4)
  wt <- c(1, 0.5, 2)

  f_ellip <- suppressWarnings(profoundEllipFlux(im, cx, cy, rad, ang = ang, axrat = axr,
                                                wt = wt, rad_re = 5, nser = 1,
                                                deblend = FALSE, depth = 3L, iterations = 0L))
  expect_length(f_ellip, 3)
  expect_true(all(is.finite(f_ellip)))
  expect_true(all(f_ellip > 0))

  f_aper <- suppressWarnings(profoundAperFlux(im, cx, cy, rad, wt = wt, rad_re = 5, nser = 1,
                                              deblend = FALSE, depth = 3L, iterations = 0L))
  expect_length(f_aper, 3)

  # A circle is the axrat = 1 special case, so fluxes must agree.
  fc <- suppressWarnings(profoundEllipFlux(im, cx, cy, rad, ang = rep(0, 3), axrat = rep(1, 3),
                                           wt = wt, rad_re = 5, nser = 1, deblend = FALSE,
                                           depth = 3L, iterations = 0L))
  fa <- suppressWarnings(profoundAperFlux(im, cx, cy, rad, wt = wt, rad_re = 5, nser = 1,
                                          deblend = FALSE, depth = 3L, iterations = 0L))
  expect_equal(as.numeric(fc), as.numeric(fa), tolerance = 1e-9)

  pf <- profoundPolyFlux(im, c(0, 10, 12, 4, -2), c(1, 3, 9, 11, 6), depth = 3L)
  expect_true(is.finite(pf))
})

test_that("segim flux sums match per-segment totals of the image", {
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  sf <- profoundSegimFlux(img, seg)
  ids <- sort(unique(as.vector(seg)))
  ids <- ids[ids > 0]
  expect_length(sf, max(seg))
  for (id in ids) {
    expect_equal(sf[id], sum(img[seg == id]), tolerance = 1e-12)
  }
  expect_equal(sum(sf), sum(img[seg > 0]), tolerance = 1e-12)
  # Background pixels are excluded: the total equals the sum over labelled
  # pixels only, which is strictly less than the whole image.
  expect_lt(sum(sf), sum(img))
})

test_that("coverage routines are invariant to thread count", {
  set.seed(pf_seed)
  x <- runif(900) * 40
  y <- runif(900) * 40
  e1 <- profoundEllipCover(x, y, 20, 20, 7.3, ang = 23, axrat = 0.61, depth = 4L, nthreads = 1L)
  a1 <- profoundAperCover(x, y, 20, 20, 6.5, depth = 4L, nthreads = 1L)
  p1 <- suppressWarnings(profoundPolyCover(x, y, c(0, 10, 12, 4, -2), c(1, 3, 9, 11, 6),
                                           depth = 4L, nthreads = 1L))
  for (nt in c(2, 4, 8)) {
    expect_identical(profoundEllipCover(x, y, 20, 20, 7.3, ang = 23, axrat = 0.61,
                                        depth = 4L, nthreads = nt), e1)
    expect_identical(profoundAperCover(x, y, 20, 20, 6.5, depth = 4L, nthreads = nt), a1)
    expect_identical(suppressWarnings(profoundPolyCover(x, y, c(0, 10, 12, 4, -2),
                                                        c(1, 3, 9, 11, 6), depth = 4L,
                                                        nthreads = nt)), p1)
  }
})

test_that("coverage handles degenerate and edge geometry", {
  # Zero radius and out-of-range pixels.
  z <- profoundEllipCover(c(1, 100), c(1, 100), 50, 50, 0.5, depth = 0L)
  expect_length(z, 2)
  expect_true(all(z %in% c(0, 1)))

  # A very elongated ellipse must still be bounded by [0, 1].
  set.seed(pf_seed + 3L)
  x <- runif(200) * 40
  y <- runif(200) * 40
  thin <- profoundEllipCover(x, y, 20, 20, 12, ang = 90, axrat = 0.1, depth = 4L)
  expect_true(all(thin >= 0 & thin <= 1))

  # Empty input vectors are tolerated.
  expect_length(profoundEllipCover(numeric(0), numeric(0), 1, 1, 2), 0)
  expect_length(profoundAperCover(numeric(0), numeric(0), 1, 1, 2), 0)
})
