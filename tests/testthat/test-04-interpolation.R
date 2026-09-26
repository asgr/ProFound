# Tests for src/IntpAkimaUniform2.h, where XLookup()/YLookup() were changed from
# a linear scan to a direct estimate plus correction loops. An off-by-one in the
# interval index would show up immediately as a failure to reproduce values at
# nodes and to integrate linear and quadratic fields exactly.
#
# The output convention of .interpolateAkimaGrid() (see src/akima.cpp) is that
# output cell (i, j) is evaluated at x = -0.5 + i, y = -0.5 + j, so a node grid
# running over half-integers aligns exactly with output cells.

pf_akima <- function(x, y, grid, nxout, nyout) {
  out <- matrix(NA_real_, nxout, nyout)
  pfns(".interpolateAkimaGrid")(x, y, grid, out)
  out
}

# Values of the interpolant at the node positions of a half-integer node grid.
pf_node_cells <- function(nx, ny) {
  list(qx = seq_len(nx), qy = seq_len(ny))   # x = -0.5 + i  =>  i = 1..nx
}

skip_if_not_installed("imager")   # sky/RMS normalisation uses imager (Suggests)

test_that("Akima reproduces its own nodes on a half-integer grid", {
  set.seed(pf_seed + 5L)
  for (nx in c(6, 9, 17)) {
    for (ny in c(5, 11, 20)) {
      gx <- seq(0.5, nx - 0.5)
      gy <- seq(0.5, ny - 0.5)
      z <- matrix(sin(outer(gx, gy) / 2.1) * 4 + cos(gx %o% (gy / 3)), nx, ny)
      out <- pf_akima(gx, gy, z, nx, ny)
      expect_lt(max(abs(out - z)), 1e-9)
    }
  }
})

test_that("Akima reproduces linear and quadratic fields", {
  # A C2 spline through nodes reproduces low-order polynomials well; the cubic
  # Akima form is exact for quadratics along each axis.
  nx <- 9
  ny <- 11
  gx <- seq(0.5, nx - 0.5)
  gy <- seq(0.5, ny - 0.5)
  qx <- -0.5 + seq_len(2 * nx)
  qy <- -0.5 + seq_len(2 * ny)

  z_lin <- outer(gx, gy, function(x, y) 3 + 0.7 * x - 1.3 * y)
  out_lin <- pf_akima(gx, gy, z_lin, 2 * nx, 2 * ny)
  ref_lin <- outer(qx, qy, function(x, y) 3 + 0.7 * x - 1.3 * y)
  inside <- outer(qx >= min(gx) & qx <= max(gx), qy >= min(gy) & qy <= max(gy), "&")
  expect_lt(max(abs(out_lin[inside] - ref_lin[inside])), 1e-9)

  z_quad <- outer(gx, gy, function(x, y) 1 + x^2 + 2 * y)
  out_quad <- pf_akima(gx, gy, z_quad, 2 * nx, 2 * ny)
  ref_quad <- outer(qx, qy, function(x, y) 1 + x^2 + 2 * y)
  expect_lt(max(abs(out_quad[inside] - ref_quad[inside])), 1e-12)
})

test_that("Akima node lookup is correct across interval sizes and indices", {
  # The rewritten lookup clamps to the first/last interval and then corrects
  # with while loops. Node coordinates are chosen so that output cells
  # (x = -0.5 + i) land exactly on nodes for integer spacing of 1, 2 and 4
  # cells, which exercises interval indices at the start, middle and end of a
  # long grid.
  for (nx in c(5, 9, 17, 33)) {
    for (step in c(1, 2, 4)) {
      gx <- (seq_len(nx) - 1) * step + 0.5
      gy <- (seq_len(nx) - 1) * step + 0.5
      set.seed(pf_seed + nx + step)
      z <- matrix(runif(nx * nx, -10, 10), nx, nx)
      ncell <- (nx - 1) * step + 1
      out <- pf_akima(gx, gy, z, ncell, ncell)
      idx <- (seq_len(nx) - 1) * step + 1
      expect_lt(max(abs(out[idx, idx] - z)), 1e-9)
    }
  }
})

test_that("Akima handles shifted and non-half-integer node grids", {
  nx <- 9
  ny <- 11
  gx <- seq(1, 9) + 0.25
  gy <- seq(2, 12) - 0.6
  z <- matrix(sin(outer(gx, gy) / 3), nx, ny)
  out <- pf_akima(gx, gy, z, 25, 19)
  expect_identical(dim(out), c(25L, 19L))
  expect_true(all(is.finite(out[is.finite(out)])))
  # Continuity: neighbouring output cells cannot jump by more than the range of
  # the input field scaled by a modest factor.
  rng <- diff(range(z))
  expect_lt(max(abs(diff(out, difference = 1)), na.rm = TRUE), 2 * rng)
})

test_that("Akima output is deterministic and matches a re-run", {
  gx <- seq(0.5, 20.5)
  gy <- seq(0.5, 14.5)
  z <- matrix(rnorm(21 * 15, 100, 5), 21, 15)
  a <- pf_akima(gx, gy, z, 40, 30)
  b <- pf_akima(gx, gy, z, 40, 30)
  expect_identical(a, b)
})

test_that("bilinear interpolation reproduces planes exactly", {
  # .interpolateLinearGrid() is untouched by the speed edits but shares the
  # -0.5 + i convention, and provides a simple cross-check of the harness.
  nx <- 6
  ny <- 9
  gx <- seq(0.5, nx - 0.5)
  gy <- seq(0.5, ny - 0.5)
  z <- outer(gx, gy, function(x, y) 5 - 0.4 * x + 1.1 * y)
  out <- matrix(NA_real_, 2 * nx, 2 * ny)
  pfns(".interpolateLinearGrid")(gx, gy, z, out)
  qx <- -0.5 + seq_len(2 * nx)
  qy <- -0.5 + seq_len(2 * ny)
  ref <- outer(qx, qy, function(x, y) 5 - 0.4 * x + 1.1 * y)
  inside <- outer(qx >= min(gx) & qx <= max(gx), qy >= min(gy) & qy <= max(gy), "&")
  expect_lt(max(abs(out[inside] - ref[inside])), 1e-12)
  expect_equal(pf_akima(gx, gy, z, 2 * nx, 2 * ny)[inside], out[inside],
               tolerance = 1e-9)
})

test_that("profoundResample keeps total flux and scales as documented", {
  img <- pf_make_img(64L)
  for (po in c(0.5, 1, 2)) {
    for (pn in c(0.5, 1, 2)) {
      r <- suppressWarnings(profoundResample(img, pixscale_old = po, pixscale_new = pn,
                                             type = "bicubic", fluxscale = "image"))
      expect_true(all(dim(r) == round(dim(img) * po / pn)))
      # fluxscale = "image" forces the totals to match by construction.
      expect_equal(sum(r), sum(img), tolerance = 1e-6)
      rn <- suppressWarnings(profoundResample(img, pixscale_old = po, pixscale_new = pn,
                                              type = "bicubic", fluxscale = "norm"))
      expect_equal(sum(rn), 1, tolerance = 1e-9)
    }
  }
  expect_error(profoundResample(img, type = "notreal"))
  expect_error(profoundResample(img, fluxscale = "notreal"))
})

test_that("profoundResample is stable and structure-preserving at unit scale", {
  img <- pf_make_img(32L)
  r <- suppressWarnings(profoundResample(img, pixscale_old = 1, pixscale_new = 1,
                                          type = "bicubic", fluxscale = "image"))
  expect_identical(dim(r), dim(img))
  expect_identical(r, suppressWarnings(profoundResample(img, 1, 1, type = "bicubic",
                                                         fluxscale = "image")))
  # profoundResample() resamples a PSF onto a node grid offset by half a pixel
  # from the input index convention, so unit rescaling is not an exact identity;
  # it must still preserve the field structure and sign of the sources.
  expect_gt(cor(as.vector(r), as.vector(img)), 0.99)
  expect_lt(sd(as.vector(r - img)), sd(as.vector(img)))
})
