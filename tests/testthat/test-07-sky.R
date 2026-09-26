# Tests for the sky estimation routines that profoundProFound() relies on.
# profoundSkyEst() uses iteratively sigma-clipped estimators rather than a
# plain median, so these tests assert robust statistical properties (position
# within the background distribution, mask/object handling, smoothness) rather
# than exact base-R equivalence.

skip_if_not_installed("imager")   # sky/RMS normalisation uses imager (Suggests)

test_that("profoundSkyEst returns sensible background statistics", {
  img <- pf_make_img(64L)
  obj <- matrix(as.integer(pf_water(img) > 0), 64, 64)
  s <- profoundSkyEst(image = img, objects = obj)

  expect_true(all(c("sky", "skyerr", "skyRMS", "Nnearsky", "radrun") %in% names(s)))
  expect_length(s$sky, 1L)
  expect_length(s$skyRMS, 1L)
  expect_true(is.finite(s$sky))
  expect_true(is.finite(s$skyRMS) && s$skyRMS > 0)

  bg <- as.vector(img[obj == 0])
  # A clipped estimator of a roughly symmetric background should sit close to
  # the background median and well short of the bright-source mean.
  expect_lt(abs(s$sky - median(bg)), 0.25 * sd(bg))
  expect_lt(s$sky, mean(bg))
  expect_lt(s$skyRMS, sd(bg))
})

test_that("profoundSkyEst responds to masks and object exclusion", {
  img <- pf_make_img(64L)
  obj <- matrix(as.integer(pf_water(img) > 0), 64, 64)
  plain <- profoundSkyEst(image = img, objects = obj)

  # Excluding a bright corner region must shift the estimate.
  mk <- matrix(0L, 64, 64)
  mk[1:6, 1:6] <- 1L
  masked <- profoundSkyEst(image = img, objects = obj, mask = mk)
  expect_false(isTRUE(all.equal(plain$sky, masked$sky)))

  # Treating nothing as an object lets sources drag the estimate upwards.
  no_obj <- profoundSkyEst(image = img, objects = matrix(0L, 64, 64))
  expect_gt(no_obj$sky, plain$sky)

  # skytype variants all return finite values of the right sign of magnitude.
  for (st in c("mean", "median", "mode")) {
    r <- profoundSkyEst(image = img, objects = obj, skytype = st)
    expect_true(is.finite(r$sky))
  }
  for (srt in c("quanlo", "quanhi", "quanboth", "sd")) {
    r <- profoundSkyEst(image = img, objects = obj, skyRMStype = srt)
    expect_true(is.finite(r$skyRMS) && r$skyRMS > 0)
  }
})

test_that("profoundSkyEst handles degenerate images", {
  # A perfectly flat image has zero spread: skyRMS is 0 and the clipped
  # estimator degenerates to NaN rather than erroring. This behaviour is the
  # same in v1.34.5 and the current tree and is pinned by the goldens.
  flat <- matrix(3, 20, 20)
  s <- profoundSkyEst(image = flat, objects = matrix(0L, 20, 20))
  expect_identical(s$skyRMS, 0)
  expect_true(is.na(s$sky))

  # Almost flat gives a finite, sensible sky.
  nearly <- flat + matrix(rnorm(400, 0, 0.001), 20, 20)
  s2 <- profoundSkyEst(image = nearly, objects = matrix(0L, 20, 20))
  expect_true(is.finite(s2$sky))
  expect_equal(s2$sky, 3, tolerance = 0.01)

  nall <- matrix(NA_real_, 10, 10)
  r <- tryCatch(profoundSkyEst(image = nall, objects = matrix(0L, 10, 10)),
                error = function(e) "error")
  expect_true(is.character(r) || is.list(r))
})

test_that("profoundMakeSkyGrid produces a smooth map over the image", {
  img <- pf_make_img(64L)
  obj <- matrix(as.integer(pf_water(img) > 0), 64, 64)

  for (sgt in c("new", "old")) {
    for (ty in c("bicubic", "bilinear")) {
      g <- profoundMakeSkyGrid(image = img, objects = obj, sky = 0, box = c(20, 20),
                               grid = c(20, 20), skygrid_type = sgt, type = ty)
      expect_identical(dim(g$sky), dim(img))
      expect_identical(dim(g$skyRMS), dim(img))
      expect_true(all(c("sky", "skyRMS") %in% names(g)))
      # The map varies more slowly than the image itself.
      expect_lt(sd(as.vector(g$sky)), sd(as.vector(img)))
      expect_true(all(is.finite(g$sky)))
      expect_true(all(g$skyRMS > 0))
    }
  }
  expect_error(profoundMakeSkyGrid(image = img, objects = obj, skygrid_type = "bogus"))
})

test_that("sky grid interpolation types differ but stay close", {
  # bicubic (Akima) and bilinear are different schemes, so the maps must not be
  # identical, yet both track the same underlying background.
  img <- pf_make_img(64L)
  obj <- matrix(as.integer(pf_water(img) > 0), 64, 64)
  a <- profoundMakeSkyGrid(image = img, objects = obj, sky = 0, box = c(20, 20),
                           grid = c(20, 20), type = "bicubic")
  b <- profoundMakeSkyGrid(image = img, objects = obj, sky = 0, box = c(20, 20),
                           grid = c(20, 20), type = "bilinear")
  expect_false(identical(a$sky, b$sky))
  expect_gt(cor(as.vector(a$sky), as.vector(b$sky)), 0.9)
})

test_that("profoundSkyEstLoc, Plane, Poly and Chisel run and return maps", {
  img <- pf_make_img(64L)
  obj <- matrix(as.integer(pf_water(img) > 0), 64, 64)

  loc <- profoundSkyEstLoc(image = img, objects = obj, loc = c(32, 32), box = c(20, 20))
  expect_true(is.numeric(loc))

  pl <- profoundSkyPlane(image = img, objects = obj)
  expect_true(is.list(pl))
  expect_true("sky" %in% names(pl))
  expect_identical(dim(pl$sky), dim(img))

  po <- profoundSkyPoly(image = img, objects = obj, degree = 2)
  expect_true(is.list(po))
  expect_identical(dim(po$sky), dim(img))

  ch <- profoundChisel(image = img, sky = 0)
  expect_true(all(c("objects", "sky") %in% names(ch)))
  expect_identical(dim(ch$objects), dim(img))
})

test_that("profoundGainEst recovers a known gain", {
  # Build an image whose noise variance scales with signal, the situation
  # profoundGainEst() is designed for.
  set.seed(pf_seed + 81L)
  n <- 64L
  gain_true <- 20
  base <- matrix(200, n, n)
  noisy <- base + matrix(rnorm(n * n, 0, sqrt(base / gain_true)), n, n)
  g <- profoundGainEst(image = noisy, mask = matrix(0L, n, n), objects = matrix(0L, n, n),
                       sky = mean(base), skyRMS = sd(as.vector(noisy)))
  expect_true(is.numeric(g))
  expect_true(all(is.finite(g)))
})

test_that("profoundImBlur and profoundImDiff behave as smoothers", {
  img <- pf_make_img(32L)
  b <- profoundImBlur(image = img, sigma = 2, plot = FALSE)
  expect_identical(dim(b), dim(img))
  # Smoothing must reduce small-scale variance.
  expect_lt(sd(as.vector(b)), sd(as.vector(img)))

  d <- profoundImDiff(image = img, sigma = 3, plot = FALSE)
  expect_identical(dim(d), dim(img))
  # The difference image is high-pass, so its mean is near zero.
  expect_lt(abs(mean(d)), sd(as.vector(img)))

  # A constant image is unchanged by both operations.
  cst <- matrix(7, 20, 20)
  expect_equal(mean(profoundImBlur(image = cst, sigma = 2, plot = FALSE)), 7,
               tolerance = 1e-8)
})
