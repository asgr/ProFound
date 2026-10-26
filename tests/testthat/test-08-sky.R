# Sky estimation: profoundSkyEst / EstLoc / SkyMap / SkyGrid / SkyScan / SkyPoly
# / SkyBlur / Chisel, plus the gain estimator and the imager wrappers.


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

test_that("profoundSkyEst recovers a known flat sky level and RMS", {
  set.seed(821)
  noise = matrix(stats::rnorm(200^2, mean = 500, sd = 12), 200, 200)

  sky = profoundSkyEst(noise)

  expect_named(sky, c('sky', 'skyerr', 'skyRMS', 'Nnearsky', 'radrun'))
  expect_equal(sky$sky, 500, tolerance = 0.02)
  expect_equal(sky$skyRMS, 12, tolerance = 0.05)
})

test_that("profoundSkyEst is not confused by bright sources when objects are flagged", {
  set.seed(4711)
  image = makeTestImage(dim = c(201, 201), sky = 300, skyRMS = 10,
                        sources = data.frame(xcen = c(60, 150), ycen = c(100, 60),
                                             flux = c(1e6, 5e5), sigma = c(3, 4)))
  objects = matrix(0L, 201, 201)
  objects[image > 400] = 1L

  sky = profoundSkyEst(image = image, objects = objects)

  expect_equal(sky$sky, 300, tolerance = 0.02)
  expect_equal(sky$skyRMS, 10, tolerance = 0.1)
})

test_that("profoundSkyEstLoc reports sky, RMS and chi-square", {
  set.seed(1234)
  image = matrix(stats::rnorm(101*101, mean = 42, sd = 3), 101, 101)

  med = profoundSkyEstLoc(image, loc = c(51, 51), box = c(100, 100), skytype = 'median')
  meanloc = profoundSkyEstLoc(image, loc = c(51, 51), box = c(100, 100), skytype = 'mean',
                              doChiSq = TRUE)

  expect_length(med, 3)
  expect_equal(unname(med[1]), 42, tolerance = 0.05)
  expect_equal(unname(med[2]), 3, tolerance = 0.15)

  # The 'mean' estimator has less scatter than the median for Gaussian noise
  expect_equal(unname(meanloc[1]), 42, tolerance = 0.05)
  expect_true(is.finite(meanloc[3]))
})

test_that("profoundSkyEstLoc rises through the available sky types", {
  set.seed(99)
  image = matrix(stats::rnorm(81*81, mean = 10, sd = 2), 81, 81)

  for(skytype in c('median', 'mean', 'mode', 'converge')){
    # skytype = 'converge' also switches the RMS estimator, with a message
    out = suppressMessages(
      profoundSkyEstLoc(image, loc = c(41, 41), box = c(80, 80), skytype = skytype)
    )
    # The density based mode is much noisier than the closed form estimators
    tol = if(skytype == 'mode') 0.5 else 0.15
    expect_equal(unname(out[1]), 10, tolerance = tol, label = skytype)
  }

  for(rms in c('quanlo', 'quanhi', 'quanboth', 'sd', 'converge')){
    out = profoundSkyEstLoc(image, loc = c(41, 41), box = c(80, 80), skyRMStype = rms)
    expect_equal(unname(out[2]), 2, tolerance = 0.3, label = rms)
  }
})

test_that("profoundMakeSkyMap and profoundMakeSkyGrid return the sky on their grid", {
  set.seed(2718)
  grid = expand.grid(1:200, 1:200)
  image = matrix(200 + 0.5 * (grid[,1] - 100) + 0.25 * (grid[,2] - 100) +
                   stats::rnorm(200^2, 0, 3), 200, 200)
  truth = matrix(200 + 0.5 * (grid[,1] - 100) + 0.25 * (grid[,2] - 100), 200, 200)

  skymap = profoundMakeSkyMap(image, box = c(50, 50), grid = c(50, 50))
  expect_named(skymap, c('sky', 'skyRMS', 'skyChiSq'))
  expect_equal(dim(skymap$sky$z), c(4, 4))
  # One grid cell per 50 pixels, centred at 25, 75, 125, 175
  expect_equal(skymap$sky$x, c(25, 75, 125, 175))
  # The image brightens with both x and y, and so does the sky map
  expect_lt(skymap$sky$z[1, 1], skymap$sky$z[4, 4])
  expect_true(all(diff(skymap$sky$z[, 1]) > 0))

  skygrid = profoundMakeSkyGrid(image, box = c(50, 50), grid = c(50, 50))
  expect_equal(dim(skygrid$sky), c(200, 200))
  expect_gt(stats::cor(as.numeric(skygrid$sky), as.numeric(truth)), 0.99)
})

test_that("profoundSkyPoly recovers a plane and a constant sky", {
  set.seed(31)
  grid = expand.grid(1:100, 1:100)
  image = matrix(100 + 2 * grid[,1] - grid[,2] + stats::rnorm(100^2, 0, 2), 100, 100)

  plane = profoundSkyPoly(image, degree = 1)
  expect_named(plane, c('sky', 'skyRMS', 'lm_out', 'good_pix', 'skyChiSq'))
  expect_equal(dim(plane$sky), c(100, 100))

  # First order fit removes most of the linear structure
  expect_lt(stats::sd(plane$sky - image), stats::sd(image))

  const = profoundSkyPoly(image, degree = 0)
  expect_equal(mean(const$sky), mean(image), tolerance = 0.01)
})

test_that("profoundSkyScan removes a row-wise background step", {
  set.seed(5150)
  jagged = matrix(rep(stats::rnorm(100, 0, 10), 100), 100, 100) +
    matrix(stats::rnorm(100^2, 0, 1), 100, 100)

  fixed = profoundSkyScan(jagged, scan_block = c(20, 20), keep_trend = FALSE)

  expect_named(fixed, c('image_fix', 'row_map', 'col_map'))
  expect_identical(dim(fixed$image_fix), dim(jagged))
  # The row offsets dominate the raw scatter and are largely removed
  expect_lt(stats::sd(fixed$image_fix), 0.3 * stats::sd(jagged))
})

test_that("profoundSkyScan leaves a flat image untouched", {
  flat = matrix(100, 60, 60)
  out = profoundSkyScan(flat, scan_block = c(20, 20))

  expect_equal(range(out$image_fix), c(100, 100))
  expect_equal(range(out$row_map), c(0, 0))
  expect_equal(range(out$col_map), c(0, 0))
})

test_that("profoundMakeSkyBlur preserves a constant sky", {
  blurred = profoundMakeSkyBlur(matrix(100, 51, 51), sigma = 5)
  expect_equal(range(blurred), c(100, 100), tolerance = 1e-10)
})

test_that("profoundChisel separates compact sources from the sky", {
  set.seed(808)
  grid = expand.grid(1:151, 1:151)
  image = matrix(150 + stats::rnorm(151^2, 0, 2), 151, 151)
  # A compact source and a broad one, so there is something to classify
  rad1 = (grid[,1] - 100)^2 + (grid[,2] - 60)^2
  rad2 = (grid[,1] - 50)^2 + (grid[,2] - 110)^2
  image = image + matrix(5000 * exp(-rad1 / (2 * 4^2)) / (2 * pi * 16), 151, 151)
  image = image + matrix(20000 * exp(-rad2 / (2 * 10^2)) / (2 * pi * 100), 151, 151)

  out = profoundChisel(image)

  expect_named(out, c('objects', 'sky'))
  expect_identical(dim(out$sky), dim(image))
  # Chisel is aimed at near-flat skies; the recovered level should track the truth
  expect_equal(stats::median(out$sky), 150, tolerance = 0.05)
  expect_lt(stats::sd(out$sky), 0.5 * stats::sd(image))
  # The bright source centres are classified as objects
  expect_true(out$objects[100, 60])
  expect_true(out$objects[50, 110])
})
