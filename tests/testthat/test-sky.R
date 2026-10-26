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
