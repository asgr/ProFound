# Image manipulation helpers: smoothing, differences, gradients, sigma maps and
# resampling. (The imager wrappers are also exercised from test-sky.R, where
# they feed the sky estimator.)


test_that("profoundMakeSigma combines noise sources in quadrature", {
  image = matrix(100, 11, 11)

  sigma = profoundMakeSigma(image = image, sky = 100, skyRMS = 5, gain = 1)

  # Pure sky: sigma is just skyRMS away from the sources
  expect_true(all(sigma == 5))

  # read noise adds in quadrature (0 counts above sky)
  sigma_read = profoundMakeSigma(image = image, sky = 100, skyRMS = 3, readRMS = 4)
  expect_true(all(abs(sigma_read - 5) < 1e-10))
})

test_that("profoundMakeSigma rejects units it does not understand", {
  image = matrix(100, 11, 11)

  expect_error(profoundMakeSigma(image = image, image_units = 'junk'), "not recognised")
  expect_error(profoundMakeSigma(image = image, sky_units = 'junk'), "not recognised")
  expect_error(profoundMakeSigma(image = image, image_units = 'ADU/s'),
               "exptime cannot be NULL")
})

test_that("profoundMakeSigma converts between electrons and ADU", {
  image = matrix(100, 11, 11)

  # In electrons the Poisson term is sqrt(image); in ADU it is divided by the gain
  elec = profoundMakeSigma(image = image, sky = 100, gain = 1, image_units = 'elec')
  adu = profoundMakeSigma(image = image, sky = 100, gain = 2, image_units = 'ADU',
                          output_units = 'ADU')

  expect_true(all(abs(adu - elec / 2) < 1e-10))
})

test_that("profoundResample conserves total flux by default", {
  set.seed(301)
  image = makeTestImage(dim = c(100, 100))

  coarser = profoundResample(image, pixscale_old = 1, pixscale_new = 2)
  finer = profoundResample(image, pixscale_old = 1, pixscale_new = 0.5)

  expect_identical(dim(coarser), c(50L, 50L))
  expect_identical(dim(finer), c(200L, 200L))
  expect_equal(sum(coarser), sum(image), tolerance = 1e-6)
  expect_equal(sum(finer), sum(image), tolerance = 1e-6)
})

test_that("profoundResample fluxscale options behave as documented", {
  set.seed(302)
  image = makeTestImage(dim = c(100, 100))

  norm = profoundResample(image, 1, 2, fluxscale = 'norm')
  expect_equal(sum(norm), 1, tolerance = 1e-8)

  expect_error(profoundResample(image, 1, 2, fluxscale = 'junk'), "fluxscale must be")
  expect_error(profoundResample(image, 1, 2, type = 'junk'), "type must be")
})
