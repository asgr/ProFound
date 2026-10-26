test_that("profoundFitMagPSF recovers the fluxes of isolated point sources", {
  set.seed(606)
  psf = makeTestPSF()
  trueflux = c(1000, 500, 200)
  image = makeTestImage(dim = c(101, 101), sky = 0, skyRMS = 1,
                        sources = data.frame(xcen = 1, ycen = 1, flux = 0, sigma = 1))
  image = addTestSources(image, xcen = c(30.5, 60.5, 75.5), ycen = c(40.5, 50.5, 25.5),
                         flux = trueflux, psf = psf)

  fit = profoundFitMagPSF(xcen = c(30.5, 60.5, 75.5), ycen = c(40.5, 50.5, 25.5),
                          image = image, im_sigma = matrix(1, 101, 101),
                          psf = psf, magzero = 25)

  expect_named(fit, c('psfstats', 'origLL', 'finalLL', 'origmodel', 'finalmodel', 'image',
                      'keyvalues', 'psfstats_extra', 'profound', 'mask', 'call', 'date',
                      'time', 'ProFound.version', 'R.version'))

  expect_equal(nrow(fit$psfstats), 3)
  # Sources are sorted brightest first
  expect_equal(fit$psfstats$xcen, c(30.5, 60.5, 75.5))
  # Fluxes are recovered to within the noise (bright end is most accurate)
  expect_equal(fit$psfstats$flux, trueflux, tolerance = 0.25)
  expect_true(all(diff(fit$psfstats$flux) < 0))
  expect_equal(fit$psfstats$mag, profoundFlux2Mag(fit$psfstats$flux, magzero = 25))
  expect_true(all(fit$psfstats$flux_err > 0))
  expect_gt(fit$finalLL, fit$origLL)
})

test_that("profoundFitMagPSF handles a single source", {
  set.seed(607)
  psf = makeTestPSF()
  image = makeTestImage(dim = c(51, 51), sky = 0, skyRMS = 1,
                        sources = data.frame(xcen = 1, ycen = 1, flux = 0, sigma = 1))
  image = addTestSources(image, xcen = 26.5, ycen = 26.5, flux = 2000, psf = psf)

  fit = profoundFitMagPSF(xcen = 26.5, ycen = 26.5, image = image,
                          im_sigma = matrix(1, 51, 51), psf = psf)

  expect_equal(nrow(fit$psfstats), 1)
  expect_equal(fit$psfstats$flux, 2000, tolerance = 0.1)
})

test_that("profoundFitMagPSF validates its inputs", {
  psf = makeTestPSF()
  image = matrix(0, 51, 51)

  expect_error(profoundFitMagPSF(xcen = 26, ycen = 26, im_sigma = 1, psf = psf),
               "Missing image")
  expect_error(profoundFitMagPSF(image = image, xcen = 26, ycen = 26, im_sigma = 1),
               "Missing psf")
  expect_error(profoundFitMagPSF(image = image, xcen = 26, ycen = 26, psf = psf),
               "Need im_sigma")
  expect_error(profoundFitMagPSF(image = image, xcen = 26, im_sigma = 1, psf = psf),
               "Need xcen/ycen pair")
  expect_error(profoundFitMagPSF(image = image, ycen = 26, im_sigma = 1, psf = psf),
               "Need xcen/ycen pair")
  expect_error(profoundFitMagPSF(image = image, im_sigma = 1, psf = psf),
               "Need RAcen/Decen pair")
  expect_error(profoundFitMagPSF(image = image, xcen = c(26, 30), ycen = 26,
                                 im_sigma = 1, psf = psf),
               "same length")
})

test_that("plot.fitmagpsf draws without error", {
  skip_on_cran()
  set.seed(608)
  psf = makeTestPSF()
  image = makeTestImage(dim = c(51, 51), sky = 0, skyRMS = 1,
                        sources = data.frame(xcen = 1, ycen = 1, flux = 0, sigma = 1))
  image = addTestSources(image, xcen = 26.5, ycen = 26.5, flux = 2000, psf = psf)
  fit = profoundFitMagPSF(xcen = 26.5, ycen = 26.5, image = image,
                          im_sigma = matrix(1, 51, 51), psf = psf, fit_iters = 1)

  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())

  expect_error(plot(fit), NA)
})
