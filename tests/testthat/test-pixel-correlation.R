test_that("profoundPixelCorrelation returns a correlation table and FFT", {
  set.seed(306)
  image = makeTestImage(dim = c(101, 101))
  pro = runProFound(image, skycut = 2, box = 50, grid = 50)

  cor = profoundPixelCorrelation(image = image, profound = pro)

  expect_named(cor, c('cortab', 'fft', 'image_sky', 'cor_err_func'))
  expect_true(all(c('lag', 'corx', 'cory', 'corx_neg', 'corx_pos') %in% names(cor$cortab)))
  # Uncorrelated Gaussian noise has near zero pixel-to-pixel correlation
  expect_lt(max(abs(cor$cortab$corx[cor$cortab$lag < 5])), 0.4)
  expect_true(is.function(cor$cor_err_func))
  expect_equal(cor$cor_err_func(1), 1)
})

test_that("profoundCovMat returns a covariance matrix", {
  set.seed(307)
  image = makeTestImage(dim = c(101, 101))
  pro = runProFound(image, skycut = 2, box = 50, grid = 50)

  covmat = profoundCovMat(image = image, profound = pro, xsel = 1:50, ysel = 1:50)

  expect_named(covmat, c('cov_mat', 'var2cov'))
  expect_identical(dim(covmat$cov_mat), c(50L, 50L))
  expect_false(anyNA(covmat$var2cov))
})

test_that("profoundSkySplitFFT splits a sky into low and high frequency parts", {
  set.seed(308)
  image = makeTestImage(dim = c(101, 101), sources = NULL)
  # Pure noise, so nothing has to be padded with random values
  split = profoundSkySplitFFT(image = image, sky = 100, skyRMS = 5, skyscale = 20)

  expect_named(split, c('sky', 'sky_lo', 'sky_hi'))
  expect_identical(dim(split$sky), dim(image))
  # The low and high frequency parts reconstruct the sky subtracted image
  expect_equal(split$sky_lo + split$sky_hi, image - 100, tolerance = 1e-6)
  # sky_lo is the smooth background and sky_hi the small scale structure
  expect_lt(stats::sd(split$sky_lo), stats::sd(split$sky_hi))
})
