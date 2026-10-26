test_that("profoundProFound returns a complete profound object", {
  set.seed(2201)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(30, 70, 50), ycen = c(50, 50, 25), flux = c(3000, 3000, 1500), sigma = c(2, 2, 3)))

  pro = runProFound(image, skycut = 2, box = 50, grid = 50)

  expect_s3_class(pro, 'profound')
  expect_named(pro, c('segim', 'segim_orig', 'objects', 'objects_redo', 'sky', 'skyRMS',
                      'image', 'mask', 'segstats', 'Nseg', 'near', 'group', 'groupstats',
                      'haralick', 'keyvalues', 'SBlim', 'magzero', 'dim', 'pixscale',
                      'imarea', 'skyLL', 'skyChiSq', 'skyChiSqMap', 'gain', 'call', 'date',
                      'time', 'ProFound.version', 'R.version'))

  expect_identical(dim(pro$segim), dim(image))
  expect_identical(dim(pro$sky), dim(image))
  expect_identical(dim(pro$segim), dim(pro$image))
  expect_equal(pro$Nseg, 3)
  expect_equal(nrow(pro$segstats), 3)
  expect_equal(sort(pro$segstats$xcen), c(30, 50, 70), tolerance = 1)
  expect_equal(range(pro$sky), c(100, 100), tolerance = 0.5)
  expect_equal(mean(pro$skyRMS), 5, tolerance = 0.5)
  expect_true(pro$time > 0)
})

test_that("profoundProFound finds the sources it was given", {
  set.seed(2202)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(25, 75), ycen = c(30, 70), flux = c(5000, 5000), sigma = c(2, 2)))

  pro = runProFound(image, skycut = 2, box = 50, grid = 50)

  expect_equal(pro$Nseg, 2)
  expect_equal(sort(pro$segstats$xcen), c(25, 75), tolerance = 1)
  expect_equal(sort(pro$segstats$ycen), c(30, 70), tolerance = 1)
  # Recovered flux is close to the truth but isophotally truncated
  expect_equal(sort(pro$segstats$flux), c(5000, 5000), tolerance = 0.2)
})

test_that("profoundProFound segim_orig lies inside segim", {
  set.seed(2203)
  image = makeTestImage()

  pro = runProFound(image, skycut = 2, box = 50, grid = 50)

  expect_true(all(pro$segim[pro$segim_orig > 0] == pro$segim_orig[pro$segim_orig > 0]))
  expect_true(all(tabulate(pro$segim) >= tabulate(pro$segim_orig)))
  expect_identical(pro$objects, (pro$segim > 0) * 1L)
})

test_that("profoundProFound honour a user supplied segim", {
  set.seed(2204)
  image = makeTestImage()
  segim = matrix(0L, 101, 101)
  segim[20:30, 20:30] = 1L

  pro = runProFound(image, segim = segim, static_photom = TRUE,
                           sky = 100, skyRMS = 5)

  expect_equal(pro$segstats$segID, 1L)
  expect_equal(pro$segstats$N100, 121L)
  # Static photometry is measured on exactly the supplied map
  expect_equal(pro$segim, segim)
})

test_that("profoundProFound magnitude zero point scales the magnitudes", {
  set.seed(2205)
  image = makeTestImage()

  pro0 = runProFound(image, skycut = 2, box = 50, grid = 50)
  pro25 = runProFound(image, skycut = 2, box = 50, grid = 50, magzero = 25)

  expect_equal(pro25$magzero, 25)
  expect_equal(pro25$segstats$mag - pro0$segstats$mag, rep(25, nrow(pro0$segstats)))
})

test_that("profoundProFound group statistics describe touching segments", {
  set.seed(2206)
  # A close pair that the watershed splits but that share a footprint
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(45, 55), ycen = c(50, 50), flux = c(4000, 3000), sigma = c(3, 3)))

  pro = runProFound(image, skycut = 2, box = 50, grid = 50,
                           tolerance = 1, groupstats = TRUE, nearstats = TRUE)

  expect_equal(pro$Nseg, 2)
  expect_named(pro$group, c('groupim', 'groupsegID'))
  expect_equal(nrow(pro$group$groupsegID), 1)
  expect_equal(pro$group$groupsegID$Ngroup, 2)
  expect_equal(length(unlist(pro$group$groupsegID$segID)), 2)
  expect_true(is.data.frame(pro$groupstats))
  expect_true(is.data.frame(pro$near))
})

test_that("profoundProFound boundstats adds boundary columns", {
  set.seed(2207)
  image = makeTestImage()
  # Put a bright source half off the edge of the image
  image[1:15, 1:15] = image[1:15, 1:15] + 300

  pro = runProFound(image, skycut = 2, box = 50, grid = 50, boundstats = TRUE)

  expect_true(all(c('Nborder', 'flag_border', 'edge_frac', 'Nedge') %in% names(pro$segstats)))
  expect_true(any(pro$segstats$flag_border > 0))
})

test_that("profoundProFound keepim = FALSE drops the image", {
  set.seed(2208)
  image = makeTestImage()

  pro = runProFound(image, skycut = 2, box = 50, grid = 50, keepim = FALSE)

  expect_null(pro$image)
  expect_equal(pro$Nseg, 2)
})

test_that("profoundProFound respects a mask and NA pixels", {
  set.seed(2209)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(25, 75), ycen = c(50, 50), flux = c(4000, 4000), sigma = c(2, 2)))

  image_masked = image
  image_masked[1:30, ] = NA

  pro = runProFound(image_masked, skycut = 2, box = 50, grid = 50)

  expect_equal(pro$Nseg, 1)
  expect_equal(pro$segstats$xcen, 75, tolerance = 1)
  expect_true(all(is.finite(pro$segstats$flux)))
})

test_that("profoundProFound rejects bad input", {
  expect_error(profoundProFound(), "Missing image")
  expect_error(profoundProFound(image = 'junk'), "only Rfits_image FITS inputs are allowed")
})

test_that("plot.profound draws without error", {
  set.seed(2210)
  image = makeTestImage()
  pro = runProFound(image, skycut = 2, box = 50, grid = 50)

  grDevices::pdf(file = NULL)
  on.exit(grDevices::dev.off())

  expect_invisible(plot(pro))
})

test_that("profoundProFound runs on a real FITS image", {
  skip_if_not_installed('Rfits')

  image = Rfits::Rfits_read_image(
    system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package = "ProFound"))

  pro = runProFound(image, skycut = 3, box = 100, grid = 100)

  expect_s3_class(pro, 'profound')
  expect_gt(pro$Nseg, 10)
  expect_true(all(c('RAcen', 'Deccen', 'R50', 'axrat') %in% names(pro$segstats)))
  # The WCS comes through from the FITS header
  expect_true(all(pro$segstats$RAcen > 0 & pro$segstats$RAcen < 360))
  expect_true(pro$magzero == 0)
})
