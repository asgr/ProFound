test_that("profoundSegimStats recovers the flux and geometry of a uniform square", {
  image = matrix(0, 51, 51)
  image[11:20, 11:20] = 10
  segim = matrix(0L, 51, 51)
  segim[11:20, 11:20] = 1L

  stats = profoundSegimStats(image = image, segim = segim, sky = 0, skyRMS = 1, pixscale = 1)

  expect_equal(nrow(stats), 1)
  expect_equal(stats$segID, 1L)
  expect_equal(stats$flux, 1000)
  expect_equal(stats$N100, 100L)
  expect_equal(stats$xcen, 15)
  expect_equal(stats$ycen, 15)
  expect_equal(stats$axrat, 1, tolerance = 1e-6)

  # Uniform square: variance = (N^2 - 1)/12 with a 10 pixel side
  expect_equal(stats$xsd, sqrt((10^2 - 1)/12))
  expect_equal(stats$ysd, sqrt((10^2 - 1)/12))
  expect_equal(stats$corxy, 0)
})

test_that("profoundSegimStats magnitudes agree with profoundFlux2Mag", {
  image = matrix(0, 51, 51)
  image[11:20, 11:20] = 10
  image[35:38, 30:41] = 5
  segim = matrix(0L, 51, 51)
  segim[11:20, 11:20] = 1L
  segim[35:38, 30:41] = 2L

  stats = profoundSegimStats(image = image, segim = segim, sky = 0, skyRMS = 1, magzero = 25)

  expect_equal(stats$mag, profoundFlux2Mag(stats$flux, magzero = 25))
  expect_equal(sort(stats$flux), c(240, 1000))
})

test_that("profoundSegimStats subtraction of the sky changes the flux", {
  image = matrix(10, 31, 31)
  segim = matrix(0L, 31, 31)
  segim[10:14, 10:14] = 1L

  raw = profoundSegimStats(image = image, segim = segim, sky = 0, skyRMS = 1)
  sub = profoundSegimStats(image = image, segim = segim, sky = 10, skyRMS = 1)

  expect_equal(raw$flux, 250)
  expect_equal(sub$flux, 0)
  expect_equal(sub$sky_mean, 10)
  expect_equal(sub$sky_sum, 250)
})

test_that("profoundSegimStats reports NA magnitudes for non-positive flux", {
  image = matrix(10, 31, 31)
  segim = matrix(0L, 31, 31)
  segim[10:14, 10:14] = 1L

  stats = profoundSegimStats(image = image, segim = segim, sky = 10, skyRMS = 1)

  expect_true(is.na(stats$mag))
  # N50/N90 collapse onto N100 when there is no flux to rank
  expect_equal(stats$N50, stats$N100)
  expect_equal(stats$N90, stats$N100)
})

test_that("profoundSegimStats elongation and angle track an elongated bar", {
  image = matrix(0, 41, 41)
  segim = matrix(0L, 41, 41)
  segim[16:20, 6:35] = 1L
  image[segim > 0] = 5

  stats = profoundSegimStats(image = image, segim = segim, sky = 0, skyRMS = 1)

  expect_lt(stats$axrat, 0.35)
  expect_equal(stats$axrat, 5/30, tolerance = 0.02)
  # semimaj runs along y (the long axis), ang near 0 degrees
  expect_equal(stats$ang, 0, tolerance = 1)
})

test_that("profoundSegimStats flags segments touching the image border", {
  image = matrix(1, 21, 21)
  segim = matrix(0L, 21, 21)
  segim[1:5, 1:5] = 1L
  segim[9:12, 9:12] = 2L

  stats = profoundSegimStats(image = image, segim = segim, sky = 0, skyRMS = 1, boundstats = TRUE)

  expect_equal(stats$Nborder[stats$segID == 1], 10)
  expect_equal(stats$Nborder[stats$segID == 2], 0)
  expect_equal(stats$flag_border[stats$segID == 1], 3) # left + bottom
  expect_equal(stats$flag_border[stats$segID == 2], 0)
  # edge_frac is the fraction of the segment edge bordering sky, not the image
  # boundary, so the isolated interior segment is 1 and the truncated one is not
  expect_equal(stats$edge_frac[stats$segID == 2], 1)
  expect_lt(stats$edge_frac[stats$segID == 1], 1)
  expect_gt(stats$Nedge[stats$segID == 1], 0)
})

test_that("profoundSegimStats omits boundary statistics they are not requested", {
  image = matrix(1, 21, 21)
  segim = matrix(0L, 21, 21)
  segim[9:12, 9:12] = 1L

  stats = profoundSegimStats(image = image, segim = segim, sky = 0, skyRMS = 1)

  expect_true(all(is.na(stats$Nedge)))
  expect_true(all(is.na(stats$flag_border)))
})

test_that("profoundSegimStats rotstats measures reflection symmetry", {
  segim = matrix(0L, 41, 41)
  segim[15:25, 15:25] = 1L

  uniform = matrix(0, 41, 41)
  uniform[segim > 0] = 2
  symmetric = profoundSegimStats(uniform, segim, sky = 0, skyRMS = 1, rotstats = TRUE)

  expect_equal(symmetric$asymm, 0)
  # Reflecting a symmetric object onto itself reproduces its flux
  expect_equal(symmetric$flux_reflect, symmetric$flux)

  # Now make it lopsided by brightening one half
  lopsided = uniform
  lopsided[15:20, 15:25] = 10
  asymmetric = profoundSegimStats(lopsided, segim, sky = 0, skyRMS = 1, rotstats = TRUE)

  expect_gt(asymmetric$asymm, 0.1)
  expect_lt(asymmetric$flux_reflect, asymmetric$flux)

  # rotstats = FALSE reports NA rather than computing them
  plain = profoundSegimStats(uniform, segim, sky = 0, skyRMS = 1)
  expect_true(is.na(plain$asymm))
  expect_true(is.na(plain$flux_reflect))
})

test_that("profoundSegimStats propagates errors from skyRMS and gain", {
  set.seed(501)
  image = makeTestImage(dim = c(101, 101))
  segim = matrix(0L, 101, 101)
  segim[20:40, 20:40] = 1L
  segim[60:80, 55:75] = 2L

  stats = profoundSegimStats(image = image, segim = segim, sky = 100, skyRMS = 5,
                             magzero = 25, gain = 2)

  expect_true(all(stats$flux_err > 0))
  expect_equal(stats$flux_err, sqrt(stats$flux_err_sky^2 + stats$flux_err_skyRMS^2 +
                                      stats$flux_err_shot^2))
  expect_equal(stats$flux_err_shot, sqrt(stats$flux) / 2)
  expect_equal(stats$mag_err, (2.5/log(10)) * abs(stats$flux_err / stats$flux))
  # skyRMS contribution scales with sqrt(N100)
  expect_equal(stats$flux_err_skyRMS, sqrt(stats$N100) * 5)
})

test_that("profoundSegimStats sortcol controls the row order", {
  image = matrix(0, 51, 51)
  image[11:12, 11:12] = 1
  image[31:38, 31:38] = 2
  segim = matrix(0L, 51, 51)
  segim[11:12, 11:12] = 1L
  segim[31:38, 31:38] = 2L

  expect_equal(profoundSegimStats(image, segim, sky = 0, skyRMS = 1)$segID, c(1L, 2L))
  expect_equal(
    profoundSegimStats(image, segim, sky = 0, skyRMS = 1, sortcol = 'flux', decreasing = TRUE)$segID,
    c(2L, 1L)
  )
})

test_that("profoundSegimStats works without any sky information", {
  image = matrix(5, 31, 31)
  segim = matrix(0L, 31, 31)
  segim[10:14, 10:14] = 1L

  stats = profoundSegimStats(image = image, segim = segim)

  expect_equal(stats$flux, 125)
  expect_equal(stats$flux_err_sky, 0)
  expect_equal(stats$flux_err_skyRMS, 0)
})

test_that("profoundSegimStats maps to sky coordinates with keyvalues", {
  keyvalues = Rfits::Rfits_read_image(
    system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package = "ProFound"),
    ext = 1)$keyvalues

  image = matrix(0, 51, 51)
  image[20:30, 20:30] = 1
  segim = matrix(0L, 51, 51)
  segim[20:30, 20:30] = 1L

  stats = profoundSegimStats(image = image, segim = segim, sky = 0, skyRMS = 1,
                             keyvalues = keyvalues)

  expect_true(is.finite(stats$RAcen))
  expect_true(is.finite(stats$Deccen))
  # The peak sits just inside the segment, so it should be close to the centroid
  expect_equal(stats$RAmax, stats$RAcen, tolerance = 1e-3)
  expect_equal(stats$Decmax, stats$Deccen, tolerance = 1e-3)
})
