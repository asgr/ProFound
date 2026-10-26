test_that("profoundAperPhot returns flux in the expected columns", {
  set.seed(1234)
  image = makeTestImage(dim = c(101, 101))
  seg = profoundMakeSegim(image, skycut = 3)

  aper = profoundAperPhot(image = image - seg$sky, segim = seg$segim, app_diam = 30,
                          magzero = 25)

  expect_equal(nrow(aper), nrow(seg$segstats))
  expect_equal(aper$segID, sort(seg$segstats$segID))
  expect_true(all(c('flux_app_1', 'mag_app_1', 'SB_app_1', 'N_app_1', 'frac_app_1') %in% names(aper)))
  # magzero is applied
  expect_equal(aper$mag_app_1, profoundFlux2Mag(aper$flux_app_1, magzero = 25))
})

test_that("profoundAperPhot on a fully enclosed segment matches segstats flux", {
  set.seed(1235)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = 51, ycen = 51, flux = 4000, sigma = 3))
  seg = profoundMakeSegim(image, skycut = 3)

  # The aperture radius in pixels is app_diam/2/pixscale, so 40 arcsec at a
  # pixel scale of 1 is a radius of 20 pixels around the 94 pixel segment
  aper = profoundAperPhot(image = image - seg$sky, segim = seg$segim, app_diam = 40)

  expect_equal(aper$flux_app_1, seg$segstats$flux, tolerance = 1e-6)
  expect_equal(aper$frac_app_1, seg$segstats$N100 / (pi * 20^2), tolerance = 1e-6)
  expect_equal(aper$N_app_1, seg$segstats$N100)
})

test_that("profoundAperPhot aperture flux grows with aperture size", {
  set.seed(1236)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = 51, ycen = 51, flux = 4000, sigma = 4))
  seg = profoundMakeSegim(image, skycut = 3)

  aper = profoundAperPhot(image = image - seg$sky, segim = seg$segim,
                          app_diam = c(10, 20, 40))

  # Each aperture returns its own set of columns
  expect_true(all(c('flux_app_1', 'flux_app_2', 'flux_app_3') %in% names(aper)))
  expect_true(all(diff(unlist(aper[1, c('flux_app_1', 'flux_app_2', 'flux_app_3')])) >= 0))

  # A 2 arcsec aperture on a sigma = 4 pixel source captures very little flux
  expect_lt(aper$flux_app_1, aper$flux_app_3)
})

test_that("profoundAperPhot centre types both work", {
  set.seed(1237)
  image = makeTestImage(dim = c(101, 101))
  seg = profoundMakeSegim(image, skycut = 3)

  for(centype in c('mean', 'wt', 'max')){
    aper = profoundAperPhot(image = image, segim = seg$segim, app_diam = 20, centype = centype)
    expect_equal(nrow(aper), nrow(seg$segstats), label = centype)
    expect_true(all(aper$flux_app_1 > 0), label = centype)
  }

  expect_equal(
    profoundAperPhot(image, seg$segim, centype = 'mean'),
    profoundAperPhot(image, seg$segim, centype = 'wt')
  )

  # An unknown centype fails later, when the aperture centres cannot be built
  expect_error(profoundAperPhot(image, seg$segim, centype = 'junk'))
})

test_that("profoundAperPhot accepts a target table of positions", {
  set.seed(1238)
  image = makeTestImage(dim = c(101, 101))
  seg = profoundMakeSegim(image, skycut = 3)

  tar = data.frame(xcen = seg$segstats$xcen, ycen = seg$segstats$ycen)
  aper = profoundAperPhot(image = image, segim = seg$segim, app_diam = 20, tar = tar)

  expect_equal(aper$xcen, tar$xcen, tolerance = 1e-6)
  expect_equal(aper$segID, seg$segstats$segID)
})

test_that("profoundAperPhot accepts fluxtype and units", {
  set.seed(1239)
  image = makeTestImage(dim = c(101, 101))
  seg = profoundMakeSegim(image, skycut = 3)

  raw = profoundAperPhot(image, seg$segim, app_diam = 20, magzero = 25, fluxtype = 'Raw')
  jy = profoundAperPhot(image, seg$segim, app_diam = 20, magzero = 25, fluxtype = 'Jansky')

  expect_equal(jy$flux_app_1 / raw$flux_app_1, rep(10^(-0.4 * (25 - 8.9)), nrow(raw)))
  # Magnitudes are unaffected by the flux units
  expect_equal(jy$mag_app_1, raw$mag_app_1)
  expect_error(profoundAperPhot(image, seg$segim, fluxtype = 'junk'), "fluxtype must be")
})

test_that("profoundAperPhot requires an image and segim", {
  expect_error(profoundAperPhot(segim = makeTestSegim()), "Missing image")
  expect_error(profoundAperPhot(image = makeTestImage()), "Need segim")
})

test_that("profoundAperRan measures the noise of empty apertures", {
  set.seed(1240)
  image = makeTestImage(dim = c(101, 101))
  seg = profoundMakeSegim(image, skycut = 3)

  ran = profoundAperRan(image = image - seg$sky, segim = seg$segim, app_diam = 20, Nran = 20)

  expect_named(ran, c('AperPhot', 'errors', 'segim_ran'))
  expect_equal(nrow(ran$AperPhot), 20)
  expect_length(ran$errors, 1)
  # Empty apertures on a sky subtracted image scatter around zero
  expect_lt(abs(mean(ran$AperPhot$flux_app_1)), 3 * ran$errors)
  # Random positions never sit on a real segment
  expect_equal(sum(ran$segim_ran > 0 & seg$segim > 0), 0)
})

test_that("profoundAperCover integrates the fraction of each pixel in an aperture", {
  grid = expand.grid(1:20 - 0.5, 1:20 - 0.5)

  frac = profoundAperCover(grid[,1], grid[,2], cx = 10.5, cy = 10.5, rad = 5)

  expect_true(all(frac >= 0 & frac <= 1))
  # depth = 0 is a hard pixel-centre test (pixels centred on the circle count)
  centres = profoundAperCover(grid[,1], grid[,2], 10.5, 10.5, 5, depth = 0)
  expect_equal(sum(centres), sum((grid[,1] - 10.5)^2 + (grid[,2] - 10.5)^2 <= 25))

  # depth > 0 allows partial coverage, and the total area approaches pi r^2
  expect_true(any(frac > 0 & frac < 1))
  expect_equal(sum(frac), pi * 25, tolerance = 0.1)
  # Every pixel counted by depth = 0 is at least partly covered
  expect_true(all(frac[centres == 1] > 0))
})

test_that("profoundPolyCover and profoundEllipCover produce sensible areas", {
  grid = expand.grid(1:40 - 0.5, 1:40 - 0.5)

  square = profoundPolyCover(grid[,1], grid[,2], poly_x = c(5, 25, 25, 5), poly_y = c(5, 5, 25, 25))
  expect_equal(sum(square), 400)

  ellipse = profoundEllipCover(grid[,1], grid[,2], cx = 20.5, cy = 20.5, rad = 10,
                               ang = 0, axrat = 0.5)
  expect_equal(sum(ellipse), pi * 10 * 5, tolerance = 0.05)
})

test_that("profoundEllipseSeg fills a circle and an ellipse", {
  circle = profoundEllipseSeg(dim = c(21, 21), rad = 8)
  ellipse = profoundEllipseSeg(dim = c(21, 21), rad = 8, axrat = 0.5)

  expect_equal(sum(circle > 0), sum(profoundAperCover(
    expand.grid(1:21 - 0.5, 1:21 - 0.5)[,1], expand.grid(1:21 - 0.5, 1:21 - 0.5)[,2],
    10.5, 10.5, 8, depth = 0)))
  expect_lt(sum(ellipse > 0), sum(circle > 0))

  # segID labels the output, and an existing segim is added to
  labelled = profoundEllipseSeg(dim = c(21, 21), rad = 5, segID = 7L)
  expect_equal(max(labelled), 7L)

  base = matrix(0L, 21, 21)
  base[1:3, 1:3] = 3L
  combined = profoundEllipseSeg(segim = base, rad = 5)
  expect_equal(sort(unique(as.integer(combined))), c(0L, 1L, 3L))
})
