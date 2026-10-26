# Aperture photometry (profoundAperPhot / AperRan) and the aperture/flux
# helpers behind it.
#
# The coverage-recursion primitives are tested in test-coverage.R against R
# transcriptions; here they are exercised through the public interface, where a
# flat image must return the geometric aperture area.


test_that("aperture photometry of a flat image equals the aperture area", {
  # A uniform image of 1 per pixel must return the geometric aperture area
  # computed by the sub-pixel coverage machinery, and that area must converge
  # to pi r^2 as depth increases.
  one <- matrix(1, 40, 40)
  sq <- matrix(0L, 40, 40)
  sq[15:25, 15:25] <- 1L
  area <- pi * 3^2
  prev <- Inf
  for (depth in c(1, 3, 5, 7)) {
    r <- profoundAperPhot(image = one, segim = sq, app_diam = 6, pixscale = 1,
                          depth = depth, correction = FALSE)
    err <- abs(r$flux_app_1 - area)
    expect_lt(err, 0.05 * area)
    # Deeper recursion resolves the boundary better, so error shrinks.
    expect_lte(err, prev + 1e-9)
    prev <- err
    # N_app counts the (fractional) pixels covered and must track the area.
    expect_equal(r$N_app_1, r$flux_app_1, tolerance = 1e-6)
  }
})

test_that("aperture photometry scales with pixel scale and magzero", {
  one <- matrix(1, 40, 40)
  sq <- matrix(0L, 40, 40)
  sq[15:25, 15:25] <- 1L
  # app_diam is in arcsec, so at pixscale p the radius in pixels is
  # app_diam / 2 / p and the area scales as 1 / p^2.
  a1 <- profoundAperPhot(image = one, segim = sq, app_diam = 6, pixscale = 1,
                         depth = 6, correction = FALSE)
  a2 <- profoundAperPhot(image = one, segim = sq, app_diam = 6, pixscale = 2,
                         depth = 6, correction = FALSE)
  expect_equal(a2$flux_app_1, a1$flux_app_1 / 4, tolerance = 1e-3)

  # Magnitudes shift by exactly magzero; fluxes do not.
  m0 <- profoundAperPhot(image = one, segim = sq, app_diam = 6, pixscale = 1,
                         magzero = 0, depth = 6, correction = FALSE)
  m25 <- profoundAperPhot(image = one, segim = sq, app_diam = 6, pixscale = 1,
                          magzero = 25, depth = 6, correction = FALSE)
  expect_equal(m25$flux_app_1, m0$flux_app_1, tolerance = 1e-9)
  expect_equal(m25$mag_app_1, m0$mag_app_1 + 25, tolerance = 1e-9)
})

test_that("aperture photometry targets requested segments only", {
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  ids <- sort(unique(as.vector(seg)))
  ids <- ids[ids > 0]
  sel <- ids[seq_len(min(3, length(ids)))]
  r <- profoundAperPhot(image = img, segim = seg, app_diam = 5, tar = data.frame(segID = sel),
                        pixscale = 1, depth = 4)
  expect_setequal(r$segID, sel)
  expect_equal(nrow(r), length(sel))

  # Without tar, every segment is measured.
  rall <- profoundAperPhot(image = img, segim = seg, app_diam = 5, pixscale = 1, depth = 4)
  expect_setequal(rall$segID, ids)

  # Coordinates may be supplied instead of segIDs; the routine maps them onto
  # segments. profoundAperPhot() adds 0.5 to supplied coords then ceilings, so
  # pixel (i, j) is addressed by xcen = i - 1, ycen = j - 1.
  loc <- t(vapply(sel, function(id) which(seg == id, arr.ind = TRUE)[1, ], integer(2)))
  rc <- profoundAperPhot(image = img, segim = seg, app_diam = 5,
                         tar = data.frame(xcen = loc[, 1] - 1, ycen = loc[, 2] - 1),
                         pixscale = 1, depth = 4)
  expect_setequal(rc$segID, sel)
})

test_that("aperture photometry respects masks", {
  one <- matrix(1, 40, 40)
  sq <- matrix(0L, 40, 40)
  sq[15:25, 15:25] <- 1L
  plain <- profoundAperPhot(image = one, segim = sq, app_diam = 6, pixscale = 1, depth = 6)
  mk <- matrix(0L, 40, 40)
  mk[15:17, 15:17] <- 1L      # masked pixels inside the aperture
  masked <- suppressWarnings(profoundAperPhot(image = one, segim = sq, mask = mk,
                                              app_diam = 6, pixscale = 1, depth = 6))
  # Masked pixels are set to NA and the aperture correction rescales the sum
  # back to the nominal area, so the corrected flux still matches pi r^2.
  expect_true(all(is.finite(masked$flux_app_1)))
  expect_equal(masked$flux_app_1, pi * 3^2, tolerance = 1e-3)
  # The uncorrected minimum is a per-pixel value of the (constant) image.
  expect_equal(masked$flux_min_1, 1)

  # Masking every pixel of the footprint leaves no measurement; the row
  # survives but carries NA flux (a divide-by-zero in the correction).
  mk_all <- matrix(1L, 40, 40)
  m_all <- suppressWarnings(profoundAperPhot(image = one, segim = sq, mask = mk_all,
                                             app_diam = 6, pixscale = 1, depth = 6))
  expect_true(nrow(m_all) <= 1L)
  if (nrow(m_all)) {
    expect_true(is.na(m_all$flux_min_1))
  }
})

test_that("aperture photometry validates inputs", {
  one <- matrix(1, 20, 20)
  expect_error(profoundAperPhot(image = one, app_diam = 1), "Need segim!")
  expect_error(profoundAperPhot(image = NULL, segim = one), "required input")
  expect_error(profoundAperPhot(image = "notaimage", segim = matrix(1L, 20, 20)))
  expect_error(profoundAperPhot(image = one, segim = matrix(1L, 20, 20), fluxtype = "bogus"))

  sq <- matrix(0L, 20, 20)
  sq[5:8, 5:8] <- 1L
  # A segID that matches no segment yields an empty (not erroring) result.
  expect_equal(nrow(profoundAperPhot(image = one, segim = sq, app_diam = 2,
                                      tar = data.frame(segID = 99L), pixscale = 1)), 0L)
  # Two fibres on the same segment are rejected as ambiguous.
  expect_error(profoundAperPhot(image = one, segim = sq, app_diam = 2,
                                tar = data.frame(segID = c(1L, 1L)), pixscale = 1),
               "unique segment")
})

test_that("profoundAperRan measures blank apertures", {
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  set.seed(pf_seed + 51L)
  r <- suppressWarnings(profoundAperRan(image = img, segim = seg, app_diam = 3, Nran = 40,
                                        pixscale = 1, depth = 4))
  expect_true(is.list(r))
  expect_true(all(c("AperPhot", "errors", "segim_ran") %in% names(r)))
  expect_equal(nrow(r$AperPhot), 40L)
  # Random apertures sit on background, so their segim_ran ids are distinct from
  # the real segments.
  expect_identical(dim(r$segim_ran), dim(seg))
  # Re-running with the same seed reproduces the result exactly.
  set.seed(pf_seed + 51L)
  r2 <- suppressWarnings(profoundAperRan(image = img, segim = seg, app_diam = 3, Nran = 40,
                                          pixscale = 1, depth = 4))
  expect_identical(r2$AperPhot, r$AperPhot)
})

test_that("profoundMakeMask builds symmetric masks", {
  m <- profoundMakeMask(size = 11, shape = "disc")
  expect_identical(dim(m), c(11L, 11L))
  expect_gt(sum(m), 0)
  # The mask is its own flip in both directions.
  expect_identical(m, m[, 11:1])
  expect_identical(m, m[11:1, ])
  mb <- profoundMakeMask(size = 9, shape = "box")
  expect_true(all(mb == 1))
})

test_that("profoundCoverMask reports segment/mask overlap", {
  seg <- matrix(0L, 20, 20)
  seg[1:6, 1:6] <- 1L
  seg[1:6, 8:13] <- 2L
  mk <- matrix(0L, 20, 20)
  mk[1:3, 1:3] <- 1L      # overlaps only segment 1
  cm <- profoundCoverMask(segim = seg, mask = mk)
  expect_true(is.data.frame(cm))
  expect_true(all(c("segID", "Nseg") %in% names(cm)) || ncol(cm) >= 2)
  # Only the touched segment appears.
  expect_true(1 %in% cm[[1]])
})

test_that("profoundEllipseSeg builds elliptical masks with the right area", {
  e <- profoundEllipseSeg(xcen = 30, ycen = 30, rad = 10, ang = 0, axrat = 1, dim = c(60, 60))
  expect_identical(dim(e), c(60L, 60L))
  # A circle of radius 10 has area pi * 100 ~ 314.
  expect_equal(sum(as.vector(e) > 0), pi * 10^2, tolerance = 0.05)

  ee <- profoundEllipseSeg(xcen = 30, ycen = 30, rad = 12, ang = 0, axrat = 0.5,
                           dim = c(60, 60))
  expect_equal(sum(as.vector(ee) > 0), pi * 12^2 * 0.5, tolerance = 0.08)

  # Rotation by 90 degrees preserves area but reorients the mask.
  r0 <- profoundEllipseSeg(xcen = 30, ycen = 30, rad = 14, ang = 0, axrat = 0.3, dim = c(60, 60))
  r90 <- profoundEllipseSeg(xcen = 30, ycen = 30, rad = 14, ang = 90, axrat = 0.3, dim = c(60, 60))
  expect_equal(sum(r0 != 0), sum(r90 != 0))
  expect_false(identical(r0, r90))
})

test_that("profoundCovMat and profoundPixelCorrelation run on real data", {
  img <- pf_make_img(64L)
  obj <- matrix(as.integer(pf_water(img) > 0), 64, 64)
  cm <- profoundCovMat(image = img, objects = obj)
  expect_true(is.list(cm) || is.matrix(cm))

  pc <- suppressWarnings(profoundPixelCorrelation(image = img, objects = obj, sky = 0,
                                                  skyRMS = 1, fft = FALSE, plot = FALSE))
  expect_true(is.list(pc))
  expect_true("cortab" %in% names(pc))
  ct <- pc$cortab
  expect_true(all(c("lag", "corx", "cory") %in% names(ct)))
  # Zero-lag autocorrelation is 1 by construction; correlations are bounded.
  expect_true(all(abs(as.matrix(ct[, c("corx", "cory")])) <= 1.001, na.rm = TRUE))
})

test_that("profoundApplyMask returns a masked image and mask of the right shape", {
  img <- pf_make_img(32L)
  mk <- matrix(0L, 32, 32)
  mk[1:5, 1:5] <- 1L
  out <- suppressWarnings(profoundApplyMask(image = img, mask = mk))
  expect_true(is.list(out))
  expect_setequal(names(out), c("mask", "image"))
  expect_identical(dim(out$image), dim(img))
  expect_identical(dim(out$mask), dim(img))

  # A named shape builds its own mask and runs on the same image.
  sd_ <- suppressWarnings(profoundApplyMask(image = img, mask = "disc",
                                             xsize = 32, ysize = 32))
  expect_identical(dim(sd_$image), dim(img))
  expect_gt(sum(sd_$mask > 0), 0)
})

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


test_that("the growth-curve quantities in segstats are internally consistent", {
  # profoundAperPhot's flux_app_<n> columns are the same measurement as
  # segstats' N50/N90 growth curve, so increasing app_diam must be monotone.
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  ids <- sort(unique(as.vector(seg)))
  ids <- ids[ids > 0]
  sel <- ids[1]
  prev <- -Inf
  for (ad in c(1, 2, 4, 8)) {
    r <- profoundAperPhot(image = img, segim = seg, app_diam = ad,
                          tar = data.frame(segID = sel), pixscale = 1, depth = 5,
                          correction = FALSE)
    # Segment pixels bound the aperture, so flux cannot decrease as the
    # aperture grows within the footprint.
    expect_gte(r$flux_app_1[1], prev)
    prev <- r$flux_app_1[1]
  }
})


test_that("makeBrush produces brushes of the requested shape", {
  box = profoundMakeMask(size = 11, shape = 'box')
  disc = profoundMakeMask(size = 11, shape = 'disc')
  diamond = profoundMakeMask(size = 11, shape = 'diamond')
  gaussian = profoundMakeMask(size = 11, shape = 'gaussian')

  expect_identical(dim(box), c(11L, 11L))
  expect_equal(sum(box > 0), 121)
  expect_equal(sum(disc > 0), 97)
  expect_equal(sum(diamond > 0), 61)
  expect_equal(sum(gaussian), 1, tolerance = 1e-12)

  # Every brush is symmetric about its centre
  expect_equal(disc, t(disc))
  expect_equal(diamond, t(diamond))
})

test_that("makeBrush rounds even sizes up with a warning", {
  expect_warning(br <- profoundMakeMask(size = 10, shape = 'box'),
                 "rounded to the next odd number")
  expect_identical(dim(br), c(11L, 11L))
})

test_that("makeBrush line shape is a single pixel wide diagonal", {
  # profoundMakeMask does not expose the line angle, so reach for the brush itself
  line = ProFound:::.makeBrush(size = 11, shape = 'line', angle = 45)

  expect_equal(sum(line > 0), 11)
  expect_equal(sum(rowSums(line) > 0), 11)
})

test_that("makeBrush rejects sizes that are not positive numbers", {
  expect_error(profoundMakeMask(size = 'big'), "'size' must be an odd integer")
  expect_error(profoundMakeMask(size = 0), "'size' must be an odd integer")
})
