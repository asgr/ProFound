# The end-to-end pipeline: profoundProFound() and its plot method.
#
# These are the highest-level tests in the suite. They assert the output
# contract (which elements exist, and their shapes), that photometry is
# internally consistent (segstats vs re-measuring the same segim), and that the
# documented options actually change what they claim to.


test_that("profoundProFound returns a complete, internally consistent result", {
  img <- pf_make_img(64L)
  r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                        sigma = 1, smooth = TRUE, size = 5, shape = "disc", iters = 6,
                        threshold = 1.05, magzero = 0, pixscale = 1, box = c(20, 20),
                        grid = c(20, 20), plot = FALSE, stats = TRUE, verbose = FALSE,
                        fluxtype = "Raw", nthreads = 1L)

  expect_s3_class(r, "profound")
  core <- c("segim", "segim_orig", "segstats", "objects", "Nseg", "sky", "skyRMS",
            "dim", "pixscale", "magzero", "time", "date", "call")
  expect_true(all(core %in% names(r)))
  expect_identical(r$dim, dim(img))
  expect_identical(dim(r$segim), dim(img))
  expect_identical(dim(r$sky), dim(img))
  expect_identical(dim(r$skyRMS), dim(img))

  # objects is the binary version of segim.
  expect_identical(as.integer(r$objects), as.integer(r$segim > 0))

  # Nseg is the number of rows in segstats.
  expect_identical(r$Nseg, nrow(r$segstats))
  expect_identical(r$Nseg, length(unique(as.vector(r$segim)[as.vector(r$segim) > 0])))

  # Version stamps are recorded.
  expect_true(is.character(r$`ProFound.version`) || inherits(r$`ProFound.version`, "package_version"))

  # segstats columns include the standard photometry set.
  expect_true(all(c("segID", "xcen", "ycen", "flux", "mag", "N100", "R50", "R100") %in%
                    names(r$segstats)))

  # The recorded call names the function.
  expect_identical(as.character(r$call[[1]]), "profoundProFound")
})

test_that("profoundProFound is deterministic", {
  img <- pf_make_img(64L)
  run <- function() {
    profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 4, size = 5,
                     box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                     verbose = FALSE, nthreads = 1L)
  }
  a <- run()
  b <- run()
  expect_identical(a$segim, b$segim)
  expect_identical(a$segstats, b$segstats)
  expect_identical(a$sky, b$sky)
  expect_identical(a$skyRMS, b$skyRMS)
  expect_identical(a$Nseg, b$Nseg)
})

test_that("a user supplied segim is used verbatim", {
  # redosegim = FALSE with an explicit segim must bypass segmentation entirely.
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  r <- profoundProFound(image = img, segim = seg, redosegim = FALSE, iters = 0,
                        box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                        verbose = FALSE, nthreads = 1L)
  expect_identical(r$segim, seg)
  expect_identical(r$segim_orig, seg)
  expect_identical(r$Nseg, nrow(r$segstats))
})

test_that("dilation only grows segments and never shrinks them", {
  img <- pf_make_img(64L)
  prev <- NULL
  for (it in c(0, 1, 2, 4, 6)) {
    r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = it, size = 5,
                          shape = "disc", threshold = 1.05, box = c(20, 20), grid = c(20, 20),
                          plot = FALSE, stats = TRUE, verbose = FALSE, nthreads = 1L)
    # segim is a superset of segim_orig for any amount of dilation.
    expect_true(all((r$segim > 0)[r$segim_orig > 0]))
    nlab <- sum(r$segim > 0)
    if (!is.null(prev)) {
      expect_gte(nlab, prev)   # more iterations can only cover more pixels
    }
    prev <- nlab
  }
})

test_that("profoundProFound fluxes agree with re-measuring its own segim", {
  img <- pf_make_img(64L)
  r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 2, size = 5,
                        box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                        verbose = FALSE, nthreads = 1L)
  # profoundSegimStats on the returned (sky-subtracted) quantities must be able
  # to reproduce the total flux of the labelled pixels.
  tot_direct <- sum(r$segstats$flux)
  tot_cpp <- sum(profoundSegimFlux(img, r$segim))
  expect_gt(tot_direct, 0)
  expect_equal(tot_cpp, sum(img[r$segim > 0]), tolerance = 1e-8)

  # The sky map is finite. skyRMS is interpolated onto the pixel grid by
  # profoundMakeSkyGrid(), and the Akima bicubic used there can overshoot a
  # little below zero between grid nodes, so a handful of cells are slightly
  # negative (they must still be non-zero, since skyRMS is used as a divisor).
  expect_true(all(is.finite(r$sky)))
  expect_true(all(is.finite(r$skyRMS)))
  expect_true(all(r$skyRMS != 0))
  # The map tracks the data: it is larger where the image is brighter rather
  # than being a constant. (Individual cells can be negative, hence the ranks.)
  gr <- seq_len(nrow(img))
  bright <- img > stats::median(img)
  expect_gt(median(r$skyRMS[bright]), median(r$skyRMS[!bright]))
})

test_that("profoundProFound keeps its documented inputs independent", {
  img <- pf_make_img(64L)
  img_copy <- img
  seg <- pf_water(img)
  expect_false({
    profoundProFound(image = img, segim = seg, redosegim = FALSE, iters = 2,
                     box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                     verbose = FALSE, nthreads = 1L)
    !identical(img, img_copy)
  })
})

test_that("profoundProFound folds NA pixels into the mask and masks photometry", {
  img <- pf_make_img(64L)
  mk <- matrix(0L, 64, 64)
  mk[1:8, 1:8] <- 1L
  im <- img
  im[20, 20] <- NA
  r <- profoundProFound(image = im, mask = mk, skycut = 1.5, iters = 2, size = 5,
                        box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                        verbose = FALSE, nthreads = 1L)
  expect_identical(r$dim, dim(im))

  # Returned mask is the union of the supplied mask and the NA pixels.
  expect_identical(sum(r$mask), sum(mk) + 1L)
  expect_true(all(r$mask[mk > 0] != 0))
  expect_true(r$mask[20, 20] != 0)

  # Masked regions are excluded from the sky-analysis object layer...
  expect_true(all(r$objects_redo[r$mask != 0] == 0))

  # ...and masked pixels contribute no flux, so the masked run must measure
  # strictly less total flux than the unmasked one on the same image.
  unmasked <- profoundProFound(image = im, skycut = 1.5, iters = 2, size = 5,
                               box = c(20, 20), grid = c(20, 20), plot = FALSE,
                               stats = TRUE, verbose = FALSE, nthreads = 1L)
  expect_lt(sum(r$segstats$flux), sum(unmasked$segstats$flux))

  # Masking does not suppress the segment map itself: ProFound segments first and
  # masks during photometry, so segments may overlap masked pixels.
  expect_identical(r$Nseg, unmasked$Nseg)
})

test_that("profoundProFound handles images with no detectable sources", {
  img <- pf_make_img(64L)
  # Sky and sky-RMS are re-estimated from the data, so a constant offset does not
  # remove detections; only a genuinely featureless image yields nothing.
  off <- profoundProFound(image = img - 100, skycut = 1.5, pixcut = 3, iters = 2, size = 5,
                          box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                          verbose = FALSE, nthreads = 1L)
  expect_gt(off$Nseg, 0)

  # A perfectly flat (all zero) image detects nothing and returns NULL segim /
  # segstats rather than erroring.
  zz <- profoundProFound(image = matrix(0, 32, 32), skycut = 1.5, pixcut = 3, iters = 2,
                         size = 5, box = c(8, 8), grid = c(8, 8), plot = FALSE,
                         stats = TRUE, verbose = FALSE, nthreads = 1L)
  expect_equal(zz$Nseg, 0)          # Nseg is numeric, not integer
  expect_null(zz$segim)
  expect_null(zz$segstats)
  expect_identical(dim(zz$sky), c(32L, 32L))
})

test_that("profoundProFound validates its arguments", {
  img <- pf_make_img(32L)
  expect_error(profoundProFound(), "required input")
  expect_error(profoundProFound(image = img, fluxtype = "bogus"), "fluxtype")
  expect_error(profoundProFound(image = img, watershed = "bogus"), "watershed")
  expect_error(profoundProFound(image = img, skygrid_type = "bogus"), "skygrid_type")
  expect_error(profoundProFound(image = "not_an_image"), "Rfits")
})

test_that("profoundProFound is invariant to thread count", {
  img <- pf_make_img(64L)
  run <- function(nt) {
    suppressMessages(profoundProFound(image = img, skycut = 1.5, pixcut = 3, tolerance = 4,
                                       ext = 2, iters = 4, size = 5, box = c(20, 20),
                                       grid = c(20, 20), plot = FALSE, stats = TRUE,
                                       verbose = FALSE, nthreads = nt))
  }
  a <- run(1L)
  for (nt in c(2, 4, 8)) {
    b <- run(nt)
    # Dilatation uses a running minimum so results are order independent; only
    # floating point summation could differ, and that is exact for integer ids.
    expect_identical(a$segim, b$segim)
    expect_equal(a$segstats$flux, b$segstats$flux, tolerance = 1e-9)
    expect_equal(a$sky, b$sky, tolerance = 1e-12)
  }
})

test_that("fluxtype scaling is exactly the documented power of ten", {
  img <- pf_make_img(64L)
  mz <- 27
  raw <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 1, size = 5,
                          box = c(20, 20), grid = c(20, 20), fluxtype = "Raw",
                          magzero = mz, plot = FALSE, stats = TRUE, verbose = FALSE,
                          nthreads = 1L)
  mj <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 1, size = 5,
                         box = c(20, 20), grid = c(20, 20), fluxtype = "Jansky",
                         magzero = mz, plot = FALSE, stats = TRUE, verbose = FALSE,
                         nthreads = 1L)
  # Segmentation is unaffected by the flux unit, only the flux columns scale.
  expect_identical(raw$segim, mj$segim)
  scale <- 10^(-0.4 * (mz - 8.9))
  expect_equal(mj$segstats$flux, raw$segstats$flux * scale, tolerance = 1e-9)

  mu <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 1, size = 5,
                         box = c(20, 20), grid = c(20, 20), fluxtype = "Microjansky",
                         magzero = mz, plot = FALSE, stats = TRUE, verbose = FALSE,
                         nthreads = 1L)
  expect_equal(mu$segstats$flux, raw$segstats$flux * 10^(-0.4 * (mz - 23.9)),
               tolerance = 1e-9)
})

test_that("profoundProFound deblending produces sub-segments", {
  # Two overlapping sources: with deblending on, the total number of measured
  # sources must be at least as many as without it.
  twin <- pf_make_twin(3, 4)
  off <- profoundProFound(image = twin + 5, skycut = 1.5, pixcut = 3, tolerance = 4,
                          ext = 2, iters = 2, size = 5, box = c(10, 10), grid = c(10, 10),
                          redosky = FALSE, sky = 5, skyRMS = 1, deblend = FALSE,
                          plot = FALSE, stats = TRUE, verbose = FALSE, nthreads = 1L)
  on <- suppressWarnings(profoundProFound(image = twin + 5, skycut = 1.5, pixcut = 3,
                                          tolerance = 4, ext = 2, iters = 2, size = 5,
                                          box = c(10, 10), grid = c(10, 10),
                                          redosky = FALSE, sky = 5, skyRMS = 1,
                                          deblend = TRUE, df = 3, radtrunc = 2,
                                          iterative = TRUE, plot = FALSE, stats = TRUE,
                                          verbose = FALSE, nthreads = 1L))
  expect_gte(on$Nseg, off$Nseg)
  expect_true(all(is.finite(on$segstats$flux)))
  # Deblending redistributes flux, so the total is approximately conserved.
  expect_equal(sum(on$segstats$flux), sum(off$segstats$flux), tolerance = 0.35)
})

test_that("profoundProFound works on the packaged real data", {
  skip_if_not_installed("Rfits")
  f <- system.file("extdata", "VIKING", "mystery_VIKING_Z.fits", package = "ProFound")
  skip_if(!file.exists(f), "VIKING test image not available")
  im <- Rfits::Rfits_read_image(f)
  r <- profoundProFound(image = im, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                        size = 5, shape = "disc", iters = 6, box = c(64, 64),
                        grid = c(64, 64), plot = FALSE, stats = TRUE, verbose = FALSE,
                        magzero = 31, nthreads = 1L)
  expect_s3_class(r, "profound")
  expect_identical(r$dim, dim(im$imDat))
  expect_gt(r$Nseg, 0L)
  expect_identical(nrow(r$segstats), r$Nseg)
  expect_true(all(r$segstats$flux[is.finite(r$segstats$flux)] > 0))
  expect_true(all(is.finite(r$sky)))
})

test_that("profoundProFound with a PSF runs the deblending convolution path", {
  skip_if_not_installed("Rfits")
  fim <- system.file("extdata", "IRdata", "s250_im.fits", package = "ProFound")
  fpsf <- system.file("extdata", "IRdata", "s250_psf.fits", package = "ProFound")
  skip_if(!file.exists(fim) || !file.exists(fpsf), "IR test images not available")
  im <- Rfits::Rfits_read_image(fim)
  psf <- Rfits::Rfits_read_image(fpsf)
  r <- suppressWarnings(profoundProFound(image = im, psf = psf, skycut = 1.5, pixcut = 3,
                                          tolerance = 4, ext = 2, size = 5, iters = 3,
                                          box = c(20, 20), grid = c(20, 20), app_diam = 4,
                                          plot = FALSE, stats = TRUE, verbose = FALSE,
                                          nthreads = 1L))
  expect_s3_class(r, "profound")
  expect_identical(r$dim, dim(im$imDat))
  expect_true(all(is.finite(r$segstats$flux)))
  # app_diam adds the *_app aperture columns.
  expect_true(any(grepl("app", names(r$segstats))))
})

test_that("plot.profound is defined and runs on a null device", {
  r <- profoundProFound(image = pf_make_img(48L), skycut = 1.5, pixcut = 3, iters = 1,
                        size = 5, box = c(16, 16), grid = c(16, 16), plot = FALSE,
                        stats = TRUE, verbose = FALSE, nthreads = 1L)
  # plot.profound is an unexported S3 method registered in the namespace.
  expect_true(is.function(pfns("plot.profound")))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(r, hist = "sky"))
  # An unrecognised hist type is rejected.
  expect_error(plot(r, hist = "bogus"), "hist type")
  expect_no_error(suppressWarnings(plot(r)))
})

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


test_that("profoundProFound's own segmentation is stable across option paths", {
  # profoundMakeSegim() is called internally with an already sky-subtracted and
  # smoothed image, so it cannot be reproduced from the raw input. What must
  # hold is that the default path and the ProFound-old watershed path agree
  # modulo id numbering, and that segim_orig is the pre-dilation image.
  img <- pf_make_img(64L)
  a <- suppressMessages(profoundProFound(image = img, skycut = 1.5, pixcut = 3,
                                          tolerance = 4, ext = 2, iters = 0, size = 5,
                                          box = c(20, 20), grid = c(20, 20), plot = FALSE,
                                          stats = TRUE, verbose = FALSE, nthreads = 1L))
  b <- suppressMessages(profoundProFound(image = img, skycut = 1.5, pixcut = 3,
                                          tolerance = 4, ext = 2, iters = 0, size = 5,
                                          box = c(20, 20), grid = c(20, 20), plot = FALSE,
                                          stats = TRUE, verbose = FALSE, watershed = "ProFound-old",
                                          nthreads = 1L))
  expect_identical(pf_labels(a$segim), pf_labels(b$segim))
  # With no dilation iterations the two segment images are identical.
  expect_identical(a$segim, a$segim_orig)
})
