test_that("profoundMakeSegim finds the expected number of sources at the right places", {
  set.seed(4493)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(30, 70, 50), ycen = c(50, 50, 25), flux = c(3000, 3000, 1500), sigma = c(2, 2, 3)))

  seg = profoundMakeSegim(image, skycut = 3, pixcut = 3)

  expect_named(seg, c('segim', 'objects', 'sky', 'skyRMS', 'segstats', 'keyvalues', 'call'))
  expect_equal(max(seg$segim), 3)
  expect_equal(nrow(seg$segstats), 3)
  expect_equal(sort(seg$segstats$xcen), c(30, 50, 70), tolerance = 1)
  expect_equal(sort(seg$segstats$ycen), c(25, 50, 50), tolerance = 1)
  expect_equal(max(seg$sky), 100, tolerance = 0.3)
  expect_equal(seg$skyRMS, 5, tolerance = 0.5)
})

test_that("profoundMakeSegim segim and objects agree", {
  set.seed(5)
  image = makeTestImage()

  seg = profoundMakeSegim(image, skycut = 3)

  expect_identical(dim(seg$segim), dim(image))
  expect_identical(seg$objects, (seg$segim > 0) * 1L)
  expect_true(is.integer(seg$segim))
})

test_that("profoundMakeSegim skycut and pixcut control the number of detections", {
  set.seed(613)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(25, 50, 75), ycen = c(30, 60, 35), flux = c(4000, 800, 300), sigma = c(2, 2, 2)))

  nsig = vapply(c(1, 2, 3, 5), function(sc) max(profoundMakeSegim(image, skycut = sc, pixcut = 3)$segim), 0)
  expect_true(all(diff(nsig) <= 0)) # higher skycut -> fewer or equal detections

  # A very high pixcut removes everything
  expect_equal(max(profoundMakeSegim(image, skycut = 3, pixcut = 1e4)$segim), 0)
})

test_that("profoundMakeSegim returns an empty result when nothing is detected", {
  set.seed(991)
  # Pure noise with a skycut too high for anything to survive
  image = matrix(stats::rnorm(51*51, mean = 100, sd = 5), 51, 51)

  seg = profoundMakeSegim(image, skycut = 10)

  expect_equal(max(seg$segim), 0)
  expect_equal(sum(seg$objects), 0)
  # An empty (zero row) catalogue, not NULL, is returned when segstats are asked for
  expect_equal(nrow(seg$segstats), 0)
})

test_that("profoundMakeSegim handles an all-constant image with a known sky", {
  # A perfectly flat image has no sky scatter for profoundSkyEst to work with, so
  # the sky has to be supplied. This is the degenerate zero-signal case.
  seg = profoundMakeSegim(matrix(0, 31, 31), sky = 0, skyRMS = 1)

  expect_equal(max(seg$segim), 0)
  expect_null(seg$segstats)
})

test_that("profoundMakeSegim respects a mask", {
  set.seed(7)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(30, 70), ycen = c(50, 40), flux = c(3000, 3000), sigma = c(2, 2)))
  # Cover the source at (30, 50) but leave the one at (70, 40) alone
  mask = matrix(0L, 101, 101)
  mask[20:40, 40:60] = 1L

  image_masked = image
  image_masked[mask > 0] = NA

  seg_all = profoundMakeSegim(image, skycut = 3)
  seg_mask = profoundMakeSegim(image_masked, mask = mask, skycut = 3)

  expect_equal(max(seg_all$segim), 2)
  expect_equal(max(seg_mask$segim), 1)
  expect_equal(seg_mask$segstats$xcen, 70, tolerance = 1)
  expect_equal(max(seg_mask$segim[mask > 0]), 0)
})

test_that("profoundMakeSegim treats NA pixels as masked", {
  set.seed(8)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(30, 70), ycen = c(50, 40), flux = c(3000, 3000), sigma = c(2, 2)))

  image_na = image
  image_na[55:101, ] = NA # swallows the source at (70, 40)

  seg = profoundMakeSegim(image_na, skycut = 3)

  expect_identical(dim(seg$segim), dim(image_na))
  expect_equal(max(seg$segim), 1)
  expect_equal(seg$segstats$xcen, 30, tolerance = 1)
  expect_equal(max(seg$segim[55:101, ]), 0)
})

test_that("profoundMakeSegim stats = FALSE skips the segstats", {
  set.seed(9)
  image = makeTestImage()

  seg = profoundMakeSegim(image, skycut = 3, stats = FALSE)

  expect_null(seg$segstats)
  expect_equal(max(seg$segim), 2)
})

test_that("profoundMakeSegim accepts a user supplied sky and skyRMS", {
  set.seed(10)
  image = makeTestImage()

  seg = profoundMakeSegim(image, skycut = 3, sky = 100, skyRMS = 5,
                          skytype = 'median')

  expect_equal(max(seg$segim), 2)
  # The supplied values are passed through unaltered
  expect_equal(range(seg$sky), c(100, 100))
  expect_equal(range(seg$skyRMS), c(5, 5))
})

test_that("profoundMakeSegim rejects unknown watershed methods", {
  expect_error(profoundMakeSegim(makeTestImage(), watershed = 'junk'),
               "watershed option must either be")
})

test_that("profoundMakeSegim agrees with the EBImage watershed", {
  set.seed(11)
  image = makeTestImage()

  native = profoundMakeSegim(image, skycut = 3)
  ebimage = profoundMakeSegim(image, skycut = 3, watershed = 'EBImage')

  expect_equal(max(ebimage$segim), max(native$segim))
  expect_equal(sort(ebimage$segstats$xcen), sort(native$segstats$xcen), tolerance = 0.5)
})

test_that("profoundMakeSegimExpand grows segments but keeps their identity", {
  set.seed(12)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(30, 70, 50), ycen = c(50, 50, 25), flux = c(3000, 3000, 1500), sigma = c(2, 2, 3)))
  seg = profoundMakeSegim(image, skycut = 3)

  expanded = profoundMakeSegimExpand(image, segim = seg$segim, skycut = 2, expandsigma = 5)

  expect_equal(max(expanded$segim), max(seg$segim))
  expect_true(all(tabulate(expanded$segim) >= tabulate(seg$segim)))
  # Every core pixel is retained
  expect_true(all(expanded$segim[seg$segim > 0] == seg$segim[seg$segim > 0]))
  expect_equal(nrow(expanded$segstats), nrow(seg$segstats))
})

test_that("profoundMakeSegimExpand with expand=NULL is a no-op", {
  set.seed(13)
  image = makeTestImage()
  seg = profoundMakeSegim(image, skycut = 3)

  same = profoundMakeSegimExpand(image, segim = seg$segim, expand = NULL)

  expect_identical(same$segim, seg$segim)
})

test_that("profoundMakeSegimExpand can expand a subset of segments", {
  set.seed(14)
  image = makeTestImage()
  seg = profoundMakeSegim(image, skycut = 3)

  expand_one = profoundMakeSegimExpand(image, segim = seg$segim, expand = 1, skycut = 2)
  expand_all = profoundMakeSegimExpand(image, segim = seg$segim, expand = 'all', skycut = 2)

  expect_gt(tabulate(expand_one$segim)[1], tabulate(seg$segim)[1])
  # Segment 2 is not in the expand list, so it keeps its original footprint
  expect_equal(tabulate(expand_one$segim)[2], tabulate(seg$segim)[2])
  expect_gt(tabulate(expand_all$segim)[2], tabulate(expand_one$segim)[2])
  # Core pixels of both segments are preserved in every case
  expect_true(all(expand_one$segim[seg$segim > 0] == seg$segim[seg$segim > 0]))
  expect_true(all(expand_all$segim[seg$segim > 0] == seg$segim[seg$segim > 0]))
})

test_that("profoundDilate grows all segments but never loses pixels", {
  segim = makeTestSegim()

  dilated = profoundDilate(segim, size = 3)

  expect_equal(sort(unique(as.integer(dilated))), sort(unique(as.integer(segim))))
  expect_true(all(tabulate(dilated) >= tabulate(segim)))
  expect_true(all(dilated[segim > 0] == segim[segim > 0]))
})

test_that("profoundDilate can dilate only selected segments", {
  segim = makeTestSegim()

  dilated = profoundDilate(segim, size = 3, expand = 3)

  expect_equal(tabulate(dilated)[3] > tabulate(segim)[3], TRUE)
  expect_equal(tabulate(dilated)[4], tabulate(segim)[4])
  expect_equal(tabulate(dilated)[1], tabulate(segim)[1])
})

test_that("profoundMakeSegimDilate returns segstats for the dilated segments", {
  set.seed(15)
  image = makeTestImage()
  seg = profoundMakeSegim(image, skycut = 3)

  dilated = profoundMakeSegimDilate(image = image, segim = seg$segim, size = 5,
                                    sky = seg$sky, skyRMS = seg$skyRMS)

  expect_named(dilated, c('segim', 'objects', 'segstats', 'keyvalues', 'call'))
  expect_true(all(dilated$segstats$N100 >= seg$segstats$N100))
  expect_equal(nrow(dilated$segstats), nrow(seg$segstats))
})

test_that("profoundMakeSegimPropagate assigns every pixel to a seed", {
  set.seed(16)
  image = makeTestImage()
  seg = profoundMakeSegim(image, skycut = 3)

  propagated = profoundMakeSegimPropagate(image = image, segim = seg$segim, sky = seg$sky)

  expect_named(propagated, c('propim', 'propim_sky'))
  expect_equal(max(propagated$propim), max(seg$segim))
  expect_true(all(propagated$propim_sky[seg$segim > 0] == 0))
})
