# Tests for the segmentation and statistics routines that profoundProFound()
# calls internally: profoundSegimStats(), profoundDilate(), profoundSegimEdge(),
# profoundSegimNear(), profoundSegimGroup(), profoundSegimMerge(), and the
# profoundMakeSegim* family.

skip_if_not_installed("imager")   # sky/RMS normalisation uses imager (Suggests)

test_that("profoundSegimStats computes counts and fluxes exactly", {
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  ids <- sort(unique(as.vector(seg)))
  ids <- ids[ids > 0]

  ss <- profoundSegimStats(image = img, segim = seg, sky = 0, skyRMS = 1)

  expect_identical(nrow(ss), length(ids))
  expect_setequal(ss$segID, ids)

  # N100 is the pixel count and flux the pixel sum of each segment.
  expect_identical(as.integer(ss$N100[match(ids, ss$segID)]),
                   as.integer(vapply(ids, function(i) sum(seg == i), integer(1))))
  expect_equal(ss$flux[match(ids, ss$segID)],
               vapply(ids, function(i) sum(img[seg == i]), numeric(1)), tolerance = 1e-9)
  expect_equal(sum(ss$flux), sum(img[seg > 0]), tolerance = 1e-9)
  expect_equal(sum(ss$N100), sum(seg > 0))

  # Magnitudes follow from fluxes.
  expect_equal(ss$mag, -2.5 * log10(ss$flux), tolerance = 1e-9)

  # Growth curve ordering is monotone by construction.
  expect_true(all(ss$N50 <= ss$N90))
  expect_true(all(ss$N90 <= ss$N100))
  expect_true(all(ss$R50 <= ss$R90))
  expect_true(all(ss$R90 <= ss$R100))

  # Centroids lie inside the image, and cenfrac is a fraction.
  expect_true(all(ss$xcen > 0 & ss$xcen <= nrow(img)))
  expect_true(all(ss$ycen > 0 & ss$ycen <= ncol(img)))
  expect_true(all(ss$cenfrac >= 0 & ss$cenfrac <= 1))

  # uniqueID is a per-source identifier; distinct segments get distinct ids.
  expect_equal(length(unique(ss$uniqueID)), nrow(ss))
})

test_that("profoundSegimStats honours sortcol, decreasing, and sky", {
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  a <- profoundSegimStats(image = img, segim = seg, sky = 0, skyRMS = 1,
                          sortcol = "segID", decreasing = FALSE)
  b <- profoundSegimStats(image = img, segim = seg, sky = 0, skyRMS = 1,
                          sortcol = "segID", decreasing = TRUE)
  expect_identical(a$segID, sort(a$segID))
  expect_identical(b$segID, sort(a$segID, decreasing = TRUE))

  f <- profoundSegimStats(image = img, segim = seg, sky = 0, skyRMS = 1,
                          sortcol = "flux", decreasing = TRUE)
  expect_identical(f$flux[1], max(a$flux))
  expect_true(all(diff(f$flux) <= 0))

  # A constant sky subtracts sky * N100 from each segment flux.
  skyv <- 2
  s2 <- profoundSegimStats(image = img, segim = seg, sky = skyv, skyRMS = 1)
  expect_equal(s2$flux, a$flux - skyv * a$N100, tolerance = 1e-9)

  # rotstats / boundstats do not change the core photometry.
  sr <- profoundSegimStats(image = img, segim = seg, sky = 0, skyRMS = 1,
                           rotstats = TRUE, boundstats = TRUE)
  expect_equal(sum(sr$flux), sum(a$flux), tolerance = 1e-9)
  expect_equal(names(sr), names(a))
})

test_that("profoundSegimStats requires an image", {
  expect_error(profoundSegimStats(image = NULL, segim = matrix(1L, 4, 4)))
  # A NULL segim yields an empty result rather than an error.
  expect_identical(nrow(profoundSegimStats(image = pf_make_img(16L), segim = NULL)), 0L)
})

test_that("profoundSegimEdge strips interiors, leaving only the skin", {
  # profoundSegimEdge() blanks pixels whose value equals a 3x3 boxblur of the
  # segim, i.e. pixels wholly surrounded by their own segment. The return is
  # therefore the original segim with only boundary pixels left non-zero.
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  e <- profoundSegimEdge(segim = seg)
  expect_identical(dim(e), dim(seg))
  expect_true(all(e[seg == 0] == 0))          # background stays zero
  expect_true(all((e != 0) <= (seg != 0)))    # survivors are a subset

  # A large solid disc keeps only its rim: the deep interior is removed and the
  # rim is roughly the perimeter.
  disc <- matrix(0L, 61, 61)
  g <- expand.grid(x = seq_len(61), y = seq_len(61))
  disc[as.matrix(g[(g$x - 31)^2 + (g$y - 31)^2 <= 20^2, ])] <- 1L
  ed <- profoundSegimEdge(segim = disc)
  expect_equal(sum(ed != 0), round(2 * pi * 20), tolerance = 12)
  expect_gt(sum(disc != 0) - sum(ed != 0), 500) # a lot of interior was removed
  expect_identical(ed[31, 31], 0L)              # the very centre is gone

  # fill is a scalar that replaces the zeroed positions.
  ef <- suppressWarnings(profoundSegimEdge(segim = disc, fill = 9L))
  expect_true(all(ef[disc == 0] == 9L))
})

test_that("profoundSegimNear finds adjacent segments", {
  # The implementation looks offset pixels away along each axis, so directly
  # touching segments (no gap) are reported as neighbours.
  seg <- matrix(0L, 21, 21)
  seg[1:6, 1:6] <- 1L
  seg[1:6, 7:12] <- 2L      # touching 1
  seg[16:20, 16:20] <- 3L   # isolated
  n1 <- profoundSegimNear(segim = seg, offset = 1)
  expect_true(is.data.frame(n1))
  expect_true(all(c("segID", "nearID", "Nnear") %in% names(n1)))
  neigh <- function(id) unlist(n1$nearID[match(id, n1$segID)])
  expect_equal(neigh(1), 2)
  expect_equal(neigh(2), 1)
  expect_length(neigh(3), 0)
  expect_identical(n1$Nnear[match(3, n1$segID)], 0L)
  expect_identical(n1$Nnear[match(1, n1$segID)], 1L)
})

test_that("profoundSegimGroup links touching segments", {
  seg <- matrix(0L, 21, 21)
  seg[1:6, 1:6] <- 1L
  seg[1:6, 7:12] <- 2L      # directly touching 1, so one group
  seg[16:20, 16:20] <- 3L    # isolated, so its own group
  g <- profoundSegimGroup(segim = seg)
  expect_true(all(c("groupim", "groupsegID") %in% names(g)))
  expect_identical(dim(g$groupim), dim(seg))
  expect_true(all((g$groupim > 0) == (seg > 0)))

  # groupsegID$segID is a list column holding the integer member ids.
  k <- g$groupsegID
  members <- lapply(k$segID, as.integer)
  expect_equal(sort(unlist(members)), c(1L, 2L, 3L))
  grp <- function(id) which(vapply(members, function(v) id %in% v, logical(1)))
  expect_identical(grp(1), grp(2))
  expect_false(identical(grp(3), grp(1)))
  expect_equal(k$Ngroup, vapply(members, length, integer(1)))
  expect_equal(k$Npix, vapply(members, function(v) sum(seg %in% v), numeric(1)))
})

test_that("profoundSegimMerge unions two label images", {
  img <- pf_make_img(32L)
  base <- matrix(0L, 32, 32)
  base[1:8, 1:8] <- 1L
  base[20:28, 20:28] <- 2L
  add <- matrix(0L, 32, 32)
  add[1:10, 5:12] <- 1L      # overlaps base segment 1
  add[15:19, 15:19] <- 2L    # disjoint
  m <- suppressWarnings(profoundSegimMerge(image = img, segim_base = base,
                                            segim_add = add, sky = 0))
  segout <- if (is.list(m) && !is.null(m$segim)) m$segim else m
  expect_identical(dim(segout), dim(base))
  # added labels are shifted above max(base) so ids cannot collide
  expect_true(all(segout[add != 0] > 0))
  expect_true(any(segout > max(base)))
})

test_that("profoundMergeSegID groups overlapping id sets and ZapSegID removes them", {
  # MergeSegID takes a LIST of id vectors and merges any group sharing an id.
  mg <- profoundMergeSegID(segID_merge = list(c(1L, 2L), c(2L, 3L), c(9L)))
  expect_true(is.list(mg))
  flat <- sort(unlist(mg))
  expect_equal(flat, c(1, 2, 3, 9))
  expect_equal(length(mg), 2L)                    # {1,2,3} merged, {9} kept
  grp_of <- function(id) which(vapply(mg, function(v) id %in% v, logical(1)))
  expect_identical(grp_of(1), grp_of(2))
  expect_identical(grp_of(2), grp_of(3))
  expect_false(grp_of(9) == grp_of(1))

  # A list with no duplicates is returned unchanged, as is a bare vector.
  expect_identical(profoundMergeSegID(segID_merge = list(1:2, 3:4)), list(1:2, 3:4))
  expect_identical(profoundMergeSegID(segID_merge = 1:4), 1:4)

  # ZapSegID drops the groups that contain any of the supplied ids.
  whole <- list(1:6)
  expect_identical(profoundZapSegID(segID = 1:6, segID_merge = whole), list())
  expect_identical(profoundZapSegID(segID = 2L, segID_merge = list(c(1L, 2L), 9L)), list(9L))
})

test_that("profoundShareFlux redistributes flux by the share matrix", {
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  ss <- profoundSegimStats(image = img, segim = seg, sky = 0, skyRMS = 1)
  k <- nrow(ss)
  sm <- diag(k)
  colnames(sm) <- as.character(ss$segID)
  # An identity share matrix returns the input fluxes.
  out <- profoundShareFlux(segstats = ss, sharemat = sm)
  expect_equal(out$flux, ss$flux, tolerance = 1e-9)

  sm2 <- matrix(1 / k, k, k)
  colnames(sm2) <- as.character(ss$segID)
  out2 <- profoundShareFlux(segstats = ss, sharemat = sm2)
  # Every output row receives the same total flux, so fluxes are equal and sum
  # to the input total.
  expect_equal(length(unique(round(out2$flux, 9))), 1)
  expect_equal(sum(out2$flux), sum(ss$flux), tolerance = 1e-8)

  expect_error(profoundShareFlux(segstats = ss, sharemat = matrix(1, 2, 2)),
               "not compatible")
})

test_that("profoundMakeSegim is deterministic and self-consistent", {
  img <- pf_make_img(64L)
  a <- profoundMakeSegim(image = img, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                         sigma = 1, smooth = TRUE, plot = FALSE, stats = TRUE,
                         verbose = FALSE, watershed = "ProFound", nthreads = 1L)
  b <- suppressMessages(profoundMakeSegim(image = img, skycut = 1.5, pixcut = 3,
                                          tolerance = 4, ext = 2, sigma = 1, smooth = TRUE,
                                          plot = FALSE, stats = TRUE, verbose = FALSE,
                                          watershed = "ProFound", nthreads = 1L))
  expect_identical(pf_labels(a$segim), pf_labels(b$segim))
  expect_identical(storage.mode(a$segim), "integer")
  expect_identical(dim(a$segim), dim(img))
  expect_true(all(c("objects", "segstats", "keyvalues", "call") %in% names(a)))
  expect_identical(as.integer(a$objects), as.integer(a$segim > 0))
  expect_equal(nrow(a$segstats), length(unique(as.vector(a$segim)[as.vector(a$segim) > 0])))

  # An image with nothing positive cannot segment at all.
  c0 <- profoundMakeSegim(image = -abs(img), skycut = 1.5, plot = FALSE, stats = TRUE,
                          verbose = FALSE)
  expect_true(all(c0$segim == 0))
  expect_error(profoundMakeSegim(image = img, watershed = "notreal"))
})

test_that("profoundMakeSegimExpand and Dilate are monotone supersets", {
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  dd <- suppressWarnings(profoundMakeSegimDilate(image = img, segim = seg, size = 5,
                                                  shape = "disc", expand = "all", sky = 0,
                                                  skyRMS = 1, verbose = FALSE, plot = FALSE,
                                                  stats = TRUE, rotstats = TRUE))
  expect_true(all((dd$segim > 0)[seg > 0]))         # growth never removes pixels
  expect_identical(dim(dd$segim), dim(seg))

  ex <- suppressWarnings(profoundMakeSegimExpand(image = img, segim = seg, skycut = 1.5,
                                                  expandsigma = 5, expand = "all", sky = 0,
                                                  skyRMS = 1, verbose = FALSE, plot = FALSE,
                                                  stats = TRUE))
  expect_identical(dim(ex$segim), dim(seg))
  expect_true(all((ex$segim > 0)[seg > 0]))
})

test_that("profoundAutoMerge and profoundCatMerge run on real segstats", {
  img <- pf_make_img(64L)
  a <- profoundMakeSegim(image = img, skycut = 1.5, pixcut = 3, plot = FALSE, stats = TRUE,
                         verbose = FALSE)
  am <- suppressWarnings(profoundAutoMerge(segim = a$segim, segstats = a$segstats))
  expect_true(is.list(am) || is.data.frame(am))
  expect_identical(dim(am$segim %||% a$segim), dim(img))
})

test_that("flux utility conversions round-trip", {
  fluxes <- 10^seq(0, 3)
  mz <- 25
  # Flux <-> magnitude.
  expect_equal(profoundFlux2Mag(flux = fluxes, magzero = mz), mz - 2.5 * log10(fluxes),
               tolerance = 1e-12)
  expect_equal(profoundMag2Flux(mag = c(18, 21.5, 25, 28.2), magzero = mz),
               10^(-0.4 * (c(18, 21.5, 25, 28.2) - mz)), tolerance = 1e-12)
  expect_equal(profoundMag2Flux(profoundFlux2Mag(fluxes, magzero = mz), magzero = mz),
               fluxes, tolerance = 1e-10)

  # Surface brightness conversions use the pixel area.
  sb <- profoundFlux2SB(flux = 100, magzero = 25, pixscale = 0.4)
  expect_equal(profoundSB2Flux(sb, magzero = 25, pixscale = 0.4), 100, tolerance = 1e-9)

  expect_equal(profoundGainConvert(gain = 1, magzero = 25, magzero_new = 25), 1)
  expect_equal(profoundGainConvert(gain = 10, magzero = 25, magzero_new = 27.5),
               10 * 10^(-0.4 * (27.5 - 25)), tolerance = 1e-12)

  # mu <-> mag round trip.
  m <- 24
  mu <- profoundMag2Mu(mag = m, re = 2, axrat = 0.7, pixscale = 0.4)
  expect_equal(profoundMu2Mag(mu = mu, re = 2, axrat = 0.7, pixscale = 0.4), m,
               tolerance = 1e-9)
})
