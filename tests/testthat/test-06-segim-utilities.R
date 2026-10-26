# Segmentation map post-processing: grouping, adjacency, borders, relabelling
# (Keep / Zap), merging, comparison, masked-fraction and masking.
#
# profoundSegimGroup()/SegimNear() define what "touching" means, which is what
# the deblender and the grouped photometry depend on.

test_that("profoundSegimGroup groups touching segments and keeps isolated ones apart", {
  segim = makeTestSegim()

  group = profoundSegimGroup(segim)

  expect_named(group, c('groupim', 'groupsegID'))

  groupsegID = group$groupsegID
  expect_equal(groupsegID$groupID, c(1, 3, 4))
  expect_equal(groupsegID$Ngroup, c(2, 1, 1))
  expect_equal(groupsegID$Npix, c(32, 16, 16))
  expect_equal(unlist(groupsegID$segID[groupsegID$groupID == 1]), c(1L, 2L))

  # The group image uses the minimum segID of the group as its label
  expect_equal(sort(unique(as.integer(group$groupim))), c(0L, 1L, 3L, 4L))
  expect_equal(sum(group$groupim == 1), 32)
})

test_that("profoundSegimNear lists adjacent segments", {
  segim = makeTestSegim()

  near = profoundSegimNear(segim)

  expect_equal(near$segID, 1:4)
  expect_equal(near$Nnear, c(1, 1, 0, 0))
  expect_equal(near$nearID[[1]], 2L)
  expect_equal(near$nearID[[2]], 1L)
})

test_that("profoundSegimNear uses the requested offset", {
  segim = matrix(0L, 21, 21)
  segim[5:8, 5:8] = 1L
  segim[5:8, 9:12] = 2L
  segim[5:8, 13:16] = 3L

  near1 = profoundSegimNear(segim, offset = 1)
  expect_equal(near1$Nnear, c(1, 2, 1))
  expect_equal(sapply(near1$nearID, paste, collapse = ','), c('2', '1,3', '2'))

  # At a wider offset every segment is 'near' both of the others
  expect_equal(profoundSegimNear(segim, offset = 5)$Nnear, c(2, 2, 2))
})

test_that("profoundSegimEdge keeps only the bordering layer of each segment", {
  segim = makeTestSegim()

  edge = profoundSegimEdge(segim)

  expect_identical(dim(edge), dim(segim))
  # Only the pixels with a fully interior 3x3 neighbourhood are stripped, so a
  # 4x4 block keeps a 2x2 hole of zeros in the middle and 12 pixels of edge
  expect_equal(tabulate(edge), rep(12L, 4))
  expect_true(all(edge <= segim))
  expect_equal(sum(edge > 0), 48)
  expect_equal(sum(edge[4:5, 4:5]), 0)
  # The original interior is what gets dropped
  expect_equal(sum(segim > 0) - sum(edge > 0), 4 * 4)
})

test_that("profoundSegimEdge replaces NAs and can fill the sky", {
  segim = makeTestSegim()
  segim[1, 1] = NA

  expect_false(anyNA(profoundSegimEdge(segim)))
  filled = profoundSegimEdge(segim, fill = 9L)
  # The fill value replaces everything that is not part of the retained edge
  expect_equal(sum(filled == 9L), prod(dim(segim)) - sum(profoundSegimEdge(segim) > 0))
  expect_equal(sort(unique(as.integer(filled))), c(1L, 2L, 3L, 4L, 9L))
})

test_that("profoundSegimKeep merges segments onto the minimum segID", {
  segim = makeTestSegim()

  merged = profoundSegimKeep(segim, segID_merge = list(c(1L, 2L)))

  expect_equal(sort(unique(as.integer(merged))), c(0L, 1L, 3L, 4L))
  expect_equal(tabulate(merged)[1], 32)
  expect_equal(tabulate(merged)[2], 0)
})

test_that("profoundSegimKeep applies several merge groups at once", {
  segim = makeTestSegim()

  merged = profoundSegimKeep(segim, segID_merge = list(c(1L, 2L), c(3L, 4L)))

  expect_equal(sort(unique(as.integer(merged))), c(0L, 1L, 3L))
  expect_equal(tabulate(merged)[1], 32)
  expect_equal(tabulate(merged)[3], 32)
})

test_that("profoundSegimKeep requires a list of segments", {
  expect_error(profoundSegimKeep(makeTestSegim(), segID_merge = c(1L, 2L)),
               "segID_merge must be a list")
})

test_that("profoundSegimKeep can merge whole groups", {
  segim = makeTestSegim()
  group = profoundSegimGroup(segim)

  # Group 1 holds segments 1 and 2, which collapse onto the group label
  merged = profoundSegimKeep(segim, groupim = group$groupim, groupID_merge = 1)

  expect_equal(sort(unique(as.integer(merged))), c(0L, 1L, 3L, 4L))
  expect_equal(tabulate(merged)[1], 32)
  expect_equal(tabulate(merged)[3], 16)

  # Merging a single-segment group is a no-op: group 3 is already segment 3
  single = profoundSegimKeep(segim, groupim = group$groupim, groupID_merge = 3)
  expect_equal(tabulate(single), tabulate(segim))
})

test_that("profoundMergeSegID merges overlapping groups transitively", {
  expect_equal(profoundMergeSegID(list(c(1, 2), c(2, 3), c(9))), list(c(1, 2, 3), 9))
  # Non-overlapping groups are returned untouched
  expect_equal(profoundMergeSegID(list(c(1, 2), c(3, 4))), list(c(1, 2), c(3, 4)))
})

test_that("profoundZapSegID drops the groups containing a segID", {
  expect_equal(profoundZapSegID(segID = 2, segID_merge = list(c(1, 2), c(3, 4))), list(c(3, 4)))
  expect_equal(profoundZapSegID(segID = c(1, 2), segID_merge = list(c(1, 2), c(3, 4))), list(c(3, 4)))
  # A segID that appears nowhere leaves the merge list alone
  expect_equal(profoundZapSegID(segID = 99, segID_merge = list(c(1, 2))), list(c(1, 2)))
})

test_that("profoundSegimCompare matches segments between two segmentations", {
  segim_1 = matrix(0L, 20, 20)
  segim_1[2:6, 2:6] = 1L
  segim_1[12:16, 12:16] = 2L

  segim_2 = segim_1
  segim_2[2:6, 2:6] = 3L # same object, different ID

  compare = profoundSegimCompare(segim_1, segim_2)

  expect_named(compare, c('seg_bij', 'tab_match', 'seg1_in_seg2', 'seg2_in_seg1'))
  expect_equal(nrow(compare$tab_match), 2)

  match_pair = compare$tab_match[compare$tab_match$segID_1 == 1, ]
  expect_equal(match_pair$segID_2, 3)
  # A perfect match sits at a bij value of 1
  expect_equal(compare$seg_bij[1, 3], 1)
})

test_that("profoundSegimCompare reports nothing for disjoint segmentations", {
  segim_1 = matrix(0L, 20, 20)
  segim_1[1:5, 1:5] = 1L
  segim_2 = matrix(0L, 20, 20)
  segim_2[15:19, 15:19] = 1L

  compare = profoundSegimCompare(segim_1, segim_2)

  expect_equal(nrow(compare$tab_match), 0)
})

test_that("profoundSegimMerge keeps the segment containing the image peak", {
  segim_base = matrix(0L, 21, 21)
  segim_base[3:8, 3:8] = 1L
  segim_add = matrix(0L, 21, 21)
  segim_add[5:9, 5:9] = 1L

  # Peak in the overlap: both segments share a uniqueID, so base wins outright
  peak_shared = matrix(1, 21, 21)
  peak_shared[segim_base > 0] = 3
  peak_shared[6, 6] = 500
  merged_shared = profoundSegimMerge(peak_shared, segim_base, segim_add)
  expect_equal(sort(unique(as.integer(merged_shared))), c(0L, 1L))

  # Peak only inside the add segment: base is dropped, add survives as ID 2
  peak_add = matrix(1, 21, 21)
  peak_add[segim_add > 0] = 3
  peak_add[9, 9] = 500
  merged_add = profoundSegimMerge(peak_add, segim_base, segim_add)
  expect_equal(sort(unique(as.integer(merged_add))), c(0L, 1L, 2L))
  # Add segments are renumbered above the base maximum
  expect_true(any(merged_add == 2L))
  expect_equal(sum(merged_add == 2L), sum(segim_add > 0))
})

test_that("profoundCoverMask reports the masked fraction of each segment", {
  segim = makeTestSegim()
  mask = matrix(0L, 21, 21)
  mask[3:4, 3:6] = 1L # half of segment 1

  cover = profoundCoverMask(segim, mask)

  expect_named(cover, c('segim', 'orig_N', 'mask_N', 'frac'))
  expect_equal(cover$frac[cover$segim == 1], 0.5)
  expect_equal(cover$frac[cover$segim == 2], 1)
  expect_equal(cover$orig_N[cover$segim == 1], 16)
})

test_that("profoundShareFlux redistributes flux through a share matrix", {
  stats = data.frame(segID = c(1, 2), flux = c(100, 50), flux_err = c(10, 5),
                     mag = c(20, 21), N50 = c(10, 5), N90 = c(20, 10), N100 = c(30, 15))
  sharemat = matrix(c(0.5, 0.5, 0, 1), 2, 2, dimnames = list(NULL, c('1', '2')))

  shared = profoundShareFlux(segstats = stats, sharemat = sharemat)

  expect_equal(shared$segID, 1:2)
  expect_equal(shared$flux, c(75, 50))
  expect_equal(shared$flux_err, as.numeric(sqrt(c(100, 25) %*% sharemat)))
})

test_that("profoundShareFlux rejects mismatched inputs", {
  stats = data.frame(segID = 1:2, flux = 1:2, flux_err = 1:2, mag = 1:2,
                     N50 = 1:2, N90 = 1:2, N100 = 1:2)
  expect_error(profoundShareFlux(segstats = stats, sharemat = matrix(1, 3, 3)),
               "not compatible")
})

test_that("profoundAutoMerge merges spur-like segments onto their neighbours", {
  image = matrix(1, 21, 21)
  segim = makeTestSegim()

  stats = profoundSegimStats(image = image, segim = segim, sky = 0, skyRMS = 1)
  stats = stats[, c('segID', 'xmax', 'ymax', 'cenfrac', 'sep')]
  # Segment 2 is a low central-fraction spur of the group containing 1
  stats$cenfrac = c(0.4, 0.005, 0.3, 0.02)
  stats$sep = c(1, 1, 1, 1)

  merged = profoundAutoMerge(segim, stats)
  expect_equal(merged$mergeID, 2)
  expect_equal(merged$groupID, 1)
})

test_that("profoundAutoMerge returns NULL when nothing looks like a spur", {
  image = matrix(1, 21, 21)
  segim = makeTestSegim()
  stats = profoundSegimStats(image = image, segim = segim, sky = 0, skyRMS = 1)

  expect_null(profoundAutoMerge(segim, stats))
})

test_that("profoundApplyMask writes a mask into an image", {
  skip_if_not_installed("imager")

  image = matrix(1, 51, 51)
  out = profoundApplyMask(image = image, mask = 'disc', xcen = 26, ycen = 26,
                          xsize = 21, ysize = 21)

  expect_named(out, c('mask', 'image'))
  expect_identical(dim(out$mask), dim(image))
  expect_identical(dim(out$image), dim(image))
  expect_setequal(which(is.na(out$image)), which(out$mask > 0))
  expect_gt(sum(out$mask > 0), 0)
})

test_that("profoundApplyMask rotates a mask and supports several positions", {
  skip_if_not_installed("imager")

  image = matrix(1, 101, 101)
  out = profoundApplyMask(image = image, mask = 'disc', xcen = c(30, 70), ycen = c(30, 70),
                          xsize = 21, ysize = 21)

  # The two mask positions get distinct integer labels
  expect_setequal(unique(as.integer(out$mask)), c(0L, 1L, 2L))
  expect_equal(sum(out$mask == 1), sum(out$mask == 2), tolerance = 2)
})

test_that("profoundSegimExtend adds new segments above the existing maximum", {
  set.seed(17)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(30, 70), ycen = c(50, 40), flux = c(3000, 3000), sigma = c(2, 2)))
  seg = profoundMakeSegim(image, skycut = 3)

  # Start from a segim that only knows about the first source
  segim = seg$segim
  segim[segim > 1] = 0L

  # SegimExtend internally runs profoundProFound, which messages about
  # clipping box/grid to a third of the image size
  extended = suppressMessages(profoundSegimExtend(image = image, segim = segim, skycut = 3))

  expect_equal(max(extended), 2)
  expect_true(all(extended[segim > 0] == segim[segim > 0]))
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
