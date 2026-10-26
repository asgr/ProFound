# Catalogue merging: profoundCatMerge / profoundAutoMerge.


test_that("profoundCatMerge replaces a merged group with the group photometry", {
  pro = makeGroupedProFound()

  expect_equal(nrow(pro$segstats), 2)
  expect_equal(pro$group$groupsegID$Ngroup, 2)

  merged = profoundCatMerge(segstats = pro$segstats, groupstats = pro$groupstats,
                            groupsegID = pro$group$groupsegID, groupID_merge = 1)

  expect_equal(nrow(merged), 1)
  expect_equal(merged$segID, 1)
  expect_equal(merged$flux, pro$groupstats$flux)
  expect_gt(merged$flux, max(pro$segstats$flux))
  expect_equal(merged$origin, 'group')
})

test_that("profoundCatMerge only merges the requested groups", {
  pro = makeGroupedProFound()

  merged = profoundCatMerge(segstats = pro$segstats, groupstats = pro$groupstats,
                            groupsegID = pro$group$groupsegID, groupID_merge = 99)

  expect_equal(merged$segID, pro$segstats$segID)
  expect_equal(merged$flux, pro$segstats$flux)
  expect_equal(merged$origin, rep('seg', nrow(pro$segstats)))
})

test_that("profoundCatMerge flags non-merged segments", {
  set.seed(305)
  seg = profoundMakeSegim(makeTestImage(dim = c(101, 101)), skycut = 3)

  flagged = profoundCatMerge(segstats = seg$segstats, flag = TRUE, rowreset = TRUE)

  expect_equal(flagged$origin, rep('seg', nrow(seg$segstats)))
  expect_identical(row.names(flagged), as.character(seq_len(nrow(flagged))))
  expect_equal(flagged$flux, seg$segstats$flux)
})

test_that("profoundCatMerge can leave the frame unflagged", {
  set.seed(306)
  seg = profoundMakeSegim(makeTestImage(dim = c(101, 101)), skycut = 3)

  unflagged = profoundCatMerge(segstats = seg$segstats, flag = FALSE)

  expect_false('origin' %in% names(unflagged))
  expect_equal(unflagged, seg$segstats)
})

test_that("profoundImBlur smooths and profoundImDiff is the residual", {
  skip_if_not_installed('imager')

  set.seed(305)
  image = makeTestImage(dim = c(101, 101))

  blurred = profoundImBlur(image, sigma = 3)
  expect_identical(dim(blurred), dim(image))
  expect_lt(stats::sd(blurred), stats::sd(image))

  diff = profoundImDiff(image, sigma = 3)
  expect_equal(diff, image - blurred, tolerance = 1e-10)
  # It is a difference image, so it straddles zero and carries almost no total flux
  expect_lt(abs(sum(diff)), 1e-4 * sum(image))
  expect_lt(min(diff), 0)
  expect_gt(max(diff), 0)

  grad = profoundImGrad(image, sigma = 3)
  expect_identical(dim(grad), dim(image))
  expect_true(all(grad >= 0))
})

test_that("profoundAutoMerge and profoundCatMerge run on real segstats", {
  img <- pf_make_img(64L)
  a <- profoundMakeSegim(image = img, skycut = 1.5, pixcut = 3, plot = FALSE, stats = TRUE,
                         verbose = FALSE)
  am <- suppressWarnings(profoundAutoMerge(segim = a$segim, segstats = a$segstats))
  expect_true(is.list(am) || is.data.frame(am))
  expect_identical(dim(am$segim %||% a$segim), dim(img))
})
