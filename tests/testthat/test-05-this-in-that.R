# Tests for src/this_in_that.cpp and its R wrapper this_in_that(). The
# non-fastmatch branch of profoundProFound() was switched from %in% to
# .mat_this_in_vec_that(), so this routine's exact semantics matter.

skip_if_not_installed("imager")   # sky/RMS normalisation uses imager (Suggests)

test_that("this_in_that agrees with base R %in% on valid input", {
  set.seed(pf_seed + 61L)
  a <- as.integer(sample(0:30, 4000, replace = TRUE))
  b <- as.integer(sort(unique(sample(0:30, 8))))
  expect_identical(this_in_that(a, b), a %in% b)
  expect_identical(this_in_that(a, b, invert = TRUE), !(a %in% b))

  m <- matrix(a, 100, 40)
  expect_identical(this_in_that(m, b), matrix(a %in% b, 100, 40))

  # Duplicated values in `that` must not change the answer.
  expect_identical(this_in_that(a, c(b, b, b[1])), a %in% b)
})

test_that("the ProFound operators match their definitions", {
  set.seed(pf_seed + 62L)
  a <- as.integer(sample(0:9, 500, replace = TRUE))
  b <- as.integer(c(2, 4, 6))
  expect_identical(a %fin% b, this_in_that(a, b))
  expect_identical(a %nin% b, this_in_that(a, b, invert = TRUE))
  expect_identical(a %fin% b, a %in% b)
  expect_identical(a %nin% b, !(a %in% b))
})

test_that("this_in_that returns NA for NA and negative values of `this`", {
  a <- as.integer(c(-5, -1, 0, 1, NA, 3, 100))
  r <- this_in_that(a, 0:3L)
  expect_identical(r, c(NA, NA, TRUE, TRUE, NA, TRUE, FALSE))
  ri <- this_in_that(a, 0:3L, invert = TRUE)
  expect_identical(ri, c(NA, NA, FALSE, FALSE, NA, FALSE, TRUE))

  # A value larger than max(that) is FALSE, not an error: the bounds check in
  # the C++ keeps out-of-range lookups out of the table.
  expect_identical(this_in_that(c(1L, 99L), 1:2L), c(TRUE, FALSE))
})

test_that("this_in_that validates its input types", {
  expect_error(this_in_that(c(1.5, 2), 1:2), "this must be integer")
  expect_error(this_in_that(1:3, c(1.5, 2)), "that must be integer")
  # An unrecognised `type` is not validated and falls through to NULL.
  expect_null(this_in_that(1:3, 1:3, type = "bogus"))
})

test_that("this_in_that type = 'which' and arr.ind behave like which()", {
  a <- as.integer(c(3, 1, 3, 2, 3))
  expect_identical(this_in_that(a, 3L, type = "which"), which(a %in% 3))
  m <- matrix(a, 5, 1)
  expect_identical(this_in_that(m, 3L, type = "which", arr.ind = TRUE),
                   which(matrix(a %in% 3, 5, 1), arr.ind = TRUE))
  # The logical form keeps matrix dimensions.
  expect_identical(dim(this_in_that(m, 3L)), dim(m))
})

test_that("this_in_that is invariant to thread count", {
  set.seed(pf_seed + 63L)
  a <- as.integer(sample(0:20, 5000, replace = TRUE))
  b <- as.integer(c(0, 5, 11, 20))
  ref <- this_in_that(a, b)
  refi <- this_in_that(a, b, invert = TRUE)
  for (nt in c(2, 4, 8)) {
    expect_identical(this_in_that(a, b, nthreads = nt), ref)
    expect_identical(this_in_that(a, b, invert = TRUE, nthreads = nt), refi)
  }
  m <- matrix(a, 100, 50)
  refm <- this_in_that(m, b)
  for (nt in c(2, 4, 8)) {
    expect_identical(this_in_that(m, b, nthreads = nt), refm)
  }
})

test_that("mat_this_in_vec_that matches this_in_that on matrices", {
  # profoundProFound() calls .mat_this_in_vec_that() directly in its non-
  # fastmatch dilation branch, so check it against the R wrapper and base R.
  set.seed(pf_seed + 64L)
  a <- as.integer(sample(0:15, 1200, replace = TRUE))
  m <- matrix(a, 40, 30)
  b <- as.integer(c(2, 4, 9, 15))
  # The C++ version keeps the matrix shape; base R %in% strips it.
  expect_identical(pfns(".mat_this_in_vec_that")(m, b), this_in_that(m, b))
  expect_identical(dim(pfns(".mat_this_in_vec_that")(m, b)), dim(m))
  expect_equal(as.vector(pfns(".mat_this_in_vec_that")(m, b)), as.vector(m %in% b))
})

test_that("degenerate `that` vectors error rather than silently misbehave", {
  # max() of these is negative or NA, so the C++ lookup table would need a
  # negative length and Rcpp raises. Identical in v1.34.5 and the current tree.
  a <- as.integer(c(1, 2, 3))
  expect_error(this_in_that(a, integer(0)), "negative length vectors")
  expect_error(this_in_that(a, NA_integer_), "negative length vectors")
  expect_error(this_in_that(a, c(3L, NA_integer_)), "negative length vectors")
  expect_error(this_in_that(a, -5L), "negative length vectors")
})

test_that("known defects are pinned (documented, not asserted as correct)", {
  # DEFECT 1: .mat_this_in_vec_that(invert = TRUE) seeds its lookup table with
  # `invert` rather than FALSE, so every in-range entry reads TRUE and the
  # inverted result is always FALSE. The R wrapper this_in_that() avoids this
  # by always calling the vector version, so no package code path is affected.
  set.seed(pf_seed + 65L)
  a <- as.integer(sample(0:9, 200, replace = TRUE))
  m <- matrix(a, 20, 10)
  bad <- pfns(".mat_this_in_vec_that")(m, 3:7L, invert = TRUE)
  expect_false(any(bad))
  expect_true(any(this_in_that(m, 3:7L, invert = TRUE)))

  # DEFECT 2 (out of bounds): a `that` vector mixing a negative with a positive
  # value passes the table length check (max(that) + 1 is positive) but then
  # writes ref_ID at a negative index, so the result is heap-dependent and
  # differs between runs of the *same* build. It is therefore not exercised
  # here and the regression case set deliberately omits it; the all-negative
  # variant above is the case that is caught cleanly by Rcpp.
})

test_that("this_in_that handles a whole-image segim the way ProFound needs", {
  # The ProFound dilation branch needs exactly `segim_new %in% expand_segID`.
  img <- pf_make_img(64L)
  seg <- pf_water(img, pixcut = 1L)
  ids <- sort(unique(as.vector(seg)))
  expand <- as.integer(ids[ids > 0][seq_len(min(3, sum(ids > 0)))])
  via_cpp <- which(pfns(".mat_this_in_vec_that")(seg, expand))
  via_base <- which(seg %in% expand)
  expect_identical(via_cpp, via_base)
})

test_that("the non-fastmatch ProFound branch matches the fastmatch branch", {
  # profoundProFound() picks .mat_this_in_vec_that() when fastmatch is not
  # attached, and fastmatch::fmatch() when it is. Both must select exactly the
  # same pixels during segment dilation, which is the behaviour this edit
  # depends on.
  img <- pf_make_img(64L)
  run <- function() {
    profoundProFound(image = img, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                     size = 5, shape = "disc", iters = 6, threshold = 1.05,
                     box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                     verbose = FALSE, nthreads = 1L)
  }
  skip_if_not_installed("fastmatch")
  without <- run()
  expect_false("fastmatch" %in% .packages())

  suppressPackageStartupMessages(library(fastmatch))
  on.exit(detach("package:fastmatch", character.only = TRUE), add = TRUE)
  expect_true("fastmatch" %in% .packages())
  with <- run()

  expect_identical(without$segim, with$segim)
  expect_identical(without$segim_orig, with$segim_orig)
  expect_identical(without$segstats, with$segstats)
  expect_identical(without$Nseg, with$Nseg)
})
