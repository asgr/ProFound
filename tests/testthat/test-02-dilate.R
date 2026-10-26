skip_if_not_installed("imager")   # sky/RMS normalisation uses imager (Suggests)

test_that("dilate matches the pure-R transcription of the original algorithm", {
  # dilate_reference() reproduces the pre-rewrite dilate_cpp() loop exactly, so
  # this is an independent oracle for the rewritten two-phase scan/scatter code.
  kerns <- list(
    box3 = pf_brush(3, "box"), box5 = pf_brush(5, "box"), box9 = pf_brush(9, "box"),
    disc3 = pf_brush(3, "disc"), disc5 = pf_brush(5, "disc"), disc7 = pf_brush(7, "disc"),
    diamond3 = pf_brush(3, "diamond"), diamond5 = pf_brush(5, "diamond"),
    gaussian5 = pf_brush(5, "Gaussian"), line3 = pf_brush(3, "line"),
    asym = { k <- matrix(0L, 5, 3); k[, 2] <- 1L; k[2, 1] <- 1L; k },
    offcentre = { k <- pf_brush(5, "box"); k[3, 3] <- 0L; k },
    corners = { k <- matrix(0L, 3, 3); k[1, 1] <- 1L; k[3, 3] <- 1L; k },
    even = matrix(1L, 4, 4)
  )
  set.seed(pf_seed + 11L)
  nfail <- 0L
  for (t in seq_len(60)) {
    nr <- sample(3:26, 1)
    nc <- sample(3:26, 1)
    dens <- sample(c(0.05, 0.15, 0.35, 0.6), 1)
    segim <- matrix(as.integer(sample(0:6, nr * nc, replace = TRUE,
                                      prob = c(1 - dens, rep(dens / 6, 6)))), nr, nc)
    kname <- sample(names(kerns), 1)
    kern <- kerns[[kname]]
    expand <- if (t %% 4 == 0) 0L else as.integer(sample(1:6, sample(1:3, 1), replace = TRUE))
    got <- pf_dilate(segim, kern, expand)
    ref <- dilate_reference(segim, kern, expand)
    if (!identical(got, ref)) {
      nfail <- nfail + 1L
      if (nfail <= 3L) {
        cat("mismatch: kern =", kname, " dim =", nr, "x", nc,
            " expand =", paste(expand, collapse = ","), "\n")
      }
    }
  }
  expect_identical(nfail, 0L)
})

test_that("dilate reproduces the documented single-source behaviour", {
  # One 3x3 block, boxed: grows by one pixel in every direction.
  s <- matrix(0L, 11, 11)
  s[5:7, 5:7] <- 1L
  out <- pf_dilate(s, pf_brush(3, "box"))
  expect_identical(sum(out > 0), 5L * 5L)      # 3x3 grown by 1 -> 5x5
  expect_identical(out[4, 4], 1L)
  expect_identical(out[8, 8], 1L)
  expect_identical(out[3, 3], 0L)
  # Original block keeps its own id.
  expect_true(all(out[5:7, 5:7] == 1L))
})

test_that("dilate keeps original labels and only fills empty cells", {
  s <- pf_make_sparse()
  out <- pf_dilate(s, pf_brush(3, "box"))
  original <- which(s != 0, arr.ind = TRUE)
  for (r in seq_len(nrow(original))) {
    i <- original[r, 1]
    j <- original[r, 2]
    expect_identical(out[i, j], s[i, j])
  }
  # No pixel that was empty can take a label it did not touch.
  grown <- which(out != 0 & s == 0, arr.ind = TRUE)
  if (nrow(grown)) {
    expect_true(all(apply(grown, 1, function(p) any(s[max(1, p[1] - 1):min(nrow(s), p[1] + 1),
                                                    max(1, p[2] - 1):min(ncol(s), p[2] + 1)] != 0))))
  }
})

test_that("dilate resolves contention deterministically to the smallest id", {
  # Two labels equidistant from an empty cell; the smaller id must win, and the
  # answer must not depend on scan order (which the rewrite changed).
  s <- matrix(0L, 9, 9)
  s[2, 5] <- 7L
  s[8, 5] <- 3L
  out <- pf_dilate(s, pf_brush(3, "box"))
  expect_identical(out[5, 5], 0L)              # equidistant: neither reaches
  s2 <- matrix(0L, 7, 7)
  s2[1, 4] <- 4L
  s2[7, 4] <- 2L
  o2 <- pf_dilate(s2, pf_brush(3, "box"))
  expect_identical(o2[4, 4], 0L)
  s3 <- matrix(0L, 5, 5)
  s3[1, 3] <- 5L
  s3[3, 3] <- 0L
  s3[5, 3] <- 2L
  o3 <- pf_dilate(s3, pf_brush(3, "box"))
  # cell (2,3) is adjacent to id 5 only; cell (4,3) adjacent to id 2 only
  expect_identical(o3[2, 3], 5L)
  expect_identical(o3[4, 3], 2L)
  # A cell touching both takes the minimum.
  s4 <- matrix(0L, 3, 3)
  s4[1, 1] <- 9L
  s4[3, 3] <- 2L
  o4 <- pf_dilate(s4, pf_brush(3, "box"))
  expect_identical(o4[1, 2], 9L)
  expect_identical(o4[2, 1], 9L)
  expect_identical(o4[2, 3], 2L)
  expect_identical(o4[3, 2], 2L)
  expect_identical(o4[2, 2], 2L)               # touches both -> minimum
})

test_that("dilate expand only grows the requested segments", {
  s <- pf_make_sparse()
  ids <- sort(unique(as.vector(s)))
  ids <- ids[ids > 0]
  keep <- ids[c(1, 3)]
  out <- pf_dilate(s, pf_brush(3, "box"), expand = as.integer(keep))
  # Non-expanded segments keep exactly their original pixels...
  for (id in setdiff(ids, keep)) {
    expect_identical(sum(out == id), sum(s == id))
    expect_identical(which(out == id, arr.ind = TRUE), which(s == id, arr.ind = TRUE))
  }
  # ...and expanded ones grow.
  for (id in keep) {
    expect_gte(sum(out == id), sum(s == id))
  }
  # expand = 0 grows everything.
  allgrow <- pf_dilate(s, pf_brush(3, "box"), expand = 0L)
  expect_gte(sum(allgrow > 0), sum(out > 0))
  # An expand list with no matching ids leaves the image unchanged.
  expect_identical(pf_dilate(s, pf_brush(5, "disc"), expand = 999L), s)
  # A leading 0 in expand means "expand everything" (the original semantics).
  expect_identical(pf_dilate(s, pf_brush(3, "box"), expand = c(0L, as.integer(keep))),
                   pf_dilate(s, pf_brush(3, "box"), expand = 0L))
})

test_that("dilate handles borders without wrapping or crashing", {
  # Labels at the very edge must not bleed to the opposite side.
  frame <- matrix(0L, 25, 25)
  frame[c(1, 25), ] <- 1L
  frame[, c(1, 25)] <- 2L
  out <- pf_dilate(frame, pf_brush(7, "box"))
  expect_identical(dim(out), dim(frame))
  expect_false(anyNA(out))
  expect_identical(out[13, 13], 0L)             # centre untouched by a 7x7 kernel at edges? (see below)

  # A single label in a corner spreads at most 3 pixels away for a 7x7 kernel.
  corner <- matrix(0L, 20, 20)
  corner[1, 1] <- 5L
  oc <- pf_dilate(corner, pf_brush(7, "box"))
  nz <- which(oc != 0, arr.ind = TRUE)
  expect_true(all(nz[, 1] <= 4 & nz[, 2] <= 4)) # clipped at the border
  expect_true(all(oc[nz] == 5L))

  # Full-width stripe: dilation must not wrap from column 25 to column 1.
  stripe <- matrix(0L, 40, 40)
  stripe[, 1] <- 1L
  os <- pf_dilate(stripe, pf_brush(5, "box"))
  expect_true(all(os[, 40] == 0L))
  expect_true(all(os[, 1:3] == 1L))

  # Degenerate inputs.
  expect_identical(pf_dilate(matrix(0L, 8, 8), pf_brush(3, "disc")), matrix(0L, 8, 8))
  tiny <- matrix(c(0L, 3L, 0L, 0L), 2, 2)
  ot <- pf_dilate(tiny, pf_brush(3, "box"))
  expect_identical(dim(ot), c(2L, 2L))
  expect_true(all(ot[tiny != 0] == tiny[tiny != 0]))
})

test_that("dilate interior fast path agrees with the boundary path", {
  # The rewritten code uses a branch-free interior loop when the kernel fits
  # away from the border, and a bounds-checked loop otherwise. Both must give
  # the same answer, so a small image (all boundary) and the same pattern
  # embedded in a large one (with interior pixels) behave consistently.
  small <- pf_make_sparse(21L)
  big <- matrix(0L, 61, 61)
  big[21:41, 21:41] <- small
  k <- pf_brush(5, "box")
  a <- pf_dilate(small, k)
  b <- pf_dilate(big, k)
  expect_identical(pf_labels(a), pf_labels(b[21:41, 21:41]))
  # The embedded copy additionally grows outwards into the surrounding empty
  # space, which is where the interior fast path does its work.
  expect_gt(sum(b != 0) - sum(b[21:41, 21:41] != 0), 0)
  expect_identical(sum(b[21:41, 21:41] != 0), sum(a != 0))
})

test_that("dilate is invariant to thread count", {
  big <- pf_make_big_sparse()
  ids <- sort(unique(as.vector(big)))
  ids <- ids[ids > 0]
  k <- pf_brush(5, "box")
  ref <- pf_dilate(big, k)
  ref_e <- pf_dilate(big, k, expand = as.integer(ids[1:3]))
  for (nt in c(2, 4, 8)) {
    expect_identical(pf_dilate(big, k, nthreads = nt), ref)
    expect_identical(pf_dilate(big, k, expand = as.integer(ids[1:3]), nthreads = nt), ref_e)
  }
})

test_that("dilate output type and NA-freedom are stable", {
  s <- pf_make_sparse()
  out <- pf_dilate(s, pf_brush(5, "disc"))
  expect_identical(storage.mode(out), "integer")
  expect_identical(dim(out), dim(s))
  expect_false(anyNA(out))
  expect_true(all(out >= 0))
})

test_that("profoundDilate wraps .dilate_cpp consistently", {
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  d <- profoundDilate(segim = seg, size = 5, shape = "disc", expand = "all", iters = 1)
  expect_identical(d, pf_dilate(seg, pf_brush(5, "disc")))
  # iters = 3 equals three single applications.
  d3 <- profoundDilate(segim = seg, size = 3, shape = "box", expand = "all", iters = 3)
  man <- pf_dilate(pf_dilate(pf_dilate(seg, pf_brush(3, "box")), pf_brush(3, "box")), pf_brush(3, "box"))
  expect_identical(d3, man)
  # Original pixels are never lost by dilation.
  expect_true(all((d > 0)[seg > 0]))
  # expand = NULL is an unsupported input: length(expand) == 0 makes the
  # internal expand[1] == "all" test fail. Same in v1.34.5 and now.
  expect_error(profoundDilate(segim = seg, expand = NULL))
  # Expanding only a subset leaves the other segments untouched.
  keep <- as.integer(sort(unique(as.vector(seg)))[1])
  ds <- profoundDilate(segim = seg, size = 3, shape = "box", expand = keep, iters = 1)
  expect_identical(sum(ds != 0), sum(pf_dilate(seg, pf_brush(3, "box"), keep) != 0))
})
