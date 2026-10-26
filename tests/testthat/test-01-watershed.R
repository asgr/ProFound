# Watershed segmentation (src/water.h) and segmentation dilation (src/dilate.cpp).
#
# The C++ in src/ was rewritten for speed after v1.34.5. These tests pin its
# *behaviour* against independent oracles: a pure-R transcription of the original
# dilate loop, the untouched water_cpp_old() legacy build, and monotonicity /
# invariance properties that must hold for any correct implementation. The
# exact v1.34.5 numbers are pinned separately by the goldens (test-goldens.R).


test_that("watershed basics: shapes, dtypes, and background handling", {
  img <- pf_make_img(64L)
  seg <- pf_water(img)

  expect_identical(dim(seg), dim(img))
  expect_identical(storage.mode(seg), "integer")
  expect_false(anyNA(seg))
  expect_gte(min(seg), 0L)

  # Nothing above the cut means nothing segmented.
  cold <- pf_water(matrix(0, 20, 20), skycut = 1.5)
  expect_true(all(cold == 0L))
  expect_identical(dim(cold), c(20L, 20L))

  # All-negative images cannot segment.
  expect_true(all(pf_water(-abs(img)) == 0L))

  # A single pixel image is handled without crashing.
  expect_identical(pf_water(matrix(7, 1, 1), skycut = 0.5, pixcut = 1L), matrix(1L, 1, 1))
  expect_true(all(pf_water(matrix(-5, 12, 12)) == 0L))
})

test_that("watershed output is a subset of pixels above the cut, and exact at pixcut = 1", {
  img <- pf_make_img(64L)
  above <- as.vector(img) > 1.5
  labelled <- as.vector(pf_water(img, pixcut = 1L)) > 0
  expect_true(all(labelled <= above))
  expect_identical(sum(labelled), sum(above))

  # Larger pixcut can only remove pixels, never add them.
  expect_lte(sum(as.vector(pf_water(img, pixcut = 10L)) > 0), sum(labelled))
})

test_that("watershed pixcut enforces a minimum segment size", {
  img <- pf_make_img(64L)
  for (pc in c(2, 3, 5, 10, 25)) {
    tab <- tabulate(as.vector(pf_water(img, skycut = 1, pixcut = pc)))
    tab <- tab[tab > 0]
    if (length(tab)) {
      expect_gte(min(tab), pc)
    }
  }
})

test_that("watershed pixcut drops only genuinely small segments", {
  # A flat image above the cut becomes one giant segment plus a few plateau
  # seeds, so the sizes lost at each pixcut are explicit and predictable.
  tiny <- matrix(0L, 30, 30)
  tiny[1, 1] <- 1L
  tiny[2:3, 5] <- 2L
  tiny[4:6, 9] <- 3L
  tiny[8:12, 14] <- 4L
  im <- tiny + 10
  n_above <- sum(im > 1.5)

  sizes_1 <- sort(tabulate(as.vector(pf_water(im, skycut = 1.5, pixcut = 1L))))
  sizes_1 <- sizes_1[sizes_1 > 0]
  expect_equal(sum(sizes_1), n_above)

  for (pc in 2:6) {
    tab <- tabulate(as.vector(pf_water(im, skycut = 1.5, pixcut = pc)))
    tab <- tab[tab > 0]
    if (length(tab)) {
      expect_gte(min(tab), pc)
    }
    expect_equal(sum(tab), n_above - sum(sizes_1[sizes_1 < pc]))
  }
})

test_that("watershed segment counts are monotone in ext and abstol", {
  img <- pf_make_img(64L)

  n_ext <- vapply(0:10, function(e) max(pf_water(img, ext = e)), integer(1))
  expect_lte(max(diff(n_ext)), 0)

  n_ab <- vapply(c(0.5, 1, 2, 4, 8, 16, 32, 64),
                 function(a) max(pf_water(img, abstol = a, pixcut = 1L)), integer(1))
  expect_lte(max(diff(n_ab)), 0)

  # With no merging tolerance every above-cut pixel becomes its own seed,
  # bounded above by the number of above-cut pixels.
  expect_lte(max(pf_water(img, abstol = 0, pixcut = 1L)), sum(img > 1.5))

  # A very large merging tolerance collapses everything into one segment.
  expect_identical(max(pf_water(img, abstol = 1e6, pixcut = 1L)), 1L)
})

test_that("watershed reltol shortcut reproduces the general pow() branch", {
  img <- pf_make_img(64L)
  twin <- pf_make_twin(3, 4)

  # within_merge_tolerance() skips pow() when reltol == 0 because pow(x, 0) is
  # exactly 1. A tiny non-zero reltol must take the pow() branch yet return the
  # same decision, since abstol * pow(., 1e-12) == abstol to within 1e-12.
  expect_identical(pf_water(img, abstol = 4, reltol = 0, pixcut = 1L),
                   pf_water(img, abstol = 4, reltol = 1e-12, pixcut = 1L))
  expect_identical(pf_water(twin, abstol = 4, reltol = 0, pixcut = 1L),
                   pf_water(twin, abstol = 4, reltol = 1e-12, pixcut = 1L))

  # reltol enters as abstol * (pixel/central)^reltol with ratio >= 1, so
  # increasing it can only merge more aggressively.
  n_rel <- vapply(c(0, 0.25, 0.5, 1, 2, 5),
                  function(r) max(pf_water(img, abstol = 4, reltol = r, pixcut = 1L)),
                  integer(1))
  expect_lte(max(diff(n_rel)), 0)
})

test_that("watershed cliptol branch behaves as documented", {
  img <- pf_make_img(64L)
  # central_pixel > cliptol forces a merge, so cliptol = 0 maximises merging
  # and cliptol = Inf disables the shortcut.
  expect_lte(max(pf_water(img, abstol = 4, cliptol = 0, pixcut = 1L)),
             max(pf_water(img, abstol = 4, cliptol = Inf, pixcut = 1L)))
  expect_identical(max(pf_water(img, abstol = 1e6, cliptol = 0, pixcut = 1L)), 1L)
})

test_that("watershed splitting and merging behave on two-peak sources", {
  expect_identical(max(pf_water(pf_make_twin(3, 4), abstol = 50, pixcut = 1L)), 1L)

  # The brightest pixel of a two-peak source is always labelled.
  tw <- pf_make_twin(3, 8)
  s <- pf_water(tw, pixcut = 1L)
  pk <- which(tw == max(tw), arr.ind = TRUE)[1, ]
  expect_gt(s[pk[1], pk[2]], 0L)

  # Well separated peaks yield at least as many segments as overlapping ones.
  expect_gte(max(pf_water(pf_make_twin(1.5, 8), pixcut = 1L)),
             max(pf_water(pf_make_twin(3, 2), pixcut = 1L)))
})

test_that("watershed ids are non-negative, gap-tolerant, and canonicalise cleanly", {
  img <- pf_make_img(64L)
  seg <- pf_water(img)
  ids <- sort(unique(as.vector(seg)))
  expect_identical(ids[1], 0L)   # 0 always denotes background
  expect_true(all(ids >= 0))
  # Merging can fully consume a segment, leaving its id unused, so gaps in the
  # labelling are legitimate; ids remain strictly increasing.
  expect_true(all(diff(ids) >= 1))
  expect_equal(max(ids), max(seg))

  # Canonical relabelling is idempotent and label-preserving.
  can <- pf_labels(seg)
  expect_identical(pf_labels(can), can)
  expect_equal(as.vector(can > 0), as.vector(seg > 0))
  expect_equal(sort(unique(as.vector(can))), 0:max(can))
})

test_that("watershed is invariant to thread count", {
  img128 <- pf_make_img(128L, seed = pf_seed + 32L)
  ref <- pf_water(img128, nthreads = 1L)
  for (nt in c(2, 4, 8)) {
    expect_identical(pf_water(img128, nthreads = nt), ref)
  }
})

test_that("watershed agrees with the untouched legacy implementation", {
  # water_cpp_old() compiles the original watershed (src/water_old.cpp), which
  # the speed-focused rewrite does not touch, making it an independent oracle
  # for the rewritten water_cpp(). Compared modulo id numbering.
  #
  # Restricted to images without exactly-equal bright pixels: the two
  # implementations resolve ties differently because the C++ sort that orders
  # pixels by brightness is not stable. That divergence is identical in
  # v1.34.5 and in the current tree (and so is covered by the goldens), rather
  # than being introduced by the rewrite.
  set.seed(pf_seed + 77L)
  mism <- 0L
  detail <- character(0)
  for (t in seq_len(20)) {
    nr <- sample(12:36, 1)
    kind <- sample(c("rand", "ramp", "twin", "blob"), 1)
    im <- switch(kind,
      rand = matrix(runif(nr * nr, 0, 100), nr, nr),
      ramp = outer(seq_len(nr), seq_len(nr), "+") + matrix(rnorm(nr * nr), nr, nr),
      twin = pf_make_twin(sample(c(1.5, 3), 1), sample(c(4, 8), 1), nr),
      blob = {
        g <- expand.grid(x = seq_len(nr), y = seq_len(nr))
        z <- rep(0, nrow(g))
        for (k in seq_len(sample(3:7, 1))) {
          xc <- runif(1, 3, nr - 3); yc <- runif(1, 3, nr - 3); w <- runif(1, 1, 4)
          z <- z + runif(1, 20, 90) * exp(-((g$x - xc)^2 + (g$y - yc)^2) / (2 * w^2))
        }
        matrix(z, nr, nr) + 2 + rnorm(nr * nr, 0, 0.3)
      })
    ab <- sample(c(0.5, 2, 4, 10), 1)
    ex <- sample(1:5, 1)
    sk <- sample(c(0.5, 1.5, 3), 1)
    pc <- sample(c(1L, 3L, 8L), 1)
    rl <- sample(c(0, 0.3, 1), 1)
    new <- pf_water(im, abstol = ab, reltol = rl, ext = ex, skycut = sk, pixcut = pc)
    old <- pfns("water_cpp_old")(image = as.vector(im), nx = as.integer(nr),
                                 ny = as.integer(nr), abstol = ab, reltol = rl,
                                 cliptol = Inf, ext = as.integer(ex), skycut = sk,
                                 pixcut = pc, verbose = FALSE, Ncheck = 1e6)
    if (!identical(pf_labels(new), pf_labels(old))) {
      mism <- mism + 1L
      detail <- c(detail, paste0("t", t, " ", kind, " ext=", ex, " abstol=", ab))
    }
  }
  expect_identical(mism, 0L)
})

test_that("watershed never labels pixels below the cut", {
  img <- pf_make_img(64L)
  seg <- pf_water(img, skycut = 3, pixcut = 1L)
  expect_true(all(as.vector(img)[as.vector(seg) > 0] > 3))
})

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

test_that("the C++ coverage helpers reproduce their R areas", {
  grid = expand.grid(1:40 - 0.5, 1:40 - 0.5)

  square = profoundPolyCover(grid[,1], grid[,2], poly_x = c(5, 25, 25, 5),
                             poly_y = c(5, 5, 25, 25), depth = 0)
  expect_equal(sum(square), 400)

  # A disjoint pixel is not covered; a pixel at the centre is
  expect_equal(profoundPolyCover(2, 2, poly_x = c(5, 25, 25, 5), poly_y = c(5, 5, 25, 25), depth = 0), 0)
  expect_equal(profoundPolyCover(10, 10, poly_x = c(5, 25, 25, 5), poly_y = c(5, 5, 25, 25), depth = 0), 1)
})

test_that("profoundAperFlux, profoundPolyFlux and profoundEllipFlux sum an image", {
  image = matrix(1, 101, 101)

  # Uniform unit image: the flux is the covered area
  expect_equal(profoundAperFlux(image, cx = 50, cy = 50, rad = 10), pi * 100, tolerance = 1e-2)
  expect_equal(profoundPolyFlux(image, poly_x = c(20, 80, 80, 20), poly_y = c(20, 20, 80, 80)), 3600)
  expect_equal(profoundEllipFlux(image, cx = 50, cy = 50, rad = 10, ang = 0, axrat = 0.5),
               pi * 10 * 5, tolerance = 1e-2)

  # Scaling the image scales the flux linearly
  expect_equal(profoundAperFlux(image * 3, 50, 50, 10), 3 * profoundAperFlux(image, 50, 50, 10))
})

test_that("profoundSegimFlux sums per segment", {
  image = matrix(0, 21, 21)
  image[3:6, 3:6] = 2  # segment 1
  image[15:18, 15:18] = 3 # segment 3

  flux = profoundSegimFlux(image, makeTestSegim())

  expect_length(flux, 4)
  expect_equal(flux[1], 32) # 16 pixels at 2
  expect_equal(flux[2], 0)
  expect_equal(flux[3], 48) # 16 pixels at 3
  expect_equal(flux[4], 0)
})
