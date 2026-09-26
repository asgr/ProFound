skip_if_not_installed("imager")   # sky/RMS normalisation uses imager (Suggests)

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
