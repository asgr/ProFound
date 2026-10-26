# The shared regression case set.
#
# Each entry is a named zero-argument thunk returning a value that depends only
# on ProFound itself (all inputs are seeded). This one file drives both:
#   - tests/scripts/make-goldens.R  (produces golden/goldens.csv)
#   - tests/testthat/test-11-goldens.R (checks the build under test against it)
#
# The checked-in golden table was generated from the committed v1.34.5 source,
# whose DESCRIPTION already declares Version: 1.34.5 and whose only later
# changes are the uncommitted speed-focused edits in src/. That is therefore the
# release whose numerical behaviour the suite pins; it is NOT a GitHub release
# tag, so the reference is identified by its DESCRIPTION version and its recorded
# commit (see golden/README.txt), not by a tag.

pf_cases <- function() {
  img <- pf_make_img(64L)
  img128 <- pf_make_img(128L, seed = pf_seed + 32L)
  n <- nrow(img)
  sparse <- pf_make_sparse()
  big_sparse <- pf_make_big_sparse()
  twin_1_5_2 <- pf_make_twin(1.5, 2)
  twin_3_4 <- pf_make_twin(3, 4)
  twin_3_8 <- pf_make_twin(3, 8)

  water <- function(image, nx, ny, abstol = 4, reltol = 0, cliptol = Inf, ext = 2,
                    skycut = 1.5, pixcut = 3, nthreads = 1L) {
    pfns("water_cpp")(image = as.vector(image), nx = as.integer(nx), ny = as.integer(ny),
                      abstol = abstol, reltol = reltol, cliptol = cliptol, ext = as.integer(ext),
                      skycut = skycut, pixcut = as.integer(pixcut), verbose = FALSE,
                      Ncheck = 1e6, nthreads = as.integer(nthreads))
  }
  dilate <- function(segim, kern, expand = 0L, nthreads = 1L) {
    pfns(".dilate_cpp")(as.matrix(segim), as.matrix(kern), as.integer(expand), as.integer(nthreads))
  }

  # Normalised segim used by the dilation / stats cases.
  segim <- water(img, n, n)
  segids <- sort(unique(as.vector(segim)))
  segids <- segids[segids > 0]

  set.seed(pf_seed)
  rxy <- runif(500) * 40
  rxy2 <- runif(500) * 40
  set.seed(pf_seed + 3)
  tiT <- as.integer(sample(0:12, 5000, replace = TRUE))

  im45 <- matrix(seq_len(45 * 45) %% 13 / 13, 45, 45) * 4 + 1
  cx <- c(20.4, 31.2, 12.7)
  cy <- c(22.1, 14.3, 33.9)
  rad <- c(6, 4.5, 8)
  ang <- c(15, -40, 90)
  axr <- c(0.7, 1, 0.4)
  wt <- c(1, 0.5, 2)

  list(
    # ---- watershed: the heuristically rewritten core (water.h) ----
    water_basic = function() water(img, n, n),
    water_pixcut0 = function() water(img, n, n, skycut = 1, pixcut = 0),
    water_pixcut1 = function() water(img, n, n, skycut = 1, pixcut = 1),
    water_pixcut20 = function() water(img, n, n, skycut = 1, pixcut = 20),
    water_pixcut500 = function() water(img, n, n, skycut = 1, pixcut = 500),
    water_ext0 = function() water(img, n, n, ext = 0),
    water_ext1 = function() water(img, n, n, ext = 1),
    water_ext5 = function() water(img, n, n, ext = 5),
    water_ext12 = function() water(img, n, n, ext = 12),
    water_abstol0 = function() water(img, n, n, abstol = 0, skycut = 1, pixcut = 1),
    water_abstolBig = function() water(img, n, n, abstol = 100, skycut = 1, pixcut = 1),
    water_reltol0 = function() water(img, n, n, reltol = 0),
    water_reltol = function() water(img, n, n, abstol = 4, reltol = 0.5, cliptol = 50,
                                    skycut = 1, pixcut = 1),
    water_reltol2 = function() water(img, n, n, abstol = 2, reltol = 2, cliptol = 80,
                                     ext = 3, skycut = 0.5, pixcut = 2),
    water_reltolNeg = function() water(img, n, n, abstol = 2, reltol = -1, cliptol = 80,
                                       skycut = 1, pixcut = 1),
    water_cliptol0 = function() water(img, n, n, cliptol = 0),
    water_cliptol30 = function() water(img, n, n, cliptol = 30),
    water_neg = function() {
      v <- img
      v[v > 30] <- -v[v > 30]
      water(v, n, n)
    },
    water_allneg = function() water(-abs(img), n, n),
    water_flat = function() water(matrix(1, 10, 10), 10, 10, skycut = 0.5),
    water_single = function() water(matrix(7, 1, 1), 1, 1, skycut = 0.5, pixcut = 1),
    water_tiny = function() water(matrix(c(1, 2, 3, 4), 2, 2), 2, 2, skycut = 0.5, pixcut = 1),
    water_zero = function() water(matrix(0, 8, 8), 8, 8, skycut = 0.5),
    water_twin_1_5_2 = function() water(twin_1_5_2, 40, 40, pixcut = 1),
    water_twin_3_4 = function() water(twin_3_4, 40, 40, pixcut = 1),
    water_twin_3_8 = function() water(twin_3_8, 40, 40, pixcut = 1),
    water_twin_merge = function() water(twin_3_4, 40, 40, abstol = 50, pixcut = 1),
    water_twin_merge2 = function() water(twin_3_8, 40, 40, abstol = 50, pixcut = 1),
    water_nthreads4 = function() water(img, n, n, nthreads = 4),
    water_128 = function() water(img128, 128, 128),
    water_128_nt = function() {
      lapply(c(1, 2, 4, 8), function(t) water(img128, 128, 128, nthreads = t))
    },
    water_old_match = function() {
      a <- water(img128, 128, 128)
      b <- pfns("water_cpp_old")(image = as.vector(img128), nx = 128L, ny = 128L,
                                 abstol = 4, reltol = 0, cliptol = Inf, ext = 2L,
                                 skycut = 1.5, pixcut = 3L, verbose = FALSE, Ncheck = 1e6)
      list(identical = identical(a, b), n_new = max(a), n_old = max(b),
           sum_new = sum(a), sum_old = sum(b))
    },
    water_sweep = function() {
      out <- list()
      for (ab in c(0.5, 2, 4, 10, 25)) {
        for (ex in c(0, 1, 2, 4, 8)) {
          for (pc in c(1, 5, 15)) {
            key <- paste0("a", ab, "_e", ex, "_p", pc)
            out[[key]] <- water(img, n, n, abstol = ab, ext = ex, pixcut = pc)
          }
        }
      }
      out
    },
    water_sweep_relclip = function() {
      out <- list()
      for (rl in c(0, 0.25, 1, 3)) {
        for (ct in c(1, 20, 60, 100, Inf)) {
          out[[paste0("r", rl, "_c", ct)]] <- water(img128, 128, 128, reltol = rl, cliptol = ct)
        }
      }
      out
    },
    water_pixcut_boundary = function() {
      tiny <- matrix(0L, 30, 30)
      tiny[1, 1] <- 1L
      tiny[2:3, 5] <- 2L
      tiny[4:6, 9] <- 3L
      tiny[8:12, 14] <- 4L
      lapply(1:6, function(pc) water(tiny + 10, 30, 30, skycut = 1.5, pixcut = pc))
    },

    # ---- dilate: rewritten scan/scatter loop ----
    dilate_all_brushes = function() {
      out <- list()
      for (shape in c("box", "disc", "diamond", "Gaussian", "line")) {
        for (size in c(3, 5, 7)) {
          key <- paste0(shape, size)
          out[[paste0("sparse_", key)]] <- dilate(sparse, pf_brush(size, shape))
          out[[paste0("segim_", key)]] <- dilate(segim, pf_brush(size, shape))
        }
      }
      out
    },
    dilate_big = function() {
      out <- list()
      for (shape in c("box", "disc", "diamond")) {
        for (size in c(3, 5, 7, 9, 11)) {
          out[[paste0(shape, size)]] <- dilate(big_sparse, pf_brush(size, shape))
        }
      }
      out
    },
    dilate_expand = function() {
      out <- list()
      out[["all0"]] <- dilate(segim, pf_brush(3, "box"), expand = 0L)
      out[["subset"]] <- dilate(segim, pf_brush(3, "box"), expand = segids[c(1, 3)])
      out[["first0"]] <- dilate(segim, pf_brush(3, "box"), expand = c(0L, segids[1:2]))
      out[["none"]] <- dilate(segim, pf_brush(5, "disc"), expand = 999L)
      out[["single"]] <- dilate(segim, pf_brush(3, "box"), expand = segids[1])
      out[["all"]] <- dilate(segim, pf_brush(3, "box"), expand = as.integer(segids))
      out
    },
    dilate_expand_big = function() {
      ids <- sort(unique(as.vector(big_sparse)))
      ids <- ids[ids > 0]
      out <- list()
      for (sub in list(ids[1:3], ids[c(1, 4, 7)], ids)) {
        out[[paste(sub, collapse = "_")]] <- dilate(big_sparse, pf_brush(5, "box"), expand = as.integer(sub))
      }
      out
    },
    dilate_centre_off = function() {
      k <- pf_brush(5, "box")
      k[3, 3] <- 0L
      k2 <- matrix(0L, 5, 5)
      k2[3, ] <- 1L
      k2[, 3] <- 1L
      k2[3, 3] <- 0L
      list(box5_offcentre = dilate(big_sparse, k),
           cross_offcentre = dilate(big_sparse, k2))
    },
    dilate_asymmetric = function() {
      k <- matrix(0L, 5, 3)
      k[, 2] <- 1L
      k[2, 1] <- 1L
      list(asym = dilate(big_sparse, k), even = dilate(big_sparse, matrix(1L, 4, 4)))
    },
    dilate_borders = function() {
      stripe <- matrix(0L, 40, 40)
      stripe[, seq(1, 40, by = 7)] <- 1L
      stripe[seq(1, 40, by = 5), ] <- 2L
      frame <- matrix(0L, 25, 25)
      frame[c(1, 25), ] <- 1L
      frame[, c(1, 25)] <- 2L
      list(stripe = dilate(stripe, pf_brush(5, "box")),
           frame = dilate(frame, pf_brush(7, "box")),
           empty = dilate(matrix(0L, 8, 8), pf_brush(3, "disc")),
           single = dilate(matrix(c(0L, 3L, 0L, 0L), 2, 2), pf_brush(3, "box")))
    },
    dilate_nthreads = function() {
      k <- pf_brush(5, "box")
      lapply(c(1, 2, 4, 8), function(t) {
        list(segim = dilate(segim, k, nthreads = t),
             big = dilate(big_sparse, k, nthreads = t),
             exp_ = dilate(big_sparse, k, expand = as.integer(sort(unique(as.vector(big_sparse)))[1:3]),
                           nthreads = t))
      })
    },

    # ---- pixel coverage / weighting / flux: quad-bound and pow shortcuts ----
    ellipcover_basic = function() {
      profoundEllipCover(rxy, rxy2, 20, 20, 7.3, ang = 23, axrat = 0.61, depth = 3L)
    },
    ellipcover_depths = function() {
      lapply(c(0, 1, 2, 3, 4, 5, 6, 8), function(d) {
        profoundEllipCover(rxy, rxy2, 20, 20, 7.3, ang = 23, axrat = 0.61, depth = as.integer(d))
      })
    },
    ellipcover_sweep = function() {
      out <- list()
      for (r in c(1.5, 3.7, 8, 15)) {
        for (a in c(0.15, 0.5, 1, 2)) {
          for (g in c(0, 17, 45, 90, 133)) {
            for (d in c(0, 2, 4)) {
              key <- paste0("r", r, "_a", a, "_g", g, "_d", d)
              out[[key]] <- profoundEllipCover(rxy, rxy2, 20.3, 19.8, r, ang = g,
                                               axrat = a, depth = as.integer(d))
            }
          }
        }
      }
      out
    },
    ellipcover_grid = function() {
      grid <- expand.grid(x = seq_len(60), y = seq_len(60))
      out <- list()
      for (d in c(0, 2, 4, 6)) {
        cv <- profoundEllipCover(grid$x, grid$y, 30.5, 30.5, 9.3, ang = 33, axrat = 0.42,
                                 depth = as.integer(d))
        out[[paste0("depth", d)]] <- list(total = sum(cv), n1 = sum(cv == 1), n0 = sum(cv == 0),
                                          analytic = pi * 9.3^2 * 0.42)
      }
      out
    },
    ellipcover_vs_weight = function() {
      grid <- expand.grid(x = seq_len(30), y = seq_len(30))
      cv <- profoundEllipCover(grid$x, grid$y, 15.5, 15.5, 7, ang = 25, axrat = 0.5, depth = 4L)
      w <- profoundEllipWeight(15.5, 15.5, 7, ang = 25, axrat = 0.5, dimx = 30L, dimy = 30L, depth = 4L)
      list(cover = as.numeric(cv), weight = as.numeric(w),
           max_abs = max(abs(as.numeric(cv) - as.numeric(w))))
    },
    apercover_basic = function() profoundAperCover(rxy, rxy2, 20, 20, 6.5, depth = 3L),
    apercover_depths = function() {
      lapply(c(0, 1, 2, 3, 4, 5, 6, 8), function(d) {
        profoundAperCover(rxy, rxy2, 20, 20, 6.5, depth = as.integer(d))
      })
    },
    apercover_grid = function() {
      grid <- expand.grid(x = seq_len(60), y = seq_len(60))
      out <- list()
      for (d in c(0, 2, 4, 6)) {
        cv <- profoundAperCover(grid$x, grid$y, 30.5, 30.5, 9.3, depth = as.integer(d))
        out[[paste0("depth", d)]] <- list(total = sum(cv), analytic = pi * 9.3^2)
      }
      out
    },
    polycover_basic = function() {
      profoundPolyCover(rxy, rxy2, c(0, 10, 12, 4, -2), c(1, 3, 9, 11, 6), depth = 3L)
    },
    polycover_grid = function() {
      grid <- expand.grid(x = seq_len(60), y = seq_len(60))
      profoundPolyCover(grid$x, grid$y, c(20, 45, 40, 25, 15), c(20, 25, 45, 40, 30), depth = 6L)
    },
    weight_sweep = function() {
      out <- list()
      specs <- list(
        exp1 = list(rad_re = rep(5, 3), nser = rep(1, 3)),
        dev = list(rad_re = rep(5, 3), nser = rep(4, 3)),
        off = list(rad_re = rep(0, 3), nser = rep(1, 3)),
        mixed = list(rad_re = c(3, 6, 2), nser = c(1, 4, 1)),
        nser_odd = list(rad_re = c(3, 6, 2), nser = c(2, 0.5, 1))
      )
      for (nm in names(specs)) {
        sp <- specs[[nm]]
        out[[paste0("ellip_", nm)]] <- profoundEllipWeight(cx, cy, rad, ang = ang, axrat = axr,
                                                           dimx = 45L, dimy = 45L, wt = wt,
                                                           rad_re = sp$rad_re, nser = sp$nser,
                                                           depth = 3L)
        out[[paste0("aper_", nm)]] <- profoundAperWeight(cx, cy, rad, dimx = 45L, dimy = 45L,
                                                          wt = wt, rad_re = sp$rad_re,
                                                          nser = sp$nser, depth = 3L)
      }
      out[["ellip_plain"]] <- profoundEllipWeight(cx, cy, rad, ang = ang, axrat = axr,
                                                   dimx = 45L, dimy = 45L, depth = 3L)
      out[["aper_plain"]] <- profoundAperWeight(cx, cy, rad, dimx = 45L, dimy = 45L, depth = 3L)
      out
    },
    radial_weight_nser = function() {
      list(
        ellip = lapply(c(0.5, 1, 2, 4), function(ns) {
          profoundEllipWeight(15.5, 15.5, 7, ang = 20, axrat = 0.6, dimx = 30L, dimy = 30L,
                              wt = 1, rad_re = 3, nser = ns, depth = 4L)
        }),
        aper = lapply(c(0.5, 1, 2, 4), function(ns) {
          profoundAperWeight(15.5, 15.5, 7, dimx = 30L, dimy = 30L, wt = 1, rad_re = 3,
                             nser = ns, depth = 4L)
        })
      )
    },
    flux_fits = function() {
      list(
        ellip = profoundEllipFlux(im45, cx, cy, rad, ang = ang, axrat = axr, wt = wt,
                                  rad_re = 5, nser = 1, deblend = FALSE, depth = 3L,
                                  iterations = 0L),
        ellip_deblend = profoundEllipFlux(im45, cx, cy, rad, ang = ang, axrat = axr, wt = wt,
                                          rad_re = 5, nser = 1, deblend = TRUE, depth = 3L,
                                          iterations = 5L),
        ellip_dev = profoundEllipFlux(im45, cx, cy, rad, ang = ang, axrat = axr, wt = wt,
                                      rad_re = 5, nser = 4, deblend = TRUE, depth = 3L,
                                      iterations = 5L),
        aper = profoundAperFlux(im45, cx, cy, rad, wt = wt, rad_re = 5, nser = 1,
                                deblend = FALSE, depth = 3L, iterations = 0L),
        aper_deblend = profoundAperFlux(im45, cx, cy, rad, wt = wt, rad_re = 5, nser = 1,
                                         deblend = TRUE, depth = 3L, iterations = 5L),
        poly = profoundPolyFlux(im45, c(0, 10, 12, 4, -2), c(1, 3, 9, 11, 6), depth = 3L)
      )
    },
    segimflux_basic = function() profoundSegimFlux(img, img),
    segimflux_seg = function() {
      profoundSegimFlux(img, matrix(as.integer(segim > 0), n, n) * 1L)
    },
    segimflux_real = function() profoundSegimFlux(img, matrix(as.integer(segim), n, n)),

    # ---- Akima interpolation: rewritten O(1) node lookup ----
    akima_grid = function() {
      gx <- seq(1, 31, length.out = 8)
      gy <- seq(2, 40, length.out = 10)
      z <- matrix(sin(outer(gx, gy) / 7) * 3 + outer(seq_along(gx), seq_along(gy)), 8, 10)
      out <- matrix(0, 24, 30)
      pfns(".interpolateAkimaGrid")(gx, gy, z, out)
      out
    },
    akima_edge = function() {
      gx <- seq(0, 7, length.out = 8)
      gy <- seq(-3, 11, length.out = 15)
      z <- matrix(outer(seq_along(gx) %% 5, seq_along(gy) %% 7), 8, 15)
      out <- matrix(0, 40, 60)
      pfns(".interpolateAkimaGrid")(gx, gy, z, out)
      out
    },
    akima_sweep = function() {
      out <- list()
      for (nx in c(5, 8, 17)) {
        for (ny in c(4, 11, 20)) {
          gx <- seq(-2, 13, length.out = nx)
          gy <- seq(0.5, 22.5, length.out = ny)
          z <- matrix(sin(outer(gx, gy) / 3.3) * 5 + outer(gx, gy / 20), nx, ny)
          o <- matrix(NA_real_, 33, 41)
          pfns(".interpolateAkimaGrid")(gx, gy, z, o)
          out[[paste0(nx, "x", ny)]] <- o
        }
      }
      out
    },
    akima_last_interval = function() {
      # The rewritten lookup clamps the interval index to mXBound - 2, so the
      # last interval of the grid is reachable. Queries right up to the final
      # node must reproduce the node values; an index clamped one interval early
      # would silently extrapolate instead.
      gx <- seq(0.5, 12.5)                 # 13 nodes, spacing 1
      gy <- seq(0.5, 9.5)                  # 10 nodes
      set.seed(pf_seed + 91L)
      z <- matrix(rnorm(13 * 10, 50, 8), 13, 10)
      out <- matrix(NA_real_, 13, 10)      # cells land exactly on every node
      ProFound:::.interpolateAkimaGrid(gx, gy, z, out)
      list(at_nodes = out, truth = z, last_col = out[, 10], last_row = out[13, ])
    },
    akima_offsets = function() {
      out <- list()
      for (off in c(0, 0.37, 1.5, -1.5, 0.999999)) {
        gx <- seq(1, 9, by = 1) + off
        gy <- seq(2, 12, by = 1) - off
        z <- matrix(outer(seq_along(gx), seq_along(gy)) / 7, 9, 11)
        o <- matrix(NA_real_, 25, 19)
        pfns(".interpolateAkimaGrid")(gx, gy, z, o)
        out[[as.character(off)]] <- o
      }
      out
    },
    linear_grid = function() {
      gx <- seq(1, 20, length.out = 6)
      gy <- seq(1, 30, length.out = 9)
      set.seed(pf_seed + 4)
      z <- matrix(rnorm(54, 10, 2), 6, 9)
      out <- matrix(0, 15, 25)
      pfns(".interpolateLinearGrid")(gx, gy, z, out)
      out
    },
    resample_sweep = function() {
      out <- list()
      for (po in c(0.5, 1, 2)) {
        for (pn in c(0.5, 1, 2)) {
          for (ty in c("bilinear", "bicubic")) {
            for (fs in c("image", "pixscale", "norm")) {
              key <- paste0(po, "_", pn, "_", ty, "_", fs)
              out[[key]] <- tryCatch(
                profoundResample(img, pixscale_old = po, pixscale_new = pn, type = ty,
                                 fluxscale = fs),
                error = function(e) pf_error_val(e))
            }
          }
        }
      }
      out
    },

    # ---- this_in_that: used by the rewritten non-fastmatch branch of ProFound ----
    tit_vec = function() this_in_that(tiT, 3:6L),
    tit_vec_inv = function() this_in_that(tiT, 3:6L, invert = TRUE),
    tit_vec_all = function() this_in_that(tiT, 0:12),
    tit_vec_none = function() this_in_that(tiT, 99:120),
    tit_mat = function() this_in_that(matrix(tiT, 100, 50), as.integer(c(2, 5, 9))),
    tit_which = function() this_in_that(tiT, as.integer(c(1, 4, 11)), type = "which"),
    tit_which_ai = function() {
      this_in_that(matrix(tiT, 100, 50), as.integer(c(1, 4, 11)), type = "which", arr.ind = TRUE)
    },
    tit_ops = function() list(fin = tiT %fin% as.integer(2:4), nin = tiT %nin% as.integer(2:4)),
    tit_mat_cpp = function() {
      list(plain = pfns(".mat_this_in_vec_that")(matrix(tiT, 100, 50), as.integer(c(2, 5, 9))),
           inv = pfns(".mat_this_in_vec_that")(matrix(tiT, 100, 50), as.integer(c(2, 5, 9)),
                                               invert = TRUE),
           nt4 = pfns(".mat_this_in_vec_that")(matrix(tiT, 100, 50), as.integer(c(2, 5, 9)),
                                               nthreads = 4L))
    },
    tit_sweep = function() {
      # Note: a `that` vector containing negative values is deliberately absent.
      # vec/mat_this_in_vec_that() size its lookup table as max(that) + 1 and
      # then writes at each value of that, so a negative value writes out of
      # bounds and yields heap-dependent (nondeterministic) results. That is a
      # pre-existing defect in both v1.34.5 and the current tree; see
      # test-07-this-in-that.R for the deterministic subset of this behaviour.
      variants <- list(integer(0), 0:3, 5L, 1:30, c(0L, 12L), c(12L, 0L),
                       NA_integer_, c(3L, NA_integer_))
      out <- lapply(seq_along(variants), function(i) {
        tryCatch(this_in_that(tiT, variants[[i]]), error = function(e) pf_error_val(e))
      })
      names(out) <- NULL
      out$mat <- this_in_that(matrix(tiT, 100, 50), as.integer(c(0, 5, 12)))
      out$which <- this_in_that(tiT, as.integer(c(0, 5, 12)), type = "which")
      out$dup <- this_in_that(as.integer(c(1, 2, 2, 3, 3, 3)), as.integer(c(2, 3)))
      out$single <- this_in_that(5L, 5L)
      out$mat_cpp <- pfns(".mat_this_in_vec_that")(matrix(tiT, 100, 50), as.integer(c(0, 5, 12)))
      out$vec_cpp <- pfns(".vec_this_in_vec_that")(tiT, as.integer(c(0, 5, 12)))
      out
    },
    tit_vs_base = function() {
      set.seed(pf_seed + 61)
      a <- as.integer(sample(0:30, 4000, replace = TRUE))
      b <- as.integer(sort(unique(a[sample(length(a), 8)])))
      list(
        logical = identical(this_in_that(a, b), a %in% b),
        inverted = identical(this_in_that(a, b, invert = TRUE), !(a %in% b)),
        mat = identical(this_in_that(matrix(a, 100, 40), b), matrix(a %in% b, 100, 40))
      )
    },
    tit_nthreads = function() {
      m <- matrix(as.integer(sample(0:20, 2500, replace = TRUE)), 50, 50)
      lapply(c(1, 2, 4, 8), function(t) pfns(".mat_this_in_vec_that")(m, as.integer(c(3, 7, 19)),
                                                                      nthreads = t))
    },

    # ---- segim sub-routines ----
    segimstats_basic = function() {
      profoundSegimStats(image = img, segim = segim, sky = 0, skyRMS = 1, sortcol = "segID",
                         decreasing = FALSE, rotstats = FALSE, boundstats = FALSE)
    },
    segimstats_rich = function() {
      profoundSegimStats(image = img, segim = segim, sky = 0, skyRMS = 1, sortcol = "flux",
                         decreasing = TRUE, rotstats = TRUE, boundstats = TRUE, offset = 2)
    },
    segimedge = function() profoundSegimEdge(segim = segim),
    segimnear = function() profoundSegimNear(segim = segim, offset = 1),
    segimgroup = function() profoundSegimGroup(segim = segim),
    segimcompare = function() {
      profoundSegimCompare(segim_1 = segim,
                           segim_2 = dilate(segim, pf_brush(3, "disc")), threshold = 0.5)
    },
    segimmerge = function() {
      add <- matrix(0L, n, n)
      add[c(3, 12, 30, 45), c(3, 20, 33, 50)] <- 1:4
      profoundSegimMerge(image = img, segim_base = segim, segim_add = add, sky = 5)
    },
    segimkeep = function() {
      gg <- profoundSegimGroup(segim = segim)
      profoundSegimKeep(segim = segim, groupim = gg$groupim,
                        groupID_merge = gg$groupsegID$groupID[1:2])
    },
    mergsegid = function() profoundMergeSegID(segID_merge = segids),
    zapsegid = function() profoundZapSegID(segID = segids, segID_merge = segids[1:2]),
    shareflux = function() {
      ss <- profoundSegimStats(image = img, segim = segim, sky = 0, skyRMS = 1)
      idsel <- ss$segID
      k <- length(idsel)
      sm <- matrix(1 / k, k, k)
      colnames(sm) <- as.character(idsel)
      list(
        uniform = profoundShareFlux(segstats = ss, sharemat = sm),
        weighted = profoundShareFlux(segstats = ss, sharemat = sm, weights = seq_len(k))
      )
    },
    dilate_public = function() {
      list(all = profoundDilate(segim = segim, size = 5, shape = "disc", expand = "all", iters = 1),
           subset = profoundDilate(segim = segim, size = 3, shape = "box",
                                   expand = as.integer(segids[c(1, 3)]), iters = 1),
           iters3 = profoundDilate(segim = segim, size = 3, shape = "box", expand = "all", iters = 3))
    },
    makesegim = function() {
      r <- profoundMakeSegim(image = img, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                             sigma = 1, smooth = TRUE, plot = FALSE, stats = TRUE,
                             rotstats = FALSE, verbose = FALSE, watershed = "ProFound",
                             nthreads = 1L)
      list(segim = r$segim, segstats = r$segstats)
    },
    makesegim_nosmooth = function() {
      r <- profoundMakeSegim(image = img, skycut = 2, pixcut = 5, tolerance = 6, ext = 3,
                             smooth = FALSE, plot = FALSE, stats = TRUE, rotstats = TRUE,
                             verbose = FALSE, nthreads = 1L)
      list(segim = r$segim, segstats = r$segstats)
    },
    makesegim_old = function() {
      r <- profoundMakeSegim(image = img, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                             sigma = 1, smooth = TRUE, plot = FALSE, stats = TRUE,
                             verbose = FALSE, watershed = "ProFound-old", nthreads = 1L)
      list(segim = r$segim, segstats = r$segstats)
    },
    makesegimdilate = function() {
      suppressWarnings(profoundMakeSegimDilate(image = img, segim = segim, size = 5,
                                                shape = "disc", expand = "all", sky = 0,
                                                skyRMS = 1, verbose = FALSE, plot = FALSE,
                                                stats = TRUE, rotstats = TRUE))
    },
    makesegimexpand = function() {
      suppressWarnings(profoundMakeSegimExpand(image = img, segim = segim, skycut = 1.5,
                                                expandsigma = 5, expand = "all", sky = 0,
                                                skyRMS = 1, verbose = FALSE, plot = FALSE,
                                                stats = TRUE))
    },

    # ---- sky ----
    sky_est = function() profoundSkyEst(image = img, objects = matrix(as.integer(segim > 0), n, n)),
    sky_est_loc = function() {
      profoundSkyEstLoc(image = img, objects = matrix(as.integer(segim > 0), n, n),
                        loc = c(32, 32))
    },
    sky_grid_new = function() {
      profoundMakeSkyGrid(image = img, objects = matrix(as.integer(segim > 0), n, n), sky = 0,
                          box = c(20, 20), grid = c(20, 20), skygrid_type = "new",
                          type = "bicubic")
    },
    sky_grid_old = function() {
      profoundMakeSkyGrid(image = img, objects = matrix(as.integer(segim > 0), n, n), sky = 0,
                          box = c(20, 20), grid = c(20, 20), skygrid_type = "old",
                          type = "bicubic")
    },
    sky_grid_bilinear = function() {
      profoundMakeSkyGrid(image = img, objects = matrix(as.integer(segim > 0), n, n), sky = 0,
                          box = c(20, 20), grid = c(20, 20), skygrid_type = "new",
                          type = "bilinear")
    },
    sky_types = function() {
      obj <- matrix(as.integer(segim > 0), n, n)
      out <- list()
      for (st in c("mean", "median", "mode", "quanlo", "quanhi", "quan50")) {
        for (srt in c("quanlo", "mean")) {
          out[[paste0(st, "-", srt)]] <- tryCatch({
            r <- profoundMakeSkyGrid(image = img, objects = obj, sky = 0, box = c(20, 20),
                                     grid = c(20, 20), skytype = st, skyRMStype = srt)
            list(sky = r$sky, skyRMS = r$skyRMS)
          }, error = function(e) pf_error_val(e))
        }
      }
      out
    },
    sky_plane = function() profoundSkyPlane(image = img, objects = matrix(as.integer(segim > 0), n, n)),
    sky_poly = function() profoundSkyPoly(image = img, objects = matrix(as.integer(segim > 0), n, n),
                                          degree = 2),
    sky_chisel = function() profoundChisel(image = img, sky = 0, skythresh = 0.005),
    sky_scan = function() profoundSkyScan(image = img, mask = NULL, clip = c(0, 1)),

    # ---- aperphot ----
    aperphot_tar = function() {
      tar <- data.frame(segID = segids[seq_len(min(8, length(segids)))])
      profoundAperPhot(image = img, segim = segim, app_diam = 5, tar = tar, pixscale = 1,
                       magzero = 0, correction = TRUE, centype = "mean", depth = 4,
                       verbose = FALSE)
    },
    aperphot_all = function() {
      profoundAperPhot(image = img, segim = segim, app_diam = 3, pixscale = 1,
                       correction = FALSE, depth = 4, verbose = FALSE)
    },
    aperphot_depths = function() {
      lapply(c(0, 2, 4, 6), function(d) {
        suppressWarnings(profoundAperPhot(image = img, segim = segim, app_diam = 4,
                                          pixscale = 1, depth = d, verbose = FALSE))
      })
    },
    aperphot_mask = function() {
      mk <- matrix(0L, n, n)
      mk[1:10, 1:10] <- 1L
      profoundAperPhot(image = img, segim = segim, app_diam = 5, mask = mk, pixscale = 1,
                       depth = 4, verbose = FALSE)
    },
    aperphot_coord = function() {
      # profoundAperPhot() adds 0.5 to supplied coordinates then takes ceiling(),
      # so pixel (i, j) is addressed by xcen = i - 1, ycen = j - 1.
      sel <- segids[seq_len(min(6, length(segids)))]
      loc <- t(vapply(sel, function(id) which(segim == id, arr.ind = TRUE)[1, ], integer(2)))
      tar <- data.frame(xcen = loc[, 1] - 1, ycen = loc[, 2] - 1)
      profoundAperPhot(image = img, segim = segim, app_diam = 5, tar = tar, pixscale = 1,
                       depth = 4, verbose = FALSE)
    },

    # ---- utility ----
    flux_magenta = function() {
      list(
        f2m = profoundFlux2Mag(flux = c(1, 10, 100, 1e3), magzero = 25),
        m2f = profoundMag2Flux(mag = c(20, 22.5, 25, 27.5), magzero = 25),
        f2sb = profoundFlux2SB(flux = c(1, 10, 100), magzero = 25, pixscale = 0.4),
        sb2f = profoundSB2Flux(SB = c(26, 27, 28), magzero = 25, pixscale = 0.4),
        mu2m = profoundMu2Mag(mu = c(26, 27, 28), re = 2, axrat = 0.7, pixscale = 0.4),
        m2mu = profoundMag2Mu(mag = 25, re = 2, axrat = 0.7, pixscale = 0.4),
        gain = profoundGainConvert(gain = 10, magzero = 25, magzero_new = 30)
      )
    },
    ellipse_seg = function() {
      list(
        disc = profoundEllipseSeg(xcen = 30, ycen = 30, rad = 12, ang = 0, axrat = 1, dim = c(60, 60)),
        ell = profoundEllipseSeg(xcen = 30, ycen = 30, rad = 15, ang = 30, axrat = 0.5,
                                 dim = c(60, 60)),
        box = profoundEllipseSeg(xcen = 30, ycen = 30, rad = 10, ang = 0, axrat = 1, dim = c(60, 60))
      )
    },
    brush_shapes = function() {
      out <- list()
      for (shape in c("box", "disc", "diamond", "Gaussian", "line")) {
        for (size in c(1, 3, 5, 7, 9)) {
          out[[paste0(shape, size)]] <- pf_brush(size, shape)
        }
      }
      out
    },
    covmat_basic = function() {
      profoundCovMat(image = img, objects = matrix(as.integer(segim > 0), n, n))
    },

    # ---- the main entry point ----
    ppf_basic = function() {
      r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                            sigma = 1, smooth = TRUE, size = 5, shape = "disc", iters = 6,
                            threshold = 1.05, magzero = 0, pixscale = 1, redosky = TRUE,
                            redoskysize = 21, box = c(20, 20), grid = c(20, 20),
                            skygrid_type = "new", type = "bicubic", skytype = "median",
                            skyRMStype = "quanlo", sigmasel = 1, skypixmin = 200,
                            boxadd = c(10, 10), conviters = 100, deblend = FALSE,
                            doclip = TRUE, verbose = FALSE, plot = FALSE, stats = TRUE,
                            rotstats = FALSE, boundstats = FALSE, sortcol = "segID",
                            decreasing = FALSE, pixelcov = FALSE, convtype = "brute",
                            convmode = "extended", fluxtype = "Raw", nthreads = 1L)
      pf_core_result(r)
    },
    ppf_dilate3 = function() {
      r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                            size = 3, shape = "box", iters = 3, threshold = 1.02,
                            box = c(16, 16), grid = c(16, 16), plot = FALSE, stats = TRUE,
                            verbose = FALSE, rotstats = TRUE, boundstats = TRUE, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_shapes = function() {
      out <- list()
      for (sh in c("box", "disc", "diamond")) {
        r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 5, size = 7,
                              shape = sh, box = c(20, 20), grid = c(20, 20), plot = FALSE,
                              stats = TRUE, verbose = FALSE, nthreads = 1L)
        out[[sh]] <- pf_core_result(r)
      }
      out
    },
    ppf_no_dilate = function() {
      r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 0, size = 5,
                            box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                            verbose = FALSE, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_iters_sweep = function() {
      lapply(c(0, 1, 2, 3, 4, 6, 10), function(it) {
        r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = it, size = 5,
                              threshold = 1.05, box = c(20, 20), grid = c(20, 20),
                              plot = FALSE, stats = TRUE, verbose = FALSE, nthreads = 1L)
        pf_core_result(r)
      })
    },
    ppf_thresholds = function() {
      lapply(c(1.01, 1.05, 1.2, 2), function(th) {
        r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 6, size = 5,
                              threshold = th, box = c(20, 20), grid = c(20, 20), plot = FALSE,
                              stats = TRUE, verbose = FALSE, nthreads = 1L)
        pf_core_result(r)
      })
    },
    ppf_sky_old = function() {
      r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 2, size = 3,
                            box = c(20, 20), grid = c(20, 20), skygrid_type = "old",
                            type = "bicubic", plot = FALSE, stats = TRUE, verbose = FALSE,
                            nthreads = 1L)
      pf_core_result(r)
    },
    ppf_pixelcov = function() {
      r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 2, size = 5,
                            box = c(20, 20), grid = c(20, 20), pixelcov = TRUE, plot = FALSE,
                            stats = TRUE, verbose = FALSE, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_deblend = function() {
      im <- twin_3_4 + 5
      r <- profoundProFound(image = im, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                            iters = 3, size = 5, box = c(10, 10), grid = c(10, 10),
                            redosky = FALSE, sky = 5, skyRMS = 1, deblend = TRUE, df = 3,
                            radtrunc = 2, iterative = TRUE, plot = FALSE, stats = TRUE,
                            verbose = FALSE, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_water_old = function() {
      r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                            iters = 3, size = 5, box = c(20, 20), grid = c(20, 20),
                            watershed = "ProFound-old", plot = FALSE, stats = TRUE,
                            verbose = FALSE, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_relclip = function() {
      r <- profoundProFound(image = img, skycut = 1, pixcut = 2, tolerance = 3, reltol = 0.7,
                            cliptol = 60, ext = 3, iters = 4, size = 5, box = c(20, 20),
                            grid = c(20, 20), plot = FALSE, stats = TRUE, verbose = FALSE,
                            nthreads = 1L)
      pf_core_result(r)
    },
    ppf_segim_given = function() {
      r <- profoundProFound(image = img, segim = segim, redosegim = FALSE, skycut = 1.5,
                            pixcut = 3, iters = 4, size = 5, box = c(20, 20), grid = c(20, 20),
                            plot = FALSE, stats = TRUE, verbose = FALSE, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_mask_na = function() {
      mk <- matrix(0L, n, n)
      mk[1:10, 1:10] <- 1L
      mk[50:60, 30:40] <- 1L
      im <- img
      im[15, 15] <- NA
      r <- profoundProFound(image = im, mask = mk, skycut = 1.5, pixcut = 3, iters = 4,
                            size = 5, box = c(20, 20), grid = c(20, 20), plot = FALSE,
                            stats = TRUE, verbose = FALSE, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_lowmemory = function() {
      r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 4, size = 5,
                            box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                            verbose = FALSE, lowmemory = TRUE, keepim = FALSE, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_jansky = function() {
      r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 2, size = 5,
                            box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                            verbose = FALSE, fluxtype = "Jansky", magzero = 30, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_stats_rich = function() {
      r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 3, size = 5,
                            box = c(20, 20), grid = c(20, 20), plot = FALSE, stats = TRUE,
                            verbose = FALSE, boundstats = TRUE, nearstats = TRUE,
                            groupstats = TRUE, haralickstats = TRUE, groupby = "segim_orig",
                            offset = 2, nthreads = 1L)
      pf_core_result(r, extra = TRUE)
    },
    ppf_nthreads = function() {
      lapply(c(1, 2, 4, 8), function(t) {
        r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                              iters = 4, size = 5, box = c(20, 20), grid = c(20, 20),
                              plot = FALSE, stats = TRUE, verbose = FALSE, nthreads = t)
        pf_core_result(r)
      })
    },
    ppf_grid_types = function() {
      out <- list()
      for (ty in c("bicubic", "bilinear")) {
        for (sgt in c("new", "old")) {
          key <- paste0(sgt, "-", ty)
          out[[key]] <- tryCatch({
            r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 2, size = 5,
                                  box = c(20, 20), grid = c(20, 20), skygrid_type = sgt,
                                  type = ty, plot = FALSE, stats = TRUE, verbose = FALSE,
                                  nthreads = 1L)
            pf_core_result(r)
          }, error = function(e) pf_error_val(e))
        }
      }
      out
    },
    ppf_fluxtypes = function() {
      out <- list()
      for (ft in c("Raw", "Jansky", "Microjansky")) {
        out[[ft]] <- tryCatch({
          r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 1, size = 5,
                                box = c(20, 20), grid = c(20, 20), fluxtype = ft, magzero = 27,
                                plot = FALSE, stats = TRUE, verbose = FALSE, nthreads = 1L)
          pf_core_result(r)
        }, error = function(e) pf_error_val(e))
      }
      out[["bogus"]] <- tryCatch(profoundProFound(image = img, fluxtype = "notreal"),
                                 error = function(e) pf_error_val(e))
      out
    },
    ppf_sbdilate = function() {
      out <- list()
      for (sbd in c(1, 2, 5)) {
        r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 6, size = 5,
                              SBdilate = sbd, SBN100 = 100, box = c(20, 20), grid = c(20, 20),
                              plot = FALSE, stats = TRUE, verbose = FALSE, nthreads = 1L)
        out[[as.character(sbd)]] <- pf_core_result(r)
      }
      out
    },
    ppf_aperphot = function() {
      r <- profoundProFound(image = img, skycut = 1.5, pixcut = 3, iters = 2, size = 5,
                            box = c(20, 20), grid = c(20, 20), app_diam = 5, plot = FALSE,
                            stats = TRUE, verbose = FALSE, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_viking = function() {
      f <- system.file("extdata", "VIKING", "mystery_VIKING_Z.fits", package = "ProFound")
      if (!file.exists(f)) {
        return("MISSING")
      }
      im <- Rfits::Rfits_read_image(f)
      r <- profoundProFound(image = im, skycut = 1.5, pixcut = 3, tolerance = 4, ext = 2,
                            sigma = 1, smooth = TRUE, size = 5, shape = "disc", iters = 6,
                            threshold = 1.05, box = c(64, 64), grid = c(64, 64), plot = FALSE,
                            stats = TRUE, rotstats = FALSE, boundstats = FALSE, verbose = FALSE,
                            nthreads = 1L, magzero = 31)
      pf_core_result(r)
    },
    ppf_ir_psf = function() {
      fim <- system.file("extdata", "IRdata", "s250_im.fits", package = "ProFound")
      fpsf <- system.file("extdata", "IRdata", "s250_psf.fits", package = "ProFound")
      if (!file.exists(fim) || !file.exists(fpsf)) {
        return("MISSING")
      }
      im <- Rfits::Rfits_read_image(fim)
      psf <- Rfits::Rfits_read_image(fpsf)
      r <- profoundProFound(image = im, psf = psf, skycut = 1.5, pixcut = 3, tolerance = 4,
                            ext = 2, size = 5, iters = 3, box = c(20, 20), grid = c(20, 20),
                            plot = FALSE, stats = TRUE, verbose = FALSE, app_diam = 4,
                            nthreads = 1L)
      pf_core_result(r)
    },
    ppf_all_zero = function() {
      r <- profoundProFound(image = matrix(0, 32, 32), skycut = 1.5, pixcut = 3, iters = 2,
                            size = 5, box = c(8, 8), grid = c(8, 8), plot = FALSE,
                            stats = TRUE, verbose = FALSE, nthreads = 1L)
      pf_core_result(r)
    },
    ppf_errors = function() {
      list(
        no_image = tryCatch(profoundProFound(), error = function(e) pf_error_val(e)),
        bad_fluxtype = tryCatch(profoundProFound(image = img, fluxtype = "x"),
                                error = function(e) pf_error_val(e)),
        bad_watershed = tryCatch(profoundProFound(image = img, watershed = "x"),
                                 error = function(e) pf_error_val(e)),
        bad_skygrid = tryCatch(profoundProFound(image = img, skygrid_type = "x"),
                               error = function(e) pf_error_val(e)),
        matrix_ok = tryCatch({
          r <- profoundProFound(image = img, iters = 0, box = c(20, 20), grid = c(20, 20),
                                plot = FALSE, stats = FALSE, verbose = FALSE)
          "ok"
        }, error = function(e) pf_error_val(e))
      )
    }
  )
}

# profoundProFound returns wall-clock fields and an environment-heavy call; keep
# only the scientifically meaningful components for the goldens.
pf_core_result <- function(r, extra = FALSE) {
  keep <- c("segim", "segim_orig", "segstats", "Nseg", "sky", "skyRMS", "objects",
             "imarea", "pixscale", "magzero", "gain", "SBlim")
  if (extra) {
    keep <- c(keep, "groupstats", "group", "near", "haralick", "objects_redo", "skyLL",
              "skyChiSq", "skyChiSqMap")
  }
  out <- keep[keep %in% names(r)]
  res <- r[out]
  if ("call" %in% names(r)) {
    res$call_name <- as.character(r$call[[1]])
  }
  res
}

