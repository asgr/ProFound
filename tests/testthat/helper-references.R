# Reference implementations used as oracles by the unit tests. These are
# deliberately direct transcriptions of the ORIGINAL (pre speed-rewrite) C++
# logic, so the rewritten src/ code is checked against the algorithm it is
# supposed to implement rather than against itself.

# Faithful R transcription of dilate_cpp() before the two-phase scan/scatter
# rewrite: for every labelled pixel, scatter its id into the empty cells of the
# kernel footprint, keeping the smallest id that reaches each cell.
dilate_reference <- function(segim, kern, expand = 0L) {
  # .dilate_cpp() takes an IntegerMatrix, so Rcpp truncates a numeric kernel
  # with as.integer() before the C++ loop ever sees it. Mirror that here (it is
  # what makes a Gaussian brush, whose values are all < 1, a no-op).
  kern <- matrix(as.integer(kern), nrow(kern), ncol(kern))
  srow <- nrow(segim)
  scol <- ncol(segim)
  krow <- nrow(kern)
  kcol <- ncol(kern)
  ro <- (krow - 1) %/% 2
  co <- (kcol - 1) %/% 2
  max_segim <- max(segim)

  seglogic <- rep(FALSE, max(max_segim, 0L))
  if (length(expand) > 0 && expand[1] > 0) {
    for (k in seq_along(expand)) {
      if (expand[k] <= max_segim) {
        seglogic[expand[k]] <- TRUE
      }
    }
  }

  segim_new <- matrix(0L, srow, scol)
  for (i in seq_len(srow)) {
    for (j in seq_len(scol)) {
      if (segim[i, j] > 0) {
        checkseg <- TRUE
        if (expand[1] > 0) {
          checkseg <- seglogic[segim[i, j]]
        }
        if (checkseg) {
          for (m in seq(krow)) {
            for (n in seq(kcol)) {
              if (kern[m, n] > 0) {
                if (m - 1 == ro && n - 1 == co) {
                  segim_new[i, j] <- segim[i, j]
                } else {
                  x <- i + m - 1 - ro
                  y <- j + n - 1 - co
                  if (x >= 1 && x <= srow && y >= 1 && y <= scol && segim[x, y] == 0) {
                    if (segim[i, j] < segim_new[x, y] || segim_new[x, y] == 0) {
                      segim_new[x, y] <- segim[i, j]
                    }
                  }
                }
              }
            }
          }
        } else {
          segim_new[i, j] <- segim[i, j]
        }
      }
    }
  }
  segim_new
}

# Watershed segmentation with the same signature style as water_cpp().
pf_water <- function(image, abstol = 4, reltol = 0, cliptol = Inf, ext = 2L,
                     skycut = 1.5, pixcut = 3L, nthreads = 1L) {
  pfns("water_cpp")(image = as.vector(image), nx = as.integer(nrow(image)),
                    ny = as.integer(ncol(image)), abstol = abstol, reltol = reltol,
                    cliptol = cliptol, ext = as.integer(ext), skycut = skycut,
                    pixcut = as.integer(pixcut), verbose = FALSE, Ncheck = 1e6,
                    nthreads = as.integer(nthreads))
}

pf_dilate <- function(segim, kern, expand = 0L, nthreads = 1L) {
  pfns(".dilate_cpp")(as.matrix(segim), as.matrix(kern), as.integer(expand),
                     as.integer(nthreads))
}

# Relabel by first appearance so segment id numbering cannot affect an
# assertion about segmentation geometry.
pf_labels <- function(x) pf_canon_labels(x)

# Random labelled images, for oracle comparison across many shapes.
pf_random_segim <- function(n, seed) {
  set.seed(seed)
  out <- matrix(0L, n, n)
  k <- sample.int(6, 1) + 1L
  for (i in seq_len(k)) {
    cx <- sample.int(n, 1)
    cy <- sample.int(n, 1)
    w <- sample.int(max(2L, n %/% 4), 1)
    grid <- expand.grid(x = seq_len(n), y = seq_len(n))
    hot <- (grid$x - cx)^2 + (grid$y - cy)^2 <= w^2
    out[as.matrix(grid[hot, ])][out[as.matrix(grid[hot, ])] == 0] <- i
  }
  out
}

# ---------------------------------------------------------------------------
# Coverage-recursion oracles for src/ellip_cover.cpp and src/aper_cover.cpp.
# The speed edits added a quadratic bounding-box shortcut to the ellipse
# recursion, so these transcribe the *unsimplified* recursion and the code is
# checked against the algorithm it is meant to implement.
# ---------------------------------------------------------------------------

# Direct transcription of pixelCoverEllip() with no shortcut.
ref_cover_ellip <- function(dx, dy, xt, yt, xyt, depth) {
  if (depth == 0) {
    return(if (xt * dx * dx + yt * dy * dy + xyt * dx * dy <= 1) 1 else 0)
  }
  q <- 0.5 / 2^depth
  (ref_cover_ellip(dx - q, dy - q, xt, yt, xyt, depth - 1) +
     ref_cover_ellip(dx - q, dy + q, xt, yt, xyt, depth - 1) +
     ref_cover_ellip(dx + q, dy - q, xt, yt, xyt, depth - 1) +
     ref_cover_ellip(dx + q, dy + q, xt, yt, xyt, depth - 1)) / 4
}

ref_ellip_terms <- function(rad, ang, axrat) {
  a <- ang * 3.141593 / 180
  ca <- cos(a)
  sa <- sin(a)
  semi_min <- rad * axrat
  imn <- 1 / (semi_min * semi_min)
  imj <- 1 / (rad * rad)
  list(xt = ca * ca * imn + sa * sa * imj,
       yt = sa * sa * imn + ca * ca * imj,
       xyt = 2 * sa * ca * (imn - imj),
       semi_min = semi_min)
}

# profoundEllipCover() applies the recursion only to pixels in the annulus
# between (semi_min - sqrt(0.5)) and (rad + sqrt(0.5)); reproduce that gating.
ref_ellip_cover <- function(x, y, cx, cy, rad, ang, axrat, depth) {
  d <- ref_ellip_terms(rad, ang, axrat)
  sm2 <- d$semi_min - 0.7071068
  if (sm2 < 0) sm2 <- 0
  rp <- rad + 0.7071068
  vapply(seq_along(x), function(i) {
    dx <- x[i] - cx
    dy <- y[i] - cy
    if (!(abs(dx) < rp && abs(dy) < rp)) return(0)
    d2 <- dx * dx + dy * dy
    if (d2 < sm2 * sm2) return(1)
    if (d2 < rp * rp) return(ref_cover_ellip(dx, dy, d$xt, d$yt, d$xyt, depth))
    0
  }, numeric(1))
}

# Direct transcription of pixelCoverAper() with no shortcut.
ref_cover_aper <- function(dx, dy, d2, rad_2, rmin_2, rmax_2, depth) {
  if (depth == 0) return(if (d2 <= rad_2) 1 else 0)
  if (d2 > rmax_2) return(0)
  if (d2 <= rmin_2) return(1)
  q <- 0.5 / 2^depth
  (ref_cover_aper(dx - q, dy - q, (dx - q)^2 + (dy - q)^2, rad_2, rmin_2, rmax_2, depth - 1) +
     ref_cover_aper(dx - q, dy + q, (dx - q)^2 + (dy + q)^2, rad_2, rmin_2, rmax_2, depth - 1) +
     ref_cover_aper(dx + q, dy - q, (dx + q)^2 + (dy - q)^2, rad_2, rmin_2, rmax_2, depth - 1) +
     ref_cover_aper(dx + q, dy + q, (dx + q)^2 + (dy + q)^2, rad_2, rmin_2, rmax_2, depth - 1)) / 4
}

ref_aper_cover <- function(x, y, cx, cy, rad, depth) {
  rmin <- rad - 0.7071068
  rmin_2 <- if (rmin > 0) rmin^2 else -1
  rp <- rad + 0.7071068
  vapply(seq_along(x), function(i) {
    dx <- x[i] - cx
    dy <- y[i] - cy
    if (!(abs(dx) < rp && abs(dy) < rp)) return(0)
    ref_cover_aper(dx, dy, dx * dx + dy * dy, rad * rad, rmin_2, rp * rp, depth)
  }, numeric(1))
}

# Akima grid interpolation wrappers (src/akima.cpp), used by
# test-03-interpolation.R.
pf_akima <- function(x, y, grid, nxout, nyout) {
  out <- matrix(NA_real_, nxout, nyout)
  pfns(".interpolateAkimaGrid")(x, y, grid, out)
  out
}

# Values of the interpolant at the node positions of a half-integer node grid.
pf_node_cells <- function(nx, ny) {
  list(qx = seq_len(nx), qy = seq_len(ny))   # x = -0.5 + i  =>  i = 1..nx
}
