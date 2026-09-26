# Helpers shared by every ProFound test file.
#
# The suite has two layers:
#   1. Focused unit tests that assert mathematical properties and invariants.
#   2. A golden-value regression layer (see helper-cases.R and
#      golden/goldens.csv) whose reference values were produced by the
#      committed ProFound v1.34.5 (2026-09-24) source, so it pins the speed
#      focused C++ edits in src/ to the pre-edit numerical behaviour.

pf_seed <- 4242L

# Internal (non-exported) routines are addressed explicitly so tests keep
# working whether or not the package is attached with its internals visible.
pfns <- function(name) get(name, envir = asNamespace("ProFound"), inherits = FALSE)

# Deterministic, self-contained synthetic inputs -----------------------------

# A Gaussian-source image plus noise, returned in sky/RMS units so that
# watershed segmentation genuinely produces several segments and background.
pf_make_img <- function(n = 64L, seed = pf_seed) {
  set.seed(seed)
  grid <- expand.grid(x = seq_len(n), y = seq_len(n))
  amp <- runif(9, 30, 200)
  xc <- runif(9, 0.2, 0.8) * n
  yc <- runif(9, 0.2, 0.8) * n
  w <- runif(9, 2, 8)
  z <- rep(0, nrow(grid))
  for (k in seq_along(amp)) {
    z <- z + amp[k] * exp(-((grid$x - xc[k])^2 + (grid$y - yc[k])^2) / (2 * w[k]^2))
  }
  raw <- matrix(z, n, n) + matrix(rnorm(n * n, sd = 1), n, n) + 5
  no_obj <- matrix(0L, n, n)
  sky <- profoundSkyEst(image = raw, objects = no_obj)$sky
  skyRMS <- profoundSkyEst(image = profoundImDiff(raw - sky, 3), objects = no_obj)$skyRMS
  out <- (raw - sky) / skyRMS
  out[!is.finite(out)] <- 0
  out
}

# Two overlapping Gaussians; the classic de-blending stress case.
pf_make_twin <- function(sd = 3, gap = 4, n = 40L) {
  grid <- expand.grid(x = seq_len(n), y = seq_len(n))
  z <- 100 * exp(-((grid$x - (n / 2 - gap / 2))^2 + (grid$y - n / 2)^2) / (2 * sd^2)) +
    120 * exp(-((grid$x - (n / 2 + gap / 2))^2 + (grid$y - n / 2)^2) / (2 * sd^2))
  matrix(z, n, n)
}

# Sparse label image with segments away from and right at the borders.
pf_make_sparse <- function(n = 21L) {
  out <- matrix(0L, n, n)
  out[c(2, 5, 10, 15, 20), c(2, 5, 11, 20)] <- 1:4
  out[1, 1] <- 5L
  out[n, n] <- 6L
  out
}

# Scattered label image large enough to exercise the interior fast path of
# dilate alongside its boundary path.
pf_make_big_sparse <- function(n = 60L, nseg = 12L, seed = pf_seed + 21L) {
  set.seed(seed)
  out <- matrix(0L, n, n)
  loc <- cbind(sample.int(n, 30), sample.int(n, 30))
  out[loc] <- sample.int(nseg, 30, replace = TRUE)
  out
}

pf_brush <- function(size, shape) pfns(".makeBrush")(size, shape)

# Canonicalisation and digests -----------------------------------------------

# Watershed segment ids are an artefact of pixel processing order, and the C++
# sort that produces that order is not stable for exactly-equal pixel values.
# Integer matrices are therefore relabelled into 1..k by first appearance
# (column-major, i.e. R's storage order). This pins segmentation geometry
# without pinning an arbitrary label.
pf_canon_labels <- function(x) {
  vals <- as.vector(x)
  ids <- unique(vals[vals > 0])
  if (!length(ids)) {
    return(x)
  }
  map <- seq_along(ids)
  names(map) <- as.character(ids)
  fresh <- vals
  hit <- vals > 0
  fresh[hit] <- as.integer(map[as.character(vals[hit])])
  matrix(fresh, nrow(x), ncol(x))
}

# Round to 12 significant digits so last-bit platform noise in pow/exp/sqrt
# cannot break a golden, while any real difference (>= ~1e-11 relative) still
# does. NaN and NA are kept distinct.
pf_sig <- function(x) {
  if (is.integer(x)) {
    return(x)
  }
  if (!is.numeric(x)) {
    return(as.character(x))
  }
  out <- signif(as.numeric(x), digits = 12)
  out[is.nan(x)] <- NaN
  out[is.na(x) & !is.nan(x)] <- NA_real_
  out
}

# serialize() stores environments by memory address, so any object holding one
# (an lm fit, a closure's body, ...) would produce a digest that changes between
# sessions. Such pieces are reduced to their address-free content.
pf_canon <- function(x) {
  if (is.function(x)) {
    return("function")
  }
  if (is.environment(x)) {
    return(NULL)
  }
  if (inherits(x, "lm")) {
    return(list(coef = x$coef, sigma = x$sigma, df_residual = x$df.residual,
                rank = x$rank, residuals = as.vector(x$residuals),
                fitted = as.vector(x$fitted.values), weights = x$weights))
  }
  if (inherits(x, "data.frame")) {
    if (!ncol(x) || !nrow(x)) {
      return(x)
    }
    cols <- lapply(x, pf_sig)
    nms <- names(cols)
    # segID is arbitrary (see pf_canon_labels), so rows are ordered by the
    # remaining physical columns and segID itself is replaced by a rank.
    keycols <- setdiff(nms, "segID")
    if (length(keycols) && length(cols)) {
      # unname() because some ProFound data frames have a column named "sep"
      keylist <- unname(lapply(cols[keycols], format))
      keylist <- c(keylist, list(sep = "|"))
      keys <- do.call(paste, keylist)
      ord <- order(keys, seq_along(keys), method = "radix")
      cols <- lapply(cols, function(v) v[ord])
    }
    if ("segID" %in% nms) {
      cols[["segID"]] <- seq_along(cols[[1]])
    }
    attr(cols, "class") <- "data.frame"
    attr(cols, "row.names") <- seq_along(cols[[1]])
    attr(cols, "names") <- nms
    return(cols)
  }
  if (is.integer(x) && is.matrix(x)) {
    return(pf_canon_labels(x))
  }
  if (is.numeric(x)) {
    return(pf_sig(x))
  }
  if (is.logical(x) || is.character(x) || is.raw(x)) {
    return(x)
  }
  if (is.list(x)) {
    out <- lapply(x, pf_canon)
    names(out) <- names(x)
    return(out)
  }
  x
}

pf_digest <- function(x) {
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeBin(serialize(pf_canon(x), connection = NULL, version = 2), tmp)
  unname(tools::md5sum(tmp))
}

# A case that fails is recorded as its message so the golden pins the failure
# behaviour too, rather than aborting the whole run.
pf_error_val <- function(e) paste("ERROR:", conditionMessage(e))

# Locate the testthat directory whether running under R CMD check, test_dir(),
# or interactively with the working directory already at tests/testthat.
pf_test_dir <- function() {
  cand <- character(0)
  tp <- tryCatch(testthat::test_path(), error = function(e) NULL)
  if (!is.null(tp)) {
    cand <- c(cand, tp)
  }
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  cand <- c(cand, cwd, file.path(cwd, "tests", "testthat"),
            file.path(dirname(cwd), "testthat"))
  hit <- cand[file.exists(file.path(cand, "helper-cases.R"))]
  if (!length(hit)) {
    stop("could not locate the tests/testthat directory")
  }
  hit[[1]]
}

golden_dir <- function() file.path(pf_test_dir(), "golden")

pf_golden_table <- function() {
  path <- file.path(golden_dir(), "goldens.csv")
  if (!file.exists(path)) {
    return(NULL)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

pf_fixture <- function(name) {
  path <- file.path(golden_dir(), "fixtures", name)
  if (!file.exists(path)) {
    return(NULL)
  }
  readRDS(path)
}

# The fixture .rds files hold full-precision outputs captured from the
# committed v1.34.5 build. Their checksums are pinned here so a fixture can
# never be silently swapped for one produced by a different build.
pf_fixture_manifest <- function() {
  c(
    "akima_v1345.rds"      = "6a4de4201981af9d7eef05ff537d2cd2",
    "cover_v1345.rds"      = "0c426c972c40d8edf56790357dde940a",
    "dilate_v1345.rds"      = "4b76bc09c6903b3e68629795a30c2d15",
    "thisinthat_v1345.rds"      = "1aeb1c76a346d1e22c6faafd24ad499b",
    "watershed_v1345.rds"      = "5f809772412276586faa035421db1032",
    "weights_v1345.rds"      = "02839fbfa5a93efbbb09f0d0669c72da"
  )
}

# Returns the fixture, or NULL (with the caller expected to skip) when missing
# or not matching the pinned digest.
pf_fixture_checked <- function(name) {
  path <- file.path(golden_dir(), "fixtures", name)
  if (!file.exists(path)) {
    return(NULL)
  }
  man <- pf_fixture_manifest()
  if (!isTRUE(unname(tools::md5sum(path)) == man[[name]])) {
    return(NULL)
  }
  readRDS(path)
}
