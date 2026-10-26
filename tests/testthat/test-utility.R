test_that("magnitude and flux conversions are inverses", {
  mag = c(10, 15, 22.5)
  magzero = 25

  expect_equal(profoundFlux2Mag(profoundMag2Flux(mag, magzero = magzero), magzero = magzero), mag)
  expect_equal(profoundMag2Flux(profoundFlux2Mag(c(1, 10, 1000), magzero = 0), magzero = 0), c(1, 10, 1000))
})

test_that("flux to magnitude handles zero, negative and missing flux", {
  expect_equal(profoundFlux2Mag(c(1, 0, NA, -1), magzero = 0), c(0, NA, NA, NA))

  # magzero shifts magnitudes but not their validity
  expect_equal(profoundFlux2Mag(10, magzero = 30), -2.5 * log10(10) + 30)
})

test_that("surface brightness conversions round trip", {
  flux = c(1, 100, 1e4)
  magzero = 26
  pixscale = 0.2

  expect_equal(
    profoundSB2Flux(profoundFlux2SB(flux, magzero = magzero, pixscale = pixscale),
                    magzero = magzero, pixscale = pixscale),
    flux
  )
})

test_that("surface brightness is flux spread over a pixel", {
  # A flux of 1 in a 1 arcsec pixel is the same surface brightness as 1 mag
  expect_equal(profoundFlux2SB(flux = 1, magzero = 0, pixscale = 1), 0)

  # Smaller pixels pile the same flux into less area, so surface brightness rises
  expect_lt(profoundFlux2SB(1, pixscale = 0.5), profoundFlux2SB(1, pixscale = 1))
  expect_equal(profoundFlux2SB(1, pixscale = 0.5), 5 * log10(0.5))
})

test_that("Mag2Mu / Mu2Mag round trip", {
  args = list(re = 2, axrat = 0.6, pixscale = 0.3)
  expect_equal(
    do.call(profoundMu2Mag, c(list(mu = do.call(profoundMag2Mu, c(list(mag = 18), args))), args)),
    18
  )
})

test_that("gain conversion follows the magnitude zero point ratio", {
  # Reducing magzero by 5 magnitudes means each ADU holds 100x more flux
  expect_equal(profoundGainConvert(gain = 2, magzero = 30, magzero_new = 25), 200)
  expect_equal(profoundGainConvert(gain = 1), 1)
})

test_that("this_in_that gives set membership with logical and which output", {
  this = c(1L, 5L, 7L, 9L)
  that = c(5L, 9L, 11L)

  expect_equal(this_in_that(this, that), c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(this_in_that(this, that, type = 'which'), c(2L, 4L))
  expect_equal(this_in_that(this, that, invert = TRUE), c(TRUE, FALSE, TRUE, FALSE))
})

test_that("%fin% and %nin% are complementary", {
  this = c(1L, 5L, 7L, 9L)
  that = c(5L, 9L, 11L)

  expect_equal(this %fin% that, this_in_that(this, that))
  expect_false(any((this %fin% that) & (this %nin% that)))
})

test_that("this_in_that preserves the shape of matrix input", {
  mat = matrix(c(1L, 2L, 3L, 4L), 2, 2)
  out = this_in_that(mat, c(2L, 4L))

  expect_true(is.matrix(out))
  expect_identical(dim(out), dim(mat))
  # Column-major fill: mat is 1, 2 in column 1 and 3, 4 in column 2
  expect_equal(out, matrix(c(FALSE, TRUE, FALSE, TRUE), 2, 2))
  expect_equal(this_in_that(mat, c(2L, 4L), type = 'which'), c(2L, 4L))
})

test_that("this_in_that rejects non-integer input", {
  expect_error(this_in_that(c(1, 2), c(1L, 2L)), "this must be integer")
  expect_error(this_in_that(c(1L, 2L), c(1, 2)), "that must be integer")
})

test_that("makeBrush produces brushes of the requested shape", {
  box = profoundMakeMask(size = 11, shape = 'box')
  disc = profoundMakeMask(size = 11, shape = 'disc')
  diamond = profoundMakeMask(size = 11, shape = 'diamond')
  gaussian = profoundMakeMask(size = 11, shape = 'gaussian')

  expect_identical(dim(box), c(11L, 11L))
  expect_equal(sum(box > 0), 121)
  expect_equal(sum(disc > 0), 97)
  expect_equal(sum(diamond > 0), 61)
  expect_equal(sum(gaussian), 1, tolerance = 1e-12)

  # Every brush is symmetric about its centre
  expect_equal(disc, t(disc))
  expect_equal(diamond, t(diamond))
})

test_that("makeBrush rounds even sizes up with a warning", {
  expect_warning(br <- profoundMakeMask(size = 10, shape = 'box'),
                 "rounded to the next odd number")
  expect_identical(dim(br), c(11L, 11L))
})

test_that("makeBrush line shape is a single pixel wide diagonal", {
  # profoundMakeMask does not expose the line angle, so reach for the brush itself
  line = ProFound:::.makeBrush(size = 11, shape = 'line', angle = 45)

  expect_equal(sum(line > 0), 11)
  expect_equal(sum(rowSums(line) > 0), 11)
})

test_that("makeBrush rejects sizes that are not positive numbers", {
  expect_error(profoundMakeMask(size = 'big'), "'size' must be an odd integer")
  expect_error(profoundMakeMask(size = 0), "'size' must be an odd integer")
})
