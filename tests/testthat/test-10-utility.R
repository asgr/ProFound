# Unit conversions. These are pure arithmetic, so they are the cheapest and most
# exact tests in the suite: every conversion is checked against its analytic
# inverse or its defining formula.


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

test_that("flux utility conversions round-trip", {
  fluxes <- 10^seq(0, 3)
  mz <- 25
  # Flux <-> magnitude.
  expect_equal(profoundFlux2Mag(flux = fluxes, magzero = mz), mz - 2.5 * log10(fluxes),
               tolerance = 1e-12)
  expect_equal(profoundMag2Flux(mag = c(18, 21.5, 25, 28.2), magzero = mz),
               10^(-0.4 * (c(18, 21.5, 25, 28.2) - mz)), tolerance = 1e-12)
  expect_equal(profoundMag2Flux(profoundFlux2Mag(fluxes, magzero = mz), magzero = mz),
               fluxes, tolerance = 1e-10)

  # Surface brightness conversions use the pixel area.
  sb <- profoundFlux2SB(flux = 100, magzero = 25, pixscale = 0.4)
  expect_equal(profoundSB2Flux(sb, magzero = 25, pixscale = 0.4), 100, tolerance = 1e-9)

  expect_equal(profoundGainConvert(gain = 1, magzero = 25, magzero_new = 25), 1)
  expect_equal(profoundGainConvert(gain = 10, magzero = 25, magzero_new = 27.5),
               10 * 10^(-0.4 * (27.5 - 25)), tolerance = 1e-12)

  # mu <-> mag round trip.
  m <- 24
  mu <- profoundMag2Mu(mag = m, re = 2, axrat = 0.7, pixscale = 0.4)
  expect_equal(profoundMu2Mag(mu = mu, re = 2, axrat = 0.7, pixscale = 0.4), m,
               tolerance = 1e-9)
})
