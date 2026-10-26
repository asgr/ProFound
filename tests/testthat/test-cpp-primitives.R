test_that("water_cpp segments a single isolated hill", {
  grid = expand.grid(1:41 - 0.5, 1:41 - 0.5)
  image = matrix(exp(-((grid[,1] - 20)^2 + (grid[,2] - 20)^2) / 18), 41, 41)

  segim = water_cpp(image, nx = 41L, ny = 41L, abstol = 1, ext = 2L, skycut = 0.5, pixcut = 1L)

  expect_true(is.matrix(segim))
  expect_equal(sort(unique(as.integer(segim))), c(0L, 1L))
  expect_equal(sum(segim > 0), 32)
})

test_that("water_cpp separates two well separated hills", {
  image = makeTwoHills()

  segim = water_cpp(image, nx = 41L, ny = 41L, abstol = 1, ext = 2L, skycut = 0.5, pixcut = 1L)

  expect_equal(max(segim), 2)
  # tabulate drops the zeros, so this is the pixel count of each segment
  expect_equal(unname(tabulate(segim)), c(32L, 32L))
  # Each segment is centred on one of the two input hills
  rows = which(segim > 0, arr.ind = TRUE)
  expect_equal(
    as.numeric(tapply(rows[, 1], segim[rows], mean)),
    c(12, 30),
    tolerance = 1.5
  )
})

test_that("water_cpp pixcut drops islands smaller than the cut", {
  image = makeTwoHills()

  expect_equal(max(water_cpp(image, 41L, 41L, abstol = 1, ext = 2L, skycut = 0.5, pixcut = 1L)), 2)
  expect_equal(max(water_cpp(image, 41L, 41L, abstol = 1, ext = 2L, skycut = 0.5, pixcut = 100L)), 0)
})

test_that("water_cpp skycut controls what counts as sky", {
  image = makeTwoHills()

  # A skycut above the peak value leaves everything as sky
  expect_equal(max(water_cpp(image, 41L, 41L, abstol = 1, ext = 2L, skycut = 2, pixcut = 1L)), 0)
})

test_that("water_cpp returns all sky for a featureless image", {
  image = matrix(-1, 11, 11)
  segim = water_cpp(image, nx = 11L, ny = 11L, abstol = 1, ext = 1L, skycut = 0, pixcut = 1L)

  expect_equal(max(segim), 0)
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

