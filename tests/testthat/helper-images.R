# Shared fixtures for the functional (non-golden) tests.
#
# These build small synthetic images and segmentations with *known* answers, so
# assertions can be written against the truth rather than against the code's own
# output. The golden layer uses helper-profound.R / helper-cases.R instead, since
# its inputs must be reproducible bit-for-bit rather than interpretable.
#
# Synthetic image with a flat sky, Gaussian noise and a handful of Gaussian
# "sources". Positions are given in ProFound image coordinates, i.e. the centre
# of the first pixel is (0.5, 0.5), so a source at (30, 50) peaks in the pixel
# with index (30, 50). Defaults give two well separated, clearly detected
# sources.
makeTestImage = function(dim = c(101, 101), sky = 100, skyRMS = 5, sources = NULL){
  if(is.null(sources)){
    sources = data.frame(
      xcen = c(30, 70),
      ycen = c(50, 40),
      flux = c(2000, 1200),
      sigma = c(2, 3)
    )
  }

  image = matrix(stats::rnorm(prod(dim), mean = sky, sd = skyRMS), dim[1], dim[2])
  # Pixel indices converted to ProFound coordinates
  grid = expand.grid(1:dim[1] - 0.5, 1:dim[2] - 0.5)

  for(i in seq_len(nrow(sources))){
    rad2 = (grid[,1] - sources$xcen[i])^2 + (grid[,2] - sources$ycen[i])^2
    image = image + matrix(
      sources$flux[i] / (2 * pi * sources$sigma[i]^2) *
        exp(-rad2 / (2 * sources$sigma[i]^2)),
      dim[1], dim[2]
    )
  }

  image
}

# Small integer segmentation map with a touching pair (1, 2) and two isolated
# segments (3, 4). Layout is documented in test-segim-utilities.R.
makeTestSegim = function(){
  segim = matrix(0L, 21, 21)
  segim[3:6, 3:6] = 1L
  segim[3:6, 7:10] = 2L
  segim[15:18, 15:18] = 3L
  segim[15:18, 5:8] = 4L
  segim
}

# Gaussian PSF normalised to unit sum, used by the PSF fitting tests.
makeTestPSF = function(size = 15, sigma = 1.5){
  psf = matrix(0, size, size)
  grid = expand.grid(1:size - (size + 1)/2, 1:size - (size + 1)/2)
  psf[] = exp(-(grid[,1]^2 + grid[,2]^2) / (2 * sigma^2))
  psf / sum(psf)
}

# profoundProFound is chatty (the sky estimate, and messages when box/grid are
# clipped for small images), so run it quietly in tests.
runProFound = function(...){
  suppressMessages(profoundProFound(...))
}

# Place point sources with known fluxes into an image using a PSF. Positions are
# in ProFound image coordinates (pixel centre at index - 0.5).
addTestSources = function(image, xcen, ycen, flux, psf){
  dimpsf = dim(psf)

  for(i in seq_along(xcen)){
    xstart = round(xcen[i] + 0.5) - floor(dimpsf[1]/2)
    ystart = round(ycen[i] + 0.5) - floor(dimpsf[2]/2)
    xpix = xstart:(xstart + dimpsf[1] - 1)
    ypix = ystart:(ystart + dimpsf[2] - 1)

    inx = xpix >= 1 & xpix <= dim(image)[1]
    iny = ypix >= 1 & ypix <= dim(image)[2]

    image[xpix[inx], ypix[iny]] = image[xpix[inx], ypix[iny]] +
      psf[inx, iny] * flux[i]
  }

  image
}

# A blended pair, run through ProFound with group statistics, giving a real
# segstats / groupstats / groupsegID triple for the catalogue merging tests.
makeGroupedProFound = function(){
  set.seed(2206)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(45, 55), ycen = c(50, 50), flux = c(4000, 3000), sigma = c(3, 3)))
  runProFound(image, skycut = 2, box = 50, grid = 50, tolerance = 1,
              groupstats = TRUE)
}
