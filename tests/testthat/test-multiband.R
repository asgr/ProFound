test_that("profoundMultiBand runs a two band photometry workflow", {
  skip_if_not_installed('Rfits')
  skip_if_not_installed('ProPane')

  image = Rfits::Rfits_read_image(
    system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package = "ProFound"))

  multi = suppressMessages(
    profoundMultiBand(inputlist = list(Z = image, Z2 = image), detectbands = 'Z',
                      multibands = c('Z', 'Z2'), magzero = 30, cores = 1)
  )

  expect_s3_class(multi, 'profoundmulti')
  expect_named(multi, c('pro_detect', 'cat_tot', 'cat_col', 'cat_grp', 'segimlist',
                        'detectbands', 'multibands', 'call', 'date', 'time'))

  expect_s3_class(multi$pro_detect, 'profound')
  expect_equal(multi$detectbands, 'Z')
  expect_equal(multi$multibands, c('Z', 'Z2'))

  # Each photometry block carries its own per band column suffix: t for total,
  # c for isophotal colour, g for grouped segments
  expect_true(all(c('segID', 'flux_Zt', 'mag_Zt') %in% names(multi$cat_tot)))
  expect_true(all(c('segID', 'flux_Zc', 'mag_Zc') %in% names(multi$cat_col)))
  expect_true(all(c('groupID', 'flux_Zg', 'mag_Zg') %in% names(multi$cat_grp)))

  # The two bands are the same image, so the photometry must be identical
  expect_equal(multi$cat_tot$mag_Zt, multi$cat_tot$mag_Z2t)
  expect_equal(multi$cat_col$mag_Zc, multi$cat_col$mag_Z2c)
  expect_equal(multi$cat_grp$mag_Zg, multi$cat_grp$mag_Z2g)
})

test_that("profoundMultiBand total magnitudes are brighter than colour magnitudes", {
  skip_if_not_installed('Rfits')
  skip_if_not_installed('ProPane')

  image = Rfits::Rfits_read_image(
    system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package = "ProFound"))

  multi = suppressMessages(
    profoundMultiBand(inputlist = list(Z = image), detectbands = 'Z', multibands = 'Z',
                      magzero = 30, cores = 1)
  )

  # Dilated total photometry includes at least as much flux as the isophote
  expect_true(all(multi$cat_tot$flux_Zt >= multi$cat_col$flux_Zc, na.rm = TRUE))
})

test_that("profoundMultiBand can run from a directory of FITS files", {
  skip_if_not_installed('Rfits')
  skip_if_not_installed('ProPane')

  dir = file.path(tempdir(), 'profound-multiband-test')
  dir.create(dir, showWarnings = FALSE)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)

  # The bundled MultiBand images have no usable WCS, so copy the VIKING frame.
  # Note that dir must end with a path separator: it is pasted straight onto the
  # band name when the FITS files are opened.
  from = system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package = "ProFound")
  file.copy(from, file.path(dir, 'Z.fits'), overwrite = TRUE)

  multi = suppressMessages(
    profoundMultiBand(dir = paste0(dir, '/'), detectbands = 'Z', multibands = 'Z',
                      magzero = 30, cores = 1)
  )

  expect_s3_class(multi, 'profoundmulti')
  expect_gt(nrow(multi$cat_tot), 0)
  expect_equal(nrow(multi$cat_tot), multi$pro_detect$Nseg)
})

test_that("profoundMultiBand can run without a detection band if given a segim", {
  skip_if_not_installed('Rfits')
  skip_if_not_installed('ProPane')

  image = Rfits::Rfits_read_image(
    system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package = "ProFound"))
  # profoundMakeSegim wants the bare image matrix, not the Rfits container
  segim = profoundMakeSegim(image$imDat, skycut = 3)$segim

  multi = suppressMessages(
    profoundMultiBand(inputlist = list(Z = image, Z2 = image), segim = segim,
                      detectbands = NULL, multibands = c('Z', 'Z2'), magzero = 30,
                      cores = 1)
  )

  expect_s3_class(multi, 'profoundmulti')
  expect_null(multi$pro_detect)
  # Only total and grouped photometry are possible without a detection image
  expect_true(is.data.frame(multi$cat_tot))
  expect_null(multi$cat_col)
})
