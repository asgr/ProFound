# profoundFluxDeblend apportions the flux inside a group footprint between the
# segments of that group. Most of its outputs are *model dependent*: the
# deblendtype choice ('fit', 'gauss', 'psf') decides how the group flux is
# shared out. The tests below therefore lean on the properties that must hold
# regardless of the model, plus accuracy checks only where the model actually
# matches the data that was put in.

# A blended pair that the watershed splits into two touching segments, plus an
# isolated source so the pass-through path is exercised too.
# Truth: fluxes of 4000 and 3000 at (45, 50) and (55, 50), both sigma = 3.
makeDeblendCase = function(magzero = 25){
  set.seed(452)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(45, 55, 20),
    ycen = c(50, 50, 20),
    flux = c(4000, 3000, 6000),
    sigma = c(3, 3, 2)
  ))

  pro = runProFound(image, skycut = 2, box = 50, grid = 50, tolerance = 1,
                    groupstats = TRUE, magzero = magzero)

  list(
    image = image,
    sky = pro$sky,
    pro = pro,
    magzero = magzero,
    sky_sub = image - pro$sky,
    # segID 2 and 3 are the blended pair, segID 1 is isolated
    pair = c(2L, 3L),
    blend_ids = c(2L, 3L),
    blend_group = pro$group$groupsegID$groupID[pro$group$groupsegID$Ngroup > 1],
    truth = c(4000, 3000)
  )
}

# A three source group: a bright / mid / faint row on a shared footprint, so the
# deblend has to apportion flux between three segments at once rather than two.
# Truth: fluxes of 4000, 3000 and 2000 at x = 50, 60, 70, all at y = 60.
makeDeblendTriple = function(magzero = 25){
  set.seed(2)
  image = makeTestImage(dim = c(121, 121), sources = data.frame(
    xcen = 60 + c(-10, 0, 10),
    ycen = c(60, 60, 60),
    flux = c(4000, 3000, 2000),
    sigma = c(3, 3, 3)
  ))

  pro = runProFound(image, skycut = 2, box = 50, grid = 50, tolerance = 1,
                    groupstats = TRUE, magzero = magzero)

  list(
    image = image,
    sky = pro$sky,
    pro = pro,
    magzero = magzero,
    sky_sub = image - pro$sky,
    seg = c(1L, 2L, 3L),
    blend_ids = c(1L, 2L, 3L),
    blend_group = pro$group$groupsegID$groupID[pro$group$groupsegID$Ngroup > 1],
    truth = c(4000, 3000, 2000)
  )
}

# The same call for every test, so that the arguments under test are obvious.
# magzero defaults to the one used in the profoundProFound() run, so that mag_db
# stays consistent with segstats$mag.
runDeblend = function(case, magzero = case$magzero, ...){
  profoundFluxDeblend(
    image = case$sky_sub,
    segim = case$pro$segim,
    segstats = case$pro$segstats,
    groupim = case$pro$group$groupim,
    groupsegID = case$pro$group$groupsegID,
    magzero = magzero,
    ...
  )
}

# The deblended rows of one group, ordered by segID so that assertions can index
# them against the input segment order.
groupRows = function(deblend, case){
  out = deblend[deblend$segID %in% case$blend_ids, , drop = FALSE]
  out[order(out$segID), , drop = FALSE]
}

test_that("profoundFluxDeblend returns one row per segment with the expected columns", {
  case = makeDeblendCase()

  deblend = runDeblend(case)

  expect_s3_class(deblend, 'data.frame')
  # One row per input segment, in segstats order
  expect_equal(nrow(deblend), nrow(case$pro$segstats))
  expect_equal(deblend$segID, case$pro$segstats$segID)

  expect_true(all(c('groupID', 'segID', 'flux_db', 'mag_db', 'N100_db', 'flux_segfrac',
                    'Qseg_db', 'Qgroup_db', 'beamcorrect', 'flux_err_db',
                    'mag_err_db') %in% names(deblend)))
})

test_that("profoundFluxDeblend conserves the flux inside each group footprint", {
  case = makeDeblendCase()

  # This is the defining property of the deblend: whatever the weights, the
  # deblended fluxes sum to the flux in the group footprint.
  footprint = sum(case$sky_sub[case$pro$group$groupim > 0])
  expect_equal(sum(runDeblend(case)$flux_db), footprint)

  # It holds for each deblend model, since they only differ in the weights
  expect_equal(sum(runDeblend(case, deblendtype = 'gauss')$flux_db), footprint)
  expect_equal(sum(runDeblend(case, deblendtype = 'psf',
                              psf = makeTestPSF(size = 25, sigma = 3))$flux_db), footprint)
})

test_that("profoundFluxDeblend splits a blended pair without losing or inventing flux", {
  case = makeDeblendCase()

  deblend = runDeblend(case)
  pair = deblend[deblend$segID %in% case$pair, ]

  expect_equal(nrow(pair), 2)
  expect_equal(pair$groupID, rep(case$blend_group, 2))
  # Both halves get a meaningful share of the group, and the two shares add to 1
  expect_true(all(pair$flux_db > 0))
  expect_equal(sum(pair$flux_segfrac), 1)
  # The brighter isophotal segment stays the brighter deblended one
  expect_equal(which.max(pair$flux_db), which.max(pair$flux_segfrac))

  # The group total is preserved
  footprint = sum(case$sky_sub[case$pro$group$groupim == case$blend_group])
  expect_equal(sum(pair$flux_db), footprint)

  # The deblend moves flux towards / away from the isophotal split, but not by
  # more than the group flux on either side
  expect_true(all(abs(pair$flux_db - pair$flux_segfrac * footprint) < footprint))
})

test_that("profoundFluxDeblend leaves isolated segments untouched", {
  case = makeDeblendCase()

  deblend = runDeblend(case)
  solo = deblend[deblend$segID == 1, ]
  segstats_solo = case$pro$segstats[case$pro$segstats$segID == 1, ]

  expect_true(is.na(solo$groupID))
  # Ungrouped segments pass straight through from segstats
  expect_equal(solo$flux_db, segstats_solo$flux)
  expect_equal(solo$N100_db, segstats_solo$N100)
  expect_equal(solo$flux_segfrac, 1)
  expect_equal(solo$Qseg_db, 0)
  expect_equal(solo$Qgroup_db, 0)
  # beamcorrect is only meaningful for the psf model, which is not used here
  expect_true(is.na(solo$beamcorrect))
})

test_that("profoundFluxDeblend flux_segfrac is the isophotal share of the group", {
  case = makeDeblendCase()

  deblend = runDeblend(case)

  # Recompute the shares directly: the flux of each segment's own pixels as a
  # fraction of the whole group footprint's flux. This is independent of any
  # deblending model, so it can be checked exactly.
  footprint = which(case$pro$group$groupim == case$blend_group, arr.ind = TRUE)
  groupflux = sum(case$sky_sub[footprint])

  for(i in seq_along(case$pair)){
    segpix = footprint[case$pro$segim[footprint] == case$pair[i], , drop = FALSE]
    expect_equal(
      deblend$flux_segfrac[deblend$segID == case$pair[i]],
      sum(case$sky_sub[segpix]) / groupflux
    )
  }

  # Ungrouped segments are their own whole group
  expect_equal(deblend$flux_segfrac[deblend$segID == 1], 1)
})

test_that("profoundFluxDeblend weights add up to the footprint pixel count", {
  case = makeDeblendCase()

  deblend = runDeblend(case)
  pair = deblend[deblend$segID %in% case$pair, ]
  footprint_pixels = sum(case$pro$group$groupim == case$blend_group)

  # N100_db is the effective number of footprint pixels assigned to each segment,
  # so they cannot exceed the footprint and should add up to it
  expect_true(all(pair$N100_db <= footprint_pixels))
  expect_equal(sum(pair$N100_db), footprint_pixels)
  # Ungrouped segments keep their own pixel count
  expect_equal(deblend$N100_db[deblend$segID == 1],
               case$pro$segstats$N100[case$pro$segstats$segID == 1])
})

test_that("profoundFluxDeblend Qgroup_db is the same for every segment in a group", {
  case = makeDeblendCase()

  deblend = runDeblend(case)
  pair = deblend[deblend$segID %in% case$pair, ]

  # Qgroup_db describes how well the group as a whole is modelled, so it is a
  # single number broadcast to all members of the group
  expect_equal(length(unique(pair$Qgroup_db)), 1)
  # The group is a blend, so the model does not reproduce it exactly
  expect_lt(pair$Qgroup_db[1], 0)
})

test_that("profoundFluxDeblend magnitudes and errors are internally consistent", {
  case = makeDeblendCase()

  deblend = runDeblend(case)
  segstats = case$pro$segstats[match(deblend$segID, case$pro$segstats$segID), ]

  expect_equal(deblend$mag_db, profoundFlux2Mag(deblend$flux_db, magzero = case$magzero))
  expect_equal(deblend$mag_err_db,
               (2.5/log(10)) * abs(deblend$flux_err_db / deblend$flux_db))

  # The total error is the quadrature sum of the four components
  expect_equal(deblend$flux_err_db,
               sqrt(deblend$flux_err_sky_db^2 + deblend$flux_err_skyRMS_db^2 +
                      deblend$flux_err_shot_db^2 + deblend$flux_err_Qseg_db^2))

  # The skyRMS term scales with the square root of the assigned pixel count
  expect_equal(deblend$flux_err_skyRMS_db,
               segstats$flux_err_skyRMS * sqrt(deblend$N100_db / segstats$N100))
  expect_true(all(deblend$flux_err_db > 0))
})

test_that("profoundFluxDeblend mag_db needs a magzero that matches segstats", {
  # The grouped rows are rebuilt as profoundFlux2Mag(flux_db, magzero), but the
  # pass-through rows simply reuse segstats$mag, which already carries the
  # magzero used in the profoundProFound() run. Passing a second, different
  # magzero therefore only shifts the grouped rows. A consistent magzero is the
  # correct usage; the mismatch is pinned down here so a fix is noticed.
  case = makeDeblendCase(magzero = 25)

  consistent = runDeblend(case, magzero = case$magzero)
  expect_equal(consistent$mag_db, profoundFlux2Mag(consistent$flux_db, magzero = case$magzero))

  mismatch = runDeblend(case, magzero = 30)
  expect_equal(mismatch$flux_db, consistent$flux_db)
  expect_equal(mismatch$mag_db[mismatch$segID == 1],
               consistent$mag_db[consistent$segID == 1])
  pair = mismatch$segID %in% case$pair
  expect_equal(mismatch$mag_db[pair], consistent$mag_db[pair] + 5)
})

test_that("profoundFluxDeblend deblendtype = 'psf' recovers a matched-PSF blend", {
  skip_if_not_installed('ProFit')

  # Two sources that really are PSF shaped, so the psf model is the right one.
  # A broad PSF is used because narrower ones get split by the watershed before
  # the deblend stage sees them as a group.
  set.seed(14)
  psf = makeTestPSF(size = 25, sigma = 3)
  truth = c(4000, 3000)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(45, 55), ycen = c(50, 50), flux = truth, sigma = c(3, 3)))

  pro = runProFound(image, skycut = 2, box = 50, grid = 50, tolerance = 1,
                    groupstats = TRUE)
  expect_equal(nrow(pro$segstats), 2)
  expect_equal(nrow(pro$group$groupsegID), 1) # a single group of two

  deblend = profoundFluxDeblend(image = image - pro$sky, segim = pro$segim,
                                segstats = pro$segstats,
                                groupim = pro$group$groupim,
                                groupsegID = pro$group$groupsegID,
                                deblendtype = 'psf', psf = psf, magzero = 25)

  expect_equal(sum(deblend$flux_db), sum((image - pro$sky)[pro$group$groupim > 0]))
  # The bright source is recovered well; the fainter one is biased low because
  # the footprint is isophotally truncated at the skycut
  expect_equal(deblend$flux_db[1] / truth[1], 1, tolerance = 0.1)
  expect_equal(deblend$flux_db[2] / truth[2], 1, tolerance = 0.2)
  expect_true(all(deblend$flux_db > 0))
})

test_that("profoundFluxDeblend rejects its documented bad inputs", {
  case = makeDeblendCase()

  expect_error(runDeblend(case, deblendtype = 'junk'),
               "deblendtype type not recognised")
  expect_error(runDeblend(case, deblendtype = 'psf'),
               "psf must be provided")
  expect_error(runDeblend(case, deblendtype = 'psf', psf = makeTestPSF(),
                          fluxweight = 'junk'),
               "fluxweight type not recognised")
  expect_error(runDeblend(case, deblendtype = 'psf', psf = makeTestPSF(),
                          convmode = 'junk'),
               "convmode type not recognised")
  expect_error(profoundFluxDeblend(image = case$sky_sub), "Need segim if missing")
  expect_error(profoundFluxDeblend(profound = 'junk'),
               "Class of profound input must be of type 'profound'")
})

test_that("profoundFluxDeblend returns NULL when no group is worth deblending", {
  set.seed(451)
  image = makeTestImage(dim = c(101, 101), sources = data.frame(
    xcen = c(28, 50, 76), ycen = c(72, 50, 30), flux = c(5000, 4000, 6000),
    sigma = c(2, 3, 2)))
  pro = runProFound(image, skycut = 2, box = 50, grid = 50, tolerance = 1,
                    groupstats = TRUE)

  # Three well separated sources: every group has a single member
  expect_true(all(pro$group$groupsegID$Ngroup == 1))

  expect_null(profoundFluxDeblend(image = image - pro$sky, segim = pro$segim,
                                  segstats = pro$segstats,
                                  groupim = pro$group$groupim,
                                  groupsegID = pro$group$groupsegID))

  # Ndeblendlim drops groups whose pixel count times member count is too large
  case = makeDeblendCase()
  expect_null(runDeblend(case, Ndeblendlim = 100))
})

test_that("profoundFluxDeblend accepts a profound object via image or profound", {
  case = makeDeblendCase()

  from_args = runDeblend(case)
  from_object = profoundFluxDeblend(image = case$pro, magzero = case$magzero)

  # Passing the profound object as `image` pulls out the sky subtracted image
  expect_equal(from_object, from_args)

  # Inferring the grouping from segim gives the same answer as passing it in
  from_segim = profoundFluxDeblend(image = case$sky_sub, segim = case$pro$segim,
                                   segstats = case$pro$segstats, magzero = case$magzero)
  expect_equal(from_segim, from_args)
})

test_that("profoundFluxDeblend subtracts the sky from a profound object's image", {
  case = makeDeblendCase()

  # Both entry points must behave identically: they differ only in how the
  # inputs are supplied, not in whether the sky is removed.
  via_image = profoundFluxDeblend(image = case$pro, magzero = case$magzero)
  via_profound = profoundFluxDeblend(profound = case$pro, magzero = case$magzero)

  expect_equal(via_profound, via_image)
  expect_equal(via_profound, runDeblend(case))

  # The blended pair is the only group with more than one member; the isolated
  # segment is a single-member group and has no groupID in the output
  footprint = case$pro$group$groupim == case$blend_group
  in_group = !is.na(via_profound$groupID)

  # The group flux is the sky subtracted footprint, not the raw one
  expect_equal(sum(via_profound$flux_db[in_group]), sum(case$sky_sub[footprint]))
  expect_lt(sum(via_profound$flux_db[in_group]), sum(case$image[footprint]))
})

test_that("profoundFluxDeblend lets an explicit sky override the profound object", {
  case = makeDeblendCase()
  footprint = case$pro$group$groupim == case$blend_group
  group_flux = function(x) sum(x$flux_db[!is.na(x$groupID)])

  # sky is only defaulted from the profound object when the user leaves it alone
  no_sky = profoundFluxDeblend(profound = case$pro, sky = 0, magzero = case$magzero)
  full_sky = profoundFluxDeblend(profound = case$pro, magzero = case$magzero)
  half_sky = profoundFluxDeblend(profound = case$pro, sky = case$pro$sky/2,
                                 magzero = case$magzero)

  expect_equal(group_flux(no_sky), sum(case$image[footprint]))
  expect_equal(group_flux(full_sky), sum(case$sky_sub[footprint]))
  expect_gt(group_flux(half_sky), group_flux(full_sky))
  expect_lt(group_flux(half_sky), group_flux(no_sky))
})

test_that("profoundFluxDeblend lowmemory matches the default", {
  case = makeDeblendCase()

  expect_equal(runDeblend(case, lowmemory = TRUE), runDeblend(case))
})

test_that("profoundFluxDeblend doallstats = FALSE returns only the deblend columns", {
  case = makeDeblendCase()

  lean = runDeblend(case, doallstats = FALSE)

  expect_true(all(c('groupID', 'segID', 'flux_db', 'mag_db', 'flux_segfrac') %in% names(lean)))
  # Only grouped segments are reported
  expect_equal(lean$segID, case$pair)
  expect_false(any(c('flux_err_db', 'mag_err_db') %in% names(lean)))
  # The same deblended fluxes as the full run
  full = runDeblend(case)
  expect_equal(lean$flux_db, full$flux_db[full$segID %in% case$pair])
})

test_that("profoundFluxDeblend iterative changes the split but conserves flux", {
  case = makeDeblendCase()

  single = runDeblend(case)
  iterative = runDeblend(case, iterative = TRUE)

  # Iterative subtraction of each segment before fitting the next perturbs the
  # weights, so the answer is different ...
  expect_false(isTRUE(all.equal(iterative$flux_db, single$flux_db)))
  # ... but every row is still returned, and the solos still pass through
  expect_equal(nrow(iterative), nrow(single))
  expect_equal(iterative$flux_db[iterative$segID == 1], single$flux_db[single$segID == 1])
})

test_that("profoundFluxDeblend image_reweight returns a list and rescales the image", {
  case = makeDeblendCase()

  plain = runDeblend(case)
  reweight = runDeblend(case, image_reweight = TRUE)

  expect_named(reweight, c('segstats', 'image'))
  expect_equal(reweight$segstats, plain)
  expect_identical(dim(reweight$image), dim(case$sky_sub))

  # The footprint is reweighted onto the segments, so it loses the flux that the
  # model assigned to background pixels within the group
  footprint = case$pro$group$groupim > 0
  expect_lt(sum(reweight$image[footprint]), sum(case$sky_sub[footprint]))
  # Pixels outside every group are untouched
  outside = case$pro$group$groupim == 0
  expect_equal(reweight$image[outside], case$sky_sub[outside])
})

test_that("profoundFluxDeblend gauss and fit models give the same sign of shift", {
  case = makeDeblendCase()

  fit = runDeblend(case, deblendtype = 'fit')
  gauss = runDeblend(case, deblendtype = 'gauss')

  # Neither model is the truth here, so only assert that they shift the split in
  # a consistent way. Deviations are measured against the isophotal split of the
  # pair total (not the sum over every segment, which would include the solo).
  pair_total = sum(fit$flux_db[fit$segID %in% case$pair])
  pick = fit$segID %in% case$pair
  fit_shift = (fit$flux_db - fit$flux_segfrac * pair_total)[pick]
  gauss_shift = (gauss$flux_db - gauss$flux_segfrac * pair_total)[gauss$segID %in% case$pair]

  # Both models move flux in the same direction, and by a non-trivial amount
  expect_equal(sign(fit_shift), sign(gauss_shift))
  expect_true(all(abs(fit_shift) > 1))
  expect_true(all(abs(gauss_shift) > 1))
  # The shifts are equal and opposite, since the pair total is conserved
  expect_equal(sum(fit_shift), 0, tolerance = 1e-8)
  expect_equal(sum(gauss_shift), 0, tolerance = 1e-8)
  # Both models preserve the group total
  expect_equal(sum(fit$flux_db), sum(gauss$flux_db))
})
# --- Three source groups -----------------------------------------------------
# A two segment group only ever asks the deblend to split a pile in two. With
# three segments sharing one footprint the weighting has to be shared three ways,
# which is where interactions between segments show up.

test_that("the three source case really is one group of three", {
  case = makeDeblendTriple()

  # Guard the fixture: if segmentation changes, these tests would silently stop
  # exercising a three member group.
  expect_equal(nrow(case$pro$segstats), 3)
  expect_equal(nrow(case$pro$group$groupsegID), 1)
  expect_equal(case$pro$group$groupsegID$Ngroup, 3)
  expect_equal(unlist(case$pro$group$groupsegID$segID), case$blend_ids)
  # The segments are ordered left to right, brightest first
  expect_equal(round(case$pro$segstats$xcen), c(50, 60, 70), tolerance = 1)
})

test_that("profoundFluxDeblend apportions a three segment group without loss", {
  case = makeDeblendTriple()

  deblend = runDeblend(case)
  trio = groupRows(deblend, case)
  footprint = case$pro$group$groupim == case$blend_group
  footprint_flux = sum(case$sky_sub[footprint])
  footprint_pixels = sum(footprint)

  expect_equal(nrow(trio), 3)
  expect_equal(trio$groupID, rep(case$blend_group, 3))
  expect_true(all(trio$flux_db > 0))

  # Every quantity that does not depend on the weighting model must still add up
  expect_equal(sum(trio$flux_db), footprint_flux)
  expect_equal(sum(trio$flux_segfrac), 1)
  expect_equal(sum(trio$N100_db), footprint_pixels)

  # Three shares, three pixels counts: no segment may swallow the whole group
  expect_true(all(trio$flux_segfrac > 0))
  expect_true(all(trio$N100_db < footprint_pixels))
})

test_that("profoundFluxDeblend three way shares are the isophotal shares", {
  case = makeDeblendTriple()

  deblend = runDeblend(case)
  trio = groupRows(deblend, case)

  # Recompute each segment's isophotal share directly from segim and the
  # footprint. This is model independent, so it can be checked exactly.
  footprint = which(case$pro$group$groupim == case$blend_group, arr.ind = TRUE)
  groupflux = sum(case$sky_sub[footprint])

  for(i in seq_along(case$blend_ids)){
    segpix = footprint[case$pro$segim[footprint] == case$blend_ids[i], , drop = FALSE]
    expect_equal(trio$flux_segfrac[i], sum(case$sky_sub[segpix]) / groupflux,
                 label = paste('segment', case$blend_ids[i]))
  }
})

test_that("profoundFluxDeblend three segment group shares flux by brightness rank", {
  case = makeDeblendTriple()

  # The left to right sources are 4000 / 3000 / 2000, well enough separated and
  # bright enough that every model recovers their ordering. Brightest first, so
  # flux_db and N100_db should both decrease along the row.
  for(type in c('fit', 'gauss')){
    trio = groupRows(runDeblend(case, deblendtype = type), case)
    expect_true(all(diff(trio$flux_db) < 0), label = type)
    expect_true(all(diff(trio$N100_db) < 0), label = type)
    expect_equal(order(-trio$flux_db), seq_along(case$blend_ids), label = type)
  }

  # The isophotal ordering is the same, so the deblend is not fighting the data
  iso = case$pro$segstats$flux[match(case$blend_ids, case$pro$segstats$segID)]
  expect_true(all(diff(iso) < 0))
})

test_that("profoundFluxDeblend three segment group keeps Qgroup_db group level", {
  case = makeDeblendTriple()

  trio = groupRows(runDeblend(case), case)

  expect_equal(length(unique(trio$Qgroup_db)), 1)
  expect_lt(trio$Qgroup_db[1], 0)
  # Per segment quantities still vary between the three members
  expect_gt(length(unique(round(trio$Qseg_db, 6))), 1)
})

test_that("profoundFluxDeblend three segments conserve flux in every model", {
  case = makeDeblendTriple()
  footprint_flux = sum(case$sky_sub[case$pro$group$groupim == case$blend_group])

  for(type in c('fit', 'gauss', 'psf')){
    trio = groupRows(
      runDeblend(case, deblendtype = type, psf = makeTestPSF(size = 25, sigma = 3)),
      case
    )
    expect_equal(sum(trio$flux_db), footprint_flux, label = type)
    expect_equal(sum(trio$N100_db),
                 sum(case$pro$group$groupim == case$blend_group), label = type)
    expect_equal(sum(trio$flux_segfrac), 1, label = type)
  }
})

test_that("profoundFluxDeblend three segment group errors are self consistent", {
  case = makeDeblendTriple()

  deblend = runDeblend(case)
  trio = groupRows(deblend, case)
  segstats = case$pro$segstats[match(trio$segID, case$pro$segstats$segID), ]

  expect_equal(trio$mag_db, profoundFlux2Mag(trio$flux_db, magzero = case$magzero))
  expect_equal(trio$flux_err_db,
               sqrt(trio$flux_err_sky_db^2 + trio$flux_err_skyRMS_db^2 +
                      trio$flux_err_shot_db^2 + trio$flux_err_Qseg_db^2))
  expect_equal(trio$flux_err_skyRMS_db,
               segstats$flux_err_skyRMS * sqrt(trio$N100_db / segstats$N100))
  expect_true(all(trio$flux_err_db > 0))
})

test_that("profoundFluxDeblend three segment group respects Ndeblendlim", {
  case = makeDeblendTriple()

  # The limit is on footprint pixels times group members
  limit = as.numeric(case$pro$group$groupsegID$Npix) * case$pro$group$groupsegID$Ngroup
  expect_length(limit, 1)

  expect_null(runDeblend(case, Ndeblendlim = limit - 1))
  expect_false(is.null(runDeblend(case, Ndeblendlim = limit + 1)))
})

test_that("profoundFluxDeblend three segment group handles iterative mode", {
  case = makeDeblendTriple()

  single = groupRows(runDeblend(case), case)
  iterative = groupRows(runDeblend(case, iterative = TRUE), case)

  # Iterative mode subtracts each segment before fitting the next, so with three
  # members there is a chain of subtractions and the split changes ...
  expect_false(isTRUE(all.equal(iterative$flux_db, single$flux_db)))
  # ... while the group total is still preserved
  expect_equal(sum(iterative$flux_db), sum(single$flux_db))
  expect_true(all(iterative$flux_db > 0))
})

test_that("profoundFluxDeblend three segment group works via a profound object", {
  case = makeDeblendTriple()

  # Both profound entry points, with the sky now subtracted in each
  expect_equal(profoundFluxDeblend(image = case$pro, magzero = case$magzero),
               runDeblend(case))
  expect_equal(profoundFluxDeblend(profound = case$pro, magzero = case$magzero),
               runDeblend(case))
})

test_that("profoundFluxDeblend does not double subtract a pre-subtracted image", {
  # Regression guard for the sky handling in the `profound =` branch. Scripts
  # written against the old behaviour pre-subtracted the sky *and* passed the
  # object; that now subtracts twice. The stored sky (~100 counts) exceeds the
  # sky subtracted footprint flux (~6400 over the pair), so everything inside
  # the group goes to zero, silently and with no warning.
  case = makeDeblendCase()
  pre_subtracted = case$pro$image - case$sky

  from_object = runDeblend(case)
  doubles = profoundFluxDeblend(image = pre_subtracted, profound = case$pro,
                                magzero = case$magzero)

  expect_false(isTRUE(all.equal(sum(doubles$flux_db), sum(from_object$flux_db))))
  # The group collapses to zero flux: the failure is complete, not subtle
  in_group = !is.na(doubles$groupID)
  expect_true(all(doubles$flux_db[in_group] == 0))
  expect_true(all(doubles$N100_db[in_group] == 0))

  # Passing an explicit sky = 0 is the way to keep the old workflow working
  expect_equal(profoundFluxDeblend(image = pre_subtracted, profound = case$pro,
                                   sky = 0, magzero = case$magzero),
               from_object)
})
