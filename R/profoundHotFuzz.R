profoundHotFuzz = function(profound=NULL, loc=NULL, box=c(200,200),
                           magzero=0, pixscale=1, size=21, shape="disc", select=NULL,
                           rough=TRUE, nser=2, dofit=TRUE, Niters=c(200,0)){
  
  if(!requireNamespace("ProFit", quietly = TRUE)){
    stop('The ProFit package is needed for smoothing to work. Please install from ICRAR/ProFit', call. = FALSE)
  }
  if(!requireNamespace("Highlander", quietly = TRUE)){
    stop('The Highlander package is needed for smoothing to work. Please install from asgr/Highlander', call. = FALSE)
  }
  
  image = profound$image - profound$sky
  segim_orig = profound$segim_orig
  segim = profound$segim
  sigma = profound$skyRMS
  
  segstats = profound$segstats
  magzero = profound$magzero
  pixscale = profound$pixscale
  
  if(!is.null(select)){
    segstats = segstats[select,]
    segim_orig[!segim_orig %in% segstats$segID] = 0L
    segim[!segim %in% segstats$segID] = 0L
  }
  
  segim_redo = profoundMakeSegimDilate(segim=segim, size=size, shape=shape)$segim
  
  if(!is.null(loc)){
    image = magcutout(image, loc=loc, box=box)
    loc.diff = image$loc.diff
    image = image$image
    segim_orig = magcutout(segim_orig, loc=loc, box=box)$image
    segim = magcutout(segim, loc=loc, box=box)$image
    segim_redo = magcutout(segim_redo, loc=loc, box=box)$image
    sigma = magcutout(sigma, loc=loc, box=box)$image
    
    segstats = segstats[segstats$segID %in% segim_redo,]
    
    segstats[,'xcen'] = segstats[,'xcen'] - loc.diff[1]
    segstats[,'ycen'] = segstats[,'ycen'] - loc.diff[2]
    segstats[,'xmax'] = segstats[,'xmax'] - loc.diff[1]
    segstats[,'ymax'] = segstats[,'ymax'] - loc.diff[2]
  }
  
  region = (segim_redo - segim > 0)
  
  Ncomp = dim(segstats)[1]
  
  modellist = list(
    sersic = list(
      xcen = segstats[,'xcen'],
      ycen = segstats[,'ycen'],
      mag = segstats[,'mag'],
      re = segstats[,'R50'],
      nser = rep(nser, Ncomp),
      ang = segstats[,'ang'],
      axrat = segstats[,'axrat']
    )
  )

  tofit = list(
    sersic = list(
      xcen = rep(FALSE, Ncomp),
      ycen = rep(FALSE, Ncomp),
      mag = rep(TRUE, Ncomp),
      re = rep(TRUE, Ncomp),
      nser = rep(FALSE, Ncomp),
      ang = rep(FALSE, Ncomp),
      axrat = rep(FALSE, Ncomp)
    )
  )
  
  tolog = list(
    sersic = list(
      xcen = rep(FALSE, Ncomp),
      ycen = rep(FALSE, Ncomp),
      mag = rep(FALSE, Ncomp),
      re = rep(TRUE, Ncomp),
      nser = rep(TRUE, Ncomp),
      ang = rep(FALSE, Ncomp),
      axrat = rep(TRUE, Ncomp) #axrat is best fit in log space
    )
  )
  
  maxsize = max(segstats[, 'R50']*4, na.rm=TRUE)
  
  intervals = list(
    sersic = list(
      xcen = rep(list(c(0, dim(image)[1])), Ncomp),
      ycen = rep(list(c(0, dim(image)[2])), Ncomp),
      mag = rep(list(c(0, 40)), Ncomp),
      re = rep(list(c(0.1, maxsize)), Ncomp),
      nser = rep(list(c(0.5, 5.3)), Ncomp),
      ang = rep(list(c(-180, 360)), Ncomp),
      axrat = rep(list(c(0.01, 1)), Ncomp)
    )
  )
  
  constraints = NULL
  
  Data = ProFit::profitSetupData(
    image = image,
    region = region,
    sigma = sigma,
    psf = NULL,
    modellist = modellist,
    tofit = tofit,
    tolog = tolog,
    intervals = intervals,
    constraints = constraints,
    magzero = magzero,
    algo.func = 'CMA',
    verbose = FALSE,
    rough = rough
  )
  
  # lowers = unlist(Data$intervals)[c(T, F)]
  # lowers[unlist(Data$tolog) == T] = log10(lowers[unlist(Data$tolog) == T])
  # lowers = lowers[which(unlist(Data$tofit))]
  # uppers = unlist(Data$intervals)[c(F, T)]
  # uppers[unlist(Data$tolog) == T] = log10(uppers[unlist(Data$tolog) == T])
  # uppers = uppers[which(unlist(Data$tofit))]
  
  lowers = c(segstats[,'mag'] - 5, rep(0, Ncomp))
  uppers = c(segstats[,'mag'] + 5, rep(log10(maxsize), Ncomp))
  
  Data$lowers = lowers
  Data$uppers = uppers
  
  if(dofit){
    highfit = Highlander::Highlander(
      parm = Data$init,
      Data = Data,
      likefunc = ProFit::profitLikeModel,
      likefunctype = 'CMA',
      liketype = 'max',
      Niters = Niters,
      lower = lowers,
      upper = uppers,
      applyintervals = FALSE,
      applyconstraints = FALSE
    )
    
    highfit$image = image
    highfit$segim_orig = segim_orig
    highfit$segim = segim
    highfit$segim_redo = segim_redo
    highfit$sigma = sigma
    highfit$segstats = segstats
    highfit$modellist = ProFit::profitRemakeModellist(highfit$parm, Data=Data)$modellist
    highfit$image_model = ProFit::profitMakeModel(highfit$modellist, magzero=magzero, dim=dim(image))$z
    
    highfit$Data = Data
    return(highfit)
  }else{
    return(list(Data, image=image, segim_orig = segim_orig, segim = segim, segim_redo = segim_redo, sigma = sigma, segstats = segstats))
  }
}
