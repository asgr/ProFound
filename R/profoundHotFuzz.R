profoundHotFuzz = function(profound=NULL, loc=NULL, box=c(200,200),
                           region=NULL, size=21, shape="disc", select=NULL,
                           rough=TRUE, nser=1, dofit=TRUE, Niters=c(200,0), fitRe=TRUE, axrat='profound', ...){
  
  if(!requireNamespace("ProFit", quietly = TRUE)){
    stop('The ProFit package is needed for smoothing to work. Please install from ICRAR/ProFit', call. = FALSE)
  }
  if(!requireNamespace("Highlander", quietly = TRUE)){
    stop('The Highlander package is needed for smoothing to work. Please install from asgr/Highlander', call. = FALSE)
  }
  
  image = profound$image - profound$sky
  #segim_orig = profound$segim_orig
  segim = profound$segim
  sigma = profound$skyRMS
  
  segstats = profound$segstats
  magzero = profound$magzero

  if(!is.null(select)){
    segstats = segstats[select,]
    #segim_orig[!segim_orig %in% segstats$segID] = 0L
    segim[!segim %in% segstats$segID] = 0L
  }
  
  # Make sure we are all good!
  segstats = segstats[is.finite(segstats[,'xcen']) & is.finite(segstats[,'ycen']) & is.finite(segstats[,'mag']) & is.finite(segstats[,'R50']) & is.finite(segstats[,'ang']) & is.finite(segstats[,'axrat']),]
  segim[!segim %in% segstats$segID] = 0L
  
  if(!is.null(loc)){
    image = magcutout(image, loc=loc, box=box)
    loc.diff = image$loc.diff
    image = image$image
    #segim_orig = magcutout(segim_orig, loc=loc, box=box)$image
    segim = magcutout(segim, loc=loc, box=box)$image
    sigma = magcutout(sigma, loc=loc, box=box)$image
    
    segstats = segstats[segstats$segID %in% segim,]
  }
  
  if(is.null(region)){
    segim_redo = profoundMakeSegimDilate(segim=segim, size=size, shape=shape)$segim
    if(!is.null(loc)){
      segim_redo = magcutout(segim_redo, loc=loc, box=box)$image
    }
    region = (segim_redo - segim > 0) & !is.na(image)
    image[region != 1L] = NA
  }else{
    region[is.na(image)] = FALSE
    segim_redo = NULL
  }
  
  Ncomp = dim(segstats)[1]
  
  modellist = list(
    sersic = list(
      xcen = if(is.null(loc)){segstats[,'xcen']}else{segstats[,'xcen'] - loc.diff[1]},
      ycen = if(is.null(loc)){segstats[,'ycen']}else{segstats[,'ycen'] - loc.diff[2]},
      mag = segstats[,'mag'],
      re = segstats[,'R50'],
      nser = rep(nser, Ncomp),
      ang = segstats[,'ang'],
      axrat = if(axrat=='profound'){segstats[,'axrat']}else{axrat}
    )
  )

  tofit = list(
    sersic = list(
      xcen = rep(FALSE, Ncomp),
      ycen = rep(FALSE, Ncomp),
      mag = rep(TRUE, Ncomp),
      re = rep(fitRe, Ncomp),
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
  
  if(fitRe){
    lowers = c(segstats[,'mag'] - 5, rep(-1, Ncomp))
    uppers = c(segstats[,'mag'] + 5, rep(log10(maxsize), Ncomp))
  }else{
    lowers = segstats[,'mag'] - 5
    uppers = segstats[,'mag'] + 5
  }
  
  Data$lowers = lowers
  Data$uppers = uppers

  if(dofit){
    if(requireNamespace("Highlander", quietly = TRUE)){
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
        applyconstraints = FALSE,
        ...
      )
      
      highfit$modellist = ProFit::profitRemakeModellist(highfit$parm, Data=Data)$modellist
      highfit$image_model = ProFit::profitMakeModel(highfit$modellist, magzero=magzero, dim=dim(image))$z
      
      highfit$Data = Data
      highfit$image = image
      #highfit$segim_orig = segim_orig
      highfit$segim = segim
      highfit$segim_redo = segim_redo
      highfit$sigma = sigma
      highfit$segstats = segstats

      return(highfit)
    }else{
      stop('The Highlander package is needed to fit with HotFuzz!')
    }
  }else{
    return(list(Data=Data, image=image, segim = segim, segim_redo = segim_redo, sigma = sigma, segstats = segstats))
  }
}

profoundPSFFuzz = function(profound=NULL, loc=NULL, box=c(200,200),
                           size=21, shape="disc", select=NULL,
                           rough=TRUE, nser=2, dofit=TRUE, Niters=c(200,0), ...){
  
  if(!requireNamespace("ProFit", quietly = TRUE)){
    stop('The ProFit package is needed for smoothing to work. Please install from ICRAR/ProFit', call. = FALSE)
  }
  if(!requireNamespace("Highlander", quietly = TRUE)){
    stop('The Highlander package is needed for smoothing to work. Please install from asgr/Highlander', call. = FALSE)
  }
  
  image = profound$image - profound$sky
  #segim_orig = profound$segim_orig
  segim = profound$segim
  sigma = profound$skyRMS
  
  segstats = profound$segstats
  magzero = profound$magzero

  if(!is.null(select)){
    segstats = segstats[select,]
    #segim_orig[!segim_orig %in% segstats$segID] = 0L
    segim[!segim %in% segstats$segID] = 0L
  }
  
  # Make sure we are all good!
  segstats = segstats[is.finite(segstats[,'xcen']) & is.finite(segstats[,'ycen']) & is.finite(segstats[,'mag']) & is.finite(segstats[,'R50']) & is.finite(segstats[,'ang']) & is.finite(segstats[,'axrat']),]
  segim[!segim %in% segstats$segID] = 0L
  
  segim_redo = profoundMakeSegimDilate(segim=segim, size=size, shape=shape)$segim
  
  if(!is.null(loc)){
    image = magcutout(image, loc=loc, box=box)
    loc.diff = image$loc.diff
    image = image$image
    #segim_orig = magcutout(segim_orig, loc=loc, box=box)$image
    segim = magcutout(segim, loc=loc, box=box)$image
    segim_redo = magcutout(segim_redo, loc=loc, box=box)$image
    sigma = magcutout(sigma, loc=loc, box=box)$image
    
    segstats = segstats[segstats$segID %in% segim_redo,]
  }
  
  region = (segim_redo - segim > 0)
  
  Ncomp = dim(segstats)[1]
  
  modellist = list(
    sersic = list(
      xcen = if(is.null(loc)){segstats[,'xcen']}else{segstats[,'xcen'] - loc.diff[1]},
      ycen = if(is.null(loc)){segstats[,'ycen']}else{segstats[,'ycen'] - loc.diff[2]},
      mag = segstats[,'mag'],
      re = rep(max(segstats[,'R50']), Ncomp),
      nser = rep(nser, Ncomp),
      ang = rep(0, Ncomp ),
      axrat = rep(1, Ncomp )
    )
  )
  
  tofit = list(
    sersic = list(
      xcen = rep(FALSE, Ncomp),
      ycen = rep(FALSE, Ncomp),
      mag = rep(FALSE, Ncomp),
      re = c(TRUE,rep(NA, Ncomp - 1)),
      nser = c(TRUE,rep(NA, Ncomp - 1)),
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
  
  minsize = max(segstats[, 'R50']/4, na.rm=TRUE)
  maxsize = max(segstats[, 'R90']*4, na.rm=TRUE)
  
  intervals = list(
    sersic = list(
      xcen = rep(list(c(0, dim(image)[1])), Ncomp),
      ycen = rep(list(c(0, dim(image)[2])), Ncomp),
      mag = rep(list(c(0, 40)), Ncomp),
      re = rep(list(c(minsize, maxsize)), Ncomp),
      nser = rep(list(c(0.5, 10)), Ncomp),
      ang = rep(list(c(-180, 360)), Ncomp),
      axrat = rep(list(c(0.5, 1)), Ncomp)
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
  
  lowers = c(log10(minsize), log10(0.5))
  uppers = c(log10(maxsize), log10(10))
  
  Data$lowers = lowers
  Data$uppers = uppers
  
  if(dofit){
    if(requireNamespace("Highlander", quietly = TRUE)){
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
      
      highfit$modellist = ProFit::profitRemakeModellist(highfit$parm, Data=Data)$modellist
      highfit$image_model = ProFit::profitMakeModel(highfit$modellist, magzero=magzero, dim=dim(image))$z
      
      highfit$Data = Data
      highfit$image = image
      highfit$segim = segim
      highfit$segim_redo = segim_redo
      highfit$sigma = sigma
      highfit$segstats = segstats
      
      return(highfit)
    }else{
      stop('The Highlander package is needed to fit with HotFuzz!')
    }
  }else{
    return(list(Data=Data, image=image, segim=segim, segim_redo=segim_redo, sigma=sigma, segstats=segstats))
  }
}
