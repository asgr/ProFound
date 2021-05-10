profoundHotFuzz = function(image=NULL, region=NULL, sigma=NULL, segstats=NULL, magzero=0, pixscale=1,
                            rough=TRUE, profound=NULL, nser=2, dofit=TRUE, Niters=c(200,0)){
  
  if(!requireNamespace("ProFit", quietly = TRUE)){
    stop('The ProFit package is needed for smoothing to work. Please install from ICRAR/ProFit', call. = FALSE)
  }
  if(!requireNamespace("Highlander", quietly = TRUE)){
    stop('The Highlander package is needed for smoothing to work. Please install from asgr/Highlander', call. = FALSE)
  }
  
  if(!is.null(image)){
    if(inherits(image, 'profound')){
      profound = image
      image = NULL
    }
  }
  
  if(!is.null(profound)){
    if(is.null(image)){
      image = profound$image - profound$sky
    }
    if(is.null(region)){
      region = profound$objects_redo - (profound$segim_orig > 0)
    }
    if(is.null(sigma)){
      sigma = profound$skyRMS
    }
    if(is.null(segstats)){
      segstats = profound$segstats
    }
    if(missing(magzero)){
      magzero = profound$magzero
    }
    if(missing(pixscale)){
      pixscale = profound$pixscale
    }
  }
  
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
    algo.func = 'LD',
    verbose = FALSE,
    rough = rough
  )
  
  lowers = unlist(Data$intervals)[c(T, F)]
  lowers[unlist(Data$tolog) == T] = log10(lowers[unlist(Data$tolog) == T])
  lowers = lowers[which(unlist(Data$tofit))]
  uppers = unlist(Data$intervals)[c(F, T)]
  uppers[unlist(Data$tolog) == T] = log10(uppers[unlist(Data$tolog) == T])
  uppers = uppers[which(unlist(Data$tofit))]
  
  Data$lowers = lowers
  Data$uppers = uppers
  
  if(dofit){
    highfit = Highlander::Highlander(
      parm = Data$init,
      Data = Data,
      likefunc = ProFit::profitLikeModel,
      Niters=Niters,
      lower = lowers,
      upper = uppers,
      applyintervals = FALSE,
      applyconstraints = FALSE
    )
    
    highfit$modellist = ProFit::profitRemakeModellist(highfit$parm, Data=Data)$modellist
    highfit$image_model = ProFit::profitMakeModel(highfit$modellist, magzero=magzero, dim=dim(image))$z
    
    highfit$Data = Data
    return(highfit)
  }else{
    return(Data)
  }
}
