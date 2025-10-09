.fluxcalcapp = function(x=NULL, y=NULL, rad2=NULL, flux=NULL, xcen=NA, ycen=NA, rad_app=NULL, centype='mean'){
  if(is.null(rad2)){
    if(is.na(xcen)){
      if(centype == 'wt' | centype == 'mean'){
        xcen = .meanwt(x, flux)
      }else if(centype == 'max'){
        xcen = x[which.max(flux)]
      }else{
        stop('centype must be max or mean!')
      }
    }
    
    if(is.na(ycen)){
      if(centype == 'wt' | centype == 'mean'){
        ycen = .meanwt(y, flux)
      }else if(centype == 'max'){
        ycen = y[which.max(flux)]
      }else{
        stop('centype must be max or mean!')
      }
    }
    
    if(xcen == 0 & ycen == 0){
      rad2 = x^2 + y^2
    }else{
      rad2 = (x - xcen)^2 + (y - ycen)^2 
    }
  }else{
    xcen = 0
    ycen = 0
  }
  
  if(is.na(rad_app)){
    #left here for a different mode I tested for applying rad2 selection outside of this function.
    #that wasn't any faster, so probably remove this if case in the future
    Nsel = length(rad2)
    rad_out = max(rad2)
    sel_out = which(rad2 == rad_out)
    
    flux_app = sum(flux, na.rm=TRUE)
    flux_min = mean(flux[sel_out], na.rm=TRUE)
  }else{
    sel = which(rad2 <= rad_app^2)
    Nsel = length(sel)
    
    if(length(Nsel) == 0){
      flux_app = NA_real_
      flux_min = 0
    }else{
      if(all(is.na(flux[sel]))){
        flux_app = NA_real_
        flux_min = 0
      }else{
        rad_out = max(rad2[sel])
        sel_out = which(rad2 == rad_out)
        
        flux_app = sum(flux[sel], na.rm=TRUE)
        flux_min = mean(flux[sel_out], na.rm=TRUE)
      }
    }
  }
  
  if(!isTRUE(is.finite(flux_min))){ #this catches NA, NaN, NULL, Inf events
    flux_min = 0 #don't want to penalise when masked or other weird events
  }else if(flux_min < 0){
    flux_min = 0 #don't want to penalise when in the sky noise 
  }
  
  return(list(flux_app=flux_app, flux_min=flux_min, N=Nsel))
}

profoundAperPhot = function(image=NULL, segim=NULL, app_diam=1, mask=NULL, keyvalues=NULL, tar=NULL,
                           pixscale=1, magzero=0, correction=TRUE, centype='mean', fluxtype='Raw',
                           verbose=FALSE){
  if(!is.null(image)){
    if(inherits(image, 'Rfits_image')){
      keyvalues = image$keyvalues
      image = image$imDat
    }else if(inherits(image, 'Rfits_pointer')){
      keyvalues = image$keyvalues
      image = image[,]$imDat
    }else if(inherits(image, 'matrix')){
      'Do nothing'
    }else{
      stop('As of ProFound v1.21.0 only Rfits_image FITS inputs are allowed. Please install from GitHub asgr/Rfits')
    }
  }else{
    stop('Missing image - this is a required input!')
  }
  
  if(!is.null(keyvalues)){
    if(!inherits(keyvalues, 'Rfits_keylist')){
      if(is.list(keyvalues)){
        class(keyvalues) = 'Rfits_keylist'
      }else{
        stop('keyvalues is the wrong format- should be a list!')
      }
    }
  }
  
  if(missing(pixscale) & !is.null(keyvalues)){
    pixscale = pixscale(keyvalues)
    if(verbose){message(paste('Extracted pixel scale from keyvalues provided:',signif(pixscale,4),'asec/pixel'))}
  }else{
    if(verbose){message(paste('Using suggested pixel scale:',signif(pixscale,4),'asec/pixel'))}
  }
  
  if(is.null(segim)){
    stop('Need segim!')
  }
  
  if(!is.null(mask)){
    image[mask > 0] = NA #means we will ignore the masked bits when doing the LL
  }
  
  fluxtype = tolower(fluxtype)
  
  if(fluxtype=='raw' | fluxtype=='adu' | fluxtype=='adus'){
    if(verbose){message('Using raw flux units')}
    fluxscale=1
  }else if (fluxtype=='jansky'){
    if(verbose){message('Using Jansky flux units (WARNING: magzero must take system to AB)')}
    fluxscale=10^(-0.4*(magzero - 8.9))
  }else if (fluxtype=='microjansky'){
    if(verbose){message('Using Micro-Jansky flux units (WARNING: magzero must take system to AB)')}
    fluxscale=10^(-0.4*(magzero - 23.9))
  }else{
    stop('fluxtype must be Jansky / Microjansky / Raw!')
  }
  
  segID = x = y = flux = j = rad2 = NULL
  
  Rapp = (app_diam / 2 / pixscale) #in pixels
  Aapp = (pi * Rapp^2) #in pixels
  
  
  if(is.null(tar)){
    segsel = which(segim > 0, arr.ind = TRUE)
  }else{
    tar = as.data.frame(tar)
    
    if(is.null(tar$xcen) & !is.null(tar$xmax)){
      stop('Fibre x position must be called xcen (not xmax)!')
    }
    
    if(is.null(tar$ycen) & !is.null(tar$ymax)){
      stop('Fibre y position must be called ycen (not ymax)!')
    }
    
    if(!is.null(tar$xcen)){
      # + 0.5 since arr.ind is offset compared to our normal image definition
      tar[,'xcen'] = tar[,'xcen'] + 0.5
    }
    
    if(!is.null(tar$ycen)){
      # + 0.5 since arr.ind is offset compared to our normal image definition
      tar[,'ycen'] = tar[,'ycen'] + 0.5
    }
    
    if(!is.null(tar$xcen) & !is.null(tar$ycen) & is.null(tar$segID)){
      tar$segID = segim[as.matrix(tar[,c('xcen', 'ycen')])]
    }
    
    if(any(tar$segID == 0L)){
      stop('Some fibres are not on a segment!')
    }
    
    if(anyDuplicated(tar$segID)){
      stop('Not all fibre positions fall on a unique segment!')
    }
    
    segsel = segim %in% tar$segID
    dim(segsel) = dim(segim)
    segsel = which(segsel, arr.ind = TRUE)
  }
  
  segID_all = as.integer(segim[segsel])
  
  tempDT = data.table(
    segID = segID_all,
    x = as.numeric(segsel[,1]),
    y = as.numeric(segsel[,2]),
    flux = as.numeric(image[segsel]),
    key = 'segID'
  )
  
  segID_tar = unique(tempDT$segID)
  
  if(is.null(tar$xcen) | is.null(tar$ycen)){
    if(centype == 'wt' | centype == 'mean'){
      if(verbose){message('Computing mean flux weighted aperture position')}
      xcen = tempDT[, .meanwt(x, flux), by=segID]$V1 #already has +0.5 pix offset
      ycen = tempDT[, .meanwt(y, flux), by=segID]$V1 #already has +0.5 pix offset
    }else if(centype == 'max'){
      if(verbose){message('Computing max flux aperture position')}
      xcen = tempDT[,x[which.max(flux)], by=segID]$V1 #already has +0.5 pix offset
      ycen = tempDT[,y[which.max(flux)], by=segID]$V1 #already has +0.5 pix offset
    }
    
    tar = data.frame(segID=segID_tar, xcen=xcen, ycen=ycen)
  }else{
    if(verbose){message('User provided tar table of fibre positions')}
  }
  
  if(verbose){message('Placing fibres and sorting')}
  
  match_segID = match(tempDT$segID, tar$segID)
  tempDT[, x:= x - tar[match_segID, 'xcen']]
  tempDT[, y:= y - tar[match_segID, 'ycen']]
  tempDT[, rad2:= x^2 + y^2]
  
  #setkey(tempDT, segID, rad2)
  
  if(verbose){message('Computing fibre photometry...')}
  
  output = foreach(j = seq_along(Rapp), .combine='cbind')%do%{
    #newer more accurate fibre mag calculation
    if(verbose){message('  Aperture: ', app_diam[j], ' [asec]')}
    # the top one can completely remove segments in some weird edge cases, and doesn't appear to be faster. Use the second!
    #temp_app = tempDT[rad2 <= Rapp[j]^2, .fluxcalcapp(x=x, y=y, rad2=rad2, flux=flux, xcen=0, ycen=0, rad_app=NA), by=segID]
    temp_app = tempDT[, .fluxcalcapp(rad2=rad2, flux=flux, xcen=0, ycen=0, rad_app=Rapp[j]), by=segID]
    
    if(correction){
      temp_app$flux_app = temp_app$flux_app - (temp_app$N - Aapp[j])*temp_app$flux_min
    }
    
    mag_app = profoundFlux2Mag(flux = temp_app$flux_app, magzero = magzero)
    return(data.frame(flux_app = temp_app$flux_app*fluxscale,
                      mag_app = mag_app,
                      SB_app = mag_app + 2.5*log10(Aapp[j]) + 5*log10(pixscale),
                      N_app = temp_app$N,
                      frac_app = temp_app$N/Aapp[j],
                      flux_min = temp_app$flux_min*fluxscale)
    )
  }
  
  colnames(output) = paste(colnames(output), rep(1:length(Rapp), each=6), sep='_')
  
  output = cbind(segID=segID_tar, xcen=tar$xcen - 0.5, ycen=tar$ycen - 0.5, output)
  if(verbose){message('Done!')}
  return(output)
}

profoundAperRan = function(image=NULL, segim=NULL, app_diam=1, mask=NULL, Nran=100, keyvalues=NULL,
                          pixscale=1, magzero=0, correction=TRUE, fluxtype='Raw', verbose=FALSE){
  if(!is.null(image)){
    if(inherits(image, 'Rfits_image')){
      keyvalues = image$keyvalues
      image = image$imDat
    }else if(inherits(image, 'Rfits_pointer')){
      keyvalues = image$keyvalues
      image = image[,]$imDat
    }else if(inherits(image, 'matrix')){
      'Do nothing'
    }else{
      stop('As of ProFound v1.21.0 only Rfits_image FITS inputs are allowed. Please install from GitHub asgr/Rfits')
    }
  }else{
    stop('Missing image - this is a required input!')
  }
  
  if(!is.null(keyvalues)){
    if(!inherits(keyvalues, 'Rfits_keylist')){
      if(is.list(keyvalues)){
        class(keyvalues) = 'Rfits_keylist'
      }else{
        stop('keyvalues is the wrong format- should be a list!')
      }
    }
  }
  
  if(missing(pixscale) & !is.null(keyvalues)){
    pixscale = pixscale(keyvalues)
    if(verbose){message(paste('Extracted pixel scale from keyvalues provided:',signif(pixscale,4),'asec/pixel'))}
  }else{
    if(verbose){message(paste('Using suggested pixel scale:',signif(pixscale,4),'asec/pixel'))}
  }
  
  if(is.null(segim)){
    stop('Need segim!')
  }
  
  if(!is.null(mask)){
    image[mask > 0] = NA #means we will ignore the masked bits when doing the LL
  }
  
  fluxtype = tolower(fluxtype)
  
  if(fluxtype=='raw' | fluxtype=='adu' | fluxtype=='adus'){
    if(verbose){message('Using raw flux units')}
    fluxscale=1
  }else if (fluxtype=='jansky'){
    if(verbose){message('Using Jansky flux units (WARNING: magzero must take system to AB)')}
    fluxscale=10^(-0.4*(magzero - 8.9))
  }else if (fluxtype=='microjansky'){
    if(verbose){message('Using Micro-Jansky flux units (WARNING: magzero must take system to AB)')}
    fluxscale=10^(-0.4*(magzero - 23.9))
  }else{
    stop('fluxtype must be Jansky / Microjansky / Raw!')
  }
  
  temp_arr = which(segim == 0L & !is.na(image), arr.ind=TRUE)
  temp_sel = sample(dim(temp_arr)[1], Nran)
  
  segim_ran = matrix(0L, dim(image)[1], dim(image)[2])
  
  tar = data.frame(segID = 1:Nran, temp_arr[temp_sel,] - 0.5)
  colnames(tar)[2:3] = c('xcen', 'ycen')
  
  segim_ran[temp_arr[temp_sel,]] = tar$segID
  
  size = ceiling(max(app_diam)/pixscale) + 1L
  if(size %% 2 == 0){
    size = size + 1L
  }
  
  segim_ran = profoundDilate(segim_ran, size = size)
  segim_ran[segim > 0L] = 0L
  
  output = profoundAperPhot(image=image, segim=segim_ran, app_diam=app_diam, keyvalues=keyvalues,
                            tar=tar, pixscale=pixscale, magzero=magzero, correction=correction, 
                            fluxtype=fluxtype, verbose=verbose)
  
  i = NULL
  
  errors = foreach(i = 1:length(app_diam), .combine='c')%do%{
    Nmax = max(output[,paste0('N_app_',i)], na.rm=TRUE)
    sel = which(output[,paste0('N_app_',i)] == Nmax)
    as.numeric(diff(quantile(output[sel,paste0('flux_app_',i)],c(0.16,0.84), na.rm=TRUE))/2)
  }
  
  return(list(AperPhot=output, errors=errors, segim_ran=segim_ran))
}
