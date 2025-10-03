.fluxcalcapp = function(x=NULL, y=NULL, flux=NULL, xcen=NA, ycen=NA, rad_app=NULL, centype='max'){
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
  
  sel = which(rad2 < rad_app^2)
  rad_out = max(rad2[sel])
  sel_out = which(rad2 == rad_out)
  
  flux_app = sum(flux[sel], na.rm=TRUE)
  flux_min = mean(flux[sel_out], na.rm=TRUE)
  
  if(!isTRUE(is.finite(flux_min))){ #this catches NA, NaN, NULL, Inf events
    flux_min = 0 #don't want to penalise when masked or other weird events
  }else if(flux_min < 0){
    flux_min = 0 #don't want to penalise when in the sky noise 
  }
  
  return(list(flux_app=flux_app, flux_min=flux_min, N=length(sel)))
}

profoundAperPhot = function(image=NULL, segim=NULL, app_diam=1, keyvalues=NULL, tar=NULL,
                           pixscale=1, magzero=0, correction=TRUE, centype='max', fluxtype='Raw',
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
  
  segID = x = y = flux = j = NULL
  
  Rapp = (app_diam / 2 / pixscale)
  Aapp = (pi * Rapp^2)
  
  if(is.null(tar$segID)){
    segsel = which(segim > 0, arr.ind = TRUE)
  }else{
    tar = as.data.frame(tar)
    
    segsel = segim %in% tar$segID
    dim(segsel) = dim(segim)
    segsel = which(segsel, arr.ind = TRUE)
    
    if(!is.null(tar$xcen)){
      # + 0.5 since arr.ind is offset compared to our normal image definition
      tar[,'xcen'] = tar[,'xcen'] + 0.5
    }
    
    if(!is.null(tar$ycen)){
      # + 0.5 since arr.ind is offset compared to our normal image definition
      tar[,'ycen'] = tar[,'ycen'] + 0.5
    }
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
  
  # for(i in 1:dim(tar)[1]){
  #   tempDT[segID == tar[i,'segID'], x := x - tar[i,'xcen']]
  #   tempDT[segID == tar[i,'segID'], y := y - tar[i,'ycen']]
  # }
  
  match_segID = match(tempDT$segID, tar$segID)
  tempDT[, x:= x - tar[match_segID, 'xcen']]
  tempDT[, y:= y - tar[match_segID, 'ycen']]
  
  if(verbose){message('Computing fibre photometry...')}
  
  output = foreach(j = seq_along(Rapp), .combine='cbind')%do%{
    #newer more accurate fibre mag calculation
    if(verbose){message('  Aperture: ', app_diam[j], ' / asec')}
    temp_app = tempDT[, .fluxcalcapp(x=x, y=y, flux=flux, xcen=0, ycen=0, rad_app=Rapp[j]), by=segID]
    
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