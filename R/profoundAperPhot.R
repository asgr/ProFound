.fluxcalcapp = function(x=NULL, y=NULL, flux=NULL, xcen=NA, ycen=NA,
                        rad_app=NULL, centype='max'){
  if(is.na(xcen)){
    if(centype == 'wt' | centype == 'mean'){
      xcen = .meanwt(x, flux)
    }else if(centype == 'max'){
      xcen = x[which.max(flux)]
    }
  }
  
  if(is.na(ycen)){
    if(centype == 'wt' | centype == 'mean'){
      ycen = .meanwt(y, flux)
    }else if(centype == 'max'){
      ycen = y[which.max(flux)]
    }
  }
  
  if(xcen != 0 & ycen != 0){
    rad2 = (x - xcen)^2 + (y - ycen)^2
  }
  
  sel = which(rad2 < rad_app^2)
  rad_out = max(rad2[sel])
  sel_out = which(rad2 == rad_out)
  
  flux_app=sum(flux[sel], na.rm=TRUE)
  flux_min = mean(flux[sel_out])
  flux_min[flux_min < 0] = 0 #don't want to penalise when in the sky noise
  
  return(list(flux_app=flux_app, flux_min=flux_min, N=length(sel)))
}

profoundAperPhot = function(image=NULL, segim=NULL, app_diam=1, keyvalues=NULL, tar=NULL,
                           pixscale=1, magzero=0, correction=TRUE, centype='max',
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
  
  segID = x = y = flux = j = NULL
  
  Rapp = (app_diam / 2 / pixscale)
  Aapp = (pi * Rapp^2)
  
  if(is.null(tar)){
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
  
  tempDT = data.table(
    segID = as.integer(segim[segsel]),
    x = as.numeric(segsel[,1]),
    y = as.numeric(segsel[,2]),
    flux = as.numeric(image[segsel]),
    key = 'segID'
  )
  
  segID_out = unique(tempDT$segID)
  
  if(!is.null(tar$xcen) & !is.null(tar$ycen)){
    for(i in 1:dim(tar)[1]){
      if(!is.na(tar[i,'xcen']) & !is.na(tar[i,'ycen'])){
        tempDT[segID == tar[i,'segID'], x := x - tar[i,'xcen']]
        tempDT[segID == tar[i,'segID'], y := y - tar[i,'ycen']]
      }
    }
    
    output = foreach(j = seq_along(Rapp), .combine='cbind')%do%{
      #newer more accurate fibre mag calculation
      temp_app = tempDT[, .fluxcalcapp(x=x, y=y, flux=flux, xcen=0, ycen=0, rad_app=Rapp[j]), by=segID]
      if(correction){
        temp_app$flux_app = temp_app$flux_app - (temp_app$N - Aapp[j])*temp_app$flux_min
      }
      mag_app = profoundFlux2Mag(flux = temp_app$flux_app, magzero = magzero)
      return(data.frame(flux_app = temp_app$flux_app,
                        mag_app = mag_app,
                        N_app = temp_app$N,
                        frac_app = temp_app$N/Aapp[j],
                        flux_min = temp_app$flux_min)
      )
    }
  }else{
    output = foreach(j = seq_along(Rapp), .combine='cbind')%do%{
      #newer more accurate fibre mag calculation
      temp_app = tempDT[, .fluxcalcapp(x=x, y=y, flux=flux, rad_app=Rapp[j],
                                       centype=centype), by=segID]
      #Here we correct by the lowest value pixel in the outer aperture
      if(correction){
        temp_app$flux_app = temp_app$flux_app - (temp_app$N - Aapp[j])*temp_app$flux_min
      }
      mag_app = profoundFlux2Mag(flux = temp_app$flux_app, magzero = magzero)
      return(data.frame(flux_app = temp_app$flux_app,
                        mag_app = mag_app,
                        N_app = temp_app$N,
                        frac_app = temp_app$N/Aapp[j],
                        flux_min = temp_app$flux_min)
             )
    }
  }
  
  colnames(output) = paste(colnames(output), rep(1:length(Rapp), each=5), sep='_')
  
  output = cbind(segID=segID_out, output)
  
  return(output)
}