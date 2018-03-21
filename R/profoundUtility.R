profoundMag2Mu=function(mag=15, re=1, axrat=1, pixscale=1){
  return(mag+2.5*log10(pi*re^2*axrat)-2.5*log10(0.5)+5*log10(pixscale))
}

profoundMu2Mag=function(mu=17, re=1, axrat=1, pixscale=1){
  return(mu-2.5*log10(pi*re^2*axrat)+2.5*log10(0.5)-5*log10(pixscale))
}

profoundGainConvert=function(gain=1, magzero=0, magzero_new=0){
  return(gain*10^(-0.4*(magzero_new-magzero)))
}

profoundMag2Flux=function(mag=0, magzero=0){
  return(10^(-0.4*(mag-magzero)))
}

profoundFlux2Mag=function(flux=1, magzero=0){
  flux[is.na(flux)]=0
  output=flux
  output[]=NA
  output[which(flux>0)]=-2.5*log10(flux[which(flux>0)])+magzero
  return(output)
}

profoundFlux2SB=function(flux=1, magzero=0, pixscale=1){
  return(profoundFlux2Mag(flux=flux, magzero=magzero)+5*log10(pixscale))
}

profoundSB2Flux=function(SB=0, magzero=0, pixscale=1){
  mag=SB-5*log10(pixscale)
  return(profoundMag2Flux(mag=mag, magzero=magzero))
}

profoundImBlur=function(image, sigma=1, plot=FALSE, ...){
  if(requireNamespace("imager", quietly = TRUE)){
    output=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  }else{
    if(!requireNamespace("EBImage", quietly = TRUE)){
      stop('The imager or EBImage package is needed for the profoundImBlur function to work. Please install from CRAN.', call. = FALSE)
    }
    message(" - WARNING: imager package not installed, using EBImage gblur smoothing!")
    output=as.matrix(EBImage::gblur(image,sigma))
  }
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profoundImGrad=function(image, sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  output=as.matrix(imager::enorm(imager::imgradient(imager::isoblur(imager::as.cimg(image),sigma), "xy")))
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profoundImDiff=function(image,sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  blur=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  output=image-blur
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profoundMakeSigma=function(image, objects, sky=0, skyRMS=0, readRMS=0, darkRMS=0, skycut=0, gain=1, image_units='ADU', sky_units='ADU', read_units='ADU', dark_units='ADU', output_units='ADU', plot=FALSE, ...){
  if(!missing(objects)){
    if(length(objects)==length(image)){
      image[objects==0]=0
    }
  }
  if(image_units=='ADU'){
    image=gain*image
  }else if(image_units=='elec'){
    NULL
  }else{
    stop(paste('image_units unit type of',image_units,'not recognised, must be ADU or elec'))
  }
  
  if(sky_units=='ADU'){
    sky=gain*sky
    skyRMS=gain*skyRMS
  }else if(sky_units=='elec'){
    NULL
  }else{
    stop(paste('sky_units unit type of',sky_units,'not recognised, must be ADU or elec'))
  }
  
  if(read_units=='ADU'){
    readRMS=gain*readRMS
  }else if(read_units=='elec'){
    NULL
  }else{
    stop(paste('read_units unit type of',read_units,'not recognised, must be ADU or elec'))
  }
  
  if(dark_units=='ADU'){
    darkRMS=gain*darkRMS
  }else if(dark_units=='elec'){
    NULL
  }else{
    stop(paste('dark_units unit type of',dark_units,'not recognised, must be ADU or elec'))
  }
  
  image=image-sky
  image[image < skyRMS*skycut]=0
  
  if(output_units=='ADU'){
    sigma=sqrt(image+skyRMS^2+readRMS^2+darkRMS^2)/gain
  }else if(output_units=='elec'){
    sigma=sqrt(image+skyRMS^2+readRMS^2+darkRMS^2)
  }else{
    stop(paste('output_units unit type of',output_units,'not recognised, must be ADU or elec'))
  }
  
  if(plot){
    magimage(sigma, ...)
  }
  return=sigma
}

profoundGainEst=function(image, mask=0, objects=0, sky=0, skyRMS=1){
  if(missing(sky)){
    sky=profoundSkyEst(image=image, mask=mask, objects=objects,plot=FALSE)$sky
  }
  image_sky=image-sky
  if(missing(skyRMS)){
    skyRMS=profoundSkyEst(image=profoundImDiff(image_sky,3), mask=mask, objects=objects,plot=FALSE)$skyRMS
  }
  tempval=as.numeric(image_sky[mask==0 & objects==0])
  
  startgain=ceiling(log10(abs(min(tempval, na.rm=T)-sky)/(skyRMS^2)))+1
  
  tempfunc=function(gain,tempval,skyRMS){
    gain=10^gain
    floor=(skyRMS*gain)^2
    trialdata=tempval*gain+floor
    value=-sum(dpois(x=round(trialdata), lambda=floor, log=T))
    return=value
  }

  suppressWarnings({findgain=optim(par=startgain, fn=tempfunc, method="Brent", tempval=tempval, skyRMS=skyRMS, lower=startgain-2, upper=startgain+2)})
  return=10^findgain$par
}

### Deprecated Functions ###

# profoundGetPixScale=function(header, CD1_1=1, CD1_2=0, CD2_1=0, CD2_2=1){
#   if(!missing(header)){
#     if(is.data.frame(header) | is.matrix(header)){
#       locs=match(c('CD1_1','CD1_2','CD2_1','CD2_2'),header[,1])
#       headerWCS=data.frame(header[locs,1],as.numeric(header[locs,2]))
#       if('CD1_1' %in% headerWCS[,1]){
#         CD1_1=headerWCS[headerWCS[,1]=='CD1_1',2]
#         if('CD1_2' %in% headerWCS[,1]){CD1_2=headerWCS[headerWCS[,1]=='CD1_2',2]}else{message('Missing CD1_2')}
#       }else{
#         if('CDELT1' %in% headerWCS[,1]){
#           CD1_1=headerWCS[headerWCS[,1]=='CDELT1',2]
#         }else{
#           message("Missing CD1_1 and CDELT1")
#         }
#       }
#       if('CD2_2' %in% headerWCS[,1]){
#         CD2_2=headerWCS[headerWCS[,1]=='CD2_2',2]
#         if('CD2_1' %in% headerWCS[,1]){CD2_1=headerWCS[headerWCS[,1]=='CD2_1',2]}else{message('Missing CD2_1')}
#       }else{
#         if('CDELT2' %in% headerWCS[,1]){
#           CD2_2=headerWCS[headerWCS[,1]=='CDELT2',2]
#         }else{
#           message("Missing CD2_2 and CDELT2")
#         }
#       }
#     }else{
#       if('CD1_1' %in% header){
#         CD1_1=as.numeric(header[which(header=='CD1_1')+1])
#         if('CD1_2' %in% header){CD1_2=as.numeric(header[which(header=='CD1_2')+1])}else{message('Missing CD1_2')}
#       }else{
#         if('CDELT1' %in% header){
#           CD1_1=as.numeric(header[which(header=='CDELT1')+1])
#         }else{
#           message("Missing CD1_1 and CDELT1")
#         }
#       }
#       if('CD2_2' %in% header){
#         CD2_2=as.numeric(header[which(header=='CD2_2')+1])
#         if('CD2_1' %in% header){CD2_1=as.numeric(header[which(header=='CD2_1')+1])}else{message('Missing CD2_1')}
#       }else{
#         if('CDELT1' %in% header){
#           CD2_2=as.numeric(header[which(header=='CDELT2')+1])
#         }else{
#           message("Missing CD2_2 and CDELT2")
#         }
#       }
#     }
#   }
#   return(3600*(sqrt(CD1_1^2+CD1_2^2)+sqrt(CD2_1^2+CD2_2^2))/2)
# }

# profoundInterp2d=function(x,y,image){
#     scale=sum(image)
#     imagelist=list(x=seq(-dim(image)[1]/2,dim(image)[1]/2,len=dim(image)[1]),y=seq(-dim(image)[2]/2,dim(image)[2]/2,len=dim(image)[2]),z=image)
#     ximage = seq(-dim(image)[1]/2,dim(image)[1]/2,len=dim(image)[1])
#     yimage = seq(-dim(image)[2]/2,dim(image)[2]/2,len=dim(image)[2])
#     zimage = image
#     nx = length(ximage)
#     ny = length(yimage)
#     lx = approx(ximage, 1:nx, x, rule=2)$y
#     ly = approx(yimage, 1:ny, y, rule=2)$y
#     lx1 = floor(lx)
#     ly1 = floor(ly)
#     ex = lx - lx1
#     ey = ly - ly1
#     ex[lx1 == nx] = 1
#     ey[ly1 == ny] = 1
#     lx1[lx1 == nx] = nx - 1
#     ly1[ly1 == ny] = ny - 1
#     z=
# 	zimage[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) +
# 	zimage[cbind(lx1 + 1, ly1)] * ex * (1 - ey) +
# 	zimage[cbind(lx1, ly1 + 1)] * (1 - ex) * ey +
# 	zimage[cbind(lx1 + 1, ly1 + 1)] * ex * ey
#   return = cbind(X=x,Y=y,Z=z)
# }
