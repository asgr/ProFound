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

profoundImBlur=function(image=NULL, sigma=1, plot=FALSE, ...){
  if(requireNamespace("imager", quietly = TRUE)){
    output=as.matrix(imager::isoblur(imager::as.cimg(image),sigma,na.rm=TRUE))
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
  invisible(output)
}

profoundImGrad=function(image=NULL, sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  output=as.matrix(imager::enorm(imager::imgradient(imager::isoblur(imager::as.cimg(image),sigma,na.rm=TRUE), "xy")))
  if(plot){
    magimage(output, ...)
  }
  invisible(output)
}

profoundImDiff=function(image=NULL,sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  blur=as.matrix(imager::isoblur(imager::as.cimg(image),sigma,na.rm=TRUE))
  output=image-blur
  if(plot){
    magimage(output, ...)
  }
  invisible(output)
}

profoundMakeSigma=function(image=NULL, objects=NULL, sky=0, skyRMS=0, readRMS=0, darkRMS=0, skycut=0, gain=1, image_units='ADU', sky_units='ADU', read_units='ADU', dark_units='ADU', output_units='ADU', plot=FALSE, ...){
  if(!is.null(objects)){
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
  invisible(sigma)
}

profoundGainEst=function(image=NULL, mask=0, objects=0, sky=0, skyRMS=1){
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
    invisible(value)
  }

  suppressWarnings({findgain=optim(par=startgain, fn=tempfunc, method="Brent", tempval=tempval, skyRMS=skyRMS, lower=startgain-2, upper=startgain+2)})
  invisible(10^findgain$par)
}

profoundCatMerge=function(segstats=NULL, groupstats=NULL, groupsegID=NULL, groupID_merge=NULL, flag=TRUE, rowreset=FALSE){
  if(! is.null(groupID_merge)){
    remove_segIDs=unique(unlist(groupsegID[groupsegID$groupID %in% groupID_merge,'segID']))
    remove_segIDs=remove_segIDs[!remove_segIDs %in% groupID_merge]
    segstats=segstats[! segstats$segID %in% remove_segIDs,]
    segstats[segstats$segID %in% groupID_merge,2:dim(segstats)[2]]=NA
    segstats[segstats$segID %in% groupID_merge,1:dim(groupstats)[2]]=groupstats[groupstats$groupID %in% groupID_merge,]
    segstats=segstats[order(segstats$segID),]
  }
  if(flag){
    segstats=cbind(segstats, origin='seg', stringsAsFactors=FALSE)
    segstats[segstats$segID %in% groupID_merge,'origin']='group'
  }
  if(rowreset){
    row.names(segstats)=NULL
  }
  invisible(segstats)
}

profoundResample=function(image, pixscale_old=1, pixscale_new=1, type='bicubic', fluxscale='image', recentre=FALSE){
  xseq=1:dim(image)[1]-dim(image)[1]/2-0.5
  yseq=1:dim(image)[2]-dim(image)[2]/2-0.5
  
  relscale=pixscale_new/pixscale_old
  
  xout=seq(relscale,xseq[length(xseq)],by=relscale)
  xout=c(-rev(xout),0,xout)
  yout=seq(relscale,yseq[length(yseq)],by=relscale)
  yout=c(-rev(yout),0,yout)
  bigrid=expand.grid(xout,yout)
  
  output=matrix(0,length(xout),length(yout))
  
  if(type=='bilinear'){
    output[]=.interp.2d(bigrid[,1], bigrid[,2], list(x=xseq, y=yseq, z=image))
  }else if(type=='bicubic'){
    if(!requireNamespace("akima", quietly = TRUE)){
      stop('The akima package is needed for bicubic interpolation to work. Please install it from CRAN.', call. = FALSE)
    }
    output[]=akima::bicubic(xseq, yseq, image,bigrid[,1], bigrid[,2])$z
  }else{
    stop('type must be one of bilinear / bicubic !')
  }
  
  if(recentre){
    maxloc=as.numeric(bigrid[which.max(output),])
    xout=xout+maxloc[1]
    yout=yout+maxloc[2]
    bigrid=expand.grid(xout,yout)
    
    if(type=='bilinear'){
      output[]=.interp.2d(bigrid[,1], bigrid[,2], list(x=xseq, y=yseq, z=image))
    }else{
      output[]=akima::bicubic(xseq, yseq, image,bigrid[,1], bigrid[,2])$z
    }
  }
  
  if(fluxscale=='image'){
    output=output*sum(image)/sum(output)
  }else if(fluxscale=='pixscale'){
    output=output*relscale^2
  }else if(fluxscale=='norm'){
    output=output/sum(output)
  }else{
    stop('fluxscale must be one of image / pixscale / norm !')
  }
  
  return(invisible(output))
}

# Hidden utility functions

.interp.2d=function(x, y, obj){
    xobj = obj$x
    yobj = obj$y
    zobj = obj$z
    nx = length(xobj)
    ny = length(yobj)
    lx = approx(xobj, 1:nx, x, rule = 2)$y
    rm(x)
    ly = approx(yobj, 1:ny, y, rule = 2)$y
    rm(y)
    lx1 = floor(lx)
    ly1 = floor(ly)
    ex = lx - lx1
    rm(lx)
    ey = ly - ly1
    rm(ly)
    ex[lx1 == nx] = 1
    ey[ly1 == ny] = 1
    lx1[lx1 == nx] = nx - 1
    ly1[ly1 == ny] = ny - 1
    temp=rep(0,length(lx1))
    temp=zobj[cbind(lx1, ly1)] * (1 - ex) * (1 - ey)
    temp=temp+zobj[cbind(lx1 + 1, ly1)] * ex * (1 - ey)
    temp=temp+zobj[cbind(lx1, ly1 + 1)] * (1 - ex) * ey
    temp=temp+zobj[cbind(lx1 + 1, ly1 + 1)] * ex * ey
    invisible(temp)
}
