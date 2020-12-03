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
  #xseq = 1:dim(image)[1]-dim(image)[1]/2-0.5
  #yseq = 1:dim(image)[2]-dim(image)[2]/2-0.5
  
  relscale = pixscale_new/pixscale_old
  
  # xout = seq(relscale,xseq[length(xseq)],by=relscale)
  # xout = c(-rev(xout),0,xout)
  # yout = seq(relscale,yseq[length(yseq)],by=relscale)
  # yout = c(-rev(yout),0,yout)
  # bigrid = expand.grid(xout,yout)
  xin = (1:(dim(image)[1]/2)) * relscale
  xin = c(-rev(xin),0,xin) + dim(image)[1] * relscale/2

  yin = (1:(dim(image)[2]/2)) * relscale
  yin = c(-rev(yin),0,yin) + dim(image)[2] * relscale/2
  
  output = matrix(0,dim(image)[1] * relscale, dim(image)[2] * relscale)
  
  if(type=='bilinear'){
    .interpolateLinearGrid(xin, yin, image, output)
  }else if(type=='bicubic'){
    .interpolateAkimaGrid(xin, yin, image, output)
  }else{
    stop('type must be one of bilinear / bicubic !')
  }
  
  # if(type=='bilinear'){
  #   output[]=.interp.2d(bigrid[,1], bigrid[,2], list(x=xseq, y=yseq, z=image))
  # }else if(type=='bicubic'){
  #   if(!requireNamespace("akima", quietly = TRUE)){
  #     stop('The akima package is needed for bicubic interpolation to work. Please install it from CRAN.', call. = FALSE)
  #   }
  #   output[]=akima::bicubic(xseq, yseq, image,bigrid[,1], bigrid[,2])$z
  # }else{
  #   stop('type must be one of bilinear / bicubic !')
  # }
  
  if(recentre){
    bigrid = expand.grid(1:dim(output)[1]-0.5,1:dim(output)[2]-0.5)
    maxloc = as.numeric(bigrid[which.max(output),])
    xin = xin - maxloc[1] + dim(output)[1]/2
    yin = yin - maxloc[2] + dim(output)[2]/2
    #bigrid=expand.grid(xout,yout)
    
    if(type=='bilinear'){
      .interpolateLinearGrid(xin, xin, image, output)
    }else if(type=='bicubic'){
      .interpolateAkimaGrid(yin, yin, image, output)
    }
    
    # if(type=='bilinear'){
    #   output[]=.interp.2d(bigrid[,1], bigrid[,2], list(x=xseq, y=yseq, z=image))
    # }else{
    #   output[]=akima::bicubic(xseq, yseq, image,bigrid[,1], bigrid[,2])$z
    # }
  }
  
  if(fluxscale == 'image'){
    output = output*sum(image)/sum(output)
  }else if(fluxscale == 'pixscale'){
    output = output*relscale^2
  }else if(fluxscale == 'norm'){
    output = output/sum(output)
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

.genPointSource = function(xcen=50, ycen=50, flux=1, psf, dim=c(100,100)){
  dimpsf = dim(psf)
  psfcen=dim(psf)/2
  
  image=matrix(0,dim[1],dim[2])
  
  for(i in 1:length(xcen)){
    x_off = (xcen[i] - psfcen[1])
    pixshift_x = floor(x_off)
    x_off = 1L - (x_off - pixshift_x)
    
    y_off = (ycen[i] - psfcen[2])
    pixshift_y = floor(y_off)
    y_off = 1L - (y_off - pixshift_y)
    
    left_bottom = x_off * y_off
    left_top = x_off * (1-y_off)
    right_bottom = (1-x_off) * y_off
    right_top = (1-x_off) * (1-y_off)
    
    xpix_left = 1:dimpsf[1] + pixshift_x
    inim_xpix_left = xpix_left >= 1 & xpix_left <= dim[1]
    xpix_right = xpix_left + 1L
    inim_xpix_right = xpix_right >= 1 & xpix_right <= dim[1]
    
    ypix_bottom = 1:dimpsf[2] + pixshift_y
    inim_ypix_bottom = ypix_bottom >= 1 & ypix_bottom <= dim[2]
    ypix_top = ypix_bottom + 1L
    inim_ypix_top = ypix_top >= 1 & ypix_top <= dim[2]
    
    dobottom = any(inim_ypix_bottom)
    doleft = any(inim_xpix_left)
    dotop = any(inim_ypix_top)
    doright = any(inim_xpix_right)
    
    if(dobottom & doleft & left_bottom>0){
      image = .addmat_cpp(
        image,
        psf[(1:dimpsf[1])[inim_xpix_left], (1:dimpsf[2])[inim_ypix_bottom]] * left_bottom * flux[i],
        range(xpix_left[inim_xpix_left]),
        range(ypix_bottom[inim_ypix_bottom])
      )
      #image[xpix_left[inim_xpix_left],ypix_bottom[inim_ypix_bottom]] = 
        #image[xpix_left[inim_xpix_left],ypix_bottom[inim_ypix_bottom]] + psf[(1:dimpsf[1])[inim_xpix_left], (1:dimpsf[2])[inim_ypix_bottom]] * left_bottom * flux[i]
    }
    
    if(dotop & doleft & left_top>0){
      image = .addmat_cpp(
        image,
        psf[(1:dimpsf[1])[inim_xpix_left], (1:dimpsf[2])[inim_ypix_top]] * left_top * flux[i],
        range(xpix_left[inim_xpix_left]),
        range(ypix_top[inim_ypix_top])
      )
      #image[xpix_left[inim_xpix_left],ypix_top[inim_ypix_top]] =
        #image[xpix_left[inim_xpix_left],ypix_top[inim_ypix_top]] + psf[(1:dimpsf[1])[inim_xpix_left], (1:dimpsf[2])[inim_ypix_top]] * left_top * flux[i]
    }
    
    if(dobottom & doright & right_bottom>0){
      image = .addmat_cpp(
        image,
        psf[(1:dimpsf[1])[inim_xpix_right], (1:dimpsf[2])[inim_ypix_bottom]] * right_bottom * flux[i],
        range(xpix_right[inim_xpix_right]),
        range(ypix_bottom[inim_ypix_bottom])
      )
      #image[xpix_right[inim_xpix_right],ypix_bottom[inim_ypix_bottom]] = 
        #image[xpix_right[inim_xpix_right],ypix_bottom[inim_ypix_bottom]] + psf[(1:dimpsf[1])[inim_xpix_right], (1:dimpsf[2])[inim_ypix_bottom]] * right_bottom * flux[i]
    }
    
    if(dotop & doright & right_top>0){
      image = .addmat_cpp(
        image,
        psf[(1:dimpsf[1])[inim_xpix_right], (1:dimpsf[2])[inim_ypix_top]] * right_top * flux[i],
        range(xpix_right[inim_xpix_right]),
        range(ypix_top[inim_ypix_top])
      )
      #image[xpix_right[inim_xpix_right],ypix_top[inim_ypix_top]] =
        #image[xpix_right[inim_xpix_right],ypix_top[inim_ypix_top]] + psf[(1:dimpsf[1])[inim_xpix_right], (1:dimpsf[2])[inim_ypix_top]] * right_top * flux[i]
    }
  }
  return(invisible(image))
}

.makeBrush = function(size, shape=c('box', 'disc', 'diamond', 'Gaussian', 'line'), step=TRUE, sigma=0.3, angle=45) {
  #This is a direct port of EBImage::makeBrush. This reduces code dependencies, and EBImage does not appear to be well maintained.
  if(! (is.numeric(size) && (length(size)==1L) && (size>=1)) ) stop("'size' must be an odd integer.")
  shape = match.arg(arg = tolower(shape), choices = c('box', 'disc', 'diamond', 'gaussian', 'line'))
  
  if(size %% 2 == 0){
    size = size + 1
    warning(paste("'size' was rounded to the next odd number: ", size))
  }
  
  if (shape=='box') z = matrix(1L, size, size)
  else if (shape == 'line') {
    angle = angle %% 180
    angle.radians = angle * pi / 180;
    tg = tan(angle.radians)
    sizeh = (size-1)/2
    if ( angle < 45 || angle > 135) {
      z.x = sizeh
      z.y = round(sizeh*tg)
    }
    else {
      z.y = sizeh
      z.x = round(sizeh/tg)
    }
    z = array(0L, dim=2*c(z.x, z.y)+1);
    for (i in -sizeh:sizeh) {
      if ( angle < 45 || angle > 135) {
        ## scan horizontally
        i.x = i
        i.y = round(i*tg)
      }
      else {
        ## scan vertically
        i.y = i
        i.x = round(i/tg) 
      }
      z[i.x+z.x+1, i.y+z.y+1] = 1L
    }
  }
  else if (shape=='gaussian') {
    x = seq(-(size-1)/2, (size-1)/2, length=size)
    x = matrix(x, size, size)
    z = exp(- (x^2 + t(x)^2) / (2*sigma^2))
    z = z / sum(z)
  } else {
    ## pixel center coordinates
    x = 1:size -((size+1)/2)
    
    ## for each pixel, compute the distance from its center to the origin, using L1 norm ('diamond') or L2 norm ('disc')
    if (shape=='disc') {
      z = outer(x, x, FUN=function(X,Y) (X*X+Y*Y))
      mz = (size/2)^2
      z = (mz - z)/mz
      z = sqrt(ifelse(z>0, z, 0))
    } else {
      z = outer(x, x, FUN=function(X,Y) (abs(X)+abs(Y)))
      mz = (size/2)
      z = (mz - z)/mz
      z = ifelse(z>0, z, 0)
    }
    
    if (step) z = ifelse(z>0, 1L, 0L)
  }
  z
}
