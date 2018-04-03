.meanwt=function(x, wt){
  wt[wt<0]=0
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  sum(x*wt, na.rm=TRUE)/sum(wt, na.rm=TRUE)
}

.varwt=function(x, wt, xcen){
  wt[wt<0]=0
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  return=(sum((x-xcen)^2*wt, na.rm=TRUE)/sum(wt, na.rm=TRUE))
}

.covarwt=function(x, y, wt, xcen, ycen){
  wt[wt<0]=0
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  if(missing(ycen)){ycen=.meanwt(y, wt)}
  return=(sum((x-xcen)*(y-ycen)*wt, na.rm=TRUE)/sum(wt, na.rm=TRUE))
}

.cov2eigval=function(sx,sy,sxy){
  b=-sx^2-sy^2
  c=sx^2*sy^2-sxy^2
  return=(list(hi=(-b+sqrt(b^2-4*c))/2,lo=(-b-sqrt(b^2-4*c))/2))
}

.cov2eigvec=function(sx,sy,sxy){
  eigval=.cov2eigval(sx,sy,sxy)$hi
  eigvec=(sx^2-eigval)/sxy
  return=eigvec
}

.eigvec2ang=function(eigvec){
  return=(90-atan(eigvec)*180/pi)%%180
}

.asymm=function(x, y, wt, xcen, ycen){
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  if(missing(ycen)){ycen=.meanwt(y, wt)}
  relx=round(x-xcen)
  rely=round(y-ycen)
  frame1=data.frame(x=relx,y=rely,wt1=wt)
  frame2=data.frame(x=-relx,y=-rely,wt2=wt)
  comp=merge(frame1,frame2,by=c('x','y'), all=TRUE)
  overlap=which(comp$wt1>0 & comp$wt2>0)
  asymm=sum(abs(comp[overlap,'wt1']-comp[overlap,'wt2']), na.rm=TRUE)/sum(abs(comp[overlap,'wt1']+comp[overlap,'wt2']), na.rm=TRUE)
  return=asymm
}

.reflect=function(x, y, wt, xcen, ycen){
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  if(missing(ycen)){ycen=.meanwt(y, wt)}
  relx=round(x-xcen)
  rely=round(y-ycen)
  frame1=data.frame(x=relx,y=rely,wt1=wt)
  frame2=data.frame(x=-relx,y=-rely,wt2=wt)
  comp=merge(frame1,frame2,by=c('x','y'),all=TRUE)
  overlap=is.na(comp$wt1)==FALSE & is.na(comp$wt2)==FALSE
  asymm=2*sum(abs(comp[overlap,'wt1']-comp[overlap,'wt2']), na.rm=TRUE)/sum(comp[overlap,'wt1'], comp[overlap,'wt2'], na.rm=TRUE)
  len=dim(frame1)[1]
  flux_reflect=sum(comp[overlap,'wt1'], na.rm=TRUE)+2*sum(frame1[is.na(comp$wt2),'wt1'], na.rm=TRUE)
  return=flux_reflect
}

#Not currently used:
.nser2ccon=function(nser=0.5, lo=0.5, hi=0.9){
  return=(((qgamma(lo, 2 * nser)/qgamma(hi, 2 * nser))^nser)^2)
}

#Not currently used (too slow):
.match2col=function(tab1, tab2){
  return(which( outer(tab1[,1], tab2[,1], "==") & outer(tab1[,2], tab2[,2], "=="), arr.ind=TRUE))
}

.fluxcalc=function(flux){
  
  good=which(!is.na(flux))
  N100seg=length(good)
  
  if(N100seg>0){
    
    sumflux=sum(flux[good], na.rm=TRUE)
    temp=cumsum(flux[good])/sumflux
    
    # print(sumflux)
    # print(flux)
    # print(temp)
  
    loc50=min(which(temp>=0.5))
    loc50cumsumhi=temp[loc50]
    loc50cumsumlo=temp[loc50-1]
    N50seg=N100seg-(loc50-1)+(0.5-loc50cumsumlo)/(loc50cumsumhi-loc50cumsumlo)
    
    loc90=min(which(temp>=0.1))
    loc90cumsumhi=temp[loc90]
    loc90cumsumlo=temp[loc90-1]
    N90seg=N100seg-(loc90-1)+(0.1-loc90cumsumlo)/(loc90cumsumhi-loc90cumsumlo)
  
    cenfrac=temp[N100seg]-temp[N100seg-1]
    
  }else{
    sumflux=0
    N50seg=0
    N90seg=0
    cenfrac=0
  }
  
  #N50seg=sum(temp>=0.5)
  #N90seg=sum(temp>=0.1)
  
  return=list(flux=sumflux, N50seg=N50seg, N90seg=N90seg, N100seg=N100seg, cenfrac=cenfrac)
}

profoundMakeSegim=function(image, mask, objects, skycut=1, pixcut=3, tolerance=4, ext=2, sigma=1, smooth=TRUE, SBlim, magzero=0, gain=NULL, pixscale=1, sky, skyRMS, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){
  if(length(image)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  call=match.call()
  if(verbose){message(' - Running MakeSegim:')}
  timestart = proc.time()[3]
  # if(!requireNamespace("imager", quietly = TRUE)){
  #   stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  # }
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  #Treat image NAs as masked regions:
  
  if(!missing(mask)){
    mask[is.na(image)]=1
  }else{
    if(any(is.na(image))){
      mask=matrix(0,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1
    }
  }
  
  if(missing(pixscale) & !missing(header)){
    pixscale=getpixscale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  image_orig=image
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profoundSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky - User provided sky")}
  }
  
  image_sky=image-sky
  
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profoundSkyEst(image=profoundImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  
  image=image_sky/skyRMS
  image[!is.finite(image)]=0
 
  if(smooth){
    if(verbose){message(paste(" - Smoothing the image -", round(proc.time()[3]-timestart,3), "sec"))}
    if(requireNamespace("imager", quietly = TRUE)){
      image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
    }else{
      if(!requireNamespace("EBImage", quietly = TRUE)){
      stop('The imager or EBImage package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
      }
      message(" - WARNING: imager package not installed, using EBImage gblur smoothing!")
      image=as.matrix(EBImage::gblur(image,sigma))
    }
  }else{
    if(verbose){message(" - Skipping smoothing - smooth set to FALSE")}
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!missing(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profoundSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
  if(!missing(mask)){
    image[mask!=0]=0
  }
  if(verbose){message(paste(" - Watershed de-blending -", round(proc.time()[3]-timestart,3), "sec"))}
  if(any(image>0)){
    segim=EBImage::imageData(EBImage::watershed(image,tolerance=tolerance,ext=ext))
    segim=matrix(as.integer(segim),xlen,ylen)
  }else{
    segim=image
  }
  
  segtab=tabulate(segim)
  segim[segim %in% which(segtab<pixcut)]=0

  if(plot){
    if(verbose){message(paste(" - Plotting segments -", round(proc.time()[3]-timestart,3), "sec"))}
    profoundSegimPlot(image=image_orig, segim=segim, mask=mask, sky=sky, ...)
  }else{
    if(verbose){message(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  objects=segim
  objects[objects!=0]=1
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profoundSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making final local estimate of the sky - User provided sky")}
  }
  
  image_sky=image_orig-sky
  
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profoundSkyEst(image=profoundImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  
  if(stats & any(image>0)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profoundSegimStats(image=image_orig, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE or no segments")}
    segstats=NULL
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  
  if(missing(header)){header=NULL}
  
  if(verbose){message(paste(" - MakeSegim is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return=list(segim=segim, objects=objects, sky=sky, skyRMS=skyRMS, segstats=segstats, header=header, SBlim=SBlim, call=call)
}

profoundMakeSegimExpand=function(image, segim, mask, objects, skycut=1, SBlim, magzero=0, gain=NULL, pixscale=1, sigma=1, smooth=TRUE, expandsigma=5, expand='all', sky, skyRMS, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){
  
  if(verbose){message(' - Running MakeSegimExpand:')}
  timestart = proc.time()[3]
  
  if(length(image)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  
  call=match.call()
  
  # if(!requireNamespace("imager", quietly = TRUE)){
  #   stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  # }
  
  #Treat image NAs as masked regions:
  
  if(!missing(mask)){
    mask[is.na(image)]=1
  }else{
    if(any(is.na(image))){
      mask=matrix(0,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1
    }
  }
  
  if(missing(pixscale) & !missing(header)){
    pixscale=getpixscale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  image_orig=image
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profoundSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky - User provided sky")}
  }
  image_sky=image-sky
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profoundSkyEst(image=profoundImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  image=image_sky/skyRMS
  image[!is.finite(image)]=0

  if(smooth){
    if(verbose){message(paste(" - Smoothing the image -", round(proc.time()[3]-timestart,3), "sec"))}
    if(requireNamespace("imager", quietly = TRUE)){
      image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
    }else{
      if(!requireNamespace("EBImage", quietly = TRUE)){
      stop('The imager or EBImage package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
      }
      message(" - WARNING: imager package not installed, using EBImage gblur smoothing!")
      image=as.matrix(EBImage::gblur(image,sigma))
    }
  }else{
    if(verbose){message(" - Skipping smoothing - smooth set to FALSE")}
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!missing(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profoundSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
  if(!missing(mask)){
    image[mask!=0]=0
  }
  #kernel=profoundMakeGaussianPSF(fwhm = expandsigma, dim=dim)
  maxmat=matrix(min(image, na.rm=TRUE), xlen, ylen)
  segim_new=segim
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  if(expand[1]=='all'){expand=segvec}
  if(verbose){message(paste(" - Expanding segments -", round(proc.time()[3]-timestart,3), "sec"))}
  for(i in expand){
    segtemp=segim
    segtemp[segim==i]=1
    segtemp[segim!=i]=0
    segim_new[segim_new==i]=0
    if(i %in% expand){
      temp=profoundImBlur(segtemp, expandsigma)
    }else{
      temp=segtemp
    }
    tempmult=temp*image
    segsel=tempmult>maxmat & temp>0.01 & image>skycut
    segim_new[segsel]=i
    maxmat[segsel]=tempmult[segsel]
  }
  
  objects=segim_new
  objects[objects!=0]=1
  
  if(plot){
    if(verbose){message(paste(" - Plotting segments -", round(proc.time()[3]-timestart,3), "sec"))}
    profoundSegimPlot(image=image_orig, segim=segim_new, mask=mask, sky=sky, ...)
  }else{
    if(verbose){message(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profoundSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making final local estimate of the sky - User provided sky")}
  }
  
  image_sky=image_orig-sky
  
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profoundSkyEst(image=profoundImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making final local estimate of the sky RMS - User provided sky")}
  }
  
  if(stats){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profoundSegimStats(image=image_orig, segim=segim_new, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  
  if(missing(header)){header=NULL}
  
  if(verbose){message(paste(" - profoundMakeSegimExpand is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return=list(segim=segim_new, objects=objects, sky=sky, skyRMS=skyRMS, segstats=segstats, header=header, SBlim=SBlim, call=call)
}

profoundMakeSegimDilate=function(image, segim, mask, size=9, shape='disc', expand='all', magzero=0, gain=NULL, pixscale=1, sky=0, skyRMS=0, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){

  if(verbose){message(' - Running MakeSegimDilate:')}
  timestart = proc.time()[3]
  
  if(length(image)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  
  call=match.call()
  
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  #Treat image NAs as masked regions:
  
  if(!missing(mask)){
    mask[is.na(image)]=1
  }else{
    if(any(is.na(image))){
      mask=matrix(0,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1
    }
  }
  
  if(missing(pixscale) & !missing(header)){
    pixscale=getpixscale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  kern = EBImage::makeBrush(size, shape=shape)
  
  if(verbose){message(paste(" - Dilating segments -", round(proc.time()[3]-timestart,3), "sec"))}
  
  if(expand=='all'){
    segim_new=segim
    maxorig=max(segim_new, na.rm=TRUE)
    segim_new[segim_new>0]=maxorig+1-segim_new[segim_new>0]
    segim_new=EBImage::dilate(segim_new, kern)
    segim_new[segim_new>0]=maxorig+1-segim_new[segim_new>0]
    segim_new=EBImage::imageData(segim_new)
    segim_new[segim!=0]=segim[segim!=0]
  }else{
    segim_new=segim
    segim_new[!(segim_new %in% expand)]=0
    segim_new=EBImage::dilate(segim_new, kern)
    segim_new=EBImage::imageData(segim_new)
    segim_new[segim!=0]=segim[segim!=0]
  }
  
  if(rembig){
    rm(segim)
    invisible(gc())
  }
  
  if(!missing(mask)){
    image[mask!=0]=0
    segim_new[mask!=0]=0
  }
  
  if(stats & !missing(image)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profoundSegimStats(image=image, segim=segim_new, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  objects=segim_new
  objects[objects!=0]=1
  
  if(plot){
    profoundSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }
  
  if(missing(header)){header=NULL}
  
  if(verbose){message(paste(" - profoundMakeSegimDilate is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return=list(segim=segim_new, objects=objects, segstats=segstats, header=header, call=call)
}

profoundMakeSegimPropagate=function(image, segim, objects, mask, sky=0, lambda=1e-4, plot=FALSE, ...){
  
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  if(length(image)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  
  if(!missing(mask)){
    mask[is.na(image)]=1
  }else{
    if(any(is.na(image))){
      mask=matrix(0,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1
    }
  }
  
  if(missing(objects)){
    if(!missing(segim)){
      objects=segim
      objects[objects != 0] = 1
    }
  }
  
  image_sky=image-sky
  
  if(rembig){
    rm(image)
    rm(sky)
    invisible(gc())
  }
  
  if(missing(mask)){
    propim=EBImage::imageData(EBImage::propagate(image_sky, seeds=segim, lambda=lambda))
  }else{
    #Because EBImage is odd we need to use mask to mean pixels to be propagated, i.e. where the ProFound mask=0
    propim=EBImage::imageData(EBImage::propagate(image_sky, seeds=segim, mask=(mask==0), lambda=lambda))
  }
  
  if(rembig){
    rm(segim)
    rm(mask)
    invisible(gc())
  }
  
  propim_sky=propim
  propim_sky[objects>0]=0
  
  if(plot){
    profoundSegimPlot(image=image_sky, segim=propim, mask=mask, ...)
  }
  
  return=list(propim=propim, propim_sky=propim_sky)
}

profoundSegimStats=function(image, segim, mask, sky=0, skyRMS=0, magzero=0, gain=NULL, pixscale=1, header, sortcol='segID', decreasing=FALSE, rotstats=FALSE, boundstats=FALSE, offset=1){
  
  if(length(image)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  
  if(missing(pixscale) & !missing(header)){
    pixscale=getpixscale(header)
  }
  
  if(!missing(sky)){
    hassky=any(is.finite(sky))
  }else{
    hassky=FALSE
  }
  if(hassky){
    if(length(hassky)==1 & sky[1]==0){
      hassky=FALSE
    }else{
      image=image-sky
    }
  }
  if(!missing(skyRMS)){
    hasskyRMS=any(is.finite(skyRMS))
  }else{
    hasskyRMS=FALSE
  }
  if(hasskyRMS){
    if(length(hasskyRMS)==1 & skyRMS[1]==0){hasskyRMS=FALSE}
  }
  
  #Treat image NAs as masked regions:
  
  if(!missing(mask)){
    mask[is.na(image)]=1
  }else{
    if(any(is.na(image))){
      mask=matrix(0,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1
    }
  }
  
  #Set masked things to 0, to be safe:
  
  if(!missing(mask)){
    image[mask!=0]=NA
    #segim[mask!=0]=NA
  }
  
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  
  segsel=which(segim>0)
  
  xloc = rep(1:xlen, times = ylen)[segsel]
  yloc = rep(1:ylen, each = xlen)[segsel]
  
  if(hassky & hasskyRMS){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]), sky=as.numeric(sky[segsel]), skyRMS=as.numeric(skyRMS[segsel]))
    if(rembig){
      rm(sky)
      rm(skyRMS)
      invisible(gc())
    }
  }
  if(hassky & hasskyRMS==FALSE){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]), sky=as.numeric(sky[segsel]))
    if(rembig){
      rm(sky)
      invisible(gc())
    }
  }
  if(hassky==FALSE & hasskyRMS){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]), skyRMS=as.numeric(skyRMS[segsel]))
    if(rembig){
      rm(skyRMS)
      invisible(gc())
    }
  }
  if(hassky==FALSE & hasskyRMS==FALSE){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]))
  }
  
  setkey(tempDT, segID, flux)
  
  if(rembig){
    rm(xloc)
    rm(yloc)
    rm(image)
    invisible(gc())
  }
  
  #tempDT[is.na(tempDT)]=0
  segID=tempDT[,.BY,by=segID]$segID
  
  x=NULL; y=NULL; flux=NULL; sky=NULL; skyRMS=NULL
  
  fluxout=tempDT[,.fluxcalc(flux), by=segID]
  mag=profoundFlux2Mag(flux=fluxout$flux, magzero=magzero)
  
  if(any(flux==0)){
    fluxout$N50seg[flux==0]=fluxout$N100seg[flux==0]
    fluxout$N90seg[flux==0]=fluxout$N100seg[flux==0]
  }
  
  if(hassky){
    #With one version of data.table sd doesn't work when all numbers are identical (complains about gsd and negative length vectors). Fixed with explicit sd caclulation until this gets fixed.
    flux_err_sky=tempDT[,sd(sky, na.rm=FALSE)*1, by=segID]$V1*fluxout$N100seg
    #flux_err_sky=tempDT[,sqrt(sum((sky-mean(sky, na.rm=TRUE))^2)/(.N-1)), by=segID]$V1*fluxout$N100seg
  }else{
    flux_err_sky=0
  }
  
  if(hasskyRMS){
    flux_err_skyRMS=tempDT[,sqrt(sum(skyRMS^2, na.rm=FALSE)), by=segID]$V1
    pchi=pchisq(tempDT[,sum((flux/skyRMS)^2, na.rm=TRUE), by=segID]$V1, df=fluxout$N100seg, log.p=TRUE)
    signif=qnorm(pchi, log.p=TRUE)
    FPlim=qnorm(1-fluxout$N100seg/(xlen*ylen))
  }else{
    flux_err_skyRMS=0
    signif=NA
    FPlim=NA
  }
  
  if(!is.null(gain)){
    flux_err_shot=sqrt(fluxout$flux)/gain
  }else{
    flux_err_shot=0
  }
  
  flux_err=sqrt(flux_err_sky^2+flux_err_skyRMS^2+flux_err_shot^2)
  mag_err=(2.5/log(10))*abs(flux_err/fluxout$flux)
  
  if(hassky){
    sky_mean=tempDT[,mean(sky, na.rm=FALSE), by=segID]$V1
  }else{
    sky_mean=0
  }
  
  if(hasskyRMS){
    skyRMS_mean=tempDT[,mean(skyRMS, na.rm=FALSE), by=segID]$V1
  }else{
    skyRMS_mean=0
  }
  
  xcen=tempDT[,.meanwt(x-0.5, flux),by=segID]$V1
  ycen=tempDT[,.meanwt(y-0.5, flux),by=segID]$V1
  xsd=tempDT[,sqrt(.varwt(x-0.5,flux)),by=segID]$V1
  ysd=tempDT[,sqrt(.varwt(y-0.5,flux)),by=segID]$V1
  covxy=tempDT[,.covarwt(x-0.5,y-0.5,flux),by=segID]$V1
  
  xmax=xcen
  ymax=ycen
  xmax[fluxout$N100seg>0]=tempDT[,x[which.max(flux)]-0.5,by=segID]$V1
  ymax[fluxout$N100seg>0]=tempDT[,y[which.max(flux)]-0.5,by=segID]$V1
  
  sep=sqrt((xcen-xmax)^2+(ycen-ymax)^2)*pixscale
  
  pad=10^ceiling(log10(ylen+1))
  uniqueID=ceiling(xmax)*pad+ceiling(ymax)
  
  if(rotstats){
    asymm=tempDT[,.asymm(x-0.5,y-0.5,flux),by=segID]$V1
    flux_reflect=tempDT[,.reflect(x-0.5,y-0.5,flux),by=segID]$V1
    mag_reflect=profoundFlux2Mag(flux=flux_reflect, magzero=magzero)
  }else{
    asymm=NA
    flux_reflect=NA
    mag_reflect=NA
  }
  
  corxy=covxy/(xsd*ysd)
  rad=.cov2eigval(xsd, ysd, covxy)
  rad$hi=sqrt(abs(rad$hi)+0.08333333) #Added variance of uniform in quadrature (prevents zeros)
  rad$lo=sqrt(abs(rad$lo)+0.08333333) #Added variance of uniform in quadrature (prevents zeros)
  axrat=rad$lo/rad$hi
  eigvec=.cov2eigvec(xsd, ysd, covxy)
  ang=.eigvec2ang(eigvec)
  
  R50seg=sqrt(fluxout$N50seg/(axrat*pi))*pixscale
  R90seg=sqrt(fluxout$N90seg/(axrat*pi))*pixscale
  R100seg=sqrt(fluxout$N100seg/(axrat*pi))*pixscale
  
  con=R50seg/R90seg

  SB_N50=profoundFlux2SB(flux=fluxout$flux*0.5/fluxout$N50seg, magzero=magzero, pixscale=pixscale)
  SB_N90=profoundFlux2SB(flux=fluxout$flux*0.9/fluxout$N90seg, magzero=magzero, pixscale=pixscale)
  SB_N100=profoundFlux2SB(flux=fluxout$flux/fluxout$N100seg, magzero=magzero, pixscale=pixscale)
  
  if(!missing(header)){
    coord=magWCSxy2radec(xcen, ycen, header=header)
    RAcen=coord[,1]
    Deccen=coord[,2]
    coord=magWCSxy2radec(xmax, ymax, header=header)
    RAmax=coord[,1]
    Decmax=coord[,2]
  }else{
    RAcen=NA
    Deccen=NA
    RAmax=NA
    Decmax=NA
  }
  
  if(boundstats){
    
    segim_inner=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)]
    off_down=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)-offset]
    off_left=segim[(offset+1):(xlen-offset)-offset,(offset+1):(ylen-offset)]
    off_up=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)+offset]
    off_right=segim[(offset+1):(xlen-offset)+offset,(offset+1):(ylen-offset)]
    
    inner_segim=segim_inner>0 & off_down==segim_inner & off_left==segim_inner & off_up==segim_inner & off_right==segim_inner
    
    segim_edge=segim_inner
    segim_edge[inner_segim==1]=0
    tab_edge=tabulate(segim_edge)
    tab_edge=c(tab_edge,rep(0,max(segID)-length(tab_edge)))
    tab_edge=cbind(1:max(segID),tab_edge)
    Nedge=tab_edge[match(segID,tab_edge[,1]),2]
    
    outer_sky=segim_inner>0 & (off_down==0 | off_left==0 | off_up==0 | off_right==0)
    
    segim_sky=segim_inner
    segim_sky[outer_sky==0]=0
    tab_sky=tabulate(segim_sky)
    tab_sky=c(tab_sky,rep(0,max(segID)-length(tab_sky)))
    tab_sky=cbind(1:max(segID),tab_sky)
    Nsky=tab_sky[match(segID,tab_sky[,1]),2]
    
    if(rembig){
      rm(off_down)
      rm(off_left)
      rm(off_up)
      rm(off_right)
      rm(tab_edge)
      rm(tab_sky)
      invisible(gc())
    }
    
    BorderBottom=segim[segim[,1]>0,1]
    BorderLeft=segim[1,segim[1,]>0]
    BorderTop=segim[segim[,ylen]>0,ylen]
    BorderRight=segim[xlen,segim[xlen,]>0]
    tab_bottom=tabulate(BorderBottom)
    tab_left=tabulate(BorderLeft)
    tab_top=tabulate(BorderTop)
    tab_right=tabulate(BorderRight)
    
    tab_border=cbind(1:max(segID),0,0,0,0)
    tab_border[1:length(tab_bottom),2]=tab_border[1:length(tab_bottom),2]+tab_bottom
    tab_border[1:length(tab_left),3]=tab_border[1:length(tab_left),3]+tab_left
    tab_border[1:length(tab_top),4]=tab_border[1:length(tab_top),4]+tab_top
    tab_border[1:length(tab_right),5]=tab_border[1:length(tab_right),5]+tab_right
    bordersel=match(segID,tab_border[,1])
    Nborder=tab_border[bordersel,2]+tab_border[bordersel,3]+tab_border[bordersel,4]+tab_border[bordersel,5]
    flag_border=1*(tab_border[bordersel,2]>0)+2*(tab_border[bordersel,3]>0)+4*(tab_border[bordersel,4]>0)+8*(tab_border[bordersel,5]>0)
    
    if(!missing(mask)){
      outer_mask=segim_inner>0 & (mask[2:(xlen-1)+1,2:(ylen-1)]==1 | mask[2:(xlen-1)-1,2:(ylen-1)]==1 | mask[2:(xlen-1),2:(ylen-1)+1]==1 | mask[2:(xlen-1),2:(ylen-1)-1]==1)
      segim_mask=segim_inner
      segim_mask[outer_mask==0]=0
      tab_mask=tabulate(segim_mask)
      tab_mask=c(tab_mask,rep(0,max(segID)-length(tab_mask)))
      tab_mask=cbind(1:max(segID),tab_mask)
      Nmask=tab_mask[match(segID,tab_mask[,1]),2]
      
      if(rembig){
        rm(segim_inner)
        rm(tab_mask)
        invisible(gc())
      }
    }else{
      Nmask=0
    }
    
    Nedge=Nedge+Nborder
    Nsky=Nsky-Nmask #Raw Nsky-Nmask-Nborder, to correct for masked pixels
    Nobject=Nedge-Nsky-Nmask-Nborder # Nedge-Nsky
    edge_frac=Nsky/Nedge

    #Using Ramanujan approximation from Wikipedia:
    
    A=R100seg/pixscale
    B=R100seg*axrat/pixscale
    h=(A-B)^2/(A+B)^2
    C=pi*(A+B)*(1+(3*h)/(10+sqrt(4-3*h)))
    edge_excess=Nedge/C
    
  }else{
    Nedge=NA
    Nsky=NA
    Nobject=NA
    Nborder=NA
    Nmask=NA
    flag_border=NA
    edge_frac=NA
    edge_excess=NA
  }
    
  segstats=data.table(segID=segID, uniqueID=uniqueID, xcen=xcen, ycen=ycen, xmax=xmax, ymax=ymax, RAcen=RAcen, Deccen=Deccen, RAmax=RAmax, Decmax=Decmax, sep=sep, flux=fluxout$flux, mag=mag, cenfrac=fluxout$cenfrac, N50=fluxout$N50seg, N90=fluxout$N90seg, N100=fluxout$N100seg, R50=R50seg, R90=R90seg, R100=R100seg, SB_N50=SB_N50, SB_N90=SB_N90, SB_N100=SB_N100, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=con, asymm=asymm, flux_reflect=flux_reflect, mag_reflect=mag_reflect, semimaj=rad$hi, semimin=rad$lo, axrat=axrat, ang=ang, signif=signif, FPlim=FPlim, flux_err=flux_err, mag_err=mag_err, flux_err_sky=flux_err_sky, flux_err_skyRMS=flux_err_skyRMS, flux_err_shot=flux_err_shot, sky_mean=sky_mean, sky_sum=sky_mean*fluxout$N100seg, skyRMS_mean=skyRMS_mean, Nedge=Nedge, Nsky=Nsky, Nobject=Nobject, Nborder=Nborder, Nmask=Nmask, edge_frac=edge_frac, edge_excess=edge_excess, flag_border=flag_border)
  return=as.data.frame(segstats[order(segstats[[sortcol]], decreasing=decreasing),])
}

profoundSegimPlot=function(image, segim, mask, sky=0, header, col=rainbow(max(segim), end=2/3), profound, ...){
  if(!missing(image)){
    if(class(image)=='profound'){
      if(missing(segim)){segim=image$segim}
      if(missing(mask)){mask=image$mask}
      if(missing(sky)){sky=image$sky}
      if(missing(header)){header=image$header}
      image=image$image
      if(is.null(image)){stop('Need image in profound object to be non-Null')}
    }
  }
  if(!missing(profound)){
    if(class(profound) != 'profound'){
      stop('Class of profound input must be of type \'profound\'')
    }
    if(missing(image)){image=profound$image}
    if(is.null(image)){stop('Need image in profound object to be non-Null')}
    if(missing(segim)){segim=profound$segim}
    if(missing(mask)){mask=profound$mask}
    if(missing(sky)){sky=profound$sky}
    if(missing(header)){header=profound$header}
  }
  if(!missing(image)){
    if(any(names(image)=='imDat') & missing(header)){
      header=image$hdr
      image=image$imDat
    }else if(any(names(image)=='imDat') & !missing(header)){
      image=image$imDat
    }
    if(any(names(image)=='dat') & missing(header)){
      header=image$hdr[[1]]
      header=data.frame(key=header[,1],value=header[,2], stringsAsFactors = FALSE)
      image=image$dat[[1]]
    }else if(any(names(image)=='dat') & !missing(header)){
      image=image$dat[[1]]
    }
    if(any(names(image)=='image') & missing(header)){
      header=image$header
      image=image$image
    }else if(any(names(image)=='image') & !missing(header)){
      image=image$image
    }
  }
  
  image=image-sky
  
  segim[is.na(segim)]=0
  
  if(missing(header)){header=NULL}
  if(is.null(header)){
    temp=magimage(image, ...)
  }else{
    temp=magimageWCS(image, header=header, ...)
  }
  if(min(segim,na.rm=TRUE)!=0){segim=segim-min(segim,na.rm=TRUE)}
  segvec=which(tabulate(segim)>0)
  for(i in segvec){
    z=segim==i
    z=z[ceiling(temp$x), ceiling(temp$y)]
    contour(temp$x,temp$y,z,add=T,col=col[i],zlim=c(0,1),drawlabels=FALSE,nlevels=1)
  }
  if(!missing(mask)){
    if(!is.null(mask)){
      magimage(mask, locut=0, hicut=1, col=c(NA,hsv(alpha=0.3)), add=T)
    }
  }
}

profoundSegimNear=function(segim, offset=1){
  
  if(length(segim)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  
  xlen=dim(segim)[1]
  ylen=dim(segim)[2]
  
  segim_inner=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)]
  off_down=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)-offset]
  off_left=segim[(offset+1):(xlen-offset)-offset,(offset+1):(ylen-offset)]
  off_up=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)+offset]
  off_right=segim[(offset+1):(xlen-offset)+offset,(offset+1):(ylen-offset)]
  
  tabcomb=data.table(segID=as.integer(segim_inner), down=as.integer(off_down), left=as.integer(off_left), up=as.integer(off_up), right=as.integer(off_right))
  
  if(rembig){
      rm(segim_inner)
      rm(off_down)
      rm(off_left)
      rm(off_up)
      rm(off_right)
      invisible(gc())
    }
  
  tabcomb=tabcomb[segID>0,]
  setorder(tabcomb,segID)
  
  segID=NULL; down=NULL; left=NULL; up=NULL; right=NULL; nearID=NULL; Nnear=NULL
  
  #The line means for each set of segID pixels, identify unique adjacent pixels that do not have the ID of the segID of interest or 0 (sky) and sort the touching ID list and return to a listed data.frame. Ouch!
  
  tabnear=tabcomb[,list(nearID=list(sort(setdiff(unique(c(down,left,up,right)),c(0,segID))))),by=segID]
  tabnear[,Nnear:=length(unlist(nearID)),by=segID]
  return=as.data.frame(tabnear)
}

profoundSegimGroup=function(segim){
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  Ngroup=NULL; segID=NULL; Npix=NULL
  
  groupim=EBImage::bwlabel(segim)
  segimDT=data.table(segID=as.integer(segim), groupID=as.integer(groupim))
  groupID=segimDT[groupID>0,.BY,by=groupID]$groupID
  segIDmin=segimDT[groupID>0,min(segID, na.rm=FALSE),by=groupID]$V1
  remap=vector(length=max(groupID))
  remap[groupID]=segIDmin
  groupim[groupim>0]=remap[segimDT[groupID>0,groupID]]
  
  segimDT=data.table(segID=as.integer(segim), groupID=as.integer(groupim))
  groups=segimDT[groupID>0,.N,by=groupID]
  groupID=groups$groupID
  groupsegID=segimDT[groupID>0,list(segID=list(sort(unique(segID)))),by=groupID]
  groupsegID[,Ngroup:=length(unlist(segID)),by=groupID]
  groupsegID[,Npix:=groups$N]
  setkey(groupsegID,groupID)
  return=list(groupim=groupim, groupsegID=as.data.frame(groupsegID))
}

profoundSegimMerge=function(image, segim_base, segim_add, mask, sky=0){
  
  segim_add[segim_add>0]=segim_add[segim_add>0]+max(segim_base, na.rm=FALSE)
  
  image=image-sky
  
  stats_base=profoundSegimStats(image=image, segim=segim_base, mask=mask)[,c('segID','uniqueID','flux')]
  stats_add=profoundSegimStats(image=image, segim=segim_add, mask=mask)[,c('segID','uniqueID','flux')]
  
  compID=cbind(1:length(stats_base$uniqueID), match(stats_base$uniqueID, stats_add$uniqueID))
  flux_base=stats_base[compID[,1],'flux']
  flux_add=stats_add[compID[,2],'flux']
  
  segim_base[segim_base %in% stats_base[compID[flux_base<flux_add,1],'segID']]=0
  segim_add[segim_add %in% stats_add[compID[flux_base>flux_add,2],'segID']]=0
  
  segim_merge=segim_base
  segim_merge[segim_add>0]=segim_add[segim_add>0]
  
  return=segim_merge
}

profoundSegimWarp=function(segim_in, header_in, header_out){
  return=magwarp(image_in = segim_in, header_out = header_out, header_in = header_in, doscale = FALSE, interpolation = 'nearest')$image
}
