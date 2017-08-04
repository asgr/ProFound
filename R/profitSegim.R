.meanwt=function(x, wt){
  wt[wt<0]=0
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  sum(x*wt, na.rm = T)/sum(wt, na.rm = T)
}

.varwt=function(x, wt, xcen){
  wt[wt<0]=0
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  return=(sum((x-xcen)^2*wt, na.rm = T)/sum(wt, na.rm = T))
}

.covarwt=function(x, y, wt, xcen, ycen){
  wt[wt<0]=0
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  if(missing(xcen)){xcen=.meanwt(x, wt)}
  if(missing(ycen)){ycen=.meanwt(y, wt)}
  return=(sum((x-xcen)*(y-ycen)*wt, na.rm = T)/sum(wt, na.rm = T))
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

profitMakeSegim=function(image, mask, objects, tolerance=4, ext=2, sigma=1, smooth=TRUE, pixcut=5, skycut=2, SBlim, magzero=0, gain=NULL, pixscale=1, sky, skyRMS, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){
  if(length(image)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  call=match.call()
  if(verbose){message(' - Running profitMakeSegim:')}
  timestart = proc.time()[3]
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
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
    pixscale=profitGetPixScale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  image_orig=image
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky - User provided sky")}
  }
  
  image_sky=image-sky
  
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  
  image=image_sky/skyRMS
  image[!is.finite(image)]=0
 
  if(smooth){
    if(verbose){message(paste(" - Smoothing the image -", round(proc.time()[3]-timestart,3), "sec"))}
    image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  }else{
    if(verbose){message(" - Skipping smoothing - smooth set to FALSE")}
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!missing(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profitSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
  if(!missing(mask)){
    image[mask!=0]=0
  }
  if(verbose){message(paste(" - Watershed de-blending -", round(proc.time()[3]-timestart,3), "sec"))}
  if(any(image>0)){
    segim=EBImage::imageData(EBImage::watershed(image,tolerance=tolerance,ext=ext))
  }else{
    segim=image
  }
  
  segtab=tabulate(segim)
  segim[segim %in% which(segtab<pixcut)]=0

  if(plot){
    if(verbose){message(paste(" - Plotting segments -", round(proc.time()[3]-timestart,3), "sec"))}
    profitSegimPlot(image=image_orig, segim=segim, mask=mask, sky=sky, ...)
  }else{
    if(verbose){message(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  objects=segim
  objects[objects!=0]=1
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making final local estimate of the sky - User provided sky")}
  }
  
  image_sky=image_orig-sky
  
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  
  if(stats & any(image>0)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profitSegimStats(image=image_orig, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE or no segments")}
    segstats=NULL
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  
  if(missing(header)){header=NULL}
  
  if(verbose){message(paste(" - profitMakeSegim is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return=list(segim=segim, objects=objects, sky=sky, skyRMS=skyRMS, segstats=segstats, header=header, SBlim=SBlim, call=call)
}

profitMakeSegimExpand=function(image, segim, mask, objects, skycut=1, SBlim, magzero=0, gain=NULL, pixscale=1, sigma=1, smooth=TRUE, expandsigma=5, expand='all', sky, skyRMS, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){
  
  if(verbose){message(' - Running profitMakeSegimExpand:')}
  timestart = proc.time()[3]
  
  if(length(image)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  
  call=match.call()
  
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
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
    pixscale=profitGetPixScale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  image_orig=image
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky - User provided sky")}
  }
  image_sky=image-sky
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making initial local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making initial local estimate of the sky RMS - User provided sky RMS")}
  }
  image=image_sky/skyRMS
  image[!is.finite(image)]=0

  if(smooth){
    if(verbose){message(paste(" - Smoothing image -", round(proc.time()[3]-timestart,3), "sec"))}
    image=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  }else{
    if(verbose){message(" - Skipping smoothing - smooth set to FALSE")}
  }
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  if(!missing(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profitSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
  if(!missing(mask)){
    image[mask!=0]=0
  }
  #kernel=profitMakeGaussianPSF(fwhm = expandsigma, dim=dim)
  maxmat=matrix(min(image, na.rm=TRUE), xlen, ylen)
  segim_new=matrix(0,xlen,ylen)
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  if(expand[1]=='all'){expand=segvec}
  if(verbose){message(paste(" - Expanding segments -", round(proc.time()[3]-timestart,3), "sec"))}
  for(i in segvec){
    segtemp=segim
    segtemp[segim==i]=1
    segtemp[segim!=i]=0
    if(i %in% expand){
      #temp=profitConvolvePSF(segtemp, kernel)
      temp=profitImBlur(segtemp, expandsigma)
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
    profitSegimPlot(image=image_orig, segim=segim_new, mask=mask, sky=sky, ...)
  }else{
    if(verbose){message(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  if(hassky==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky -", round(proc.time()[3]-timestart,3), "sec"))}
    sky=profitSkyEst(image=image_orig, mask=mask, objects=objects, plot=FALSE)$sky
  }else{
    if(verbose){message(" - Skipping making final local estimate of the sky - User provided sky")}
  }
  
  image_sky=image_orig-sky
  
  if(hasskyRMS==FALSE){
    if(verbose){message(paste(" - Making final local estimate of the sky RMS -", round(proc.time()[3]-timestart,3), "sec"))}
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects, plot=FALSE)$skyRMS
  }else{
    if(verbose){message(" - Skipping making final local estimate of the sky RMS - User provided sky")}
  }
  
  if(stats){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profitSegimStats(image=image_orig, segim=segim_new, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  if(!missing(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(missing(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  
  if(missing(header)){header=NULL}
  
  if(verbose){message(paste(" - profitMakeSegimExpand is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return=list(segim=segim_new, objects=objects, sky=sky, skyRMS=skyRMS, segstats=segstats, header=header, SBlim=SBlim, call=call)
}

profitMakeSegimDilate=function(image, segim, mask, size=9, shape='disc', expand='all', magzero=0, gain=NULL, pixscale=1, sky=0, skyRMS=0, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){

  if(verbose){message(' - Running profitMakeSegimDilate:')}
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
    pixscale=profitGetPixScale(header)
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
    segim_new[segim_new!=expand]=0
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
    segstats=profitSegimStats(image=image, segim=segim_new, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  objects=segim_new
  objects[objects!=0]=1
  
  if(plot){
    profitSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }
  
  if(missing(header)){header=NULL}
  
  if(verbose){message(paste(" - profitMakeSegimDilate is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return=list(segim=segim_new, objects=objects, segstats=segstats, header=header, call=call)
}

profitSegimStats=function(image, segim, mask, sky=0, skyRMS=0, magzero=0, gain=NULL, pixscale=1, header, sortcol='segID', decreasing=FALSE, rotstats=FALSE, boundstats=FALSE, offset=1){
  
  if(length(image)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  
  if(missing(pixscale) & !missing(header)){
    pixscale=profitGetPixScale(header)
  }
  
  hassky=!missing(sky)
  if(hassky){
    if(length(hassky)==1 & hassky[1]==0){
      hassky=FALSE
    }else{
      image=image-sky
    }
  }
  hasskyRMS=!missing(skyRMS)
  if(hasskyRMS){
    if(length(hasskyRMS)==1 & hasskyRMS[1]==0){hasskyRMS=FALSE}
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
    image[mask!=0]=0
    segim[mask!=0]=0
  }
  
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  #locs=expand.grid(1:xlen,1:ylen)-0.5
  
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
  
  if(rembig){
    rm(xloc)
    rm(yloc)
    rm(image)
    invisible(gc())
  }
  
  tempDT[is.na(tempDT)]=0
  segID=tempDT[,.BY,by=segID]$segID
  
  x=NULL; y=NULL; flux=NULL; sky=NULL; skyRMS=NULL
  
  flux=tempDT[,sum(flux, na.rm=TRUE),by=segID]$V1
  mag=profitFlux2Mag(flux=flux, magzero=magzero)
  
  N50seg=tempDT[,length(which(cumsum(sort(flux))/sum(flux, na.rm=TRUE)>=0.5)),by=segID]$V1
  N90seg=tempDT[,length(which(cumsum(sort(flux))/sum(flux, na.rm=TRUE)>=0.1)),by=segID]$V1
  N100seg=tempDT[,.N,by=segID]$N
  
  if(any(flux==0)){
    N50seg[flux==0]=N100seg[flux==0]
    N90seg[flux==0]=N100seg[flux==0]
  }
  
  if(hassky){
    flux_err_sky=tempDT[,sd(sky, na.rm=TRUE), by=segID]$V1*N100seg
  }else{
    flux_err_sky=0
  }
  
  if(hasskyRMS){
    flux_err_skyRMS=tempDT[,sqrt(sum(skyRMS^2, na.rm=TRUE)), by=segID]$V1
  }else{
    flux_err_skyRMS=0
  }
  
  if(!is.null(gain)){
    flux_err_shot=sqrt(flux)/gain
  }else{
    flux_err_shot=0
  }
  
  flux_err=sqrt(flux_err_sky^2+flux_err_skyRMS^2+flux_err_shot^2)
  mag_err=(2.5/log(10))*abs(flux_err/flux)
  
  if(hassky){
    sky_mean=tempDT[,mean(sky, na.rm=TRUE), by=segID]$V1
  }else{
    sky_mean=0
  }
  
  if(hasskyRMS){
    skyRMS_mean=tempDT[,mean(skyRMS, na.rm=TRUE), by=segID]$V1
  }else{
    skyRMS_mean=0
  }
  
  xcen=tempDT[,.meanwt(x-0.5, flux),by=segID]$V1
  ycen=tempDT[,.meanwt(y-0.5, flux),by=segID]$V1
  xsd=tempDT[,sqrt(.varwt(x-0.5,flux)),by=segID]$V1
  ysd=tempDT[,sqrt(.varwt(y-0.5,flux)),by=segID]$V1
  covxy=tempDT[,.covarwt(x-0.5,y-0.5,flux),by=segID]$V1
  
  xmax=tempDT[,x[which.max(flux)]-0.5,by=segID]$V1
  ymax=tempDT[,y[which.max(flux)]-0.5,by=segID]$V1
  
  sep=sqrt((xcen-xmax)^2+(ycen-ymax)^2)*pixscale
  
  pad=10^ceiling(log10(ylen))
  uniqueID=ceiling(xcen)*pad+ceiling(ycen)
  
  if(rotstats){
    asymm=tempDT[,.asymm(x-0.5,y-0.5,flux),by=segID]$V1
    flux_reflect=tempDT[,.reflect(x-0.5,y-0.5,flux),by=segID]$V1
    mag_reflect=profitFlux2Mag(flux=flux_reflect, magzero=magzero)
  }else{
    asymm=NA
    flux_reflect=NA
    mag_reflect=NA
  }
  
  corxy=covxy/(xsd*ysd)
  rad=.cov2eigval(xsd, ysd, covxy)
  rad$hi=sqrt(abs(rad$hi))
  rad$lo=sqrt(abs(rad$lo))
  axrat=rad$lo/rad$hi
  eigvec=.cov2eigvec(xsd, ysd, covxy)
  ang=.eigvec2ang(eigvec)
  
  R50seg=sqrt(N50seg/(axrat*pi))*pixscale
  R90seg=sqrt(N90seg/(axrat*pi))*pixscale
  R100seg=sqrt(N100seg/(axrat*pi))*pixscale
  
  con=R50seg/R90seg

  SB_N50=profitFlux2SB(flux=flux*0.5/N50seg, magzero=magzero, pixscale=pixscale)
  SB_N90=profitFlux2SB(flux=flux*0.9/N90seg, magzero=magzero, pixscale=pixscale)
  SB_N100=profitFlux2SB(flux=flux/N100seg, magzero=magzero, pixscale=pixscale)
  
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
    Flagborder=1*(tab_border[bordersel,2]>0)+2*(tab_border[bordersel,3]>0)+4*(tab_border[bordersel,4]>0)+8*(tab_border[bordersel,5]>0)
    
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
    Flagborder=NA
    edge_frac=NA
    edge_excess=NA
  }
    
  segstats=data.table(segID=segID, uniqueID=uniqueID, xcen=xcen, ycen=ycen, xmax=xmax, ymax=ymax, RAcen=RAcen, Deccen=Deccen, RAmax=RAmax, Decmax=Decmax, sep=sep, flux=flux, mag=mag, N50=N50seg, N90=N90seg, N100=N100seg, R50=R50seg, R90=R90seg, R100=R100seg, SB_N50=SB_N50, SB_N90=SB_N90, SB_N100=SB_N100, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=con, asymm=asymm, flux_reflect=flux_reflect, mag_reflect=mag_reflect, maj=rad$hi, min=rad$lo, axrat=axrat, ang=ang, flux_err=flux_err, mag_err=mag_err, flux_err_sky=flux_err_sky, flux_err_skyRMS=flux_err_skyRMS, flux_err_shot=flux_err_shot, sky_mean=sky_mean, sky_sum=sky_mean*N100seg, skyRMS_mean=skyRMS_mean, Nedge=Nedge, Nsky=Nsky, Nobject=Nobject, Nborder=Nborder, Nmask=Nmask, edge_frac=edge_frac, edge_excess=edge_excess, Flagborder=Flagborder)
  return=as.data.frame(segstats[order(segstats[[sortcol]], decreasing=decreasing),])
}

profitSegimPlot=function(image, segim, mask, sky=0, ...){
  image=image-sky
  temp=magimage(image, ...)
  if(min(segim,na.rm=TRUE)!=0){segim=segim-min(segim,na.rm=TRUE)}
  segvec=which(tabulate(segim)>0)
  for(i in segvec){
    z=segim==i
    z=z[ceiling(temp$x), ceiling(temp$y)]
    contour(temp$x,temp$y,z,add=T,col=rainbow(1e3)[sample(1e3,1)],zlim=c(0,1),drawlabels=FALSE,nlevels=1)
  }
  if(!missing(mask)){
    magimage(mask, lo=0, hi=1, col=c(NA,hsv(alpha=0.3)), add=T)
  }
}