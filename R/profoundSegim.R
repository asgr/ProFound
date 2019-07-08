.meanwt=function(x=NULL, wt=NULL){
  wt[wt<0]=0
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  sum(x*wt, na.rm=TRUE)/sum(wt, na.rm=TRUE)
}

.varwt=function(x=NULL, wt=NULL, xcen=NULL){
  wt[wt<0]=0
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  if(is.null(xcen)){xcen=.meanwt(x, wt)}
  invisible((sum((x-xcen)^2*wt, na.rm=TRUE)/sum(wt, na.rm=TRUE)))
}

.covarwt=function(x=NULL, y=NULL, wt=NULL, xcen=NULL, ycen=NULL){
  wt[wt<0]=0
  if(all(wt==wt[1], na.rm=TRUE)){wt[]=1}
  if(is.null(xcen)){xcen=.meanwt(x, wt)}
  if(is.null(ycen)){ycen=.meanwt(y, wt)}
  invisible((sum((x-xcen)*(y-ycen)*wt, na.rm=TRUE)/sum(wt, na.rm=TRUE)))
}

.cov2eigval=function(sx=NULL, sy=NULL, sxy=NULL){
  b=-sx^2-sy^2
  c=sx^2*sy^2-sxy^2
  invisible((list(hi=(-b+sqrt(b^2-4*c))/2,lo=(-b-sqrt(b^2-4*c))/2)))
}

.cov2eigvec=function(sx=NULL, sy=NULL, sxy=NULL){
  eigval=.cov2eigval(sx,sy,sxy)$hi
  eigvec=(sx^2-eigval)/sxy
  invisible(eigvec)
}

.eigvec2ang=function(eigvec=NULL){
  invisible((90-atan(eigvec)*180/pi)%%180)
}

.asymm=function(x=NULL, y=NULL, wt=NULL, xcen=NULL, ycen=NULL){
  if(is.null(xcen)){xcen=.meanwt(x, wt)}
  if(is.null(ycen)){ycen=.meanwt(y, wt)}
  relx=round(x-xcen)
  rely=round(y-ycen)
  frame1=data.frame(x=relx,y=rely,wt1=wt)
  frame2=data.frame(x=-relx,y=-rely,wt2=wt)
  comp=merge(frame1,frame2,by=c('x','y'), all=TRUE)
  overlap=which(comp$wt1>0 & comp$wt2>0)
  asymm=sum(abs(comp[overlap,'wt1']-comp[overlap,'wt2']), na.rm=TRUE)/sum(abs(comp[overlap,'wt1']+comp[overlap,'wt2']), na.rm=TRUE)
  invisible(asymm)
}

.reflect=function(x=NULL, y=NULL, wt=NULL, xcen=NULL, ycen=NULL){
  if(is.null(xcen)){xcen=.meanwt(x, wt)}
  if(is.null(ycen)){ycen=.meanwt(y, wt)}
  relx=round(x-xcen)
  rely=round(y-ycen)
  frame1=data.frame(x=relx,y=rely,wt1=wt)
  frame2=data.frame(x=-relx,y=-rely,wt2=wt)
  comp=merge(frame1,frame2,by=c('x','y'),all=TRUE)
  overlap=is.na(comp$wt1)==FALSE & is.na(comp$wt2)==FALSE
  asymm=2*sum(abs(comp[overlap,'wt1']-comp[overlap,'wt2']), na.rm=TRUE)/sum(comp[overlap,'wt1'], comp[overlap,'wt2'], na.rm=TRUE)
  len=dim(frame1)[1]
  flux_reflect=sum(comp[overlap,'wt1'], na.rm=TRUE)+2*sum(frame1[is.na(comp$wt2),'wt1'], na.rm=TRUE)
  invisible(flux_reflect)
}

#Not currently used:
.nser2ccon=function(nser=0.5, lo=0.5, hi=0.9){
  invisible((((qgamma(lo, 2 * nser)/qgamma(hi, 2 * nser))^nser)^2))
}

#Not currently used (too slow):
.match2col=function(tab1, tab2){
  invisible(which(outer(tab1[,1], tab2[,1], "==") & outer(tab1[,2], tab2[,2], "=="), arr.ind=TRUE))
}

.fluxcalc=function(flux, Napp=0){
  
  if(anyNA(flux)){
    good=which(!is.na(flux))
    N100seg=length(good)
  }else{
    good=TRUE
    N100seg=length(flux)
  }
  
  if(Napp>0){
    Nsel=N100seg:(N100seg-Napp+1)
    Nsel=Nsel[Nsel>0]
  }else{
    Nsel=0
  }
  
  if(N100seg>0){
    
    sumflux=sum(flux[good])
    
    if(length(Nsel)>0){
      if(Nsel[1]>0){
        sumflux_app=sum(flux[good][Nsel])
      }else{
        sumflux_app=0
      }
    }else{
      sumflux_app=0
    }
    
    temp=cumsum(flux[good])/sumflux
    
    if(sumflux>0){
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
      N100seg=length(flux)
      N50seg=NA
      N90seg=NA
      cenfrac=NA
    }
    
  }else{
    sumflux=NA
    sumflux_app=NA
    N100seg=length(flux)
    N50seg=N100seg*0.5
    N90seg=N100seg*0.9
    cenfrac=NA
  }
  
  mode(sumflux)='numeric'
  mode(sumflux_app)='numeric'
  mode(N50seg)='numeric'
  mode(N90seg)='numeric'
  mode(cenfrac)='numeric'
  mode(cenfrac)='numeric'
  
  invisible(list(flux=sumflux, flux_app=sumflux_app, N50seg=N50seg, N90seg=N90seg, N100seg=N100seg, cenfrac=cenfrac))
}

.fluxcalcmin=function(flux){
  
  if(anyNA(flux)){
    good=which(!is.na(flux))
    N100seg=length(good)
  }else{
    good=TRUE
    N100seg=length(flux)
  }
  
  if(N100seg>0){
    return(list(flux=as.numeric(sum(flux[good], na.rm=TRUE)),N100=as.integer(N100seg)))
  }else{
    return(list(flux=as.numeric(0),N100=as.integer(0)))
  }
}

profoundMakeSegim=function(image=NULL, mask=NULL, objects=NULL, skycut=1, pixcut=3, tolerance=4, ext=2, reltol=0, cliptol=Inf, sigma=1, smooth=TRUE, SBlim=NULL, magzero=0, gain=NULL, pixscale=1, sky=NULL, skyRMS=NULL, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, watershed = 'ProFound', ...){
  
  call=match.call()
  if(verbose){message(' - Running MakeSegim:')}
  timestart = proc.time()[3]
  # if(!requireNamespace("imager", quietly = TRUE)){
  #   stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  # }
  
  #Treat image NAs as masked regions:
  
  if(!is.null(mask)){
    if(anyNA(image)){
      mask[is.na(image)]=1L
    }
  }else{
    if(anyNA(image)){
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1L
    }
  }
  
  if(missing(pixscale) & !is.null(header)){
    pixscale=getpixscale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  hassky=!is.null(sky)
  hasskyRMS=!is.null(skyRMS)
  
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
  if(!is.null(SBlim) & !missing(magzero)){
    #image[image<skycut | image_sky<profoundSB2Flux(SBlim, magzero, pixscale)]=0
    image[image_sky<profoundSB2Flux(SBlim, magzero, pixscale)]=0
  }
  if(!is.null(mask)){
    image[mask>0]=0
  }
  if(verbose){message(paste(" - Watershed de-blending -", round(proc.time()[3]-timestart,3), "sec"))}
  if(any(image>0)){
    if(watershed=='EBImage'){
      if(!requireNamespace("EBImage", quietly = TRUE)){
        stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
      }
      image[image<skycut]=0
      segim=EBImage::imageData(EBImage::watershed(image,tolerance=tolerance,ext=ext))
      segtab=tabulate(segim)
      segim[segim %in% which(segtab<pixcut)]=0L
      mode(segim)='integer'
    }else if(watershed=='ProFound'){
      segim=water_cpp(image=image, nx=dim(image)[1], ny=dim(image)[2], abstol=tolerance, reltol=reltol, cliptol=cliptol, ext=ext, skycut=skycut, pixcut=pixcut, verbose=verbose)
    }else if(watershed=='ProFound-old'){
      segim=water_cpp_old(image=image, nx=dim(image)[1], ny=dim(image)[2], abstol=tolerance, reltol=reltol, cliptol=cliptol, ext=ext, skycut=skycut, pixcut=pixcut, verbose=verbose)
    }else{
      stop('watershed option must either be EBImage/ProFound/ProFound-old!')
    }
  }else{
    segim=image
  }
  
  if(plot){
    if(verbose){message(paste(" - Plotting segments -", round(proc.time()[3]-timestart,3), "sec"))}
    profoundSegimPlot(image=image_orig, segim=segim, mask=mask, sky=sky, ...)
  }else{
    if(verbose){message(" - Skipping segmentation plot - plot set to FALSE")}
  }
  
  objects=matrix(0L,dim(segim)[1],dim(segim)[2])
  objects[]=as.logical(segim)
  
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
    if(verbose){message(" - Skipping making final local estimate of the sky RMS - User provided sky RMS")}
  }
  
  if(stats & any(image>0)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profoundSegimStats(image=image_orig, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE or no segments")}
    segstats=NULL
  }
  
  if(!is.null(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(is.null(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  
  if(is.null(header)){header=NULL}
  
  if(verbose){message(paste(" - MakeSegim is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  invisible(list(segim=segim, objects=objects, sky=sky, skyRMS=skyRMS, segstats=segstats, header=header, SBlim=SBlim, call=call))
}

profoundMakeSegimExpand=function(image=NULL, segim=NULL, mask=NULL, objects=NULL, skycut=1, SBlim=NULL, magzero=0, gain=NULL, pixscale=1, sigma=1, smooth=TRUE, expandsigma=5, expand='all', sky=NULL, skyRMS=NULL, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){
  
  if(verbose){message(' - Running MakeSegimExpand:')}
  timestart = proc.time()[3]
  
  call=match.call()
  
  # if(!requireNamespace("imager", quietly = TRUE)){
  #   stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  # }
  
  #Treat image NAs as masked regions:
  
  if(!is.null(mask)){
    if(anyNA(image)){
      mask[is.na(image)]=1L
    }
  }else{
    if(anyNA(image)){
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1L
    }
  }
  
  if(missing(pixscale) & !is.null(header)){
    pixscale=getpixscale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  hassky=!is.null(sky)
  hasskyRMS=!is.null(skyRMS)
  
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
  if(!is.null(SBlim) & !missing(magzero)){
    image[image<skycut | image_sky<profoundSB2Flux(SBlim, magzero, pixscale)]=0
  }else{
    image[image<skycut]=0
  }
  if(!is.null(mask)){
    image[mask!=0]=NA
  }
  #kernel=profoundMakeGaussianPSF(fwhm = expandsigma, dim=dim)
  maxmat=matrix(min(image, na.rm=TRUE), xlen, ylen)
  segim_new=segim
  segvec=which(tabulate(segim)>0)
  segvec=segvec[segvec>0]
  
  if(is.null(expand) | length(expand)==0){
    objects=matrix(0L,dim(segim)[1],dim(segim)[2])
    objects[]=as.logical(segim)
    
    if(stats){
      if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
      segstats=profoundSegimStats(image=image_orig, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
    }else{
      if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
      segstats=NULL
    }
    
    return(invisible(list(segim=segim, objects=objects, segstats=segstats, header=header, call=call)))
  }
  
  if(expand[1]=='all'){expand=segvec}
  if(verbose){message(paste(" - Expanding segments -", round(proc.time()[3]-timestart,3), "sec"))}
  for(i in expand){
    segtemp=segim
    segtemp[segim==i]=1L
    segtemp[segim!=i]=0L
    segim_new[segim_new==i]=0L
    if(i %in% expand){
      temp=profoundImBlur(segtemp, expandsigma)
    }else{
      temp=segtemp
    }
    tempmult=temp*image
    segsel=which(tempmult>maxmat & temp>0.01 & image>skycut)
    segim_new[segsel]=i
    maxmat[segsel]=tempmult[segsel]
  }
  mode(segim_new)='integer'
  segim_new[segim>0]=segim[segim>0]
  
  if(!is.null(mask)){
    segim_new[mask!=0]=0
  }
  
  objects=matrix(0L,dim(segim_new)[1],dim(segim_new)[2])
  objects[]=as.logical(segim_new)
  
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
  
  if(!is.null(SBlim) & !missing(magzero)){
    SBlim=min(SBlim, profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
  }else if(is.null(SBlim) & !missing(magzero) & skycut>0){
    SBlim=profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
  }else{
    SBlim=NULL
  }
  
  if(is.null(header)){header=NULL}
  
  if(verbose){message(paste(" - profoundMakeSegimExpand is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return(invisible(list(segim=segim_new, objects=objects, sky=sky, skyRMS=skyRMS, segstats=segstats, header=header, SBlim=SBlim, call=call)))
}

profoundMakeSegimDilate=function(image=NULL, segim=NULL, mask=NULL, size=9, shape='disc', expand='all', magzero=0, gain=NULL, pixscale=1, sky=0, skyRMS=0, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, offset=1, sortcol = "segID", decreasing = FALSE, ...){
  
  if(verbose){message(' - Running MakeSegimDilate:')}
  timestart = proc.time()[3]
  
  call=match.call()
  
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  #Treat image NAs as masked regions:
  
  if(!is.null(mask) & !is.null(image)){
    mask[is.na(image)]=1L
  }else{
    if(anyNA(image)){
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1L
    }
  }
  
  if(missing(pixscale) & !is.null(header)){
    pixscale=getpixscale(header)
    if(verbose){message(paste(' - Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  kern = EBImage::makeBrush(size, shape=shape)
  
  if(verbose){message(paste(" - Dilating segments -", round(proc.time()[3]-timestart,3), "sec"))}
  
  if(is.null(expand) | length(expand)==0){
    objects=matrix(0L,dim(segim)[1],dim(segim)[2])
    objects[]=as.logical(segim)
    
    if(stats){
      if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
      segstats=profoundSegimStats(image=image, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
    }else{
      if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
      segstats=NULL
    }
    
    return(invisible(list(segim=segim, objects=objects, segstats=segstats, header=header, call=call)))
  }
  
  if(expand[1]=='all'){
    segim_new=segim
    maxorig=max(segim_new, na.rm=TRUE)+1L
    replace=which(segim_new!=0)
    segim_new[replace]=maxorig-segim_new[replace]
    segim_new=EBImage::imageData(EBImage::dilate(segim_new, kern)) #Run Dilate
    replace=which(segim_new!=0)
    segim_new[replace]=maxorig-segim_new[replace]
    replace=which(segim!=0) #put back non-dilated segments
    segim_new[replace]=segim[replace] #put back non-dilated segments
  }else{
    segim_new=segim
    #segim_new[!(segim_new %in% expand)]=0L #remove things that will not be dilated
    if('fastmatch' %in% .packages()){ #remove things that will not be dilated
      segim_new[fastmatch::fmatch(segim_new, expand, nomatch = 0L) == 0L] = 0L
    }else{
      segim_new[!(segim_new %in% expand)] = 0L
    }
    maxorig=max(segim_new, na.rm=TRUE)+1L
    replace=which(segim_new!=0)
    segim_new[replace]=maxorig-segim_new[replace]
    segim_new=EBImage::imageData(EBImage::dilate(segim_new, kern)) #Run Dilate
    replace=which(segim_new!=0)
    segim_new[replace]=maxorig-segim_new[replace]
    replace=which(segim!=0) #put back non-dilated segments
    segim_new[replace]=segim[replace] #put back non-dilated segments
    rm(replace)
  }
  mode(segim_new)='integer'
  
  rm(segim)
  
  if(!is.null(mask) & !is.null(image)){
    segim_new[mask!=0]=0
    image[mask!=0]=NA
  }
  
  if(stats & !is.null(image)){
    if(verbose){message(paste(" - Calculating segstats -", round(proc.time()[3]-timestart,3), "sec"))}
    segstats=profoundSegimStats(image=image, segim=segim_new, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
  }else{
    if(verbose){message(" - Skipping segmentation statistics - segstats set to FALSE")}
    segstats=NULL
  }
  
  objects=matrix(0L,dim(segim_new)[1],dim(segim_new)[2])
  objects[]=as.logical(segim_new)
  
  if(plot & !is.null(image)){
    profoundSegimPlot(image=image, segim=segim_new, mask=mask, sky=sky, ...)
  }
  
  if(is.null(header)){header=NULL}
  
  if(verbose){message(paste(" - profoundMakeSegimDilate is finished! -", round(proc.time()[3]-timestart,3), "sec"))}
  
  return(invisible(list(segim=segim_new, objects=objects, segstats=segstats, header=header, call=call)))
}

profoundMakeSegimPropagate=function(image=NULL, segim=NULL, objects=NULL, mask=NULL, sky=0, lambda=1e-4, plot=FALSE, ...){
  
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  if(!is.null(mask)){
    mask[is.na(image)]=1L
  }else{
    if(anyNA(image)){
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1L
    }
  }
  
  if(is.null(objects)){
    if(!is.null(segim)){
      objects=matrix(0L,dim(segim)[1],dim(segim)[2])
      objects[]=as.logical(segim)
    }
  }
  
  image_sky=image-sky
  
  rm(image)
  rm(sky)
  
  if(is.null(mask)){
    propim=EBImage::imageData(EBImage::propagate(image_sky, seeds=segim, lambda=lambda))
  }else{
    #Because EBImage is odd we need to use mask to mean pixels to be propagated, i.e. where the ProFound mask=0
    propim=EBImage::imageData(EBImage::propagate(image_sky, seeds=segim, mask=(mask==0), lambda=lambda))
  }
  mode(propim)='integer'
  
  rm(segim)
  rm(mask)
  
  propim_sky=propim
  propim_sky[objects>0]=0L
  mode(propim_sky)='integer'
  
  if(plot){
    profoundSegimPlot(image=image_sky, segim=propim, mask=mask, ...)
  }
  
  invisible(list(propim=propim, propim_sky=propim_sky))
}

profoundSegimStats=function(image=NULL, segim=NULL, mask=NULL, sky=NULL, skyRMS=NULL, magzero=0, gain=NULL, pixscale=1, header=NULL, sortcol='segID', decreasing=FALSE, rotstats=FALSE, boundstats=FALSE, offset=1, cor_err_func=NULL, app_diam=1){
  
  if(missing(pixscale) & !is.null(header)){
    pixscale=getpixscale(header)
  }
  
  Napp=ceiling(pi*(app_diam/2/pixscale)^2)
  
  if(!is.null(sky)){
    hassky=any(is.finite(sky))
    if(hassky & length(sky)==1){
      sky=rep(sky,length(image))
    }
  }else{
    hassky=FALSE
  }
  if(hassky){
    image=image-sky
  }
  if(!is.null(skyRMS)){
    hasskyRMS=any(is.finite(skyRMS))
    if(hasskyRMS & length(skyRMS)==1){
      skyRMS=rep(skyRMS,length(image))
    }
  }else{
    hasskyRMS=FALSE
  }
  
  #Treat image NAs as masked regions:
  
  if(!is.null(mask)){
    mask[is.na(image)]=1L
  }else{
    if(anyNA(image)){
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      mask[is.na(image)]=1L
    }
  }
  
  #Set masked things to NA, to be safe:
  
  if(!is.null(mask)){
    image[mask!=0]=NA
    #segim[mask!=0]=NA
  }
  
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  #segvec=which(tabulate(segim)>0)
  #segvec=segvec[segvec>0]
  
  segsel=which(segim>0)
  
  xloc = rep(1:xlen, times = ylen)[segsel]
  yloc = rep(1:ylen, each = xlen)[segsel]
  
  if(hassky & hasskyRMS){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]), sky=as.numeric(sky[segsel]), skyRMS=as.numeric(skyRMS[segsel]))
    rm(sky)
    rm(skyRMS)
  }
  if(hassky & hasskyRMS==FALSE){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]), sky=as.numeric(sky[segsel]))
    rm(sky)
  }
  if(hassky==FALSE & hasskyRMS){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]), skyRMS=as.numeric(skyRMS[segsel]))
    rm(skyRMS)
  }
  if(hassky==FALSE & hasskyRMS==FALSE){
    tempDT=data.table(segID=as.integer(segim[segsel]), x=xloc, y=yloc, flux=as.numeric(image[segsel]))
  }
  
  setkey(tempDT, segID, flux)
  
  rm(xloc)
  rm(yloc)
  rm(image)
  
  #tempDT[is.na(tempDT)]=0
  segID=tempDT[,.BY,by=segID]$segID
  
  x=NULL; y=NULL; flux=NULL; sky=NULL; skyRMS=NULL
  
  fluxout=tempDT[,.fluxcalc(flux,Napp=Napp), by=segID]
  fluxout$flux_app[which(fluxout$flux_app>fluxout$flux)]=fluxout$flux[which(fluxout$flux_app>fluxout$flux)]
  mag=profoundFlux2Mag(flux=fluxout$flux, magzero=magzero)
  mag_app=profoundFlux2Mag(flux=fluxout$flux_app, magzero=magzero)
  
  if(any(fluxout$flux==0, na.rm=TRUE)){
    fluxout$N50seg[fluxout$flux==0]=fluxout$N100seg[fluxout$flux==0]
    fluxout$N90seg[fluxout$flux==0]=fluxout$N100seg[fluxout$flux==0]
  }
  
  if(hassky){
    #With one version of data.table sd doesn't work when all numbers are identical (complains about gsd and negative length vectors). Fixed with explicit sd caclulation until this gets fixed.
    flux_err_sky=tempDT[,sd(sky, na.rm=TRUE)*1, by=segID]$V1*fluxout$N100seg
    #flux_err_sky=tempDT[,sqrt(sum((sky-mean(sky, na.rm=TRUE))^2)/(.N-1)), by=segID]$V1*fluxout$N100seg
  }else{
    flux_err_sky=0
  }
  
  if(hasskyRMS){
    flux_err_skyRMS=tempDT[,sqrt(sum(skyRMS^2, na.rm=TRUE)), by=segID]$V1
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
  
  if(!is.null(cor_err_func)){
    cor_seg=cor_err_func(fluxout$N100seg)
    flux_err_cor=sqrt((flux_err_sky^2)/(1-cor_seg)-flux_err_sky^2)
  }else{
    cor_seg=0
    flux_err_cor=0
  }
  
  flux_err_sky[!is.finite(flux_err_sky)]=0
  flux_err_skyRMS[!is.finite(flux_err_skyRMS)]=0
  flux_err_shot[!is.finite(flux_err_shot)]=0
  flux_err_cor[!is.finite(flux_err_cor)]=0
  
  flux_err=sqrt(flux_err_sky^2+flux_err_skyRMS^2+flux_err_shot^2+flux_err_cor^2)
  mag_err=(2.5/log(10))*abs(flux_err/fluxout$flux)
  
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
  
  xmax=xcen
  ymax=ycen
  xmax[!is.na(fluxout$flux)]=tempDT[,x[which.max(flux)]-0.5,by=segID]$V1
  ymax[!is.na(fluxout$flux)]=tempDT[,y[which.max(flux)]-0.5,by=segID]$V1
  
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
  con[R90seg==0]=NA
  
  SB_N50=profoundFlux2SB(flux=fluxout$flux*0.5/fluxout$N50seg, magzero=magzero, pixscale=pixscale)
  SB_N90=profoundFlux2SB(flux=fluxout$flux*0.9/fluxout$N90seg, magzero=magzero, pixscale=pixscale)
  SB_N100=profoundFlux2SB(flux=fluxout$flux/fluxout$N100seg, magzero=magzero, pixscale=pixscale)
  
  if(!is.null(header)){
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
    
    rm(off_down)
    rm(off_left)
    rm(off_up)
    rm(off_right)
    rm(tab_edge)
    rm(tab_sky)
    
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
    
    if(!is.null(mask)){
      outer_mask=segim_inner>0 & (mask[2:(xlen-1)+1,2:(ylen-1)]==1 | mask[2:(xlen-1)-1,2:(ylen-1)]==1 | mask[2:(xlen-1),2:(ylen-1)+1]==1 | mask[2:(xlen-1),2:(ylen-1)-1]==1)
      segim_mask=segim_edge
      segim_mask[outer_mask==0]=0
      tab_mask=tabulate(segim_mask)
      tab_mask=c(tab_mask,rep(0,max(segID)-length(tab_mask)))
      tab_mask=cbind(1:max(segID),tab_mask)
      Nmask=tab_mask[match(segID,tab_mask[,1]),2]
      
      rm(segim_inner)
      rm(tab_mask)
      
    }else{
      Nmask=0
    }
    
    Nedge=Nedge+Nborder
    #Nsky=Nsky-Nmask #Raw Nsky-Nmask-Nborder, to correct for masked pixels
    Nobject=Nedge-Nsky-Nborder # Nedge-Nsky
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
  
  if(anyNA(fluxout$flux)){
    bad=is.na(fluxout$flux)
    asymm[bad]=NA
    flux_reflect[bad]=NA
    mag_reflect[bad]=NA
    signif[bad]=NA
    FPlim[bad]=NA
    edge_excess[bad]=NA
  }
  
  segstats=data.table(segID=segID, uniqueID=uniqueID, xcen=xcen, ycen=ycen, xmax=xmax, ymax=ymax, RAcen=RAcen, Deccen=Deccen, RAmax=RAmax, Decmax=Decmax, sep=sep, flux=fluxout$flux, mag=mag, flux_app=fluxout$flux_app, mag_app=mag_app, cenfrac=fluxout$cenfrac, N50=fluxout$N50seg, N90=fluxout$N90seg, N100=fluxout$N100seg, R50=R50seg, R90=R90seg, R100=R100seg, SB_N50=SB_N50, SB_N90=SB_N90, SB_N100=SB_N100, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy, con=con, asymm=asymm, flux_reflect=flux_reflect, mag_reflect=mag_reflect, semimaj=rad$hi, semimin=rad$lo, axrat=axrat, ang=ang, signif=signif, FPlim=FPlim, flux_err=flux_err, mag_err=mag_err, flux_err_sky=flux_err_sky, flux_err_skyRMS=flux_err_skyRMS, flux_err_shot=flux_err_shot, flux_err_cor=flux_err_cor, cor_seg=cor_seg, sky_mean=sky_mean, sky_sum=sky_mean*fluxout$N100seg, skyRMS_mean=skyRMS_mean, Nedge=Nedge, Nsky=Nsky, Nobject=Nobject, Nborder=Nborder, Nmask=Nmask, edge_frac=edge_frac, edge_excess=edge_excess, flag_border=flag_border)
  invisible(as.data.frame(segstats[order(segstats[[sortcol]], decreasing=decreasing),]))
}

profoundSegimPlot=function(image=NULL, segim=NULL, mask=NULL, sky=NULL, header=NULL, col=rainbow(max(segim), end=2/3), profound=NULL, ...){
  if(!is.null(image)){
    if(class(image)=='profound'){
      if(is.null(segim)){segim=image$segim}
      if(is.null(mask)){mask=image$mask}
      if(is.null(sky)){sky=image$sky}
      if(is.null(header)){header=image$header}
      image=image$image
      if(is.null(image)){stop('Need image in profound object to be non-Null')}
    }
  }
  if(!is.null(profound)){
    if(class(profound) != 'profound'){
      stop('Class of profound input must be of type \'profound\'')
    }
    if(is.null(image)){image=profound$image}
    if(is.null(image)){stop('Need image in profound object to be non-Null')}
    if(is.null(segim)){segim=profound$segim}
    if(is.null(mask)){mask=profound$mask}
    if(is.null(sky)){sky=profound$sky}
    if(is.null(header)){header=profound$header}
  }
  if(!is.null(image)){
    if(any(names(image)=='imDat') & is.null(header)){
      header=image$hdr
      image=image$imDat
    }else if(any(names(image)=='imDat') & !is.null(header)){
      image=image$imDat
    }
    if(any(names(image)=='dat') & is.null(header)){
      header=image$hdr[[1]]
      header=data.frame(key=header[,1],value=header[,2], stringsAsFactors = FALSE)
      image=image$dat[[1]]
    }else if(any(names(image)=='dat') & !is.null(header)){
      image=image$dat[[1]]
    }
    if(any(names(image)=='image') & is.null(header)){
      header=image$header
      image=image$image
    }else if(any(names(image)=='image') & !is.null(header)){
      image=image$image
    }
  }
  
  if(!is.null(sky)){
    image=image-sky
  }
  
  segim[is.na(segim)]=0L
  
  if(is.null(header)){header=NULL}
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
  if(!is.null(mask)){
    if(!is.null(mask)){
      magimage(mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
  }
}

profoundSegimNear=function(segim=NULL, offset=1){
  
  xlen=dim(segim)[1]
  ylen=dim(segim)[2]
  
  segim_inner=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)]
  off_down=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)-offset]
  off_left=segim[(offset+1):(xlen-offset)-offset,(offset+1):(ylen-offset)]
  off_up=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)+offset]
  off_right=segim[(offset+1):(xlen-offset)+offset,(offset+1):(ylen-offset)]
  
  tabcomb=data.table(segID=as.integer(segim_inner), down=as.integer(off_down), left=as.integer(off_left), up=as.integer(off_up), right=as.integer(off_right))
  
  rm(segim_inner)
  rm(off_down)
  rm(off_left)
  rm(off_up)
  rm(off_right)
  
  tabcomb=tabcomb[segID>0,]
  setorder(tabcomb,segID)
  
  segID=NULL; down=NULL; left=NULL; up=NULL; right=NULL; nearID=NULL; Nnear=NULL
  
  #The line means for each set of segID pixels, identify unique adjacent pixels that do not have the ID of the segID of interest or 0 (sky) and sort the touching ID list and return to a listed data.frame. Ouch!
  
  tabnear=tabcomb[,list(nearID=list(sort(setdiff(unique(c(down,left,up,right)),c(0,segID))))),by=segID]
  tabnear[,Nnear:=length(unlist(nearID)),by=segID]
  invisible(as.data.frame(tabnear))
}

profoundSegimGroup=function(segim=NULL){
  if(!requireNamespace("EBImage", quietly = TRUE)){
    stop('The EBImage package is needed for this function to work. Please install it from Bioconductor.', call. = FALSE)
  }
  
  Ngroup=NULL; segID=NULL; Npix=NULL
  
  groupim=EBImage::bwlabel(segim)
  segimDT=data.table(segID=as.integer(segim), groupID=as.integer(groupim))
  segimDT[groupID>0,groupID:=which.max(tabulate(groupID)),by=segID]
  groupID=segimDT[groupID>0,.BY,by=groupID]$groupID
  segIDmin=segimDT[groupID>0,min(segID, na.rm=FALSE),by=groupID]$V1
  remap=vector(length=max(groupID))
  remap[groupID]=segIDmin
  groupim[groupim>0]=remap[segimDT[groupID>0,groupID]]
  
  segimDT=data.table(segID=as.integer(segim), groupID=as.integer(groupim))
  segimDT=segimDT[groupID>0,]
  groups=segimDT[,.N,by=groupID]
  groupsegID=segimDT[,list(segID=list(sort(unique(segID)))),by=groupID]
  groupsegID[,Npix:=groups$N]
  setkey(groupsegID,groupID)
  groupsegID=groupsegID[groupID %in% which(tabulate(unlist(groupsegID$segID))==1),]
  groupsegID[,Ngroup:=length(unlist(segID)),by=groupID]
  groupim[!groupim %in% groupsegID$groupID]=0
  invisible(list(groupim=groupim, groupsegID=as.data.frame(groupsegID[,list(groupID, segID, Ngroup, Npix)])))
}

profoundSegimMerge=function(image=NULL, segim_base=NULL, segim_add=NULL, mask=NULL, sky=0){
  
  segim_add[segim_add>0]=segim_add[segim_add>0]+max(segim_base, na.rm=TRUE)
  
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
  
  invisible(segim_merge)
}

profoundSegimWarp=function(segim_in=NULL, header_in=NULL, header_out=NULL){
  invisible(magwarp(image_in = segim_in, header_out = header_out, header_in = header_in, doscale = FALSE, interpolation = 'nearest')$image)
}

profoundSegimShare=function(segim_in=NULL, header_in=NULL, header_out=NULL, pixcut=1, weights=NULL){
  segID_in=sort(unique(as.integer(segim_in[segim_in>0])))
  segim_warp=profoundSegimWarp(segim_in=segim_in, header_in=header_in, header_out=header_out)
  segim_warp_tab=tabulate(segim_warp)
  segID_warp=which(segim_warp_tab>=pixcut)
  segim_warp[!segim_warp %in% segID_warp]=0
  segim_unwarp=profoundSegimWarp(segim_in=segim_warp, header_in=header_out, header_out=header_in)
  segim_out=NULL
  segimDT=data.table(segim_in=as.integer(segim_in), segim_out=as.integer(segim_unwarp))
  segimDT=segimDT[segim_in>0,]
  segim_groups=segimDT[,list(segim_out=list(tabulate(segim_out))),keyby=segim_in]
  sharemat=matrix(0,max(segim_warp),dim(segim_groups)[1])
  for(i in 1:(dim(sharemat)[2])){sharemat[1:length(unlist(segim_groups$segim_out[i])),i]=unlist(segim_groups$segim_out[i])}
  sharemat=sharemat/rowSums(sharemat)
  sharemat=sharemat[is.finite(sharemat[,1]),,drop=FALSE]
  colnames(sharemat)=segID_in
  rownames(sharemat)=segID_warp
  if(! is.null(weights)){
    t(t(sharemat)*weights)
    sharemat=sharemat/rowSums(sharemat)
  }
  shareseg=diag(sharemat[segID_warp %in% segID_in,segID_in %in% segID_warp])
  invisible(list(segID_in=segID_in, segID_warp=segID_warp, segim_warp=segim_warp, sharemat=sharemat, shareseg=shareseg))
}

profoundShareFlux=function(segstats=NULL, sharemat=NULL, weights=NULL){
  if(dim(segstats)[1] != dim(sharemat)[1]){
    stop('Input segstats and sharemat are not compatible!')
  }
  if(! is.null(weights)){
    t(t(sharemat)*weights)
    sharemat=sharemat/rowSums(sharemat)
  }
  segstats[is.na(segstats)]=0
  flux_out = as.numeric(segstats$flux %*% sharemat)
  flux_err_out = as.numeric(sqrt(segstats$flux_err^2 %*% sharemat))
  mag_out = as.numeric(-2.5*log10(10^(-0.4*segstats$mag) %*% sharemat))
  mag_err_out = as.numeric((2.5/log(10))*abs(flux_err_out/flux_out))
  
  N50_out = as.numeric(segstats$N50 %*% sharemat)
  N90_out = as.numeric(segstats$N90 %*% sharemat)
  N100_out = as.numeric(segstats$N100 %*% sharemat)
  invisible(data.frame(segID=as.integer(colnames(sharemat)), flux=flux_out, flux_err=flux_err_out, mag=mag_out, mag_err=mag_err_out, N50=N50_out, N90=N90_out, N100=N100_out))
}

profoundSegimKeep=function(segim=NULL, groupim=NULL, groupID_merge=NULL, segID_merge=NULL, clean=FALSE){
  segim_out=segim
  
  if(! is.null(groupID_merge)){
    groupID_merge=groupID_merge[groupID_merge %in% groupim]
    if(clean){
      removeID=segim[groupim %in% groupID_merge]
      segim_out[segim_out %in% removeID]=0
    }
    segim_out[groupim %in% groupID_merge]=groupim[groupim %in% groupID_merge]
  }
  
  if(! is.null(segID_merge)){
    if(! is.list(segID_merge)){
      stop('segID_merge must be a list of segments to be merged!')
    }
    whichpix=which(segim_out %in% unlist(segID_merge))
    pixsel=segim_out[whichpix]
    for(i in 1:length(segID_merge)){
      tempID=segID_merge[[i]]
      if('fastmatch' %in% .packages()){
        pixsel[fastmatch::fmatch(pixsel, tempID, nomatch = 0L) > 0L]=min(tempID)
      }else{
        pixsel[pixsel %in% tempID]=min(tempID)
      }
    }
    segim_out[whichpix]=pixsel
  }
  
  invisible(segim_out)
}

profoundSegimExtend=function(image=NULL, segim=NULL, mask=segim, ...){
  if(is.null(image)){stop('Missing image - this is a required input!')}
  if(is.null(segim)){stop('Missing segim - this is a required input!')}
  
  segimadd=profoundProFound(image=image, mask=mask, ...)$segim
  newloc=which(segimadd>0)
  segimadd[newloc]=segimadd[newloc]+max(segim)
  segim=segim+segimadd
  
  invisible(segim)
}

.profoundFluxCalcMin=function(image=NULL, segim=NULL, mask=NULL){
  
  #Set masked things to NA, to be safe:
  
  if(!is.null(mask)){
    image[mask!=0]=NA
  }
  
  segsel=which(segim>0)
  segID=flux=NULL
  tempDT=data.table(segID=as.integer(segim[segsel]),flux=as.numeric(image[segsel]))
  
  output=tempDT[,.fluxcalcmin(flux), by=segID]
  setkey(output, segID)
  
  return(as.data.frame(output))
}
