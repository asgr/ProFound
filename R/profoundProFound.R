profoundProFound=function(image=NULL, segim=NULL, objects=NULL, mask=NULL, skycut=1, pixcut=3, tolerance=4, ext=2, reltol=0, cliptol=Inf, sigma=1, smooth=TRUE, SBlim=NULL, size=5, shape='disc', iters=6, threshold=1.05, magzero=0, gain=NULL, pixscale=1, sky=NULL, skyRMS=NULL, redosegim=FALSE, redosky=TRUE, redoskysize=21, box=c(100,100), grid=box, type='bicubic', skytype='median', skyRMStype='quanlo', roughpedestal=FALSE, sigmasel=1, skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, iterskyloc=TRUE, deblend=FALSE, df=3, radtrunc=2, iterative=FALSE, doclip=TRUE, shiftloc = FALSE, paddim = TRUE, header=NULL, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, nearstats=boundstats, groupstats=boundstats, group=NULL, groupby='segim_orig', offset=1, haralickstats=FALSE, sortcol="segID", decreasing=FALSE, lowmemory=FALSE, keepim=TRUE, watershed='ProFound', pixelcov=FALSE, deblendtype='fit', psf=NULL, fluxweight='sum', convtype = 'brute', convmode = 'extended', fluxtype='Raw', app_diam=1, Ndeblendlim=Inf, ...){
  if(verbose){message('Running ProFound:')}
  timestart=proc.time()[3]
  
  call=match.call()
  
  if(length(box)==1){
    box=rep(box,2)
    if(missing(grid)){grid=box}
    if(missing(boxadd)){boxadd=box/2}
    if(missing(skypixmin)){skypixmin=prod(box)/2}
  }
  if(length(grid)==1){
    grid=rep(grid,2)
  }
  if(length(boxadd)==1){
    boxadd=rep(boxadd,2)
  }
  
  fluxtype=tolower(fluxtype)
  
  if(fluxtype=='raw' | fluxtype=='adu' | fluxtype=='adus'){
    if(verbose){message('Using raw flux units')}
    fluxscale=1
  }else if (fluxtype=='jansky'){
    if(verbose){message('Using Jansky flux units (WARNING: magzero must take system to AB)')}
    fluxscale=10^(-0.4*(magzero-8.9))
  }else{
    stop('fluxtype must be Jansky / Raw!')
  }
  
  #Split out image and header parts of input:
  
  if(!is.null(image)){
    if(any(names(image)=='imDat') & is.null(header)){
      if(verbose){message('Supplied image contains image and header components')}
      header=image$hdr
      image=image$imDat
    }else if(any(names(image)=='imDat') & !is.null(header)){
      if(verbose){message('Supplied image contains image and header but using specified header')}
      image=image$imDat
    }
    if(any(names(image)=='dat') & is.null(header)){
      if(verbose){message('Supplied image contains image and header components')}
      header=image$hdr[[1]]
      header=data.frame(key=header[,1],value=header[,2], stringsAsFactors = FALSE)
      image=image$dat[[1]]
    }else if(any(names(image)=='dat') & !is.null(header)){
      if(verbose){message('Supplied image contains image and header but using specified header')}
      image=image$dat[[1]]
    }
    if(any(names(image)=='image') & is.null(header)){
      if(verbose){message('Supplied image contains image and header components')}
      header=image$header
      image=image$image
    }else if(any(names(image)=='image') & !is.null(header)){
      if(verbose){message('Supplied image contains image and header but using specified header')}
      image=image$image
    }
  }else{
    stop('Missing image - this is a required input!')
  }
  
  if(verbose){message(paste('Supplied image is',dim(image)[1],'x',dim(image)[2],'pixels'))}
  
  #Treat image NAs as masked regions:
  
  badpix=NULL
  if(!is.null(mask)){
    mask=mask*1L #Looks silly, but this ensures a logical mask becomes integer.
    if(length(mask)==1 & !is.na(mask[1])){
      maskflag=mask
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      mask[image==maskflag]=1L
    }
    if(anyNA(image)){
      badpix=which(is.na(image))
      mask[badpix]=1L
      image[badpix]=0
    }
  }else{
    if(anyNA(image)){
      mask=matrix(0L,dim(image)[1],dim(image)[2])
      badpix=which(is.na(image))
      mask[badpix]=1L
      image[badpix]=0
    }
  }
  
  #if(!is.null(segim) & !is.null(mask)){
  #  segim=segim*(1-mask) #I don't think we actually need this
  #}
  
  #Get the pixel scale, if possible and not provided:
  
  if(missing(pixscale) & !is.null(header)){
    pixscale=getpixscale(header)
    if(verbose){message(paste('Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel'))}
  }else{
    if(verbose){message(paste('Using suggested pixel scale:',round(pixscale,3),'asec/pixel'))}
  }
  
  skyarea=prod(dim(image))*pixscale^2/(3600^2)
  if(verbose){message(paste('Supplied image is',round(dim(image)[1]*pixscale/60,3),'x',round(dim(image)[2]*pixscale/60,3),'amin, ', round(skyarea,3),'deg-sq'))}
  
  if(is.null(objects)){
    if(!is.null(segim)){
      objects=matrix(0L,dim(segim)[1],dim(segim)[2])
      objects[]=as.logical(segim)
    }
  }else{
    objects=objects*1 #Looks silly, but this ensures a logical mask becomes integer.
  }
  
  #Check for user provided sky, and compute if missing:
  
  hassky=!is.null(sky)
  hasskyRMS=!is.null(skyRMS)
  
  if((hassky==FALSE | hasskyRMS==FALSE) & is.null(segim)){
    if(verbose){message(paste('Making initial sky map -',round(proc.time()[3]-timestart,3),'sec'))}
    roughsky=profoundMakeSkyGrid(image=image, objects=objects, mask=mask, box=box, grid=grid, type=type, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin, boxadd=boxadd, boxiters=0, doclip=doclip, shiftloc=shiftloc, paddim=paddim)
    if(roughpedestal){
      roughsky$sky=median(roughsky$sky,na.rm=TRUE)
      roughsky$skyRMS=median(roughsky$skyRMS,na.rm=TRUE)
    }
    if(hassky==FALSE){
      sky=roughsky$sky
      if(verbose){message(' - Sky statistics :')}
      if(verbose){print(summary(as.numeric(sky)))}
    }
    if(hasskyRMS==FALSE){
      skyRMS=roughsky$skyRMS
      if(verbose){message(' - Sky-RMS statistics :')}
      if(verbose){print(summary(as.numeric(skyRMS)))}
    }
      rm(roughsky)
  }else{
    if(verbose){message("Skipping making initial sky map - User provided sky and sky RMS, or user provided segim")}
  }
  
  #Make the initial segmentation map, if not provided.
  
  if(is.null(segim)){
    if(verbose){message(paste('Making initial segmentation image -',round(proc.time()[3]-timestart,3),'sec'))}
    segim=profoundMakeSegim(image=image, objects=objects, mask=mask, tolerance=tolerance, ext=ext, reltol=reltol, cliptol=cliptol, sigma=sigma, smooth=smooth, pixcut=pixcut, skycut=skycut, SBlim=SBlim,  sky=sky, skyRMS=skyRMS, magzero=magzero, pixscale=pixscale, verbose=verbose, watershed=watershed, plot=FALSE, stats=FALSE)
    objects=segim$objects
    segim=segim$segim
  }else{
    redosegim=FALSE
    # Commented out. New profoundExtendSegim function made instead.
    # if(extendsegim){
    #   if(verbose){message("Finding additional sources - User provided segim")}
    #   segimadd=profoundMakeSegim(image=image, mask=objects, tolerance=tolerance, ext=ext, reltol=reltol, cliptol=cliptol, sigma=sigma, smooth=smooth, pixcut=pixcut, skycut=skycut, SBlim=SBlim,  sky=sky, skyRMS=skyRMS, magzero=magzero, pixscale=pixscale, verbose=verbose, watershed=watershed, plot=FALSE, stats=FALSE)
    #   newloc=which(segimadd$segim>0)
    #   segimadd$segim[newloc]=segimadd$segim[newloc]+max(segim)
    #   segim=segim+segimadd$segim
    #   objects[newloc]=1
    #   rm(newloc)
    #   rm(segimadd)
    # }else{
    if(verbose){message("Skipping making an initial segmentation image - User provided segim")}
    #}
  }
  
  if(any(segim>0)){
    if((hassky==FALSE | hasskyRMS==FALSE)){
      if(redosky){
        if(verbose){message(paste('Doing initial aggressive dilation -',round(proc.time()[3]-timestart,3),'sec'))}
        objects_redo=profoundMakeSegimDilate(segim=objects, size=redoskysize, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$objects
      }else{
        objects_redo=objects
      }
      if(verbose){message(paste('Making better sky map -',round(proc.time()[3]-timestart,3),'sec'))}
      bettersky=profoundMakeSkyGrid(image=image, objects=objects_redo, mask=mask, box=box, grid=grid, type=type, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim)
      if(hassky==FALSE){
        sky=bettersky$sky
        if(verbose){message(' - Sky statistics :')}
        if(verbose){print(summary(as.numeric(sky)))}
      }
      if(hasskyRMS==FALSE){
        skyRMS=bettersky$skyRMS
        if(verbose){message(' - Sky-RMS statistics :')}
        if(verbose){print(summary(as.numeric(skyRMS)))}
      }
      if(redosegim){
        if(verbose){message(paste('Making better segmentation image -',round(proc.time()[3]-timestart,3),'sec'))}
        imagescale=(image-sky)/skyRMS
        imagescale[!is.finite(imagescale)]=0
        if(!is.null(SBlim) & !missing(magzero)){
          imagescale[imagescale<skycut | sky<profoundSB2Flux(SBlim, magzero, pixscale)]=0
        }else{
          imagescale[imagescale<skycut]=0
        }
        if(!is.null(mask)){
          imagescale[mask!=0]=0
        }
        segim[imagescale==0]=0
        objects[segim==0]=0
      }
    }else{
      if(verbose){message("Skipping making better sky map - User provided sky and sky RMS")}
    }
    
    if(iters>0 | iterskyloc){
      if(verbose){message(paste('Calculating initial segstats -',round(proc.time()[3]-timestart,3),'sec'))}
      segstats=.profoundFluxCalcMin(image=image, segim=segim, mask=mask)
      skystats=.profoundFluxCalcMin(image=sky, segim=segim, mask=mask)
      skystats=skystats$flux/skystats$N100
      skymed=median(skystats, na.rm=TRUE)
      origfrac=segstats$flux - (skystats*segstats$N100)
      
      if(iterskyloc){
        localadd=1
      }else{
        localadd=0
      }
      
      #compmat=matrix(0,nrow = dim(segstats)[1], ncol = iters+1+localadd)
      #Nmat=compmat
      #compmat[,1]=segstats[,'flux']
      #Nmat[,1]=segstats[,'N100']
      #flux_old=segstats[,'flux']
      #N100_old=segstats[,'N100']
      
      #segim_array=array(0L, dim=c(dim(segim),iters+1+localadd))
      #segim_array[,,1]=segim
      
      segim_orig=segim
      expand_segID=segstats[,'segID']
      SBlast=rep(Inf,length(expand_segID))
      selseg=rep(0,length(expand_segID))
      
      if(verbose){message('Doing dilations:')}
        
      for(i in 1:(iters)){
        if(verbose){message(paste('Iteration',i,'of',iters,'-',round(proc.time()[3]-timestart,3),'sec'))}
        segim_new=profoundMakeSegimDilate(segim=segim, expand=expand_segID, size=size, shape=shape, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$segim
        segstats_new=.profoundFluxCalcMin(image=image, segim=segim_new, mask=mask)
        SBnew=(segstats_new$flux - segstats$flux) / (segstats_new$N100 - segstats$N100)
        fluxgrowth = (segstats_new$flux - skystats * segstats_new$N100) / (segstats$flux - skystats * segstats$N100) #account for flux growth
        skyfrac = abs(((skystats-skymed) * (segstats_new$N100-segstats$N100)) / (segstats_new$flux - segstats$flux)) #account for sky growth
        expand_segID=segstats[which(segstats_new$flux>0 & fluxgrowth > threshold & SBnew < (SBlast/threshold) & skyfrac < 0.5 & selseg==(i-1)),'segID']
        expand_segID=expand_segID[is.finite(expand_segID)]
        if(length(expand_segID)==0){break}
        updateID=which(segstats$segID %in% expand_segID)
        selseg[updateID] = i
        segstats[updateID,] = segstats_new[updateID,]
        SBlast = SBnew
        if('fastmatch' %in% .packages()){ #dilate segments that pass tests
          selpix = which(fastmatch::fmatch(segim_new, expand_segID, nomatch = 0L) > 0) 
        }else{
          selpix = which(segim_new %in% expand_segID)
        }
        segim[selpix]=segim_new[selpix]
      }
      
      if(iterskyloc){
        segim_skyloc=profoundMakeSegimDilate(segim=segim, size=size, shape=shape, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$segim
        segstats_new=.profoundFluxCalcMin(image=image, segim=segim_skyloc, mask=mask)
        skyseg_mean=(segstats_new$flux-segstats$flux)/(segstats_new$N100-segstats$N100)
        skyseg_mean[!is.finite(skyseg_mean)]=0
      }else{
        skyseg_mean=NA
      }
      
      objects=matrix(0L,dim(segim)[1],dim(segim)[2])
      objects[]=as.logical(segim)
      
      origfrac = origfrac / (segstats$flux - (skystats * segstats$N100))
    }else{
      if(verbose){message('Iters set to 0 - keeping segim un-dilated')}
      segim_orig=segim
      selseg=0
      origfrac=1
    }
    
    if(redosky){
      if(redoskysize %% 2 == 0){redoskysize=redoskysize+1}
      if(verbose){message(paste('Doing final aggressive dilation -',round(proc.time()[3]-timestart,3),'sec'))}
      objects_redo=profoundMakeSegimDilate(segim=objects, mask=mask, size=redoskysize, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$objects
      if(verbose){message(paste('Making final sky map -',round(proc.time()[3]-timestart,3),'sec'))}
      sky=profoundMakeSkyGrid(image=image, objects=objects_redo, mask=mask, box=box, grid=grid, type=type, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim)
      skyRMS=sky$skyRMS
      sky=sky$sky
      if(verbose){message(' - Sky statistics :')}
      if(verbose){print(summary(as.numeric(sky)))}
      if(verbose){message(' - Sky-RMS statistics :')}
      if(verbose){print(summary(as.numeric(skyRMS)))}
    }else{
      if(verbose){message("Skipping making final sky map - redosky set to FALSE")}
      objects_redo=NULL
    }
    
    Norig=tabulate(segim_orig)
    
    if(pixelcov){
      if(verbose){message(paste('Calculating pixel covariance -',round(proc.time()[3]-timestart,3),'sec'))}
      cor_err_func=profoundPixelCorrelation(image=image, objects=objects, mask=mask, sky=sky, skyRMS=skyRMS, fft=FALSE, lag=apply(expand.grid(c(1,2,4),c(1,10,100,1000,1e4)),MARGIN=1,FUN=prod))$cor_err_func
    }else{
      cor_err_func=NULL
    }
    
    if(lowmemory){
      image=image-sky
      sky=0
      skyRMS=0
      segim_orig=NULL
      objects=NULL
      objects_redo=NULL
    }
    
    if(stats & !is.null(image)){
      if(verbose){message(paste('Calculating final segstats for',length(which(tabulate(segim)>0)),'objects -',round(proc.time()[3]-timestart,3),'sec'))}
      if(verbose){message(paste(' - magzero =', round(magzero,3)))}
      if(verbose){
        if(is.null(gain)){
          message(paste(' - gain = NULL (ignored)'))
        }else{
          message(paste(' - gain =', round(gain,3)))
        }
      }
      if(verbose){message(paste(' - pixscale =', round(pixscale,3)))}
      if(verbose){message(paste(' - rotstats =', rotstats))}
      if(verbose){message(paste(' - boundstats =', boundstats))}
      segstats=profoundSegimStats(image=image, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset, cor_err_func=cor_err_func, app_diam=app_diam)
      segstats=cbind(segstats, iter=selseg, origfrac=origfrac, Norig=Norig[segstats$segID], skyseg_mean=skyseg_mean)
      segstats=cbind(segstats, flag_keep=segstats$origfrac>= median(segstats$origfrac[segstats$iter==iters]) | segstats$iter<iters)
    }else{
      if(verbose){message("Skipping segmentation statistics - segstats set to FALSE")}
      segstats=NULL
    }
    
    if(nearstats){
      near=profoundSegimNear(segim=segim, offset=offset)
    }else{
      near=NULL
    }
    
    if(deblend){
      groupstats=TRUE
    }
    
    if(groupstats){
      if(verbose){message(paste(' - groupstats = TRUE - ',round(proc.time()[3]-timestart,3),'sec'))}
      if(groupby=='segim'){
        if(is.null(group)){
          group=profoundSegimGroup(segim)
        }
      }else if(groupby=='segim_orig'){
        if(is.null(group)){
          #message(round(proc.time()[3]-timestart,3))
          group=profoundSegimGroup(segim_orig)
          if(any(group$groupsegID$Ngroup>1)){
            #message(round(proc.time()[3]-timestart,3))
            group$groupim=profoundSegimKeep(segim=segim, segID_merge=group$groupsegID[group$groupsegID$Ngroup>1,'segID'])
            #message(round(proc.time()[3]-timestart,3))
            group$groupsegID$Npix=tabulate(group$groupim)[group$groupsegID$groupID]
          }
        }
      }else{
        stop('Non legal groupby option, must be segim or segim_orig!')
      }

      if(stats & !is.null(image) & !is.null(group)){
        groupstats=profoundSegimStats(image=image, segim=group$groupim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset, cor_err_func=cor_err_func, app_diam=app_diam)
        colnames(groupstats)[1]='groupID'
      }else{
        groupstats=NULL
      }
    }else{
      if(verbose){message(' - groupstats = FALSE')}
      group=NULL
      groupstats=NULL
    }
    
    if(deblend & stats & !is.null(image) & any(group$groupsegID$Ngroup>1)){
      if(verbose){message(paste(' - deblend = TRUE - ',round(proc.time()[3]-timestart,3),'sec'))}
      tempblend=profoundFluxDeblend(image=image-sky, segim=segim, segstats=segstats, groupim=group$groupim, groupsegID=group$groupsegID, magzero=magzero, df=df, radtrunc=radtrunc, iterative=iterative, doallstats=TRUE, deblendtype=deblendtype, psf=psf, fluxweight=fluxweight, convtype=convtype, convmode=convmode, Ndeblendlim = Ndeblendlim)
      if(!is.null(tempblend)){
        segstats=cbind(segstats,tempblend[,-2])
      }
    }else{
      if(verbose){message(' - deblend = FALSE')}
    }
    
    if(haralickstats){
      if(requireNamespace("EBImage", quietly = TRUE)){
        scale=10^(0.4*(30-magzero))
        temphara=(image-sky)*scale
        if(!is.null(mask)){
          temphara[mask!=0]=0
        }
        temphara[!is.finite(temphara)]=0
        haralick=as.data.frame(EBImage::computeFeatures.haralick(segim,temphara))
        haralick=haralick[segstats$segID,]
      }else{
        if(verbose){
          message('The EBImage package is needed to compute Haralick statistics.')
          haralick=NULL
        }
      }
    }else{
      haralick=NULL
    }
    
    if(plot){
      if(verbose){message(paste('Plotting segments -',round(proc.time()[3]-timestart,3),'sec'))}
      if(any(is.finite(sky))){
        profoundSegimPlot(image=image-sky, segim=segim, mask=mask, header=header, ...)
      }else{
        profoundSegimPlot(image=image, segim=segim, mask=mask, header=header, ...)
      }
    }else{
      if(verbose){message("Skipping segmentation plot - plot = FALSE")}
    }
    
    if(is.null(SBlim)){
      SBlim=NULL
    }else if(is.numeric(SBlim)){
      SBlimtemp=profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
      SBlimtemp=matrix(SBlimtemp,dim(skyRMS)[1],dim(skyRMS)[2])
      SBlimtemp[which(SBlimtemp>SBlim)]=SBlim
      SBlim=SBlimtemp
    }else if(SBlim[1]=='get' & skycut> -Inf){
      SBlim=profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
    }
    
    if(is.null(header)){header=NULL}
    if(keepim==FALSE){image=NULL; mask=NULL}
    if(is.null(mask)){mask=NULL}
    if(!is.null(badpix)){image[badpix]=NA}
    row.names(segstats)=NULL
    
    segstats[,grep('flux',colnames(segstats))]=fluxscale*segstats[,grep('flux',colnames(segstats))]
    
    if(verbose){message(paste('ProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
    output=list(segim=segim, segim_orig=segim_orig, objects=objects, objects_redo=objects_redo, sky=sky, skyRMS=skyRMS, image=image, mask=mask, segstats=segstats, Nseg=dim(segstats)[1], near=near, group=group, groupstats=groupstats, haralick=haralick, header=header, SBlim=SBlim, magzero=magzero, dim=dim(segim), pixscale=pixscale, skyarea=skyarea, gain=gain, call=call, date=date(), time=proc.time()[3]-timestart, ProFound.version=packageVersion('ProFound'), R.version=R.version)
  }else{
    if(is.null(header)){header=NULL}
    if(keepim==FALSE){image=NULL; mask=NULL}
    if(is.null(mask)){mask=NULL}
    if(!is.null(badpix)){image[badpix]=NA}
    if(verbose){message('No objects in segmentation map - skipping dilations and CoG')}
    if(verbose){message(paste('ProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
    output=list(segim=NULL, segim_orig=NULL, objects=NULL, objects_redo=NULL, sky=sky, skyRMS=skyRMS, image=image, mask=mask, segstats=NULL, Nseg=0, near=NULL, group=NULL, groupstats=NULL, haralick=NULL, header=header, SBlim=NULL,  magzero=magzero, dim=dim(segim), pixscale=pixscale, skyarea=skyarea, gain=gain, call=call, date=date(), time=proc.time()[3]-timestart, ProFound.version=packageVersion('ProFound'), R.version=R.version)
  }
  class(output)='profound'
  invisible(output)
}

plot.profound=function(x, logR50=TRUE, dmag=0.5, hist='sky', ...){
  
  if(class(x)!='profound'){
    stop('Object class is not of type profound!')
  }
  
  if(is.null(x$image)){
    stop('Missing image!')
  }
  
  if(is.null(x$segim)){
    stop('Missing segmentation map!')
  }
  
  if(is.null(x$sky)){
    x$sky=matrix(0, x$dim[1], x$dim[2])
  }
  if(length(x$sky)==1){
    x$sky=matrix(x$sky, x$dim[1], x$dim[2])
  }
  
  if(is.null(x$skyRMS)){
    x$skyRMS=matrix(1, x$dim[1], x$dim[2])
  }
  if(length(x$skyRMS)==1){
    x$skyRMS=matrix(x$skyRMS, x$dim[1], x$dim[2])
  }
  
  segdiff=x$segim-x$segim_orig
  segdiff[segdiff<0]=0
  
  if(all(x$skyRMS>0, na.rm=TRUE)){
    image = (x$image-x$sky)/x$skyRMS
  }else{
    image = (x$image-x$sky)
  }
  
  if(!is.null(x$mask)){
    image[x$mask==1]=NA
  }
  
  cmap = rev(colorRampPalette(brewer.pal(9,'RdYlBu'))(100))
  maximg = quantile(abs(image[is.finite(image)]), 0.995, na.rm=TRUE)
  stretchscale = 1/median(abs(image), na.rm=TRUE)
  
  layout(matrix(1:9, 3, byrow=TRUE))
  
  if(!is.null(x$header)){
  
    par(mar=c(3.5,3.5,0.5,0.5))
    magimageWCS(image, x$header, stretchscale=stretchscale, locut=-maximg, hicut=maximg, range=c(-1,1), type='num', zlim=c(-1,1), col=cmap)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    
    par(mar=c(3.5,3.5,0.5,0.5))
    magimageWCS(x$segim, x$header, col=c(NA, rainbow(max(x$segim,na.rm=TRUE), end=2/3)), magmap=FALSE)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    abline(v=c(0,dim(x$image)[1]))
    abline(h=c(0,dim(x$image)[2]))
    
    par(mar=c(3.5,3.5,0.5,0.5))
    magimageWCS(image, x$header)
    magimage(segdiff, col=c(NA, rainbow(max(x$segim,na.rm=TRUE), end=2/3)), magmap=FALSE, add=TRUE)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    
    par(mar=c(3.5,3.5,0.5,0.5))
    if(is.null(x$skyarea)){
      skyarea=prod(x$dim)*x$pixscale^2/(3600^2)
    }else{
      skyarea=x$skyarea
    }
    temphist=maghist(x$segstats$mag, log='y', scale=(2*dmag)/x$skyarea, breaks=seq(floor(min(x$segstats$mag, na.rm = TRUE)), ceiling(max(x$segstats$mag, na.rm = TRUE)),by=0.5), xlab='mag', ylab=paste('#/deg-sq/d',dmag,'mag',sep=''), grid=TRUE, verbose=FALSE)
    #magplot(temphist, log='y', xlab='mag', ylab=expression('#'/'deg-sq'/'dmag'), grid=TRUE)
    ymax=log10(max(temphist$counts,na.rm = T))
    xmax=temphist$mids[which.max(temphist$counts)]
    abline(ymax - xmax*0.4, 0.4, col='red')
    abline(v=xmax+0.25, col='red')
    axis(side=1, at=xmax+0.25, labels=xmax+0.25, tick=FALSE, line=-1, col.axis='red')
      
    par(mar=c(3.5,3.5,0.5,0.5))
    #stretchscale = 1/median(abs(x$sky), na.rm=TRUE)
    #maxsky = quantile(abs(x$sky[is.finite(x$sky)]), 0.995, na.rm=TRUE)
    #magimageWCS(x$sky, x$header, locut=-maxsky, hicut=maxsky, range=c(-1,1), type='num', zlim=c(-1,1), stretchscale=stretchscale, col=cmap)
    magimageWCS(x$sky, x$header, qdiff=TRUE)
    legend('topleft',legend='sky',bg='white')
    
    par(mar=c(3.5,3.5,0.5,0.5))
    magimageWCS(x$skyRMS, x$header)
    legend('topleft',legend='skyRMS',bg='white')
    
    if(hist=='iters'){
      maghist(x$segstats$iter, breaks=seq(-0.5,max(x$segstats$iter, na.rm=TRUE)+0.5,by=1), majorn=max(x$segstats$iter, na.rm=TRUE)+1, xlab='Number of Dilations', ylab='#', verbose=FALSE)
    }else if(hist=='sky'){
      try({
        if(!is.null(x$objects_redo)){
          tempsky=image[x$objects_redo==0]
        }else{
          tempsky=image[x$objects==0]
        }
        tempsky=tempsky[tempsky> -6 & tempsky<6]
        magplot(density(tempsky[is.finite(tempsky)], bw=0.1), grid=TRUE, xlim=c(-5,5), xlab='(image - sky) / skyRMS', ylab='PDF', log='y', ylim=c(1e-5,0.5))
        curve(dnorm(x, mean=0, sd=1), add=TRUE, col='red', lty=2)
        legend('topleft',legend='sky pixels',bg='white')
        })
    }else{stop('Not a recognised hist type! Must be iters / sky.')}
    
    par(mar=c(3.5,3.5,0.5,0.5))
    if(logR50){
      magplot(x$segstats$mag, x$segstats$R50, pch='.', col=hsv(alpha=0.5), ylim=c(min(x$segstats$R50, 0.1, na.rm = TRUE), max(x$segstats$R50, 1, na.rm = TRUE)), cex=3, xlab='mag', ylab='R50 / asec', grid=TRUE, log='y')
    }else{
      magplot(x$segstats$mag, x$segstats$R50, pch='.', col=hsv(alpha=0.5), ylim=c(0, max(x$segstats$R50, 1, na.rm = TRUE)), cex=3, xlab='mag', ylab='R50 / asec', grid=TRUE)
    }
    
    par(mar=c(3.5,3.5,0.5,0.5))
    fluxrat=x$segstats$flux/x$segstats$flux_err
    magplot(x$segstats$SB_N90, fluxrat, pch='.', col=hsv(alpha=0.5), ylim=c(0.5,max(fluxrat, 1, na.rm=TRUE)), cex=3, xlab='SB90 / mag/asec-sq', ylab='Flux/Flux-Error', grid=TRUE, log='y')
  
  }else{
    
    par(mar=c(3.5,3.5,0.5,0.5))
    magimage(image, stretchscale=stretchscale, locut=-maximg, hicut=maximg, range=c(-1,1), type='num', zlim=c(-1,1), col=cmap)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    
    par(mar=c(3.5,3.5,0.5,0.5))
    magimage(x$segim, col=c(NA, rainbow(max(x$segim,na.rm=TRUE), end=2/3)), magmap=FALSE)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    abline(v=c(0,dim(image)[1]))
    abline(h=c(0,dim(image)[2]))
    
    par(mar=c(3.5,3.5,0.5,0.5))
    magimage(image)
    magimage(segdiff, col=c(NA, rainbow(max(x$segim,na.rm=TRUE), end=2/3)), magmap=FALSE, add=TRUE)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}

    par(mar=c(3.5,3.5,0.5,0.5))
    temphist=maghist(x$segstats$mag, log='y', scale=(2*dmag), breaks=seq(floor(min(x$segstats$mag, na.rm = TRUE)), ceiling(max(x$segstats$mag, na.rm = TRUE)),by=0.5), xlab='mag', ylab=paste('#/d',dmag,'mag',sep=''), grid=TRUE, verbose=FALSE)
    ymax=log10(max(temphist$counts,na.rm = T))
    xmax=temphist$mids[which.max(temphist$counts)]
    abline(ymax - xmax*0.4, 0.4, col='red')
    abline(v=xmax+0.25, col='red')
    axis(side=1, at=xmax+0.25, labels=xmax+0.25, tick=FALSE, line=-1, col.axis='red')
    
    par(mar=c(3.5,3.5,0.5,0.5))
    #stretchscale = 1/median(abs(x$sky), na.rm=TRUE)
    #maxsky = quantile(abs(x$sky[is.finite(x$sky)]), 0.995, na.rm=TRUE)
    #magimage(x$sky, locut=-maxsky, hicut=maxsky, range=c(-1,1), type='num', zlim=c(-1,1), col=cmap)
    magimage(x$sky, qdiff=TRUE)
    legend('topleft',legend='sky',bg='white')
    
    par(mar=c(3.5,3.5,0.5,0.5))
    magimage(x$skyRMS)
    legend('topleft',legend='skyRMS',bg='white')
    
    if(hist=='iters'){
      maghist(x$segstats$iter, breaks=seq(-0.5,max(x$segstats$iter, na.rm=TRUE)+0.5,by=1), majorn=max(x$segstats$iter, na.rm=TRUE)+1, xlab='Number of Dilations', ylab='#', verbose=FALSE)
    }else if(hist=='sky'){
      try({
        if(!is.null(x$objects_redo)){
          tempsky=image[x$objects_redo==0]
        }else{
          tempsky=image[x$objects==0]
        }
        tempsky=tempsky[tempsky> -6 & tempsky<6]
        magplot(density(tempsky[is.finite(tempsky)], bw=0.1), grid=TRUE, xlim=c(-5,5), xlab='(image - sky) / skyRMS', ylab='PDF', log='y', ylim=c(1e-5,0.5))
        curve(dnorm(x, mean=0, sd=1), add=TRUE, col='red', lty=2)
        legend('topleft',legend='sky pixels',bg='white')
        })
    }else{stop('Not a recognised hist type! Must be iters / sky.')}
    
    par(mar=c(3.5,3.5,0.5,0.5))
    if(logR50){
      magplot(x$segstats$mag, x$segstats$R50, pch='.', col=hsv(alpha=0.5), ylim=c(min(x$segstats$R50, 0.1, na.rm = TRUE), max(x$segstats$R50, 1, na.rm = TRUE)), cex=3, xlab='mag', ylab='R50 / asec', grid=TRUE, log='y')
    }else{
      magplot(x$segstats$mag, x$segstats$R50, pch='.', col=hsv(alpha=0.5), ylim=c(0, max(x$segstats$R50, 1, na.rm = TRUE)), cex=3, xlab='mag', ylab='R50 / asec', grid=TRUE)
    }
    
    par(mar=c(3.5,3.5,0.5,0.5))
    fluxrat=x$segstats$flux/x$segstats$flux_err
    magplot(x$segstats$SB_N90, fluxrat, pch='.', col=hsv(alpha=0.5), ylim=c(0.5,max(fluxrat, 1, na.rm=TRUE)), cex=3, xlab='SB90 / mag/pix-sq', ylab='Flux/Flux-Error', grid=TRUE, log='y')
  }
  
}

.selectCoG=function(diffmat, threshold=1.05){
  IDmat=matrix(rep(1:dim(diffmat)[2],each=dim(diffmat)[1]),nrow=dim(diffmat)[1])
  logmat=diffmat>1 & diffmat<threshold
  IDfin=IDmat
  IDfin[logmat==FALSE]=NA
  NegFlux=which(diffmat<threshold^0.2,arr.ind=TRUE)
  if(length(NegFlux)>0){
    NegFlux[,2]=NegFlux[,2]-1
    IDfin[NegFlux]=IDmat[NegFlux]
    IDfin[NegFlux[NegFlux[,2]==0,1],1]=0
  }
  tempout=suppressWarnings(apply(IDfin,1,min,na.rm=TRUE))
  tempout[is.infinite(tempout)]=dim(diffmat)[2]
  tempout+1
}
