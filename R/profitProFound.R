.selectCoG=function(diffmat, threshold=1.05){
  tempout={}
  for(i in 1:dim(diffmat)[1]){
    tempsel=which(diffmat[i,]>1 & diffmat[i,]<threshold)+1
    if(length(tempsel)==0){
      if(any(diffmat[i,]<1, na.rm=TRUE)){
        tempsel=min(which(diffmat[i,]<1))
      }else{
        tempsel=which.min(diffmat[i,])+1
        if(length(tempsel)==0){
          tempsel=1
        }
      }
    }else{
      tempsel=min(tempsel)
    }
    tempout=c(tempout, tempsel)
  }
  return=tempout
}

profitProFound=function(image, segim, objects, mask, tolerance=4, ext=2, sigma=1, smooth=TRUE, pixcut=5, skycut=2, SBlim, size=5, shape='disc', iters=6, threshold=1.05, converge='flux', magzero=0, gain=NULL, pixscale=1, sky, skyRMS, redosky=TRUE, redoskysize=21, box=c(100,100), grid=box, type='bilinear', skytype='median', skyRMStype='quanlo', sigmasel=1, doclip=TRUE, shiftloc = TRUE, paddim = TRUE, header, verbose=FALSE, plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE, nearstats=boundstats, offset=1, sortcol="segID", decreasing=FALSE, lowmemory=FALSE, ...){
  if(verbose){message('Running profitProFound:')}
  timestart=proc.time()[3]
  call=match.call()
  if(length(image)>1e6){rembig=TRUE}else{rembig=FALSE}
  
  #Split out image and header parts of input:
  
  if(!missing(image)){
    if(any(names(image)=='imDat') & missing(header)){
      if(verbose){message('Supplied image contains image and header components.')}
      header=image$hdr
      image=image$imDat
    }else if(any(names(image)=='imDat') & !missing(header)){
      if(verbose){message('Supplied image contains image and header but using specified header.')}
      image=image$imDat
    }
    if(any(names(image)=='dat') & missing(header)){
      if(verbose){message('Supplied image contains image and header components.')}
      header=image$hdr[[1]]
      header=data.frame(key=header[,1],value=header[,2], stringsAsFactors = FALSE)
      image=image$dat[[1]]
    }else if(any(names(image)=='dat') & !missing(header)){
      if(verbose){message('Supplied image contains image and header but using specified header.')}
      image=image$dat[[1]]
    }
    if(any(names(image)=='image') & missing(header)){
      if(verbose){message('Supplied image contains image and header components.')}
      header=image$header
      image=image$image
    }else if(any(names(image)=='image') & !missing(header)){
      if(verbose){message('Supplied image contains image and header but using specified header.')}
      image=image$image
    }
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
  
  #Get the pixel scale, if possible and not provided:
  
  if(missing(pixscale) & !missing(header)){
    pixscale=profitGetPixScale(header)
    if(verbose){message(paste('Extracted pixel scale from header provided:',round(pixscale,3),'asec/pixel.'))}
  }
  
  if(missing(objects)){
    if(!missing(segim)){
      objects=segim
      objects[objects != 0] = 1
    }
  }
  
  #Check for user provided sky, and compute if missing:
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  if((hassky==FALSE | hasskyRMS==FALSE) & iters>0){
    if(verbose){message(paste('Making initial sky map -',round(proc.time()[3]-timestart,3),'sec'))}
    roughsky=profitMakeSkyGrid(image=image, objects=objects, mask=mask, box=box, grid=grid, type=type, shiftloc = shiftloc, paddim = paddim)
    if(hassky==FALSE){
      sky=roughsky$sky
    }
    if(hasskyRMS==FALSE){
      skyRMS=roughsky$skyRMS
    }
  }else{
    if(verbose){message("Skipping making initial sky map - User provided sky and sky RMS or iters=0")}
  }
  
  #Make the initial segmentation map, if not provided.
  
  if(missing(segim)){
    if(verbose){message(paste('Making initial segmentation image -',round(proc.time()[3]-timestart,3),'sec'))}
    segim=profitMakeSegim(image=image, objects=objects, mask=mask, tolerance=tolerance, ext=ext, sigma=sigma, smooth=smooth, pixcut=pixcut, skycut=skycut, SBlim=SBlim,  sky=sky, skyRMS=skyRMS, verbose=verbose, plot=FALSE, stats=FALSE)
    objects=segim$objects
    segim=segim$segim
  }else{
    if(verbose){message("Skipping making an initial segmentation image - User provided segim")}
  }
  
  segim_orig=segim
  
  if(any(segim>0)){
    if((hassky==FALSE | hasskyRMS==FALSE) & iters>0){
      if(verbose){message(paste('Doing initial aggressive dilation -',round(proc.time()[3]-timestart,3),'sec'))}
      objects_redo=profitMakeSegimDilate(image=image, segim=objects, mask=mask, size=redoskysize, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$objects
      if(verbose){message(paste('Making better sky map -',round(proc.time()[3]-timestart,3),'sec'))}
      bettersky=profitMakeSkyGrid(image=image, objects=objects_redo, mask=mask, box=box, grid=grid, type=type, shiftloc = shiftloc, paddim = paddim)
      if(hassky==FALSE){
        sky=bettersky$sky
      }
      if(hasskyRMS==FALSE){
        skyRMS=bettersky$skyRMS
      }
    }else{
      if(verbose){message("Skipping making better sky map - User provided sky and sky RMS or iters=0")}
    }
    
    if(iters>0){
      if(verbose){message(paste('Calculating initial segstats -',round(proc.time()[3]-timestart,3),'sec'))}
      segstats=profitSegimStats(image=image, segim=segim, mask=mask, sky=sky)
      compmat=cbind(segstats[,converge])
      segim_array=array(0, dim=c(dim(segim),iters+1))
      segim_array[,,1]=segim
      
      if(verbose){message('Doing dilations:')}
      
      for(i in 1:iters){
        if(verbose){message(paste('Iteration',i,'of',iters,'-',round(proc.time()[3]-timestart,3),'sec'))}
        segim=profitMakeSegimDilate(image=image, segim=segim_array[,,i], mask=mask, size=size, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=TRUE, rotstats=FALSE)
        compmat=cbind(compmat, segim$segstats[,converge])
        segim_array[,,i+1]=segim$segim
      }
      
      if(verbose){message(paste('Finding CoG convergence -',round(proc.time()[3]-timestart,3),'sec'))}
      
      diffmat=rbind(compmat[,2:iters]/compmat[,1:(iters-1)])
      selseg=.selectCoG(diffmat, threshold)
      
      segim=segim$segim
      segim[]=0
      
      if(verbose){message(paste('Constructing final segim -',round(proc.time()[3]-timestart,3),'sec'))}
      for(i in 1:(iters+1)){
        select=segim_array[,,i] %in% segstats[selseg==i,'segID']
        segim[select]=segim_array[,,i][select]
      }
      
      if(rembig){
        rm(select)
        rm(segim_array)
        invisible(gc())
      }
      
      origfrac=compmat[,1]/compmat[cbind(1:length(selseg),selseg)]
      
      objects=segim
      objects[objects!=0]=1
      
    }else{
      if(verbose){message('Iters set to 0 - keeping segim un-dilated')}
      selseg=0
      origfrac=1
    }
    
    if(redosky){
      if(redoskysize %% 2 == 0){redoskysize=redoskysize+1}
      if(verbose){message(paste('Doing final aggressive dilation -',round(proc.time()[3]-timestart,3),'sec'))}
      objects_redo=profitMakeSegimDilate(image=image, segim=objects, mask=mask, size=redoskysize, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$objects
      if(verbose){message(paste('Making final sky map -',round(proc.time()[3]-timestart,3),'sec'))}
      sky=profitMakeSkyGrid(image=image, objects=objects_redo, mask=mask, box=box, grid=grid, type=type, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, doclip=doclip, shiftloc = shiftloc, paddim = paddim)
      skyRMS=sky$skyRMS
      sky=sky$sky
    }else{
      if(verbose){message("Skipping making final sky map - redosky set to FALSE")}
      objects_redo=NULL
    }
    
    if(lowmemory){
      image=image-sky
      sky=0
      skyRMS=0
      segim_orig=NULL
      objects=NULL
      objects_redo=NULL
      invisible(gc())
    }
    
    if(stats & !missing(image)){
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
      segstats=profitSegimStats(image=image, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, pixscale=pixscale, header=header, sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, offset=offset)
      segstats=cbind(segstats, iter=selseg, origfrac=origfrac)
    }else{
      if(verbose){message("Skipping sementation statistics - segstats set to FALSE")}
      segstats=NULL
    }
    
    if(nearstats){
      near=profitSegimNear(segim=segim, offset=offset)
    }else{
      near=NULL
    }
    
    if(plot){
      if(verbose){message(paste('Plotting segments -',round(proc.time()[3]-timestart,3),'sec'))}
      profitSegimPlot(image=image, segim=segim, mask=mask, sky=sky, ...)
    }else{
      if(verbose){message("Skipping segmentation plot - plot set to FALSE")}
    }
    
    if(!missing(SBlim) & !missing(magzero)){
      SBlim=min(SBlim, profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale), na.rm=TRUE)
    }else if(missing(SBlim) & skycut>0){
      SBlim=profitFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
    }else{
      SBlim=NULL
    }
    if(missing(header)){header=NULL}
    if(verbose){message(paste('profitProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
    return=list(segim=segim, segim_orig=segim_orig, objects=objects, objects_redo=objects_redo, sky=sky, skyRMS=skyRMS, segstats=segstats, near=near, header=header, SBlim=SBlim, magzero=magzero, gain=gain, pixscale=pixscale, call=call)
  }else{
    if(missing(header)){header=NULL}
    if(verbose){message('No objects in segmentation map - skipping dilations and CoG')}
    if(verbose){message(paste('profitProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
    return=list(segim=segim, segim_orig=segim_orig, objects=objects, objects_redo=segim, sky=sky, skyRMS=skyRMS, segstats=NULL, near=NULL, header=header, SBlim=NULL,  magzero=magzero, gain=gain, pixscale=pixscale, call=call)
  }
}
