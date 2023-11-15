profoundProFound=function(image=NULL, segim=NULL, objects=NULL, mask=NULL, skycut=1.5, pixcut=3,
                          tolerance=4, ext=2, reltol=0, cliptol=Inf, sigma=1, smooth=TRUE,
                          SBlim=NULL, SBdilate=2, SBN100=100, size=5, shape='disc', iters=6,
                          threshold=1.05, magzero=0, gain=NULL, pixscale=1, sky=NULL, skyRMS=NULL,
                          redosegim=FALSE, redosky=TRUE, redoskysize=21, box=c(100,100), grid=box,
                          skygrid_type = 'new', type='bicubic', skytype='median', skyRMStype='quanlo', roughpedestal=FALSE,
                          sigmasel=1, skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, conviters=100, iterskyloc=FALSE,
                          deblend=FALSE, df=3, radtrunc=2, iterative=FALSE, doChiSq=FALSE, doclip=TRUE,
                          shiftloc = FALSE, paddim = TRUE, keyvalues=NULL, verbose=FALSE,
                          plot=FALSE, stats=TRUE, rotstats=FALSE, boundstats=FALSE,
                          nearstats=boundstats, groupstats=boundstats, group=NULL,
                          groupby='segim_orig', offset=1, haralickstats=FALSE, sortcol="segID",
                          decreasing=FALSE, lowmemory=FALSE, keepim=TRUE, watershed='ProFound',
                          pixelcov=FALSE, deblendtype='fit', psf=NULL, fluxweight='sum',
                          convtype = 'brute', convmode = 'extended', fluxtype='Raw',
                          app_diam=1, Ndeblendlim=Inf, static_photom = FALSE, rem_mask=FALSE, ...){

  if(verbose){message('Running ProFound:')}
  timestart=proc.time()[3]
  
  call=match.call(expand.dots=TRUE)
  
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
  
  fluxtype = tolower(fluxtype)
  
  if(fluxtype=='raw' | fluxtype=='adu' | fluxtype=='adus'){
    if(verbose){message('Using raw flux units')}
    fluxscale=1
  }else if (fluxtype=='jansky'){
    if(verbose){message('Using Jansky flux units (WARNING: magzero must take system to AB)')}
    fluxscale=10^(-0.4*(magzero-8.9))
  }else if (fluxtype=='microjansky'){
    if(verbose){message('Using Micro-Jansky flux units (WARNING: magzero must take system to AB)')}
    fluxscale=10^(-0.4*(magzero-23.9))
  }else{
    stop('fluxtype must be Jansky / Microjansky / Raw!')
  }
  
  #Split out image and keyvalues parts of input:
  
  if(!is.null(image)){
    if(inherits(image, 'Rfits_image')){
      keyvalues = image$keyvalues
      image = image$imDat
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
  
  if(box[1] > ceiling(dim(image)[1]/3)){
    box[1] = ceiling(dim(image)[1]/3)
    message('dim(image)[1]/box[1] must be >=3, box[1] modified to ',box[1])
  }
  if(box[2] > ceiling(dim(image)[2]/3)){
    box[2] = ceiling(dim(image)[2]/3)
    message('dim(image)[2]/box[2] must be >=3, box[2] modified to ',box[2])
  }
  
  if(grid[1] > ceiling(dim(image)[1]/3)){
    grid[1] = ceiling(dim(image)[1]/3)
    message('dim(image)[1]/grid[1] must be >=3, grid[1] modified to ',grid[1])
  }
  if(grid[2] > ceiling(dim(image)[2]/3)){
    grid[2] = ceiling(dim(image)[2]/3)
    message('dim(image)[2]/grid[2] must be >=3, grid[2] modified to ',grid[2])
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
    if(is.numeric(mask)){
      mask=as.integer(mask)
      attributes(mask)$dim = dim(image)
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
  
  if(missing(pixscale) & !is.null(keyvalues)){
    pixscale = pixscale(keyvalues)
    if(verbose){message(paste('Extracted pixel scale from keyvalues provided:',signif(pixscale,4),'asec/pixel'))}
  }else{
    if(verbose){message(paste('Using suggested pixel scale:',signif(pixscale,4),'asec/pixel'))}
  }
  
  imarea = prod(dim(image))*pixscale^2/(3600^2)
  if(verbose){message(paste('Supplied image is',signif(dim(image)[1]*pixscale/60,4),'x',signif(dim(image)[2]*pixscale/60,4),'amin, ', signif(imarea,4),'deg-sq'))}
  
  if(is.null(objects)){
    if(!is.null(segim)){
      objects = matrix(0L,dim(segim)[1],dim(segim)[2])
      objects[] = as.logical(segim)
    }
  }#else{
  #   objects = objects*1L #Looks silly, but this ensures a logical mask becomes integer.
  # }
  
  if(static_photom & !is.null(segim)){
    #set all the flags to make static photom
    sky = 0
    iters = 0
    redosky = FALSE
    redosegim = FALSE
  }
  
  #Check for user provided sky, and compute if missing:
  
  hassky = !is.null(sky)
  hasskyRMS = !is.null(skyRMS)
  
  if((hassky==FALSE | hasskyRMS==FALSE) & is.null(segim)){
    if(hassky==FALSE){
      sky = 0
    }
    if(verbose){message(paste('Making initial sky map -',round(proc.time()[3]-timestart,3),'sec'))}
    roughsky = profoundMakeSkyGrid(image=image, objects=objects, mask=mask, sky=sky, box=box, 
                                 grid=grid, skygrid_type=skygrid_type, type=type, 
                                 skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, 
                                 skypixmin=skypixmin, boxadd=boxadd, boxiters=0, conviters=conviters,
                                 doclip=doclip, shiftloc=shiftloc, paddim=paddim, rem_mask=rem_mask)
    if(roughpedestal){
      roughsky$sky = median(roughsky$sky,na.rm=TRUE)
      roughsky$skyRMS = median(roughsky$skyRMS,na.rm=TRUE)
    }
    if(hassky==FALSE){
      sky = roughsky$sky
      if(verbose){message(' - Sky statistics :')}
      if(verbose){print(summary(as.numeric(sky[!is.na(sky)])))}
    }
    if(hasskyRMS==FALSE){
      skyRMS = roughsky$skyRMS
      if(verbose){message(' - Sky-RMS statistics :')}
      if(verbose){print(summary(as.numeric(skyRMS[!is.na(skyRMS)])))}
    }
    rm(roughsky)
  }else{
    if(verbose){message("Skipping making initial sky map - User provided sky and sky RMS, or user provided segim")}
  }
  
  #Make the initial segmentation map, if not provided.
  
  if(is.null(segim)){
    if(verbose){message(paste('Making initial segmentation image -',round(proc.time()[3]-timestart,3),'sec'))}
    segim = profoundMakeSegim(image=image, objects=objects, mask=mask, tolerance=tolerance, 
                            ext=ext, reltol=reltol, cliptol=cliptol, sigma=sigma, 
                            smooth=smooth, pixcut=pixcut, skycut=skycut, SBlim=SBlim,  
                            sky=sky, skyRMS=skyRMS, magzero=magzero, pixscale=pixscale, 
                            verbose=verbose, watershed=watershed, plot=FALSE, stats=FALSE)
    objects = segim$objects
    segim = segim$segim
  }else{
    redosegim = FALSE
    if(verbose){message("Skipping making an initial segmentation image - User provided segim")}
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
      if(hasskyRMS==FALSE){rm(skyRMS)}
      bettersky = profoundMakeSkyGrid(image=image, objects=objects_redo, mask=mask, sky=sky,
                                    box=box, skygrid_type=skygrid_type, grid=grid, 
                                    type=type, skytype=skytype, skyRMStype=skyRMStype, 
                                    sigmasel=sigmasel, skypixmin=skypixmin, boxadd=boxadd, 
                                    boxiters=boxiters, conviters=conviters, doclip=doclip, shiftloc=shiftloc, 
                                    paddim=paddim, rem_mask=rem_mask)
      if(hassky==FALSE){
        sky = bettersky$sky
        if(verbose){message(' - Sky statistics :')}
        if(verbose){print(summary(as.numeric(sky[!is.na(sky)])))}
      }
      if(hasskyRMS==FALSE){
        skyRMS = bettersky$skyRMS
        if(verbose){message(' - Sky-RMS statistics :')}
        if(verbose){print(summary(as.numeric(skyRMS[!is.na(skyRMS)])))}
      }
      rm(bettersky)
      if(redosegim){
        if(verbose){message(paste('Making better segmentation image -',round(proc.time()[3]-timestart,3),'sec'))}
        imagescale = (image-sky)/skyRMS
        imagescale[!is.finite(imagescale)] = 0
        if(!is.null(SBlim) & !missing(magzero)){
          imagescale[imagescale<skycut | sky<profoundSB2Flux(SBlim, magzero, pixscale)]=0
        }else{
          imagescale[imagescale<skycut] = 0
        }
        if(!is.null(mask)){
          imagescale[mask!=0] = 0
        }
        segim[imagescale==0] = 0
        objects[segim==0] = 0
        rm(imagescale)
      }
    }else{
      if(verbose){message("Skipping making better sky map - User provided sky and sky RMS")}
    }
    
    if(iters>0 | iterskyloc){
      if(verbose){message(paste('Calculating initial segstats -',round(proc.time()[3]-timestart,3),'sec'))}
      # if(length(sky)>1){
      #   skystats = .profoundFluxCalcMin(image = sky - median(sky,na.rm=TRUE), segim=segim, mask=mask) #run on sky
      #   skystats = skystats$flux / skystats$N100 #get per pixel mean flux per segment
      #   #skymed = median(skystats, na.rm=TRUE) #median per pixel mean flux per segment for the whole image
      # }else{
      #   skystats = 0 #specified
      #   #skymed = sky #specified
      # }
      
      image_sky = image - sky
      segstats = .profoundFluxCalcMin(image=image_sky, segim=segim, mask=mask) #run on initial image
      #segstats$flux = segstats$flux - (skystats * segstats$N100) #remove the local sky component
      origfrac = segstats$flux #set to initial fluxes
      
      segim_orig = segim
      expand_segID = segstats[,'segID']
      SBlast = rep(Inf,length(expand_segID))
      selseg = rep(0,length(expand_segID))
      
      if(is.null(SBdilate)){
        SBdilate = Inf
      }else{
        skyRMSstats = .profoundFluxCalcMin(image=skyRMS, segim=segim, mask=mask) #run on sky
        skyRMSstats = skyRMSstats$flux/skyRMSstats$N100 #get per pixel mean flux per segment
        SBdilate = skyRMSstats * 10^(-0.4*SBdilate)
        SBdilate[!is.finite(SBdilate)]=Inf
      }
      
      if(iters>0){
        if(verbose){message('Doing dilations:')}
        for(i in 1:(iters)){
          if(verbose){message(paste('Iteration',i,'of',iters,'-',round(proc.time()[3]-timestart,3),'sec'))}
          segim_new = profoundDilate(segim=segim, size=size, shape=shape, expand=expand_segID, iters=1)
          #segim_new = profoundMakeSegimDilate(segim=segim, expand=expand_segID, size=size, shape=shape, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$segim #dilate
          segstats_new = .profoundFluxCalcMin(image=image_sky, segim=segim_new, mask=mask) #run on image with dilated segments
          N100diff = (segstats_new$N100 - segstats$N100)
          SBnew = (segstats_new$flux - segstats$flux) / N100diff #calculate the surface brightness of the new grown anulus only
          fluxgrowthratio = segstats_new$flux / segstats$flux #calculate flux growth ratio
          #skyfrac = abs((skystats * N100diff) / (segstats_new$flux - segstats$flux)) #estimate how much of the flux growth might be coming from unusual local sky
          expand_segID = segstats[which((fluxgrowthratio > threshold | (SBnew > SBdilate & N100diff>SBN100)) & segstats_new$flux>0 & SBnew < SBlast & selseg==(i-1)),'segID']
          expand_segID = expand_segID[is.finite(expand_segID)] #safety first
          if(length(expand_segID)==0){break}
          updateID = which(segstats$segID %in% expand_segID)
          selseg[updateID] = i #iteration number for segments flagged for growth
          segstats[updateID,] = segstats_new[updateID,] #update only the segstats for things that are growing, the rest can be ignored
          SBlast = SBnew
          if('fastmatch' %in% .packages()){ #dilate segments that pass tests
            selpix = which(fastmatch::fmatch(segim_new, expand_segID, nomatch = 0L) > 0) 
          }else{
            selpix = which(segim_new %in% expand_segID)
          }
          segim[selpix] = segim_new[selpix]
        }
      }
      
      origfrac = origfrac / segstats$flux
      
      if(iterskyloc){
        segim_skyloc = profoundMakeSegimDilate(segim=segim, size=size, shape=shape, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$segim
        segstats_sky = .profoundFluxCalcMin(image=image, segim=segim_skyloc-segim, mask=mask)
        #segstats$flux = segstats$flux + (skystats * segstats$N100) #add back the local sky component
        skyseg_temp = segstats_sky$flux / segstats_sky$N100
        skyseg_temp[!is.finite(skyseg_temp)] = 0
        skyseg_mean = rep(NA,length(segstats$segID))
        skyseg_mean[match(segstats_sky$segID,segstats$segID)] = skyseg_temp
      }else{
        skyseg_mean = NA
      }
      
      objects = matrix(0L,dim(segim)[1],dim(segim)[2])
      objects[] = as.logical(segim)
      
    }else{
      if(verbose){message('Iters set to 0 - keeping segim un-dilated')}
      segim_orig = segim
      selseg = 0
      origfrac = 1
      skyseg_mean = NA
    }
    
    objects_redo = profoundMakeSegimDilate(segim=objects, mask=mask, size=redoskysize, shape=shape, sky=sky, verbose=verbose, plot=FALSE, stats=FALSE, rotstats=FALSE)$objects
    
    if(redosky){
      if(redoskysize %% 2 == 0){redoskysize=redoskysize+1}
      if(verbose){message(paste('Doing final aggressive dilation -',round(proc.time()[3]-timestart,3),'sec'))}
      
      if(verbose){message(paste('Making final sky map -',round(proc.time()[3]-timestart,3),'sec'))}
      rm(skyRMS)
      sky = profoundMakeSkyGrid(image=image, objects=objects_redo, mask=mask, sky=sky, box=box, 
                              skygrid_type=skygrid_type, grid=grid, type=type, skytype=skytype, 
                              skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin, 
                              boxadd=boxadd, boxiters=boxiters, conviters=conviters, doclip=doclip, shiftloc=shiftloc, 
                              paddim=paddim, rem_mask=rem_mask)$sky
      if(verbose){message(paste('Making final sky RMS map -',round(proc.time()[3]-timestart,3),'sec'))}
      skyRMS = profoundMakeSkyGrid(image=image, objects=objects_redo, mask=mask, sky=sky, box=box, 
                                skygrid_type=skygrid_type, grid=grid, type=type, skytype=skytype, 
                                skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin, 
                                boxadd=boxadd, boxiters=boxiters, conviters=conviters, doclip=doclip, shiftloc=shiftloc, 
                                paddim=paddim, rem_mask=rem_mask)$skyRMS
      if(verbose){message(' - Sky statistics :')}
      if(verbose){print(summary(as.numeric(sky[!is.na(sky)])))}
      if(verbose){message(' - Sky-RMS statistics :')}
      if(verbose){print(summary(as.numeric(skyRMS[!is.na(skyRMS)])))}
    }else{
      if(verbose){message("Skipping making final sky map - redosky set to FALSE")}
    }
    
    if(doChiSq){
      if(verbose){message(paste('Making sky LL ChiSq map -',round(proc.time()[3]-timestart,3),'sec'))}
      skyChiSqMap = profoundMakeSkyGrid(image=image, objects=objects_redo, mask=mask, sky=sky, box=box, 
                                        skygrid_type='old', grid=grid, type=type, skytype=skytype, 
                                        skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin, 
                                        boxadd=boxadd, boxiters=boxiters, conviters=conviters, doChiSq=TRUE,
                                        doclip=doclip, shiftloc=shiftloc, paddim=paddim, rem_mask=rem_mask)$skyChiSq
    }else{
      skyChiSqMap = NA
    }
    
    Norig = tabulate(segim_orig)
    
    if(is.function(pixelcov)){
      cor_err_func = pixelcov
    }else{
      if(pixelcov){
        if(verbose){message(paste('Calculating pixel covariance -',round(proc.time()[3]-timestart,3),'sec'))}
        if(is.null(profoundMakeSkyGrid)){
          cor_err_func = profoundPixelCorrelation(image=image, objects=objects, mask=mask, sky=sky, skyRMS=skyRMS, 
                                              fft=FALSE, lag=apply(expand.grid(c(1,2,4),c(1,10,100,1000,1e4)),MARGIN=1,FUN=prod))$cor_err_func
        }else{
          cor_err_func = profoundPixelCorrelation(image=image, objects=objects_redo, mask=mask, sky=sky, skyRMS=skyRMS, 
                                                  fft=FALSE, lag=apply(expand.grid(c(1,2,4),c(1,10,100,1000,1e4)),MARGIN=1,FUN=prod))$cor_err_func
        }
      }else{
        cor_err_func = NULL
      }
    }
    
    if(lowmemory){
      image = image - sky
      sky = 0
      skyRMS = 0
      segim_orig = NULL
      objects = NULL
      objects_redo = NULL
    }
    
    if(stats & !is.null(image)){
      if(verbose){message(paste('Calculating final segstats for',length(which(tabulate(segim)>0)),'objects -',round(proc.time()[3]-timestart,3),'sec'))}
      if(verbose){message(paste(' - magzero =', signif(magzero,4)))}
      if(verbose){
        if(is.null(gain)){
          message(paste(' - gain = NULL (ignored)'))
        }else{
          message(paste(' - gain =', signif(gain,4)))
        }
      }
      if(verbose){message(paste(' - pixscale =', signif(pixscale,4)))}
      if(verbose){message(paste(' - rotstats =', rotstats))}
      if(verbose){message(paste(' - boundstats =', boundstats))}
      segstats = profoundSegimStats(image=image, segim=segim, mask=mask, sky=sky, skyRMS=skyRMS, 
                                  magzero=magzero, gain=gain, pixscale=pixscale, keyvalues=keyvalues, 
                                  sortcol=sortcol, decreasing=decreasing, rotstats=rotstats, 
                                  boundstats=boundstats, offset=offset, cor_err_func=cor_err_func, 
                                  app_diam=app_diam)
      segstats = cbind(segstats, iter=selseg, origfrac=origfrac, Norig=Norig[segstats$segID], skyseg_mean=skyseg_mean)
      segstats = cbind(segstats, flag_keep=segstats$origfrac >= median(segstats$origfrac[segstats$iter==iters]) | segstats$iter<iters)
    }else{
      if(verbose){message("Skipping segmentation statistics - segstats set to FALSE")}
      segstats = NULL
    }
    
    if(nearstats){
      near = profoundSegimNear(segim=segim, offset=offset)
    }else{
      near = NULL
    }
    
    if(deblend){
      groupstats = TRUE
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
          group = profoundSegimGroup(segim_orig)
          if(any(group$groupsegID$Ngroup>1)){
            #message(round(proc.time()[3]-timestart,3))
            group$groupim = profoundSegimKeep(segim=segim, segID_merge=group$groupsegID[group$groupsegID$Ngroup>1,'segID'])
            #message(round(proc.time()[3]-timestart,3))
            group$groupsegID$Npix = tabulate(group$groupim)[group$groupsegID$groupID]
          }
        }
      }else{
        stop('Non legal groupby option, must be segim or segim_orig!')
      }

      if(stats & !is.null(image) & !is.null(group)){
        groupstats = profoundSegimStats(image=image, segim=group$groupim, mask=mask, 
                                      sky=sky, skyRMS=skyRMS, magzero=magzero, gain=gain, 
                                      pixscale=pixscale, keyvalues=keyvalues, sortcol=sortcol, 
                                      decreasing=decreasing, rotstats=rotstats, boundstats=boundstats, 
                                      offset=offset, cor_err_func=cor_err_func, app_diam=app_diam)
        colnames(groupstats)[1] = 'groupID'
      }else{
        groupstats = NULL
      }
    }else{
      if(verbose){message(' - groupstats = FALSE')}
      group = NULL
      groupstats = NULL
    }
    
    if(deblend & stats & !is.null(image) & any(group$groupsegID$Ngroup>1)){
      if(verbose){message(paste(' - deblend = TRUE - ',round(proc.time()[3]-timestart,3),'sec'))}
      tempblend = profoundFluxDeblend(image=image-sky, segim=segim, segstats=segstats, 
                                    groupim=group$groupim, groupsegID=group$groupsegID, 
                                    magzero=magzero, df=df, radtrunc=radtrunc, iterative=iterative, 
                                    doallstats=TRUE, deblendtype=deblendtype, psf=psf, 
                                    fluxweight=fluxweight, convtype=convtype, convmode=convmode, 
                                    Ndeblendlim = Ndeblendlim)
      if(!is.null(tempblend)){
        segstats = cbind(segstats,tempblend[,-2])
      }
    }else{
      if(verbose){message(' - deblend = FALSE')}
    }
    
    if(haralickstats){
      if(requireNamespace("EBImage", quietly = TRUE)){
        scale = 10^(0.4*(30-magzero))
        temphara = (image-sky)*scale
        if(!is.null(mask)){
          temphara[mask!=0] = 0
        }
        temphara[!is.finite(temphara)] = 0
        haralick = as.data.frame(EBImage::computeFeatures.haralick(segim,temphara))
        haralick = haralick[segstats$segID,]
      }else{
        if(verbose){
          message('The EBImage package is needed to compute Haralick statistics.')
          haralick = NULL
        }
      }
    }else{
      haralick = NULL
    }
    
    if(plot){
      if(verbose){message(paste('Plotting segments -',round(proc.time()[3]-timestart,3),'sec'))}
      if(any(is.finite(sky))){
        profoundSegimPlot(image=image-sky, segim=segim, mask=mask, keyvalues=keyvalues, ...)
      }else{
        profoundSegimPlot(image=image, segim=segim, mask=mask, keyvalues=keyvalues, ...)
      }
    }else{
      if(verbose){message("Skipping segmentation plot - plot = FALSE")}
    }
    
    if(is.null(SBlim)){
      SBlim = NULL
    }else if(is.numeric(SBlim)){
      SBlimtemp = profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
      SBlimtemp = matrix(SBlimtemp,dim(skyRMS)[1],dim(skyRMS)[2])
      SBlimtemp[which(SBlimtemp>SBlim)] = SBlim
      SBlim = SBlimtemp
    }else if(SBlim[1]=='get' & skycut> -Inf){
      SBlim = profoundFlux2SB(flux=skyRMS*skycut, magzero=magzero, pixscale=pixscale)
    }
    
    if(is.null(keyvalues)){keyvalues = NULL}
    if(keepim==FALSE){image = NULL; mask = NULL}
    if(is.null(mask)){mask = NULL}
    if(!is.null(badpix)){image[badpix] = NA}
    row.names(segstats) = NULL
    
    segstats[,grep('flux',colnames(segstats))] = fluxscale*segstats[,grep('flux',colnames(segstats))]
    
    if(!is.null(groupstats)){
      groupstats[,grep('flux',colnames(groupstats))] = fluxscale*groupstats[,grep('flux',colnames(groupstats))]
    }
    
    if(stats){
      cutsky = (image - sky) / skyRMS
      cutsky[!is.finite(cutsky)] = NA
      if(!is.null(mask)){
        cutsky[mask>0] = NA
      }
      if(!is.null(objects_redo)){
        cutsky[objects_redo == 1] = NA
      }else{
        cutsky[objects == 1] = NA
      }
      cutsky = cutsky[which(cutsky<=0)] #Only consider sky pixels below the mode
      
      if(length(cutsky) > 0){
        df = length(cutsky) - 1
        skyChiSq = sum(cutsky^2, na.rm =TRUE)
        skyLL = dchisq(skyChiSq, df=df, log=TRUE)
        skyChiSq = skyChiSq / df
      }else{
        skyLL = NA
        skyChiSq = NA
      }
    }else{
      skyLL = NULL
      skyChiSq = NULL
    }
    
    if(verbose){message(paste('ProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
    output = list(segim=segim, segim_orig=segim_orig, objects=objects, objects_redo=objects_redo, 
                sky=sky, skyRMS=skyRMS, image=image, mask=mask, segstats=segstats, 
                Nseg=dim(segstats)[1], near=near, group=group, groupstats=groupstats, 
                haralick=haralick, keyvalues=keyvalues, SBlim=SBlim, magzero=magzero, dim=dim(segim), 
                pixscale=pixscale, imarea=imarea, skyLL=skyLL, skyChiSq=skyChiSq, skyChiSqMap=skyChiSqMap,
                gain=gain, call=call, date=date(), time=proc.time()[3]-timestart, ProFound.version=packageVersion('ProFound'), 
                R.version=R.version)
  }else{
    if(is.null(keyvalues)){keyvalues = NULL}
    if(keepim==FALSE){image = NULL; mask = NULL}
    if(is.null(mask)){mask = NULL}
    if(!is.null(badpix)){image[badpix] = NA}
    if(verbose){message('No objects in segmentation map - skipping dilations and CoG')}
    if(verbose){message(paste('ProFound is finished! -',round(proc.time()[3]-timestart,3),'sec'))}
    output = list(segim=NULL, segim_orig=NULL, objects=NULL, objects_redo=NULL, sky=sky, 
                skyRMS=skyRMS, image=image, mask=mask, segstats=NULL, Nseg=0, near=NULL, 
                group=NULL, groupstats=NULL, haralick=NULL, keyvalues=keyvalues, SBlim=NULL,  
                magzero=magzero, dim=dim(segim), pixscale=pixscale, imarea=imarea, skyLL=NULL, skyChiSq=NULL,
                skyChiSqMap=NULL, gain=gain, call=call, date=date(), time=proc.time()[3]-timestart, 
                ProFound.version=packageVersion('ProFound'), R.version=R.version)
  }
  class(output) = 'profound'
  return(invisible(output))
}

plot.profound=function(x, logR50=TRUE, dmag=0.5, hist='sky', ...){
  suppressWarnings({
  if(!inherits(x, 'profound')){
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
  
  segdiff = x$segim - x$segim_orig
  segdiff[segdiff<0]=0
  
  if(!is.null(x$mask)){
    x$image[x$mask>0] = NA
    x$mask[is.na(x$image)] = 1
  }
  
  image = (x$image-x$sky)/x$skyRMS
  image[x$skyRMS < 0] = NA
  
  if(is.null(x$objects)){
    x$objects = matrix(0L,dim(x$segim)[1],dim(x$segim)[2])
    x$objects[] = as.logical(x$segim)
  }
  
  #Masked and unmasked stats:
  
  if(!is.null(x$mask)){
    stat_mean_sky = signif(mean(x$sky[x$mask==0], na.rm=TRUE),4)
    stat_sd_sky = signif(sd(x$sky[x$mask==0], na.rm=TRUE),4)
    stat_mean_skyRMS = signif(mean(x$skyRMS[x$mask==0], na.rm=TRUE),4)
    stat_sd_skyRMS = signif(sd(x$skyRMS[x$mask==0], na.rm=TRUE),4)
    stat_cor_sky_skyRMS = signif(cor(as.numeric(x$sky[x$mask==0]), as.numeric(x$skyRMS[x$mask==0])^2, use="pairwise.complete.obs"),4)
    if(!is.null(x$objects_redo)){
      stat_cor_sky_image = signif(cor(as.numeric(x$sky[x$mask==0 & x$objects_redo>0]), as.numeric(x$image[x$mask==0 & x$objects_redo>0] - as.numeric(x$sky[x$mask==0 & x$objects_redo>0])), use="pairwise.complete.obs"),4)
    }else{
      stat_cor_sky_image = signif(cor(as.numeric(x$sky[x$mask==0 & x$objects>0]), as.numeric(x$image[x$mask==0 & x$objects>0] - as.numeric(x$sky[x$mask==0 & x$objects>0])), use="pairwise.complete.obs"),4)
    }
  }else{
    stat_mean_sky = signif(mean(x$sky, na.rm=TRUE),4)
    stat_sd_sky = signif(sd(x$sky, na.rm=TRUE),4)
    stat_mean_skyRMS = signif(mean(x$skyRMS, na.rm=TRUE),4)
    stat_sd_skyRMS = signif(sd(x$skyRMS, na.rm=TRUE),4)
    stat_cor_sky_skyRMS = signif(cor(as.numeric(x$sky), as.numeric(x$skyRMS)^2, use="pairwise.complete.obs"),4)
    if(!is.null(x$objects_redo)){
      stat_cor_sky_image = signif(cor(as.numeric(x$sky[x$objects_redo>0]), as.numeric(x$image[x$objects_redo>0] - as.numeric(x$sky[x$objects_redo>0])), use="pairwise.complete.obs"),4)
    }else{
      stat_cor_sky_image = signif(cor(as.numeric(x$sky[x$objects>0]), as.numeric(x$image[x$objects>0] - as.numeric(x$sky[x$objects>0])), use="pairwise.complete.obs"),4)
    }
  }
  
  if(requireNamespace("imager", quietly = TRUE)){
    outline = x$objects - as.matrix(imager::erode(imager::as.cimg(x$objects),imager::imfill(3,3,val=1)))
  }else{
    outline = NULL
  }
  
  cmap = hcl.colors(100, 'RdYlBu', rev=TRUE)
  maximg = quantile(abs(image[is.finite(image)]), 0.995, na.rm=TRUE)
  stretchscale = 1/median(abs(image), na.rm=TRUE)
  
  layout(matrix(1:9, 3, byrow=TRUE))
  
  if(!is.null(x$keyvalues)){
  
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 1 image
    if(requireNamespace("Rwcs", quietly = TRUE)){
      Rwcs::Rwcs_image(image=image, keyvalues=x$keyvalues, stretchscale=stretchscale, locut=-maximg, hicut=maximg, range=c(-1,1), type='num', zlim=c(-1,1), col=cmap)
    }else{
      stop("The Rwcs package is need to process the keyvalues. Install from GitHub asgr/Rfits.")
    }
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    if(!is.null(outline)){
      magimage(outline, col=c(NA,'black'), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 2 seg
    profoundSegimPlot(profound=x)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 3 dilate
    if(requireNamespace("Rwcs", quietly = TRUE)){
      Rwcs::Rwcs_image(image=image, keyvalues=x$keyvalues)
    }else{
      stop("The Rwcs package is need to process the keyvalues. Install from GitHub asgr/Rfits.")
    }
    magimage(segdiff, col=c(NA, rainbow(max(x$segim,na.rm=TRUE), end=2/3)), magmap=FALSE, add=TRUE)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 4 mag hist
    if(is.null(x$imarea)){
      imarea=prod(x$dim)*x$pixscale^2/(3600^2)
    }else{
      imarea=x$imarea
    }
    temphist=maghist(x$segstats$mag, log='y', scale=(2*dmag)/x$imarea, breaks=seq(floor(min(x$segstats$mag, na.rm = TRUE)), ceiling(max(x$segstats$mag, na.rm = TRUE)),by=0.5),
                     xlab='mag', ylab=paste('#/deg-sq/d',dmag,'mag',sep=''), grid=TRUE, verbose=FALSE)
    #magplot(temphist, log='y', xlab='mag', ylab=expression('#'/'deg-sq'/'dmag'), grid=TRUE)
    ymax=log10(max(temphist$counts,na.rm = T))
    xmax=temphist$mids[which.max(temphist$counts)]
    abline(ymax - xmax*0.4, 0.4, col='red')
    abline(v=xmax+0.25, col='red')
    axis(side=1, at=xmax+0.25, labels=xmax+0.25, tick=FALSE, line=-1, col.axis='red')
    legend('topleft', legend=paste('N:',length(x$segstats$mag)), bty='n', text.col='black')
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 5 sky
    if(requireNamespace("Rwcs", quietly = TRUE)){
      Rwcs::Rwcs_image(image=x$sky-median(x$sky,na.rm=TRUE), keyvalues=x$keyvalues, qdiff=TRUE)
    }else{
      stop("The Rwcs package is need to process the keyvalues. Install from GitHub asgr/Rfits.")
    }
    if(!is.null(x$mask)){
      magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
    if(!is.null(outline)){
      magimage(outline, col=c(NA,'black'), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
    legend('topleft',legend=c('Sky',paste0('Mean: ',stat_mean_sky),paste0('SD: ',stat_sd_sky)), bty='n', text.col='magenta')
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 6 skyRMS
    if(requireNamespace("Rwcs", quietly = TRUE)){
      Rwcs::Rwcs_image(image=x$skyRMS, keyvalues=x$keyvalues)
    }else{
      stop("The Rwcs package is need to process the keyvalues. Install from GitHub asgr/Rfits.")
    }
    if(!is.null(x$mask)){
      magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
    if(!is.null(outline)){
      magimage(outline, col=c(NA,'orange'), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
    legend('topleft',legend=c('Sky RMS',paste0('Mean: ',stat_mean_skyRMS),paste0('SD: ',stat_sd_skyRMS)), bty='n', text.col='magenta')
    
    #plot 7 sky den
    if(hist=='iters'){
      maghist(x$segstats$iter, breaks=seq(-0.5,max(x$segstats$iter, na.rm=TRUE)+0.5,by=1), majorn=max(x$segstats$iter, na.rm=TRUE)+1, xlab='Number of Dilations', ylab='#', verbose=FALSE)
    }else if(hist=='sky'){
      try({
        if(!is.null(x$objects_redo)){
          tempsky = image[x$objects_redo==0]
        }else{
          tempsky = image[x$objects==0]
        }
        tempsky = tempsky[tempsky > -8 & tempsky < 8 & !is.na(tempsky)]
        if(!is.null(x$skyLL)){stat_LL = signif(x$skyLL,3)}else{stat_LL = NA}
        if(!is.null(x$skyChiSq)){stat_ChiSq = signif(x$skyChiSq,3)}else{stat_ChiSq = NA}
        magplot(function(x){dnorm(x,mean=0, sd=1)}, grid=TRUE, xlim=c(-6,6), xlab='(image - sky) / skyRMS', ylab='PDF', log='y', ylim=c(1e-8,0.5), lty=2, type='l', col='green4')
        lines(density(tempsky[is.finite(tempsky)], bw=0.1, na.rm=TRUE, from=-8, to=8), col='black')
        if(stat_sd_sky > 0){
          lines(density((x$sky[x$objects==1] - stat_mean_sky) / stat_sd_sky, na.rm=TRUE, from=-8, to=8), col='red')
          lines(density((x$sky[x$objects==0] - stat_mean_sky) / stat_sd_sky, na.rm=TRUE, from=-8, to=8), col='blue')
          legend('bottom', legend=c('Sky pixels','Normal Dist','Sky in objects','Sky in sky'), lty=c(1,2,1,1), col=c('black','green4', 'red','blue'))
          legend('topright',legend=c(paste0('Cor RMS: ',stat_cor_sky_skyRMS), paste0('Cor Image: ',stat_cor_sky_image)), bty='n', text.col='black')
        }else{
          legend('bottom', legend=c('Sky pixels','Normal Dist'), lty=c(1,2), col=c('black','green4'))
        }
        legend('topleft',legend=c(paste0('LL: ',stat_LL), paste0('Chi-Sq: ',stat_ChiSq)), bty='n', text.col='black')
      })
    }else{stop('Not a recognised hist type! Must be iters / sky.')}
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 8 size v mag
    if(logR50){
      magplot(x$segstats$mag, x$segstats$R50, pch='.', col=hsv(alpha=0.5), ylim=c(min(x$segstats$R50, 0.3, na.rm = TRUE), max(x$segstats$R50, 1, na.rm = TRUE)), cex=3, xlab='mag', ylab='R50 / asec', grid=TRUE, log='y')
    }else{
      magplot(x$segstats$mag, x$segstats$R50, pch='.', col=hsv(alpha=0.5), ylim=c(0, max(x$segstats$R50, 1, na.rm = TRUE)), cex=3, xlab='mag', ylab='R50 / asec', grid=TRUE)
    }
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 9 SB v rel_err
    fluxrat=x$segstats$flux/x$segstats$flux_err
    magplot(x$segstats$SB_N90, fluxrat, pch='.', col=hsv(alpha=0.5), ylim=c(0.5,max(fluxrat, 1, na.rm=TRUE)), cex=3, xlab='SB90 / mag/asec-sq', ylab='Flux/Flux-Error', grid=TRUE, log='y')
  
  }else{
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 1 image
    magimage(image, stretchscale=stretchscale, locut=-maximg, hicut=maximg, range=c(-1,1), type='num', zlim=c(-1,1), col=cmap)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    if(!is.null(outline)){
      magimage(outline, col=c(NA,'black'), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 2 seg
    # magimage(x$segim, col=c(NA, rainbow(max(x$segim,na.rm=TRUE), end=2/3)), magmap=FALSE)
    # if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    profoundSegimPlot(profound=x)
    #abline(v=c(0,dim(image)[1]))
    #abline(h=c(0,dim(image)[2]))
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 3 dilate
    magimage(image)
    magimage(segdiff, col=c(NA, rainbow(max(x$segim,na.rm=TRUE), end=2/3)), magmap=FALSE, add=TRUE)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}

    par(mar=c(3.5,3.5,0.5,0.5)) #plot 4 mag hist
    temphist=maghist(x$segstats$mag, log='y', scale=(2*dmag), breaks=seq(floor(min(x$segstats$mag, na.rm = TRUE)), ceiling(max(x$segstats$mag, na.rm = TRUE)),by=0.5), xlab='mag', ylab=paste('#/d',dmag,'mag',sep=''), grid=TRUE, verbose=FALSE)
    ymax=log10(max(temphist$counts,na.rm = T))
    xmax=temphist$mids[which.max(temphist$counts)]
    abline(ymax - xmax*0.4, 0.4, col='red')
    abline(v=xmax+0.25, col='red')
    axis(side=1, at=xmax+0.25, labels=xmax+0.25, tick=FALSE, line=-1, col.axis='red')
    legend('topleft', legend=paste('N:',length(x$segstats$mag)), bty='n', text.col='black')
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 5 sky
    magimage(x$sky-median(x$sky,na.rm=TRUE), qdiff=TRUE)
    if(!is.null(x$mask)){
      magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
    if(!is.null(outline)){
      magimage(outline, col=c(NA,'black'), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
    legend('topleft',legend=c('Sky',paste0('Mean: ',stat_mean_sky),paste0('SD: ',stat_sd_sky)), bty='n', text.col='magenta')
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 6 skyRMS
    magimage(x$skyRMS)
    if(!is.null(x$mask)){
      magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
    if(!is.null(outline)){
      magimage(outline, col=c(NA,'orange'), add=TRUE, magmap=FALSE, zlim=c(0,1))
    }
    legend('topleft',legend=c('Sky RMS',paste0('Mean: ',stat_mean_skyRMS),paste0('SD: ',stat_sd_skyRMS)), bty='n', text.col='magenta')
    
    #plot 7 sky den
    if(hist=='iters'){
      maghist(x$segstats$iter, breaks=seq(-0.5,max(x$segstats$iter, na.rm=TRUE)+0.5,by=1), majorn=max(x$segstats$iter, na.rm=TRUE)+1, xlab='Number of Dilations', ylab='#', verbose=FALSE)
    }else if(hist=='sky'){
      try({
        if(!is.null(x$objects_redo)){
          tempsky = image[x$objects_redo==0]
        }else{
          tempsky = image[x$objects==0]
        }
        tempsky = tempsky[tempsky > -8 & tempsky < 8 & !is.na(tempsky)]
        if(!is.null(x$skyLL)){stat_LL = signif(x$skyLL,3)}else{stat_LL = NA}
        if(!is.null(x$skyChiSq)){stat_ChiSq = signif(x$skyChiSq,3)}else{stat_ChiSq = NA}
        magplot(function(x){dnorm(x,mean=0, sd=1)}, grid=TRUE, xlim=c(-6,6), xlab='(image - sky) / skyRMS', ylab='PDF', log='y', ylim=c(1e-8,0.5), lty=2, type='l', col='green4')
        lines(density(tempsky[is.finite(tempsky)], bw=0.1, na.rm=TRUE, from=-8, to=8), col='black')
        if(stat_sd_sky > 0){
          lines(density((x$sky[x$objects==1] - stat_mean_sky) / stat_sd_sky, na.rm=TRUE, from=-8, to=8), col='red')
          lines(density((x$sky[x$objects==0] - stat_mean_sky) / stat_sd_sky, na.rm=TRUE, from=-8, to=8), col='blue')
          legend('bottom', legend=c('Sky pixels','Normal Dist','Sky in objects','Sky in sky'), lty=c(1,2,1,1), col=c('black','green4', 'red','blue'))
          legend('topright',legend=c(paste0('Cor RMS: ',stat_cor_sky_skyRMS), paste0('Cor Image: ',stat_cor_sky_image)), bty='n', text.col='black')
        }else{
          legend('bottom', legend=c('Sky pixels','Normal Dist'), lty=c(1,2), col=c('black','green4'))
        }
        legend('topleft',legend=c(paste0('LL: ',stat_LL), paste0('Chi-Sq: ',stat_ChiSq)), bty='n', text.col='black')
      })
    }else{stop('Not a recognised hist type! Must be iters / sky.')}
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 8 size v mag
    if(logR50){
      magplot(x$segstats$mag, x$segstats$R50, pch='.', col=hsv(alpha=0.5), ylim=c(min(x$segstats$R50, 0.3, na.rm = TRUE), max(x$segstats$R50, 1, na.rm = TRUE)), cex=3, xlab='mag', ylab='R50 / asec', grid=TRUE, log='y')
    }else{
      magplot(x$segstats$mag, x$segstats$R50, pch='.', col=hsv(alpha=0.5), ylim=c(0, max(x$segstats$R50, 1, na.rm = TRUE)), cex=3, xlab='mag', ylab='R50 / asec', grid=TRUE)
    }
    
    par(mar=c(3.5,3.5,0.5,0.5)) #plot 9 SB v rel_err
    fluxrat=x$segstats$flux/x$segstats$flux_err
    magplot(x$segstats$SB_N90, fluxrat, pch='.', col=hsv(alpha=0.5), ylim=c(0.5,max(fluxrat, 1, na.rm=TRUE)), cex=3, xlab='SB90 / mag/pix-sq', ylab='Flux/Flux-Error', grid=TRUE, log='y')
  }
  })
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
