profoundFitMagPSF=function(xcen=NULL, ycen=NULL, RAcen=NULL, Deccen=NULL, mag=NULL, image=NULL, im_sigma=NULL, mask=NULL, psf=NULL, fit_iters=5, magdiff=1, modxy=FALSE, sigthresh=0, itersub=TRUE, magzero=0, modelout=TRUE, fluxtype='Raw', header=NULL, doProFound=FALSE, findextra=FALSE, verbose=FALSE, ...){
  
  timestart=proc.time()[3]
  
  call=match.call()
  
  if(!requireNamespace("ProFit", quietly = TRUE)){
    hasProFit=TRUE
    message('ProFit package will run things faster! Please install it from CRAN.', call. = FALSE)
  }else{
    hasProFit=FALSE
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
  
  if(is.null(psf)){stop('Missing psf - this is a required input!')}
  
  fluxtype=tolower(fluxtype)
  
  if(fluxtype=='raw' | fluxtype=='adu' | fluxtype=='adus'){
    fluxscale=1
  }else if (fluxtype=='jansky'){
    fluxscale=10^(-0.4*(magzero-8.9))
  }else{
    stop('fluxtype must be Jansky / Raw!')
  }
  
  if((is.null(xcen) & !is.null(ycen))){stop('Need xcen/ycen pair!')}
  if((!is.null(xcen) & is.null(ycen))){stop('Need xcen/ycen pair!')}
  if((is.null(xcen) & is.null(ycen)) & (is.null(RAcen) | is.null(Deccen))){stop('Need RAcen/Decen pair!')}
  
  if((!is.null(xcen) & !is.null(ycen)) & (is.null(RAcen) & is.null(Deccen)) & !is.null(header)){
    RAdeccoords=magWCSxy2radec(x = xcen, y = ycen, header=header)
    RAcen=RAdeccoords[,'RA']
    Deccen=RAdeccoords[,'Dec']
  }
  
  if((is.null(xcen) & is.null(ycen)) & (!is.null(RAcen) & !is.null(Deccen))){
    if(is.null(header)){stop('Need header if using RAcen and Deccec input')}
    xycoords=magWCSradec2xy(RA = RAcen, Dec = Deccen, header=header)
    xcen=xycoords[,'x']
    ycen=xycoords[,'y']
  }
  
  if(doProFound){
    if(verbose){message('Running initial ProFound')}
    protemp=profoundProFound(image, mask=mask, magzero=magzero, header=header, verbose=FALSE, ...)
    if(verbose){message('- subtracting ProFound sky model from input image')}
    image=image-protemp$sky
    if(is.null(xcen) & is.null(ycen) & is.null(mag)){
      if(verbose){message('- using ProFound xcen / ycen / mag for source properties')}
      xcen=protemp$segstats$xcen
      ycen=protemp$segstats$ycen
      mag=protemp$segstats$mag
    }
    if(is.null(im_sigma)){
      if(verbose){message('- using ProFound skyRMS model for im_sigma')}
      im_sigma=protemp$skyRMS
    }
  }else{
    protemp=NULL
  }
  
  if(is.null(mag)){
    tempflux=image[cbind(ceiling(xcen),ceiling(ycen))]
    tempflux=tempflux*sum(image, na.rm=TRUE)/sum(tempflux, na.rm=TRUE)
    mag=profoundFlux2Mag(flux=tempflux, magzero=magzero)
  }
  if(is.null(im_sigma)){stop('Need im_sigma!')}
  
  if(is.null(RAcen)){RAcen=rep(NA,length(mag))}
  if(is.null(Deccen)){Deccen=rep(NA,length(mag))}
  
  mag[is.na(mag)]=Inf
  
  if(!is.null(mask)){
    image[mask!=0]=NA #means we will ignore the masked bits when doing the LL
  }
  
  image_orig=image
  
  Nmodels=length(mag)
  
  if(verbose){
    Nprint=signif(Nmodels/5,1)
  }else{
    Nprint=Inf
  }
  
  if(Nmodels != length(xcen)){
    stop('xcen must be the same length as mag!')
  }
  if(Nmodels != length(ycen)){
    stop('ycen must be the same length as mag!')
  }
  
  #sort vectors by flux
  
  fluxorder=order(mag)
  
  xcen=xcen[fluxorder]
  ycen=ycen[fluxorder]
  RAcen=RAcen[fluxorder]
  Deccen=Deccen[fluxorder]
  mag=mag[fluxorder]
  
  if(modxy){
    xygrid=expand.grid(1:dim(psf)[1],1:dim(psf)[2])-0.5
  }
  
  psf=psf/sum(psf,na.rm=TRUE)
  beamscale=sqrt(1/sum(psf^2,na.rm=TRUE))
  
  diffmag=rep(0,Nmodels)
  psfLL=rep(0,Nmodels)
  flux=rep(0,Nmodels)
  flux_err=rep(0,Nmodels)
  beam_err=rep(0,Nmodels)
  
  if(verbose){message('Iterating over source model')}
  
  for(j in 1:fit_iters){
    if(verbose){message(paste('- iteration',j,'of',fit_iters))}
    
    if(itersub){
      if(j==1){
        dofull=TRUE
      }else{
        dofull=FALSE
      }
    }else{
      dofull=TRUE
    }
    
    if(dofull){
      if(hasProFit){
        modellist = list(
          pointsource = list(
            xcen = xcen[is.finite(mag)],
            ycen = ycen[is.finite(mag)],
            mag = mag[is.finite(mag)]
          )
        )
        fullmodel=ProFit::profitMakeModel(modellist=modellist, dim=dim(image), psf=psf, magzero=magzero)$z
      }else{
        fullmodel=.genPointSource(
          xcen = xcen[is.finite(mag)],
          ycen = ycen[is.finite(mag)],
          flux=profoundMag2Flux(mag=mag[is.finite(mag)], magzero=magzero),
          psf=psf,
          dim=dim(image)
        )
      }
      
      image=image_orig-fullmodel
    }
    
    if(j==1){
      if(modelout){
        origmodel=fullmodel
        origLL= -0.5*sum((image_orig-fullmodel)^2/(im_sigma^2))
      }else{
        origmodel=NULL
        origLL=NULL
      }
    }
    
    for(i in 1:Nmodels){
      if(verbose & (i %% Nprint == 0)  & Nmodels>=1e3){
        message(paste(' - model',i,'of',Nmodels))
      }
      if(is.finite(mag[i])){
        image_cut=magcutout(image, loc=c(xcen[i],ycen[i]), box=dim(psf))
        sigma_cut=magcutout(im_sigma, loc=c(xcen[i],ycen[i]), box=dim(psf))$image
        
        if(hasProFit){
          singlist = list(
            pointsource = list(
              xcen = image_cut$loc[1],
              ycen = image_cut$loc[2],
              mag = mag[i]
            )
          )
          singmodel=ProFit::profitMakeModel(modellist=singlist, dim=dim(psf), psf=psf, magzero=magzero)$z
        }else{
          singmodel=.genPointSource(
            xcen = image_cut$loc[1],
            ycen = image_cut$loc[2],
            flux=profoundMag2Flux(mag=mag[i], magzero=magzero),
            psf=psf,
            dim=dim(psf)
          )
        }
        
        if(modxy){
          if(anyNA(image_cut$image)){
            selNA=which(is.na(image_cut$image))
            image_cut$image[selNA]=image_cut$image[length(image_cut$image)+1-selNA]
            sigma_cut[selNA]=sigma_cut[length(image_cut$image)+1-selNA]
          }
          
          select=which(image_cut$image+singmodel > sigma_cut*sigthresh)
          
          if(length(select)>0){
            weights=(singmodel*(image_cut$image+singmodel))[select]
            newx=.meanwt(xygrid[select,1], weights)
            newy=.meanwt(xygrid[select,2], weights)
            if(newx-image_cut$loc[1] > 0.5){newx=image_cut$loc[1]+0.5}
            if(newx-image_cut$loc[1] < -0.5){newx=image_cut$loc[1]-0.5}
            if(newy-image_cut$loc[2] > 0.5){newy=image_cut$loc[2]+0.5}
            if(newy-image_cut$loc[2] < -0.5){newy=image_cut$loc[2]-0.5}
            
            image_cut$loc=c(newx,newy)
            xcen[i]=image_cut$loc[1]+image_cut$loc.diff[1]
            ycen[i]=image_cut$loc[2]+image_cut$loc.diff[2]
            
            if(hasProFit){
              singlist = list(
                pointsource = list(
                  xcen = image_cut$loc[1],
                  ycen = image_cut$loc[2],
                  mag = mag[i]
                )
              )
              singmodel=ProFit::profitMakeModel(modellist=singlist, dim=dim(psf), psf=psf, magzero=magzero)$z
            }else{
              singmodel=.genPointSource(
                xcen = image_cut$loc[1],
                ycen = image_cut$loc[2],
                flux=profoundMag2Flux(mag=mag[i], magzero=magzero),
                psf=psf,
                dim=dim(psf)
              )
            }
          }
        }
        
        if(j<fit_iters){
          diffmag[i] = optim(par=0, fn=.minlike_mag, method='Brent', singmodel = singmodel, image=image_cut$image, im_sigma = sigma_cut, lower=-magdiff, upper=magdiff)$par
        }else{
          finaloptim=optim(par=0, fn=.minlike_mag, method='Brent', singmodel = singmodel, image=image_cut$image, im_sigma = sigma_cut, lower=-magdiff, upper=magdiff)
          diffmag[i]=finaloptim$par
          psfLL[i]=-finaloptim$value
          flux[i]=profoundMag2Flux(mag=mag[i]+diffmag[i], magzero=magzero)
          fluxhess=optimHess(par=flux[i]*(10^(-0.4*diffmag[i])-1), fn=.minlike_flux, singmodel = singmodel, image=image_cut$image, im_sigma = sigma_cut)
          flux_err[i]=sqrt(1/abs(as.numeric(fluxhess)))
          beam_err[i]=mean(sigma_cut,na.rm =TRUE)*beamscale
        }
        if(itersub){
          rescale=10^(-0.4*diffmag[i])-1
          #Now this is a weird thing. Much more efficient memory wise to use this Rcpp function to allocate matrix subsets. Not clear why since it should not be copying on modification in this case...?
          image=.addmat_cpp(image, -singmodel[image_cut$xsel-image_cut$loc.diff[1], image_cut$ysel-image_cut$loc.diff[2]]*rescale, range(image_cut$xsel), range(image_cut$ysel))
          #image[image_cut$xsel,image_cut$ysel]=image[image_cut$xsel,image_cut$ysel]-(singmodel[image_cut$xsel-image_cut$loc.diff[1], image_cut$ysel-image_cut$loc.diff[2]]*rescale)
        }
      }else{
        if(j==fit_iters){
          psfLL[i]=NA
          flux[i]=0
          sigma_cut=magcutout(im_sigma, loc=c(xcen[i],ycen[i]), box=dim(psf))$image
          flux_err[i]=mean(sigma_cut,na.rm =TRUE)*beamscale
        }
      }
    }
    mag=mag+diffmag
  }
  
  
  if(modelout){
    if(hasProFit){
      modellist = list(
        pointsource = list(
          xcen = xcen[is.finite(mag)],
          ycen = ycen[is.finite(mag)],
          mag = mag[is.finite(mag)]
        )
      )
      fullmodel=ProFit::profitMakeModel(modellist=modellist, dim=dim(image), psf=psf, magzero=magzero)$z
    }else{
      fullmodel=.genPointSource(
        xcen = xcen[is.finite(mag)],
        ycen = ycen[is.finite(mag)],
        flux=profoundMag2Flux(mag=mag[is.finite(mag)], magzero=magzero),
        psf=psf,
        dim=dim(image)
      )
    }
    
    finalLL= -0.5*sum((image_orig-fullmodel)^2/(im_sigma^2))
  }else{
    fullmodel=NULL
    finalLL=NULL
  }
  flux_err=sqrt(flux_err^2 + (diffmag*flux/(2.5/log(10)))^2)
  flux_err[!is.finite(flux_err)]=beam_err[!is.finite(flux_err)]
  flux_err[flux_err<beam_err]=beam_err[flux_err<beam_err]
  mag_err=(2.5/log(10))*flux_err/flux
  
  pchi=pchisq(prod(dim(psf))*(flux/flux_err)^2,prod(dim(psf))-1,log.p=TRUE)
  signif=qnorm(pchi, log.p=TRUE)
  
  if(findextra){
    if(verbose){message('Finding extra sources with ProFound')}
    protemp=profoundProFound(image_orig-fullmodel+protemp$sky, mask=mask, magzero=magzero, header=header, verbose=FALSE, ...)
    image_orig=image_orig-protemp$sky
    if(!is.null(protemp$segstats)){
      xcen=c(xcen, protemp$segstats$xcen)
      ycen=c(ycen, protemp$segstats$ycen)
      mag=c(mag, protemp$segstats$mag)
      
      if(verbose){message(paste('- gained',protemp$Nseg,'additional sources'))}
      
      if(verbose){message('Rerunning source model with extra sources')}
      
      rerun=profoundFitMagPSF(xcen=xcen, ycen=ycen, mag=mag, image=image_orig, im_sigma=im_sigma, mask=mask, psf=psf, fit_iters=fit_iters, magdiff=magdiff, modxy=modxy, sigthresh=sigthresh, itersub=itersub, magzero=magzero, modelout=modelout, fluxtype='Raw', header=header, doProFound=FALSE, findextra=FALSE, verbose=verbose)
      
      seltarget=1:Nmodels
      xcen=rerun$psfstats$xcen[seltarget]
      ycen=rerun$psfstats$ycen[seltarget]
      RAcen=rerun$psfstats$RAcen[seltarget]
      Deccen=rerun$psfstats$Deccen[seltarget]
      flux=rerun$psfstats$flux[seltarget]
      flux_err=rerun$psfstats$flux_err[seltarget]
      mag=rerun$psfstats$mag[seltarget]
      mag_err=rerun$psfstats$mag_err[seltarget]
      psfLL=rerun$psfstats$psfLL[seltarget]
      signif=rerun$psfstats$signif[seltarget]
      
      selextra=(Nmodels+1):(dim(rerun$psfstats)[1])
      xcen_extra=rerun$psfstats$xcen[selextra]
      ycen_extra=rerun$psfstats$ycen[selextra]
      RAcen_extra=rerun$psfstats$RAcen[selextra]
      Deccen_extra=rerun$psfstats$Deccen[selextra]
      flux_extra=rerun$psfstats$flux[selextra]
      flux_err_extra=rerun$psfstats$flux_err[selextra]
      mag_extra=rerun$psfstats$mag[selextra]
      mag_err_extra=rerun$psfstats$mag_err[selextra]
      psfLL_extra=rerun$psfstats$psfLL[selextra]
      signif_extra=rerun$psfstats$signif[selextra]
      psfstats_extra=data.frame(xcen=xcen_extra, ycen=ycen_extra, RAcen=RAcen_extra, Deccen=Deccen_extra, flux=flux_extra, flux_err=flux_err_extra, mag=mag_extra, mag_err=mag_err_extra, psfLL=psfLL_extra, signif=signif_extra)
      
      fullmodel=rerun$finalmodel
    }else{
      psfstats_extra=NULL
    }
  }else{
    psfstats_extra=NULL
  }
  
  flux=flux*fluxscale
  flux_err=flux_err*fluxscale
  
  if(findextra & !is.null(psfstats_extra)){
    psfstats_extra$flux=psfstats_extra$flux*fluxscale
    psfstats_extra$flux_err=psfstats_extra$flux_err*fluxscale
  }
  
  psfstats=data.frame(xcen=xcen, ycen=ycen, RAcen=RAcen, Deccen=Deccen, flux=flux, flux_err=flux_err, mag=mag, mag_err=mag_err, psfLL=psfLL, signif=signif)
  psfstats=psfstats[match(1:Nmodels,fluxorder),] # Reorder back to the input
  
  return(invisible(list(psfstats=psfstats, origLL=origLL, finalLL=finalLL, origmodel=origmodel, finalmodel=fullmodel, image=image_orig, header=header, psfstats_extra=psfstats_extra, profound=protemp, call=call, date=date(), time=proc.time()[3]-timestart, ProFound.version=packageVersion('ProFound'), R.version=R.version)))
}

.minlike_mag=function(par,singmodel,image,im_sigma){
  0.5*sum((image-(singmodel*(10^(-0.4*par)-1)))^2/im_sigma^2, na.rm=TRUE)
}

.minlike_flux=function(par,singmodel,image,im_sigma){
  flux=sum(singmodel, na.rm=TRUE)
  0.5*sum((image-(singmodel*(1+par/flux)))^2/im_sigma^2, na.rm=TRUE)
}
