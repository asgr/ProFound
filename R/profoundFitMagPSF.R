profoundFitMagPSF=function(xcen=NULL, ycen=NULL, RAcen=NULL, Deccen=NULL, mag=NULL, 
                           image=NULL, im_sigma=NULL, mask=NULL, psf=NULL, fit_iters=5, 
                           magdiff=1, modxy=FALSE, sigthresh=0, itersub=TRUE, magzero=0, 
                           modelout=TRUE, fluxtype='Raw', psf_redosky=FALSE, fluxext=FALSE, 
                           header=NULL, doProFound=FALSE, findextra=FALSE, verbose=FALSE, ...){
  
  timestart=proc.time()[3]
  
  call=match.call()
  
  if(!requireNamespace("ProFit", quietly = TRUE)){
    hasProFit=FALSE
    message('ProFit package will run things faster! Please install it from CRAN.', call. = FALSE)
  }else{
    hasProFit=TRUE
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
  
  if(is.null(psf)){stop('Missing psf - this is a required input!')}
  
  fluxtype=tolower(fluxtype)
  
  if(fluxtype=='raw' | fluxtype=='adu' | fluxtype=='adus'){
    fluxscale=1
  }else if (fluxtype=='jansky'){
    fluxscale=10^(-0.4*(magzero-8.9))
  }else if (fluxtype=='microjansky'){
    fluxscale=10^(-0.4*(magzero-23.9))
  }else{
    stop('fluxtype must be Jansky / Microjansky / Raw!')
  }
  
  if((is.null(xcen) & !is.null(ycen))){stop('Need xcen/ycen pair!')}
  if((!is.null(xcen) & is.null(ycen))){stop('Need xcen/ycen pair!')}
  if((is.null(xcen) & is.null(ycen)) & (is.null(RAcen) | is.null(Deccen))){stop('Need RAcen/Decen pair!')}
  
  if((!is.null(xcen) & !is.null(ycen)) & (is.null(RAcen) & is.null(Deccen)) & !is.null(header)){
    #if(requireNamespace("Rwcs", quietly = TRUE)){
    #  RAdeccoords=Rwcs::Rwcs_p2s(x = xcen, y = ycen, pixcen = 'R', header=header)
    #}else{
      RAdeccoords=magWCSxy2radec(x = xcen, y = ycen, header=header)
    #}
    RAcen=RAdeccoords[,'RA']
    Deccen=RAdeccoords[,'Dec']
  }
  
  if((is.null(xcen) & is.null(ycen)) & (!is.null(RAcen) & !is.null(Deccen))){
    if(is.null(header)){stop('Need header if using RAcen and Deccec input')}
    #if(requireNamespace("Rwcs", quietly = TRUE)){
    #  xycoords=Rwcs::Rwcs_s2p(RA = RAcen, Dec = Deccen, pixcen = 'R', header=header)
    #}else{
      xycoords=magWCSradec2xy(RA = RAcen, Dec = Deccen, header=header)
    #}
    xcen=xycoords[,'x']
    ycen=xycoords[,'y']
  }
  
  if(doProFound){
    if(verbose){message('Running initial ProFound')}
    if(psf_redosky){
      protemp=profoundProFound(image, mask=mask, magzero=magzero, header=header, verbose=FALSE, ...)
      if(verbose){message('- subtracting ProFound sky model from input image')}
      image=image-protemp$sky
    }else{
      protemp=profoundProFound(image, mask=mask, sky=0, magzero=magzero, header=header, verbose=FALSE, ...)
    }
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
  
  if(!is.null(mask)){
    image[mask!=0]=NA #means we will ignore the masked bits when doing the LL
  }
  
  goodlocs=ceiling(xcen)>=1 & ceiling(xcen)<= dim(image)[1] & ceiling(ycen)>=1 & ceiling(ycen)<= dim(image)[2]
  
  if(is.null(mag)){
    tempflux=rep(NA,length(xcen))
    tempflux[goodlocs]=image[cbind(ceiling(xcen[goodlocs]),ceiling(ycen[goodlocs]))] #only if source is in image
    tempflux=tempflux*sum(magclip(image)$x, na.rm=TRUE)/sum(tempflux, na.rm=TRUE) # clipping to be safe from bad pixels!
    mag=profoundFlux2Mag(flux=tempflux, magzero=magzero)
  }
  if(is.null(im_sigma)){stop('Need im_sigma!')}
  
  if(is.null(RAcen)){RAcen=rep(NA,length(mag))}
  if(is.null(Deccen)){Deccen=rep(NA,length(mag))}
  
  mag[is.na(mag)]=Inf
  
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
  goodlocs=goodlocs[fluxorder]
  
  if(modxy){
    xygrid=expand.grid(1:dim(psf)[1],1:dim(psf)[2])-0.5
  }
  
  psf=psf/sum(psf,na.rm=TRUE)
  beamscale=sqrt(1/sum(psf^2,na.rm=TRUE))
  
  diffmag=rep(NA,Nmodels)
  psfLL=rep(NA,Nmodels)
  flux=rep(NA,Nmodels)
  flux_err=rep(NA,Nmodels)
  beam_err=rep(mean(im_sigma, na.rm=TRUE)*beamscale,Nmodels)
  
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
        origLL= -0.5*sum((image_orig-fullmodel)^2/(im_sigma^2), na.rm = TRUE)
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
          flux[i]=NA
          if(goodlocs[i]){
            sigma_cut=magcutout(im_sigma, loc=c(xcen[i],ycen[i]), box=dim(psf))$image
            flux_err[i]=mean(sigma_cut, na.rm=TRUE)*beamscale
          }else{
            flux_err[i]=NA
          }
        }
      }
    }
    mag=mag+diffmag
  }
  
  if(modelout | findextra | fluxext){
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
    
    finalLL= -0.5*sum((image_orig-fullmodel)^2/(im_sigma^2), na.rm = TRUE)
  }else{
    fullmodel=NULL
    finalLL=NULL
  }
  
  if(fluxext){
    for(i in 1:Nmodels){
      if(verbose & (i %% Nprint == 0)  & Nmodels>=1e3){
        message(paste(' - model',i,'of',Nmodels))
      }
      if(is.finite(mag[i])){
        image_cut = magcutout(image_orig, loc=c(xcen[i],ycen[i]), box=dim(psf))
        model_cut = magcutout(fullmodel, loc=c(xcen[i],ycen[i]), box=dim(psf))

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
        
        sing_rel = singmodel / model_cut$image
        sing_rel[is.na(sing_rel)] = 0
        sing_rel[is.nan(sing_rel)] = 0
        sing_rel[sing_rel<0] = 0
        sing_rel[sing_rel>1] = 1
        
        flux[i] = sum(image_cut$image*sing_rel, na.rm=TRUE)
        mag[i] = profoundFlux2Mag(flux[i], magzero=magzero)
      }
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
      
      finalLL= -0.5*sum((image_orig-fullmodel)^2/(im_sigma^2), na.rm = TRUE)
    }else{
      fullmodel=NULL
      finalLL=NULL
    }
  }
  
  flux_err=sqrt(flux_err^2 + (diffmag*flux/(2.5/log(10)))^2)
  flux_err[!is.finite(flux_err)]=beam_err[!is.finite(flux_err)]
  flux_err[which(flux_err<beam_err)]=beam_err[which(flux_err<beam_err)]
  mag_err=(2.5/log(10))*flux_err/flux
  
  pchi=pchisq(prod(dim(psf))*(flux/flux_err)^2,prod(dim(psf))-1,log.p=TRUE)
  signif=qnorm(pchi, log.p=TRUE)
  
  if(findextra){
    if(verbose){message('Finding extra sources with ProFound')}
    if(psf_redosky){
      image_orig=image_orig+protemp$sky
      protemp=profoundProFound(image_orig-fullmodel, mask=mask, magzero=magzero, header=header, verbose=FALSE, ...)
      image_orig=image_orig-protemp$sky
    }else{
      protemp=profoundProFound(image_orig-fullmodel, mask=mask, sky=0, magzero=magzero, header=header, verbose=FALSE, ...)
    }
    if(!is.null(protemp$segstats)){
      addxcen = protemp$segstats$xcen
      addycen = protemp$segstats$ycen
      addmag = protemp$segstats$mag
      goodadd = which(is.finite(addxcen) & is.finite(addycen) & is.finite(addmag))
      if(length(goodadd)>0){
        xcen=c(xcen, addxcen[goodadd])
        ycen=c(ycen, addycen[goodadd])
        mag=c(mag, addmag[goodadd])
      
        if(verbose){message(paste('- gained',protemp$Nseg,'additional sources'))}
      
        if(verbose){message('Rerunning source model with extra sources')}
      
        rerun=profoundFitMagPSF(xcen=xcen, ycen=ycen, mag=mag, image=image_orig,
                                im_sigma=im_sigma, mask=mask, psf=psf, fit_iters=fit_iters,
                                magdiff=magdiff, modxy=modxy, sigthresh=sigthresh,
                                itersub=itersub, magzero=magzero, modelout=modelout,
                                fluxtype='Raw', header=header, doProFound=FALSE,
                                findextra=FALSE, fluxext=fluxext, verbose=verbose)
      
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
        beam_err=rerun$psfstats$beam_err[seltarget]
        
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
        beam_err_extra=rerun$psfstats$beam_err[selextra]
        psfstats_extra=data.frame(xcen=xcen_extra, ycen=ycen_extra, RAcen=RAcen_extra, 
                                  Deccen=Deccen_extra, flux=flux_extra, flux_err=flux_err_extra, 
                                  mag=mag_extra, mag_err=mag_err_extra, psfLL=psfLL_extra, 
                                  signif=signif_extra, beam_err=beam_err_extra)
        
        fullmodel=rerun$finalmodel
      }else{
        psfstats_extra=NULL
      }
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
    row.names(psfstats_extra)=NULL
  }
  
  psfstats=data.frame(xcen=xcen, ycen=ycen, RAcen=RAcen, Deccen=Deccen, flux=flux,
                      flux_err=flux_err, mag=mag, mag_err=mag_err, psfLL=psfLL, 
                      signif=signif, beam_err=beam_err)
  psfstats=psfstats[match(1:Nmodels,fluxorder),] # Reorder back to the input
  row.names(psfstats)=NULL
  
  output=list(psfstats=psfstats, origLL=origLL, finalLL=finalLL, origmodel=origmodel, 
              finalmodel=fullmodel, image=image_orig, header=header, psfstats_extra=psfstats_extra, 
              profound=protemp, mask=mask, call=call, date=date(), time=proc.time()[3]-timestart, 
              ProFound.version=packageVersion('ProFound'), R.version=R.version)
  class(output)='fitmagpsf'
  return(invisible(output))
}

plot.fitmagpsf=function(x, ...){
  
  if(!is.null(x$mask)){
    x$image[x$mask!=0]=NA #means we will ignore the masked bits when doing the LL
    x$finalmodel[x$mask!=0]=NA #means we will ignore the masked bits when doing the LL
  }
  
  layout(rbind(1:3))
  if(is.null(x$header)){
    magimage(x$image, qdiff=TRUE, ...)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    legend('topleft', legend='image', bty='n', cex=3, pch='')
    
    magimage(x$finalmodel, qdiff=TRUE, ...)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    legend('topleft', legend='model', bty='n', cex=3, pch='')
    
    magimage(x$image - x$finalmodel, qdiff=TRUE, ...)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    legend('topleft', legend='image - model', bty='n', cex=3, pch='')
  }else{
    magimageWCS(x$image, header=x$header, qdiff=TRUE, ...)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    legend('topleft', legend='image', bty='n', cex=3, pch='')
    
    magimageWCS(x$finalmodel, header=x$header, qdiff=TRUE, ...)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    legend('topleft', legend='model', bty='n', cex=3, pch='')
    
    magimageWCS(x$image - x$finalmodel, header=x$header, qdiff=TRUE, ...)
    if(!is.null(x$mask)){magimage(x$mask!=0, col=c(NA,hsv(alpha=0.2)), add=TRUE, magmap=FALSE, zlim=c(0,1))}
    legend('topleft', legend='image - model', bty='n', cex=3, pch='')
  }
}

.minlike_mag=function(par,singmodel,image,im_sigma){
  0.5*sum((image-(singmodel*(10^(-0.4*par)-1)))^2/im_sigma^2, na.rm=TRUE)
}

.minlike_flux=function(par,singmodel,image,im_sigma){
  flux=sum(singmodel, na.rm=TRUE)
  0.5*sum((image-(singmodel*(1+par/flux)))^2/im_sigma^2, na.rm=TRUE)
}
