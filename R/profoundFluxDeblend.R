profoundFluxDeblend=function(image=NULL, segim=NULL, segstats=NULL, groupim=NULL, groupsegID=NULL, sky=0, profound=NULL, magzero=0, df=3, radtrunc=2, iterative=FALSE, doallstats=TRUE, lowmemory=FALSE, deblendtype='fit', psf=NULL, fluxweight='sum', convtype='brute', convmode='extended', Ndeblendlim=Inf){
  
  if(!is.null(image)){
    if(class(image)=='profound'){
      if(is.null(segim)){segim=image$segim}
      if(is.null(segstats)){segstats=image$segstats}
      if(is.null(groupim)){groupim=image$group$groupim}
      if(is.null(groupsegID)){groupsegID=image$group$groupsegID}
      if(missing(sky)){sky=image$sky}
      if(missing(magzero)){magzero=image$magzero}
      image=image$image
      if(is.null(image)){stop('Need image in profound object to be non-Null')}
    }
  }
  
  if(!is.null(profound)){
    if(class(profound) != 'profound'){
      stop('Class of profound input must be of type \'profound\'')
    }
    if(is.null(image)){image=profound$image}
    if(is.null(segim)){segim=profound$segim}
    if(is.null(segstats)){segstats=profound$segstats}
    if(is.null(groupim)){groupim=profound$group$groupim}
    if(is.null(groupsegID)){groupsegID=profound$group$groupsegID}
    if(missing(magzero)){magzero=profound$magzero}
  }
  
  if(is.null(groupim) | is.null(groupsegID)){
    if(is.null(segim)){
      stop('Need segim if missing groupim or groupsegID!')
    }
    group=profoundSegimGroup(segim)
    groupim=group$groupim
    groupsegID=group$groupsegID
  }
  
  image=image-sky
  
  groupsegID=groupsegID[groupsegID$Ngroup>1 & as.numeric(groupsegID$Npix)*groupsegID$Ngroup<Ndeblendlim,,drop=FALSE]
  
  if(dim(groupsegID)[1]==0){
    return(NULL)
  }
  
  if(lowmemory){
    rm(sky)
  }else{
    groupID=Var1=Var2=NULL
    groupref=data.table(groupID=as.integer(groupim), expand.grid(1:dim(groupim)[1],1:dim(groupim)[2]), keyby='groupID')
    setkey(groupref, groupID)
    groupref[groupID %in% groupsegID$groupID]
  }
  
  output=data.frame(groupID=rep(groupsegID$groupID,groupsegID$Ngroup), segID=unlist(groupsegID$segID), flux_db=NA, mag_db=NA, N100_db=NA)
  
  if(deblendtype=='fit' & iterative){
    output[,"flux_db"]=segstats[match(output$segID, segstats$segID),"flux"]
    output=output[order(output[,"groupID"],-output[,"flux_db"]),,drop=FALSE]
  }
  Npad=1
  image_temp=image
  
  if(deblendtype=='psf'){
    if(!requireNamespace("ProFit", quietly = TRUE)){
     stop('The ProFit package is needed for the deblendtype=psf option to work. Please install it from CRAN.', call. = FALSE)
    }
    if(convmode=='extended'){
      convolver_brute=ProFit::profitMakeConvolver(type=convtype, dim(image_temp), psf)
    }
    segim_temp=matrix(FALSE,dim(image_temp)[1],dim(image_temp)[2])
    groupim_temp=matrix(FALSE,dim(image_temp)[1],dim(image_temp)[2])
    beamsum=sum(psf)
  }
  
  for(i in 1:dim(groupsegID)[1]){
    #print(paste(i,'in',dim(groupsegID)[1]))
    Ngroup=groupsegID[i,"Ngroup"]
    segIDlist=output[output[,'groupID']==groupsegID[i,"groupID"],"segID"]
    segIDlist=segIDlist[segIDlist>0]
  
    if(lowmemory){
      tempgridgroup=which(groupim==groupsegID[i,"groupID"], arr.ind=TRUE)
    }else{
      tempgridgroup=groupref[groupID==groupsegID[i,"groupID"],matrix(c(Var1,Var2),ncol=2)]
    }
    
    weightmatrix=matrix(0,length(tempgridgroup[,1]),length(segIDlist))
    Qseg_db=rep(0,length(segIDlist))
    fluxfrac=rep(0,length(segIDlist))
    beamcorrect=rep(1,length(segIDlist))
    
    groupsum=sum(image[tempgridgroup],na.rm=TRUE)
    
    for(j in 1:length(segIDlist)){
      
      tempgridseg=which(segim[tempgridgroup]==segIDlist[j])
      segsum=sum(image_temp[tempgridgroup[tempgridseg,,drop=FALSE]],na.rm=TRUE)
      if(segsum<0){segsum=0}
      fluxfrac[j]=segsum/groupsum
        
      if(deblendtype=='fit'){
        groupellip=.profoundEllipse(x=tempgridgroup[,1],y=tempgridgroup[,2],flux=image_temp[tempgridgroup],xcen=segstats[segstats$segID==segIDlist[j],"xmax"]+0.5,ycen=segstats[segstats$segID==segIDlist[j],"ymax"]+0.5,ang=segstats[segstats$segID==segIDlist[j],"ang"],axrat=segstats[segstats$segID==segIDlist[j],"axrat"])
        segout=.profoundEllipse(x=tempgridgroup[tempgridseg,1],y=tempgridgroup[tempgridseg,2],flux=image_temp[tempgridgroup[tempgridseg,,drop=FALSE]],xcen=segstats[segstats$segID==segIDlist[j],"xmax"]+0.5,ycen=segstats[segstats$segID==segIDlist[j],"ymax"]+0.5,ang=segstats[segstats$segID==segIDlist[j],"ang"],axrat=segstats[segstats$segID==segIDlist[j],"axrat"])
        
        select=which(segout[,2]>0 & is.finite(segout[,1]) & is.finite(segout[,2]))
        #Try various levels of segment fitting: spline, linear, flat
        uniq_xlen=length(unique(round(segout[select,1],6)))
        if(uniq_xlen>max(3,df) & df>=2){
          tempsafe=try(10^predict(smooth.spline(segout[select,1],log10(segout[select,2]), df=df)$fit, x=groupellip[,1])$y)
          if(class(tempsafe)=="try-error"){
            weightmatrix[,j]=10^predict(lm(y~x, data=list(x=segout[select,1],y=log10(segout[select,2]))), newdata=list(x=groupellip[,1]))
          }else{
            weightmatrix[,j]=tempsafe
          }
          weightmatrix[groupellip[,1]>radtrunc*max(segout[,1]),j]=0
          Qseg_db[j]=(sum(weightmatrix[,j])-sum(segout[select,2]))/sum(segout[select,2])
        }else if(uniq_xlen>1){
          weightmatrix[,j]=10^predict(lm(y~x, data=list(x=segout[select,1],y=log10(segout[select,2]))), newdata=list(x=groupellip[,1]))
          weightmatrix[groupellip[,1]>radtrunc*max(segout[,1]),j]=0
          Qseg_db[j]=(sum(weightmatrix[,j])-sum(segout[select,2]))/sum(segout[select,2])
        }else{
          if(segsum>0){
            weightmatrix[groupellip[,1]<=radtrunc*max(segout[,1]),j]=segsum/length(which(groupellip[,1]<=radtrunc*max(segout[,1])))
            weightmatrix[groupellip[,1]>radtrunc*max(segout[,1]),j]=0
          }else{
            weightmatrix[,j]=0
          }
          Qseg_db[j]=0
        }
        
        if(iterative){
          image_temp[tempgridgroup]=image_temp[tempgridgroup]-weightmatrix[,j]
        }
      }else if(deblendtype=='psf'){
        if(is.null(psf)){
          stop("psf must be provided if deblendtype=psf !")
        }
        if(convmode=='extended'){
          segim_temp[tempgridgroup[tempgridseg,,drop=FALSE]]=TRUE
          groupim_temp[tempgridgroup]=TRUE
          groupim_conv=ProFit::profitConvolve(convolver=convolver_brute, image=segim_temp, kernel=psf, mask=groupim_temp)
          beamcorrect[j]=beamsum/sum(groupim_conv[tempgridgroup[tempgridseg,,drop=FALSE]])
          weightmatrix[,j]=groupim_conv[tempgridgroup]/length(tempgridseg)
        }else if(convmode=='psf'){
          maxloc=which.max(image_temp[tempgridgroup[tempgridseg,,drop=FALSE]])
          addloc=as.integer(tempgridgroup[tempgridseg[maxloc],,drop=FALSE])
          addloc=addloc-c(dim(psf)[1]/2, dim(psf)[2]/2) #Offset by PSF centre
          groupim_conv=ProFit::profitAddMats(matbase=groupim_temp, matadd=psf, addloc=addloc)
          beamcorrect[j]=beamsum/sum(groupim_conv[tempgridgroup[tempgridseg,,drop=FALSE]])
          weightmatrix[,j]=groupim_conv[tempgridgroup]
        }else{
          stop('convmode type not recognised, must be one of extended / psf !')
        }
        # else if(convmode=='psf_interp'){ # No need since we can only have integer pixel source positions right now
        #   maxloc=which.max(image_temp[tempgridgroup[tempgridseg,]])
        #   addloc=as.integer(tempgridgroup[tempgridseg[maxloc],]) #True centre
        #   groupim_conv=ProFit::profitMakePointSource(xcen=addloc[1], ycen=addloc[2], matbase=groupim_temp, psf=psf, image=groupim_temp, add=TRUE)
        #   weightmatrix[,j]=groupim_conv[tempgridgroup]
        # }
        
        if(fluxweight=='sum'){
          weightmatrix[,j]=weightmatrix[,j]*segsum
        }else if(fluxweight=='peak'){
          weightmatrix[,j]=weightmatrix[,j]*max(image_temp[tempgridgroup[tempgridseg,,drop=FALSE]])
        }else if(fluxweight=='none'){
          # Nothing to do here
        }else{
          stop('fluxweight type not recognised, must be one of sum / peak / none !')
        }
          
        segim_temp[tempgridgroup[tempgridseg,,drop=FALSE]]=FALSE
        if(convmode=='extended'){
          groupim_temp[tempgridgroup]=FALSE
        }
      }else{
        stop('deblendtype type not recognised, must be one of fit / psf !')
      }
    }
    
    normmat=.rowSums(weightmatrix, dim(weightmatrix)[1], dim(weightmatrix)[2])
    weightmatrix=weightmatrix/normmat
    weightmatrix[weightmatrix<0.1]=0
    normmat=.rowSums(weightmatrix, dim(weightmatrix)[1], dim(weightmatrix)[2])
    weightmatrix=weightmatrix/normmat
    Qgroup_db=sum(normmat-image[tempgridgroup])/groupsum
    
    beamcorrect[beamcorrect<1]=1
    output[Npad:(Npad+Ngroup-1),"flux_db"]=.colSums(weightmatrix*image[tempgridgroup], dim(weightmatrix)[1], dim(weightmatrix)[2], na.rm=TRUE)
    output[Npad:(Npad+Ngroup-1),"N100_db"]=.colSums(weightmatrix, dim(weightmatrix)[1], dim(weightmatrix)[2], na.rm=TRUE)
    output[Npad:(Npad+Ngroup-1),"flux_segfrac"]=fluxfrac
    output[Npad:(Npad+Ngroup-1),"Qseg_db"]=Qseg_db
    output[Npad:(Npad+Ngroup-1),"Qgroup_db"]=Qgroup_db
    output[Npad:(Npad+Ngroup-1),"beamcorrect"]=beamcorrect
    Npad=Npad+Ngroup
  }
  
  output[,"mag_db"]=profoundFlux2Mag(flux=output[,'flux_db'], magzero=magzero)
  
  if(doallstats){
    output=output[match(segstats$segID,output$segID),,drop=FALSE]
    output[is.na(output[,"flux_segfrac"]),"flux_segfrac"]=1
    output[is.na(output[,"Qseg_db"]),"Qseg_db"]=0
    output[is.na(output[,"Qgroup_db"]),"Qgroup_db"]=0
    output[is.na(output[,"flux_db"]),c("segID", "flux_db", "mag_db", "N100_db")]=segstats[is.na(output[,"flux_db"]),c("segID", "flux","mag","N100")]
    output=cbind(output, flux_err_sky_db=segstats[,"flux_err_sky"]*sqrt(output[,'N100_db']/segstats[,'N100']))
    output=cbind(output, flux_err_skyRMS_db=segstats[,"flux_err_skyRMS"]*sqrt(output[,'N100_db']/segstats[,'N100']))
    output=cbind(output, flux_err_shot_db=segstats[,"flux_err_shot"]*suppressWarnings(sqrt(max(output[,'flux_db']/segstats[,'flux'],0,na.rm=TRUE))))
    output=cbind(output, flux_err_Qseg_db=abs(segstats[,'flux']-output[,'flux_db'])*output[,"Qseg_db"])
    output=cbind(output, flux_err_db=sqrt(output[,'flux_err_sky_db']^2+output[,'flux_err_skyRMS_db']^2+output[,'flux_err_shot_db']^2+output[,'flux_err_Qseg_db']^2))
    output=cbind(output, mag_err_db=(2.5/log(10))*abs(output[,'flux_err_db']/output[,'flux_db']))

    if(deblendtype=='psf'){
      #beam correct things not in a group too
      for(i in which(is.na(output$groupID))){
        tempgridseg=which(segim==output$segID[i],arr.ind = TRUE)
        maxloc=which.max(image_temp[tempgridseg])
        addloc=as.integer(tempgridseg[maxloc,])
        addloc=addloc-c(dim(psf)[1]/2, dim(psf)[2]/2) #Offset by PSF centre
        segim_conv=ProFit::profitAddMats(matbase=groupim_temp, matadd=psf, addloc=addloc)
        output[i,'beamcorrect']=beamsum/sum(segim_conv[tempgridseg])
      } 
    }
  }else if(iterative){
    output=output[order(output[,'groupID'],output[,'segID']),,drop=FALSE]
  }
  
  row.names(output)=NULL
  
  invisible(output)
}

profoundFitMagPSF=function(xcen=NULL, ycen=NULL, RAcen=NULL, Deccen=NULL, mag=NULL, image=NULL, im_sigma=NULL, mask=NULL, psf=NULL, fit_iters=5, magdiff=1, modxy=FALSE, sigthresh=0, itersub=TRUE, magzero=0, modelout=TRUE, fluxtype='Raw', header=NULL, doProFound=FALSE, findextra=FALSE, verbose=FALSE, ...){
  
  if(!requireNamespace("ProFit", quietly = TRUE)){
    stop('The ProFit package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  
  if((is.null(xcen) & !is.null(ycen))){stop('Need xcen/ycen pair!')}
  if((!is.null(xcen) & is.null(ycen))){stop('Need xcen/ycen pair!')}
  if((is.null(xcen) & is.null(ycen)) & (is.null(RAcen) | is.null(Deccen))){stop('Need RAcen/Decen pair!')}
  if(is.null(image)){stop('Need image!')}
  if(is.null(psf)){stop('Need psf!')}
  
  fluxtype=tolower(fluxtype)
  
  if(fluxtype=='raw' | fluxtype=='adu' | fluxtype=='adus'){
    fluxscale=1
  }else if (fluxtype=='jansky'){
    fluxscale=10^(-0.4*(magzero-8.9))
  }else{
    stop('fluxtype must be Jansky / Raw!')
  }
  
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
      modellist = list(
        pointsource = list(
          xcen = xcen,
          ycen = ycen,
          mag = mag
        )
      )
      fullmodel=ProFit::profitMakeModel(modellist=modellist, dim=dim(image), psf=psf, magzero=magzero)$z
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
      image_cut=magcutout(image, loc=c(xcen[i],ycen[i]), box=dim(psf))
      sigma_cut=magcutout(im_sigma, loc=c(xcen[i],ycen[i]), box=dim(psf))$image
      
      singlist = list(
        pointsource = list(
          xcen = image_cut$loc[1],
          ycen = image_cut$loc[2],
          mag = mag[i]
        )
      )
      
      singmodel=ProFit::profitMakeModel(modellist=singlist, dim=dim(psf), psf=psf, magzero=magzero)$z
      
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
          
          singlist = list(
          pointsource = list(
            xcen = image_cut$loc[1],
            ycen = image_cut$loc[2],
            mag = mag[i]
          )
        )
        
        if(anyNA(image_cut)){
          image_cut[selNA]=NA
          sigma_cut[selNA]=NA
        }
        
        singmodel=ProFit::profitMakeModel(modellist=singlist, dim=dim(psf), psf=psf, magzero=magzero)$z
        
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
        image[image_cut$xsel,image_cut$ysel]=image[image_cut$xsel,image_cut$ysel]-(singmodel[image_cut$xsel-image_cut$loc.diff[1], image_cut$ysel-image_cut$loc.diff[2]]*rescale)
      }
    }
    mag=mag+diffmag
  }
  
  if(modelout){
    modellist = list(
      pointsource = list(
        xcen = xcen,
         ycen = ycen,
        mag = mag
      )
    )
    fullmodel=ProFit::profitMakeModel(modellist=modellist, dim=dim(image), psf=psf, magzero=magzero)$z
    finalLL= -0.5*sum((image_orig-fullmodel)^2/(im_sigma^2))
  }else{
    fullmodel=NULL
    finalLL=NULL
  }
  flux_err=sqrt(flux_err^2 + (diffmag*flux/(2.5/log(10)))^2)
  flux_err[flux_err<beam_err]=beam_err[flux_err<beam_err]
  mag_err=(2.5/log(10))*flux_err/flux
  
  pchi=pchisq(prod(dim(psf))*(flux/flux_err)^2,prod(dim(psf))-1,log.p=TRUE)
  signif=qnorm(pchi, log.p=TRUE)
  
  if(findextra){
    if(verbose){message('Finding extra sources with ProFound')}
    protemp=profoundProFound(image_orig-fullmodel, mask=mask, magzero=magzero, header=header, verbose=FALSE, ...)
    if(!is.null(protemp$segstats)){
      xcen=c(xcen, protemp$segstats$xcen)
      ycen=c(ycen, protemp$segstats$ycen)
      mag=c(mag, protemp$segstats$mag)
      
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
  
  return(list(psfstats=psfstats, origLL=origLL, finalLL=finalLL, origmodel=origmodel, finalmodel=fullmodel, image=image_orig, psfstats_extra=psfstats_extra))
}

.minlike_mag=function(par,singmodel,image,im_sigma){
  0.5*sum((image-(singmodel*(10^(-0.4*par)-1)))^2/im_sigma^2, na.rm=TRUE)
}

.minlike_flux=function(par,singmodel,image,im_sigma){
  flux=sum(singmodel, na.rm=TRUE)
  0.5*sum((image-(singmodel*(1+par/flux)))^2/im_sigma^2, na.rm=TRUE)
}
