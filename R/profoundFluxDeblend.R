profoundFluxDeblend=function(image=NULL, segim=NULL, segstats=NULL, groupim=NULL, groupsegID=NULL, sky=0, profound=NULL, magzero=0, df=3, radtrunc=2, iterative=FALSE, doallstats=TRUE, lowmemory=FALSE, deblendtype='fit', psf=NULL, fluxweight='sum', convtype='brute', convmode='extended'){
  
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
    group=profoundSegimGroup(segim)
    groupim=group$groupim
    groupsegID=group$groupsegID
  }
  
  image=image-sky
  
  groupsegID=groupsegID[groupsegID$Ngroup>1,,drop=FALSE]
  
  if(lowmemory==FALSE){
    groupID=Var1=Var2=NULL
    groupref=data.table(groupID=as.integer(groupim), expand.grid(1:dim(groupim)[1],1:dim(groupim)[2]), keyby='groupID')
    setkey(groupref, groupID)
    groupref[groupID %in% groupsegID$groupID]
  }else{
    rm(sky)
    invisible(gc())
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
      tempgridgroup=as.matrix(cbind(groupref[groupID==groupsegID[i,"groupID"],list(Var1,Var2)])) 
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
        if(uniq_xlen>max(3,df)){
          weightmatrix[,j]=10^predict(smooth.spline(segout[select,1],log10(segout[select,2]), df=df)$fit, x=groupellip[,1])$y
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
    Qgroup_db=sum(normmat-image[tempgridgroup])/groupsum
    weightmatrix=weightmatrix/normmat
    beamcorrect[beamcorrect<1]=1
    output[Npad:(Npad+Ngroup-1),"flux_db"]=.colSums(weightmatrix*image[tempgridgroup], dim(weightmatrix)[1], dim(weightmatrix)[2])
    output[Npad:(Npad+Ngroup-1),"N100_db"]=.colSums(weightmatrix, dim(weightmatrix)[1], dim(weightmatrix)[2])
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

profoundFitMagPSF=function(xcen=NULL, ycen=NULL, mag=NULL, image=NULL, sigma=NULL, mask=NULL, psf=NULL, iters=5, magdiff=1, modxy=FALSE, sigthresh=0, itersub=TRUE, magzero=0, modelout=TRUE, verbose=FALSE){
  
  if(!requireNamespace("ProFit", quietly = TRUE)){
    stop('The ProFit package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  
  if(is.null(xcen)){stop('Need xcen!')}
  if(is.null(ycen)){stop('Need ycen!')}
  if(is.null(mag)){stop('Need mag!')}
  if(is.null(image)){stop('Need image!')}
  if(is.null(sigma)){stop('Need sigma!')}
  if(is.null(psf)){stop('Need psf!')}
  
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
  
  if(modxy){
    xygrid=expand.grid(1:dim(psf)[1],1:dim(psf)[2])-0.5
  }
  
  diffmag=rep(0,Nmodels)
  mag_err=rep(0,Nmodels)
  psfLL=rep(0,Nmodels)
  
  for(j in 1:iters){
    if(verbose){message(paste('Iteration',j,'of',iters))}
    modellist = list(
      pointsource = list(
        xcen = xcen,
        ycen = ycen,
        mag = mag
      )
    )
    
    if(itersub){
      if(j==1){
        fullmodel=ProFit::profitMakeModel(modellist=modellist, dim=dim(image), psf=psf, magzero=magzero)$z
        image=image_orig-fullmodel
      }
    }else{
      fullmodel=ProFit::profitMakeModel(modellist=modellist, dim=dim(image), psf=psf, magzero=magzero)$z
      image=image_orig-fullmodel
    }
    
    if(j==1){
      if(modelout){
        origmodel=fullmodel
        origLL= sum(-0.5*(image_orig-fullmodel)^2/(sigma^2))
      }else{
        origmodel=NA
        origLL=NA
      }
    }
    
    for(i in 1:Nmodels){
      #if(i==417){browser()}
      image_cut=magcutout(image, loc=c(xcen[i],ycen[i]), box=dim(psf))
      sigma_cut=magcutout(sigma, loc=c(xcen[i],ycen[i]), box=dim(psf))$image
      
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
      
      if(j<iters){
        diffmag[i] = optim(0, .minlike, method='Brent', singmodel = singmodel, image=image_cut$image, sigma = sigma_cut, lower=-magdiff, upper=magdiff)$par
      }else{
        finaloptim=optim(0, .minlike, method='Brent', singmodel = singmodel, image=image_cut$image, sigma = sigma_cut, lower=-magdiff, upper=magdiff, hessian=TRUE)
        diffmag[i]=finaloptim$par
        mag_err[i]=sqrt(1/abs(as.numeric(finaloptim$hessian)))
        psfLL[i]=-0.5*finaloptim$value
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
    finalLL= -0.5*sum((image_orig-fullmodel)^2/(sigma^2))
  }else{
    fullmodel=NA
    finalLL=NA
  }
  
  mag_err=sqrt(diffmag^2+mag_err^2)
  flux=profoundMag2Flux(mag=mag, magzero=magzero)
  flux_err=flux*mag_err/(2.5/log(10))
  
  pchi=pchisq(prod(dim(psf))*(flux/flux_err)^2,prod(dim(psf))-1,log.p=TRUE)
  signif=qnorm(pchi, log.p=TRUE)
  
  psfstats=data.frame(xcen=xcen, ycen=ycen, flux=flux, flux_err=flux_err, mag=mag, mag_err=mag_err, psfLL=psfLL, signif=signif)
  
  return(list(psfstats=psfstats, origLL=origLL, finalLL=finalLL, origmodel=origmodel, finalmodel=fullmodel))
}

.minlike=function(par,singmodel,image,sigma){
  sum((image-(singmodel*(10^(-0.4*par)-1)))^2/sigma^2,na.rm = TRUE)
}
