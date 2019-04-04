profoundMag2Mu=function(mag=15, re=1, axrat=1, pixscale=1){
  return(mag+2.5*log10(pi*re^2*axrat)-2.5*log10(0.5)+5*log10(pixscale))
}

profoundMu2Mag=function(mu=17, re=1, axrat=1, pixscale=1){
  return(mu-2.5*log10(pi*re^2*axrat)+2.5*log10(0.5)-5*log10(pixscale))
}

profoundGainConvert=function(gain=1, magzero=0, magzero_new=0){
  return(gain*10^(-0.4*(magzero_new-magzero)))
}

profoundMag2Flux=function(mag=0, magzero=0){
  return(10^(-0.4*(mag-magzero)))
}

profoundFlux2Mag=function(flux=1, magzero=0){
  flux[is.na(flux)]=0
  output=flux
  output[]=NA
  output[which(flux>0)]=-2.5*log10(flux[which(flux>0)])+magzero
  return(output)
}

profoundFlux2SB=function(flux=1, magzero=0, pixscale=1){
  return(profoundFlux2Mag(flux=flux, magzero=magzero)+5*log10(pixscale))
}

profoundSB2Flux=function(SB=0, magzero=0, pixscale=1){
  mag=SB-5*log10(pixscale)
  return(profoundMag2Flux(mag=mag, magzero=magzero))
}

profoundImBlur=function(image=NULL, sigma=1, plot=FALSE, ...){
  if(requireNamespace("imager", quietly = TRUE)){
    output=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  }else{
    if(!requireNamespace("EBImage", quietly = TRUE)){
      stop('The imager or EBImage package is needed for the profoundImBlur function to work. Please install from CRAN.', call. = FALSE)
    }
    message(" - WARNING: imager package not installed, using EBImage gblur smoothing!")
    output=as.matrix(EBImage::gblur(image,sigma))
  }
  if(plot){
    magimage(output, ...)
  }
  invisible(output)
}

profoundImGrad=function(image=NULL, sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  output=as.matrix(imager::enorm(imager::imgradient(imager::isoblur(imager::as.cimg(image),sigma), "xy")))
  if(plot){
    magimage(output, ...)
  }
  invisible(output)
}

profoundImDiff=function(image=NULL,sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  blur=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  output=image-blur
  if(plot){
    magimage(output, ...)
  }
  invisible(output)
}

profoundMakeSigma=function(image=NULL, objects=NULL, sky=0, skyRMS=0, readRMS=0, darkRMS=0, skycut=0, gain=1, image_units='ADU', sky_units='ADU', read_units='ADU', dark_units='ADU', output_units='ADU', plot=FALSE, ...){
  if(!is.null(objects)){
    if(length(objects)==length(image)){
      image[objects==0]=0
    }
  }
  if(image_units=='ADU'){
    image=gain*image
  }else if(image_units=='elec'){
    NULL
  }else{
    stop(paste('image_units unit type of',image_units,'not recognised, must be ADU or elec'))
  }
  
  if(sky_units=='ADU'){
    sky=gain*sky
    skyRMS=gain*skyRMS
  }else if(sky_units=='elec'){
    NULL
  }else{
    stop(paste('sky_units unit type of',sky_units,'not recognised, must be ADU or elec'))
  }
  
  if(read_units=='ADU'){
    readRMS=gain*readRMS
  }else if(read_units=='elec'){
    NULL
  }else{
    stop(paste('read_units unit type of',read_units,'not recognised, must be ADU or elec'))
  }
  
  if(dark_units=='ADU'){
    darkRMS=gain*darkRMS
  }else if(dark_units=='elec'){
    NULL
  }else{
    stop(paste('dark_units unit type of',dark_units,'not recognised, must be ADU or elec'))
  }
  
  image=image-sky
  image[image < skyRMS*skycut]=0
  
  if(output_units=='ADU'){
    sigma=sqrt(image+skyRMS^2+readRMS^2+darkRMS^2)/gain
  }else if(output_units=='elec'){
    sigma=sqrt(image+skyRMS^2+readRMS^2+darkRMS^2)
  }else{
    stop(paste('output_units unit type of',output_units,'not recognised, must be ADU or elec'))
  }
  
  if(plot){
    magimage(sigma, ...)
  }
  invisible(sigma)
}

profoundGainEst=function(image=NULL, mask=0, objects=0, sky=0, skyRMS=1){
  if(missing(sky)){
    sky=profoundSkyEst(image=image, mask=mask, objects=objects,plot=FALSE)$sky
  }
  image_sky=image-sky
  if(missing(skyRMS)){
    skyRMS=profoundSkyEst(image=profoundImDiff(image_sky,3), mask=mask, objects=objects,plot=FALSE)$skyRMS
  }
  tempval=as.numeric(image_sky[mask==0 & objects==0])
  
  startgain=ceiling(log10(abs(min(tempval, na.rm=T)-sky)/(skyRMS^2)))+1
  
  tempfunc=function(gain,tempval,skyRMS){
    gain=10^gain
    floor=(skyRMS*gain)^2
    trialdata=tempval*gain+floor
    value=-sum(dpois(x=round(trialdata), lambda=floor, log=T))
    invisible(value)
  }

  suppressWarnings({findgain=optim(par=startgain, fn=tempfunc, method="Brent", tempval=tempval, skyRMS=skyRMS, lower=startgain-2, upper=startgain+2)})
  invisible(10^findgain$par)
}

profoundCatMerge=function(segstats=NULL, groupstats=NULL, groupsegID=NULL, groupID_merge=NULL, flag=TRUE, rowreset=FALSE){
  if(! is.null(groupID_merge)){
    remove_segIDs=unique(unlist(groupsegID[groupsegID$groupID %in% groupID_merge,'segID']))
    remove_segIDs=remove_segIDs[!remove_segIDs %in% groupID_merge]
    segstats=segstats[! segstats$segID %in% remove_segIDs,]
    segstats[segstats$segID %in% groupID_merge,2:dim(segstats)[2]]=NA
    segstats[segstats$segID %in% groupID_merge,1:dim(groupstats)[2]]=groupstats[groupstats$groupID %in% groupID_merge,]
    segstats=segstats[order(segstats$segID),]
  }
  if(flag){
    segstats=cbind(segstats, origin='seg', stringsAsFactors=FALSE)
    segstats[segstats$segID %in% groupID_merge,'origin']='group'
  }
  if(rowreset){
    row.names(segstats)=NULL
  }
  invisible(segstats)
}

profoundFluxDeblend=function(image=NULL, segim=NULL, segstats=NULL, groupim=NULL, groupsegID=NULL, magzero=0, df=3, radtrunc=2, iterative=FALSE, doallstats=TRUE, lowmemory=FALSE){
  if(class(image)=='profound'){
    if(is.null(segim)){segim=image$segim}
    if(is.null(segstats)){segstats=image$segstats}
    if(is.null(groupim)){groupim=image$group$groupim}
    if(is.null(groupsegID)){groupsegID=image$group$groupsegID}
    if(is.null(magzero)){magzero=image$magzero}
    image=image$image-image$sky
  }
  groupsegID=groupsegID[groupsegID$Ngroup>1,,drop=FALSE]
  if(lowmemory==FALSE){
    groupID=Var1=Var2=NULL
    groupref=data.table(groupID=as.integer(groupim), expand.grid(1:dim(groupim)[1],1:dim(groupim)[2]))
    groupref[groupID %in% groupsegID$groupID]
  }
  output=data.frame(groupID=rep(groupsegID$groupID,groupsegID$Ngroup), segID=unlist(groupsegID$segID), flux_db=NA, mag_db=NA, N100_db=NA)
  if(iterative){
    output[,"flux_db"]=segstats[match(output$segID, segstats$segID),"flux"]
    output=output[order(output[,"groupID"],-output[,"flux_db"]),]
  }
  Npad=1
  image_temp=image
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
    
    groupsum=sum(image[tempgridgroup],na.rm=TRUE)
    
    for(j in 1:length(segIDlist)){
      tempgridseg=which(segim[tempgridgroup]==segIDlist[j])
      segsum=sum(image_temp[tempgridgroup[tempgridseg,]],na.rm=TRUE)
      fluxfrac[j]=segsum/groupsum
      groupellip=.profoundEllipse(x=tempgridgroup[,1],y=tempgridgroup[,2],flux=image_temp[tempgridgroup],xcen=segstats[segstats$segID==segIDlist[j],"xmax"]+0.5,ycen=segstats[segstats$segID==segIDlist[j],"ymax"]+0.5,ang=segstats[segstats$segID==segIDlist[j],"ang"],axrat=segstats[segstats$segID==segIDlist[j],"axrat"])
      segout=.profoundEllipse(x=tempgridgroup[tempgridseg,1],y=tempgridgroup[tempgridseg,2],flux=image_temp[tempgridgroup[tempgridseg,]],xcen=segstats[segstats$segID==segIDlist[j],"xmax"]+0.5,ycen=segstats[segstats$segID==segIDlist[j],"ymax"]+0.5,ang=segstats[segstats$segID==segIDlist[j],"ang"],axrat=segstats[segstats$segID==segIDlist[j],"axrat"])
      
      select=which(segout[,2]>0 & is.finite(segout[,1] & is.finite(segout[,2])))
      #Try various levels of segment fitting: spline, linear, flat
      if(length(unique(segout[select,1]))>max(3,df)){
        weightmatrix[,j]=10^predict(smooth.spline(segout[select,1],log10(segout[select,2]), df=df)$fit, x=groupellip[,1])$y
        weightmatrix[groupellip[,1]>radtrunc*max(segout[,1]),j]=0
        Qseg_db[j]=(sum(weightmatrix[,j])-sum(segout[select,2]))/sum(segout[select,2])
      }else if(length(unique(segout[select,1]))>1){
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
    }
    
    normmat=.rowSums(weightmatrix, dim(weightmatrix)[1], dim(weightmatrix)[2])
    Qgroup_db=sum(normmat-image[tempgridgroup])/groupsum
    weightmatrix=weightmatrix/normmat
    output[Npad:(Npad+Ngroup-1),"flux_db"]=.colSums(weightmatrix*image[tempgridgroup], dim(weightmatrix)[1], dim(weightmatrix)[2])
    output[Npad:(Npad+Ngroup-1),"N100_db"]=.colSums(weightmatrix, dim(weightmatrix)[1], dim(weightmatrix)[2])
    output[Npad:(Npad+Ngroup-1),"flux_segfrac"]=fluxfrac
    output[Npad:(Npad+Ngroup-1),"Qseg_db"]=Qseg_db
    output[Npad:(Npad+Ngroup-1),"Qgroup_db"]=Qgroup_db
    Npad=Npad+Ngroup
  }
  output[,"mag_db"]=profoundFlux2Mag(flux=output[,'flux_db'], magzero=magzero)
  
  if(doallstats){
    output=output[match(segstats$segID,output$segID),]
    output[is.na(output[,"flux_db"]),c("segID", "flux_db", "mag_db", "N100_db")]=segstats[is.na(output[,"flux_db"]),c("segID", "flux","mag","N100")]
    output=cbind(output, flux_err_sky_db=segstats[,"flux_err_sky"]*sqrt(output[,'N100_db']/segstats[,'N100']))
    output=cbind(output, flux_err_skyRMS_db=segstats[,"flux_err_skyRMS"]*sqrt(output[,'N100_db']/segstats[,'N100']))
    output=cbind(output, flux_err_shot_db=segstats[,"flux_err_shot"]*suppressWarnings(sqrt(output[,'flux_db']/segstats[,'flux'])))
    output=cbind(output, flux_err_db=sqrt(output[,'flux_err_sky_db']^2+output[,'flux_err_skyRMS_db']^2+output[,'flux_err_shot_db']^2))
    output=cbind(output, mag_err_db=(2.5/log(10))*abs(output[,'flux_err_db']/output[,'flux_db']))
  }else if(iterative){
    output=output[order(output[,'groupID'],output[,'segID']),]
  }
  
  invisible(output)
}

### Deprecated Functions ###

# profoundGetPixScale=function(header, CD1_1=1, CD1_2=0, CD2_1=0, CD2_2=1){
#   if(!is.null(header)){
#     if(is.data.frame(header) | is.matrix(header)){
#       locs=match(c('CD1_1','CD1_2','CD2_1','CD2_2'),header[,1])
#       headerWCS=data.frame(header[locs,1],as.numeric(header[locs,2]))
#       if('CD1_1' %in% headerWCS[,1]){
#         CD1_1=headerWCS[headerWCS[,1]=='CD1_1',2]
#         if('CD1_2' %in% headerWCS[,1]){CD1_2=headerWCS[headerWCS[,1]=='CD1_2',2]}else{message('Missing CD1_2')}
#       }else{
#         if('CDELT1' %in% headerWCS[,1]){
#           CD1_1=headerWCS[headerWCS[,1]=='CDELT1',2]
#         }else{
#           message("Missing CD1_1 and CDELT1")
#         }
#       }
#       if('CD2_2' %in% headerWCS[,1]){
#         CD2_2=headerWCS[headerWCS[,1]=='CD2_2',2]
#         if('CD2_1' %in% headerWCS[,1]){CD2_1=headerWCS[headerWCS[,1]=='CD2_1',2]}else{message('Missing CD2_1')}
#       }else{
#         if('CDELT2' %in% headerWCS[,1]){
#           CD2_2=headerWCS[headerWCS[,1]=='CDELT2',2]
#         }else{
#           message("Missing CD2_2 and CDELT2")
#         }
#       }
#     }else{
#       if('CD1_1' %in% header){
#         CD1_1=as.numeric(header[which(header=='CD1_1')+1])
#         if('CD1_2' %in% header){CD1_2=as.numeric(header[which(header=='CD1_2')+1])}else{message('Missing CD1_2')}
#       }else{
#         if('CDELT1' %in% header){
#           CD1_1=as.numeric(header[which(header=='CDELT1')+1])
#         }else{
#           message("Missing CD1_1 and CDELT1")
#         }
#       }
#       if('CD2_2' %in% header){
#         CD2_2=as.numeric(header[which(header=='CD2_2')+1])
#         if('CD2_1' %in% header){CD2_1=as.numeric(header[which(header=='CD2_1')+1])}else{message('Missing CD2_1')}
#       }else{
#         if('CDELT1' %in% header){
#           CD2_2=as.numeric(header[which(header=='CDELT2')+1])
#         }else{
#           message("Missing CD2_2 and CDELT2")
#         }
#       }
#     }
#   }
#   return(3600*(sqrt(CD1_1^2+CD1_2^2)+sqrt(CD2_1^2+CD2_2^2))/2)
# }

# profoundInterp2d=function(x,y,image){
#     scale=sum(image)
#     imagelist=list(x=seq(-dim(image)[1]/2,dim(image)[1]/2,len=dim(image)[1]),y=seq(-dim(image)[2]/2,dim(image)[2]/2,len=dim(image)[2]),z=image)
#     ximage = seq(-dim(image)[1]/2,dim(image)[1]/2,len=dim(image)[1])
#     yimage = seq(-dim(image)[2]/2,dim(image)[2]/2,len=dim(image)[2])
#     zimage = image
#     nx = length(ximage)
#     ny = length(yimage)
#     lx = approx(ximage, 1:nx, x, rule=2)$y
#     ly = approx(yimage, 1:ny, y, rule=2)$y
#     lx1 = floor(lx)
#     ly1 = floor(ly)
#     ex = lx - lx1
#     ey = ly - ly1
#     ex[lx1 == nx] = 1
#     ey[ly1 == ny] = 1
#     lx1[lx1 == nx] = nx - 1
#     ly1[ly1 == ny] = ny - 1
#     z=
# 	zimage[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) +
# 	zimage[cbind(lx1 + 1, ly1)] * ex * (1 - ey) +
# 	zimage[cbind(lx1, ly1 + 1)] * (1 - ex) * ey +
# 	zimage[cbind(lx1 + 1, ly1 + 1)] * ex * ey
#   return = cbind(X=x,Y=y,Z=z)
# }
