.ellipsesd=function(par=c(1,0,0,0), x, y, xcen, ycen, wt=1){
  x=(x-xcen)/par[1]
  y=(y-ycen)/par[1]
  return=.varwt(x=profitEllipse(x=x, y=y, flux=1, xcen=0, ycen=0, ang=par[2], axrat=10^par[3], box=par[4])[,1], xcen=1, wt=wt)
}

profitGetEllipse=function(x, y, z, xcen, ycen, scale=sqrt(2), pixscale=1, dobox=FALSE, plot=FALSE, ...){
  if(is.matrix(x)){
    if(dim(x)[2]==3){
      y=x[,2]
      z=x[,3]
      x=x[,1]
    }
  }
  if(missing(xcen)){xcen=.meanwt(x, wt=z)}
  if(missing(ycen)){ycen=.meanwt(y, wt=z)}
  xsd=sqrt(.varwt(x, wt=z))
  ysd=sqrt(.varwt(y, wt=z))
  covxy=.covarwt(x, y, wt=z)
  corxy=covxy/(xsd*ysd)
  rad=.cov2eigval(xsd, ysd, covxy)
  rad$hi=sqrt(abs(rad$hi))
  rad$lo=sqrt(abs(rad$lo))
  radav=.meanwt(sqrt((x-xcen)^2+(y-ycen)^2), z)
  axrat=rad$lo/rad$hi
  eigvec=.cov2eigvec(xsd, ysd, covxy)
  ang=.eigvec2ang(eigvec)
  if(dobox){
    lower=c(rad$hi/10,0,-2,-1)
    upper=c(rad$hi*10,180,0,1)
    tempoptim=optim(par=c(rad$hi, ang, log10(axrat), box=0), fn=.ellipsesd, x=x, y=y, xcen=xcen, ycen=ycen, wt=z, method='L-BFGS-B', lower=lower, upper=upper)$par
    if(tempoptim[1]>lower[1] & tempoptim[1]<upper[1] & tempoptim[2]>lower[2] & tempoptim[2]<upper[2] & tempoptim[3]>lower[3] & tempoptim[3]<upper[3] & tempoptim[4]>lower[4] & tempoptim[4]<upper[4]){
      rad$hi=tempoptim[1]
      ang=tempoptim[2]
      axrat=10^tempoptim[3]
      box=tempoptim[4]
      rad$lo=rad$hi*axrat
      scale=1
    }else{
      box=0
    }
  }else{
    box=0
  }
  if(plot){
    profitDrawEllipse(xcen=xcen, ycen=ycen, rad=rad$hi*scale, axrat=axrat, ang=ang, box=box, ...)
  }
  return=c(xcen=xcen, ycen=ycen, radhi=rad$hi*scale*pixscale, radlo=rad$lo*scale*pixscale, radav=radav*pixscale, axrat=axrat, ang=ang, box=box, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy)
}

profitGetEllipses=function(image, segim, segID=1, levels=10, magzero=0, pixscale=1, fixcen=TRUE, dobox=FALSE, plot=TRUE, ...){
  if(missing(segim)){segim=segID}
  tempxy=which(segim==segID, arr.ind = T)-0.5
  tempxy=cbind(tempxy,image[segim==segID])
  tempxy=tempxy[order(tempxy[,3],decreasing = T),]
  tempxy=cbind(tempxy,cumsum(tempxy[,3])/sum(tempxy[,3],na.rm=T))
  tempellipses={}
  segelllipses=matrix(0,dim(segim)[1],dim(segim)[2])
  isolevels=seq(0,1-1/levels,by=1/levels)
  if(fixcen){
    tempellipse=profitGetEllipse(tempxy[,1:3])
    xcen=as.numeric(tempellipse['xcen'])
    ycen=as.numeric(tempellipse['ycen'])
    for(i in isolevels){
      segelllipses[ceiling(tempxy[tempxy[,4]>i & tempxy[,4]<i+1/levels,1:2])]=round(i*levels+1,0)
      tempellipses=rbind(tempellipses,
                        c(profitGetEllipse(tempxy[tempxy[,4]>i & tempxy[,4]<i+1/levels,1:3], xcen=xcen, ycen=ycen, pixscale=pixscale, dobox=dobox), flux=sum(tempxy[tempxy[,4]>i & tempxy[,4]<i+1/levels,3], na.rm=T), N=length(which(tempxy[,4]>i & tempxy[,4]<i+1/levels)))
                        )
    }
  }else{
    for(i in isolevels){
      segelllipses[ceiling(tempxy[tempxy[,4]>i & tempxy[,4]<i+1/levels,1:2])]=round(i*levels+1,0)
      tempellipses=rbind(tempellipses,
                        c(profitGetEllipse(tempxy[tempxy[,4]>i & tempxy[,4]<i+1/levels,1:3], pixscale=pixscale, dobox=dobox), flux=sum(tempxy[tempxy[,4]>i & tempxy[,4]<i+0.05,3], na.rm=T), N=length(which(tempxy[,4]>i & tempxy[,4]<i+1/levels)))
                        )
    }
  }
  SB=profitFlux2SB(tempellipses[,'flux']/tempellipses[,'N'], magzero=magzero, pixscale=pixscale)
  tempellipses=cbind(segellipseID=1:length(tempellipses[,1]), fluxfrac=isolevels+1/levels, tempellipses, SB=SB)
  tempellipses=as.data.frame(tempellipses)
  if(plot){
    profitGetEllipsesPlot(image=image, ellipses=tempellipses, pixscale=pixscale, ...)
  }
  return=list(ellipses=tempellipses, segellipses=segelllipses)
}

profitEllipse=function(x, y, flux, xcen=0, ycen=0, ang=0, axrat=1, box=0){
  if(is.matrix(x)){
    z=x
    x = seq(0.5, dim(z)[1] - 0.5)
    y = seq(0.5, dim(z)[2] - 0.5)
    temp=expand.grid(x,y)
    x=temp[,1]
    y=temp[,2]
    flux=as.numeric(z)
  }
  if(!is.numeric(box)) box = 0
  rad=sqrt((x-xcen)^2+(y-ycen)^2)
  angrad=-ang*pi/180
  angmod=atan2((x-xcen),(y-ycen))-angrad
  xmod=rad*sin(angmod)
  ymod=rad*cos(angmod)
  xmod=xmod/axrat
  radmod=(abs(xmod)^(2+box)+abs(ymod)^(2+box))^(1/(2+box))
  output=cbind(rad=radmod, flux=flux)
  output=output[order(radmod),]
  return(output)
}

profitEllipsePlot=function(Data, modellist, bulgeloc=1, diskloc=2, pixscale=1, FWHM=1, SBlim=26, df=100, raw=FALSE){
  if(missing(Data)){stop('Data object of class profit.data must be provided!')}
  if(class(Data)!="profit.data"){stop("Data must be of class profit.data, as output by profitSetupData!")}
  if(missing(modellist)){modellist=Data$modellist}
  
  bulge=profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list(sersic = bulgeloc), dim = Data$imagedim, psf=Data$psf)
  disk=profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list(sersic = diskloc), dim = Data$imagedim, psf=Data$psf)
  total=profitMakeModel(modellist, magzero = Data$magzero, whichcomponents=list(sersic = 'all'), dim = Data$imagedim, psf=Data$psf)

  region = Data$region
  if(is.null(region)) region = numeric(length(Data$image)) + 1
  
  imageellipse=profitEllipse(Data$image*region, xcen=modellist$sersic$xcen[diskloc], ycen=modellist$sersic$ycen[diskloc], ang=modellist$sersic$ang[diskloc], axrat=modellist$sersic$axrat[diskloc], box=modellist$sersic$box[diskloc])
  sigmaellipse=profitEllipse(Data$sigma*region, xcen=modellist$sersic$xcen[diskloc], ycen=modellist$sersic$ycen[diskloc], ang=modellist$sersic$ang[diskloc], axrat=modellist$sersic$axrat[diskloc], box=modellist$sersic$box[diskloc])
  bulgeellipse=profitEllipse(bulge$z, xcen=modellist$sersic$xcen[bulgeloc], ycen=modellist$sersic$ycen[bulgeloc], ang=modellist$sersic$ang[bulgeloc], axrat=modellist$sersic$axrat[bulgeloc], box=modellist$sersic$box[bulgeloc])
  diskellipse=profitEllipse(disk$z, xcen=modellist$sersic$xcen[diskloc], ycen=modellist$sersic$ycen[diskloc], ang=modellist$sersic$ang[diskloc], axrat=modellist$sersic$axrat[diskloc], box=modellist$sersic$box[diskloc])
  totalellipse=profitEllipse(total$z, xcen=modellist$sersic$xcen[diskloc], ycen=modellist$sersic$ycen[diskloc], ang=modellist$sersic$ang[diskloc], axrat=modellist$sersic$axrat[diskloc], box=modellist$sersic$box[diskloc])
  psfellipse=cbind(seq(0,10*FWHM,len=1e3), dnorm(seq(0,10*FWHM,len=1e3),sd=FWHM/(2*sqrt(2*log(2)))))
  
  imageellipse[,1]=imageellipse[,1]*pixscale
  sigmaellipse[,1]=sigmaellipse[,1]*pixscale
  bulgeellipse[,1]=bulgeellipse[,1]*pixscale
  diskellipse[,1]=diskellipse[,1]*pixscale
  totalellipse[,1]=totalellipse[,1]*pixscale
  
  sigmaellipse=cbind(sigmaellipse,imageellipse[,2]-sigmaellipse[,2])
  sigmaellipse=cbind(sigmaellipse,imageellipse[,2]+sigmaellipse[,2])
  imageellipse[,2]=-2.5*suppressWarnings(log10(imageellipse[,2]))+5*log10(pixscale)+Data$magzero
  sigmaellipse[,2:4]=-2.5*suppressWarnings(log10(sigmaellipse[,2:4]))+5*log10(pixscale)+Data$magzero
  bulgeellipse[,2]=-2.5*suppressWarnings(log10(bulgeellipse[,2]))+5*log10(pixscale)+Data$magzero
  diskellipse[,2]=-2.5*suppressWarnings(log10(diskellipse[,2]))+5*log10(pixscale)+Data$magzero
  totalellipse[,2]=-2.5*suppressWarnings(log10(totalellipse[,2]))+5*log10(pixscale)+Data$magzero
  psfellipse[,2]=-2.5*suppressWarnings(log10(psfellipse[,2]))+5*log10(pixscale)+Data$magzero
  
  imageellipse[is.na(imageellipse)]=SBlim
  imageellipse[is.infinite(imageellipse)]=SBlim
  sigmaellipse[is.na(sigmaellipse)]=SBlim
  sigmaellipse[is.infinite(sigmaellipse)]=SBlim
  bulgeellipse[is.na(bulgeellipse)]=SBlim+5
  bulgeellipse[is.infinite(bulgeellipse)]=SBlim+5
  diskellipse[is.na(diskellipse)]=SBlim
  diskellipse[is.infinite(diskellipse)]=SBlim
  totalellipse[is.na(totalellipse)]=SBlim
  totalellipse[is.infinite(totalellipse)]=SBlim
  psfellipse[is.na(psfellipse)]=SBlim+5
  psfellipse[is.infinite(psfellipse)]=SBlim+5
    
  if(raw){
    predict.image=list(x=imageellipse[,1],y=imageellipse[,2])
    predict.sigma=list(x=sigmaellipse[,1],y=sigmaellipse[,2])
    predict.sigma.lo=list(x=sigmaellipse[,1],y=sigmaellipse[,3])
    predict.sigma.hi=list(x=sigmaellipse[,1],y=sigmaellipse[,4])
    predict.bulge=list(x=bulgeellipse[,1],y=bulgeellipse[,2])
    predict.disk=list(x=diskellipse[,1],y=diskellipse[,2])
    predict.total=list(x=totalellipse[,1],y=totalellipse[,2])
    predict.psf=list(x=psfellipse[,1],y=psfellipse[,2])
  }else{
    
    smooth.image=smooth.spline(imageellipse,df=df)
    smooth.sigma.mid=smooth.spline(sigmaellipse[,c(1,2)],df=df)
    smooth.sigma.lo=smooth.spline(sigmaellipse[,c(1,3)],df=df)
    smooth.sigma.hi=smooth.spline(sigmaellipse[,c(1,4)],df=df)
    smooth.bulge=smooth.spline(bulgeellipse,df=df)
    smooth.disk=smooth.spline(diskellipse,df=df)
    smooth.total=smooth.spline(totalellipse,df=df)
  
    predict.image=predict(smooth.image,imageellipse[,1])
    predict.sigma.mid=predict(smooth.sigma.mid,imageellipse[,1])
    predict.sigma.lo=predict(smooth.sigma.lo,imageellipse[,1])
    predict.sigma.hi=predict(smooth.sigma.hi,imageellipse[,1])
    predict.bulge=predict(smooth.bulge,imageellipse[,1])
    predict.disk=predict(smooth.disk,imageellipse[,1])
    predict.total=predict(smooth.total,imageellipse[,1])
  }
  
  psfellipse[,2]=psfellipse[,2]+min(predict.bulge$y)-min(psfellipse[1,2])
  
  smooth.total=smooth.spline(totalellipse,df=df)
  refpredict=predict(smooth.total,imageellipse[,1])
  
  xhi=refpredict$x[min(which(refpredict$y>SBlim))]
  yhi=min(predict.image$y,predict.total$y)
  
  sigma.polygon=rbind(cbind(predict.sigma.lo$x,predict.sigma.lo$y),
                      cbind(rev(predict.sigma.hi$x),rev(predict.sigma.hi$y))
  )

  layout(rbind(1,2), heights=c(0.7,0.4))
  par(oma=c(3.1,3.6,1.1,1.1))
  par(mar=c(0,0,0,0))
  
  if(pixscale==1){
    xlab='Project Major Axis / pix'
    ylab1=expression(mu*' (mag/pix'*''^2*')')
    ylab2=expression(Delta*mu*' (mag/pix'*''^2*')')
  }else{
    xlab='Project Major Axis / asec'
    ylab1=expression(mu*' (mag/asec'*''^2*')')
    ylab2=expression(Delta*mu*' (mag/asec'*''^2*')')
  }
  
  magplot(0,0,type='n',ylim=c(SBlim+1,yhi),xlim=c(0,xhi),col='black',xlab='', ylab=ylab1, grid=T,labels = c(F,T))
  polygon(sigma.polygon[,1],sigma.polygon[,2],col = hsv(v=0,alpha = 0.1), border=NA)
  lines(predict.image,col='black')
  lines(predict.bulge,col='red')
  lines(predict.disk,col='blue')
  lines(predict.total,col='darkgreen')
  lines(psfellipse,col='purple')
  abline(h=SBlim,lty=3)
  abline(v=xhi,lty=3)
  abline(v=FWHM/2,lty=3)
  legend('topright',legend=c('Data','ProFit Total','ProFit Bulge','ProFit Disk','PSF'),lty=1,col=c('black','darkgreen','red','blue','purple'),bg='white')
  
  magplot(predict.image$x, predict.image$y-predict.total$y, type='l', xlim=c(0,xhi), ylim=c(-0.5,0.5), col='darkgreen', xlab=xlab, ylab=ylab2, grid=T, labels=c(T,T))
  polygon(sigma.polygon[,1], sigma.polygon[,2]-c(predict.total$y,rev(predict.total$y)), col = hsv(v=0,alpha = 0.1), border=NA)
  abline(h=0)
  abline(v=FWHM/2,lty=3)
}