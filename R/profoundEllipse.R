.ellipsesd=function(par=c(1,0,0,0), x, y, xcen, ycen, wt=1){
  x=(x-xcen)/par[1]
  y=(y-ycen)/par[1]
  invisible(.varwt(x=.profoundEllipse(x=x, y=y, flux=1, xcen=0, ycen=0, ang=par[2], axrat=10^par[3], box=par[4])[,1], xcen=1, wt=wt))
}

.profoundEllipse=function(x, y, flux, xcen=0, ycen=0, ang=0, axrat=1, box=0){
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
  if(box==0){
    radmod=sqrt(xmod^2 + ymod^2)
  }else{
    radmod=(abs(xmod)^(2+box)+abs(ymod)^(2+box))^(1/(2+box))
  }
  output=cbind(rad=radmod, flux=flux)
  invisible(output)
}

profoundGetEllipse=function(x, y, z, xcen=NULL, ycen=NULL, scale=sqrt(2), pixscale=1, dobox=FALSE, plot=FALSE, ...){
  if(is.matrix(x)){
    if(dim(x)[2]==3){
      y=x[,2]
      z=x[,3]
      x=x[,1]
    }
  }
  if(is.null(xcen)){xcen=.meanwt(x, wt=z)}
  if(is.null(ycen)){ycen=.meanwt(y, wt=z)}
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
    profoundDrawEllipse(xcen=xcen, ycen=ycen, rad=rad$hi*scale, axrat=axrat, ang=ang, box=box, ...)
  }
  invisible(c(xcen=xcen, ycen=ycen, radhi=rad$hi*scale*pixscale, radlo=rad$lo*scale*pixscale, radav=radav*pixscale, axrat=axrat, ang=ang, box=box, xsd=xsd, ysd=ysd, covxy=covxy, corxy=corxy))
}

profoundGetEllipses=function(image=NULL, segim=NULL, segID=1, levels=10, magzero=0, pixscale=1, fixcen=TRUE, dobox=FALSE, plot=TRUE, ...){
  if(is.null(segim)){segim=segID}
  tempxy=which(segim==segID, arr.ind = T)-0.5
  tempxy=cbind(tempxy,image[segim==segID])
  tempxy=tempxy[order(tempxy[,3],decreasing = T),]
  tempxy=cbind(tempxy,cumsum(tempxy[,3])/sum(tempxy[,3],na.rm=T))
  tempellipses={}
  segelllipses=matrix(0,dim(segim)[1],dim(segim)[2])
  if(length(levels)>1){
    isolevels=levels
    difflevels=diff(isolevels)
    levels=length(levels)
  }else{
    isolevels=seq(0,1,by=1/levels)
    difflevels=diff(isolevels)
  }
  if(fixcen){
    tempellipse=profoundGetEllipse(tempxy[,1:3])
    xcen=as.numeric(tempellipse['xcen'])
    ycen=as.numeric(tempellipse['ycen'])
    for(i in 1:(length(isolevels)-1)){
      segelllipses[ceiling(tempxy[tempxy[,4]>isolevels[i] & tempxy[,4]<isolevels[i]+difflevels[i],1:2])]=i
      tempellipses=rbind(tempellipses,
                        c(profoundGetEllipse(tempxy[tempxy[,4]>isolevels[i]  & tempxy[,4]<isolevels[i]+difflevels[i],1:3], xcen=xcen, ycen=ycen, pixscale=pixscale, dobox=dobox), flux=sum(tempxy[tempxy[,4]>isolevels[i]  & tempxy[,4]<isolevels[i]+difflevels[i],3], na.rm=T), N=length(which(tempxy[,4]>isolevels[i] & tempxy[,4]<isolevels[i]+difflevels[i])))
                        )
    }
  }else{
    for(i in 1:(length(isolevels)-1)){
      segelllipses[ceiling(tempxy[tempxy[,4]>isolevels[i]  & tempxy[,4]<isolevels[i]+difflevels[i],1:2])]=i
      tempellipses=rbind(tempellipses,
                        c(profoundGetEllipse(tempxy[tempxy[,4]>isolevels[i]  & tempxy[,4]<isolevels[i]+difflevels[i],1:3], pixscale=pixscale, dobox=dobox), flux=sum(tempxy[tempxy[,4]>isolevels[i] & tempxy[,4]<isolevels[i]+difflevels[i],3], na.rm=T), N=length(which(tempxy[,4]>isolevels[i] & tempxy[,4]<isolevels[i]+difflevels[i])))
                        )
    }
  }
  SB=profoundFlux2SB(tempellipses[,'flux']/tempellipses[,'N'], magzero=magzero, pixscale=pixscale)
  tempellipses=cbind(segellipseID=1:length(tempellipses[,1]), fluxfrac=isolevels[2:length(isolevels)], tempellipses, SB=SB)
  tempellipses=as.data.frame(tempellipses)
  if(plot){
    profoundGetEllipsesPlot(image=image, ellipses=tempellipses, segim=segim, segID=segID, pixscale=pixscale, ...)
  }
  invisible(list(ellipses=tempellipses, segellipses=segelllipses))
}

profoundGetEllipsesPlot=function(image=NULL, ellipses=NULL, segim=NULL, segID=1, segellipseID='all', pixscale=1, col=rep(rainbow(10,s=0.5),4), border='auto', lty='auto', lwd='auto', ...){
  tempcon = magimage(image, col=col, ...)
  if(segellipseID[1]=='all'){segellipseID=1:length(ellipses[,1])}
  for(i in segellipseID){
    if(round(ellipses[ellipses$segellipseID==i,'fluxfrac'],2)<0.5){
      if(border=='auto'){tempborder='black'}else{tempborder=border}
      if(lty=='auto'){templty=1}else{templty=lty}
      if(lwd=='auto'){templwd=0.5}else{templwd=lwd}
    }else if(round(ellipses[ellipses$segellipseID==i,'fluxfrac'],2)==0.5){
      if(border=='auto'){tempborder='black'}else{tempborder=border}
      if(lty=='auto'){templty=1}else{templty=lty}
      if(lwd=='auto'){templwd=2}else{templwd=lwd}
    }else if(round(ellipses[ellipses$segellipseID==i,'fluxfrac'],2)>0.5 & round(ellipses[ellipses$segellipseID==i,'fluxfrac'],2)<0.9){
      if(border=='auto'){tempborder='black'}else{tempborder=border}
      if(lty=='auto'){templty=1}else{templty=lty}
      if(lwd=='auto'){templwd=1}else{templwd=lwd}
    }else if(round(ellipses[ellipses$segellipseID==i,'fluxfrac'],2)==0.9){
      if(border=='auto'){tempborder='black'}else{tempborder=border}
      if(lty=='auto'){templty=1}else{templty=lty}
      if(lwd=='auto'){templwd=2}else{templwd=lwd}
    }else if(round(ellipses[ellipses$segellipseID==i,'fluxfrac'],2)>0.9){
      if(border=='auto'){tempborder='black'}else{tempborder=border}
      if(lty=='auto'){templty=2}else{templty=lty}
      if(lwd=='auto'){templwd=1}else{templwd=lwd}
    }
    if(!is.null(segim)){
      tempcon = magimage(1-(segim==segID), add=T, magmap=F, zlim=c(0,1), col=NA)
      contour(tempcon, add=T, drawlabels=F, levels=1, col = "darkgreen")
    }
    profoundDrawEllipse(xcen=ellipses[ellipses$segellipseID==i,'xcen'], ycen=ellipses[ellipses$segellipseID==i,'ycen'], rad=ellipses[ellipses$segellipseID==i,'radhi']/pixscale, axrat=ellipses[ellipses$segellipseID==i,'axrat'], ang=ellipses[ellipses$segellipseID==i,'ang'], box=ellipses[ellipses$segellipseID==i,'box'], col=tempborder, lty=templty, lwd=templwd)
  }
}

profoundDrawEllipse=function(xcen=0, ycen=0, rad=1, axrat=1, ang=0, box=0, ...){
  tempellipse=data.frame(xcen=xcen, ycen=ycen, rad=rad, axrat=axrat, ang=ang, box=box)
  tempcirc=seq(0, 2*pi, length=100)
  xtemp=cos(tempcirc)
  ytemp=sin(tempcirc)
  for(i in 1:length(xcen)){
    boxscale=1/((abs(xtemp)^(2+tempellipse$box[i])+abs(ytemp)^(2+tempellipse$box[i]))^(1/(2+tempellipse$box[i])))
    angtemp = (pi/180)*(tempellipse$ang[i]+90)
    radlo = tempellipse$rad[i]*tempellipse$axrat[i]
    radhi = tempellipse$rad[i]
    x = tempellipse$xcen[i] + boxscale*(radhi * cos(tempcirc) * cos(angtemp) - radlo * sin(tempcirc) * sin(angtemp))
    y = tempellipse$ycen[i] + boxscale*(radhi * cos(tempcirc) * sin(angtemp) + radlo * sin(tempcirc) * cos(angtemp))
    lines(x=x, y=y, ...)
  }
}
