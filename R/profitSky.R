.interp.2d=function (x, y, obj) 
{
    if(length(x)>1e6){rembig=TRUE}else{rembig=FALSE}
    xobj = obj$x
    yobj = obj$y
    zobj = obj$z
    nx = length(xobj)
    ny = length(yobj)
    lx = approx(xobj, 1:nx, x, rule = 2)$y
    if(rembig){
      rm(x)
      invisible(gc())
    }
    ly = approx(yobj, 1:ny, y, rule = 2)$y
    if(rembig){
      rm(y)
      invisible(gc())
    }
    lx1 = floor(lx)
    ly1 = floor(ly)
    ex = lx - lx1
    if(rembig){
      rm(lx)
      invisible(gc())
    }
    ey = ly - ly1
    if(rembig){
      rm(ly)
      invisible(gc())
    }
    ex[lx1 == nx] = 1
    ey[ly1 == ny] = 1
    lx1[lx1 == nx] = nx - 1
    ly1[ly1 == ny] = ny - 1
    temp=rep(0,length(lx1))
    temp=zobj[cbind(lx1, ly1)] * (1 - ex) * (1 - ey)
    temp=temp+zobj[cbind(lx1 + 1, ly1)] * ex * (1 - ey)
    temp=temp+zobj[cbind(lx1, ly1 + 1)] * (1 - ex) * ey
    temp=temp+zobj[cbind(lx1 + 1, ly1 + 1)] * ex * ey
    return = temp
}

profitSkyEst=function(image, objects, mask, cutlo=cuthi/2, cuthi=sqrt(sum((dim(image)/2)^2)), skycut='auto', clipiters=5, radweight=0, plot=FALSE, ...){
  radweight=-radweight
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  tempref=as.matrix(expand.grid(1:xlen,1:ylen))
  xcen=xlen/2; ycen=ylen/2
  temprad=sqrt((tempref[,1]-xcen)^2+(tempref[,2]-ycen)^2)
  #Keep only pixels inside the radius bounds given by cutlo and cuthi
  keep=temprad>=cutlo & temprad<=cuthi
  #Trim
  if(!missing(mask)){
    keep=keep & mask==0
  }
  if(!missing(objects)){
    keep=keep & objects==0
  }
  tempref=tempref[keep,]
  tempval=image[tempref]
  temprad=temprad[keep]
  clip=magclip(tempval, sigma=skycut, estimate='lo')
  tempval=tempval[clip$clip]
  temprad=temprad[clip$clip]
  #Find the running medians for the data
  tempmedian=magrun(x=temprad,y=tempval,ranges=NULL,binaxis='x',Nscale=T)
  if(plot){magplot(density(tempval),...)}
  tempylims=tempmedian$ysd
  tempy=tempmedian$y
  #Calculate worst case sky error- the sd of the medians calculated
  skyerr=sd(tempy, na.rm=TRUE)
  #Gen weights to use for weighted mean sky finding. This also weights by separation from the object of interest via radweight
  weights=1/((tempmedian$x^radweight)*(tempylims[,2]-tempylims[,1])/2)^2
  #Find the weighted mean of the medians
  sky=sum(tempy*weights)/(sum(weights))
  #Now we iterate until no running medians are outside the 1-sigma bound of the sky
  select=tempylims[,1]<=sky & tempylims[,2]>=sky
  Nselect=length(which(select))
  Nselect_old=0
  while(Nselect!=Nselect_old){
    Nselect_old=length(which(select))
    newtempy=tempy[select]
    newweights=weights[select]
    sky=sum(newtempy*newweights)/(sum(newweights))
    select=tempylims[,1]<=sky & tempylims[,2]>=sky
    Nselect=length(which(select))
  }
  #Find the number of running medians that agree with the final sky within error bounds (max=10)
  Nnearsky=length(which(select))
  if(Nnearsky>=1){
    skyRMS=mean((tempylims[select,2]-tempylims[select,1])/2)*sqrt(mean(tempmedian$Nbins[select]))
  }else{
    skyRMS=mean((tempylims[,2]-tempylims[,1])/2)*sqrt(mean(tempmedian$Nbins))
  }
  if(plot){
    lines(seq(sky-5*skyRMS, sky+5*skyRMS, len=1e3), dnorm(seq(sky-5*skyRMS, sky+5*skyRMS, len=1e3), mean=sky, sd=skyRMS), col='red')
    abline(v=c(sky-skyerr,sky,sky+skyerr),lty=c(3,1,3),col='blue')
    abline(v=c(sky-skyRMS,sky+skyRMS),lty=2,col='red')
    legend('topleft', legend=c('Sky Data', 'Sky Level', 'Sky RMS'), lty=1, col=c('black','blue','red'))
  }
  return=list(sky=sky,skyerr=skyerr,skyRMS=skyRMS,Nnearsky=Nnearsky,radrun=tempmedian)
}

profitSkyEstLoc=function(image, objects, mask, loc=dim(image)/2, box=c(100,100), skytype='median', skyRMStype='quanlo', sigmasel=1, doclip=TRUE, shiftloc = TRUE, paddim = TRUE, plot=FALSE, ...){
  if(!missing(objects) | !missing(mask)){
    if(!missing(objects)){
      tempobj=magcutout(image=objects, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0
      tempobj[is.na(tempobj)]=0
    }else{
      tempobj=TRUE
    }
    if(!missing(mask)){
      tempmask=magcutout(image=mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0
      tempmask[is.na(tempmask)]=1
    }else{
      tempmask=TRUE
    }
    select=magcutout(image, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image[tempobj & tempmask]
  }else{
    select=magcutout(image, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image
  }
  if(plot){
    image=magcutout(image, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image
    imout=magimage(image, ...)
    if(!missing(mask)){
      contour(x=imout$x, y=imout$y, magcutout(mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image, add=T, col='red', drawlabels = FALSE, zlim=c(0,1), nlevels = 1)
    }
    if(!missing(objects)){
      contour(x=imout$x, y=imout$y, magcutout(objects, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image, add=T, col='blue', drawlabels = FALSE, zlim=c(0,1), nlevels = 1)
    }
  }
  if(doclip){
    suppressWarnings({clip=magclip(select, sigmasel=sigmasel, estimate = 'lo')$x})
  }else{
    clip=select
  }
  
  if(skytype=='median'){
    skyloc=median(clip, na.rm=TRUE)
  }else if(skytype=='mean'){
    skyloc=mean(clip, na.rm=TRUE)
  }else if(skytype=='mode'){
    temp=density(clip, na.rm=TRUE)
    skyloc=temp$x[which.max(temp$y)]
  }
  
  if(skyRMStype=='quanlo'){
    temp=clip-skyloc
    temp=temp[temp<0]
    skyRMSloc=abs(as.numeric(quantile(temp, pnorm(-sigmasel)*2, na.rm=TRUE)))/sigmasel
  }else if(skyRMStype=='quanhi'){
    temp=clip-skyloc
    temp=temp[temp>0]
    skyRMSloc=abs(as.numeric(quantile(temp, (pnorm(sigmasel)-0.5)*2, na.rm=TRUE)))/sigmasel
  }else if(skyRMStype=='quanboth'){
    temp=clip-skyloc
    templo=temp[temp<0]
    temphi=temp[temp>0]
    skyRMSloclo=abs(as.numeric(quantile(templo, pnorm(-sigmasel)*2, na.rm=TRUE)))/sigmasel
    skyRMSlochi=abs(as.numeric(quantile(temphi, (pnorm(sigmasel)-0.5)*2, na.rm=TRUE)))/sigmasel
    skyRMSloc=(skyRMSloclo+skyRMSlochi)/2
  }else if(skyRMStype=='sd'){
    skyRMSloc=sqrt(.varwt(clip, wt=1, xcen=skyloc))
  }
  
  return=list(val=c(skyloc, skyRMSloc), clip=clip)
}

profitMakeSkyMap=function(image, objects, mask, box=c(100,100), grid=box, skytype='median', skyRMStype='quanlo', sigmasel=1, doclip=TRUE, shiftloc = TRUE, paddim = TRUE){
  xseq=seq(grid[1]/2,dim(image)[1],by=grid[1])
  yseq=seq(grid[2]/2,dim(image)[2],by=grid[2])
  tempgrid=expand.grid(xseq, yseq)
  tempsky=matrix(0,dim(tempgrid)[1],2)
  for(i in 1:dim(tempgrid)[1]){
    tempsky[i,]=profitSkyEstLoc(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, doclip=doclip, shiftloc=shiftloc, paddim=paddim)$val
  }
  tempmat_sky=matrix(tempsky[,1],length(xseq))
  tempmat_skyRMS=matrix(tempsky[,2],length(xseq))
  tempmat_sky[is.na(tempmat_sky)]=mean(tempmat_sky, na.rm = TRUE)
  tempmat_skyRMS[is.na(tempmat_skyRMS)]=mean(tempmat_skyRMS, na.rm = TRUE)
  return=list(sky=list(x=xseq, y=yseq, z=tempmat_sky), skyRMS=list(x=xseq, y=yseq, z=tempmat_skyRMS))
}

profitMakeSkyGrid=function(image, objects, mask, box=c(100,100), grid=box, type='bilinear', skytype='median', skyRMStype='quanlo', sigmasel=1, doclip=TRUE, shiftloc = TRUE, paddim = TRUE){
  if(length(image)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  if(!requireNamespace("akima", quietly = TRUE)){
    if(type=='bicubic'){
      stop('The akima package is needed for bicubic interpolation to work. Please install it from CRAN.', call. = FALSE)
    }
  }
  xseq=seq(grid[1]/2,dim(image)[1],by=grid[1])
  yseq=seq(grid[2]/2,dim(image)[2],by=grid[2])
  tempgrid=expand.grid(xseq, yseq)
  tempsky=matrix(0,dim(tempgrid)[1],2)
  for(i in 1:dim(tempgrid)[1]){
    tempsky[i,]=profitSkyEstLoc(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, doclip=doclip, shiftloc=shiftloc, paddim=paddim)$val
  }
  tempmat_sky=matrix(tempsky[,1],length(xseq))
  tempmat_skyRMS=matrix(tempsky[,2],length(xseq))
  tempmat_sky[is.na(tempmat_sky)]=mean(tempmat_sky, na.rm = TRUE)
  tempmat_skyRMS[is.na(tempmat_skyRMS)]=mean(tempmat_skyRMS, na.rm = TRUE)
  
  if(rembig){
    invisible(gc())
  }
  
  if(dim(tempmat_sky)[1]>1){
    bigridx=rep(1:dim(image)[1]-0.5,times=dim(image)[2])
    bigridy=rep(1:dim(image)[2]-0.5,each=dim(image)[1])
    
    if(type=='bilinear'){
      temp_bi_sky=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_sky))
      temp_bi_skyRMS=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_skyRMS))
    }else if(type=='bicubic'){
      temp_bi_sky=akima::bicubic(xseq, yseq, tempmat_sky, bigridx, bigridy)$z
      temp_bi_skyRMS=akima::bicubic(xseq, yseq, tempmat_skyRMS, bigridx, bigridy)$z
    }
    
    if(rembig){
      rm(bigridx)
      rm(bigridy)
      invisible(gc())
    }
  
    temp_bi_sky=matrix(temp_bi_sky, dim(image)[1])
    temp_bi_skyRMS=matrix(temp_bi_skyRMS, dim(image)[1])
  }else{
    temp_bi_sky=matrix(tempmat_sky[1,1], dim(image)[1], dim(image)[2])
    temp_bi_skyRMS=matrix(tempmat_skyRMS[1,1], dim(image)[1], dim(image)[2])
  }
  
  return=list(sky=temp_bi_sky, skyRMS=temp_bi_skyRMS)
}