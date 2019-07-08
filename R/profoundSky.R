# .subgrid=function(dim=c(100,100), grid=c(10,10)){
#   xhimult=ceiling(dim[1]/grid[1])
#   yhimult=ceiling(dim[2]/grid[2])
#   xhipix=xhimult*grid[1]
#   yhipix=yhimult*grid[2]
#   expandlen=xhipix*yhipix
#   gridlen=prod(grid)
#   tempgrid=matrix(0,expandlen,3)
#   tempgrid[,1]=rep(rep(1:grid[1],each=grid[2]),times=expandlen/gridlen)+rep(rep(seq(0,(xhimult-1)*grid[1],by=grid[1]),each=gridlen),times=yhimult)
#   tempgrid[,2]=rep(rep(1:grid[2],times=grid[1]),times=expandlen/gridlen)+rep(rep(seq(0,(yhimult-1)*grid[2],by=grid[2]),each=gridlen),each=xhimult)
#   tempgrid[,3]=rep(1:(xhimult*yhimult),each=gridlen)
#   tempgrid=tempgrid[tempgrid[,1]<=dim[1] & tempgrid[,2]<=dim[2],]
#   #tempmat=matrix(0,dim[1],dim[2])
#   tempgrid=tempgrid[order(tempgrid[,2],tempgrid[,1]),3]
#   #tempmat[]=tempgrid[,3]
#   invisible(tempgrid)
# }

.interp.2d.akima=function (x, y, z, xo, yo, ncp = 0, extrap = FALSE, duplicate = "error",
    dupfun = NULL){
    if (!(all(is.finite(x)) && all(is.finite(y)) && all(is.finite(z))))
        stop("missing values and Infs not allowed")
    if (is.null(xo))
        stop("xo missing")
    if (is.null(yo))
        stop("yo missing")
    if (ncp > 25) {
        ncp <- 25
        cat("ncp too large, using ncp=25\n")
    }
    drx <- diff(range(x))
    dry <- diff(range(y))
    if (drx == 0 || dry == 0)
        stop("all data collinear")
    if (drx/dry > 10000 || drx/dry < 1e-04)
        stop("scales of x and y are too dissimilar")
    n <- length(x)
    np <- length(xo)
    if (length(yo) != np)
        stop("length of xo and yo differ")
    if (length(y) != n || length(z) != n)
        stop("Lengths of x, y, and z do not match")
    xy <- paste(x, y, sep = ",")
    i <- match(xy, xy)
    if (duplicate == "user" && !is.function(dupfun))
        stop("duplicate=\"user\" requires dupfun to be set to a function")
    if (duplicate != "error") {
        centre <- function(x) {
            switch(duplicate, mean = mean(x), median = median(x),
                user = dupfun(x))
        }
        if (duplicate != "strip") {
            z <- unlist(lapply(split(z, i), centre))
            ord <- !duplicated(xy)
            x <- x[ord]
            y <- y[ord]
            n <- length(x)
        }
        else {
            ord <- (hist(i, plot = FALSE, freq = TRUE, breaks = seq(0.5,
                max(i) + 0.5, 1))$counts == 1)
            x <- x[ord]
            y <- y[ord]
            z <- z[ord]
            n <- length(x)
        }
    }
    else if (any(duplicated(xy)))
        stop("duplicate data points")
    zo <- rep(0, np)
    storage.mode(zo) <- "double"
    miss <- !extrap
    misso <- seq(miss, np)
    if (extrap & ncp == 0)
        warning("Cannot extrapolate with linear option")
    ans <- .Fortran("idbvip", as.integer(1), as.integer(ncp),
        as.integer(n), as.double(x), as.double(y), as.double(z),
        as.integer(np), x = as.double(xo), y = as.double(yo),
        z = zo, integer((31 + ncp) * n + np), double(8 * n),
        misso = as.logical(misso), PACKAGE = "akima")
    temp <- ans[c("x", "y", "z", "misso")]
    temp$z[temp$misso] <- NA
    temp[c("x", "y", "z")]
}

# .quickclip=function(flux){
#   sel=magclip(flux, estimate='lo')$x
#   invisible(list(sky=median(sel, na.rm=TRUE), skyRMS=sd(sel, na.rm=TRUE)))
# }

profoundSkyEst=function(image=NULL, objects=NULL, mask=NULL, cutlo=cuthi/2, cuthi=sqrt(sum((dim(image)/2)^2)), skycut='auto', clipiters=5, radweight=0, plot=FALSE, ...){
  radweight=-radweight
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  tempref=as.matrix(expand.grid(1:xlen,1:ylen))
  xcen=xlen/2; ycen=ylen/2
  temprad=sqrt((tempref[,1]-xcen)^2+(tempref[,2]-ycen)^2)
  #Keep only pixels inside the radius bounds given by cutlo and cuthi
  if(!is.null(mask)){
    keep=temprad>=cutlo & temprad<=cuthi & mask==0
  }else{
    keep=TRUE
  }
  if(!is.null(objects)){
    keep=temprad>=cutlo & temprad<=cuthi & objects==0 & keep
  }
  tempref=tempref[keep,]
  tempval=image[tempref]
  temprad=temprad[keep]
  clip=magclip(tempval, sigma=skycut, estimate='lo')$clip
  tempval=tempval[clip]
  temprad=temprad[clip]
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
  invisible(list(sky=sky,skyerr=skyerr,skyRMS=skyRMS,Nnearsky=Nnearsky,radrun=tempmedian))
}

profoundSkyEstLoc=function(image=NULL, objects=NULL, mask=NULL, loc=dim(image)/2, box=c(100,100), skytype='median', skyRMStype='quanlo', sigmasel=1, skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, doclip=TRUE, shiftloc = FALSE, paddim = TRUE, plot=FALSE, ...){
  if(!is.null(objects) | !is.null(mask)){
    # if(!is.null(objects)){
    #   tempobj=magcutout(image=objects, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0
    #   tempobj[is.na(tempobj)]=0
    # }else{
    #   tempobj=TRUE
    # }
    # if(!is.null(mask)){
    #   tempmask=magcutout(image=mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0
    #   tempmask[is.na(tempmask)]=1
    # }else{
    #   tempmask=TRUE
    # }
    skyN=0
    iterN=0
    tempcomb={}
    while(skyN<skypixmin & iterN<=boxiters){
      if(!is.null(objects)){
        tempcomb=magcutout(image=objects, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0
        if(!is.null(mask)){
          tempcomb=tempcomb & (magcutout(image=mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0)
        }
      }else{
        tempcomb=magcutout(image=mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image==0
      }
      tempcomb[is.na(tempcomb)]=FALSE
      if(!is.null(tempcomb)){
        tempcomb=which(tempcomb)
        skyN=length(tempcomb)
      }else{
        skyN=0
      }
      box=box+boxadd
      iterN=iterN+1
    }
    box=box-boxadd #since one too many boxadds will have occurred when it terminates

    if(skyN>0){
      select=magcutout(image, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image[tempcomb]
    }else{
      select=NA
    }
  }else{
    select=magcutout(image, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image
  }
  if(plot){
    image=magcutout(image, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image
    imout=magimage(image, ...)
    if(!is.null(mask)){
      contour(x=imout$x, y=imout$y, magcutout(mask, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image, add=T, col='red', drawlabels = FALSE, zlim=c(0,1), nlevels = 1)
    }
    if(!is.null(objects)){
      contour(x=imout$x, y=imout$y, magcutout(objects, loc=loc, box=box, shiftloc=shiftloc, paddim=paddim)$image, add=T, col='blue', drawlabels = FALSE, zlim=c(0,1), nlevels = 1)
    }
  }
  if(doclip){
    suppressWarnings({clip=magclip(select, sigmasel=sigmasel, estimate = 'lo', extra=FALSE)$x})
  }else{
    clip=select
  }
  
  if(length(clip)==1){
    if(is.na(clip)){
      return(invisible(list(val=c(NA, NA), clip=NA)))
    }
  }
  
  if(skytype=='median'){
    if('Rfast' %in% .packages()){
      skyloc=try(Rfast::med(clip, na.rm=TRUE), silent=TRUE)
      if(class(skyloc)=='try-error'){skyloc=NA}
    }else{
      skyloc=stats::median(clip, na.rm=TRUE)
    }
  }else if(skytype=='mean'){
    skyloc=mean(clip, na.rm=TRUE)
  }else if(skytype=='mode'){
    temp=density(clip, na.rm=TRUE)
    skyloc=temp$x[which.max(temp$y)]
  }
  
  if(skyRMStype=='quanlo'){
    temp=clip-skyloc
    temp=temp[temp<0]
    skyRMSloc=abs(quantile(temp, pnorm(-sigmasel)*2,na.rm=TRUE))/sigmasel
  }else if(skyRMStype=='quanhi'){
    temp=clip-skyloc
    temp=temp[temp>0]
    skyRMSloc=abs(quantile(temp, (pnorm(sigmasel)-0.5)*2,na.rm=TRUE))/sigmasel
  }else if(skyRMStype=='quanboth'){
    temp=clip-skyloc
    templo=temp[temp<0]
    temphi=temp[temp>0]
    skyRMSloclo=abs(quantile(templo, pnorm(-sigmasel)*2,na.rm=TRUE))/sigmasel
    skyRMSlochi=abs(quantile(temphi, (pnorm(sigmasel)-0.5)*2,na.rm=TRUE))/sigmasel
    skyRMSloc=(skyRMSloclo+skyRMSlochi)/2
  }else if(skyRMStype=='sd'){
    skyRMSloc=sqrt(.varwt(clip, wt=1, xcen=skyloc))
  }
  
  return(invisible(list(val=c(skyloc, skyRMSloc), clip=clip)))
}

profoundMakeSkyMap=function(image=NULL, objects=NULL, mask=NULL, box=c(100,100), grid=box, skytype='median', skyRMStype='quanlo', sigmasel=1, skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, doclip=TRUE, shiftloc = FALSE, paddim = TRUE, cores=1){
  xseq=seq(grid[1]/2,dim(image)[1],by=grid[1])
  yseq=seq(grid[2]/2,dim(image)[2],by=grid[2])
  tempgrid=expand.grid(xseq, yseq)
  registerDoParallel(cores=cores)
  
  if(cores>1){
    registerDoParallel(cores=cores)
    i=NULL
    tempsky=foreach(i = 1:dim(tempgrid)[1], .combine='rbind')%dopar%{
      profoundSkyEstLoc(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim)$val
    }
    tempsky=rbind(tempsky)
  }else{
    tempsky=matrix(0,dim(tempgrid)[1],2)
    for(i in 1:dim(tempgrid)[1]){
      tempsky[i,]=profoundSkyEstLoc(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim)$val
    }
  }
  
  tempmat_sky=matrix(tempsky[,1],length(xseq))
  tempmat_skyRMS=matrix(tempsky[,2],length(xseq))
  #tempmat_sky[is.na(tempmat_sky)]=median(tempmat_sky, na.rm = TRUE)
  #tempmat_skyRMS[is.na(tempmat_skyRMS)]=median(tempmat_skyRMS, na.rm = TRUE)
  invisible(list(sky=list(x=xseq, y=yseq, z=tempmat_sky), skyRMS=list(x=xseq, y=yseq, z=tempmat_skyRMS)))
}

profoundMakeSkyGrid=function(image=NULL, objects=NULL, mask=NULL, box=c(100,100), grid=box, type='bicubic', skytype='median', skyRMStype='quanlo', sigmasel=1, skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, doclip=TRUE, shiftloc = FALSE, paddim = TRUE, cores=1){
  if(!requireNamespace("akima", quietly = TRUE)){
    if(type=='bicubic'){
      stop('The akima package is needed for bicubic interpolation to work. Please install it from CRAN.', call. = FALSE)
    }
    if(type=='bilinear'){
      useakima=FALSE
    }
  }else{
    useakima=TRUE
  }
  
  if(box[1]>dim(image)[1]){box[1]=dim(image)[1]}
  if(box[2]>dim(image)[2]){box[2]=dim(image)[2]}
  if(grid[1]>dim(image)[1]){grid[1]=dim(image)[1]}
  if(grid[2]>dim(image)[2]){grid[2]=dim(image)[2]}
  
  xseq=seq(grid[1]/2,dim(image)[1],by=grid[1])
  yseq=seq(grid[2]/2,dim(image)[2],by=grid[2])
  tempgrid=expand.grid(xseq, yseq)
  
  if(cores>1){
    registerDoParallel(cores=cores)
    i=NULL
    tempsky=foreach(i = 1:dim(tempgrid)[1], .combine='rbind')%dopar%{
      profoundSkyEstLoc(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim)$val
    }
    tempsky=rbind(tempsky)
  }else{
    tempsky=matrix(0,dim(tempgrid)[1],2)
    for(i in 1:dim(tempgrid)[1]){
      tempsky[i,]=profoundSkyEstLoc(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, doclip=doclip, shiftloc=shiftloc, paddim=paddim)$val
    }
  }
  
  xseq=c(-grid[1]/2,xseq,max(xseq)+grid[1]/2)
  yseq=c(-grid[2]/2,yseq,max(yseq)+grid[2]/2)
  
  tempmat_sky=matrix(0,length(xseq),length(yseq))
  tempmat_sky[2:(length(xseq)-1),2:(length(yseq)-1)]=tempsky[,1]
  tempmat_sky[is.na(tempmat_sky)]= stats::median(tempmat_sky, na.rm = TRUE)
  
  tempmat_skyRMS=matrix(0,length(xseq),length(yseq))
  tempmat_skyRMS[2:(length(xseq)-1),2:(length(yseq)-1)]=tempsky[,2]
  tempmat_skyRMS[is.na(tempmat_skyRMS)]=stats::median(tempmat_skyRMS, na.rm = TRUE)
  
  xstart=min(3,dim(tempmat_sky)[1]-1)
  ystart=min(3,dim(tempmat_sky)[2]-1)
  xend=max(length(xseq)-2,2)
  yend=max(length(yseq)-2,2)
  
  tempmat_sky[1,]=tempmat_sky[2,]*2-tempmat_sky[xstart,]
  tempmat_sky[length(xseq),]=tempmat_sky[length(xseq)-1,]*2-tempmat_sky[xend,]
  tempmat_sky[,1]=tempmat_sky[,2]*2-tempmat_sky[,ystart]
  tempmat_sky[,length(yseq)]=tempmat_sky[,length(yseq)-1]*2-tempmat_sky[,yend]
  
  tempmat_skyRMS[1,]=tempmat_skyRMS[2,]*2-tempmat_skyRMS[xstart,]
  tempmat_skyRMS[length(xseq),]=tempmat_skyRMS[length(xseq)-1,]*2-tempmat_skyRMS[xend,]
  tempmat_skyRMS[,1]=tempmat_skyRMS[,2]*2-tempmat_skyRMS[,ystart]
  tempmat_skyRMS[,length(yseq)]=tempmat_skyRMS[,length(yseq)-1]*2-tempmat_skyRMS[,yend]
  
  if(dim(tempmat_sky)[1]>1){
    
    #expand out map here!! and then use akima::bilinear function
    
    bigridx=rep(1:dim(image)[1]-0.5,times=dim(image)[2])
    bigridy=rep(1:dim(image)[2]-0.5,each=dim(image)[1])
    
    if(type=='bilinear'){
      if(useakima){
        tempgrid=expand.grid(xseq, yseq)
        temp_bi_sky=.interp.2d.akima(x=tempgrid[,1], y=tempgrid[,2], z=as.numeric(tempmat_sky),xo=bigridx, yo=bigridy)$z
        temp_bi_skyRMS=.interp.2d.akima(x=tempgrid[,1], y=tempgrid[,2], z=as.numeric(tempmat_skyRMS),xo=bigridx, yo=bigridy)$z
      }else{
        temp_bi_sky=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_sky))
        temp_bi_skyRMS=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_skyRMS))
      }
    }else if(type=='bicubic'){
      temp_bi_sky=akima::bicubic(xseq, yseq, tempmat_sky, bigridx, bigridy)$z
      temp_bi_skyRMS=akima::bicubic(xseq, yseq, tempmat_skyRMS, bigridx, bigridy)$z
    }else{
      stop('type must be one of bilinear / bicubic !')
    }
    
    rm(bigridx)
    rm(bigridy)
  
    temp_bi_sky=matrix(temp_bi_sky, dim(image)[1], dim(image)[2])
    temp_bi_skyRMS=matrix(temp_bi_skyRMS, dim(image)[1], dim(image)[2])
  }else{
    temp_bi_sky=matrix(tempmat_sky[1,1], dim(image)[1], dim(image)[2])
    temp_bi_skyRMS=matrix(tempmat_skyRMS[1,1], dim(image)[1], dim(image)[2])
  }
  
  if(!is.null(mask)){
    temp_bi_sky[mask>0]=NA
    temp_bi_skyRMS[mask>0]=NA
  }
  
  invisible(list(sky=temp_bi_sky, skyRMS=temp_bi_skyRMS))
}

profoundMakeSkyBlur=function(image=NULL, objects=NULL, box=100, sigma=box*(4/pi)/sqrt(12)){
  if(!is.null(objects)){
    image[objects==1]=NA
  }
  if(requireNamespace("imager", quietly = TRUE)){
    invisible(as.matrix(imager::isoblur(imager::as.cimg(image), sigma=sigma, na.rm=TRUE)))
  }else{
    stop('The imager package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
  }
}

#Alas, the quick function does not appear to be quicker than the current profoundMakeSkyGrid function. Oh well, worth try.

# .profoundQuickSky=function(image, box=c(100,100)){
#   tempIDs=.subgrid(dim(image), grid=box)
#   tempDT=data.table(flux=as.numeric(image), subset=tempIDs)
#   setkey(tempDT, subset)
#   flux=NULL
#   invisible(tempDT[,.quickclip(flux),by=subset])
# }
