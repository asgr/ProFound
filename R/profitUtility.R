profitMag2Mu=function(mag=15, re=1, axrat=1, pixscale=1){
  return(mag+2.5*log10(pi*re^2*axrat)-2.5*log10(0.5)+5*log10(pixscale))
}

profitMu2Mag=function(mu=17, re=1, axrat=1, pixscale=1){
  return(mu-2.5*log10(pi*re^2*axrat)+2.5*log10(0.5)-5*log10(pixscale))
}

profitGainConvert=function(gain=1, magzero=0, magzero_new=0){
  return(gain*10^(-0.4*(magzero_new-magzero)))
}

profitMag2Flux=function(mag=0, magzero=0){
  return(10^(-0.4*(mag-magzero)))
}

profitFlux2Mag=function(flux=1, magzero=0){
  flux[is.na(flux)]=0
  output=rep(NA,length(flux))
  output[which(flux>0)]=-2.5*log10(flux[which(flux>0)])+magzero
  return(output)
}

profitFlux2SB=function(flux=1, magzero=0, pixscale=1){
  return(profitFlux2Mag(flux=flux, magzero=magzero)+5*log10(pixscale))
}

profitSB2Flux=function(SB=0, magzero=0, pixscale=1){
  mag=SB-5*log10(pixscale)
  return(profitMag2Flux(mag=mag, magzero=magzero))
}

profitGetPixScale=function(header, CD1_1=1, CD1_2=0, CD2_1=0, CD2_2=1){
  if(!missing(header)){
    if(is.data.frame(header) | is.matrix(header)){
      locs=match(c('CD1_1','CD1_2','CD2_1','CD2_2'),header[,1])
      headerWCS=data.frame(header[locs,1],as.numeric(header[locs,2]))
      if('CD1_1' %in% headerWCS[,1]){
        CD1_1=headerWCS[headerWCS[,1]=='CD1_1',2]
        if('CD1_2' %in% headerWCS[,1]){CD1_2=headerWCS[headerWCS[,1]=='CD1_2',2]}else{message('Missing CD1_2')}
      }else{
        if('CDELT1' %in% headerWCS[,1]){
          CD1_1=headerWCS[headerWCS[,1]=='CDELT1',2]
        }else{
          message("Missing CD1_1 and CDELT1")
        }
      }
      if('CD2_2' %in% headerWCS[,1]){
        CD2_2=headerWCS[headerWCS[,1]=='CD2_2',2]
        if('CD2_1' %in% headerWCS[,1]){CD2_1=headerWCS[headerWCS[,1]=='CD2_1',2]}else{message('Missing CD2_1')}
      }else{
        if('CDELT2' %in% headerWCS[,1]){
          CD2_2=headerWCS[headerWCS[,1]=='CDELT2',2]
        }else{
          message("Missing CD2_2 and CDELT2")
        }
      }
    }else{
      if('CD1_1' %in% header){
        CD1_1=as.numeric(header[which(header=='CD1_1')+1])
        if('CD1_2' %in% header){CD1_2=as.numeric(header[which(header=='CD1_2')+1])}else{message('Missing CD1_2')}
      }else{
        if('CDELT1' %in% header){
          CD1_1=header[which(header=='CDELT1')+1,2]
        }else{
          message("Missing CD1_1 and CDELT1")
        }
      }
      if('CD2_2' %in% header){
        CD2_2=as.numeric(header[which(header=='CD2_2')+1])
        if('CD2_1' %in% header){CD2_1=as.numeric(header[which(header=='CD2_1')+1])}else{message('Missing CD2_1')}
      }else{
        if('CDELT1' %in% header){
          CD2_2=header[which(header=='CDELT2')+1,2]
        }else{
          message("Missing CD2_2 and CDELT2")
        }
      }
    }
  }
  return(3600*(sqrt(CD1_1^2+CD1_2^2)+sqrt(CD2_1^2+CD2_2^2))/2)
}

profitInterp2d=function(x,y,image){
    scale=sum(image)
    imagelist=list(x=seq(-dim(image)[1]/2,dim(image)[1]/2,len=dim(image)[1]),y=seq(-dim(image)[2]/2,dim(image)[2]/2,len=dim(image)[2]),z=image)
    ximage = seq(-dim(image)[1]/2,dim(image)[1]/2,len=dim(image)[1])
    yimage = seq(-dim(image)[2]/2,dim(image)[2]/2,len=dim(image)[2])
    zimage = image
    nx = length(ximage)
    ny = length(yimage)
    lx = approx(ximage, 1:nx, x, rule=2)$y
    ly = approx(yimage, 1:ny, y, rule=2)$y
    lx1 = floor(lx)
    ly1 = floor(ly)
    ex = lx - lx1
    ey = ly - ly1
    ex[lx1 == nx] = 1
    ey[ly1 == ny] = 1
    lx1[lx1 == nx] = nx - 1
    ly1[ly1 == ny] = ny - 1
    z=
	zimage[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) +
	zimage[cbind(lx1 + 1, ly1)] * ex * (1 - ey) +
	zimage[cbind(lx1, ly1 + 1)] * (1 - ex) * ey +
	zimage[cbind(lx1 + 1, ly1 + 1)] * ex * ey
  return = cbind(X=x,Y=y,Z=z)
}

profitAddMats=function(matbase, matadd, addloc=c(1,1), plot=FALSE, ...){
  newmat=matrix(0,dim(matbase)[1]+dim(matadd)[1]*2,dim(matbase)[2]+dim(matadd)[2]*2)
  xrangebase=(dim(matadd)[1]+1):(dim(matadd)[1]+dim(matbase)[1])
  yrangebase=(dim(matadd)[2]+1):(dim(matadd)[2]+dim(matbase)[2])
  newmat[xrangebase,yrangebase]=matbase
  xrangeadd=(addloc[1]+dim(matadd)[1]):(addloc[1]+2*dim(matadd)[1]-1)
  yrangeadd=(addloc[2]+dim(matadd)[2]):(addloc[2]+2*dim(matadd)[2]-1)
  if(min(xrangeadd)>=1 & max(xrangeadd)<=dim(newmat)[1] & min(yrangeadd)>=1 & max(yrangeadd)<=dim(newmat)[2]){
    newmat[xrangeadd,yrangeadd]=newmat[xrangeadd,yrangeadd]+matadd
  }
  output=(newmat[xrangebase,yrangebase])
  
  
	if(plot){
	  magimage(output, ...)
  }
  
  return=output
}

profitCheckFinesample <- function(finesample)
{
  stopifnot(is.integer(finesample) && finesample >= 1L)
}

profitParseLikefunc <- function(funcname)
{
  funcname=tolower(funcname)
  if(funcname=="norm" | funcname=="normal")
  {
    return("norm")
  }
  else if(funcname=="chisq" | funcname=="chi-sq")
  {
    return("chisq")
  }
  else if(funcname=="t" | funcname=='student' | funcname=='student-t') {
    return("t")
  }
  else if(funcname=="pois" | funcname=="poisson" | funcname=="cash" | funcname=="c") {
    return("pois")
  }
  else {
    stop(paste0("Error: unknown likelihood function: '",funcname,"'"))
  }
}

profitImBlur=function(image, sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  output=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profitImGrad=function(image, sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  output=as.matrix(imager::enorm(imager::imgradient(imager::isoblur(imager::as.cimg(image),sigma), "xy")))
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profitImDiff=function(image,sigma=1, plot=FALSE, ...){
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for this function to work. Please install it from CRAN.', call. = FALSE)
  }
  blur=as.matrix(imager::isoblur(imager::as.cimg(image),sigma))
  output=image-blur
  if(plot){
    magimage(output, ...)
  }
  return=output
}

profitMakePriors <- function(modellist, sigmas, tolog, means=NULL, tofit=NULL, allowflat=FALSE)
{
  # Sanity checks
  stopifnot(is.logical(allowflat))
  stopifnot(all(is.list(sigmas),is.list(tolog)))
  if(!is.null(means)) stopifnot(is.list(means))
  if(!is.null(tofit)) stopifnot(is.list(tofit))
  
  model = unlist(modellist)
  nparams = length(model)
  stopifnot(all(is.numeric(model) & is.finite(model)))
  
  pformals = list(
    sigmas = unlist(sigmas),
    tolog = unlist(tolog)
  )
  stopifnot(all(pformals$sigmas >0))
  if(!allowflat) stopifnot(all(is.finite(pformals$sigmas)))
  if(!is.null(means)) pformals$means = unlist(means)
  if(!is.null(tofit)) pformals$tofit = unlist(tofit)
  for(formal in names(pformals)) stopifnot(length(pformals[[formal]]) == nparams)
      
  # Define a valid prior function. 
  # tofit will only calculate the prior for fitted values
  # if not otherwise specified, the means will be taken from init
  priors <- function(new, init, sigmas=NULL, tolog=NULL, tofit=NULL, means=unlist(init), allowflat=FALSE)
  {
  	LL = 0
  	parms = unlist(new)
  	if(!is.null(tofit)) ps = which(tofit)
  	else ps = 1:length(parms)
  	for(p in ps)
  	{
  	  if(!(allowflat && (sigmas[p] == Inf)))
  	  {
    		parm = parms[[p]]
    		mean = means[[p]]
    		if(tolog[p])
    		{
    		  parm = log10(parm)
    		  mean = log10(mean)
    		}
    		LL = LL + dnorm(parm,mean,sigmas[p],log=TRUE)
  	  }
  	}
  	return=LL
  }
  for(formal in names(pformals)) formals(priors)[[formal]] = pformals[[formal]]
  formals(priors)$allowflat = allowflat
  stopifnot(is.numeric(priors(modellist,modellist)))
  return=priors
}

profitMakeSigma=function(image, objects=0, sky=0, skyRMS=1, skycut=0, gain=1, readRMS=0, darkRMS=0, image_units='ADU', sky_units='ADU', read_units='ADU', dark_units='ADU', output_units='ADU', plot=FALSE, ...){
  if(!missing(objects) & length(objects)==length(image)){
    image[objects==0]=0
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
  return=sigma
}

profitGainEst=function(image, mask=0, objects=0, sky, skyRMS){
  if(missing(sky)){
    sky=profitSkyEst(image=image, mask=mask, objects=objects,plot=FALSE)$sky
  }
  image_sky=image-sky
  if(missing(skyRMS)){
    skyRMS=profitSkyEst(image=profitImDiff(image_sky,3), mask=mask, objects=objects,plot=FALSE)$skyRMS
  }
  tempval=as.numeric(image_sky[mask==0 & objects==0])
  
  startgain=ceiling(log10(abs(min(tempval, na.rm=T)-sky)/(skyRMS^2)))+1
  
  tempfunc=function(gain,tempval,skyRMS){
    gain=10^gain
    floor=(skyRMS*gain)^2
    trialdata=tempval*gain+floor
    value=-sum(dpois(x=round(trialdata), lambda=floor, log=T))
    return=value
  }

  suppressWarnings({findgain=optim(par=startgain, fn=tempfunc, method="Brent", tempval=tempval, skyRMS=skyRMS, lower=startgain-2, upper=startgain+2)})
  return=10^findgain$par
}

profitDeprojectImageEllipse <- function(image, xcen, ycen, axrat, ang, upsample=5L)
{
  stopifnot(is.integer(upsample))
  if(axrat == 1) return(image)
  stopifnot(axrat > 0 && axrat < 1)
  if(!is.list(image)) image = list(img=image)
  dimorig = dim(image[[1]])
  dimimg = upsample*dimorig
  nimages = length(image)
  for(i in 1:nimages)
  {
    if(!identical(dimorig*upsample,dim(image[[i]])))
    {
      stopifnot(identical(dimorig,dim(image[[i]])))
      image[[i]] = profitUpsample(image[[i]],upsample)
      dim(image[[i]]) = dimimg
    }
  }
  xcen = xcen*upsample
  ycen = ycen*upsample
  
  ang = (ang-90)*pi/180
  maj = c(cos(ang),sin(ang))
  min = c(-maj[2],maj[1])
  x = matrix(rep(0:(dimimg[1] - 1), times=dimimg[2]), nrow=dimimg[1], ncol=dimimg[2])
  y = matrix(rep(0:(dimimg[2] - 1), each=dimimg[1]), nrow=dimimg[1], ncol=dimimg[2])
  idx = 1 + x + dimimg[1]*y
  x = x - 0.5 - xcen
  y = y - 0.5 - xcen
  rmaj = maj[1]*x + maj[2]*y
  rmin = (min[1]*x + min[2]*y)/axrat
  x = ceiling(rmaj*maj[1] + rmin*min[1] + xcen)
  y = ceiling(rmaj*maj[2] + rmin*min[2] + ycen)
  cond = which((x>=1) & (x<=dimimg[1]) & (y>=1) & (y<=dimimg[2]))
  for(j in 1:nimages)
  {
    new = matrix(0,dimimg[1],dimimg[2])
    for(i in cond) new[x[i],y[i]] = new[x[i],y[i]] + image[[j]][idx[i]]
    image[[j]] = profitDownsample(new,upsample)/upsample^2
    dim(image[[j]]) = dimorig
  }
  return(image)
}

profitPoissonMonteCarlo <- function(x)
{
  dimx = dim(x)
  x = rpois(length(x), x)
  dim(x) = dimx
  return(x)
}
