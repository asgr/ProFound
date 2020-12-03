
# .interp.2d.akima=function (x, y, z, xo, yo, ncp = 0, extrap = FALSE, duplicate = "error",
#     dupfun = NULL){
#     if (!(all(is.finite(x)) && all(is.finite(y)) && all(is.finite(z))))
#         stop("missing values and Infs not allowed")
#     if (is.null(xo))
#         stop("xo missing")
#     if (is.null(yo))
#         stop("yo missing")
#     if (ncp > 25) {
#         ncp <- 25
#         cat("ncp too large, using ncp=25\n")
#     }
#     drx <- diff(range(x))
#     dry <- diff(range(y))
#     if (drx == 0 || dry == 0)
#         stop("all data collinear")
#     if (drx/dry > 10000 || drx/dry < 1e-04)
#         stop("scales of x and y are too dissimilar")
#     n <- length(x)
#     np <- length(xo)
#     if (length(yo) != np)
#         stop("length of xo and yo differ")
#     if (length(y) != n || length(z) != n)
#         stop("Lengths of x, y, and z do not match")
#     xy <- paste(x, y, sep = ",")
#     i <- match(xy, xy)
#     if (duplicate == "user" && !is.function(dupfun))
#         stop("duplicate=\"user\" requires dupfun to be set to a function")
#     if (duplicate != "error") {
#         centre <- function(x) {
#             switch(duplicate, mean = mean(x), median = median(x),
#                 user = dupfun(x))
#         }
#         if (duplicate != "strip") {
#             z <- unlist(lapply(split(z, i), centre))
#             ord <- !duplicated(xy)
#             x <- x[ord]
#             y <- y[ord]
#             n <- length(x)
#         }
#         else {
#             ord <- (hist(i, plot = FALSE, freq = TRUE, breaks = seq(0.5,
#                 max(i) + 0.5, 1))$counts == 1)
#             x <- x[ord]
#             y <- y[ord]
#             z <- z[ord]
#             n <- length(x)
#         }
#     }
#     else if (any(duplicated(xy)))
#         stop("duplicate data points")
#     zo <- rep(0, np)
#     storage.mode(zo) <- "double"
#     miss <- !extrap
#     misso <- seq(miss, np)
#     if (extrap & ncp == 0)
#         warning("Cannot extrapolate with linear option")
#     ans <- .Fortran("idbvip", as.integer(1), as.integer(ncp),
#         as.integer(n), as.double(x), as.double(y), as.double(z),
#         as.integer(np), x = as.double(xo), y = as.double(yo),
#         z = zo, integer((31 + ncp) * n + np), double(8 * n),
#         misso = as.logical(misso), PACKAGE = "akima")
#     temp <- ans[c("x", "y", "z", "misso")]
#     temp$z[temp$misso] <- NA
#     temp[c("x", "y", "z")]
# }

profoundSkyEst=function(image=NULL, objects=NULL, mask=NULL, cutlo=cuthi/2, 
                        cuthi=sqrt(sum((dim(image)/2)^2)), skycut='auto', 
                        clipiters=5, radweight=0, plot=FALSE, ...){
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
  return(invisible(list(sky=sky,skyerr=skyerr,skyRMS=skyRMS,Nnearsky=Nnearsky,radrun=tempmedian)))
}

profoundSkyEstLoc=function(image=NULL, objects=NULL, mask=NULL, loc=dim(image)/2, 
                           box=c(100,100), skytype='median', skyRMStype='quanlo', 
                           sigmasel=1, skypixmin=prod(box)/2, boxadd=box/2, boxiters=0, 
                           conviters = 100, doclip=TRUE, shiftloc = FALSE, paddim = TRUE, plot=FALSE, ...){
  if(!is.null(objects) | !is.null(mask)){
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
      return(invisible(c(NA, NA)))
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
  }else if(skytype=='converge'){
    stats = .converge_sky(clip, conviters = conviters)
    skyloc = stats[1]
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
    skyRMSloc=sqrt(.varwt(clip, xcen=skyloc))
  }else if(skytype=='converge' & skyRMStype=='converge'){
    skyRMSloc = stats[2]
  }else if(skytype!='converge' & skyRMStype=='converge'){
    skyRMSloc = .converge_sky(clip, conviters = conviters)[2]
  }
  
  return(invisible(c(skyloc, skyRMSloc)))
}

profoundMakeSkyMap=function(image=NULL, objects=NULL, mask=NULL, sky=0, box=c(100,100), 
                            grid=box, skytype='median', skyRMStype='quanlo', sigmasel=1, 
                            skypixmin=prod(box)/2, boxadd=box/2, boxiters=0,
                            conviters = 100, doclip=TRUE, shiftloc = FALSE, paddim = TRUE, cores=1){
  xseq=seq(grid[1]/2,dim(image)[1],by=grid[1])
  yseq=seq(grid[2]/2,dim(image)[2],by=grid[2])
  tempgrid=expand.grid(xseq, yseq)

  image = image - sky
  
  if(cores>1 & 'foreach' %in% .packages() & 'snow' %in% .packages() & 'doSNOW' %in% .packages() & 'bigmemory'  %in% .packages()){
    cl = snow::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    
    image_big = bigmemory::as.big.matrix(image)
    image_desc = bigmemory::describe(image_big)
    i=NULL
    
    tempsky=foreach(i = 1:dim(tempgrid)[1], .combine='rbind')%dopar%{
      image_loop = bigmemory::attach.big.matrix(image_desc)
      profoundSkyEstLoc(image=image_loop, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]),
                        box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel,
                        skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, conviters=conviters,
                        doclip=doclip, shiftloc=shiftloc, paddim=paddim)
    }
    closeAllConnections()
    tempsky=rbind(tempsky)
  }else{
    if(cores>1){
      message('Missing parallel backend packages (need foreach, snow, doSNOW and bigmemory)')
    }
    tempsky=matrix(0,dim(tempgrid)[1],2)
    for(i in 1:dim(tempgrid)[1]){
      tempsky[i,]=profoundSkyEstLoc(image=image, objects=objects, mask=mask,
                                    loc=as.numeric(tempgrid[i,]), box=box, skytype=skytype,
                                    skyRMStype=skyRMStype, sigmasel=sigmasel, skypixmin=skypixmin,
                                    boxadd=boxadd, boxiters=boxiters, conviters=conviters,
                                    doclip=doclip, shiftloc=shiftloc, paddim=paddim)
    }
  }
  
  tempmat_sky=matrix(tempsky[,1],length(xseq))
  tempmat_skyRMS=matrix(tempsky[,2],length(xseq))
  #tempmat_sky[is.na(tempmat_sky)]=median(tempmat_sky, na.rm = TRUE)
  #tempmat_skyRMS[is.na(tempmat_skyRMS)]=median(tempmat_skyRMS, na.rm = TRUE)
  return(invisible(list(sky=list(x=xseq, y=yseq, z=tempmat_sky+sky), skyRMS=list(x=xseq, y=yseq, z=tempmat_skyRMS))))
}

profoundMakeSkyGrid=function(image=NULL, objects=NULL, mask=NULL, sky=0, box=c(100,100),
                             grid=box, skygrid_type='new', type='bicubic', 
                             skytype='median', skyRMStype='quanlo', sigmasel=1,
                             skypixmin=prod(box)/2, boxadd=box/2, boxiters=0,
                             conviters = 100, doclip=TRUE, shiftloc = FALSE, paddim = TRUE, cores=1){
  
  if(length(box)==1){
    box=rep(box,2)
    if(missing(grid)){grid=box}
    if(missing(boxadd)){boxadd=box/2}
    if(missing(skypixmin)){skypixmin=prod(box)/2}
  }
  if(length(grid)==1){
    grid=rep(grid,2)
  }
  if(length(boxadd)==1){
    boxadd=rep(boxadd,2)
  }
  
  if(box[1] > ceiling(dim(image)[1]/3)){
    box[1] = ceiling(dim(image)[1]/3)
    message('dim(image)[1]/box[1] must be >=3, box[1] modified to ',box[1])
  }
  if(box[2] > ceiling(dim(image)[1]/3)){
    box[2] = ceiling(dim(image)[2]/3)
    message('dim(image)[2]/box[2] must be >=3, box[2] modified to ',box[2])
  }
  
  if(grid[1] > ceiling(dim(image)[1]/3)){
    grid[1] = ceiling(dim(image)[1]/3)
    message('dim(image)[1]/grid[1] must be >=3, grid[1] modified to ',grid[1])
  }
  if(grid[2] > ceiling(dim(image)[1]/3)){
    grid[2] = ceiling(dim(image)[2]/3)
    message('dim(image)[2]/grid[2] must be >=3, grid[2] modified to ',grid[2])
  }
  
  image = image - sky
  
  if(skygrid_type=='new'){
    # void .Cadacs_MakeSkyGrid(Rcpp::NumericMatrix image,
    #                         Rcpp::NumericMatrix sky, Rcpp::NumericMatrix skyRMS,
    #                         Rcpp::Function Fquantile,
    #                         Rcpp::Nullable<Rcpp::IntegerMatrix> objects = R_NilValue,
    #                         Rcpp::Nullable<Rcpp::IntegerMatrix> mask = R_NilValue,
    #                         const int box1 = 100, const int box2 = 100,
    #                         const int grid1 =100, const int grid2 = 100,
    #                         const int boxadd1 = 50, const int boxadd2 = 50,
    #                         const int type = 2, const int skypixmin = 5000,
    #                         const int boxiters = 0, const int doclip = 1,
    #                         const int skytype = 1, const int skyRMStype = 2,
    #                         const double sigmasel = 1
    # )
    #define adacs_BOTH 1
    #define adacs_LO 2
    #define adacs_HI 3
    #define adacs_SD 4
    #define adacs_RBOTH 5
    #define adacs_RLO 6
    #define adacs_RHI 7
    #define adacs_RSD 8
    
    #define adacs_AUTO 1
    #define adacs_SET 2
    
    #define adacs_MEDIAN 1
    #define adacs_MEAN 2
    #define adacs_MODE 3
    #define adacs_RMEDIAN 4
    #define adacs_RMEAN 5
    #define adacs_RMODE 6
    
    type = switch(type, bilinear=1L, bilinear_new=1L, bicubic=2L, bicubic_new=2L)
    skytype = switch(skytype, median=1L, mean=2L, mode=3L)
    skyRMStype = switch(skyRMStype, quanboth=1L, quanlo=2L, quanhi=3L, sd=4L)
    
    temp_bi_sky = matrix(0, nrow=dim(image)[1], ncol=dim(image)[2])
    temp_bi_skyRMS = matrix(0, nrow=dim(image)[1], ncol=dim(image)[2])
    .Cadacs_MakeSkyGrid(image = image,
                        sky = temp_bi_sky,
                        skyRMS = temp_bi_skyRMS,
                        objects = objects,
                        mask = mask,
                        box1 = box[1], box2 = box[2],
                        grid1 = grid[1], grid2 = grid[2],
                        boxadd1 = boxadd[1], boxadd2 = boxadd[2],
                        type = type,
                        skypixmin = skypixmin,
                        boxiters = boxiters,
                        doclip = doclip,
                        skytype = skytype,
                        skyRMStype = skyRMStype,
                        sigmasel = sigmasel
                        )
  }else if(skygrid_type=='old'){
  
    # if(!requireNamespace("akima", quietly = TRUE)){
    #   if(type=='bicubic'){
    #     stop('The akima package is needed for bicubic interpolation to work. Please install it from CRAN.', call. = FALSE)
    #   }
    #   if(type=='bilinear'){
    #     useakima=FALSE
    #   }
    # }else{
    #   useakima=TRUE
    # }
    
    if(box[1]>dim(image)[1]){box[1]=dim(image)[1]}
    if(box[2]>dim(image)[2]){box[2]=dim(image)[2]}
    if(grid[1]>dim(image)[1]){grid[1]=dim(image)[1]}
    if(grid[2]>dim(image)[2]){grid[2]=dim(image)[2]}
    
    xseq = seq(grid[1]/2, ceiling(dim(image)[1]/grid[1])*grid[1], by=grid[1])
    yseq = seq(grid[2]/2, ceiling(dim(image)[2]/grid[2])*grid[2], by=grid[2])
    tempgrid = expand.grid(xseq, yseq)
    
    if(cores>1 & 'foreach' %in% .packages() & 'snow' %in% .packages() & 'doSNOW' %in% .packages() & 'bigmemory'  %in% .packages()){
      cl = snow::makeCluster(cores)
      doSNOW::registerDoSNOW(cl)
      
      image_big = bigmemory::as.big.matrix(image)
      image_desc = bigmemory::describe(image_big)
      i=NULL
      
      tempsky=foreach(i = 1:dim(tempgrid)[1], .combine='rbind')%dopar%{
        image_loop = bigmemory::attach.big.matrix(image_desc)
        profoundSkyEstLoc(image=image_loop, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]),
                          box=box, skytype=skytype, skyRMStype=skyRMStype, sigmasel=sigmasel,
                          skypixmin=skypixmin, boxadd=boxadd, boxiters=boxiters, conviters=conviters,
                          doclip=doclip, shiftloc=shiftloc, paddim=paddim)
      }
      closeAllConnections()
      tempsky=rbind(tempsky)
    }else{
      if(cores>1){
        message('Missing parallel backend packages (need foreach, snow, doSNOW and bigmemory)')
      }
      tempsky=matrix(0,dim(tempgrid)[1],2)
      for(i in 1:dim(tempgrid)[1]){
        tempsky[i,]=profoundSkyEstLoc(image=image, objects=objects, mask=mask, loc=as.numeric(tempgrid[i,]),
                                      box=box, skytype=skytype, skyRMStype=skyRMStype,
                                      sigmasel=sigmasel, skypixmin=skypixmin, boxadd=boxadd,
                                      boxiters=boxiters, conviters=conviters, doclip=doclip, shiftloc=shiftloc, paddim=paddim)
      }
    }
    
    xseq=c(-grid[1]/2,xseq,max(xseq)+grid[1]/2)
    yseq=c(-grid[2]/2,yseq,max(yseq)+grid[2]/2)
    
    tempmat_sky=matrix(0,length(xseq),length(yseq))
    tempmat_sky[2:(length(xseq)-1),2:(length(yseq)-1)]=tempsky[,1]
    tempmat_sky[is.na(tempmat_sky)]= stats::median(tempsky[,1], na.rm = TRUE)
    
    tempmat_skyRMS=matrix(0,length(xseq),length(yseq))
    tempmat_skyRMS[2:(length(xseq)-1),2:(length(yseq)-1)]=tempsky[,2]
    tempmat_skyRMS[is.na(tempmat_skyRMS)]=stats::median(tempsky[,2], na.rm = TRUE)
    
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
      
      # if(type=='bilinear'){
      #   bigridx=rep(1:dim(image)[1]-0.5,times=dim(image)[2])
      #   bigridy=rep(1:dim(image)[2]-0.5,each=dim(image)[1])
      #   if(useakima){
      #     tempgrid=expand.grid(xseq, yseq)
      #     temp_bi_sky=.interp.2d.akima(x=tempgrid[,1], y=tempgrid[,2], z=as.numeric(tempmat_sky),xo=bigridx, yo=bigridy)$z
      #     temp_bi_skyRMS=.interp.2d.akima(x=tempgrid[,1], y=tempgrid[,2], z=as.numeric(tempmat_skyRMS),xo=bigridx, yo=bigridy)$z
      #   }else{
      #     temp_bi_sky=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_sky))
      #     temp_bi_skyRMS=.interp.2d(bigridx, bigridy, list(x=xseq, y=yseq, z=tempmat_skyRMS))
      #   }
      # }else if(type=='bilinear_new'){
      #   temp_bi_sky=matrix(0, dim(image)[1], dim(image)[2])
      #   .interpolateLinearGrid(xseq, yseq, tempmat_sky, temp_bi_sky)
      #   temp_bi_skyRMS=matrix(0, dim(image)[1], dim(image)[2])
      #   .interpolateLinearGrid(xseq, yseq, tempmat_skyRMS, temp_bi_skyRMS)
      # }else if(type=='bicubic'){
      #   bigridx=rep(1:dim(image)[1]-0.5,times=dim(image)[2])
      #   bigridy=rep(1:dim(image)[2]-0.5,each=dim(image)[1])
      #   temp_bi_sky=akima::bicubic(xseq, yseq, tempmat_sky, bigridx, bigridy)$z
      #   temp_bi_skyRMS=akima::bicubic(xseq, yseq, tempmat_skyRMS, bigridx, bigridy)$z
      # }else if(type=='bicubic_new'){
      #   temp_bi_sky=matrix(0, dim(image)[1], dim(image)[2])
      #   .interpolateAkimaGrid(xseq, yseq, tempmat_sky, temp_bi_sky)
      #   temp_bi_skyRMS=matrix(0, dim(image)[1], dim(image)[2])
      #   .interpolateAkimaGrid(xseq, yseq, tempmat_skyRMS, temp_bi_skyRMS)
      # }else{
      #   stop('type must be one of bilinear / bicubic !')
      # }
      
      if(type=='bilinear'){
        temp_bi_sky=matrix(0, dim(image)[1], dim(image)[2])
        .interpolateLinearGrid(xseq, yseq, tempmat_sky, temp_bi_sky)
        temp_bi_skyRMS=matrix(0, dim(image)[1], dim(image)[2])
        .interpolateLinearGrid(xseq, yseq, tempmat_skyRMS, temp_bi_skyRMS)
      }else if(type=='bicubic'){
        temp_bi_sky=matrix(0, dim(image)[1], dim(image)[2])
        .interpolateAkimaGrid(xseq, yseq, tempmat_sky, temp_bi_sky)
        temp_bi_skyRMS=matrix(0, dim(image)[1], dim(image)[2])
        .interpolateAkimaGrid(xseq, yseq, tempmat_skyRMS, temp_bi_skyRMS)
      }else{
        stop('type must be one of bilinear / bicubic !')
      }
    
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
  }else{
    stop('skygrid_type must be new/old!')
  }
  
  return(invisible(list(sky=temp_bi_sky+sky, skyRMS=temp_bi_skyRMS)))
}

profoundMakeSkyBlur=function(image=NULL, objects=NULL, box=100, sigma=mean(box)*(4/pi)/sqrt(12)){
  if(!is.null(objects)){
    image[objects==1]=NA
  }
  if(requireNamespace("imager", quietly = TRUE)){
    return(invisible(as.matrix(imager::isoblur(imager::as.cimg(image), sigma=sigma, na.rm=TRUE))))
  }else{
    stop('The imager package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
  }
}

profoundChisel=function(image=NULL, sky=NULL, skythresh=0.005, blurcut=0.01, objthresh=1-skythresh,
                        sigma=2, iterchisel=5, itersky=5, objbias=0.5, box=100, skyconv=1e-2){
  if(is.null(sky)){
    sky = matrix(0, dim(image)[1], dim(image)[2])
  }
  
  origimage = image
  objcomp = matrix(1, dim(image)[1], dim(image)[2])
  
  for(i in 1:itersky){
    image = origimage - sky
    
    objcut = quantile(image, objthresh, na.rm = TRUE)
    objloc = image > objcut
    objmask = profoundImBlur(objloc, sigma=sigma) > blurcut
    objects = matrix(1, dim(image)[1], dim(image)[2])
    
    skyloc = matrix(0, dim(image)[1], dim(image)[2])
    skymask = matrix(1, dim(image)[1], dim(image)[2])
    
    for(j in 1:iterchisel){
      skycut = quantile(image[objects==1], skythresh, na.rm = TRUE)
      skyloc = skymask==0 | (image < skycut)
      skymask = profoundImBlur(skyloc, sigma=sigma) < blurcut
      if(j > 1){
        objcut = quantile(image[objmask==0], objthresh, na.rm = TRUE)
        if(abs(skycut)*(1-objbias) > abs(objcut)*objbias){break}
        objloc = objmask | (image > objcut)
        objmask = profoundImBlur(objloc, sigma=sigma) > blurcut
      }
      objects = skymask | objmask
    }
    
    sky = profoundMakeSkyBlur(image=image, objects=objects, box=box) + sky
    
    if(i > 1){
      if(sum(abs(objcomp - objects)) / prod(dim(image)) < 1e-2){break}
    }
    
    objcomp = objects
  }
  return(list(objects=objects, sky=sky))
}

.tnorm_upper = function(mu,sd,upper){
  ##return the expectation of a truncated normal distribution
  upper.std=(upper-mu)/sd
  ratio = dnorm(upper.std) / pnorm(upper.std)
  mean = mu + -(sd * ratio)
  sd = sd * sqrt(1 + -upper.std*ratio - ratio^2)
  return(c(mean,sd))
}

.converge_sky = function(data, conviters = 100){
  data = as.numeric(data)
  if('Rfast' %in% .packages()){
    med_data=try(Rfast::med(data, na.rm=TRUE), silent=TRUE)
    if(class(med_data)=='try-error'){med_data=NA}
  }else{
    med_data = stats::median(data, na.rm=TRUE)
  }
  rough_sd = 2*sd(data[data < med_data])
  data = data[data < med_data + 5*rough_sd]
  pop_mu = mean(data)
  pop_sd = sd(data)
  uppercut = pop_mu + pop_sd
  data = data[data < pop_mu + pop_sd]
  samp_mean = mean(data)
  samp_sd = sd(data)
  
  for(i in 1:conviters){
    est_stats = .tnorm_upper(pop_mu, pop_sd, uppercut)
    pop_mu = pop_mu + (samp_mean - est_stats[1])*2
    pop_sd = pop_sd * samp_sd/est_stats[2]
  }
  #print(c(samp_mean, samp_sd))
  #print(est_stats)
  return(c(pop_mu, pop_sd))
}

