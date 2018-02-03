profoundPixelCorrelation=function(image, objects, mask, sky=0, skyRMS=1, lag=c(1:9,1:9*10,1:9*100,1:9*1e3,1:9*1e4), fft=TRUE, plot=FALSE, ylim=c(-1,1), log='x', grid=TRUE, ...){
  
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  
  lag=lag[lag<xlen & lag<ylen]
  
  image=(image-sky)/skyRMS
  
  flag_neg=image<0
  flag_pos=image>0
  
  if(!missing(objects)){
    image[objects>0]=NA
  }
  
  if(!missing(mask)){
    image[mask>0]=NA
  }
  
  image_neg=image
  image_neg[flag_neg==FALSE]=NA #image is left with negative pixels only
  image_pos=image
  image_pos[flag_pos==FALSE]=NA #image is left with positive pixels only
  
  
  corx={}; cory={}; corx_neg={}; cory_neg={}; corx_pos={}; cory_pos={}; relsdx={}; relsdy={}
  
  for(i in lag){
    corx=c(corx,cor(as.numeric(image[1:(xlen-i),]), as.numeric(image[1:(xlen-i)+i,]), use="na.or.complete"))
    cory=c(cory,cor(as.numeric(image[,1:(ylen-i)]), as.numeric(image[,1:(ylen-i)+i]), use="na.or.complete"))
    corx_neg=c(corx_neg,cor(as.numeric(image_neg[1:(xlen-i),]), as.numeric(image[1:(xlen-i)+i,]), use="na.or.complete"))
    cory_neg=c(cory_neg,cor(as.numeric(image_neg[,1:(ylen-i)]), as.numeric(image[,1:(ylen-i)+i]), use="na.or.complete"))
    corx_pos=c(corx_pos,cor(as.numeric(image_pos[1:(xlen-i),]), as.numeric(image[1:(xlen-i)+i,]), use="na.or.complete"))
    cory_pos=c(cory_pos,cor(as.numeric(image_pos[,1:(ylen-i)]), as.numeric(image[,1:(ylen-i)+i]), use="na.or.complete"))
    relsdx=c(relsdx,sd(as.numeric(image[1:(xlen-i),])-as.numeric(image[1:(xlen-i)+i,]), na.rm=TRUE)/sqrt(2))
    relsdy=c(relsdy,sd(as.numeric(image[,1:(ylen-i)])-as.numeric(image[,1:(ylen-i)+i]), na.rm=TRUE)/sqrt(2))
  }
  
  output_cortab=data.frame(lag=lag, corx=corx, cory=cory, corx_neg=corx_neg, cory_neg=cory_neg, corx_pos=corx_pos, cory_pos=cory_pos, corx_diff=corx_pos-corx_neg, cory_diff=cory_pos-cory_neg, relsdx=relsdx, relsdy=relsdy)
  
  if(fft){
    xlenpad=xlen+xlen%%2
    ylenpad=ylen+ylen%%2
    
    if(!missing(objects)){
      image[objects>0]=rnorm(length(which(objects>0)))
    }
    
    if(!missing(mask)){
      image[mask>0]=rnorm(length(which(mask>0)))
    }
    
    centre=matrix(c(rep(c(-1,1),xlenpad/2), rep(c(1,-1),ylenpad/2)),xlenpad,ylenpad)[1:xlen,1:ylen]
    FFT_Mod=Mod(fft(z=image*centre))
    dimx=dim(FFT_Mod)[1]
    dimy=dim(FFT_Mod)[2]
    if(dimx %% 2 == 0){
      x=c(-(dimx/2+0.5):(dimx/2-0.5))
    }else{
      x=c(-floor(dimx/2):floor(dimx/2))
    }
    if(dimy %% 2 == 0){
      y=c(-(dimy/2+0.5):(dimy/2-0.5))
    }else{
      y=c(-floor(dimy/2):floor(dimy/2))
    }
    output_FFT=list(x=x, y=y, z=FFT_Mod)
  }else{
    output_FFT=NULL
    image=NULL
  }
  
  if(plot){
    magplot(output_cortab[,c('lag','corx')], xlim=c(1,max(lag)), xlab='Pixel Lag', ylab='Correlation', type='l', col='blue', ylim=ylim, log=log, grid=grid, ...)
    lines(output_cortab[,c('lag','cory')], col='red')
    lines(output_cortab[,c('lag','corx_diff')], col='blue', lty=2)
    lines(output_cortab[,c('lag','cory_diff')], col='red', lty=2)
    lines(output_cortab[,c('lag','relsdx')], col='blue', lty=3)
    lines(output_cortab[,c('lag','relsdy')], col='red', lty=3)
    legend('topright', legend=c('x-cor','y-cor','x-cor-diff','y-cor-diff','x-rel-sd','y-rel-sd'), col=c('blue','red'), lty=c(1,1,2,2,3,3), bg='white')
  }
  
  return=list(cortab=output_cortab, fft=output_FFT, image_sky=image)
}

profoundSkySplitFFT=function(image, objects, mask, sky=0, skyRMS=1, skyscale=100, profound){
  if(!missing(image)){
    if(class(image)=='profound'){
      if(missing(objects)){objects=image$objects_redo}
      if(missing(mask)){mask=image$mask}
      if(missing(sky)){sky=image$sky}
      if(missing(skyRMS)){skyRMS=image$skyRMS}
      image=image$image
      if(is.null(image)){stop('Need image in profound object to be non-Null')}
    }
  }
  if(!missing(profound)){
    if(class(profound) != 'profound'){
      stop('Class of profound input must be of type \'profound\'')
    }
    if(missing(image)){image=profound$image}
    if(missing(objects)){objects=profound$objects_redo}
    if(missing(mask)){mask=profound$mask}
    if(missing(sky)){sky=profound$sky}
    if(missing(skyRMS)){skyRMS=profound$skyRMS}
  }
  
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  
  xlenpad=xlen+xlen%%2
  ylenpad=ylen+ylen%%2
  
  hassky=!missing(sky)
  hasskyRMS=!missing(skyRMS)
  
  image=image-sky
  
  if(!missing(objects)){
    if(is.null(objects)==FALSE){
      if(!hassky | !hasskyRMS){stop('Need sky and skyRMS for object padding')}
      sel_objects=objects>0
      image[sel_objects]=rnorm(length(which(sel_objects)),mean=0,sd=skyRMS[sel_objects])
    }
  }
    
  if(!missing(mask)){
    if(is.null(mask)==FALSE){
      if(!hassky | !hasskyRMS){stop('Need sky and skyRMS for mask padding')}
      sel_mask=mask>0
      image[sel_mask]=rnorm(length(which(sel_mask)),mean=0,sd=skyRMS[sel_objects])
    }
  }

  fft_orig=fft(image)
  clipfreqx=ceiling(xlen/skyscale)
  clipfreqy=ceiling(ylen/skyscale)
  
  fft_orig[1:clipfreqx, 1:clipfreqy]=0
  fft_orig[xlen+1-1:clipfreqx, 1:clipfreqy]=0
  fft_orig[1:clipfreqx, ylen+1-1:clipfreqy]=0
  fft_orig[xlen+1-1:clipfreqx, ylen+1-1:clipfreqy]=0
  
  # if(!missing(modclip)){
  #   FFT_Mod=Mod(fft_orig)
  #   fft_orig[FFT_Mod>quantile(FFT_Mod,modclip)]=0
  # }
  
  image_new=Re(fft(fft_orig,inverse=TRUE))/prod(xlen,ylen)
  return=list(sky=sky+image-image_new, sky_lo=image-image_new, sky_hi=image_new)
}
