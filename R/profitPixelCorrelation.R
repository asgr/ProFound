profitPixelCorrelation=function(image, objects, mask, sky=0, skyRMS=1, lag=c(1:9,1:9*10,1:9*100,1:9*1e3,1:9*1e4), fft=TRUE, plot=FALSE){
  
  xlen=dim(image)[1]
  ylen=dim(image)[2]
  
  lag=lag[lag<xlen & lag<ylen]
  
  image=(image-sky)/skyRMS
  
  if(!missing(objects)){
    image[objects>0]=NA
  }
  
  if(!missing(mask)){
    image[mask>0]=NA
  }
  
  corx={}; cory={}; relsdx={}; relsdy={}
  
  for(i in lag){
    corx=c(corx,cor(as.numeric(image[1:(xlen-i),]), as.numeric(image[1:(xlen-i)+i,]), use="complete.obs"))
    cory=c(cory,cor(as.numeric(image[,1:(ylen-i)]), as.numeric(image[,1:(ylen-i)+i]), use="complete.obs"))
    relsdx=c(relsdx,sd(as.numeric(image[1:(xlen-i),])-as.numeric(image[1:(xlen-i)+i,]), na.rm=TRUE)/sqrt(2))
    relsdy=c(relsdy,sd(as.numeric(image[,1:(ylen-i)])-as.numeric(image[,1:(ylen-i)+i]), na.rm=TRUE)/sqrt(2))
  }
  
  output_cortab=data.frame(lag=lag,corx=corx,cory=cory,relsdx=relsdx,relsdy=relsdy)
  
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
    output_FFT=Mod(fft(z=image*centre))
  }else{
    output_FFT=NULL
    image=NULL
  }
  
  if(plot){
    magplot(output_cortab[,c('lag','corx')], xlim=c(1,max(lag)), ylim=c(-1,1), xlab='Pixel Lag', ylab='Correlation', type='l', col='blue', log='x', grid=TRUE)
    lines(output_cortab[,c('lag','cory')], col='red')
    lines(output_cortab[,c('lag','relsdx')], col='blue', lty=2)
    lines(output_cortab[,c('lag','relsdy')], col='red', lty=2)
    legend('bottomright', legend=c('x-cor','y-cor','x-rel-sd','y-rel-sd'), col=c('blue','red'), lty=c(1,1,2,2))
  }
  
  return=list(cortab=output_cortab, fft=output_FFT, image_sky=image)
}