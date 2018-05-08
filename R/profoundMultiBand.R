profoundMultiBand=function(inputlist=NULL, dir='', mask, skycut = 1, pixcut = 3, tolerance = 4, ext = 2, sigma = 1, smooth = TRUE, iters_tot=2, detectbands='r', multibands=c('u','g','r','i','z'), magzero=0, gain=NULL, catappend=multibands, totappend='t', colappend='c', dotot=TRUE, docol=TRUE, boundstats=TRUE, haralickstats=TRUE, verbose=FALSE){
  
  # v1.0 of the multiband script
  # Written and maintained by Aaron Robotham (inspired by scripts by Soheil Koushan and Simon Driver)
  # This should only be sourced, for *any* modifications contact Aaron Robotham
  
  # The most important thing is that all of the input images must be pixel matched via SWarp or magwarp etc
  # detectbands and multibands are the names of the target bands, which should be the names of the images ignoring the .fits ending
  # By default Simon's script names them as u.fits and g.fits etc, hence the literal bands naming
  # detectbands can be a vector of target detection bands for stacking
  # multibands is by default a vector of all optical and NIR bands
  # You can select whether to do total and/or colour photometry with dotot and docol
  # magzero must be matched to the number of bands in multibands (which reflects the magzero in each respectively), or one value which is then recycled
  # catappend are the band names to use as the extra tag attached to the column names in the cat_tot and cat_col outputs (by default this is linked to multibands)
  # totappend and colapped is just the extra tag attached to the column names in the cat_tot and cat_col outputs
  # boundstats and haralickstats only apply to the detection photometry
  
  timestart=proc.time()[3]
  
  # Restrict outselves to data actually present (no matter what is asked for)
  
  if(is.null(inputlist)){
    presentbands=unlist(strsplit(list.files(dir), split='.fits'))
  }else{
    if(length(inputlist)!=length(multibands)){
      stop('inputlist is not the same length as multibands!')
    }
    presentbands=multibands
  }

  if(detectbands[1]=='get'){
    detectbands=presentbands
  }else{
    detectbands=detectbands[detectbands %in% presentbands]
  }
  
  if(multibands[1]=='get'){
    multibands=presentbands
  }
  
  if(catappend[1]=='get'){
    catappend=presentbands
  }
  
  if(length(magzero)==1){
    magzero=rep(magzero, length(multibands))
  }
  
  if(length(gain)==1){
    gain=rep(gain, length(multibands))
  }
  
  magzero=magzero[which(multibands %in% presentbands)]
  magzero=magzero[!is.na(magzero)]
  catappend=catappend[which(multibands %in% presentbands)]
  catappend=catappend[!is.na(catappend)]
  if(!is.null(gain)){
    gain=gain[which(multibands %in% presentbands)]
    gain=gain[!is.na(gain)]
  }
  multibands=multibands[which(multibands %in% presentbands)]
  multibands=multibands[!is.na(multibands)]
  
  # Some safety checks
  
  if(!is.null(gain)){
    if(length(unique(c(length(multibands), length(magzero), length(catappend), length(gain))))!=1){
      stop('multibands, magzero, catappend and gain are not all the same length - they must be!')
    }
  }else{
    if(length(unique(c(length(multibands), length(magzero), length(catappend))))!=1){
      stop('multibands, magzero and catappend are not all the same length - they must be!')
    }
  }
  if(!is.null(inputlist)){
    if(length(unique(c(length(multibands),length(inputlist))))!=1){
      stop('inputlist does not match the length of multibands - it must be!')
    }
  }
  
  message(paste('*** Will use',paste(detectbands,collapse=''),'for source detection ***'))
  
  if(dotot | docol){
    message(paste('*** Will use',paste(multibands,collapse=''),'for multi band photometry ***'))
    message(paste('*** Magzero:',paste(multibands,magzero,sep='=', collapse=' '),' ***'))
    if(!is.null(gain)){
      message(paste('*** Gain:',paste(multibands,gain,sep='=', collapse=' '),' ***'))
    }
    
    if(dotot){
      message('*** Will compute total multi band photometry ***')
    }
    if(docol){
      message('*** Will compute isophotal colour multi band photometry ***')
    }
  }
  
  if(length(detectbands)==1){
    
    # If only one band is specified for detection then we skip the stacking part
    
    message(paste('*** Currently processing single detection band',detectbands,'***'))
    if(is.null(inputlist)){
      detect=readFITS(paste0(dir,detectbands,'.fits'))
    }else{
      detect=inputlist[[which(multibands==detectbands)]]
    }
    temp_magzero=magzero[multibands==detectbands]
    pro_detect=profoundProFound(detect, mask=mask, skycut=skycut, pixcut=pixcut, tolerance=tolerance, ext=ext, sigma=sigma, smooth=smooth, magzero=temp_magzero, verbose=verbose, boundstats=boundstats, haralickstats=haralickstats)
  }else{
  
    # Multiple detection bands requested, so we prepare lists for stacking
    
    detect_image=list()
    detect_sky=list()
    detect_skyRMS=list()
    detect_magzero={}
    
    for(i in 1:length(detectbands)){
      
      # Loop around detection bands
      
      message(paste('*** Currently processing detection band',detectbands[i],'***'))
      
      if(is.null(inputlist)){
        detect=readFITS(paste0(dir,detectbands[i],'.fits'))
      }else{
        detect=inputlist[[which(multibands==detectbands[i])]]
      }
      temp_magzero=magzero[multibands==detectbands[i]]
      
      # Run ProFound on current detection band with input parameters
      
      pro_detect=profoundProFound(detect, mask=mask, skycut=skycut, pixcut=pixcut, tolerance=tolerance, ext=ext, sigma=sigma, smooth=smooth, magzero=temp_magzero, verbose=verbose)
      
      # Append to lists for stacking
      
      detect_image=c(detect_image,list(pro_detect$image))
      detect_sky=c(detect_sky,list(pro_detect$sky))
      detect_skyRMS=c(detect_skyRMS,list(pro_detect$skyRMS))
      detect_magzero=c(detect_magzero, temp_magzero)
    }
    
    # Grab the current header
    
    header=pro_detect$header
    
    # Delete and clean up
    
    rm(detect)
    rm(pro_detect)
    gc()
    
    # Stack!!!
    # We first stack the image then the sky.
    
    detect_image_stack=profoundMakeStack(image_list=detect_image, sky_list=detect_sky, skyRMS_list=detect_skyRMS, magzero_in=detect_magzero, magzero_out=detect_magzero[1])
    detect_sky_stack=profoundMakeStack(image_list=detect_sky, skyRMS_list=detect_skyRMS, magzero_in=detect_magzero, magzero_out=detect_magzero[1])
    
    # Delete and clean up
    
    rm(detect_image)
    rm(detect_skyRMS)
    gc()
    
    message('*** Currently processing stacked detection image ***')
    
    # For reference we run ProFound with the stacked sky added back in, passing it the stacked sky too.
    
    pro_detect=profoundProFound(image=detect_image_stack$image+detect_sky_stack$image, mask=mask, header=header, skycut=skycut, pixcut=pixcut, tolerance=tolerance, ext=ext, sigma=sigma,  smooth=smooth, magzero=detect_magzero[1], sky=detect_sky_stack$image, skyRMS=detect_image_stack$skyRMS, redosky=FALSE, verbose=verbose, boundstats=boundstats, haralickstats=haralickstats)
    
    # Delete and clean up
    
    rm(detect_image_stack)
    rm(detect_sky_stack)
    gc()
  }
  
  # Create base total and colour photometry catalogues
  
  if(dotot){
    cat_tot=data.frame(segID=1:pro_detect$Nseg)
  }else{
    cat_tot=NULL
  }
  if(docol){
    cat_col=data.frame(segID=1:pro_detect$Nseg)
  }else{
    cat_col=NULL
  }

  
  if(dotot | docol){
    for(i in 1:length(multibands)){
      
      # Loop around multi bands
      
      message(paste('*** Currently processing multi band',multibands[i],'***'))
      
      if(is.null(inputlist)){
        multi=readFITS(paste0(dir,multibands[i],'.fits'))
      }else{
        multi=inputlist[[i]]
      }
      
      if(dotot){
        
        # Compute total multi band photometry, allowing some extra dilation via the iters_tot argument
        
        pro_multi_tot=profoundProFound(multi, segim=pro_detect$segim, mask=mask, magzero=magzero[i], gain=gain[i], boundstats=boundstats, iters=iters_tot, verbose=verbose)
        
        # Append column names and concatenate cat_tot together
        
        setnames(pro_multi_tot$segstats, paste0(names(pro_multi_tot$segstats), '_', catappend[i], totappend))
        cat_tot=cbind(cat_tot, pro_multi_tot$segstats)
      }
      
      if(docol){
        
        # Compute colour multi band photometry, allowing some extra dilation via the iters_tot argument
        # If we have already run the total photometry then we use the sky and skyRMS computed there for speed
        
        if(dotot){
          pro_multi_col=profoundProFound(multi, segim=pro_detect$segim_orig, mask=mask, sky=pro_multi_tot$sky, skyRMS=pro_multi_tot$skyRMS, redosky=FALSE, magzero=magzero[i], gain=gain[i], objects=pro_detect$objects, iters=0, verbose=verbose)$segstats
        }else{
          pro_multi_col=profoundProFound(multi, segim=pro_detect$segim_orig, mask=mask, magzero=magzero[i], gain=gain[i], objects=pro_detect$objects, iters=0, verbose=verbose)$segstats
        }
        
        # Append column names and concatenate cat_col together
        
        setnames(pro_multi_col, paste0(names(pro_multi_col), '_', catappend[i], colappend))
        cat_col=cbind(cat_col, pro_multi_col)
      }
    }
  }
  
  # Return all of the things!
  
  output=list(pro_detect=pro_detect, cat_tot=cat_tot, cat_col=cat_col, detectbands=detectbands, multibands=multibands, call=match.call(), date=date(), time=proc.time()[3]-timestart)
  class(output)='profoundmulti'
  return=output
}
