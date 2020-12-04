profoundMultiBand=function(inputlist=NULL, dir='', segim=NULL, mask=NULL, detectbands='r',
                           multibands=c('u','g','r','i','z'), iters_det=6, iters_tot=0,
                           sizes_tot=5, magzero=0, gain=NULL, box=100, grid=box, boxadd=box/2,
                           app_diam=1, bandappend=multibands, totappend='t', colappend='c',
                           grpappend='g', dotot=TRUE, docol=TRUE, dogrp=TRUE, deblend=FALSE,
                           groupstats=FALSE, groupby_det='segim_orig', groupby_mul='segim_orig',
                           keepsegims=FALSE, masking='and', ...){
  
  # The most important thing is that all of the input images must be pixel matched via SWarp or magwarp etc
  # detectbands and multibands are the names of the target bands, which should be the names of the images ignoring the .fits ending
  # By default Simon's script names them as u.fits and g.fits etc, hence the literal bands naming
  # detectbands can be a vector of target detection bands for stacking
  # multibands is by default a vector of all SDSS optical bands
  # You can select whether to do total and/or colour photometry with dotot and docol
  # magzero must be matched to the number of bands in multibands (which reflects the magzero in each respectively), or one value which is then recycled
  # bandappend are the band names to use as the extra tag attached to the column names in the cat_tot and cat_col outputs (by default this is linked to multibands)
  # totappend and colapped is just the extra tag attached to the column names in the cat_tot and cat_col outputs
  
  timestart=proc.time()[3]
  
  call=match.call()
  
  dots=list(...)
  
  dotsignoredetect=c('iters', 'sky', 'skyRMS', 'plot', 'stats', 'haralickstats', 'groupby', 'pixelcov', 'box')
  dotsignoremulti=c('skycut', 'pixcut', 'tolerance', 'ext', 'sigma', 'smooth', 'iters', 'size', 'sky', 'skyRMS', 'plot', 'stats', 'redosegim', 'roughpedestal', 'haralickstats', 'groupby', 'box', 'groupstats', 'objects', 'redosky')
  
  if(length(dots)>0){
    dotsdetect=dots[! names(dots) %in% dotsignoredetect]
    dotsmulti=dots[! names(dots) %in% dotsignoremulti]
  }else{
    dotsdetect={}
    dotsmulti={}
  }
  
  # Restrict outselves to data actually present (no matter what is asked for)
  
  if(is.null(inputlist)){
    presentbands=unlist(strsplit(list.files(dir), split='.fits'))
  }else{
    if(length(inputlist)!=length(multibands)){
      stop('inputlist is not the same length as multibands!')
    }
    presentbands=multibands
  }
  
  if(multibands[1]=='get'){
    multibands=presentbands
  }
  
  if(detectbands[1]=='get'){
    detectbands=presentbands
  }
  if(detectbands[1]=='all'){
    detectbands=multibands
  }
  detectbands=detectbands[detectbands %in% presentbands]
  
  if(bandappend[1]=='get'){
    bandappend=presentbands
  }
  
  if(length(iters_tot)==1){
    iters_tot=rep(iters_tot, length(multibands))
  }
  
  if(length(iters_tot)!=length(iters_tot)){
    stop('Length of iters_tot must equal length of multibands!')
  }
  
  if(length(sizes_tot)==1){
    sizes_tot=rep(sizes_tot, length(multibands))
  }
  
  if(length(sizes_tot)!=length(sizes_tot)){
    stop('Length of sizes_tot must equal length of multibands!')
  }
  
  if(length(magzero)==1){
    magzero=rep(magzero, length(multibands))
  }
  
  if(length(magzero)!=length(multibands)){
    stop('Length of magzero must equal length of multibands!')
  }
  
  if(length(gain)==1){
    gain=rep(gain, length(multibands))
  }
  
  if(is.null(gain)==FALSE & length(gain)!=length(multibands)){
    stop('Length of gain must equal length of multibands!')
  }
  
  if(length(box)==1){
    box=rep(box, length(multibands))
  }
  
  if(length(box)!=length(multibands)){
    stop('Length of box must equal length of multibands!')
  }
  
  if(length(grid)==1){
    grid=rep(grid, length(multibands))
  }
  
  if(length(grid)!=length(multibands)){
    stop('Length of grid must equal length of multibands!')
  }
  
  if(length(boxadd)==1){
    boxadd=rep(boxadd, length(multibands))
  }
  
  if(length(boxadd)!=length(multibands)){
    stop('Length of boxadd must equal length of multibands!')
  }
  
  if(length(app_diam)==1){
    app_diam=rep(app_diam, length(multibands))
  }
  
  if(length(app_diam)!=length(app_diam)){
    stop('Length of app_diam must equal length of multibands!')
  }
  
  #if(dogrp & boundstats==FALSE){
  #  stop('If dogrp=TRUE than boundstats=TRUE must also be set!')
  #}
  
  magzero=magzero[which(multibands %in% presentbands)]
  magzero=magzero[!is.na(magzero)]
  
  bandappend=bandappend[which(multibands %in% presentbands)]
  bandappend=bandappend[!is.na(bandappend)]
  
  if(!is.null(gain)){
    gain=gain[which(multibands %in% presentbands)]
    gain=gain[!is.na(gain)]
  }
  
  iters_tot=iters_tot[which(multibands %in% presentbands)]
  iters_tot=iters_tot[!is.na(iters_tot)]
  
  sizes_tot=sizes_tot[which(multibands %in% presentbands)]
  sizes_tot=sizes_tot[!is.na(sizes_tot)]
  
  box=box[which(multibands %in% presentbands)]
  box=box[!is.na(box)]
  
  grid=grid[which(multibands %in% presentbands)]
  grid=grid[!is.na(grid)]
  
  boxadd=boxadd[which(multibands %in% presentbands)]
  boxadd=boxadd[!is.na(boxadd)]
  
  multibands=multibands[which(multibands %in% presentbands)]
  multibands=multibands[!is.na(multibands)]
  
  segimlist=NULL
  
  # Some safety checks
  
  if(!is.null(gain)){
    if(length(unique(c(length(multibands), length(magzero), length(bandappend), length(gain))))!=1){
      stop('multibands, magzero, bandappend and gain are not all the same length - they must be!')
    }
  }else{
    if(length(unique(c(length(multibands), length(magzero), length(bandappend))))!=1){
      stop('multibands, magzero and bandappend are not all the same length - they must be!')
    }
  }
  if(!is.null(inputlist)){
    if(length(unique(c(length(multibands),length(inputlist))))!=1){
      stop('inputlist does not match the length of multibands - it must be!')
    }
  }
  
  message(paste('*** Will use',paste(detectbands,collapse='+'),'for source detection ***'))
  if(!is.null(segim)){
    message('*** Will use provided segim for detection statistics ***')
  }
  
  if(dotot | docol | dogrp){
    message(paste('*** Will use',paste(multibands,collapse=' '),'for multi band photometry ***'))
    message(paste('*** Magzero:',paste(multibands,magzero,sep='=', collapse=' '),' ***'))
    if(!is.null(gain)){
      message(paste('*** Gain:',paste(multibands,gain,sep='=', collapse=' '),' ***'))
    }else{
      message('*** Gain: not specified for any bands, so shot-noise will be ignored ***')
    }
    
    if(dotot){
      message('*** Will compute total multi band photometry ***')
    }
    if(docol){
      message('*** Will compute isophotal colour multi band photometry ***')
    }
    if(dogrp){
      message('*** Will compute grouped segment multi band photometry ***')
    }
  }
  
  if(length(detectbands)==1){
    
    # If only one band is specified for detection then we skip the stacking part
    
    message(paste('*** Currently processing single detection band',detectbands,'***'))
    if(is.null(inputlist)){
      #if(requireNamespace("Rfits", quietly = TRUE)){
      #  detect=Rfits::Rfits_read_image(paste0(dir,detectbands,'.fits'))
      if(requireNamespace("FITSio", quietly = TRUE)){
        detect=FITSio::readFITS(paste0(dir,detectbands,'.fits'))
      }else{
        stop('One of Rfits / FITSio is required to read in FITS images. Get from GitHub asgr/Rfits / CRAN.')
      }
    }else{
      detect=inputlist[[which(multibands==detectbands)]]
    }
    temp_magzero=magzero[multibands==detectbands]
    temp_box=box[multibands==detectbands]
    temp_grid=grid[multibands==detectbands]
    temp_boxadd=boxadd[multibands==detectbands]
    
    # pro_detect=profoundProFound(image=detect, segim=segim, mask=mask, skycut=skycut, pixcut=pixcut, tolerance=tolerance, ext=ext, sigma=sigma, smooth=smooth, iters=iters_det, magzero=temp_magzero, verbose=verbose, boundstats=boundstats, groupstats=(groupstats | dogrp), groupby=groupby, haralickstats=haralickstats, ...)
    
    pro_detect=do.call("profoundProFound", c(list(image=quote(detect), segim=quote(segim), mask=quote(mask), iters=iters_det, magzero=temp_magzero, box=temp_box, grid=temp_grid, boxadd=temp_boxadd, deblend=FALSE, groupstats=(groupstats | dogrp), pixelcov=FALSE, group_by=groupby_det), dotsdetect))
    
  }else{
  
    # Multiple detection bands requested, so we prepare lists for stacking
    
    detect_image=list()
    detect_sky=list()
    detect_skyRMS=list()
    detect_magzero={}
    detect_box={}
    detect_grid={}
    detect_boxadd={}
    
    for(i in 1:length(detectbands)){
      
      # Loop around detection bands
      
      message(paste('*** Currently processing detection band',detectbands[i],'***'))
      
      if(is.null(inputlist)){
        detect=readFITS(paste0(dir,detectbands[i],'.fits'))
      }else{
        detect=inputlist[[which(multibands==detectbands[i])]]
      }
      
      temp_magzero=magzero[multibands==detectbands[i]]
      temp_box=box[multibands==detectbands[i]]
      temp_grid=grid[multibands==detectbands[i]]
      temp_boxadd=boxadd[multibands==detectbands[i]]
      
      # Run ProFound on current detection band with input parameters
      
      # pro_detect=profoundProFound(image=detect, segim=segim, mask=mask, skycut=skycut, pixcut=pixcut, tolerance=tolerance, ext=ext, sigma=sigma, smooth=smooth, iters=iters_det, magzero=temp_magzero, verbose=verbose, ...)
      
      pro_detect=do.call("profoundProFound", c(list(image=quote(detect), segim=quote(segim), mask=quote(mask), iters=iters_det, magzero=temp_magzero, box=temp_box, grid=temp_grid, boxadd=temp_boxadd, deblend=FALSE, groupstats=FALSE, pixelcov=FALSE), dotsdetect))
      
      # mask things per band
      
      if(!is.null(pro_detect$mask)){
        pro_detect$image[pro_detect$mask==1]=NA
      }
      
      # Append to lists for stacking
      
      detect_image=c(detect_image,list(pro_detect$image))
      detect_sky=c(detect_sky,list(pro_detect$sky))
      detect_skyRMS=c(detect_skyRMS,list(pro_detect$skyRMS))
      detect_magzero=c(detect_magzero, temp_magzero)
      detect_box=c(detect_box, temp_box)
      detect_grid=c(detect_grid, temp_grid)
      detect_boxadd=c(detect_boxadd, temp_boxadd)
    }
    
    # Grab the current header
    
    header=pro_detect$header
    
    # Delete and clean up
    
    rm(detect)
    rm(pro_detect)
    
    # Stack!!!
    # We first stack the image then the sky.
    
    detect_image_stack=profoundMakeStack(image_list=detect_image, sky_list=detect_sky, skyRMS_list=detect_skyRMS, magzero_in=detect_magzero, magzero_out=detect_magzero[1], masking=masking)
    detect_sky_stack=profoundMakeStack(image_list=detect_sky, skyRMS_list=detect_skyRMS, magzero_in=detect_magzero, magzero_out=detect_magzero[1])
    
    # Delete and clean up
    
    rm(detect_image)
    rm(detect_skyRMS)
    
    message('*** Currently processing stacked detection image ***')
    
    # For reference we run ProFound with the stacked sky added back in, passing it the stacked sky too.
    # Mask is NULL since any masked pixels will already be NA in the image.
    
    #pro_detect=profoundProFound(image=detect_image_stack$image+detect_sky_stack$image, segim=segim, mask=NULL, header=header, skycut=skycut, pixcut=pixcut, tolerance=tolerance, ext=ext, sigma=sigma,  smooth=smooth, iters=iters_det, magzero=detect_magzero[1], sky=detect_sky_stack$image, skyRMS=detect_image_stack$skyRMS, redosky=FALSE, verbose=verbose, boundstats=boundstats, groupstats=(groupstats | dogrp), groupby=groupby_det, haralickstats=haralickstats, ...)
    
    pro_detect=do.call("profoundProFound", c(list(image=quote(detect_image_stack$image+detect_sky_stack$image), segim=quote(segim), mask=NULL, header=header, iters=iters_det, magzero=detect_magzero[1], sky=quote(detect_sky_stack$image), skyRMS=quote(detect_image_stack$skyRMS), box=max(detect_box), grid=min(detect_grid), boxadd=max(detect_boxadd), deblend=FALSE, groupstats=(groupstats | dogrp), groupby=groupby_det, pixelcov=FALSE), dotsdetect))
    
    # Reset image pixels to NA if masked.
    
    if(!is.null(pro_detect$mask)){
      pro_detect$image[pro_detect$mask==1]=NA
    }
    
    if(is.null(pro_detect$segim)){
      stop('STOPPING: detection segmentation map appears to be empty! Please check input image data.')
    }
    
    # Delete and clean up
    
    rm(detect_image_stack)
    rm(detect_sky_stack)
  }
  
  #pro_detect$call=NULL
  
  # Create base total and colour photometry catalogues

  cat_tot=NULL
  cat_col=NULL
  cat_grp=NULL
  
  if(dotot | docol | dogrp){
    
    if(dotot){
      cat_tot=data.frame(segID=pro_detect$segstats$segID)
    }
    
    if(docol){
      cat_col=data.frame(segID=pro_detect$segstats$segID)
    }
    
    if(dogrp){
      cat_grp=data.frame(groupID=pro_detect$group$groupsegID$groupID)
    }
    
    for(i in 1:length(multibands)){
      
      # Loop around multi bands
      
      message(paste('*** Currently processing multi band',multibands[i],'***'))
      
      if(is.null(inputlist)){
        multi=readFITS(paste0(dir,multibands[i],'.fits'))
      }else{
        multi=inputlist[[i]]
      }
      
      if(dotot){
        
        message('* Processing total photometry *')
        
        # Compute total multi band photometry, allowing some extra dilation via the iters_tot argument
        
        # pro_multi_tot=profoundProFound(image=multi, segim=pro_detect$segim, mask=mask, magzero=magzero[i], gain=gain[i], groupstats=FALSE, iters=iters_tot[i], size=sizes_tot[i], deblend=deblend, redosegim=FALSE, roughpedestal=FALSE, ...)
        
        pro_multi_tot=do.call("profoundProFound", c(list(image=quote(multi), segim=quote(pro_detect$segim), mask=quote(mask), magzero=magzero[i], gain=gain[i], box=box[i], grid=grid[i], boxadd=boxadd[i], groupstats=FALSE, iters=iters_tot[i], size=sizes_tot[i], deblend=deblend, groupby=groupby_mul, redosegim=FALSE, roughpedestal=FALSE, app_diam=app_diam[i]), dotsmulti))
        
        # Append column names and concatenate cat_tot together
        
        setnames(pro_multi_tot$segstats, paste0(names(pro_multi_tot$segstats), '_', bandappend[i], totappend))
        cat_tot=cbind(cat_tot, pro_multi_tot$segstats)
        
        if(keepsegims){
          segimlist=c(segimlist, list(pro_multi_tot$segim))
        }
        
      }
      
      if(docol){
        
        message('* Processing colour photometry *')
        
        # Compute colour multi band photometry
        # If we have already run the total photometry then we use the sky and skyRMS computed there for speed
        
        if(dotot){
          #pro_multi_col=profoundProFound(image=multi, segim=pro_detect$segim_orig, mask=mask, sky=pro_multi_tot$sky, skyRMS=pro_multi_tot$skyRMS, redosky=FALSE, magzero=magzero[i], gain=gain[i], objects=pro_detect$objects, boundstats=boundstats, groupstats=FALSE, iters=0, verbose=verbose, ...)$segstats
          
          pro_multi_col=do.call("profoundProFound", c(list(image=quote(multi), segim=quote(pro_detect$segim_orig), mask=quote(mask), sky=quote(pro_multi_tot$sky), skyRMS=quote(pro_multi_tot$skyRMS), redosky=FALSE, magzero=magzero[i], gain=gain[i], box=box[i], grid=grid[i], boxadd=boxadd[i], objects=quote(pro_detect$objects), groupstats=FALSE, iters=0, deblend=FALSE, redosegim=FALSE, roughpedestal=FALSE), dotsmulti))$segstats
        }else{
          #pro_multi_col=profoundProFound(image=multi, segim=pro_detect$segim_orig, mask=mask, magzero=magzero[i], gain=gain[i], objects=pro_detect$objects, boundstats=boundstats, groupstats=FALSE, iters=0, verbose=verbose, ...)$segstats
          
          pro_multi_col=do.call("profoundProFound", c(list(image=quote(multi), segim=quote(pro_detect$segim_orig), mask=quote(mask), magzero=magzero[i], gain=gain[i], box=box[i], grid=grid[i], boxadd=boxadd[i], objects=quote(pro_detect$objects), groupstats=FALSE, iters=0, deblend=FALSE, redosegim=FALSE, roughpedestal=FALSE), dotsmulti))$segstats
        }
        
        # Append column names and concatenate cat_col together
        
        setnames(pro_multi_col, paste0(names(pro_multi_col), '_', bandappend[i], colappend))
        cat_col=cbind(cat_col, pro_multi_col)
      }
      
      if(dogrp){
        
        message('* Processing group photometry *')
        
        # Compute group multi band photometry
        # If we have already run the total photometry then we use the sky and skyRMS computed there for speed
        
        if(dotot){
        #   # pro_multi_grp=profoundProFound(image=multi, segim=pro_detect$group$groupim, mask=mask, sky=pro_multi_tot$sky, skyRMS=pro_multi_tot$skyRMS, redosky=FALSE, magzero=magzero[i], gain=gain[i], objects=pro_detect$objects, groupstats=FALSE, iters=iters_tot[i], size=sizes_tot[i], deblend=FALSE, redosegim=FALSE, roughpedestal=FALSE, ...)$segstats
        #   
          pro_multi_grp=do.call("profoundProFound", c(list(image=quote(multi), segim=quote(pro_detect$group$groupim), mask=quote(mask), sky=quote(pro_multi_tot$sky), skyRMS=quote(pro_multi_tot$skyRMS), redosky=FALSE, magzero=magzero[i], gain=gain[i], box=box[i], grid=grid[i], boxadd=boxadd[i], objects=quote(pro_detect$objects), groupstats=FALSE, iters=iters_tot[i], size=sizes_tot[i], deblend=FALSE, redosegim=FALSE, roughpedestal=FALSE), dotsmulti))$segstats
        }else{
        #   # pro_multi_grp=profoundProFound(image=multi, segim=pro_detect$group$groupim, mask=mask, magzero=magzero[i], gain=gain[i], objects=pro_detect$objects, boundstats=boundstats, groupstats=FALSE, iters=0, verbose=verbose, ...)$segstats
        #   
          pro_multi_grp=do.call("profoundProFound", c(list(image=quote(multi), segim=quote(pro_detect$group$groupim), mask=quote(mask), magzero=magzero[i], gain=gain[i], box=box[i], grid=grid[i], boxadd=boxadd[i], objects=quote(pro_detect$objects), groupstats=FALSE, iters=iters_tot[i], deblend=FALSE, redosegim=FALSE, roughpedestal=FALSE), dotsmulti))$segstats
        }
        
        # Append column names and concatenate cat_grp together
        names(pro_multi_grp)[1]='groupID'
        setnames(pro_multi_grp, paste0(names(pro_multi_grp), '_', bandappend[i], grpappend))
        cat_grp=cbind(cat_grp, pro_multi_grp)
      }
    }
  }
  
  # Return all of the things!
  
  if(dotot & keepsegims){
    names(segimlist)=multibands
  }
  
  output=list(pro_detect=pro_detect, cat_tot=cat_tot, cat_col=cat_col, cat_grp=cat_grp, segimlist=segimlist, detectbands=detectbands, multibands=multibands, call=call, date=date(), time=proc.time()[3]-timestart)
  class(output)='profoundmulti'
  invisible(output)
}
