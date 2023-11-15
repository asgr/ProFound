profoundApplyMask = function(image = NULL, mask='disc',
                             xcen=xsize/2, ycen=ysize/2, xsize=101, ysize=101, rot=0,
                             direction='backward', dim = c(101, 101)){
  
  if(!requireNamespace("imager", quietly = TRUE)){
    stop('The imager package is needed for smoothing to work. Please install from CRAN.', call. = FALSE)
  }
  
  if(!is.matrix(mask)){
    mask = profoundMakeMask(size=max(xsize, ysize), shape=mask)
  }
  
  if(!is.null(image)){
    dim = dim(image)
  }
  
  if(direction == 'backward'){
    .map.tran = .map.tran_backward
  }else if(direction == 'forward'){
    .map.tran = .map.tran_forward
  }
  
  N_xcen = length(xcen)
  
  if(N_xcen > 1){
    if(length(ycen) == 1){
      ycen = rep(ycen, N_xcen)
    }
    
    if(length(xsize) == 1){
      xsize = rep(xsize, N_xcen)
    }
    
    if(length(ysize) == 1){
      ysize = rep(ysize, N_xcen)
    }
    
    if(length(rot) == 1){
      rot = rep(rot, N_xcen)
    }
    
    if(length(direction) == 1){
      direction = rep(direction, N_xcen)
    }
  }
  
  temp_mask = matrix(0L, max(dim[1], xsize, dim(mask)[1]), max(dim[2], ysize, dim(mask)[2]))
  temp_mask_out = matrix(0, dim[1], dim[2])
  temp_mask[1:dim(mask)[1],1:dim(mask)[2]] = mask
  
  for(i in 1:N_xcen){
    formals(.map.tran)$in_size = dim(mask)
    formals(.map.tran)$out_size = c(xsize[i], ysize[i])
    formals(.map.tran)$xcen = xcen[i]
    formals(.map.tran)$ycen = ycen[i]
    formals(.map.tran)$rot = rot[i]
    
    local_fun = function(x, y){
      .map.tran(x = x, y = y)
    }
    
    temp_mask_loc = as.matrix(imager::imwarp(imager::as.cimg(temp_mask),
                              map=local_fun,
                              direction=direction[i],
                              coordinates='absolute',
                              interpolation='nearest'))[1:dim[1],1:dim[2]]
    temp_mask_out[temp_mask_out == 0L & temp_mask_loc > 0] = i
  }
  
  if(!is.null(image)){
    image[temp_mask_out > 0] = NA
  }
  
  return(list(mask=temp_mask_out, image=image))
}

.map.tran_forward = function(x=0, y=0, xcen, ycen, in_size, out_size, rot=0){
  x = x - in_size[1]/2
  y = y - in_size[2]/2
  x = x*(out_size[1]/in_size[1])
  y = y*(out_size[2]/in_size[2])
  if(rot != 0){
    newx = x*cos(-rot*pi/180) + y*sin(-rot*pi/180)
    newy = -x*sin(-rot*pi/180) + y*cos(-rot*pi/180)
    x = newx
    y = newy
  }
  x = x + xcen
  y = y + ycen
  list(x = x, y = y)
}

.map.tran_backward = function(x=0, y=0, xcen, ycen, in_size, out_size, rot=0){
  x = x - xcen
  y = y - ycen
  
  if(rot != 0){
    newx = x*cos(rot*pi/180) + y*sin(rot*pi/180)
    newy = -x*sin(rot*pi/180) + y*cos(rot*pi/180)
    x = newx
    y = newy
  }
  
  x = x/(out_size[1]/in_size[1])
  y = y/(out_size[2]/in_size[2])
  
  x = x + in_size[1]/2
  y = y + in_size[2]/2
  list(x = x, y = y)
}

profoundMakeMask = function(size=101, shape='disc'){
  return(.makeBrush(size, shape=shape))
}

profoundDrawMask = function(image, poly=NULL, invert_mask = FALSE, mode='draw', type='pix', poly.col='red', ...){
  
  isRfits = inherits(image, c('Rfits_image'))
  
  if(mode == 'draw'){
    if(isRfits){
      if(requireNamespace("Rfits", quietly = TRUE)){
        dev.new(noRStudioGD = TRUE)
        plot(image, ...)
      }else{
        stop("The Rfits package is need to process the header. Install from GitHub asgr/Rfits.")
      }
    }else{
      dev.new(noRStudioGD = TRUE)
      magimage(image, ...)
    }
    
    temploc = locator(type='l', col=poly.col)
    lines(temploc$x[c(length(temploc$x),1)], temploc$y[c(length(temploc$y),1)], col=poly.col)
    output = data.frame(x=temploc$x, y=temploc$y)
    
    if(isRfits){
      tempWCS = Rwcs::Rwcs_p2s(output, keyvalues=image$keyvalues, header=image$raw, pixcen='R')
      output$RA = tempWCS[,'RA']
      output$Dec = tempWCS[,'Dec']
    }
    
  }else if(mode == 'apply'){
    if(!is.null(poly)){
      output = poly
    }else{
      stop('You must provide poly!')
    }
    
    if(type == 'coord'){
      temploc = Rwcs::Rwcs_s2p(output$RA, output$Dec, keyvalues=image$keyvalues, header=image$raw, pixcen='R')
      output = data.frame(x=temploc[,'x'], y=temploc[,'y'], RA=output$RA, Dec=output$Dec)
    }else if(type == 'pix'){
      output = data.frame(x = output[,1], y = output[,2])
    }
    
  }
  
  image_grid = expand.grid(1:dim(image)[1] - 0.5, 1:dim(image)[2] - 0.5)
  
  if(requireNamespace("Rwcs", quietly = TRUE)){
    in_poly = Rwcs::Rwcs_in_poly(image_grid[,1], image_grid[,2], output$x, output$y)
    mask = matrix(0L, dim(image)[1], dim(image)[2])
    
    if(invert_mask){
      in_poly = (in_poly == FALSE)
    }
    
    mask[in_poly] = 1L
    
    if(mode == 'draw'){
      magimage(mask, col=c(NA,hsv(alpha=0.2)), add=T, magmap=FALSE)
    }
    
    if(isRfits){
      image$imDat[in_poly] = NA
    }else{
      image[in_poly] = NA
    }
    
  }else{
    mask = NULL
  }
  
  return(list(mask=mask, image=image, poly=output))
}
