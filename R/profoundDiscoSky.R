profoundDiscoSky = function(image, disco, sky_arg_list=NULL, roughpedestal=FALSE, mask=NULL, ...){
  dots=list(...)
  
  sky = matrix(NA, dim(image)[1], dim(image)[2])
  skyRMS = matrix(NA, dim(image)[1], dim(image)[2])
  
  disoc_IDs = unique(as.integer(na.omit(disco)))
  
  for(i in disoc_IDs){
    if(is.null(mask)){
      temp_mask = ! disco == i
    }else{
      temp_mask = (! disco == i) | mask > 0L
    }
    replace = which(!temp_mask)
    if(is.null(sky_arg_list[[i]])){
      temp_sky = profoundMakeSkyGrid(image=image, mask=temp_mask*1L, ...)
    }else{
      temp_sky = do.call("profoundMakeSkyGrid", c(list(image=image, mask=temp_mask*1L), sky_arg_list[[i]], dots))
    }
    if(roughpedestal){
      sky[replace] = median(temp_sky$sky[replace], na.rm=TRUE)
      skyRMS[replace] = median(temp_sky$skyRMS[replace], na.rm=TRUE)
    }else{
      sky[replace] = temp_sky$sky[replace]
      skyRMS[replace] = temp_sky$skyRMS[replace]
    }
  }
  
  return(invisible(list(sky=sky, skyRMS=skyRMS)))
}