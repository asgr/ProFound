profoundMakeStack=function(image_list, sky_list=NULL, skyRMS_list=NULL, magzero_in=0, magzero_out=0){
  if(is.list(image_list)){
    if(is.null(sky_list)==FALSE){
      if(length(image_list)!=length(sky_list)){stop('sky_list length does not match image_list length!')}
    }
    if(is.null(skyRMS_list)==FALSE){
      if(length(image_list)!=length(skyRMS_list)){stop('skyRMS_list length does not match image_list length!')}
    }
    if(length(magzero_in)!=length(image_list)){
      if(length(magzero_in)!=1){stop('Length of magzero_in should be the number of images in image_list or 1')}
      magzero_in=rep(magzero_in,length(image_list))
    }
    stack=matrix(0,dim(image_list[[1]])[1],dim(image_list[[1]])[2])
    inv_var=0
    for(i in 1:length(image_list)){
      if(is.list(sky_list) & is.list(skyRMS_list)){
        image_list[[i]]=image_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
        sky_list[[i]]=sky_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
        skyRMS_list[[i]]=skyRMS_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
        stack=stack+(image_list[[i]]-sky_list[[i]])/(skyRMS_list[[i]]^2)
        inv_var=inv_var+(1/skyRMS_list[[i]]^2)
      }
      if(is.list(sky_list) & is.list(skyRMS_list)==FALSE){
        image_list[[i]]=image_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
        sky_list[[i]]=sky_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
        stack=stack+(image_list[[i]]-sky_list[[i]])
        inv_var=1
      }
      if(is.list(sky_list)==FALSE & is.list(skyRMS_list)){
        image_list[[i]]=image_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
        skyRMS_list[[i]]=skyRMS_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
        stack=stack+image_list[[i]]/(skyRMS_list[[i]]^2)
        inv_var=inv_var+(1/skyRMS_list[[i]]^2)
      }
      if(is.list(sky_list)==FALSE & is.list(skyRMS_list)==FALSE){
        stack=stack+image_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
        inv_var=inv_var+1
      }
    }
    stack=stack/inv_var
    skyRMS=sqrt(1/inv_var)
  }else{
    if(!missing(sky_list)){
      stack=image_list-sky_list
    }else{
      stack=image_list
    }
    if(!missing(skyRMS_list)){
      skyRMS=skyRMS_list
    }else{
      skyRMS=1
    }
  }
return=list(image=stack, skyRMS=skyRMS, magzero=magzero_out)
}