profoundMakeStack=function(image_list, sky_list=NULL, skyRMS_list=NULL){
  if(is.list(image_list)){
    if(!missing(sky_list)){
      if(length(image_list)!=length(sky_list)){stop('sky_list length does not match image_list length!')}
    }
    if(!missing(skyRMS_list)){
      if(length(image_list)!=length(skyRMS_list)){stop('skyRMS_list length does not match image_list length!')}
    }
    stack=matrix(0,dim(image_list[[1]])[1],dim(image_list[[1]])[2])
    inv_var=0
    for(i in 1:length(image_list)){
      if(is.list(sky_list) & is.list(skyRMS_list)){
        stack=stack+(image_list[[i]]-sky_list[[i]])/(skyRMS_list[[i]]^2)
        inv_var=inv_var+(1/skyRMS_list[[i]]^2)
      }
      if(is.list(sky_list) & is.list(skyRMS_list)==FALSE){
        stack=stack+(image_list[[i]]-sky_list[[i]])
        inv_var=1
      }
      if(is.list(sky_list)==FALSE & is.list(skyRMS_list)){
        stack=stack+image_list[[i]]/(skyRMS_list[[i]]^2)
        inv_var=inv_var+(1/skyRMS_list[[i]]^2)
      }
      if(is.list(sky_list)==FALSE & is.list(skyRMS_list)==FALSE){
        stack=stack+image_list[[i]]
        inv_var=1
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
return=list(image=stack, skyRMS=skyRMS)
}