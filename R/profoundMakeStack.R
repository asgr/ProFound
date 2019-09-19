profoundMakeStack=function(image_list=NULL, sky_list=NULL, skyRMS_list=NULL, magzero_in=0, magzero_out=0, masking='and'){
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
    if(!masking %in% c('and', 'or')){
      stop('masking needs to be one of and/or.')
    }
    stack=matrix(0,dim(image_list[[1]])[1],dim(image_list[[1]])[2])
    inv_var=0
    masked_and_master={}
    masked_and_initial=TRUE
    for(i in 1:length(image_list)){
      masked={}
      if(is.list(sky_list) & is.list(skyRMS_list)){
        if(is.null(image_list[[i]])==FALSE & is.null(sky_list[[i]])==FALSE & is.null(skyRMS_list[[i]])==FALSE){
          if(masking=='and'){
            if(anyNA(image_list[[i]])){
              masked=which(is.na(image_list[[i]]))
              image_list[[i]][masked]=0
              sky_list[[i]][masked]=0
              skyRMS_list[[i]][masked]=Inf
            }
          }
          image_list[[i]]=image_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
          sky_list[[i]]=sky_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
          skyRMS_list[[i]]=skyRMS_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
          stack=stack+(image_list[[i]]-sky_list[[i]])/(skyRMS_list[[i]]^2)
          inv_var=inv_var+(1/skyRMS_list[[i]]^2)
        }else{
          message(paste('Missing data in image_list element',i,'so will skip for stacking!'))
        }
      }
      if(is.list(sky_list) & is.list(skyRMS_list)==FALSE){
        if(is.null(image_list[[i]])==FALSE & is.null(sky_list[[i]])==FALSE){
          if(masking=='and'){
            if(anyNA(image_list[[i]])){
              masked=which(is.na(image_list[[i]]))
              image_list[[i]][masked]=0
              sky_list[[i]][masked]=0
            }
          }
          image_list[[i]]=image_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
          sky_list[[i]]=sky_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
          stack=stack+(image_list[[i]]-sky_list[[i]])
          inv_var=1
        }else{
          message(paste('Missing data in image_list element',i,'so will skip for stacking!'))
        }
      }
      if(is.list(sky_list)==FALSE & is.list(skyRMS_list)){
        if(is.null(image_list[[i]])==FALSE & is.null(skyRMS_list[[i]])==FALSE){
          if(masking=='and'){
            if(anyNA(image_list[[i]])){
              masked=which(is.na(image_list[[i]]))
              image_list[[i]][masked]=0
              skyRMS_list[[i]][masked]=Inf
            }
          }
          image_list[[i]]=image_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
          skyRMS_list[[i]]=skyRMS_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
          stack=stack+image_list[[i]]/(skyRMS_list[[i]]^2)
          inv_var=inv_var+(1/skyRMS_list[[i]]^2)
        }else{
          message(paste('Missing data in image_list element',i,'so will skip for stacking!'))
        }
      }
      if(is.list(sky_list)==FALSE & is.list(skyRMS_list)==FALSE){
        if(is.null(image_list[[i]])==FALSE){
          if(masking=='and'){
            if(anyNA(image_list[[i]])){
              masked=which(is.na(image_list[[i]]))
              image_list[[i]][masked]=0
            }
          }
          stack=stack+image_list[[i]]*profoundMag2Flux(magzero_in[i],magzero_out)
          inv_var=inv_var+1
        }else{
          message(paste('Missing data in image_list element',i,'so will skip for stacking!'))
        }
      }
      
      if(masking=='and'){
        if(masked_and_initial){
          masked_and_master = masked
          masked_and_initial = FALSE
        }else{
          masked_and_master = masked_and_master[masked_and_master %in% masked]
        }
      }
    }
    stack=stack/inv_var
    skyRMS=sqrt(1/inv_var)
    
    if(masking=='and' & length(masked_and_master)>0){
      stack[masked_and_master] = NA
    }
  }else{
    if(!is.null(sky_list)){
      stack=image_list-sky_list
    }else{
      stack=image_list
    }
    if(!is.null(skyRMS_list)){
      skyRMS=skyRMS_list
    }else{
      skyRMS=NULL
    }
  }
invisible(list(image=stack, skyRMS=skyRMS, magzero=magzero_out))
}
