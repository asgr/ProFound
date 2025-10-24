this_in_that = function(this, that, invert=FALSE, type='logical', arr.ind=FALSE, nthreads=1L){
  if(!is.integer(this)){
    stop('this must be integer')
  }
  
  if(!is.integer(that)){
    stop('that must be integer')
  }
  
  output = .vec_this_in_vec_that(vec_this=this, vec_that=that, invert=invert, nthreads=nthreads)
  
  # if(is.vector(this)){
  #   output = .vec_this_in_vec_that(vec_this=this, vec_that=that, invert=invert, nthreads=nthreads)
  # }else if(is.matrix(this)){
  #   output = .mat_this_in_vec_that(mat_this=this, vec_that=that, invert=invert, nthreads=nthreads)
  # }
  
  if(is.matrix(this)){
    dim(output) = dim(this)
  }
  
  if(type == 'logical'){
    return(output)
  }else if(type == 'which'){
    return(which(output, arr.ind=arr.ind))
  }
}

`%fin%` = function(this, that){
  return(this_in_that(this, that))
}

`%nin%` = function(this, that){
  return(this_in_that(this, that, invert=TRUE))
}
