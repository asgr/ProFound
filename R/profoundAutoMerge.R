profoundAutoMerge = function(segim, segstats, spur_lim=4e-3, col_lim=NULL){
  
  if(is.null(segstats$col) | is.null(col_lim)){
    sel = which(segstats$cenfrac / (1 + segstats$sep) < spur_lim
    )
  }else{
    sel = which(segstats$col > col_lim[1] &
                segstats$col < col_lim[2] &
                segstats$cenfrac / (1 + segstats$sep) < spur_lim
    )
  }
  
  group = NULL
  
  if(length(sel) > 0){
    merge_segIDS = segstats[sel,"segID"]
    segim_temp = segim
    segim_temp[!segim_temp %in% merge_segIDS] = 0
    group = profoundSegimGroup(segim_temp)$groupsegID
    if(any(group$Ngroup==1)){
      merge_segIDS = merge_segIDS[!merge_segIDS %in% group[group$Ngroup==1,"groupID"]]
      segim_temp[!segim_temp %in% merge_segIDS] = 0
      if(any(segim_temp > 0)){
        group = profoundSegimGroup(segim_temp)$groupsegID
      }
    }
  }
  return(group)
}
