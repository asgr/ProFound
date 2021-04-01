profoundAutoMerge = function(segim, segstats, spur_lim=4e-3, col_lim=NULL, Ncut=1){
  
  segstats = as.data.frame(segstats)
  
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
    group_old = profoundSegimGroup(segim)
    merge_segIDS = merge_segIDS[!merge_segIDS %in% group_old$groupsegID[group_old$groupsegID$Ngroup == 1,'groupID']]
    if(length(merge_segIDS) > 0){
      segim_temp = segim
      segim_temp[!segim_temp %in% merge_segIDS] = 0
      group = profoundSegimGroup(segim_temp)$groupsegID
      if(any(group$Ngroup < Ncut)){
        merge_segIDS = merge_segIDS[!merge_segIDS %in% group[group$Ngroup < Ncut,"groupID"]]
        segim_temp[!segim_temp %in% merge_segIDS] = 0
        if(any(segim_temp > 0)){
          group = profoundSegimGroup(segim_temp)$groupsegID
        }
      }
      names(group)[1] = 'mergeID'
      names(group)[3] = 'Nmerge'
      group = cbind(group, groupID = group_old$groupim[ceiling(as.matrix(segstats[match(group$mergeID, segstats$segID),c('xmax','ymax')]))])
    }
  }
  return(group)
}
