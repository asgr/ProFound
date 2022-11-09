profoundSegimCompare = function(segim_1, segim_2, threshold=0.5, cores=1, ignoresky=TRUE){
  
  registerDoParallel(cores=cores)
  
  segim_max = max(segim_1, segim_2)
  
  segID = NULL
  
  seg1_in_seg2 = foreach(segID = 1:segim_max, .combine='rbind')%dopar%{
    if(! segID %in% segim_1){
      return(rep(NA, segim_max))
    }
    sel = which(segim_1 == segID)
    temp_tab = tabulate(segim_2[sel])
    if(ignoresky){
      temp_sum = sum(temp_tab, na.rm=TRUE)
    }else{
      temp_sum = length(sel)
    }
    output = rep(0, segim_max)
    if(temp_sum == 0){
      return(output)
    }else{
      output[1:length(temp_tab)] = temp_tab/temp_sum
      return(output)
    }
  }
  
  seg2_in_seg1 = foreach(segID = 1:segim_max, .combine='rbind')%dopar%{
    if(! segID %in% segim_2){
      return(rep(NA, segim_max))
    }
    sel = which(segim_2 == segID)
    temp_tab = tabulate(segim_1[sel])
    if(ignoresky){
      temp_sum = sum(temp_tab, na.rm=TRUE)
    }else{
      temp_sum = length(sel)
    }
    output = rep(0, segim_max)
    if(temp_sum == 0){
      return(output)
    }else{
      output[1:length(temp_tab)] = temp_tab/temp_sum
      return(output)
    }
  }
  
  colnames(seg1_in_seg2) = NULL
  rownames(seg1_in_seg2) = NULL
  colnames(seg2_in_seg1) = NULL
  rownames(seg2_in_seg1) = NULL
  
  seg_bij = seg1_in_seg2 * t(seg2_in_seg1)
  
  tab_match = which(seg_bij > threshold, arr.ind = TRUE)
  
  segim_1_unique = sort(unique(as.integer(segim_1)))[-1]
  segim_2_unique = sort(unique(as.integer(segim_2)))[-1]
  
  tab_match = cbind(tab_match, match(tab_match[,1], segim_1_unique), match(tab_match[,2], segim_2_unique))
  
  colnames(tab_match) = c('segID_1', 'segID_2', 'rowID_1', 'rowID_2')
  
  return(list(seg_bij=seg_bij, tab_match=as.data.frame(tab_match), seg1_in_seg2=seg1_in_seg2, seg2_in_seg1=seg2_in_seg1))
}
