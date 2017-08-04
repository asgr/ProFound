profitSegimNear=function(segim, offset=1){
  
  if(length(segim)>1e6){rembig=TRUE}else{rembig=FALSE}
  if(rembig){
    invisible(gc())
  }
  
  xlen=dim(segim)[1]
  ylen=dim(segim)[2]
  
  segim_inner=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)]
  off_down=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)-offset]
  off_left=segim[(offset+1):(xlen-offset)-offset,(offset+1):(ylen-offset)]
  off_up=segim[(offset+1):(xlen-offset),(offset+1):(ylen-offset)+offset]
  off_right=segim[(offset+1):(xlen-offset)+offset,(offset+1):(ylen-offset)]
  
  tabcomb=data.table(segID=as.integer(segim_inner), down=as.integer(off_down), left=as.integer(off_left), up=as.integer(off_up), right=as.integer(off_right))
  
  if(rembig){
      rm(segim_inner)
      rm(off_down)
      rm(off_left)
      rm(off_up)
      rm(off_right)
      invisible(gc())
    }
  
  tabcomb=tabcomb[segID>0,]
  setorder(tabcomb,segID)
  
  segID=NULL; down=NULL; left=NULL; up=NULL; right=NULL; nearID=NULL; Nnear=NULL
  
  #The line means for each set of segID pixels, identify unique adjacent pixels that do not have the ID of the segID of interest or 0 (sky) and sort the touching ID list and return to a listed data.frame. Ouch!
  
  tabnear=tabcomb[,list(nearID=list(sort(setdiff(unique(c(down,left,up,right)),c(0,segID))))),by=segID]
  tabnear[,Nnear:=length(unlist(nearID)),by=segID]
  return=tabnear
}