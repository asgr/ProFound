profoundSegimFix=function(image=NULL, segim=NULL, mask=NULL, sky=NULL, profound=NULL, 
                          loc=NULL, box=400, segID_merge=list(), col='magenta', pch=4, 
                          cex=2, crosshair=FALSE, crosscex=5, alpha_seg=0.3, happy_default=TRUE, 
                          continue_default=TRUE, open_window=TRUE, allow_seg_modify=FALSE, 
                          segID_max=NULL, legend_extra=NULL, group_limit=TRUE, ...){
  if(open_window){
    dev.new(noRStudioGD = TRUE)
  }
  
  if(!is.null(image)){
    if(inherits(image, 'profound')){
      if(is.null(segim)){segim=image$segim}
      if(is.null(mask)){mask=image$mask}
      if(is.null(sky)){sky=image$sky}
      image=image$image
      if(is.null(image)){stop('Need image in profound object to be non-Null')}
    }
  }
  if(!is.null(profound)){
    if(!inherits(image, 'profound')){
      stop('Class of profound input must be of type \'profound\'')
    }
    if(is.null(image)){image=profound$image}
    if(is.null(image)){stop('Need image in profound object to be non-Null')}
    if(is.null(segim)){segim=profound$segim}
    if(is.null(mask)){mask=profound$mask}
    if(is.null(sky)){sky=profound$sky}
  }
  
  if(!is.null(segID_max)){
    if(segID_max=='auto'){
      segID_max = max(segim, na.rm=TRUE)
    }
  }
  
  if(!is.null(loc)){
    if(is.list(image)){
      imageR = magcutout(image$R, loc=loc, box=box)$image
      imageG = magcutout(image$G, loc=loc, box=box)$image
      imageB = magcutout(image$B, loc=loc, box=box)$image
      image = list(R=imageR, G=imageG, B=imageB)
    }else{
      image = magcutout(image, loc=loc, box=box)$image
    }
    segim = magcutout(segim, loc=loc, box=box)$image
    if(!is.null(mask)){
      mask = magcutout(mask, loc=loc, box=box)$image
    }else{
      mask = NULL
    }
    if(!is.null(sky)){
      sky = magcutout(sky, loc=loc, box=box)$image
    }else{
      sky = NULL
    }
  }
  
  if(group_limit){
    if(!requireNamespace("imager", quietly = TRUE)){
      stop('The imager package is needed for this function to work. Please install it from CRAN', call. = FALSE)
    }
    groupim = as.matrix(imager::label(imager::as.cimg(segim>0)))
  }
  
  segim_start = segim
  segim_progress=segim_start
  if(length(segID_merge)>0){
    segim = profoundSegimKeep(segim_start, segID_merge=segID_merge)
    segim_progress[!segim_progress %in% unlist(segID_merge)] = 0
  }else{
    segim_progress[]=0
  }
  
  continue = TRUE
  quit = FALSE
  goback = FALSE
  
  par(mar=c(0.1,0.1,0.1,0.1))
  if(is.list(image)){
    magimageRGB(R=image$R, G=image$G, B=image$B, axes=FALSE, labels=FALSE, ...)
  }else{
    magimage(image, axes=FALSE, labels=FALSE, ...)
  }
  profoundSegimPlot(image=image, segim=segim, mask=mask, sky=sky, axes=FALSE, labels=FALSE, add=TRUE)
  
  while(continue){
    magimage(segim_progress, magmap=FALSE, col=c(NA, hsv(seq(0,2/3,len=max(segim_progress, na.rm=TRUE)), alpha=alpha_seg)), add=TRUE)
    if(crosshair){
      points(dim(segim)[1]/2,dim(segim)[2]/2, col='magenta', pch=5, cex=crosscex)
    }
    if(allow_seg_modify){
      legend('topleft', legend=c('ESC to stop','Seg 1: merge/mod','Seg 2: ignore','Seg 3: skip check','Grp 1: rem seg','Grp 2: ungroup','Sky 1: go back', 'Sky 2: quit', 'Sky 3: make seg'), text.col='magenta', bg='black')
    }else{
      legend('topleft', legend=c('ESC to stop','Seg 1: merge','Seg 2: ignore','Seg 3: skip check','Grp 1: rem seg','Grp 2: ungroup','Sky 1: go back', 'Sky 2: quit'), text.col='magenta', bg='black')
    }
    if(!is.null(legend_extra)){
      legend('topright', legend=legend_extra, text.col='magenta', bg='black')
    }
    cat('Click on contiguous segments to merge, and hit ESC in the plot window (not this one) when done.\n')
    
    check = NULL #this will be the tabulation of clicked segments, mainly used to count multi-clicks on segments and sky etc
    seggrid = NULL #used to store in-concave hull pixel IDs when drawing segments
    mergeIDs = NULL
    temploc = locator(type='p', col=col, pch=pch, cex=cex)
    
    if(is.null(temploc)){ #if nothing is selected move on
      legend('bottomleft', legend='Done!', text.col='magenta', bg='black')
      break
    }else{ #otherwise prepare for trouble...
      if(any(ceiling(temploc$x)<0) | any(ceiling(temploc$y)<0) | any(ceiling(temploc$x)>dim(segim)[1]) | any(ceiling(temploc$y)>dim(segim)[2])){
        legend('bottomleft', legend='Selected outside of image!', text.col='magenta', bg='black')
      }else{
        mergeIDs = segim_start[cbind(ceiling(temploc$x),ceiling(temploc$y))]
        check = tabulate(mergeIDs)
        if(all(mergeIDs>0)){
          mergeIDs = which(check %% 2 == 1)
          if(group_limit){
            groupIDs = groupim[cbind(ceiling(temploc$x),ceiling(temploc$y))]
            if(any(groupIDs==0)){
              mergeIDs = NULL
              check = NULL
              legend('bottomleft', legend='Will not merge with sky (will ignore)!', text.col='magenta', bg='black')
            }else if(length(unique(groupIDs))>1){
              mergeIDs = NULL
              check = NULL
              legend('bottomleft', legend='Discontinuous segments (will ignore)!', text.col='magenta', bg='black')
            }
          }
        }
      }
      
      if(!is.null(check)){
        if(length(check)==1 & length(mergeIDs) == 1 & all(mergeIDs==0)){ #flag to go back
          #Will set to go back if exactly one sky pixel and nothing else has been selected
          mergeIDs = NULL
          check = NULL
          legend('bottomleft', legend='Go back', text.col='magenta', bg='black')
          goback = TRUE
        }else if(length(check)==1 & length(mergeIDs) == 2 & all(mergeIDs==0)){ #quit!
          #Will set to quit if exactly two sky pixels and nothing else have been selected
          mergeIDs = NULL
          legend('bottomleft', legend='Quit', text.col='magenta', bg='black')
          quit = TRUE
        }else if(length(check)==1 & length(mergeIDs) == 3 & all(mergeIDs==0) & allow_seg_modify){ #entering segment polygon mode
          #Will set to enter segment making mode if exactly two sky pixels and nothing else have been selected
          mergeIDs = NULL
          legend('bottomleft', legend='Segment polygon mode', text.col='magenta', bg='black')
          temploc = locator(type='o', col=col, pch=pch, cex=cex)
          if(!is.null(temploc)){
            lines(temploc$x[c(1,length(temploc$x))], temploc$y[c(1,length(temploc$y))], col='magenta')
            seggrid = .inpoly_pix(temploc$x, temploc$y) #new c++ code to find pixels in segment
          }
        }else if(any(mergeIDs==0)){
          mergeIDs = NULL
          legend('bottomleft', legend='Ambiguous segment/sky click/s (will ignore)!', text.col='magenta', bg='black')
        }else if(length(which(check>0))==1){
          #If only one segment appears in the check vector
          mergeIDs = NULL
          if(max(check)==1 & allow_seg_modify==FALSE){ #select segment for group removal
            legend('bottomleft', legend='Will remove segment from group', text.col='magenta', bg='black')
            mergeIDs = NULL
            delete_seg = which(check==1)
            for(i in 1:length(segID_merge)){
              segID_merge[[i]] = segID_merge[[i]][!(segID_merge[[i]] %in% delete_seg)]
            }
            segim = profoundSegimKeep(segim_start, segID_merge=segID_merge)
            segim_progress = segim_start
            segim_progress[!segim_progress %in% unlist(segID_merge)] = 0
          }else if(max(check)==1 & allow_seg_modify){ #select segment for segID mod
            mergeIDs = NULL
            #If that object has been clicked twice flag for seg ID modification
            seggrid = which(segim==which(check==1), arr.ind=TRUE)
          }else if(max(check)==2){ #group de-merge
            legend('bottomleft', legend='Will de-merge group', text.col='magenta', bg='black')
            zapIDs = which(check>0)
            segID_merge = profoundZapSegID(zapIDs, segID_merge)
            segim = profoundSegimKeep(segim_start, segID_merge=segID_merge)
            segim_progress = segim_start
            segim_progress[!segim_progress %in% unlist(segID_merge)] = 0
          }
        }else{
          legend('bottomleft', legend='Merge segments', text.col='magenta', bg='black')
        }
      }
    }
    
    if(any(check==3)){
      happy=TRUE
    }else{
      if(happy_default){
        legend('topright', legend='Happy? [y]/n', text.col='magenta', bg='black')
        cat('HAPPY with your solution? [y]/n: ')
        happy = readLines(n=1L)
        happy = tolower(happy)
        happy = happy == "" | happy == 'yes' | happy == 'y' | happy == 't' | happy == 'true'
      }else{
        legend('topright', legend='Happy? y/[n]', text.col='magenta', bg='black')
        cat('HAPPY with your solution? y/[n]: ')
        happy = readLines(n=1L)
        happy = tolower(happy)
        happy = happy == 'yes' | happy == 'y' | happy == 't' | happy == 'true'
      }
    }
    
    if(happy){
      if(!is.null(seggrid) & allow_seg_modify){
        legend('topright', legend='segID: [auto]/#: ', text.col='magenta', bg='black')
        cat('segID: [auto]/#')
        newsegID = readLines(n=1L)
        newsegID = tolower(newsegID)
        if(newsegID=='auto' | newsegID==''){
          if(is.null(segID_max)){
            segID_max = max(segim, na.rm=TRUE)
          }
          newsegID = segID_max + 1L
          segID_max = segID_max + 1L
        }else{
          newsegID = as.integer(newsegID)
        }
        segim[seggrid]=newsegID
        segim_progress[seggrid]=newsegID
        segim_start[seggrid]=newsegID
      }
      if(length(mergeIDs)>0){
        segID_merge = c(segID_merge,list(mergeIDs))
        segID_merge = profoundMergeSegID(segID_merge)
      }
      if(length(segID_merge)>0){
        segim = profoundSegimKeep(segim_start, segID_merge=segID_merge)
        segim_progress = segim_start
        segim_progress[!segim_progress %in% unlist(segID_merge)] = 0
      }else{
        segim = segim_start
        segim_progress[] = 0
      }
    }
    
    if(happy){
      if(any(check==3) | goback){
        continue = FALSE
      }else{
        par(mar=c(0.1,0.1,0.1,0.1))
        if(is.list(image)){
          magimageRGB(R=image$R, G=image$G, B=image$B, axes=FALSE, labels=FALSE, ...)
        }else{
          magimage(image, axes=FALSE, labels=FALSE, ...)
        }
        profoundSegimPlot(image=image, segim=segim, mask=mask, sky=sky, axes=FALSE, labels=FALSE, add=TRUE) 
        
        if(quit == FALSE){
          if(continue_default){
            legend('topleft', legend='Cont? [y]/n/q', text.col='magenta', bg='black')
            cat('CONTINUE fixing segments on current image? [y]/n/q: ')
            continue = readLines(n=1L)
            continue = tolower(continue)
            quit = continue == 'q' | continue == 'quit'
            continue = continue == "" | continue == 'yes' | continue == 'y' | continue == 't' | continue == 'true'
          }else{
            legend('topleft', legend='Cont? y/[n]/q', text.col='magenta', bg='black')
            cat('CONTINUE fixing segments on current image? y/[n]/q: ')
            continue = readLines(n=1L)
            continue = tolower(continue)
            quit = continue == 'q' | continue == 'quit'
            continue = continue == 'yes' | continue == 'y' | continue == 't' | continue == 'true'
          }
        }else{
          continue = FALSE
        }
      }
    }else{
      continue = TRUE
      quit = FALSE
      par(mar=c(0.1,0.1,0.1,0.1))
      if(is.list(image)){
        magimageRGB(R=image$R, G=image$G, B=image$B, axes=FALSE, labels=FALSE, ...)
      }else{
        magimage(image, axes=FALSE, labels=FALSE, ...)
      }
      profoundSegimPlot(image=image, segim=segim, mask=mask, sky=sky, axes=FALSE, labels=FALSE, add=TRUE) 
    }
  }
  
  if(open_window){
    dev.off()
  }
  
  return(invisible(list(segim=segim, segim_start=segim_start, segID_merge=segID_merge, segID_max=segID_max, goback=goback, quit=quit)))
}

profoundSegimKeep = function(segim=NULL, groupim=NULL, groupID_merge=NULL, segID_merge=NULL, clean=FALSE){
  segim_out=segim
  
  if(! is.null(groupID_merge)){
    groupID_merge = groupID_merge[groupID_merge %in% groupim]
    if(clean){
      removeID = segim[groupim %in% groupID_merge]
      segim_out[segim_out %in% removeID] = 0
    }
    segim_out[groupim %in% groupID_merge] = groupim[groupim %in% groupID_merge]
  }
  
  if(! is.null(segID_merge)){
    if(! is.list(segID_merge)){
      stop('segID_merge must be a list of segments to be merged!')
    }
    whichpix = which(segim_out %in% unlist(segID_merge))
    pixsel = segim_out[whichpix]
    for(i in 1:length(segID_merge)){
      if(length(segID_merge[[i]])>0){
        if('fastmatch' %in% .packages()){
          pixsel[fastmatch::fmatch(pixsel, segID_merge[[i]], nomatch = 0L) > 0L] = min(segID_merge[[i]], na.rm=TRUE)
        }else{
          pixsel[pixsel %in% segID_merge[[i]]] = min(segID_merge[[i]], na.rm=TRUE)
        }
      }
    }
    segim_out[whichpix]=pixsel
  }
  
  return(invisible(segim_out))
}

profoundMergeSegID = function(segID_merge){
  if(is.list(segID_merge) & length(segID_merge)>1){
    unlist_segID = unlist(segID_merge)
    
    if(any(duplicated(unlist_segID))){
      
      groupID = unlist(lapply(segID_merge, min))
      times = unlist(lapply(segID_merge, length))
      
      segtab = data.frame(segID=unlist_segID, groupID=rep(groupID, times))
      
      for(i in 1:(length(groupID)-1)){
        for(j in (i+1):length(groupID)){
          if(i != j){
            if(any(segtab[segtab$groupID==groupID[i],'segID'] %in% segtab[segtab$groupID==groupID[j],'segID'])){
              segtab[segtab$groupID==groupID[j],'groupID'] = groupID[i]
            }
          }
        }
      }
      
      segID_merge_new = list()
      
      unique_groups = unique(segtab$groupID)
      
      for(i in unique_groups){
        segID_merge_new = c(segID_merge_new, list(unique(segtab[segtab$groupID==i,'segID'])))
      }
      
      return(invisible(segID_merge_new))
    }else{
      return(invisible(segID_merge))
    }
  }else{
    return(invisible(segID_merge))
  }
}

profoundZapSegID = function(segID, segID_merge){
  times = unlist(lapply(segID_merge, length))
  if(!is.null(times)){
    refs = rep(1:length(segID_merge), times=times)
    logic = unlist(segID_merge) %in% segID
    refs = unique(refs[logic])
    if(length(refs)>0){
      segID_merge = segID_merge[-refs]
    }
  }
  return(invisible(segID_merge))
}

.inpoly_pix=function(poly_x,poly_y){
  xseq=floor(min(poly_x)):ceiling(max(poly_x)) - 0.5
  yseq=floor(min(poly_y)):ceiling(max(poly_y)) - 0.5
  tempgrid=as.matrix(expand.grid(xseq, yseq))

  select=.point_in_polygon_cpp(tempgrid[,1], tempgrid[,2], poly_x, poly_y)
  
  return(invisible(ceiling(tempgrid[select,])))
}
