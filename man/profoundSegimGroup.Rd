\name{profoundSegimGroup}
\alias{profoundSegimGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Segmentation Groups
}
\description{
Given an input segmentation map, returns a map of groups of touching segments as well as the IDs of segments within each group.
}
\usage{
profoundSegimGroup(segim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segim}{
Integer matrix; required, the segmentation map.
}
}
\details{
To use this function you will need to have EBImage installed. Since this can be a bit cumbersome on some platforms (given its dependencies) this is only listed as a suggested package. You can have a go at installing it by running:

> source("http://bioconductor.org/biocLite.R")

> biocLite("EBImage")

Linux users might also need to install some non-standard graphics libraries (depending on your install). If you do not have them already, you should look to install **jpeg** and **tiff** libraries (these are apparently technically not entirely free, hence not coming by default on some strictly open source Linux variants).

\code{profoundSegimGroup} uses the \code{bwlabel} function from EBImage.
}
\value{
A list containting the following structures:

  \item{groupim}{An map of the unique groups identified in the input \option{segim}, where the groupID is the same as the lowest valued segID in the group.}
  \item{groupsegID}{A data.frame of lists giving the segIDs of segments in each group.}

The data.frame returned by \option{groupsegID} is a slightly unusal structure to see in R, but it allows for a compact manner of storing uneven vectors of grouped segments. E.g. you might have a massive group containing 30 other segments and many groups containing a single segment. Padding a normal matrix out to accommodate the larger figure would be quite inefficient. It contains the following:

  \item{groupID}{Group ID, which can be matched against values in \option{groupim}}
  \item{segID}{An embedded list of segmentation IDs for segments in the group. I.e. each list element of \option{segID} is a vector (see Examples for clarity).}
  \item{Ngroup}{The total number of segments that are in the group.}
  \item{Npix}{The total number of pixels that are in the group.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundSegimNear}}, ~~~
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE)

#Look for nearby (in this case touching) neighbours

group=profoundSegimGroup(profound$segim)

#Look at the first few rows (groups 1:5):

group$groupsegID[1:5,]

#To access the embedded vectors you have to use unlist:

unlist(group$groupsegID[1,2])

#We can check to see which segments are in group number 1:

profoundSegimPlot(image$imDat, profound$segim)
magimage(group$groupim==1, col=c(NA,'red'), add=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segments }% use one of  RShowDoc("KEYWORDS")
