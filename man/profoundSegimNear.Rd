\name{profoundSegimNear}
\alias{profoundSegimNear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segment Neighbour IDs
}
\description{
Returns a data.frame of all nearby (default is touching) segments surrounding every segment in a provided segim.
}
\usage{
profoundSegimNear(segim = NULL, offset = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segim}{
Integer matrix; a specified segmentation map of the image (required).
}
  \item{offset}{
Integer scalar; the distance to offset when searching for nearby segments.
}
}
\details{
This function can be run by the user directly, but usually it is called from within a higher routine in the \code{ProFound} suite of objects detection functions.
}
\value{
A data.frame of lists giving the segIDs of nearby segments for every segment. This is a slightly unusal structure to see in R, but it allows for a compact manner of storing uneven vectors of touching segmentss. E.g. you might have a massive segment touching 30 other segments and many segments touching none. Padding a normal matrix out to accommodate the larger figure would be quite inefficient.

  \item{segID}{Segmentation ID, which can be matched against values in \option{segim}}
  \item{nearID}{An embedded list of segmentation IDs for nearby segments. I.e. each list element of \option{nearID} is a vector (see Examples for clarity).}
  \item{Nnear}{The total number of segments that are considered to be nearby.}

}
\author{
Aaron Robotham
}
\note{
Due to the construction of the segmented curve-of-growth in ProFound you may have cases where the separation between segments is two or three pixels. Since these are very close to touching you might want to catch these close neighbours rather than strictly touching. By increasing \option{offset} to a larger number (2 or 3 in the cases above) you can flag these events.
}

\seealso{
\code{\link{profoundProFound}}, \code{\link{profoundMakeSegim}}, \code{\link{profoundMakeSegimDilate}}, \code{\link{profoundMakeSegimExpand}}, \code{\link{profoundSegimStats}}, \code{\link{profoundSegimPlot}}
}

\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE)

#Look for nearby (in this case touching) neighbours

near=profoundSegimNear(profound$segim)

#Look at the first few rows (segIDs 1:5):

near[1:5,]

#To access the embedded vectors you have to use unlist:

unlist(near[3,2])

#We can check to see which segments are touching segID number 3:

profoundSegimPlot(image$imDat, profound$segim)
magimage(profound$segim==3, col=c(NA,'red'), add=TRUE)
magimage(matrix(profound$segim \%in\% unlist(near[3,2]), dim(profound$segim)[1]),
col=c(NA,'blue'), add=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segments }
