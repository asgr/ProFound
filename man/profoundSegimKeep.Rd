\name{profoundSegimKeep}
\alias{profoundSegimKeep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge Segmentation Map with Grouped Segmentation map
}
\description{
Allows users to safely merge a standard segim with a groupim, where you can specify segments to be newly merged together, or groups to be merged.
}
\usage{
profoundSegimKeep(segim, groupim, groupID_merge, segID_merge, clean = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segim}{
Integer matrix; required, the segmentation map.
}
  \item{groupim}{
Integer matrix; the grouped segmentation map. This matrix *must* be the same dimensions as \option{segim} (if supplied).
}
  \item{groupID_merge}{
Integer vector; the group IDs that the user wants to persist into the final segmentation map (removing all \option{segim} segments that overlap with any of the specified group IDs).
}
  \item{segID_merge}{
Integer list; each list element should specify collections of segments to be merged.
}
  \item{clean}{
Logical; should segments partially overlapping with chosen groups be aggressively removed?
}
}
\details{
The merged segments inherit the lowest segment value, e.g. list(c(1,2,4),c(5,6)) would merge together segments 1,2,4 and to be a new segment 1, and then 5,6 to be a new segment 5.
}
\value{
Integer matrix; the merged segmentation map, where specified groups and segments have been merged.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundSegimMerge}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound=profoundProFound(image, magzero=30, groupstats=TRUE, verbose=TRUE, plot=TRUE)

segim_new=profoundSegimKeep(profound$segim, profound$group$groupim, groupID_merge=1,
segID_merge=list(c(12, 26, 62), c(13, 24)))

profoundSegimPlot(image, segim=segim_new)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ segim }% use one of  RShowDoc("KEYWORDS")
\keyword{ merge }% __ONLY ONE__ keyword per line