\name{profoundSegimMerge}
\alias{profoundSegimMerge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Merge Segmentation Maps
}
\description{
Takes two segmentation maps and merges them in a sensible manner, making sure segments representing the same object are not overlaid on each other.
}
\usage{
profoundSegimMerge(image = NULL, segim_base = NULL, segim_add = NULL, mask = NULL,
sky = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. Note, image NAs are treated as masked pixels.
}
  \item{segim_base}{
Integer matrix; required, the base segmentation map of the \option{image}. This matrix *must* be the same dimensions as \option{image}.

}
  \item{segim_add}{
Integer matrix; required, the new segmentation map of the \option{image} that is to be added. This matrix *must* be the same dimensions as \option{image}.
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{sky}{
User provided estimate of the absolute sky level. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel). This will be subtracted off the \option{image} internally, so only provide this if the sky does need to be subtracted!
}
}
\details{
The merger strategy is quite simple. Matching object segments are identified by the 'uniqueID' ID from an internal run of \code{\link{profoundSegimStats}}. Whichever segment contains more flux is determined to be the best map to use as the base segment. Unmatched segments in the \option{segim_add} map are added back in after this initial merging process, so will end up on top and potentially appear as segment islands within larger segments (which is not possible using the standard segmentation process in \code{\link{profoundMakeSegim}}).

An obvious reason to use this function is in situations where bright stars are embedded deep within an extended source. The standard watershed segmentation used in \code{\link{profoundMakeSegim}} will tend to break a large portion of the extended source off to form the segmented region. By running \code{\link{profoundProFound}} in different modes it is possible to identify the bright peaks (see Examples below), and then use \code{\link{profoundSegimMerge}} to piece the segments back together appropriately.
}
\value{
Integer matrix; the merged segmentation map matched pixel by pixel to \option{image}.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundMakeSegim}}, \code{\link{profoundSegimKeep}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))$imDat
profound=profoundProFound(image, plot=TRUE)
profound_diff=profoundProFound(profoundImDiff(image, sigma=2), plot=TRUE)
tempmerge=profoundSegimMerge(image, profound$segim, profound_diff$segim)

#Notice the new embedded blue segment near the centre:

profoundSegimPlot(image, segim=tempmerge)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segim }% use one of  RShowDoc("KEYWORDS")
