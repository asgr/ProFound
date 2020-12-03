\name{profoundMakeSegimPropagate}
\alias{profoundMakeSegimPropagate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Propagate Identified Segments
}
\description{
Propagates all identified segments across the full image, only ignoring masked regions. This serves to identify which segment every pixel is most likely to belong to using a number of image related criteria. Uses EBImage's propagate function to do the grunt work.
}
\usage{
profoundMakeSegimPropagate(image = NULL, segim = NULL, objects = NULL, mask = NULL,
sky = 0, lambda = 1e-04, plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. Note, image NAs are treated as masked pixels.
}
  \item{segim}{
Integer matrix; required, the segmentation map of the image. This matrix *must* be the same dimensions as \option{image}.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{sky}{
User provided estimate of the absolute sky level. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel). This will be subtracted off the \option{image} internally, so only provide this if the sky does need to be subtracted!
}
  \item{lambda}{
A numeric value. The regularization parameter used in the metric, determining the trade-off between the Euclidean distance in the image plane and the contribution of the gradient of x. See Details.
}
  \item{plot}{Logical; should a diagnostic plot be generated? This is useful when you only have a small number of sources (roughly a few hundred). With more than this it can start to take a long time to make the plot!
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This function propgates out the identified segments into the rest of the \option{image}, only region identified in the \option{mask} will not be assigned to a segment. To assign pixels a mixture of the Euclidian distance and the local gradient is used (as described below). The purpose of this routine is to identify all pixels in the image with their most likely segment (whether nominally object or sky pixel). The true sky pixels identified as belonging to a segment should also provide the best possible local estimate of the sky level.

For internal completeness, the below description is taken almost verbatim from the EBImage propagate function.

The method operates by computing a discretized approximation of the Voronoi regions for given seed points on a Riemann manifold with a metric controlled by local \option{image} features.

Under this metric, the infinitesimal distance d between points v and v+dv is defined by:

d^2 = ( (t(dv)*g)^2 + lambda*t(dv)*dv )/(lambda + 1) 
, where g is the gradient of \option{image} x at point v.

\option{lambda} controls the weight of the Euclidean distance term. When \option{lambda} tends to infinity, d tends to the Euclidean distance. When \option{lambda} tends to 0, d tends to the intensity gradient of the \option{image}.

The gradient is computed on a neighborhood of 3x3 pixels.

Segmentation of the Voronoi regions in the vicinity of flat areas (having a null gradient) with small values of \option{lambda} can suffer from artifacts coming from the metric approximation.
}
\value{
A list containing two images:

  \item{propim }{The propagated segmentation map including the original segments identified.}
  \item{propim_sky }{The propagated segmentation map removing the original segments identified (these pixels are set to 0).}

}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundProFound}} \code{propagate}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)

tempprop=profoundMakeSegimPropagate(image$imDat, segim=profound$segim, plot=TRUE)

tempprop_stats=profoundSegimStats(image$imDat, segim=tempprop$propim_sky,
sky=profound$sky, skyRMS=profound$skyRMS)

magplot(profound$segstats$mag, tempprop_stats$flux/tempprop_stats$N100, grid=TRUE)

#You can stop the propogation using a mask:

mask=array(0, dim=dim(image$imDat))
mask[1:50,]=1

profoundMakeSegimPropagate(image$imDat, segim=profound$segim, plot=TRUE, mask=mask)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ propagate }% use one of  RShowDoc("KEYWORDS")
