\name{profoundFluxDeblend}
\alias{profoundFluxDeblend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mid Level Image Deblender
}
\description{
Given a target image, a segmentation map, image segstats and group properties, this function will attempt a non-parametric deblend based on local fitting of B-splines to create a weight map for each segment in a group. Flux is guaranteed to be conserved, and errors are appropriately rescaled.
}
\usage{
profoundFluxDeblend(image = NULL, segim = NULL, segstats = NULL, groupim = NULL,
groupsegID = NULL, sky = 0, profound = NULL, magzero = 0,  df = 3, radtrunc = 2,
iterative = FALSE, doallstats = TRUE, lowmemory = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. As a convenience you can supply the output of \code{\link{profoundProFound}} of class profound, in which case any required input that is not explicitly set via the arguments will be inherited from the \code{\link{profoundProFound}} list.
}
  \item{segim}{
Integer matrix; a specified segmentation map of the image. This matrix *must* be the same dimensions as \option{image} if supplied.
}
  \item{segstats}{
Data.frame, segmentation catalogue output from \option{profoundProFound}.
}
  \item{groupim}{
Integer matrix; the grouped segmentation map. This matrix *must* be the same dimensions as \option{image}. If missing then this is computed using \code{\link{profoundSegimGroup}} using the supplied \option{segim}.
}
  \item{groupsegID}{
List; group information as output by \option{profoundSegimGroup} or \option{profoundProFound}. Must correspond to the supplied \option{segstats}. If missing then this is computed using \code{\link{profoundSegimGroup}} using the supplied \option{segim}.
}
  \item{sky}{
Numeric; the absolute sky level. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel).  
}
  \item{profound}{
List; object of class 'profound'. If this is provided then missing input arguments are taking directly from this structure (see Examples). As an added convenience, you can assign the profound object directly to the \option{image} input.
}
  \item{magzero}{
Numeric scalar; the magnitude zero point.
}
  \item{df}{
Integer scalar; degrees of freedom for the non-parametric spline fitting. See \code{\link{smooth.spline}}.
}
  \item{radtrunc}{
Numeric scalar; the maximum allowed radius beyond the edge-most segment pixel to consider when deblending. Keeping this low (1-3) ensures segments do not gather flux from very distant regions of the group.
}
  \item{iterative}{
Logical; should each segment profile fit be subtracted as it goes along? TRUE tends to remove the pedestal from a large galaxy that has faint objects embedded on top.
}
  \item{doallstats}{
Logical; specifies whether the output catalogue is matched against all rows of the supplied \option{segstats} (TRUE), or only the rows containing grouped (and therefore deblended) galaxies are returned and the core flux columns (see below)s.
}
\item{lowmemory}{
Logical; if TRUE then a low memory mode of FluxDeblend will be used. This is quite a bit slower, so the default of FALSE is usually preferred unless you are running into serious memory issues.
}
}
\details{
This routine only deblends with detected groups, so it is quite fast if the number of groups is quite low. If the image is more confused then this process can be quite slow.
}
\value{
A data.frame containing deblended flux information:

  \item{groupID}{The group ID reference for the deblend (as taken from groupsegID)}
  \item{segID}{The segment ID reference for the deblend (as taken from groupsegID)}
  \item{flux_db}{Total flux (calculated using \option{image}-\option{sky}) in ADUs}
  \item{mag_db}{Total flux converted to mag using \option{magzero}}
  \item{N100_db}{Total number of pixels in this segment, i.e. contains 100\% of the flux}
  \item{flux_segfrac}{Fraction of group flux in this segment. If this is very low then it is likely harder to extract good quality fluxes (with or without deblending).}
  \item{Qseg_db}{Quality flag for deblended segment. This represents what fraction of the segment image flux is in our deblend model. Negative means the model misses flux, postive means it has too much.  Nearer to 0 is better.}
  \item{Qgroup_db}{Quality flag for deblended group. This represents what fraction of the group image flux is in our deblend model. Negative means the model misses flux, postive means it has too much.  Nearer to 0 is better.}
  
  The below are only returned if \option{doallstats}=TRUE:
  
  \item{flux_err_db}{Estimated total error in the flux for the segment}
  \item{mag_err_db}{Estimated total error in the magnitude for the segment}
  \item{flux_err_sky_db}{Sky subtraction component of the flux error}
  \item{flux_err_skyRMS_db}{Sky RMS component of the flux error}
  \item{flux_err_shot_db}{Object shot-noise component of the flux error (only if \option{gain} is provided)}
}
\author{
Aaron Robotham
}
\note{
Given the large number of inputs required, this function effectively needs \code{\link{profoundProFound}} to be run first.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profoundProFound}}, \code{\link{smooth.spline}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound=profoundProFound(image, magzero=30, verbose=TRUE, plot=TRUE, groupstats=TRUE)

deblend=profoundFluxDeblend(profound)

magplot(profound$segstats$mag, profound$segstats$mag-deblend$mag_db, ylim=c(-0.3,0.3),
grid=TRUE, xlab='mag', ylab='mag_orig - mag_deblend')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ deblend }% use one of  RShowDoc("KEYWORDS")
