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
groupsegID = NULL, magzero = 0,  df = 3, radtrunc = 2, iterative = FALSE,
doallstats = TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. As a convenience you can supply the output of \code{\link{profoundProFound}} of class profound, in which case any required input that is not explicitly set via the arguments will be inherited from the \code{\link{profoundProFound}} list.
}
  \item{segim}{
Integer matrix; a specified segmentation map of the image. This matrix *must* be the same dimensions as \option{image} if supplied. If this option is used then \code{\link{profoundProFound}} will not compute its initial segmentation map using \code{\link{profoundMakeSegim}}, which is then dilated. Instead it will use the one passed through \option{segim}.
}
  \item{segstats}{
Data.frame, segmentation catalogue output from \option{profoundProFound}.
}
  \item{groupim}{
Integer matrix; the grouped segmentation map. This matrix *must* be the same dimensions as \option{image}.
}
  \item{groupsegID}{
List; group information as output by \option{profoundSegimGroup} or \option{profoundProFound}. Must correspond to the supplied \option{segstats}.
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
}
\details{
This routine only deblends with detected groups, so it is quite fast if the number of groups is quite low. If the image is more confused then this process can be quite slow.
}
\value{
A data.frame containing deblended flux information:

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

profound=profoundProFound(image, magzero=30, verbose=TRUE, plot=TRUE)

deblend=profoundFluxDeblend(profound)

plot(profound$segstats$mag, deblend$mag_db)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ deblend }% use one of  RShowDoc("KEYWORDS")
