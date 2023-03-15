\name{profoundSegimWarp}
\alias{profoundSegimWarp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remap Segmentation Map via Rwcs Warping
}
\description{
Remaps an input segmentation map WCS Tan Gnomonic projection system to a different target WCS. This uses \code{Rwcs_warp} with sensible settings \option{doscale} = FALSE, \option{interpolation} = 'nearest', \option{direction} = 'backward', but \code{Rwcs_warp} can be used more directly if the other lower level options are required (Note the \code{Rwcs} package will be needed from asgr/Rwcs GitHub.

This interface should cover most practical use cases though. Using \code{\link{profoundProFound}} with a remapped segmentation map is likely to be more sensible than remapping image flux since it will not produce flux interpolation errors.
}
\usage{
profoundSegimWarp(segim_in = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segim_in}{
Integer matrix; required, the segmentation map we want to remap. If \option{segim_in} is a list as created by \code{Rfits} then the image part of the list is passed to \option{segim_in} and the correct header part is passed to \option{header_in}.
}
  \item{\dots}{
Additional arguments passed into \code{Rwcs_warp} (note you should not set \option{doscale}, \option{interpolation}, \option{direction}).
}
}
\details{
This function uses the \option{doscale} = FALSE, \option{interpolation} = 'nearest', \option{direction} = 'backward' in \code{Rwcs_warp}.
}
\value{
Integer matrix; the remapped image using the target WCS.
}
\author{
Aaron Robotham
}

\seealso{
\code{Rwcs_warp}
}
\examples{
\dontrun{
VST_r = Rfits_read_image(system.file("extdata", 'VST_r.fits', package="magicaxis"))
GALEX_NUV = Rfits_read_image(system.file("extdata", 'GALEX_NUV.fits',
  package="magicaxis"))

profound_KiDS=profoundProFound(VST_r, sky=0, skycut=1, sigma=2, tolerance=8, plot=TRUE)

segimFUV=profoundSegimWarp(profound_KiDS$segim, profound_KiDS$header, GALEX_NUV$hdr)

profoundSegimPlot(GALEX_NUV, segim = segimFUV)

profound_GALEX=profoundProFound(GALEX_NUV, segim=segimFUV, plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ WCS }% use one of  RShowDoc("KEYWORDS")
\concept{ gnomonic }
\concept{ warp }% __ONLY ONE__ keyword per line
