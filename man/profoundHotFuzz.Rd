\name{profoundHotFuzz}
\alias{profoundHotFuzz}
\alias{profoundPSFFuzz}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hot Fuzz Fitter
}
\description{
Functions to fit the tenuous outer fuzz of extended sources.
}
\usage{
profoundHotFuzz(profound = NULL, loc = NULL, box = c(200, 200), size = 21, shape = "disc",
  select = NULL, rough = TRUE, nser = 2, dofit = TRUE, Niters = c(200, 0), fitRe = TRUE,
  axrat = 'profound', ...)
  
profoundPSFFuzz(profound = NULL, loc = NULL, box = c(200,200), size = 21, shape = "disc",
  select = NULL, rough = TRUE, nser = 2, dofit = TRUE, Niters = c(200, 0), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{profound}{
List, required; output from \code{\link{profoundProFound}}.
}
  \item{loc}{
Numeric vector; two element [x,y] position for cutout region.
}
  \item{box}{
Numeric vector; two element dimensions of the box to cut out from image centred on loc (in pixels). 
}
  \item{size}{
Integer scalar; the size (e.g. width/diameter) of the dilation kernel in pixels. Should be an odd number else will be rounded up to the nearest odd number.
}
  \item{shape}{
Character scalar; the shape of the dilation kernel. Options are 'box', 'disc', 'diamond', 'Gaussian', 'line'.
}
  \item{select}{
Vector; either logical vector of segstat rows to keep (same length as segstats), or integer vector of row numbers.
}
  \item{rough}{
Logical; should the profile be approximated roughly. Usually setting this to true is fine for the out parts of light profiles.
}
  \item{nser}{
Numeric scalar; the Sersic index of the profile to use when fitting. When fitting "Hot-Fuzz" the Magnitude, Re and Index become quite degenerate, so this can usually be left at the default.
}
  \item{dofit}{
Logical; should fitting be done (TRUE) or just produce the outputs that can be later fitted (FALSE)?
}
  \item{Niters}{
Integer vector; if \option{dofit} = TRUE, the number of CMA and MCMC steps to make with \code{Highlander}.
}
  \item{fitRe}{
Logical; if \option{dofit} = TRUE, should Re also be fitted?
}
  \item{axrat}{
Character or numeric; either set to 'profound' (so \option{axrat} in inherited from ProFound), or set to the desired value.  
}
  \item{\dots}{
Other arguments to pass to \code{Highlander}.  
}
}
\details{
These codes are probably not general user at the moment (although feel free to use them). They were produced with a mind to my (Aaron Robotham's) ARC Future Fellowship project.

In brief they construct a \code{ProFit} region map that is defined by the segmentation map from ProFound and a dilation of this (using the provided kernel arguments). This region map becomes the target for a large \code{ProFit} model that is in effect only attempting to fit the outer "Hot Fuzz" that is not inside the standard \code{ProFound} segmentation map. Because of how far out this fuzz usually easier, a computationally cheaper 'rough' model can usually be used. 
}
\value{
\code{profoundHotFuzz} fits the outer profile of target galaxies. This is used to model the "Hot Fuzz" in the outer parts of galaxy profiles. It inherits the ellipticity and orientation of the \code{\link{profoundProFound}} output, allowing variable magnitude and Re by default per source.

\code{profoundPSFFuzz} fits the outer profile of target stars. This is used to model the "Hot Fuzz" in the outer parts of star profiles. It assumes the stars are circular, allowing variable Re and Index shared for all sources (so a common PSF profile), and variable magnitude per source (different flux stars).

The either function the output is:

  If \option{dofit} = TRUE the first list elements are all the usual outputs from Highlander.
  
  \item{Modellist}{The best fit modellist (if \option{dofit} = TRUE).}
  \item{image_model}{The best fit model image (if \option{dofit} = TRUE).}
  
  \item{Data}{The fitting data for Highlander.}
  \item{image}{The sky subtracted fitting image.}
  \item{segim}{The fitting segmentation map.}
  \item{segim_redo}{The dilated fitting segmentation map.}
  \item{sigma}{The fitting segstats map.}
  \item{segstats}{The fitting segstats.}
}
\author{
Aaron Robotham
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profoundProFound}}, ~~~
}
\examples{
\dontrun{
image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))
profound = profoundProFound(image, magzero=30)

#This takes a few minutes to run
fuzz_fit = profoundHotFuzz(profound)

magimage(fuzz_fit$image, qdiff=TRUE)

no_fuzz = fuzz_fit$image
sel = profound$segim_orig == 0
no_fuzz[sel] = no_fuzz[sel] - fuzz_fit$image_model[sel] 

magimage(no_fuzz, qdiff=TRUE)
}
}

\concept{ fuzz }
