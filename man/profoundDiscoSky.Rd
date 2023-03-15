\name{profoundDiscoSky}
\alias{profoundDiscoSky}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Per Discontinuity Region Sky Map
}
\description{
This function allows the sky in an image to vary in a discontinuous manner between regions, as described by the \option{disco} discontinuity region map. The typical use case would be a CCD/CMOS sensor that has multichannel read out. This often creates extemely sharp features even within a single image that are not well captured for a smooth bicubic spline.
}
\usage{
profoundDiscoSky(image, disco, sky_arg_list = NULL, roughpedestal = FALSE, mask = NULL,
  ...)
}
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{disco}{
Integer matrix; required, the discontinuity region map. This should be numbered from 1:N for N regions, with no gaps. Regions should be contiguous, or weird things may happen.
}
  \item{sky_arg_list}{
List; a list of per discontinuity region \code{\link{profoundMakeSkyGrid}} arguments (see Examples). This allows fine control over difference regions of the image if required.
}
  \item{roughpedestal}{
Logical; should only a pedestal (based on the median of the sky/skyRMS map) be used for the sky/skyRMS? This is a good option if the image is known to contain a *very* large (many times the box size) galaxy that might otherwise be over subtracted by the initial rough sky map.
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{\dots}{
Other global arguments to pass into \code{\link{profoundMakeSkyGrid}}.
}
}
\details{
This function is not incorporated into the higher level \code{\link{profoundProFound}} function. In general it is only useful for extreme cases of big discontinuities within images due to different read out chips and/or artefacts created by image mosaicing (if this is not done optimally, you often see discrete jumps in the sky and skyRMS between regions.

Creting the \option{disco} discontinuity region map has to be done manually, but users should consider using the \option{sky} and \option{skyRMS} outputs of a standard run of \code{\link{profoundProFound}} for guidance.
}
\value{
Produces a list of 2 matrices. The first (called \option{sky}) is a matrix of values for the absolute sky; the second (called \option{skyRMS}) is a matrix of values for the absolute sky RMS.
}
\author{
Aaron Robotham
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profoundMakeSkyGrid}}
}
\examples{
\dontrun{
image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))$imDat
profound = profoundProFound(image) #We want this for the object mask.

#Create toy disco region map:
disco = matrix(1L,356,356) #firstly set everything to region 1.
disco[,1:178] = 2L #half the region is flagged as region 2.

#First we run with just an object mask and take a look:
sky = profoundDiscoSky(image, disco=disco, objects=profound$objects_redo)
#The discontinuity is obvious along the horizontal.
magimage(sky$sky)

#Different box sizes for the different disco regions (default is 100 otherwise):
sky = profoundDiscoSky(image, disco=disco, sky_arg_list=list(list(box=50),list(box=25)),
  objects=profound$objects_redo)
magimage(sky$sky)

#Just using a pedestal:
sky = profoundDiscoSky(image, disco=disco, objects=profound$objects_redo,
  roughpedestal=TRUE)
magimage(sky$sky)
table(sky$sky) #just two different sky values for the two disco regions.
}
}

\concept{ sky }
