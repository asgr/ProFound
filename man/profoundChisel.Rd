\name{profoundChisel}
\alias{profoundChisel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Noise Chisel Esque Sky
}
\description{
This is an implementation of a Noise Chisel like approach method of identifying sky pixels. It is conceptually similar on some regards, but does not attempt to produce the same solutions (i.e. it is more like "inspired by").
}
\usage{
profoundChisel(image = NULL, sky = NULL, skythresh = 0.005, blurcut = 0.01,
objthresh = 1 - skythresh, sigma = 2, iterchisel = 5, itersky = 5,
objbias = 0.5, box = 100, skyconv = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. Note, image NAs are treated as masked pixels.
}
  \item{sky}{
User provided estimate of the absolute sky level. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel). This will be subtracted off the \option{image} internally, so only provide this if the sky does need to be subtracted!
}
  \item{skythresh}{
Numeric scalar; the quantile level to use to identify sink pixels for the sky (so should be between 0 and 1).
}
  \item{blurcut}{
Numeric scalar; PDF level to cut the blurring kernel at when dilating. Lower values mean larger dilations.
}
  \item{objthresh}{
Numeric scalar; the quantile level to use to identify source pixels for real objects (so should be between 0 and 1).
}
  \item{sigma}{
Numeric scalar; standard deviation of the blur used. This should be well matched to the PSF of the image in most cases.
}
  \item{iterchisel}{
Integer scalar; how many iterations of the inner chisel routine to run for.
}
  \item{itersky}{
Integer scalar; how many iterations of the outer sky estimation routine to run for.
}
  \item{objbias}{
Numeric scalar; how much to bias the dilations towards assuming pixels belong to objects (larger values mean more object pixels will in general be created).
}
  \item{box}{
Integer scalar; the dimensions of the box car filter to estimate the sky with.
}
  \item{skyconv}{
Numeric scalar; when <= \option{skyconv} fraction of pixels change sky/object designation, the solution is considered converged and \option{itersky} iterations stop.
}
}
\details{
This is an implementation of a Noise Chisel like approach method of identifying sky pixels. It is conceptually similar on some regards, but does not attempt to produce the same solutions (i.e. it is more like "inspired by").
}
\value{
A list with two parts:

\item{objects}{Logical matrix; the object map matched pixel by pixel to \option{image}. 1 means there is an object at this pixel, 0 means it is a sky pixel. Can be used as a mask in various other functions that require objects to be masked out.}
\item{sky}{The estimated sky level of the \option{image}.}
}
\references{
Akhlaghi et al, 2015, ApJSS, 220, 1
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundMakeSkyGrid}}
}
\examples{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))$imDat
magimage(profoundChisel(image)$objects)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ sky }% use one of  RShowDoc("KEYWORDS")
