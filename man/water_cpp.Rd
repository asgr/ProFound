\name{water_cpp}
\alias{water_cpp}
\alias{water_cpp_old}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rcpp Watershed Function
}
\description{
This is a standalone implementation of a watershed deblend, with some astronomy specific tweeks. E.g. it is possible to both adapt the extent of the saddlepoint search, and it can be modified by both an absolute and relative tolerance. Defaults behave much like EBImage's watershed function. In general it is a factor of a few faster than the EBImage implementation, especially for large images with lots of deblending required.
}
\usage{
water_cpp(image = 0L, nx = 1L, ny = 1L, abstol = 1, reltol = 0, cliptol = 1e+06, ext = 1L,
skycut = 0, pixcut = 1L, verbose = FALSE, Ncheck = 1000000L)
water_cpp_old(image = 0L, nx = 1L, ny = 1L, abstol = 1, reltol = 0, cliptol = 1e+06,
ext = 1L, skycut = 0, pixcut = 1L, verbose = FALSE, Ncheck = 1000000L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. Note, image NAs are treated as masked pixels.
}
  \item{nx}{
Integer scalar; required, the dimension x of the supplied \option{image}, i.e. should be dim(image)[1].
}
  \item{ny}{
Integer scalar; required, the dimension y of the supplied \option{image}, i.e. should be dim(image)[2].
}
  \item{abstol}{
Numeric scalar; the minimum height of the object in the units of image intensity between its highest point (seed) and the point where it contacts another object (checked for every contact pixel). If the height is smaller than the tolerance, the object will be combined with its brightest neighbour. Tolerance should be chosen according to the range of \option{image}. Default works well when the \option{image} has been divided by the sky-RMS. A larger value of \option{abstol} means segments are more aggressively merged together.
}
  \item{reltol}{
Numeric scalar; a modifier to the \option{abstol}, modifying it by the ratio of the segment peak flux divided by the saddle point flux to the power \option{reltol}. The default means the \option{reltol} has no effect since this modifier becomes 1. A larger value of \option{reltol} means segments are more aggressively merged together.
}
  \item{cliptol}{
Numeric scalar; if \option{image} is above this level where segments touch then they are always merged, regardless of other criteria. When thinking in terms of sky RMS, values between 20-100 are probably appropriate for merging very bright parts of stars back together in optical data.
}
  \item{ext}{
Integer scalar; square offset of the neighborhood in pixels for the detection of neighboring objects. Higher value smoothes out small objects.
}
  \item{skycut}{
Numeric scalar; background value at or under which pixels are not considered anymore for watersheding (so with default \option{skycut}=0, 0 valued pixels and below are considered to be part of the sky).
}
  \item{pixcut}{
Integer scalar; the minimum number of pixels allowed in a segment. Below this number segments are set to 0, i.e. the background. This means they are not considered real objects in \code{\link{profoundProFound}}.
}
  \item{verbose}{
Logical; should verbose output be displayed to the user? Since big image can take a long time to run, you might want to monitor progress.  
}
  \item{Ncheck}{
Integer scalar; the pixel scanning interval to check for interupts and for printing out the verbose state.  
}
}
\details{
This was hand written from scratch by A Robotham, but in the end the approach is somewhat similar to \code{EBImage::watershed}. There do seem to be fairly large speed improvements for more sparse images though, since only pixels above the background \option{skycut} are ever looked at. This knowledge of sparcity does not exist in \code{EBImage::watershed}.

\code{water_cpp} is the newer variant re-writeen by R Tobar based on the Rcpp implentation. The older Rcpp one is still available as \code{water_cpp_old}.
}
\value{
Integer matrix; the segmentation map matched pixel by pixel to \option{image}.
}
\references{
Some aspects of Meyer's floodfill used, but not explicitly based on any published approach, so might be in detail similar by accident.
}
\author{
Aaron Robotham
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profoundMakeSegim}}, ?EBImage::watershed
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))$imDat

segim=water_cpp(im=image, nx=dim(image_smooth)[1], ny=dim(image_smooth)[2], skycut=10)
magimage(segim, col=c(0,rainbow(1e3)))
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ watershed }% use one of  RShowDoc("KEYWORDS")
