\name{profoundSegimExtend}
\alias{profoundSegimExtend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find Additional Image Segments
}
\description{
This function takes a target segmentation map, and finds additional sources on a target image. This might be useful where you want to do matched aperture photometry, but there is the possibility that the target image might have genuine additional sources that are not in your current segmentation map. This might cause issues with the sky estimation etc (even if you do not actually care about getting photometry for these additional sources).
}
\usage{
profoundSegimExtend(image = NULL, segim = NULL, mask = segim, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. Note, image NAs are treated as masked pixels.
}
  \item{segim}{
Integer matrix; required, the segmentation map of the \option{image}. This matrix *must* be the same dimensions as \option{image}.
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{\dots}{
Further arguments to be passed to \code{\link{profoundProFound}}.
}
}
\details{
This function is mostly useful if you are chaining together a script of detection images based on certain bands that are then run on data which might have very different (and previously undetected) astrophysical sources.
}
\value{
Integer matrix; the segmentation map matched pixel by pixel to \option{image} and \option{segim}. Newly identified segments are appended in number to the input \option{segim}, so if the maximum previous segment ID was 10, the new sources would start from 11 etc.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundProFound}}, \code{\link{profoundMakeSegim}}
}

\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound=profoundProFound(image, magzero=30, plot=TRUE)

#Now we remove the brightest central source:
profound$segim[profound$segim==1]=0

#And now we find it again:
segim_extend=profoundSegimExtend(image=image, segim=profound$segim, skycut=2)
profoundSegimPlot(image=image$imDat, segim=segim_extend)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segim }% use one of  RShowDoc("KEYWORDS")
