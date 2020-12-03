\name{profoundIm}
\alias{profoundImBlur}
\alias{profoundImGrad}
\alias{profoundImDiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Transformations
}
\description{
Various image transformation functions that assist in exploring data. These all require the \code{imager} package to be installed.
}
\usage{
profoundImBlur(image = NULL, sigma = 1, plot = FALSE, ...)
profoundImGrad(image = NULL, sigma = 1, plot = FALSE, ...)
profoundImDiff(image = NULL, sigma = 1, plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse.
}
  \item{sigma}{
Numeric scalar; standard deviation of the blur.
}
  \item{plot}{
Logical; should a \code{\link{magimage}} plot of the output be generated?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magimage}}. Only relevant is \option{plot}=TRUE.
}
}
\value{
Numeric matrix; a new image the same size as \option{image}, with the relevant transform applied.

For \code{profoundImBlur} the output is a smoothed version of the \option{image}.

For \code{profoundImGrad} the output is the magnitude of the gradient of the smoothed version of the \option{image}.

For \code{profoundImDiff} the output is the original \option{image} minus the smoothed version of the \option{image}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundMakeSegim}}, \code{\link{profoundMakeSegimExpand}}
}
\examples{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))$imDat
magimage(image)
profoundImBlur(image, plot=TRUE)
profoundImGrad(image, plot=TRUE)
profoundImDiff(image, plot=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ image }% use one of  RShowDoc("KEYWORDS")
