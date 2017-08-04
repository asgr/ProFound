\name{profitImBlur}
\alias{profitImBlur}
\alias{profitImGrad}
\alias{profitImDiff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Image Transformations
}
\description{
Various image transformation functions that assist in exploring data.
}
\usage{
profitImBlur(image, sigma = 1, plot = FALSE, ...)
profitImGrad(image, sigma = 1, plot = FALSE, ...)
profitImDiff(image, sigma = 1, plot = FALSE, ...)
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

For \code{profitImBlur} the output is a smoothed version of the \option{image}.

For \code{profitImGrad} the output is the magnitude of the gradient of the smoothed version of the \option{image}.

For \code{profitImDiff} the output is the original \option{image} minus the smoothed version of the \option{image}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitMakeSegim}}, \code{\link{profitMakeSegimExpand}}
}
\examples{
\dontrun{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits',
package="ProFound"))$imDat
magimage(image)
profitImBlur(image, plot=TRUE)
profitImGrad(image, plot=TRUE)
profitImDiff(image, plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ image }% use one of  RShowDoc("KEYWORDS")
