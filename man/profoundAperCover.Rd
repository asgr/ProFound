\name{profoundAperCover}
\alias{profoundAperCover}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fraction of Pixel in Aperture
}
\description{
Calculates the fraction of pixels inside an aperture using recursive descent.
}
\usage{
profoundAperCover(x, y, cx, cy, radius, depth = 4L)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector; the pixel x centre positions (where pixel centres are half integer in \code{ProFound}).
}
  \item{y}{
Numeric vector; the pixel y centre positions (where pixel centres are half integer in \code{ProFound}).
}
  \item{cx}{
Numeric scalar; the x centre of the aperture (in pixels).
}
  \item{cy}{
Numeric scalar; the y centre of the aperture (in pixels).
}
  \item{radius}{
Numeric scalar; the radius of the aperture (in pixels).
}
  \item{depth}{
Integer scalar; the recursion depth to use.
}
}
\details{
Re the recursion \option{depth}, higher numbers are more accurate, where you probably do not want to set it below 3. Above 5 will suffer notable slow-down. 4 is usually about the performance versus accuracy sweet spot.
}
\value{
Numeric vector the same length as \option{x} and \option{y}, which is what fraction of the pixel lies under the desired aperture.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundAperPhot}}
}
\examples{
library(magicaxis)
library(plotrix)

temp = expand.grid(1:10 - 0.5, 1:10 - 0.5)

output = matrix(profoundAperCover(temp[,1], temp[,2], 5, 5, 3, depth=4), 10, 10)
magimage(output, magmap=FALSE)
draw.circle(5, 5, radius=3.2, border='red')

output = matrix(profoundAperCover(temp[,1], temp[,2], 3.2, 3.2, 1.8, depth=4), 10, 10)
magimage(output, magmap=FALSE)
draw.circle(3.2, 3.2, radius=1.8, border='red')

output = matrix(profoundAperCover(temp[,1], temp[,2], 7.5, 7.5, 5.8, depth=4), 10, 10)
magimage(output, magmap=FALSE)
draw.circle(7.5, 7.5, radius=5.8, border='red')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
