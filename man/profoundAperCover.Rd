\name{profoundCover}
\alias{profoundAperCover}
\alias{profoundEllipCover}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fraction of Pixel in Aperture
}
\description{
Calculates the fraction of pixels inside an aperture using recursive descent.
}
\usage{
profoundAperCover(x, y, cx, cy, rad, depth = 4L, nthreads = 1L)

profoundEllipCover(x, y, cx, cy, rad, ang, axrat, depth = 4L, nthreads = 1L) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector; the pixel x centre positions (in pixels, where pixel centres are half integer in \code{ProFound}).
}
  \item{y}{
Numeric vector; the pixel y centre positions (in pixels, where pixel centres are half integer in \code{ProFound}).
}
  \item{cx}{
Numeric scalar; the x centre of the aperture (in pixels, where pixel centres are half integer in \code{ProFound}).
}
  \item{cy}{
Numeric scalar; the y centre of the aperture (in pixels, where pixel centres are half integer in \code{ProFound}).
}
  \item{rad}{
Numeric scalar; the radius or semi-major axis of the aperture / ellipse (in pixels).
}
  \item{ang}{
Numeric scalar; orientation of the semi-major axis in degrees. This has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). 
}
  \item{axrat}{
Numeric scalar; the semi-minor/semi-major axial ratio of the ellipse (bounded between 0 - 1).
}
  \item{depth}{
Integer scalar; the recursion depth to use. Depth 4 is a good trade off of accuracy and speed (see Examples).
}
  \item{nthreads}{
Integer scalar; the number of CPU threads to run on. Only valid if \code{ProFound} is built with OpenMP support.
}
}
\details{
Re the recursion \option{depth}, higher numbers are more accurate, where you probably do not want to set it below 3. Above 5 will suffer notable slow-down. 4 is usually about the performance versus accuracy sweet spot.
}
\value{
Numeric vector the same length as \option{x} and \option{y}, which is what fraction of the pixel lies under the desired aperture / ellipse.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundAperPhot}}, \code{\link{profoundEllipseSeg}}
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

#Testing accuracy:

\dontrun{
ap_rad = 0.5 #this is smaller than a single pixel, so hard test case
ap_area = pi*ap_rad^2

ap_grid = expand.grid(seq(0,1,by=0.01), seq(0,1,by=0.01))
pix_grid = expand.grid(seq(-0.5,1.5,by=1), seq(-0.5,1.5,by=1))

test = {}
for(i in 1:dim(ap_grid)[1]){
  test = c(test, sum(profoundAperCover(pix_grid[,1], pix_grid[,2], ap_grid[i,1],
  ap_grid[i,2], ap_rad, depth=4))/ap_area)
}

#with depth 4 we do not do worse than 3.5\% accuracy:
magimage(matrix(test,101,101), magmap=F)
maghist(test)
}

#Some ellipse examples:

output = matrix(profoundEllipCover(temp[,1], temp[,2], 5, 5, 4.5, 0.5, 30, depth=4),
  10, 10)
magimage(output)

output = matrix(profoundEllipCover(temp[,1], temp[,2], 4.2, 5.1, 2.5, 0.8, 110, depth=4),
  10, 10)
magimage(output)
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
