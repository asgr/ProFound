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
Integer scalar; the recursion depth to use. Depth 4 is a good trade off of accuracy and speed (see Examples).
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
