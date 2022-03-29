\name{profoundInPoly}
\alias{profoundInPoly}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Polygon Coordinate Checker
}
\description{
Given a convex or concave polygon, this function will check which points are inside.
}
\usage{
profoundInPoly(x, y, poly_x, poly_y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Numeric vector; x locations to be checked. If two columns then first is taken to be \option{x} and second is taken to be \option{y}.
}
  \item{y}{
Numeric vector; y locations to be checked.
}
  \item{poly_x}{
Numeric vector; x vertices of the closing polygon to be checked. This can be convex or concave, and can be arranged clockwise or counter-clockwise. If two columns then first is taken to be \option{poly_x} and second is taken to be \option{poly_y}.
}
  \item{poly_y}{
Numeric vector; y vertices of the closing polygon to be checked. This can be convex or concave, and can be arranged clockwise or counter-clockwise.
}
}
\details{
Nothing else currently.
}
\value{
Logical vector of length \option{x}. This will be a flag specifying if the point is inside or outside the polygon provided.
}
\references{
Uses code provided on \url{https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html} by W. Randolph Franklin.
}
\author{
Aaron Robotham (and C code W. Randolph Franklin)
}

\seealso{
\code{\link{profoundSegimFix}}
}
\examples{
poly_x = c(2.616268, 1.570455, 4.761987, 8.945238, 9.87955, 8.422332, 5.879925, 7.683051,
  8.800988, 6.0963, 3.35555, 3.301456, 2.616268)
poly_y = c(5.01577, 7.188526, 8.292427, 8.152249, 4.577714, 5.874359, 6.505159, 3.789214,
  2.001947, 1.406191, 2.247258, 4.38497, 5.01577)

temp_grid = expand.grid(seq(0,10,by=0.1), seq(0,10,by=0.1))

temp_sel = profoundInPoly(x=temp_grid, poly_x=poly_x, poly_y=poly_y)

magplot(poly_x, poly_y, type='l')
points(temp_grid[temp_sel,])
}
\concept{ polygon }

