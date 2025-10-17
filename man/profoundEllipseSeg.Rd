\name{profoundEllipseSeg}
\alias{profoundEllipseSeg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create Elliptical Segment
}
\description{
Creates or adds an elliptical segment. This can be useful to allow the extraction of elliptical apertures on a target image (you can just sum the multiplication of the two).
}
\usage{
profoundEllipseSeg(xcen = dim[1]/2, ycen = dim[2]/2, rad = dim[1]/2, ang = 0, axrat = 1,
  segID = 1L, dim = c(101,101), segim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xcen}{
Numeric scalar; the x centre of the segment aperture ellipse.
}
  \item{ycen}{
Numeric scalar; the y centre of the segment aperture ellipse.
}
  \item{rad}{
Numeric scalar; major axis radius for the segment aperture ellipse.
}
  \item{ang}{
Numeric scalar; major axis angle for the segment aperture ellipse in degrees. When plotted as an R image the angle (theta) has the convention that 0= | (vertical), 45= \, 90= - (horizontal), 135= /, 180= | (vertical). Values outside the range 0 <= ang <= 180 are allowed, but these get recomputed as ang = ang.
}
  \item{axrat}{
Numeric scalar; axial ratio of the segment aperture ellipse defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.
}
  \item{segID}{
Integer scalar; the segID value to use for pixels lying inside the ellipse.  
}
  \item{dim}{
Integer vector; the dimension of the integer matrix to return (only used if \option{segim} is NULL).
}
  \item{segim}{
Integer matrix; if provided the elliptical segment will be added to the \option{segim}.
}
}
\value{
Integer matrix, with the dimensions of \option{dim}, or the target \option{image} if supplied.
}
\references{
Robotham et al, 2017, MNRAS, 466, 1513 (ProFit paper)
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundDrawEllipse}}, \code{\link{profoundEllipCover}}
}
\examples{
magimage(profoundEllipseSeg())
magimage(profoundEllipseSeg(dim=c(200,200), axrat=0.5, ang=30, rad=30))
}

\concept{ ~ellipse }
