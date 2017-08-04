\name{profitMag2Mu}
\alias{profitMag2Mu}
\alias{profitMu2Mag}
\title{
Magnitude to Surface Brightness Conversions
}
\description{
Functions to convert total magnitudes to surface brightness and vica-versa. These are provided to allow models to be either specified by total magnitude or mean surface brightness within Re. The latter is a useful way of specifying a disk model since surface brightness does not span a huge range.
}
\usage{
profitMag2Mu(mag = 15, re = 1, axrat = 1, pixscale = 1)
profitMu2Mag(mu = 17, re = 1, axrat = 1, pixscale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mag}{
Total magnitude of the 2D Sersic profile.
}
  \item{mu}{
Mean surface brightness within Re of the 2D Sersic profile.
}
  \item{re}{
Effective radii of the 2D Sersic profile.
}
  \item{axrat}{
Axial ratio of Sersic profile defined as minor-axis/major-axis, i.e. 1 is a circle and 0 is a line.
}
  \item{pixscale}{
The pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). If set to 1, then the surface brightness is interpreted in terms of pixels, otherwise it is interpreted in terms of arcseconds^2.  
}
}
\value{
profitMag2Mu returns the mean surface brightness within Re of the 2D Sersic profile.

profitMag2Mu returns total magnitude of the 2D Sersic profile.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitMakeModel}}
}
\examples{
profitMag2Mu(mag=22, re=10, axrat=0.5)
profitMu2Mag(mu=28, re=10, axrat=0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ surface-brightness }
