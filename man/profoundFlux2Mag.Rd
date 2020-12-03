\name{profoundFlux2Mag}
\alias{profoundFlux2Mag}
\alias{profoundMag2Flux}
\alias{profoundFlux2SB}
\alias{profoundSB2Flux}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert between fluxes and magnitudes.
}
\description{
Simple functions to concert between magnitudes and flux given a certain magnitude zero-point.
}
\usage{
profoundFlux2Mag(flux = 1, magzero = 0)
profoundMag2Flux(mag = 0, magzero = 0)
profoundFlux2SB(flux = 1, magzero = 0, pixscale = 1)
profoundSB2Flux(SB = 0, magzero = 0, pixscale = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{flux}{
Numeric scalar/vector; flux in ADUs given the \option{magzero}.
}
  \item{mag}{
Numeric scalar/vector; magnitude given the \option{magzero}. 
}
  \item{magzero}{
Numeric scalar/vector; magnitude zero point. What this implies depends on the magnitude system being used (e.g. AB or Vega).
}
  \item{SB}{
Numeric scalar/vector; surface brightness in mag/asec^2.
}
  \item{pixscale}{
Numeric scalar/vector; the pixel scale, where pixscale=asec/pix (e.g. 0.4 for SDSS). If set to 1, then the output is in terms of pixels, otherwise it is in arcseconds.
}
}
\details{
These functions are here to prevent silly mistakes, but the conversion is almost trivial.
}
\value{
  \item{profoundFlux2Mag}{Returns the magnitude, where \option{mag} = -2.5 * log10(\option{flux}) + \option{magzero})}
  \item{profoundMag2Flux}{Returns the flux, where \option{flux} = 10^(-0.4 * (\option{mag} - \option{magzero}))}
  HERE!!!
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundGainConvert}}
}
\examples{
profoundFlux2Mag(1e5, 30)
profoundMag2Flux(17.5, 30)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ mag }% use one of  RShowDoc("KEYWORDS")
\concept{ flux }% __ONLY ONE__ keyword per line
