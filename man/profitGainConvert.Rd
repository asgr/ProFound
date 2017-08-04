\name{profitGainConvert}
\alias{profitGainConvert}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert gain between mag-zero points
}
\description{
Simple function to update the gain (electrons/ADU) when changing between magnitude zero points. These gains are what should be parsed to e.g. \code{\link{profitMakeSigma}}.
}
\usage{
profitGainConvert(gain = 1, magzero = 0, magzero_new = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gain}{
Numeric scalar or vector; the current gain/s in electrons/ADU.
}
  \item{magzero}{
Numeric scalar or vector; the current magnitude zero point/s.
}
  \item{magzero_new}{
Numeric scalar or vector; the new magnitude zero point/s.
}
}
\details{
A simple function that is mostly here to avoid silly conversion mistakes. The conversion is calculated as: gain*10^(-0.4*(magzero_new - magzero)), where an object magnitude can be calculated from ADU flux as -2.5*log10(flux_ADU)+magzero.
}
\value{
Numeric scalar or vector; the new gain/s.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profitMakeSigma}}, \code{\link{profitFlux2Mag}}, \code{\link{profitMag2Flux}}
}
\examples{
#For optical survey data typically images with gain~1 have a magzero~30:
profitGainConvert(1,30,0)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ gain }% use one of  RShowDoc("KEYWORDS")
\keyword{ magzero }% __ONLY ONE__ keyword per line
