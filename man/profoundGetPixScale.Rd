\name{profoundGetPixScale}
\alias{profoundGetPixScale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Get Pixel Scale
}
\description{
Given a \code{FITSio} of \code{astro} header, calculate the image pixel scale.
}
\usage{
profoundGetPixScale(header, CD1_1 = 1, CD1_2 = 0, CD2_1 = 0, CD2_2 = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{header}{
Full FITS header in table or vector format. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the astro package; the \option{hdr} output of \code{readFITS} in the \code{FITSio} package or the \option{header} output of \code{magcutoutWCS}. Missing header keywords are printed out and other header option arguments are used in these cases. See \code{\link{magWCSxy2radec}}.
}
  \item{CD1_1}{
FITS header CD1_1 for the Tan Gnomonic projection system. Change in RA-Tan in degrees along x-Axis.
}
  \item{CD1_2}{
FITS header CD1_2 for the Tan Gnomonic projection system. Change in RA-Tan in degrees along y-Axis.
}
  \item{CD2_1}{
FITS header CD2_1 for the Tan Gnomonic projection system. Change in Dec-Tan in degrees along x-Axis.
}
  \item{CD2_2}{
FITS header CD2_2 for the Tan Gnomonic projection system. Change in Dec-Tan in degrees along y-Axis.
}
}
\details{
In most cases users will simply provide a valid header to find the WCS, but you can enter the \option{CD} values explicitly. Calculating the pixel scale from the latter is almost trivial, but the option is there for the curious/lazy.
}
\value{
Numeric scalar; the image pixscale in asec/pixel (so typically a value of 0.1-0.5 for modern survey instruments).
}
\author{
Aaron Robotham
}

\examples{
\dontrun{
#The answer should be almost exactly 0.2 asec/pixel:

#Using FITSio package
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFound"))
profoundGetPixScale(image$hdr)
#Using astro package
image = read.fits(system.file("extdata", 'KiDS/G266035fitim.fits', package="ProFound"))
profoundGetPixScale(image$hdr[[1]])
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ pixscale }% use one of  RShowDoc("KEYWORDS")
