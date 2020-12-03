\name{profoundResample}
\alias{profoundResample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Resample Images
}
\description{
A utility function to resample input PSFs to different pixel scales on a consistent flux conserving manner.
}
\usage{
profoundResample(image, pixscale_old = 1, pixscale_new = 1, type = "bicubic",
fluxscale = "image", recentre = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to resample.
}
  \item{pixscale_old}{
Numeric scalar; the current (old) pixel scale of the supplied \option{image}.
}
  \item{pixscale_new}{
Numeric scalar; the target (new) pixel scale of the desired output.
}
  \item{type}{
Character scalar; either "bilinear" for bilinear interpolation or "bicubic" for bicubic interpolation (default, requires \code{akima} package).
}
  \item{fluxscale}{
Character scalar; how the output image should be scaled. Either 'image' (the sum of output \option{image} will exactly qual the sum of the input \option{image}), 'pixscale' (the sum of the output pixels is scaled by the ratio of pixel scales squared, which means the output \option{image} will approximately equal the sum of the input \option{image}), or 'norm' (the output \option{immge} will sum to equal exactly 1).
}
  \item{recentre}{
Logical; should the final \option{image} have its peak flux value in the centre of the image? This is useful for resampling PSFs, where we do not want small centring errors. When the input \option{image} is highly off centre, or not a simple image of a PSF, then this option might cause image artefacts.
}
}
\details{
Mostly used for resampling PSFs to different pixel scales.
}
\value{
Numeric matrix; the desired resampled \option{image}. This will have roughly dim(\option{image})[1]*\option{pixscale_old}/\option{pixscale_new} by dim(\option{image})[2]*\option{pixscale_old}/\option{pixscale_new} pixels. The interpolation 
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundFluxDeblend}}
}
\examples{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))$imDat
star = image[149 + -25:25, 74 + -25:25]
magimage(star)
magimage(profoundResample(star,1,2))
magimage(profoundResample(star,1,1/2))
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\concept{ resample }% use one of  RShowDoc("KEYWORDS")
