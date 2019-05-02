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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (image, pixscale_old = 1, pixscale_new = 1, type = "bicubic", 
    fluxscale = "image") 
{
    xseq = 1:dim(image)[1] - dim(image)[1]/2 - 0.5
    yseq = 1:dim(image)[2] - dim(image)[2]/2 - 0.5
    relscale = pixscale_new/pixscale_old
    xout = seq(xseq[1], xseq[length(xseq)], by = relscale)
    yout = seq(yseq[1], yseq[length(yseq)], by = relscale)
    bigrid = expand.grid(xout, yout)
    output = matrix(0, length(xout), length(yout))
    if (type == "bilinear") {
        output[] = .interp.2d(bigrid[, 1], bigrid[, 2], list(x = xseq, 
            y = yseq, z = image))
    }
    else if (type == "bicubic") {
        output[] = akima::bicubic(xseq, yseq, image, bigrid[, 
            1], bigrid[, 2])$z
    }
    else {
        stop("type must be one of bilinear / bicubic !")
    }
    if (fluxscale == "image") {
        output = output * sum(image)/sum(output)
    }
    else if (fluxscale == "pixscale") {
        output = output * relscale^2
    }
    else if (fluxscale == "norm") {
        output = output/sum(output)
    }
    else {
        stop("fluxscale must be one of image / pixscale / norm !")
    }
    return(invisible(output))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
