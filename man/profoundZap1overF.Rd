\name{profoundZap1overF}
\alias{profoundZap1overF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remove Sensor 1/f Noise
}
\description{
Removes classic row and column noise artefacts that are common the CCDs and CMOS due to 1/f read noise and other electronic effects.
}
\usage{
profoundZap1overF(image, mask = NULL, clip = c(0,0.9), scan_block = dim(image),
  sky_quan = 0.4, scan_direction = 'xy')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. Note, image NAs are treated as masked pixels.
}
  \item{mask}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{clip}{
Numeric scalar; the quantile clip level to apply to the internal scan blocks to remove bright real sources.
}
  \item{scan_block}{
Integer vector; the scan blocks to use for x/rows (x direction) and y/columns (y direction). Generally larger means less aggressive feature removal. These must be integer divisors of dim(\option{numeric}).
}
  \item{sky_quan}{
Numeric scalar; the quantile level to analyse to remove 1/f features. Lower is generally better, but too low is bad because you can hit negative pixel artefacts too.
}
  \item{scan_direction}{
Character scalar; how to scan image. Either x then y (default), y the x, just x or just y.
}
}
\details{
If \option{scan_direction} = 'xy' the code first scans blocks of the image in rows of x of size \option{scan_block}[1]. It clips out masked pixels, pixels below \option{clip}[1], and pixels above the quantile cut implied by \option{clip}[2]. It then detects the quantile defined by \option{sky_quan}, and the median of all blocks analysed is removed before the row wise 1/f map is reconstructed (i.e. the sky pedestal will remain intact, and only variations around this are removed).

The code next does a similar operation by in columns of y of size \option{scan_block}[2]. Note this is done with the x/row wise 1/f map subtracted already.

If \option{scan_first} = 'yy' the order of row and column wise operations is reversed. Generally you want to analyse the scan direction with strongest noise features first if doing both.
}
\value{

\item{image_fix}{A 1/f zapped version of the input \option{image}.}
\item{row_map}{The row wise (x block) 1/f map.}
\item{col_map}{The column wise (y block) 1/f map.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundMakeSkyGrid}}
}
\examples{
image = readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))$imDat
temp_fix = profoundZap1overF(image)
magimage(image)
magimage(temp_fix$image_fix)
magimage(temp_fix$row_map)
magimage(temp_fix$col_map)
}
