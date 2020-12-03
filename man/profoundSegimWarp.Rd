\name{profoundSegimWarp}
\alias{profoundSegimWarp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Remap Segmentation Map via Warping
}
\description{
Remaps an input segmentation map WCS Tan Gnomonic projection system to a different target WCS. This uses \code{\link{magwarp}} with sensible settings, but \code{\link{magwarp}} can be used more directly if the other lower level options are required. This interface should cover most practical use cases though. Using \code{\link{profoundProFound}} with a remapped segmentation map is likely to be more sensible than remapping image flux since it will not produce flux interpolation errors.
}
\usage{
profoundSegimWarp(segim_in = NULL, header_in = NULL, header_out = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segim_in}{
Integer matrix; required, the segmentation map we want to remap. If \option{segim_in} is a list as created by \code{readFITS}, \code{read.fits} of \code{\link{magcutoutWCS}} then the image part of the list is parsed to \option{segim_in} and the correct header part is passed to \option{header_in}.
}
  \item{header_in}{
Full FITS header in table or vector format. This should be the header WCS that matches \option{segim_in}. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the astro package; the \option{hdr} output of \code{readFITS} in the \code{FITSio} package or the \option{header} output of \code{magcutoutWCS}. If \option{header_in} is provided then key words will be taken from here as a priority. Missing header keywords are printed out and other header option arguments are used in these cases.
}
  \item{header_out}{
Full FITS header in table or vector format. This is the target WCS projection that \option{segim_in} will be mapped onto. Legal table format headers are provided by the \code{read.fitshdr} function or the \option{hdr} list output of \code{read.fits} in the astro package; the \option{hdr} output of \code{readFITS} in the \code{FITSio} package or the \option{header} output of \code{magcutoutWCS}. If \option{header_out} is provided then key words will be taken from here as a priority. Missing header keywords are printed out and other header option arguments are used in these cases.
}
}
\details{
This function uses the \option{interpolation}='nearest' and \option{doscale}=FALSE in \code{\link{magwarp}}.
}
\value{
Integer matrix; the remapped image using the target WCS.
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{magwarp}}
}
\examples{
\dontrun{
VST_r=readFITS(system.file("extdata", 'VST_r.fits', package="magicaxis"))
GALEX_NUV=readFITS(system.file("extdata", 'GALEX_NUV.fits', package="magicaxis"))

profound_KiDS=profoundProFound(VST_r, sky=0, skycut=1, sigma=2, tolerance=8, plot=TRUE)

segimFUV=profoundSegimWarp(profound_KiDS$segim, profound_KiDS$header, GALEX_NUV$hdr)

profoundSegimPlot(GALEX_NUV, segim = segimFUV)

profound_GALEX=profoundProFound(GALEX_NUV, segim=segimFUV, plot=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ WCS }% use one of  RShowDoc("KEYWORDS")
\concept{ gnomonic }
\concept{ warp }% __ONLY ONE__ keyword per line
