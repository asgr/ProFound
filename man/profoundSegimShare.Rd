\name{profoundSegimShare}
\alias{profoundSegimShare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Flux Share Matrix
}
\description{
Given a higher resolution input segmentation map, this function calculates how the original segments are split across a lower resolution segmentation map (\option{segim_warp}). The output of this will be used to run \code{\link{profoundProFound}} using a new target segmention map, and the outputs from this photometry run and the \option{sharemat} provided will then be passed to \code{\link{profoundShareFlux}} to create fluxes for sources that are not present in \option{segim_warp} but exist in the original \option{segim_in}. 
}
\usage{
profoundSegimShare(segim_in = NULL, header_in = NULL, header_out = NULL, pixcut = 1, 
weights = NULL)
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
  \item{pixcut}{
Integer scalar; the number of pixels required to identify an object on the new low resolution segmentation map.
}
  \item{weights}{
Numeric vector; additional weights to modify fluxes by. If provided must be the same length as dim(sharemat)[2] (i.e. the original list of segment properties, not the input \option{segstats}).
}
}
\value{
A list containing:

\item{segID_in }{Vector, the input segment IDs present in\option{segim_in}.}
\item{segID_warp }{Vector, the output segment IDs present in\option{segim_warp}.}
\item{segim_warp}{Integer matrix; the remapped image using the target WCS.}
\item{sharemat}{Numeric matrix; the sharing matrix which will have dimenstion length(segID_warp) x length(segID_in).}
\item{shareseg}{Numeric vector; the fraction of non-sky \option{segim_in} pixels that match the output \option{segim_warp}, so if near 1 then \option{segim_warp} segment should be quite clean of subpixel contaminants.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundShareFlux}}, \code{\link{profoundSegimWarp}}, \code{\link{profoundProFound}}
}
\examples{
#None yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segim }% use one of  RShowDoc("KEYWORDS")

