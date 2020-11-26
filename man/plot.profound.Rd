\name{plot.profound}
\alias{plot.profound}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ProFound Diagnostic Grid
}
\description{
A useful visual grid of ProFound diagnostics. This is useful for checking if something very odd has occurred when running the code.
}
\usage{
\method{plot}{profound}(x, logR50 = TRUE, dmag=0.5, hist='sky', ...)
}

\arguments{
  \item{x}{
Argument for the class dependent \code{plot.profound} function. An object of class profound as output by the \code{\link{profoundProFound}} function. This is the only structure that needs to be provided when executing plot(profound) class dependent plotting, which will use the plot.profound function.  
}
  \item{logR50}{
Logical; specifies whether the bottom-centre panel uses a logarithmic y-axis for R50 (default is TRUE).
}
  \item{dmag}{
Numeric scaler; the magnitude binning scale to use (default 0.5 to reflect the axis binning). The magnitude histograms always use 0.5 magnitude bins, but this controls the y-axis scaling to give the correct normalisation as if the specified binning was used. I.e. the raw counts are scaled by an additional factor of 2 if \option{dmag}=1 is specified.
}
  \item{hist}{
Character scalar; specifies the plot type for the bottom-left plot. Options are 'sky' (which is a sky pixel (image-sky)/skyRMS PDF using the objects_redo mask) or 'iters' (histogram of required iterations). Old default was 'iters', but now 'sky', since this is more useful in general.
}
  \item{\dots}{
Nothing to see here.
}
}
\details{
Run for the side effect of generating a grid of useful diagnostic plots.
}
\value{
Run for the side effect of generating a grid of useful diagnostic plots:
  
  \item{Top-left}{Sky subracted image \option{x$image}-\option{x$sky}, where blue is negative, yellow is 0, and red is positive}
  \item{Top-centre}{Output segmentation map \option{x$segim}}
  \item{Top-right}{Sky subracted and normalised image (\option{x$image}-\option{x$sky})/\option{x$skyRMS}, with segment dilation extent shown in colour}
  \item{Middle-left}{Magnitude (\option{x$segstats$mag}) counts histogram (max in red), scaled to counts per square degree if \option{x$header} is present}
  \item{Middle-centre}{Output \option{x$sky}, where blue is negative, yellow is 0, and red is positive}
  \item{Middle-right}{Output \option{x$skyRMS}, where dark is lower values and white larger values}
  \item{Bottom-left}{Sky pixel (\option{x$image}-\option{x$sky})/\option{x$skyRMS} PDF, or dilation iteration (\option{x$segstats$iter}) histogram (depends on \option{hist})}
  \item{Bottom-centre}{Output mag (\option{x$segstats$mag}) versus R50 (\option{x$segstats$R50})}
  \item{Bottom-right}{Output mag (\option{x$segstats$mag}) versus axrat (\option{x$segstats$axrat})}
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundProFound}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)

plot(profound)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diagnostic }% use one of  RShowDoc("KEYWORDS")
