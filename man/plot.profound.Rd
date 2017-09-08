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
\method{plot}{profound}(x, ...)
}

\arguments{
  \item{x}{
Argument for the class dependent \code{plot.profound} function. An object of class profound as output by the \code{\link{profoundProFound}} function. This is the only structure that needs to be provided when executing plot(profound) class dependent plotting, which will use the plot.profound function.  
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
  
  \item{Top-left}{Input \option{image}}
  \item{Top-centre}{Output segmentation map}
  \item{Top-right}{Sky subracted and normalised image with segment dilation extent shown in colour}
  \item{Middle-left}{Magnitude counts histogram, scaled to counts per sqaure degree if possible}
  \item{Middle-centre}{Output calculated sky}
  \item{Middle-right}{Output calculated skyRMS}
  \item{Bottom-left}{Dilation iteration histogram}
  \item{Bottom-centre}{Output mag versus R50}
  \item{Bottom-right}{Output mag versus axrat}
}

\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundProFound}}
}
\examples{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)

plot(profound)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Diagnostic }% use one of  RShowDoc("KEYWORDS")
