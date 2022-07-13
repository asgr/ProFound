\name{profoundSegimEdge}
\alias{profoundSegimEdge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segmentation Map Edge
}
\description{
Just show the edge-most pixels of a segmentation map. This is useful for plotting purposes.
}
\usage{
profoundSegimEdge(segim = NULL, fill = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segim}{
Integer matrix; required, the segmentation map of interest.
}
  \item{fill}{
Scalar; optional, what to fill the empty space with (default will be integer 0L if left NULL).
}
}
\details{
Simple function to set inner segmentation map values to \option{fill}.
}
\value{
Returns the edgy segmentation map.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundSegimPlot}}, ~~~
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound = profoundProFound(image, magzero=30, verbose=TRUE)
colvec = c(NA, rainbow(max(profound$segim), end=2/3))
segim_edge = profoundSegimEdge(profound$segim)

magimage(image$imDat)
magimage(segim_edge, col=colvec, magmap=FALSE, add=TRUE)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
