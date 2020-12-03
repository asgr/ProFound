\name{profoundShareFlux}
\alias{profoundShareFlux}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Redistribute Fluxes
}
\description{
Redistributes fluxes from a lower resolution segmentation map using the \option{sharemat} provided by \code{\link{profoundSegimShare}}.
}
\usage{
profoundShareFlux(segstats = NULL, sharemat = NULL, weights = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segstats}{
Data.frame, the \code{\link{profoundProFound}} \option{segstats} output that has been created using the segmentation output by \code{\link{profoundSegimShare}}.
}
  \item{sharemat}{
Numeric matrix, the \option{sharemat} output that has been created by \code{\link{profoundSegimShare}}.
}
  \item{weights}{
Numeric vector; additional weights to modify fluxes by. If provided must be the same length as dim(sharemat)[2] (i.e. the original list of segment properties, not the input \option{segstats}).
}
}
\details{
This is a high level utility to extract some useful catalogue properties for objects that have disappeared from the segmentation map when degrading the resolution. Since in this case flux is conserved, we have a few option to decide how to redistribute the flux.
}
\value{
Data.frame with minimal columns. Many properties cannot be redistributed, so we ignore those here. The few which can are flux, flux_err, mag, mag_err, N50, N90, N100 (see \code{\link{profoundSegimStats}} for a discussion of these properties).
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundSegimShare}}, \code{\link{profoundSegimWarp}}, \code{\link{profoundProFound}}
}
\examples{
#None yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segim }% use one of  RShowDoc("KEYWORDS")
