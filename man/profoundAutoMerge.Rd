\name{profoundAutoMerge}
\alias{profoundAutoMerge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Segmentation Map Auto Merging
}
\description{
Attempts to de-fragment large flocculent galaxies that tend to be broken up, even with reasonable watershed parameters. Ideally the settings should be chosen such that you do not go too far the other way, and over-merge objects again.
}
\usage{
profoundAutoMerge(segim, segstats, spur_lim = 0.004, col_lim = NULL, Ncut = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segim}{
Integer matrix; required, the segmentation map. This should usually be the \option{segim_orig} (un-dilated) segmentation map, because that usually behaves better with this function. Dilated segmentation maps can also be proved though.
}
  \item{segstats}{
Data.frame; the \option{segstats} output of ProFound. If you want to use the \option{col_lim} option, you need to add a non-standard column representing your magnitude colour of choice. A useful colour is usually optical minus NIR, e.g. (r - Z).
}
  \item{spur_lim}{
Numeric scalar; the threshold in segstats$cenfrac / (1 + segstats$sep) below which when selecting segments as potential merging candidates. The default works well for optical and NIR images, but will require some tuning if the seeing is very far from 1 arc-second.
}
  \item{col_lim}{
Numeric vector; if provided this is the lower and upper colour limit to use when selecting segments as potential merging candidates. If this is used then the 'col' column must be supplied in \option{segim}. If the colour provided is optical - NIR then a range of between 0 and 1 is usually sensible for picking out true galaxy segments.
}
  \item{Ncut}{
Integer scalar; at least this many segments need to be touching to be provided in the final output. \option{Ncut} or 1 or 2 usually makes most sense.
}
}
\details{
This function should not be applied blindly and requires some though and tuning on images obtained from a new facility. The default works well for good seeing optical imaging from e.g. KiDS.
}
\value{
The output is exactly the same as the \option{groupsegID} from \code{\link{profoundSegimGroup}}. For reference see below.

  \item{groupsegID}{A data.frame of lists giving the segIDs of segments in each group.}

The data.frame returned by \option{groupsegID} is a slightly unusal structure to see in R, but it allows for a compact manner of storing uneven vectors of grouped segments. E.g. you might have a massive group containing 30 other segments and many groups containing a single segment. Padding a normal matrix out to accommodate the larger figure would be quite inefficient. It contains the following:

  \item{mergeID}{Merging group ID. This is different to the usual \option{groupID} since it is possible to have islands that have been identified to be re-grouped that lie within superset groups.}
  \item{segID}{An embedded list of segmentation IDs for segments in the merging group. I.e. each list element of \option{segID} is a vector (see Examples for clarity).}
  \item{Nemrge}{The total number of segments that are in the merging group.}
  \item{Npix}{The total number of pixels that are in the merging group.}
  \item{groupID}{Group ID of the superset group based on the provided \option{segim}.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundSegimGroup}}
}
\examples{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound = profoundProFound(image, magzero=30)

#Below is just to see how it works, this is clearly over-merging.

groupsegID = profoundAutoMerge(segim=profound$segim_orig, segstats=profound$segstats,
  spur_lim=3e-2)

segimfix = profoundSegimKeep(segim=profound$segim, segID_merge=groupsegID$segID)
}

\concept{ segmentation }

