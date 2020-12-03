\name{profoundCatMerge}
\alias{profoundCatMerge}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Catalogue Merging Tool
}
\description{
Merges segmentation and grouped segmentation catalogues based on which groups are preferred.
}
\usage{
profoundCatMerge(segstats = NULL, groupstats = NULL, groupsegID = NULL,
groupID_merge = NULL, flag = TRUE, rowreset = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segstats}{
Data.frame, segmentation catalogue output from \option{profoundProFound}.
}
  \item{groupstats}{
Data.frame, grouped segmentation catalogue output from \option{profoundProFound}.
}
  \item{groupsegID}{
List; group information as output by \option{profoundSegimGroup} or \option{profoundProFound}. Must correspond to the supplied \option{segstats} and \option{groupstats}.
}
  \item{groupID_merge}{
Integer vector; group IDs that are preferred solutions. All segmented belonging to the corresponding group will be removed, and the new group photometry inserted instead.
}
  \item{flag}{
Logical; should an extra column be added to the end specifying the origin of the photometry (either 'seg' for the segmentation map, or 'group' for the grouped segmentation map)?  
}
  \item{rowreset}{
Logical; should the data.frame row names be reset to be 1:Nrow of the data.frame? The default leaves a trace of the group segment selection (i.e. you can see the selected row numbers from the provided \option{segstats}).
}
}
\details{
Handy tool to robustly merge catalogues based on preferred solutions.
}
\value{
Merged catalogue. This will have the same number of columns as \option{segstats}, with an additional column at the end called 'origin' that flags whether the object came from the segmentation catalogue (seg) or grouped segmentation catalogue (group).
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundSegimKeep}}
}
\examples{
\dontrun{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))
profound=profoundProFound(image, skycut=1.5, magzero=30, groupstats=TRUE, verbose=TRUE)

merge=profoundCatMerge(profound$segstats, profound$groupstats,
profound$group$groupsegID, 1)

profound$segstats[1,'mag']
merge[1,'mag'] #The merged object is brighter, as we should expect.
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segim }% use one of  RShowDoc("KEYWORDS")
