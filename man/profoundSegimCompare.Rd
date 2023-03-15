\name{profoundSegimCompare}
\alias{profoundSegimCompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare the Bijective Similarity of Segim Maps
}
\description{
Determines the fraction pixel agreement for different segim maps of the same image. The idea is this allows people to experiment with extraction settings, and easily compare the impact on properties for similar sources.
}
\usage{
profoundSegimCompare(segim_1, segim_2, threshold = 0.5, cores = 1, ignoresky = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segim_1}{
Integer matrix; required, the first segmentation map of the image.
}
  \item{segim_2}{
Integer matrix; required, the second segmentation map of the image.
}
  \item{threshold}{
Numeric scalar; the bijective threshold to use. Setting this to 0.5 ensures results have to be unambiguous, but even 0.25 is effectively a certain match since it must mean no worse than 0.5x0.5 for the input fractions (at absolute worst there could be two such cases creating ambgiuity, but in reality that will basically never happen for a real image). Higher is obviously better, but setting too high might mean an unhelpfully small number of matches (since you are basically forcing yourself to only compare segments that agree really well, and ignoring the trickier cases).
}
  \item{cores}{
Integer scalar; number of cores to use.
}
  \item{ignoresky}{
When considering similarity, should disagreement about sky by ignored (TRUE) or considered (FALSE).
}
}
\details{
The function does a bijective process where all \option{segim_1} seg IDs are iterated over and table made of associated \option{segim_2} values in those pixels. The same process is then done the other way round. The end result is two max(segID) x max(segID) matrices (\option{seg1_in_seg2} and \option{seg2_in_seg1}) with all possible matches considered (in both cases the rows/x-axis are the reference \option{segim} and the cols/y-axis are the associated \option{segim}). These two matrices are then multiplies together to make a bijective product matrix (\option{seg_bij}) where bijective products above \option{threshold} are then extracted to form a best match table (\option{tab_match})
}
\value{
Output is a list with:

\item{seg_bij}{Numeric matrix; bijective product matrix of \option{seg1_in_seg2} x transpose(\option{seg2_in_seg1}). Thus the rows/x-axis correspond to seg IDs in \option{segim_1} and cols/y-axis correspond to seg IDs in \option{segim_2}.}
\item{tab_match}{Data frame; table of best matches (i.e. those that are above \option{threshold} in \option{seg_bij}). Columns are segID_1 / segID_2 / rowID_1 / rowID_2, where the segID is the actualy segment value and rowID is the position for this segment in the segstats output (if ordered in the default manner of ascending segID).}
\item{seg1_in_seg2}{Numeric matrix; dimensions max(segID) x max(segID) where position corresponds to segID and rows/x-axis are the reference \option{segim_1} and the cols/y-axis are the associated \option{segim_2}.}
\item{seg2_in_seg1}{Numeric matrix; dimensions max(segID) x max(segID) where position corresponds to segID and rows/x-axis are the reference \option{segim_2} and the cols/y-axis are the associated \option{segim_1}.}
}
\author{
Aaron Robotham
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{profoundMakeSegim}}
}
\examples{
image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))

profound = profoundProFound(image, skycut=1, pixcut=100, plot=TRUE)
profound2 = profoundProFound(image, skycut=5, plot=TRUE)

magimage(profound$segim, col=c(NA,hcl.colors(100)), magmap=FALSE)

testcomp = profoundSegimCompare(profound$segim, profound2$segim, threshold=0.7)

magimage(testcomp$seg_bij, col=c(NA,hcl.colors(100)), xlab='segID_1', ylab='segID_2')

magplot(profound$segstats[testcomp$tab_match$rowID_1,'mag'],
 profound2$segstats[testcomp$tab_match$rowID_2,'mag'],
 xlab='mag (ProFound 1)', ylab='mag (ProFound 2)')
abline(0,1)
}

\concept{ segim }
