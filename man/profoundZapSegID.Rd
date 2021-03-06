\name{profoundMergeSegID}
\alias{profoundMergeSegID}
\alias{profoundZapSegID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Zap Elements of segID_merge
}
\description{
This function allows for safe merging or zapping of the \option{segID_merge} structure. This is important since it can be tricky to remove elements with logic if you are not an R guru.
}
\usage{
profoundMergeSegID(segID_merge)
profoundZapSegID(segID, segID_merge)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segID}{
Integer vector; segment IDs that you want to zap in \option{segID_merge}. All list entries that contain any of \option{segID} are entirely removed, so this process is very aggressive!
}
  \item{segID_merge}{
Integer list; each list element should specify collections of segments to be merged.
}
}
\details{
Merges complex lists of desired merging to their simplest form (\code{profoundMergeSegID}). Zaps unwanted list entries easily (\code{profoundZapSegID}).
}
\value{
\code{profoundMergeSegID}: Output is the simplified version of \option{segID_merge}. This means all potentially quite complex friend-of-friends groups are brought together to form a single \option{segID_merge} list item. This means each \option{segID} will occur at most 1 time in the structure.

\code{profoundZapSegID}: Output is the zapped version of \option{segID_merge}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundSegimFix}}
}
\examples{
example=list(1:5, 11:15, 21:25, 31:35)

(profoundZapSegID(3, example)) #removes first list entry since 3 appears in that vector
(profoundZapSegID(40, example)) #does nothing, since 40 does not appear anywhere
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segments }% use one of  RShowDoc("KEYWORDS")
