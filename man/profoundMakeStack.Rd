\name{profoundMakeStack}
\alias{profoundMakeStack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Stack Images
}
\description{
Stacks multiple images based on their signal-to-noise.
}
\usage{
profoundMakeStack(image_list, sky_list, skyRMS_list)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image_list}{
List; each list element is a numeric matrix representing the image to be stacked.
}
  \item{sky_list}{
List; each list element is a numeric matrix representing the sky to be subtracted.
}
  \item{skyRMS_list}{
List; each list element is a numeric matrix representing the sky-RMS to weight the stack with.
}
}
\details{
The stack is actually done based on variance weighting. In pseudo code:

stack=0
stackRMS=0
for(i in 1:length(image_list)){
  stack=stack+(image_list[[i]]-sky_list[[i]])/(skyRMS_list[[i]]^2)
  sky_stack=sky_stack+(image_list[[i]]^2)
}
stack=stack*sky_stack/(length(skyRMS_list)^2)
}
\value{
A list containing:
  \item{image}{Numeric matrix; the variance weighted stacked image.}
  \item{skyRMS}{The sky RMS of the final stacked image}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundProFound}}
}
\examples{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

stack=profoundMakeStack(list(image$imDat, image$imDat, image$imDat), skyRMS_list = list(8,8,3))
mean(stack$image==image$imDat)
stack$skyRMS
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
