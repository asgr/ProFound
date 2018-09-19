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
profoundMakeStack(image_list = NULL, sky_list = NULL, skyRMS_list = NULL, magzero_in = 0,
magzero_out = 0)
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
  \item{magzero_in}{
Numeric vector; the input mag-zero points. If length 1 then it is assumed all input frames have the same mag-zero point.
}
  \item{magzero_out}{
Numeric scalar; the output mag-zero point desired.
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

The output is explictly sky subtracted (so the sky is now 0 everywhere by definition as far as \code{\link{profoundProFound}} is concerned). The stacked sky is not returned. However, it can be computed by running \code{profoundMakeStack} again, but passing the sky list originally passed to the \option{sky_list} argument to the \option{image_list} argument instead, and not providing any input to the \option{sky_list} argument (or setting this to 0).
}
\value{
A list containing:
  \item{image}{Numeric matrix; the variance-weighted sky-subtracted stacked image.}
  \item{skyRMS}{Numeric matrix/scalar; the sky RMS image/value of the final stacked image}
  \item{magzero}{The mag-zero point of the stacked image.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundProFound}}
}
\examples{
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

stack=profoundMakeStack(list(image$imDat, image$imDat, image$imDat),
skyRMS_list = list(8,8,3))

#The new signal-to-noise weighted sky should equal sqrt(1/(1/8^2+1/8^2+1/3^2)) = 2.65

stack$skyRMS
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ stack }% use one of  RShowDoc("KEYWORDS")
