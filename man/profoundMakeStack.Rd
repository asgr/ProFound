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
magzero_out = 0, masking = 'and')
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
  \item{masking}{
Character scalar; what to do with masked pixels (NAs in the \option{image}). If 'or' a pixel is masked if *any* of the images being stacked have a masked pixel (NA in the \option{image_list}) at that location, if 'and' then a pixel is masked if *all* of the images being stacked have a masked pixel at that location.
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
  \item{image}{Numeric matrix; the variance-weighted sky-subtracted stacked image. Masked pixels are NA.}
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
im1 = im2 = im3 =readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
  package="ProFound"))$imDat

stack=profoundMakeStack(list(im1, im2, im3),
skyRMS_list = list(8,8,3))

#The new signal-to-noise weighted sky should equal sqrt(1/(1/8^2+1/8^2+1/3^2)) = 2.65

stack$skyRMS

# masking logic, here we have a wedding cake of masked regions:

im1[100:200,100:200]=NA; im2[120:180,120:180]=NA; im3[140:160,140:160]=NA

# masking='or' will conservatively mask any pixel that is masked in the stack:

magimage(profoundMakeStack(list(im1, im2, im3), masking='or')$image)

# masking='and' will optimistically only mask pixels masked in all stacked images:

magimage(profoundMakeStack(list(im1, im2, im3), masking='and')$image)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ stack }% use one of  RShowDoc("KEYWORDS")
