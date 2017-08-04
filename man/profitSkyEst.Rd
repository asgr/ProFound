\name{profitSkyEst}
\alias{profitSkyEst}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Local Sky Estimator
}
\description{
A high level utility to estimate the sky properties of a supplied \option{image}. This is closely related to the equivalent routines available in the LAMBDAR R package.
}
\usage{
profitSkyEst(image, objects, mask, cutlo = cuthi/2,
cuthi = sqrt(sum((dim(image)/2)^2)), skycut = 'auto', clipiters = 5, radweight = 0,
plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix; required, the image we want to analyse. The galaxy should be approximately central within this image since annuli weighting is done to avoid brighter central regions dominated by galaxy flux.
}
  \item{objects}{
Boolean matrix; optional, object mask where 1 is object and 0 is sky. If provided, this matrix *must* be the same dimensions as \option{image}.  
}
  \item{mask}{
Boolean matrix; optional, non galaxy parts of the image to mask out, where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{cutlo}{
Numeric scalar; radius where the code will start to calculate the sky annuli around the central object. Should be large enough to avoid significant object flux, i.e. a few times the flux 90 radius. Default is half of \option{cuthi}.
}
  \item{cuthi}{
Numeric scalar; radius where the code will stop calculating the sky annuli around the central object. Default is the corner edge of the \option{image}.
}
  \item{skycut}{
Numeric scalar; clipping threshold to make on the \option{image} in units of the skyRMS. The default scales the clipping to the number of pixels in the \option{image}, and will usually work reasonably.
}
  \item{clipiters}{
Numeric scalar; How many iterative clips of the sky will be made.
}
  \item{radweight}{
Numeric scalar; what radius power-law weighting should be used to bias the sky towards sky annuli nearer to the central object. \option{radweight}>0 weight the sky value more towards larger radii and \option{radweight}<0 weight the sky values towards the \option{image} centre. The default of 0 means there is no radial weightings. This becomes clear when plotting the \option{radrun} output (see Examples). Note this behaves differently to the similarly named option in LAMBDAR's sky.estimate.
}
  \item{plot}{
Logical; should a diagnostic plot be generated?
}
  \item{\dots}{
Further arguments to be passed to \code{\link{magplot}}. Only relevant is \option{plot}=TRUE.
}
}
\details{
This function is closely modelled on the sky.estimate function in the LAMBDAR package (the basic elements of which were written by ASGR). The defaults work well for data where the main objects (usually a galaxy) is centrally located in the \option{image} since the \option{cutlo} default will usually ignore contaminated central pixels. On top of this it does pretty aggressive object pixel rejection using the \option{skycut} and \option{clipiters} options.

The defaults should work reasonably well on modern survey data (see Examples), but should the solution not be ideal try modifying these parameters (in order of impact priority): \option{skycut}, \option{cutlo}, \option{radweight}, \option{clipiters}.

It is interesting to note that a better estimate of the sky RMS can be made by using the output of \code{\link{profitImDiff}} (see Examples).
}
\value{
Returns a list with 5 elements:

\item{sky}{The value of the estimated sky.}
\item{skyerr}{The estimated uncertainty in the sky level.}
\item{skyRMS}{The RMS of the sky pixels.}
\item{Nnearsky}{The number of sky annuli that have error bars encompassing the final sky.}
\item{radrun}{The output of \code{\link{magrun}} for radius versus sky pixels values.}
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profitMakeSegim}}, \code{\link{profitMakeSegimExpand}}
}
\examples{
\dontrun{
image = readFITS(system.file("extdata", 'KiDS/G266035fitim.fits',
package="ProFound"))$imDat
sky1 = profitSkyEst(image, plot=T)
image_sky = image-sky1$sky
sky2 = profitSkyEst(profitImDiff(image_sky), plot=T)

#You can check whether you are contaminated by the central objects by plotting the radrun
#object in the list (it should be flat for a well behaved sky):
sky = profitSkyEst(image, cutlo=0, plot=TRUE)
magplot(sky$radrun)
abline(h=sky$sky)

#The above shows heavy contamination by the central object without. We can either mask
#this out using the output of profitSegImWatershed, set cutlo to be larger or weight
#the sky towards outer annuli.

segim=profitSegImWatershed(image)
sky = profitSkyEst(image, mask=segim$objects, cutlo=0, plot=TRUE)
magplot(sky$radrun)
abline(h=sky$sky)

#The above is better, but not great. A more aggressive segmentation helps:

segim=profitSegImWatershed(image, skycut=0.5, sigma=3)
sky = profitSkyEst(image, mask=segim$objects, cutlo=0, plot=TRUE)
magplot(sky$radrun)
abline(h=sky$sky)

#Or weighting the sky to outer radii

segim=profitSegImWatershed(image)
sky = profitSkyEst(image, mask=segim$objects, cutlo=0, radweight=1, plot=TRUE)
magplot(sky$radrun)
abline(h=sky$sky)

#Finally we can leave the central cutlo mask turned on:

sky = profitSkyEst(image, mask=segim$objects, plot=TRUE)
magplot(sky$radrun)
abline(h=sky$sky)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ sky }% use one of  RShowDoc("KEYWORDS")
