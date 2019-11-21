\name{plot.fitmagpsf}
\alias{plot.fitmagpsf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
FitMagPSF diagnostic plots
}
\description{
A simple image / model / image - model grid.
}
\usage{
\method{plot}{fitmagpsf}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class 'fitmagpsf', as returned by \code{\link{profoundFitMagPSF}}.
}
  \item{\dots}{
Passed to \code{\link{magimage}} / \code{\link{magimageWCS}}
}
}
\details{
Run for the side effect of generating a grid of useful diagnostic plots.
}
\value{
Run for the side effect of generating a grid of useful diagnostic plots:
}
\author{
Aaron Robotham
}

\seealso{
\code{\link{profoundFitMagPSF}}
}
\examples{
\dontrun{
s250_im=readFITS(system.file("extdata", 'IRdata/s250_im.fits', package="ProFound"))
s250_psf=readFITS(system.file("extdata",'IRdata/s250_psf.fits', package="ProFound"))$imDat
magzero_s250=11.68

pro_s250=profoundProFound(s250_im, pixcut=1, skycut=2, ext=1, redosky=FALSE, iters=1,
tolerance=0, sigma=0, magzero=magzero_s250)
pro_s250$segstats=pro_s250$segstats[!is.na(pro_s250$segstats$mag),]

newmag=profoundFitMagPSF(RAcen=pro_s250$segstats$RAcen, Deccen=pro_s250$segstats$Deccen,
image=s250_im, psf=s250_psf, doProFound=TRUE, findextra=TRUE, verbose=TRUE, redosky=FALSE,
magzero=magzero_s250)

plot(newmag)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ plot }% use one of  RShowDoc("KEYWORDS")
