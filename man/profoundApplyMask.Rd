\name{profoundApplyMask}
\alias{profoundApplyMask}
\alias{profoundMakeMask}
\alias{profoundDrawMask}
\alias{profoundCoverMask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Apply a Mask to an Image
}
\description{
Various mask related functions. profoundApplyMask allows the application of a scale-able and rotate-able mask to an image.
}
\usage{
profoundApplyMask(image = NULL, mask = 'disc', xcen = xsize/2, ycen = ysize/2,
  xsize = 101, ysize = 101, rot = 0, direction = 'backward', dim = c(101, 101))
  
profoundMakeMask(size = 101, shape = 'disc')

profoundDrawMask(image, poly = NULL, invert_mask = FALSE, mode = 'draw', type = 'pix',
  poly.col = 'red', ...)
  
profoundCoverMask(segim, mask)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Numeric matrix or Rfits_image; the image we want to apply the \option{mask} to. If provided this will be used to compute the \option{dim} also.
}
  \item{mask}{
Numeric matrix or scalar character; the input mask. If a matrix input, any values above 0 will be considered pixels that we wisk to mask. Note for centreing to make sense the provided mask should nearly always be odd dimensions (since even dimensions will have ambiguous application solutions at pixel centres). If character then an internal mask maker is called as per \code{profoundMakeMask}, where the \option{size} argument is the maximum of \option{xsize} or \option{ysize}, and the character scalar of \option{mask} is passed into the \option{shape} argument. The default therfore creates a circular mask of diameter 101.
}
  \item{segim}{
Integer matrix; required, the segmentation map of the \option{image}. This matrix *must* be the same dimensions as \option{mask}.
}
  \item{xcen}{
Numeric scalar/vector; the x position to centre the mask at (i.e. it will be centred on dim(mask)[1]/2). If a vector then multiple masks will be applied to the image.
}
  \item{ycen}{
Numeric scalar/vector; the y position to centre the mask at (i.e. it will be centred on dim(mask)[2]/2). If a vector then multiple masks will be applied to the image. If \option{xcen} is a vector and \option{ycen} is scalar then the value of \option{ycen} will be repeated as required.
}
  \item{xsize}{
Numeric scalar/vector; the x size we wish the final mask to be. Note this operation is carried out before any rotation. This can be used to (e.g.) distort an input circle to an ellipse. If a vector then multiple masks will be applied to the image. If \option{xcen} is a vector and \option{xsize} is scalar then the value of \option{xsize} will be repeated as required.
}
  \item{ysize}{
Numeric scalar/vector; the y size we wish the final mask to be. Note this operation is carried out before any rotation. This can be used to (e.g.) distort an input circle to an ellipse. If a vector then multiple masks will be applied to the image. If \option{xcen} is a vector and \option{ysize} is scalar then the value of \option{ysize} will be repeated as required.
}
  \item{rot}{
Numeric scalar/vector; anti-clockwise rotation to apply to the \option{mask}. Note this operation is carried out after any re-sizing operation. If a vector then multiple masks will be applied to the image. If \option{xcen} is a vector and \option{rot} is scalar then the value of \option{rot} will be repeated as required.
}
  \item{direction}{
Character scalar/vector; "forward" or "backward", see \code{imwarp}. For masking purposes, 'backward' is nearly always the safer option. If a vector then multiple masks will be applied to the image. If \option{xcen} is a vector and \option{direction} is scalar then the value of \option{direction} will be repeated as required.
}
  \item{dim}{
Integer vector; the dimension of the complex mask to return. This is used if \option{image} is not provided.
}
  \item{size}{
Integer scalar; the size (e.g. width/diameter) of the dilation kernel in pixels. Should be an odd number else will be rounded up to the nearest odd number.
}
  \item{shape}{
Character scalar; the shape of the dilation kernel. Options are 'box', 'disc', 'diamond', 'Gaussian', 'line'.
}
  \item{poly}{
Data.frame; contains minimally named x/y or RA/Dec polynomial mask vertices (and possibly both). Note pixel centres are half integer, e.g. the centre of the bottom-left pixel is [0.5,0.5].
}
  \item{invert_mask}{
Logical; if \option{invert_mask} is FALSE then pixels inside the polygon are masked (e.g. to remove bright stars etc), if TRUE then pixels outside the polygon are masked (e.g. you only want to consider pixels belonging to a galaxy or cluster region etc).
}
  \item{mode}{
Character scalar; how to compute the mask. Options are 'draw', where users can then click on the image to define the mask region, or 'apply' where the user supplied \option{poly} will be used to mask pixels.  
}
  \item{type}{
Character scalar; the coordinates to use. Options are 'pix' where the named x/y columns in \option{poly} are used to define the mask (ignoring RA/Dec even if present) or 'coord' where the named RA/Dec columns in \option{poly} are used to define the mask (ignoring x/y even if present). Note pixel centres are half integer, e.g. the centre of the bottom-left pixel is [0.5,0.5]. If \option{type} = 'coord' then the input \option{image} must be type Rfits_image (i.e. have a valid WCS).
}
  \item{poly.col}{
Colour of the overlaid polygon defining the mask region.
}
  \item{\dots}{
Additional arguments to pass on to the image plotting routine.
}
}
\details{
\code{profoundApplyMask} applies a masking brush to multiple locations in an image. A typical use case would be we have a masking brush for a star and we want to scale and apply this to different regions. Details:

The masked values in the target output \option{image} are simply set to NA.

When using vector arguments to apply multiple masks, \option{xcen} is used to define the number to apply. Other scalar valued inputs are repeated to match the length of \option{xcen}. In general you minimally need to provide \option{xcen} and \option{ycen} vectors, since positions to mask will usually be distributed around the target \option{image}.

\code{profoundMakeMask} creates a simple masking brush.

\code{profoundDrawMask} allows us to directly draw or check a polygon based (vector) mask.

\code{profoundCoverMask} calculates what fraction of segment pixels are covered by a given mask.
}
\value{
For \code{profoundApplyMask} a list containing:

\item{mask}{The final complex mask. Each masked region has a unique integer ID (the position within the \option{xcen} vector), where lower valued regions take precendence in overlap regions.}
\item{image}{Same as the input \option{image} (if provided), except pixels flagged as masked based on the provided information are set to NA.}

For \code{profoundMakeMask} a matrix of the minimal masking brush that can be used as an input to \code{profoundApplyMask}.

For \code{profoundDrawMask} a list containing:

\item{mask}{The final boolean mask based on the drawn or provided (via \option{poly}) mask polygon.}
\item{image}{Same as the input \option{image} (if provided), except pixels flagged as masked based on the provided information are set to NA.}
\item{poly}{Data.frame of mask vertices x/y and possibly RA/Dec (if appropriate).}

For profoundCoverMask a data.frame containing:

\item{segID}{The segID}
\item{orig_N}{The number of segment pixels before masking.}
\item{mask_N}{The number of segment pixels after masking.}
\item{frac}{The unmasked fraction, i.e. mask_N/orig_N.}

This is ordered by segID, where segID = 0 refers to sky pixels (where it might be of some interest to see what fraction of those are masked too).
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundEllipseSeg}}
}
\examples{
image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))$imDat

magimage(image)

mask_sq = matrix(1, 101, 101)

magimage(profoundApplyMask(image, mask_sq, 100.5, 100.5)$image)
magimage(profoundApplyMask(image, mask_sq, 100.5, 100.5, xsize=51)$image)
magimage(profoundApplyMask(image, mask_sq, 100.5, 100.5, xsize=51, rot=20)$image)
magimage(profoundApplyMask(image, mask_sq, 100.5, 100.5, xsize=51, rot=90)$image)

mask_circ = profoundEllipseSeg(dim = c(101,101))
magimage(profoundApplyMask(image, mask_circ, 100.5, 100.5, xsize=51, rot=20)$image)

image = Rfits_read_image(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',
package="ProFound"))

# we will pretend we have already drawn a mask and got our vertices
poly = data.frame(RA = c(352.2922, 352.2958, 352.2954, 352.2904, 352.2870, 352.2872),
         Dec = c(-31.81325, -31.81332, -31.81657, -31.81794, -31.81763, -31.81461)
       )

mask_coord = profoundDrawMask(image, poly=poly, mode='apply', type='coord')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
