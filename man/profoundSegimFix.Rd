\name{profoundSegimFix}
\alias{profoundSegimFix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interactive Segment Fixing Tool
}
\description{
This function will launch an interactive session to let users reassemble segments of groups etc. The idea is that you might have very large galaxies that are hard to keep together with \code{\link{profoundProFound}} settings that also work well at separating close sources. Using this tool allows users to quickly put things back together in an easy fashion.
}
\usage{
profoundSegimFix(image = NULL, segim = NULL, mask = NULL, sky = NULL, profound = NULL,
loc = NULL, box = 400, segID_merge = list(), col = 'magenta', pch = 4, cex = 2,
crosshair = FALSE, crosscex = 5, alpha_seg = 0.3, happy_default = TRUE,
continue_default = TRUE, open_window = TRUE, allow_seg_modify = FALSE, segID_max = NULL,
...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{image}{
Required, the image we want to display. If a numeric matrix is passed in then \code{\link{magimage}} is used to display the image. If the \option{image} has components R/G/B then then \code{\link{magimageRGB}} is used to display a colour image.
}
  \item{segim}{
Integer matrix; required, the segmentation map of the \option{image}. This matrix *must* be the same dimensions as \option{image}.
}
  \item{mask}{
Boolean matrix; optional, parts of the image to mask out (i.e. ignore), where 1 means mask out and 0 means use for analysis. If provided, this matrix *must* be the same dimensions as \option{image}.
}
  \item{sky}{
User provided estimate of the absolute sky level. Can be a scalar or a matrix matching the dimensions of \option{image} (allows values to vary per pixel). This will be subtracted off the \option{image} internally, so only provide this if the sky does need to be subtracted!
}
  \item{profound}{
List; object of class 'profound'. If this is provided then missing input arguments are taking directly from this structure. As an added convenience, you can assign the profound object directly to the \option{image} input.
}
  \item{loc}{
Numeric vector; the [x,y] location where we want to cutout the image.
}
  \item{box}{
Numeric vector; the dimensions of the box to cut out from \option{image} centred on \option{loc}.
}
  \item{segID_merge}{
Integer list; each list element should specify collections of segments to be merged. This will be updated by the fixing session, and the updated results are returned as part of the output list.
}
  \item{col}{
Colour of the overlaid click positions.
}
  \item{pch}{
Point type of the overlaid click positions.
}
  \item{cex}{
Size of the overlaid click positions.
}
  \item{crosshair}{
Logical; should a crosshair be placed in the centre of the image to guide the eye?
}
  \item{crosscex}{
Numeric scalar; size of the overlaid crosshair.
}
  \item{alpha_seg}{
Numeric scalar; the alpha transparency of the overlaid segment polygons.
}
  \item{happy_default}{
Logical; what should the default answer be to whether you are happy with the current fix (default is TRUE, i.e. yes).
}
  \item{continue_default}{
Logical; what should the default answer be to whether you what to fix more segments (default is TRUE, i.e. yes).
}
  \item{open_window}{
Logical; should a new Quartz or X11 window be opened. When looping through different image regions this should probably be FALSE. It is here because it is a very bad idea to try to do this using the R-Studio plot window- it does not properly support interactive plotting and will not show the click positions.
}
  \item{allow_seg_modify}{
Logical; should users be allowed to add and delete segments. Since this behaviour will modify the input \option{segim} it is not always desirable. To be safe, it is off (FALSE) by default. See Details.
}
  \item{segID_max}{
Integer; the maximum segment ID to assume for the input \option{segim} when adding or replacing segments (so only relevant if \option{allow_seg_modify}=TRUE), where the default segID for a new segment becomes \option{segID_max} + 1. Providing this speeds things up since it can be slow to compute for very large images (tenths of second for 10k x 10k). If set to 'auto' then it will check the input \option{segim}, if left at the default NULL it will check the available segim when it is creating a new segment. Since this will be on a segmentation image subset when specifying a box, this might mean the segID is not unique on the larger input \option{segim}, but it will be unique on the cut down one being processed internally, and the ones returned (\option{segim} and \option{segim_start}).
}
  \item{\dots}{
Additional arguments to be passed to \code{\link{magimage}} or \code{\link{magimageRGB}} as appropriate.
}
}
\details{
Basic instructions:

Follow the prompts when you run the function (see Example).

1)  Click on the segments you want to merge together into one segment and hit ESC when done.

2)  You will then be asked if you are happy. If you are then type 'y' in the command prompt and hit enter. If not type 'n' and hit enter. The former will merge the segments, the latter will undo the most recent round of fixes and let you start again.

3)  You will then be asked if you want to continute fixing. If you want to continue fixing other broken segments in the current image then type 'y' in the command prompt and hit enter. If you want to stop type 'n' and hit enter. The former will then prompt you to fix more segments with another round of clicking, the latter will stop the fixing session and return the results as they stand.

Features:

You cannot merge segments into the sky! All clicks on the sky are ignored, and they will not produce any of the multi-click results discussed below.

Currently fixed segment complexes (as per the internal version of \option{segID_merge}, with both the input fixes and the fixes made in the current session) are shown by overlaid polygons showing the original colours of the merged segments. This stops you from fixing complexes that have already been fixed.

If you see one of these coloured polygon complexes and do not like the solution you can click once in it and hit ESC in the plot window. This will tell the code to unpick the complex and return all filled segments back to their original solution. You can then try to fix it again, if needed.

If you make a number of clicks and then change your mind, click on the unwanted segment twice and make some other this will ignore the segment when merging- this allows you to fix silly click errors by clicking again.

If you click on any segment three times and hit ESC in the plot window it will tell the code that you are both definitely happy *and* that you want to stop fixing segments, ending the session instantly. This allows you to make rapid progress when doing lots of easy fixes since you do not need to then type anything at all into the command prompt (which slows down progress).

If \option{allow_seg_modify}=TRUE and you just click on a segment twice and nothing else it will flag it for having its segID changed or potentially deleted by setting to 0 (sky). It will ask what segID you want for the segment, where the default of "auto" will make a new segID one larger than the current maximum. You can specify your own choice manually, where if you use a current segID it will merge the sources. It is possible to delete the segment by setting the segID when asked to "0", which means it will become sky.

If \option{allow_seg_modify}=TRUE and you click on the sky once and hit ESC you will enter the segment creation mode. Now all clicks made will form a concave polygon, with the contained pixels becoming a new segment with a new segID (\option{segID_max} + 1, or user specified). This new segment will be registered on both the \option{segim_start} and \option{segim} outputs to ensure internal consistency. Double clicking a single segment (see above) also enters this mode, and the segment is changed in the same manner as if it were drawn manually. In both cases it is possible to delete segments contained within the polygon simply by setting the segID when asked to "0", which means it will become sky. In this manner you can delete potentially many at once by drawing round large groups of unwanted segments.
}
\value{
List containing:
\item{segim}{Numeric matrix; fixed version of the segmentation map. If a cutout has been requested this will be for just the cutout section.}
\item{segim_start}{Numeric matrix; reference starting segmentation map.}
\item{segID_merge}{List; fixed version of the list of segIDs to be merged. Usefully, this can be directly passed into \code{\link{profoundSegimKeep}} (see Examples).}
\item{segID_max}{Integer; the current maximum segment ID in the output \option{segim}. This should only be modified compared to the input \option{segim} if new segments have be drawn on the image.}

It is always the case that \option{segim_start} combined with \option{segID_merge} will produce \option{segim} via \code{\link{profoundSegimKeep}}. We provide \option{segim_start} for reference since the segment creation mode (see Details above) will also modify the base segmentation map in order to guarantee the internal consistency. This means the output \option{segim_start} and differ from the input \option{segim_start}, but only if new segments have been manually drawn.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{profoundProFound}}; \code{\link{profoundSegimKeep}}
}
\examples{
\dontrun{
GALEX_NUV=readFITS(system.file("extdata", 'GALEX_NUV.fits', package="magicaxis"))
VST_r=readFITS(system.file("extdata", 'VST_r.fits', package="magicaxis"))
VISTA_K=readFITS(system.file("extdata", 'VISTA_K.fits', package="magicaxis"))

# Warp to common WCS:

GALEX_NUV_VST=magwarp(GALEX_NUV, VST_r$hdr)$image
VISTA_K_VST=magwarp(VISTA_K, VST_r$hdr)$image

#Good but not perfect:

profound=profoundProFound(VST_r, roughpedestal=TRUE, plot=TRUE)

#Let's fix it:

fixed=profoundSegimFix(profound)

#Colour might help:

fixedRGB=profoundSegimFix(list(R=VISTA_K_VST, G=VST_r$imDat, B=GALEX_NUV_VST),
segim=profound$segim)

#Assuming you made some fixes you can also recreate the fixed segim using the output
# \option{segID_merge})(i.e. you only really need to save \option{segID_merge}):

fixedRGB2=profoundSegimKeep(segim=profound$segim, segID_merge=fixedRGB$segID_merge)

#Check things look the same:

profoundSegimPlot(image=VST_r$imDat, segim=fixedRGB$segim)
profoundSegimPlot(image=VST_r$imDat, segim=fixedRGB2)

#We can now feed this back into ProFound for our best final effort:

profound2=profoundProFound(VST_r, segim=fixedRGB$segim, plot=TRUE, SBdilate=1)

#The object at 300, 300 is the new merged object (if you have merged the main galaxy).
#On a local run this had segID=6, old mag=15.5, new mag=15.2 (i.e. flux increase of 30\%).
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segments }% use one of  RShowDoc("KEYWORDS")
