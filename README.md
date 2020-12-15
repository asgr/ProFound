# ProFound (R package)

<!-- badges: start -->
![R-CMD-check](https://github.com/asgr/ProFound/workflows/R-CMD-check/badge.svg)
<!-- badges: end -->

## Synopsis

Core package containing all the tools for simple and advanced source extraction. This is used for source detection, extraction and photometry, and also to create inputs for **ProFit**.

## Installation

### Getting R

First things first, you will probably want to install a recent version of **R** that lets you build packages from source. The advantage of choosing this route is you can then update bleeding edge versions directly from GitHub. If you rely on the pre-built binaries on CRAN you might be waiting much longer.

#### Mac

For Mac just get the latest binaries from the **R** project pages:

<https://cloud.r-project.org/bin/macosx/>

#### Windows

For Windows just get the latest binaries from the **R** project pages:

<https://cloud.r-project.org/bin/windows/>

#### Linux

Debian:	`sudo apt-get install r-base r-base-dev`

Fedora:	`sudo yum install R`

Suse:	More of a pain, see here <https://cloud.r-project.org/bin/linux/suse/README.html>

Ubuntu:	`sudo apt-get install r-base-dev`

All the info on binaries is here: <https://cloud.r-project.org/bin/linux/>

If you have a poorly supported version of Linux (e.g. CentOS) you will need to install **R** from source with the development flags (this bit is important). You can read more here: <https://cloud.r-project.org/sources.html>

Now you have the development version of **R** installed (hopefully) I would also suggest you get yourself **R-Studio**. It is a very popular and well maintained **R** IDE that gives you a lot of helpful shortcuts to scripting and analysing with **R**. The latest version can be grabbed from <https://www.rstudio.com/products/rstudio/> where you almost certainly want the free Desktop version.

If you wish to use the command line version of **R** on Mac (why?!) then you might need to separately install **XQuartz** and set the DISPLAY system variable via something like export DISPLAY=:0 (this is not an issue for most people however).

### Getting ProFound

Source installation from GitHub should be easy:

```R
install.packages('remotes')
remotes::install_github("asgr/ProFound")
library(ProFound)
```

A few Mac people seem to have issues with the above due to the backend used to download files. A work around seems to be to either use devtools (which I do not use as the default since it has a low more dependencies, and is tricky to install on HPCs):

```R
install.packages('devtools')
devtools::install_github("asgr/ProFound")
library(ProFound)
```

Or try the following:

```R
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true”)
remotes::install_github("asgr/ProFound")
```

I also have these options set by default in my .Rprofile, which seems to help with some of the remote install issues some people face:

```R
options(download.file.method = "libcurl")
options(repos="http://cran.rstudio.com/")
options(rpubs.upload.method = "internal")
```

If all of these do not work than the nuclear option is to download (or clone) the GitHub repo, cd to where the tar.gz file is and run in the **console** (or **Terminal** on Mac):

```console
R CMD install ProFound_X.Y.Z.tar.gz
```

where X, Y and Z should be set as appropriate for the version downloaded (check the name of the file basically).

If none of the above works then you should consider burning your computer in sacrifice to the IO Gods. Then buy a newer better computer, and try all the above steps again.

Failing all of the above, please email me for help.

#### Package Dependencies

The above should also install the required packages. If you have trouble with this you can try installing the required packages manually first and then retry the installation for **ProFound**:

```R
install.packages(c('magicaxis', 'FITSio', 'data.table')) # Required packages
install.packages(c('knitr', 'rmarkdown', 'EBImage', 'akima', 'imager', 'LaplacesDemon')) # Suggested packages
install.packages('remotes')
remotes::install_github("asgr/ProFound")
```

**ProFound** now comes with its own watershed method, and can use the **imager** package to smooth the image, this means it does not strictly require the **EBImage** package any more. This is good, because it has always been a bit of a pain to install!

To use the *old* **profoundMakeSegim** and **profoundProFound** function for image segmentation you will need to have **EBImage** installed. Since this can be a bit cumbersome on some platforms (given its dependencies) this is only listed as a suggested package. You can have a go at installing it by running:

```R
install.packages("BiocManager")
BiocManager::install("EBImage")
```

Be warned, it will probably ask you to update a whole bunch of stuff. Just say **no** and only install **EBImage** for now. If something else needs updating or is missing you will get a warning when trying to use it, and you can update those packages as neccessary. If you say yes you will probably sit there for an hour whilst it builds dozens of packages from source that you do not need to really update (or need at all).

Linux users might also need to install some non-standard graphics libraries (depending on your install). If you do not have them already, you should look to install **Cairo**, **jpeg** and **tiff** libraries (these are apparently technically not entirely free, hence not coming by default on some strictly open source Linux variants). For **Cairo** you might need to install the development version, so check this if you are having issues.

Assuming you have installed all of the packages that you need/want, you should now be able to load **ProFound** within **R** with the usual:

```R
library(ProFound)
```

## Code Example

```R
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits', package="ProFound"))

profound=profoundProFound(image, skycut=1.5, magzero=30, verbose=TRUE, plot=TRUE)
```

To find more long-form examples, including complicated fitting use-cases, please check the vignettes provided. You can browse these with:

```R
browseVignettes('ProFound')
```

Or if that does not work they are all hosted externally at <http://rpubs.com/asgr/>

## Motivation

A replacement for much of the functionality of SExtractor. This was originally developed to create good quality segmentation maps to input into **ProFit**, but it serves as a handy source extraction package in its own right, hence it now exists as its own package.

## Contributors

Code:

Aaron Robotham
Rodrigo Tobar

Testing and feedback:

Sarah Casura
Luke Davies
Simon Driver
Catherine Hale

## References

Main introduction:

Robotham A.S.G., et al., 2018, MNRAS, 476, 3137

Direct usage:

Davies L.J.M, et al., 2018, MNRAS, 480, 768
Bassett R., et al., 2019, MNRAS, 487, 2354
Hale C.L., et al., 2019, MNRAS, 487, 3971

## Resources

<https://ui.adsabs.harvard.edu/abs/2018MNRAS.476.3137R/abstract>
<https://www.dropbox.com/s/xxh0c9xlaqjyjhu/Profound-Sparks.pdf?dl=0>

## Forums

Please sign up to <http://profit.freeforums.net/> if you want to ask a question (or browse the questions asked).

## License

LGPL-3+
