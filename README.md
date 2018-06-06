# ProFound (R package)

## Synopsis

Core package containing all the tools for simple and advanced source extraction. This is used to create inputs for **ProFit**, or for source detection, extraction and photometry in its own right.

## Installation

### Getting R

Firs things first, you will probably want to install a recent version of R that lets you build packages from source. The advantage of choosing this route is you can then update bleeding edge versions directly from GitHub. If you rely on the pre-build binaries on CRAN you might be waiting much longer.

Debian:	`sudo apt-get install r-base r-base-dev`

Fedora:	`sudo yum install R`

Suse:	More of a pain, see here <https://cloud.r-project.org/bin/linux/suse/README.html>

Ubuntu:	`sudo apt-get install r-base-dev`

All the info on binaries is here: <https://cloud.r-project.org/bin/linux/>

If you have a poorly supported version of Linux (e.g. CentOS) you will need to install R from source with the development flags (this bit is important). You can read more here: <https://cloud.r-project.org/sources.html>

Now you have the development version of R installed (hopefully) I would also suggest you get yourself R-Studio. It is a very popular and well maintained R IDE that gives you a lot of helpful shortcuts to scripting and analysing with R. The latest version can be grabbed from <https://www.rstudio.com/products/rstudio/> where you almost certainly want the free Desktop version.

### Getting ProFound

Source installation from GitHub should be easy:

```R
install.packages('devtools')
library(devtools)
install_github("asgr/ProFound")
library(ProFound)
```

#### Package Dependencies

The above should also install the required packages. If you have trouble with this you can try installing the required packages manually first and then retry the installation for **ProFound**:

```R
install.packages(c('magicaxis', 'FITSio', 'data.table')) # Required packages
install.packages(c('knitr', 'rmarkdown', 'EBImage', 'akima', 'imager', 'LaplacesDemon')) # Suggested packages
install.packages('devtools')
library(devtools)
install_github("asgr/ProFound")
```

To use the **profoundMakeSegim** and **profoundProFound** function for image segmentation you will need to have **EBImage** installed. Since this can be a bit cumbersome on some platforms (given its dependencies) this is only listed as a suggested package. You can have a go at installing it by running:

```R
source("http://bioconductor.org/biocLite.R")
biocLite("EBImage")
```

Be warned, it will probably ask you to update a whole bunch of stuff. Just say **no** and only install **EBImage** for now. If something else needs updating or is missing you will get a warning when trying to use it, and you can update those packages as neccessary. If you say yes you will probably sit there for an hour whilst it builds dozens of packages from source that you do not need to really update (or need at all).

Linux users might also need to install some non-standard graphics libraries (depending on your install). If you do not have them already, you should look to install **Cairo**, **jpeg** and **tiff** libraries (these are apparently technically not entirely free, hence not coming by default on some strictly open source Linux variants). For **Cairo** you might need to install the development version, so check this if you are having issues.

Assuming this has all installed successfully, you should now be able to load **ProFound** within R with the usual:

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

A.S.G. Robotham

## Reference

Robotham A.S.G., et al., 2018, MNRAS, 476, 3137

## Forums

Please sign up to <http://profit.freeforums.net/> if you want to ask a question (or browse the questions asked).

## License

GPL-3+
