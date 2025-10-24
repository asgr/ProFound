\name{this_in_that}
\alias{this_in_that}
\alias{\%fin\%}
\alias{\%nin\%}
\alias{this_in_that}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Is This in That?
}
\description{
Efficient routines to test whether integer elements of one vector exist in another. This can be used in place of the \code{R} \code{\%in\%} and \code{\link{which}} functions.
}
\usage{
this_in_that(this, that, invert = FALSE, type = "logical", arr.ind = FALSE,
  nthreads = 1L)

this \%fin\% that

this \%nin\% that
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this}{
Integer vector or matrix; the object we want to identify if elements exist in \option{that}. Only non-negative and non-NA elements will be compared to \option{that}.
}
  \item{that}{
Integer vector; the comparison vector to test \option{this} against.
}
  \item{invert}{
Logical scalar; do we want to invert the output? If set \option{invert} = TRUE, FALSE become TRUE and TRUE becomes FALSE.
}
  \item{type}{
Character scalar; either 'logical', where the output is logical (TRUE and FALSE) or 'which' where the integer indices are the output.
}
  \item{arr.ind}{
Logical scalar; should array indices be returned when \option{this} is an array? Anything other than a single true value is treated as FALSE.
}
  \item{nthreads}{
Integer scalar; the number of CPU threads to run on. Only valid if \code{ProFound} is built with OpenMP support.
}
}
\details{
These functions are more efficient that the standard \code{R} function \code{\%in\%}.
}
\value{
\code{this_in_that} returns either a logical or integer vector, depending on the setting of \option{type}. Any elements that are NA or -ve in \option{this} will be NA in the corresponding output.

\code{\%fin\%} returns a logical vector.

\code{\%nin\%} returns a logical vector and is the inverse of \code{\%fin\%}.
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{match}}
}
\examples{
this_in_that(1:9, c(1L,3L,5L), type='logical')
1:9 \%fin\% c(1L,3L,5L)
this_in_that(1:9, c(1L,3L,5L), type='which')
this_in_that(matrix(1:9,3,3), c(1L,3L,5L), type='logical')
this_in_that(matrix(1:9,3,3), c(1L,3L,5L), type='which', arr.ind=TRUE)
}
