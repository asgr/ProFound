#include <Rcpp.h>
using namespace Rcpp;

void interpolateAkimaGrid(NumericVector xseq, NumericVector yseq,
                          NumericMatrix tempmat_sky, NumericMatrix output);

void interpolateLinearGrid(NumericVector xseq, NumericVector yseq, NumericMatrix tempmat_sky,
                           NumericMatrix output);
