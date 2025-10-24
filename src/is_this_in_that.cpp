#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export(.vec_is_this_in_that)]]
LogicalVector vec_is_this_in_that(IntegerVector is_this, IntegerVector in_that, int nthreads = 1) {

  LogicalVector ref_ID(max(in_that) + 1, false);
  LogicalVector result(is_this.size(), false);
  
  #ifdef _OPENMP
    // Parallelize the main loop
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int i = 0; i < in_that.size(); ++i) {
    if (!IntegerVector::is_na(in_that[i])) {
      ref_ID[in_that[i]] = true;
    }
  }
  
  #ifdef _OPENMP
    // Parallelize the main loop
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int i = 0; i < is_this.size(); ++i) {
    int trial = is_this[i];
    if(trial < 0 || IntegerVector::is_na(trial)){
      result[i] = NA_LOGICAL;
    }else{
      if (trial < ref_ID.size()){
        result[i] = ref_ID[trial];
      }
    }
  }
  return(result);
}


// [[Rcpp::export(.mat_is_this_in_that)]]
LogicalMatrix mat_is_this_in_that(IntegerMatrix is_this, IntegerVector in_that, int nthreads = 1) {
  
  // Create a reference vector with the size of the maximum value in in_that + 1
  LogicalVector ref_ID(max(in_that) + 1, false);
  LogicalMatrix result(is_this.nrow(), is_this.ncol());
  
  // Populate the reference vector
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int i = 0; i < in_that.size(); ++i) {
    if (!IntegerVector::is_na(in_that[i])) {
      ref_ID[in_that[i]] = true;
    }
  }
  
  // Check if elements of is_this are in the reference vector
  #ifdef _OPENMP
  #pragma omp parallel for collapse(2) schedule(static) num_threads(nthreads)
  #endif
  for (int i = 0; i < is_this.nrow(); ++i) {
    for (int j = 0; j < is_this.ncol(); ++j) {
      int trial = is_this(i, j);
      if(trial < 0 || IntegerVector::is_na(trial)){
        result(i, j) = NA_LOGICAL;
      }else{
        if (trial < ref_ID.size()) {
          result(i, j) = ref_ID[trial];
        }
      }
    }
  }
  return result;
}
