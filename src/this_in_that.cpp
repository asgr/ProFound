#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export(.vec_this_in_vec_that)]]
LogicalVector vec_this_in_vec_that(IntegerVector vec_this, IntegerVector vec_that,
                                  bool invert = false, int nthreads = 1) {

  LogicalVector ref_ID(max(vec_that) + 1, false);
  LogicalVector result(vec_this.size(), invert);
  
  #ifdef _OPENMP
    // Parallelize the main loop
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int i = 0; i < vec_that.size(); ++i) {
    if (!IntegerVector::is_na(vec_that[i])) {
      ref_ID[vec_that[i]] = true;
    }
  }
  
  #ifdef _OPENMP
    // Parallelize the main loop
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int i = 0; i < vec_this.size(); ++i) {
    int trial = vec_this[i];
    if(trial < 0 || IntegerVector::is_na(trial)){
      result[i] = NA_LOGICAL;
    }else{
      if (trial < ref_ID.size()){
        if (invert){
          result[i] = !ref_ID[trial];
        }else{
          result[i] = ref_ID[trial]; 
        }
      }
    }
  }
  return result;
}


// [[Rcpp::export(.mat_this_in_vec_that)]]
LogicalMatrix mat_this_in_vec_that(IntegerMatrix mat_this, IntegerVector vec_that,
                                  bool invert = false, int nthreads = 1) {
  
  // Create a reference vector with the size of the maximum value in vec_that + 1
  LogicalVector ref_ID(max(vec_that) + 1, invert);
  LogicalMatrix result(mat_this.nrow(), mat_this.ncol());
  
  // Populate the reference vector
  #ifdef _OPENMP
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int i = 0; i < vec_that.size(); ++i) {
    if (!IntegerVector::is_na(vec_that[i])) {
      ref_ID[vec_that[i]] = true;
    }
  }
  
#ifdef _OPENMP
  // Parallelize the main loop
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
  for (int i = 0; i < mat_this.size(); ++i) {
    int trial = mat_this[i];
    if(trial < 0 || IntegerVector::is_na(trial)){
      result[i] = NA_LOGICAL;
    }else{
      if (trial < ref_ID.size()){
        if (invert){
          result[i] = !ref_ID[trial];
        }else{
          result[i] = ref_ID[trial]; 
        }
      }
    }
  }
  return result;
}
