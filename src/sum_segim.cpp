#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector profoundSegimFlux(NumericMatrix image, NumericMatrix segim, int nthreads = 1) {
  int nrow = segim.nrow();
  int ncol = segim.ncol();
  int max_seg = max(segim);
  
  NumericVector fluxes(max_seg, 0.0);
  
  #ifdef _OPENMP
    // Parallelize the main loop
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if(segim(i, j) > 0){
        if(!NumericMatrix::is_na(image(i, j))){
          fluxes(segim(i, j) - 1) += image(i, j); 
        }
      }
    }
  }
  
  return fluxes;
}
