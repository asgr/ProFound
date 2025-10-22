#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export]]
double profoundBoxFlux(NumericMatrix image, double cx, double cy, double size, int nthreads = 1) {
  int nrow = image.nrow();
  int ncol = image.ncol();
  
  cx -= 0.5;
  cy -= 0.5;
  
  double half_size = size / 2;
  int start_row = std::max(0.0, ceil(cx - half_size));
  int end_row = std::min(nrow - 1.0, floor(cx + half_size));
  int start_col = std::max(0.0, ceil(cy - half_size));
  int end_col = std::min(ncol - 1.0, floor(cy + half_size));
    
  double sum = 0.0;
  
  #ifdef _OPENMP
    // Parallelize the main loop
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int i = start_row; i <= end_row; ++i) {
    for (int j = start_col; j <= end_col; ++j) {
      if(!NumericMatrix::is_na(image(i, j))){
        sum += image(i, j);
      }
    }
  }
  
  return sum;
}
