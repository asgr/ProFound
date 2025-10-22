#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// Recursive function to determine fractional pixel coverage
double pixelCoverAper(double delta_x, double delta_y, double delta_2,
                     double rad_2, double rad_min_2, double rad_max_2, int depth) {
  if (depth == 0) {
    return delta_2 <= rad_2 ? 1.0 : 0.0;
  }
   
  if(delta_2 > rad_max_2){
    return 0.0;
  }
   
  if(delta_2 <= rad_min_2){
    return 1.0;
  }
  
  double quarter = 0.25 / (1 << (depth - 1)); // (1 << (depth - 1)) is equivalent to pow(2, depth - 1), but faster
  double coverage = 0.0;
  
  double delta_x_min  = delta_x - quarter;
  double delta_x_plus = delta_x + quarter;
  double delta_y_min  = delta_y - quarter;
  double delta_y_plus = delta_y + quarter;
  
  double deltax1_2 = delta_x_min * delta_x_min;
  double deltax2_2 = delta_x_plus * delta_x_plus;
  double deltay1_2 = delta_y_min * delta_y_min;
  double deltay2_2 = delta_y_plus * delta_y_plus;
  
  coverage += pixelCoverAper(delta_x_min, delta_y_min, deltax1_2 + deltay1_2, rad_2, rad_min_2, rad_max_2, depth - 1);
  coverage += pixelCoverAper(delta_x_min, delta_y_plus, deltax1_2 + deltay2_2, rad_2, rad_min_2, rad_max_2, depth - 1);
  coverage += pixelCoverAper(delta_x_plus, delta_y_min, deltax2_2 + deltay1_2, rad_2, rad_min_2, rad_max_2, depth - 1);
  coverage += pixelCoverAper(delta_x_plus, delta_y_plus, deltax2_2 + deltay2_2, rad_2, rad_min_2, rad_max_2, depth - 1);
  
  return coverage / 4.0;
}


// [[Rcpp::export]]
NumericVector profoundAperCover(NumericVector x,
                                NumericVector y,
                                double cx,
                                double cy,
                                double rad,
                                int depth = 4,
                                int nthreads = 1) {
  
  const int n = x.size();
  NumericVector result(n);
  
  const double rad_2 = rad * rad;
  const double rad_min = rad - 0.7071068;   // sqrt(0.5)
  const double rad_plus = rad + 0.7071068;
  double rad_min_2 = -1.0;
  if (rad_min > 0.0) {
    rad_min_2 = rad_min * rad_min;
  }
  const double rad_max_2 = rad_plus * rad_plus;
  
  #ifdef _OPENMP
    // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
  #pragma omp parallel for schedule(dynamic, 10) if(n > 100) num_threads(nthreads)
  #endif
  for (int i = 0; i < n; ++i) {
    // Make loop-local copies to avoid data races
    const double delta_x = x[i] - cx;
    if (std::abs(delta_x) < rad_plus) {
      const double delta_y = y[i] - cy;
      if (std::abs(delta_y) < rad_plus) {
        const double delta_2 = (delta_x * delta_x) + (delta_y * delta_y);
        result[i] = pixelCoverAper(delta_x, delta_y, delta_2,
                                   rad_2, rad_min_2, rad_max_2, depth);
      } else {
        result[i] = 0.0;
      }
    } else {
      result[i] = 0.0;
    }
  }
  
  return result;
}

// [[Rcpp::export]]
double profoundAperFlux(NumericMatrix image,
                        double cx,
                        double cy,
                        double rad,
                        int depth = 4,
                        int nthreads = 1) {
  int nrow = image.nrow();
  int ncol = image.ncol();
  
  cx -= 0.5;
  cy -= 0.5;
  
  const double rad_2 = rad * rad;
  const double rad_min = rad - 0.7071068;   // sqrt(0.5)
  const double rad_plus = rad + 0.7071068;
  double rad_min_2 = -1.0;
  if (rad_min > 0.0) {
    rad_min_2 = rad_min * rad_min;
  }
  const double rad_max_2 = rad_plus * rad_plus;
  
  int start_row = std::max(0.0, floor(cx - rad_plus));
  int end_row = std::min(nrow - 1.0, ceil(cx + rad_plus));
  int start_col = std::max(0.0, floor(cy - rad_plus));
  int end_col = std::min(ncol - 1.0, ceil(cy + rad_plus));
  
  double sum = 0.0;
  
  #ifdef _OPENMP
    // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
  #pragma omp parallel for schedule(dynamic, 10) num_threads(nthreads)
  #endif
  for (int i = start_row; i <= end_row; ++i) {
    for (int j = start_col; j <= end_col; ++j) {
      if(!NumericMatrix::is_na(image(i, j))){
        const double delta_x = i - cx;
        if (std::abs(delta_x) < rad_plus) {
          const double delta_y = j - cy;
          if (std::abs(delta_y) < rad_plus) {
            const double delta_2 = (delta_x * delta_x) + (delta_y * delta_y);
            if(!NumericMatrix::is_na(image(i, j))){
              sum += image(i, j)*pixelCoverAper(delta_x, delta_y, delta_2,
                                         rad_2, rad_min_2, rad_max_2, depth);
            }
          }
        }
      }
    }
  }
  
  return sum;
}
