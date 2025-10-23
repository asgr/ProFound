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
NumericMatrix profoundAperWeight(NumericMatrix image,
                                 NumericVector cx,
                                 NumericVector cy,
                                 NumericVector rad,
                                 int depth = 4,
                                 int nthreads = 1) {
  int nrow = image.nrow();
  int ncol = image.ncol();
  
  int n = cx.size();
  NumericMatrix weight(nrow, ncol);
  
  if(cy.size() != n){
    stop("Length of cx not equal to cy!");
  }
  
  if(rad.size() == 1){
    rad = NumericVector(n, rad[0]);
  }
  
  if(rad.size() != n){
    stop("Length of cx not equal to rad!");
  }
  
#ifdef _OPENMP
  // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
  for (int k = 0; k < n; ++k) {
    
    double cx_loc = cx[k] - 0.5;
    double cy_loc = cy[k] - 0.5;
    double rad_loc = rad[k];
    
    const double rad_2 = rad_loc * rad_loc;
    const double rad_min = rad_loc - 0.7071068;   // sqrt(0.5)
    const double rad_plus = rad_loc + 0.7071068;
    double rad_min_2 = -1.0;
    if (rad_min > 0.0) {
      rad_min_2 = rad_min * rad_min;
    }
    const double rad_max_2 = rad_plus * rad_plus;
    
    int start_row = std::max(0.0, floor(cx_loc - rad_plus));
    int end_row = std::min(nrow - 1.0, ceil(cx_loc + rad_plus));
    int start_col = std::max(0.0, floor(cy_loc - rad_plus));
    int end_col = std::min(ncol - 1.0, ceil(cy_loc + rad_plus));
    
    for (int i = start_row; i <= end_row; ++i) {
      for (int j = start_col; j <= end_col; ++j) {
        const double delta_x = i - cx_loc;
        if (std::abs(delta_x) < rad_plus) {
          const double delta_y = j - cy_loc;
          if (std::abs(delta_y) < rad_plus) {
            // Rcout << wt[k] << "\n";
            const double delta_2 = (delta_x * delta_x) + (delta_y * delta_y);
            weight(i,j) += pixelCoverAper(delta_x, delta_y, delta_2,
                   rad_2, rad_min_2, rad_max_2, depth);
          }
        }
      }
    }
  }
  
  // We don't want to reduce increase weight when below 1
  // for (int i = 0; i < nrow; ++i) {
  //   for (int j = 0; j < ncol; ++j) {
  //     if (weight(i, j) > 0 && weight(i, j) < 1) {
  //       weight(i, j) = 1;
  //     }
  //   }
  // }
  // this is done Flux code side now to be safe
  
  return weight;
}

// [[Rcpp::export]]
NumericVector profoundAperFlux(
                        NumericMatrix image,
                        NumericVector cx,
                        NumericVector cy,
                        NumericVector rad,
                        bool deblend = false,
                        int depth = 4,
                        int nthreads = 1) {
  int nrow = image.nrow();
  int ncol = image.ncol();
  
  int n = cx.size();
  NumericVector result(n);
  
  if(cy.size() != n){
    stop("Length of cx not equal to cy!");
  }
  
  if(rad.size() == 1){
    rad = NumericVector(n, rad[0]);
  }
  
  if(rad.size() != n){
    stop("Length of cx not equal to rad!");
  }
  
  NumericMatrix weight;
  
  if(deblend){
    weight = profoundAperWeight(image, cx, cy, rad, depth, nthreads);
  }
  
  #ifdef _OPENMP
    // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int k = 0; k < n; ++k) {
  
    double cx_loc = cx[k] - 0.5;
    double cy_loc = cy[k] - 0.5;
    double rad_loc = rad[k];
    
    const double rad_2 = rad_loc * rad_loc;
    const double rad_min = rad_loc - 0.7071068;   // sqrt(0.5)
    const double rad_plus = rad_loc + 0.7071068;
    double rad_min_2 = -1.0;
    if (rad_min > 0.0) {
      rad_min_2 = rad_min * rad_min;
    }
    const double rad_max_2 = rad_plus * rad_plus;
    
    int start_row = std::max(0.0, floor(cx_loc - rad_plus));
    int end_row = std::min(nrow - 1.0, ceil(cx_loc + rad_plus));
    int start_col = std::max(0.0, floor(cy_loc - rad_plus));
    int end_col = std::min(ncol - 1.0, ceil(cy_loc + rad_plus));
    
    double sum = 0.0;
    
    for (int i = start_row; i <= end_row; ++i) {
      for (int j = start_col; j <= end_col; ++j) {
        if(!NumericMatrix::is_na(image(i, j))){
          const double delta_x = i - cx_loc;
          if (std::abs(delta_x) < rad_plus) {
            const double delta_y = j - cy_loc;
            if (std::abs(delta_y) < rad_plus) {
              const double delta_2 = (delta_x * delta_x) + (delta_y * delta_y);
              if(deblend){
                if(weight(i,j) > 1){
                  sum += image(i, j)*pixelCoverAper(delta_x, delta_y, delta_2,
                               rad_2, rad_min_2, rad_max_2, depth)/weight(i,j);
                }else{
                  sum += image(i, j)*pixelCoverAper(delta_x, delta_y, delta_2,
                               rad_2, rad_min_2, rad_max_2, depth);
                }
              }else{
                sum += image(i, j)*pixelCoverAper(delta_x, delta_y, delta_2,
                             rad_2, rad_min_2, rad_max_2, depth);
              }
            }
          }
        }
      }
    }
    result[k] = sum;
  }
  return result;
}
