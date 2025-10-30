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
  
  // small movements first hence 0.5 / (1 << depth), or we might stop too soon
  // not clear why to me, but much faster to use this bit shift than to pre-compute and index into result
  // don't edit this!
  const double shift = 0.5 / (1 << depth); // (1 << depth) is equivalent to pow(2, depth), but faster
  double coverage = 0.0;
  
  const double delta_x_min  = delta_x - shift;
  const double delta_x_plus = delta_x + shift;
  const double delta_y_min  = delta_y - shift;
  const double delta_y_plus = delta_y + shift;
  
  const double deltax1_2 = delta_x_min * delta_x_min;
  const double deltax2_2 = delta_x_plus * delta_x_plus;
  const double deltay1_2 = delta_y_min * delta_y_min;
  const double deltay2_2 = delta_y_plus * delta_y_plus;
  
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
                                int depth = 3,
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
NumericMatrix profoundAperWeight(NumericVector cx,
                                 NumericVector cy,
                                 NumericVector rad,
                                 int dimx = 100,
                                 int dimy = 100,
                                 NumericVector wt = NumericVector::create(1),
                                 NumericVector rad_re = NumericVector::create(0),
                                 NumericVector nser = NumericVector::create(1),
                                 int depth = 3,
                                 int nthreads = 1) {
  
  int n = cx.size();
  NumericMatrix weight(dimx, dimy);
  
  if(cy.size() != n){
    stop("Length of cx not equal to cy!");
  }
  
  if(rad.size() == 1){
    rad = NumericVector(n, rad[0]);
  }
  
  if(rad.size() != n){
    stop("Length of cx not equal to rad!");
  }
  
  if(rad_re.size() == 1){
    rad_re = NumericVector(n, rad_re[0]);
  }
  
  if(rad_re.size() != n){
    stop("Length of cx not equal to rad_re!");
  }
  
  if(nser.size() == 1){
    nser = NumericVector(n, nser[0]);
  }
  
  if(nser.size() != n){
    stop("Length of cx not equal to nser!");
  }
  
  NumericVector bn(n);
  for (int k = 0; k < n; ++k) {
    bn[k] = 2*nser[k] - 1/3 + (4 / (405 * nser[k])); // use the bn approximation
  }
  
  NumericVector wt_use = NumericVector(n);
  
  if(wt.size() == 1){
    for (int k = 0; k < n; ++k) {
      wt_use[k] = wt[0];
    }
  }else{
    for (int k = 0; k < n; ++k) {
      if(rad_re[k] > 0){
        wt_use[k] = wt[k] / (rad_re[k] * rad_re[k]);
      }
      
      if(wt_use[k] < 0){
        wt_use[k] = 0;
      }
    }
  }
  
  if(wt_use.size() != n){
    stop("Length of cx not equal to wt!");
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
    
    const int start_row = std::max(0.0, floor(cx_loc - rad_plus));
    const int end_row = std::min(dimx - 1.0, ceil(cx_loc + rad_plus));
    const int start_col = std::max(0.0, floor(cy_loc - rad_plus));
    const int end_col = std::min(dimy - 1.0, ceil(cy_loc + rad_plus));
    
    for (int i = start_row; i <= end_row; ++i) {
      for (int j = start_col; j <= end_col; ++j) {
        const double delta_x = i - cx_loc;
        if (std::abs(delta_x) < rad_plus) {
          const double delta_y = j - cy_loc;
          if (std::abs(delta_y) < rad_plus) {
            const double delta_2 = (delta_x * delta_x) + (delta_y * delta_y);
            if(rad_re[k] == 0){
              weight(i,j) += wt_use[k] * pixelCoverAper(delta_x, delta_y, delta_2,
                     rad_2, rad_min_2, rad_max_2, depth);
            }else{
              weight(i,j) += wt_use[k] * pixelCoverAper(delta_x, delta_y, delta_2,
                     rad_2, rad_min_2, rad_max_2, depth) * exp(-bn[k]*pow(sqrt(delta_2) / rad_re[k], 1/nser[k]));
            }
          }
        }
      }
    }
  }
  return weight;
}

// [[Rcpp::export]]
NumericVector profoundAperFlux(
                        NumericMatrix image,
                        NumericVector cx,
                        NumericVector cy,
                        NumericVector rad,
                        NumericVector wt = NumericVector::create(1),
                        NumericVector rad_re = NumericVector::create(0),
                        NumericVector nser = NumericVector::create(1),
                        bool deblend = false,
                        int depth = 3,
                        int nthreads = 1) {
  int dimx = image.nrow();
  int dimy = image.ncol();
  
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
  
  if(rad_re.size() == 1){
    rad_re = NumericVector(n, rad_re[0]);
  }
  
  if(rad_re.size() != n){
    stop("Length of cx not equal to rad_re!");
  }
  
  NumericVector bn(n);
  for (int k = 0; k < n; ++k) {
    bn[k] = 2*nser[k] - 1/3 + (4 / (405 * nser[k])); // use the bn approximation
  }
  
  NumericVector wt_use = NumericVector(n);
  
  if(wt.size() == 1){
    for (int k = 0; k < n; ++k) {
      wt_use[k] = wt[0];
    }
  }else{
    for (int k = 0; k < n; ++k) {
      if(rad_re[k] > 0){
        wt_use[k] = wt[k] / (rad_re[k] * rad_re[k]);
      }
      
      if(wt_use[k] < 0){
        wt_use[k] = 0;
      }
    }
  }
  
  if(wt_use.size() != n){
    stop("Length of cx not equal to wt!");
  }
  
  NumericMatrix weight;
  
  if(deblend){
    weight = profoundAperWeight(cx, cy, rad, dimx, dimy, wt, rad_re, nser, depth, nthreads);
  }
  
  #ifdef _OPENMP
    // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int k = 0; k < n; ++k) {
  
    const double cx_loc = cx[k] - 0.5;
    const double cy_loc = cy[k] - 0.5;
    const double rad_loc = rad[k];
    
    const double rad_2 = rad_loc * rad_loc;
    const double rad_min = rad_loc - 0.7071068;   // sqrt(0.5)
    const double rad_plus = rad_loc + 0.7071068;
    double rad_min_2 = -1.0;
    if (rad_min > 0.0) {
      rad_min_2 = rad_min * rad_min;
    }
    const double rad_max_2 = rad_plus * rad_plus;
    
    const int start_row = std::max(0.0, floor(cx_loc - rad_plus));
    const int end_row = std::min(dimx - 1.0, ceil(cx_loc + rad_plus));
    const int start_col = std::max(0.0, floor(cy_loc - rad_plus));
    const int end_col = std::min(dimy - 1.0, ceil(cy_loc + rad_plus));
    
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
                if(weight(i,j) > 0){
                  // need the rad_pow scaling stuff here!!!
                  const double PC_temp = pixelCoverAper(delta_x, delta_y, delta_2,
                                                  rad_2, rad_min_2, rad_max_2, depth);
                  if(rad_re[k]== 0){
                    sum += image(i, j) * (PC_temp * PC_temp) * wt_use[k] / weight(i,j);
                  }else{
                    sum += image(i, j) * (PC_temp * PC_temp) * wt_use[k]  * exp(-bn[k]*pow(sqrt(delta_2) / rad_re[k], 1/nser[k])) / weight(i,j);
                  }
                }
              }else{
                sum += image(i, j)*pixelCoverAper(delta_x, delta_y, delta_2, rad_2,
                             rad_min_2, rad_max_2, depth);
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
