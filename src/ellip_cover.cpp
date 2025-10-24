#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// Function to check if a point is inside an ellipse
double in_ellip(double delta_x, double delta_y, double semi_maj, double semi_min, double cos_ang, double sin_ang) {
  // Rotate the point by the negative of the ellipse's rotation angle
  //double mod_x = (delta_x * sin_ang - delta_y * cos_ang) / semi_min;
  //double mod_y = (delta_x * cos_ang + delta_y * sin_ang) / semi_maj;
  
  double mod_x = (delta_x * cos_ang + delta_y * sin_ang) / semi_min;
  double mod_y = (-delta_x * sin_ang + delta_y * cos_ang) / semi_maj;
  
  
  // Check if the point is inside the ellipse
  return (mod_x * mod_x) + (mod_y * mod_y) <= 1.0 ? 1.0 : 0.0;
}

// Recursive function to determine fractional pixel coverage
double pixelCoverEllip(double delta_x, double delta_y, double semi_maj, double semi_min,
                      double cos_ang, double sin_ang, int depth) {
  if (depth == 0) {
    return in_ellip(delta_x, delta_y, semi_maj, semi_min, cos_ang, sin_ang);
  }
  
  double quarter = 0.25 / (1 << (depth - 1)); // (1 << (depth - 1)) is equivalent to pow(2, depth - 1), but faster
  double coverage = 0.0;
  
  double delta_x_min  = delta_x - quarter;
  double delta_x_plus = delta_x + quarter;
  double delta_y_min  = delta_y - quarter;
  double delta_y_plus = delta_y + quarter;
  
  coverage += pixelCoverEllip(delta_x_min, delta_y_min, semi_maj, semi_min, cos_ang, sin_ang, depth - 1);
  coverage += pixelCoverEllip(delta_x_min, delta_y_plus, semi_maj, semi_min, cos_ang, sin_ang, depth - 1);
  coverage += pixelCoverEllip(delta_x_plus, delta_y_min, semi_maj, semi_min, cos_ang, sin_ang, depth - 1);
  coverage += pixelCoverEllip(delta_x_plus, delta_y_plus, semi_maj, semi_min, cos_ang, sin_ang, depth - 1);
  
  return coverage / 4.0;
}

// [[Rcpp::export]]
NumericVector profoundEllipCover(NumericVector x, NumericVector y, double cx, double cy,
                                 double rad, double ang = 0, double axrat = 1, int depth = 4, int nthreads = 1) {
  int n = x.size();
  NumericVector result(n);
  
  double semi_min = rad * axrat;
  const double rad_plus = rad + 0.7071068;
  
  ang = ang*3.141593/180;
  double cos_ang = std::cos(ang);
  double sin_ang = std::sin(ang);
  
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
        result[i] = pixelCoverEllip(delta_x, delta_y, rad, semi_min, cos_ang, sin_ang, depth);
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
NumericMatrix profoundEllipWeight(NumericMatrix image,
                                  NumericVector cx,
                                  NumericVector cy,
                                  NumericVector rad,
                                  NumericVector ang = NumericVector::create(0),
                                  NumericVector axrat = NumericVector::create(1),
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
  
  if(ang.size() == 1){
    ang = NumericVector(n, ang[0]);
  }
  
  if(ang.size() != n){
    stop("Length of cx not equal to ang!");
  }
  
  if(axrat.size() == 1){
    axrat = NumericVector(n, axrat[0]);
  }
  
  if(axrat.size() != n){
    stop("Length of cx not equal to axrat!");
  }
  
#ifdef _OPENMP
  // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
  for (int k = 0; k < n; ++k) {
    
    double cx_loc = cx[k] - 0.5;
    double cy_loc = cy[k] - 0.5;
    double rad_loc = rad[k];
    double ang_loc = ang[k];
    double axrat_loc = axrat[k];
    
    const double semi_min = rad_loc * axrat_loc;
    const double rad_plus = rad_loc + 0.7071068;
    
    ang_loc = ang_loc*3.141593/180;
    const double cos_ang = std::cos(ang_loc);
    const double sin_ang = std::sin(ang_loc);
    
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
            weight(i,j) += pixelCoverEllip(delta_x, delta_y, rad_loc, semi_min, cos_ang, sin_ang, depth);
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
NumericVector profoundEllipFlux(NumericMatrix image,
                         NumericVector cx,
                         NumericVector cy,
                         NumericVector rad,
                         NumericVector ang = NumericVector::create(0),
                         NumericVector axrat = NumericVector::create(1),
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
  
  if(ang.size() == 1){
    ang = NumericVector(n, ang[0]);
  }
  
  if(ang.size() != n){
    stop("Length of cx not equal to ang!");
  }
  
  if(axrat.size() == 1){
    axrat = NumericVector(n, axrat[0]);
  }
  
  if(axrat.size() != n){
    stop("Length of cx not equal to axrat!");
  }
  
  NumericMatrix weight;
  
  if(deblend){
    weight = profoundEllipWeight(image, cx, cy, rad, ang, axrat, depth, nthreads);
  }
  
  #ifdef _OPENMP
    // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int k = 0; k < n; ++k) {
    
    double cx_loc = cx[k] - 0.5;
    double cy_loc = cy[k] - 0.5;
    double rad_loc = rad[k];
    double ang_loc = ang[k];
    double axrat_loc = axrat[k];
    
    const double semi_min = rad_loc * axrat_loc;
    const double rad_plus = rad_loc + 0.7071068;
    
    ang_loc = ang_loc*3.141593/180;
    const double cos_ang = std::cos(ang_loc);
    const double sin_ang = std::sin(ang_loc);
    
    int start_row = std::max(0.0, floor(cx_loc - rad_plus));
    int end_row = std::min(nrow - 1.0, ceil(cx_loc + rad_plus));
    int start_col = std::max(0.0, floor(cy_loc - rad_plus));
    int end_col = std::min(ncol - 1.0, ceil(cy_loc + rad_plus));
    
    double sum = 0.0;
    
    for (int i = start_row; i <= end_row; ++i) {
      for (int j = start_col; j <= end_col; ++j) {
        const double delta_x = i - cx_loc;
        if (std::abs(delta_x) < rad_plus) {
          const double delta_y = j - cy_loc;
          if (std::abs(delta_y) < rad_plus) {
            if(deblend){
              if(weight(i,j) > 1){
                sum += image(i, j)*pixelCoverEllip(delta_x, delta_y, rad_loc, semi_min, cos_ang, sin_ang, depth)/weight(i,j);
              }else{
                sum += image(i, j)*pixelCoverEllip(delta_x, delta_y, rad_loc, semi_min, cos_ang, sin_ang, depth);
              }
            }else{
              sum += image(i, j)*pixelCoverEllip(delta_x, delta_y, rad_loc, semi_min, cos_ang, sin_ang, depth);
            }
          }
        }
      }
    }
    result[k] = sum;
  }
  return result;
}
