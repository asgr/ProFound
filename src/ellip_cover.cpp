#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// Function to check if a point is inside an ellipse
double in_ellip(double delta_x, double delta_y, double semi_maj, double semi_min, double cos_ang, double sin_ang) {
  // Rotate the point by the negative of the ellipse's rotation angle
  
  const double mod_x = (delta_x * cos_ang + delta_y * sin_ang) / semi_min;
  const double mod_y = (-delta_x * sin_ang + delta_y * cos_ang) / semi_maj;
  
  
  // Check if the point is inside the ellipse
  return (mod_x * mod_x) + (mod_y * mod_y) <= 1.0 ? 1.0 : 0.0;
}

// Recursive function to determine fractional pixel coverage
double pixelCoverEllip(double delta_x, double delta_y, double x_term, double y_term,
                      double xy_term, int depth) {
  if (depth == 0) {
    // old code return in_ellip(delta_x, delta_y, x_term, y_term, xy_term);
    // new quicker results
    return x_term * (delta_x * delta_x) + y_term * (delta_y * delta_y) + xy_term * (delta_x * delta_y) <= 1.0 ? 1.0 : 0.0;

  }
  
  const double quarter = 0.5 / (1 << depth); // (1 << depth) is equivalent to pow(2, depth), but faster
  double coverage = 0.0;
  
  const double delta_x_min  = delta_x - quarter;
  const double delta_x_plus = delta_x + quarter;
  const double delta_y_min  = delta_y - quarter;
  const double delta_y_plus = delta_y + quarter;
  
  coverage += pixelCoverEllip(delta_x_min, delta_y_min, x_term, y_term, xy_term, depth - 1);
  coverage += pixelCoverEllip(delta_x_min, delta_y_plus, x_term, y_term, xy_term, depth - 1);
  coverage += pixelCoverEllip(delta_x_plus, delta_y_min, x_term, y_term, xy_term, depth - 1);
  coverage += pixelCoverEllip(delta_x_plus, delta_y_plus, x_term, y_term, xy_term, depth - 1);
  
  return coverage / 4.0;
}

// [[Rcpp::export]]
NumericVector profoundEllipCover(NumericVector x, NumericVector y, double cx, double cy,
                                 double rad, double ang = 0, double axrat = 1, int depth = 3, int nthreads = 1) {
  int n = x.size();
  NumericVector result(n);
  
  const double semi_min = rad * axrat;
  double semi_min_min = semi_min - 0.7071068;
  if(semi_min_min < 0){
    semi_min_min = 0;
  }
  const double rad_plus = rad + 0.7071068;
  
  ang = ang*3.141593/180;
  const double cos_ang = std::cos(ang);
  const double sin_ang = std::sin(ang);
  
  const double inv_semi_minor2 = 1 / (semi_min * semi_min);
  const double inv_semi_major2 = 1 / (rad * rad);
  const double sin_ang2 = sin_ang * sin_ang;
  const double cos_ang2 = cos_ang * cos_ang;
  
  const double x_term = cos_ang2 * inv_semi_minor2 + sin_ang2 * inv_semi_major2;
  const double y_term = sin_ang2 * inv_semi_minor2 + cos_ang2 * inv_semi_major2;
  const double xy_term = 2 * sin_ang * cos_ang * (inv_semi_minor2 - inv_semi_major2);
  
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
        if(delta_2 < semi_min_min * semi_min_min){
          result[i] = 1;
        }else if(delta_2 < rad_plus * rad_plus){
          result[i] = pixelCoverEllip(delta_x, delta_y, x_term, y_term, xy_term, depth);
        }
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
NumericMatrix profoundEllipWeight(NumericVector cx,
                                  NumericVector cy,
                                  NumericVector rad,
                                  NumericVector ang = NumericVector::create(0),
                                  NumericVector axrat = NumericVector::create(1),
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
        wt_use[k] = wt[k] / ((rad_re[k] * rad_re[k]) * axrat[k]);
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
    
    const double cx_loc = cx[k] - 0.5;
    const double cy_loc = cy[k] - 0.5;
    const double rad_loc = rad[k];
    const double ang_loc = ang[k] * 3.141593 / 180;
    const double axrat_loc = axrat[k];
    
    const double semi_min = rad_loc * axrat_loc;
    double semi_min_min = semi_min - 0.7071068;
    if(semi_min_min < 0){
      semi_min_min = 0;
    }
    const double rad_plus = rad_loc + 0.7071068;

    const double cos_ang = std::cos(ang_loc);
    const double sin_ang = std::sin(ang_loc);
    
    const double inv_semi_minor2 = 1 / (semi_min * semi_min);
    const double inv_semi_major2 = 1 / (rad_loc * rad_loc);
    const double sin_ang2 = sin_ang * sin_ang;
    const double cos_ang2 = cos_ang * cos_ang;
    
    const double x_term = cos_ang2 * inv_semi_minor2 + sin_ang2 * inv_semi_major2;
    const double y_term = sin_ang2 * inv_semi_minor2 + cos_ang2 * inv_semi_major2;
    const double xy_term = 2 * sin_ang * cos_ang * (inv_semi_minor2 - inv_semi_major2);
    
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
              if(delta_2 < semi_min_min * semi_min_min){
                weight(i,j) += wt_use[k];
              }else if(delta_2 < rad_plus * rad_plus){
                weight(i,j) += wt_use[k] * pixelCoverEllip(delta_x, delta_y, x_term, y_term, xy_term, depth);
              }
            }else{
              const double mod_x = (delta_x * cos_ang + delta_y * sin_ang) / axrat_loc;
              const double mod_y = (-delta_x * sin_ang + delta_y * cos_ang);
              const double mod_delta_2 = (mod_x * mod_x) + (mod_y * mod_y);
              if(delta_2 < semi_min_min * semi_min_min){
                weight(i,j) += wt_use[k] * exp(-bn[k]*pow(sqrt(mod_delta_2) / rad_re[k], 1/nser[k]));
              }else if(delta_2 < rad_plus * rad_plus){
                weight(i,j) += wt_use[k] * pixelCoverEllip(delta_x, delta_y, x_term, y_term, xy_term, depth) * exp(-bn[k]*pow(sqrt(mod_delta_2) / rad_re[k], 1/nser[k]));
              }
            }
          }
        }
      }
    }
  }
  return weight;
}

// [[Rcpp::export]]
NumericVector profoundEllipFlux(NumericMatrix image,
                         NumericVector cx,
                         NumericVector cy,
                         NumericVector rad,
                         NumericVector ang = NumericVector::create(0),
                         NumericVector axrat = NumericVector::create(1),
                         NumericVector wt = NumericVector::create(1),
                         NumericVector rad_re = NumericVector::create(0),
                         NumericVector nser = NumericVector::create(1),
                         bool deblend = false,
                         int depth = 3,
                         int iterations = 0,
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
        wt_use[k] = wt[k] / ((rad_re[k] * rad_re[k]) * axrat[k]);
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
    weight = profoundEllipWeight(cx, cy, rad, ang, axrat, dimx, dimy, wt, rad_re, nser, depth, nthreads);
  }
  
  #ifdef _OPENMP
    // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
  #pragma omp parallel for schedule(dynamic, 10) num_threads(nthreads)
  #endif
  for (int k = 0; k < n; ++k) {
    
    const double cx_loc = cx[k] - 0.5;
    const double cy_loc = cy[k] - 0.5;
    const double rad_loc = rad[k];
    const double ang_loc = ang[k] * 3.141593 / 180;
    const double axrat_loc = axrat[k];
    
    const double semi_min = rad_loc * axrat_loc;
    double semi_min_min = semi_min - 0.7071068;
    if(semi_min_min < 0){
      semi_min_min = 0;
    }
    const double rad_plus = rad_loc + 0.7071068;
    
    const double cos_ang = std::cos(ang_loc);
    const double sin_ang = std::sin(ang_loc);
    
    const double inv_semi_minor2 = 1 / (semi_min * semi_min);
    const double inv_semi_major2 = 1 / (rad_loc * rad_loc);
    const double sin_ang2 = sin_ang * sin_ang;
    const double cos_ang2 = cos_ang * cos_ang;
    
    const double x_term = cos_ang2 * inv_semi_minor2 + sin_ang2 * inv_semi_major2;
    const double y_term = sin_ang2 * inv_semi_minor2 + cos_ang2 * inv_semi_major2;
    const double xy_term = 2 * sin_ang * cos_ang * (inv_semi_minor2 - inv_semi_major2);
    
    const int start_row = std::max(0.0, floor(cx_loc - rad_plus));
    const int end_row = std::min(dimx - 1.0, ceil(cx_loc + rad_plus));
    const int start_col = std::max(0.0, floor(cy_loc - rad_plus));
    const int end_col = std::min(dimy - 1.0, ceil(cy_loc + rad_plus));
    
    double sum = 0.0;
    
    for (int i = start_row; i <= end_row; ++i) {
      for (int j = start_col; j <= end_col; ++j) {
        const double delta_x = i - cx_loc;
        if (std::abs(delta_x) < rad_plus) {
          const double delta_y = j - cy_loc;
          if (std::abs(delta_y) < rad_plus) {
            const double delta_2 = (delta_x * delta_x) + (delta_y * delta_y);
            if(deblend){
              if(weight(i,j) > 0){
                if(rad_re[k] == 0){
                  if(delta_2 < semi_min_min * semi_min_min){
                    sum += image(i, j) * wt_use[k] / weight(i,j);
                  }else if(delta_2 < rad_plus * rad_plus){
                    const double PC_temp = pixelCoverEllip(delta_x, delta_y, x_term, y_term, xy_term, depth);
                    sum += image(i, j) * (PC_temp * PC_temp) * wt_use[k] / weight(i,j);
                  }
                }else{
                  const double mod_x = (delta_x * cos_ang + delta_y * sin_ang) / axrat_loc;
                  const double mod_y = (-delta_x * sin_ang + delta_y * cos_ang);
                  const double mod_delta_2 = (mod_x * mod_x) + (mod_y * mod_y);
                  if(delta_2 < semi_min_min * semi_min_min){
                    sum += image(i, j) * wt_use[k] * exp(-bn[k]*pow(sqrt(mod_delta_2) / rad_re[k], 1/nser[k])) / weight(i,j);
                  }else if(delta_2 < rad_plus * rad_plus){
                    const double PC_temp = pixelCoverEllip(delta_x, delta_y, x_term, y_term, xy_term, depth);
                    sum += image(i, j) * (PC_temp * PC_temp) * wt_use[k] * exp(-bn[k]*pow(sqrt(mod_delta_2) / rad_re[k], 1/nser[k])) / weight(i,j);
                  }
                }
              }
            }else{
              if(delta_2 < semi_min_min * semi_min_min){
                sum += image(i, j);
              }else if(delta_2 < rad_plus * rad_plus){
                sum += image(i, j)*pixelCoverEllip(delta_x, delta_y, x_term, y_term, xy_term, depth);
              }
            }
          }
        }
      }
    }
    result[k] = sum;
  }
  
  if(deblend == false || iterations == 0){
    return result; 
  }else{
    result = profoundEllipFlux(image, cx, cy, rad, ang, axrat, result, rad_re,
      nser, deblend, depth, iterations - 1, nthreads);
    return result;
  }
}
