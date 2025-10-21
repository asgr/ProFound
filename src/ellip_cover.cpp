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
                                 double rad, double ang, double axrat, int depth = 4, int nthreads = 1) {
  int n = x.size();
  NumericVector result(n);
  
  ang = ang*3.141593/180;
  
  double cos_ang = std::cos(ang);
  double sin_ang = std::sin(ang);
  double semi_min = rad * axrat;
  const double rad_plus = rad + 0.7071068;
  
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
