#include <Rcpp.h>
using namespace Rcpp;

// Function to check if a point is inside an ellipse
bool isInsideEllip(double delta_x, double delta_y, double semi_maj, double semi_min, double cos_ang, double sin_ang) {
  // Rotate the point by the negative of the ellipse's rotation angle
  double mod_x = (delta_x * sin_ang + delta_y * cos_ang) / semi_maj;
  double mod_y = (delta_x * cos_ang - delta_y * sin_ang) / semi_min;
  
  // Check if the point is inside the ellipse
  double ellipse_eq = (mod_x * mod_x) + (mod_y * mod_y);
  return ellipse_eq <= 1.0;
}

// Recursive function to determine fractional pixel coverage
double pixelCoverEllip(double delta_x, double delta_y, double semi_maj, double semi_min,
                      double cos_ang, double sin_ang, int depth) {
  if (depth == 0) {
    return isInsideEllip(delta_x, delta_y, semi_maj, semi_min, cos_ang, sin_ang) ? 1.0 : 0.0;
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
                                 double rad, double ang, double axrat, int depth = 4) {
  int n = x.size();
  NumericVector result(n);
  
  ang = ang*3.141593/180;
  
  double cos_ang = std::cos(-ang);
  double sin_ang = std::sin(-ang);
  double semi_min = rad * axrat;
  
  for (int i = 0; i < n; ++i) {
    double delta_x = x[i] - cx;
    double delta_y = y[i] - cy;
    result[i] = pixelCoverEllip(delta_x, delta_y, rad, semi_min, cos_ang, sin_ang, depth);
  }
  
  return result;
}
