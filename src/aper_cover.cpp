#include <Rcpp.h>
using namespace Rcpp;

// Recursive function to determine fractional pixel coverage
double pixelCoverAper(double delta_x, double delta_y, double delta_2,
                     double radius2, double radius2_min, double radius2_max, int depth) {
  if (depth == 0) {
    return delta_2 <= radius2 ? 1.0 : 0.0;
  }
    
  if(delta_2 <= radius2_min){
    return 1.0;
  }
  
  if(delta_2 > radius2_max){
    return 0.0;
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
  
  coverage += pixelCoverAper(delta_x_min, delta_y_min, deltax1_2 + deltay1_2, radius2, radius2_min, radius2_max, depth - 1);
  coverage += pixelCoverAper(delta_x_min, delta_y_plus, deltax1_2 + deltay2_2, radius2, radius2_min, radius2_max, depth - 1);
  coverage += pixelCoverAper(delta_x_plus, delta_y_min, deltax2_2 + deltay1_2, radius2, radius2_min, radius2_max, depth - 1);
  coverage += pixelCoverAper(delta_x_plus, delta_y_plus, deltax2_2 + deltay2_2, radius2, radius2_min, radius2_max, depth - 1);
  
  return coverage / 4.0;
}

// [[Rcpp::export]]
NumericVector profoundAperCover(NumericVector x, NumericVector y, double cx, double cy, double rad, int depth = 4) {
  int n = x.size();
  NumericVector result(n);
  
  double radius_2 = rad * rad;
  
  double radius_min_2 = -1;
  
  if(rad > 0.7071068){
    radius_min_2 = (rad - 0.7071068) * (rad - 0.7071068);
  }
  
  double radius_max_2 = (rad + 0.7071068) * (rad + 0.7071068);
  
  for (int i = 0; i < n; ++i) {
    double delta_x = x[i] - cx;
    double delta_y = y[i] - cy;
    double delta_2 = (delta_x * delta_x) + (delta_y * delta_y);
    result[i] = pixelCoverAper(delta_x, delta_y, delta_2, radius_2, radius_min_2, radius_max_2, depth);
  }
  
  return result;
}

  