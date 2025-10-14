#include <Rcpp.h>
using namespace Rcpp;

// Function to check if a point is inside a circle
bool isInsideCircle(double deltax_2, double deltay_2, double radius2) {
  return deltax_2 + deltay_2 <= radius2;
}

// Recursive function to determine fractional pixel coverage
double pixelCoverage(double x, double y, double cx, double cy, double deltax_2, double deltay_2,
                     double radius2, double radius2_min, double radius2_max, int depth) {
  if (depth == 0) {
    return isInsideCircle(deltax_2, deltay_2, radius2) ? 1.0 : 0.0;
  }
    
  if(isInsideCircle(deltax_2, deltay_2, radius2_min)){
    return 1.0;
  }
  
  if(!isInsideCircle(deltax_2, deltay_2, radius2_max)){
    return 0.0;
  }
  
  double quarter = 0.25 / (1 << (depth - 1));
  double coverage = 0.0;
  
  double deltax1_2 = (x - quarter - cx) * (x - quarter - cx);
  double deltax2_2 = (x + quarter - cx) * (x + quarter - cx);
  double deltay1_2 = (y - quarter - cy) * (y - quarter - cy);
  double deltay2_2 = (y + quarter - cy) * (y + quarter - cy);
  coverage += pixelCoverage(x - quarter, y - quarter, cx, cy, deltax1_2, deltay1_2, radius2, radius2_min, radius2_max, depth - 1);
  coverage += pixelCoverage(x - quarter, y + quarter, cx, cy, deltax1_2, deltay2_2, radius2, radius2_min, radius2_max, depth - 1);
  coverage += pixelCoverage(x + quarter, y - quarter, cx, cy, deltax2_2, deltay1_2, radius2, radius2_min, radius2_max, depth - 1);
  coverage += pixelCoverage(x + quarter, y + quarter, cx, cy, deltax2_2, deltay2_2, radius2, radius2_min, radius2_max, depth - 1);
  
  return coverage / 4.0;
}

// [[Rcpp::export]]
NumericVector profoundAperCover(NumericVector x, NumericVector y, double cx, double cy, double radius, int depth = 4) {
  int n = x.size();
  NumericVector result(n);
  
  double radius_2 = radius * radius;
  double radius_min_2 = (radius - 0.7071068) * (radius - 0.7071068);
  double radius_max_2 = (radius + 0.7071068) * (radius + 0.7071068);
  
  for (int i = 0; i < n; ++i) {
    double deltax_2 = (x[i] - cx) * (x[i] - cx);
    double deltay_2 = (y[i] - cy) * (y[i] - cy);
    result[i] = pixelCoverage(x[i], y[i], cx, cy,
                              deltax_2, deltay_2,
                              radius_2, radius_min_2, radius_max_2,
                              depth);
  }
  
  return result;
}

  