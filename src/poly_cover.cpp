#include <Rcpp.h>
#include <vector>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export(".point_in_polygon_cpp_short")]]
double in_poly(double testx, double testy, NumericVector poly_x, NumericVector poly_y)
{
  int i, j;
  double temp = 0.0;
  
  for (i = 0, j = poly_x.size() - 1; i < poly_x.size(); j = i++) {
    if ( ((poly_y[i] > testy) != (poly_y[j] > testy)) &&
         (testx < (poly_x[j] - poly_x[i]) * (testy - poly_y[i]) / (poly_y[j] - poly_y[i]) + poly_x[i]) )
      temp = !temp;
  }
  
  return temp;
}

// Recursive function to determine fractional pixel coverage
double pixelCoverPoly(double x, double y, NumericVector poly_x, NumericVector poly_y, int depth) {
  if (depth == 0) {
    return in_poly(x, y, poly_x, poly_y);
  }
  
  double quarter = 0.25 / (1 << (depth - 1)); // (1 << (depth - 1)) is equivalent to pow(2, depth - 1), but faster
  double coverage = 0.0;
  
  double x_min  = x - quarter;
  double x_plus = x + quarter;
  double y_min  = y - quarter;
  double y_plus = y + quarter;
  
  coverage += pixelCoverPoly(x_min, y_min, poly_x, poly_y, depth - 1);
  coverage += pixelCoverPoly(x_min, y_plus, poly_x, poly_y, depth - 1);
  coverage += pixelCoverPoly(x_plus, y_min, poly_x, poly_y, depth - 1);
  coverage += pixelCoverPoly(x_plus, y_plus, poly_x, poly_y, depth - 1);
  
  return coverage / 4.0;
}

// [[Rcpp::export]]
NumericVector profoundPolyCover(NumericVector x,
                                NumericVector y,
                                NumericVector poly_x,
                                NumericVector poly_y,
                                int depth = 4,
                                int nthreads = 1) {
  
  const int n = x.size();
  NumericVector result(n);
  
  double poly_x_min = min(poly_x) - 0.5;
  double poly_x_max = max(poly_x) + 0.5;
  double poly_y_min = min(poly_y) - 0.5;
  double poly_y_max = max(poly_y) + 0.5;
  
  #ifdef _OPENMP
    // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
  #pragma omp parallel for schedule(dynamic, 10) if(n > 100) num_threads(nthreads)
  #endif
  for (int i = 0; i < n; ++i) {
    // Make loop-local copies to avoid data races
    if (x[i] > poly_x_min && x[i] < poly_x_max) {
      if (y[i] > poly_y_min && y[i] < poly_y_max) {
        result[i] = pixelCoverPoly(x[i], y[i], poly_x, poly_y, depth);
      } else {
        result[i] = 0.0;
      }
    } else {
      result[i] = 0.0;
    }
  }
  
  return result;
}
