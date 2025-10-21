#include <Rcpp.h>
#include <vector>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

static double in_poly_internal(double testx, double testy, const NumericVector &vertx, const NumericVector &verty)
{
  int i, j;
  double temp = 0.0;
  
  for (i = 0, j = vertx.size() - 1; i < vertx.size(); j = i++) {
    if ( ((verty[i] > testy) != (verty[j] > testy)) &&
         (testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]) )
      temp = !temp;
  }
  
  return temp;
}

// [[Rcpp::export(".point_in_polygon_cpp_short")]]
double in_poly(double testx, double testy, NumericVector vertx, NumericVector verty)
{
  return in_poly_internal(testx, testy, vertx, verty);
}

// Recursive function to determine fractional pixel coverage
double pixelCoverPoly(double x, double y, const NumericVector &vertx, const NumericVector &verty, int depth) {
  if (depth == 0) {
    return in_poly_internal(x, y, vertx, verty);
  }
  
  double quarter = 0.25 / (1 << (depth - 1)); // (1 << (depth - 1)) is equivalent to pow(2, depth - 1), but faster
  double coverage = 0.0;
  
  double x_min  = x - quarter;
  double x_plus = x + quarter;
  double y_min  = y - quarter;
  double y_plus = y + quarter;
  
  coverage += pixelCoverPoly(x_min, y_min, vertx, verty, depth - 1);
  coverage += pixelCoverPoly(x_min, y_plus, vertx, verty, depth - 1);
  coverage += pixelCoverPoly(x_plus, y_min, vertx, verty, depth - 1);
  coverage += pixelCoverPoly(x_plus, y_plus, vertx, verty, depth - 1);
  
  return coverage / 4.0;
}

// [[Rcpp::export]]
NumericVector profoundPolyCover(NumericVector x,
                                NumericVector y,
                                NumericVector vertx,
                                NumericVector verty,
                                int depth = 4,
                                int nthreads = 1) {
  
  const int n = x.size();
  NumericVector result(n);
  
  double vertx_min = min(vertx) - 0.5;
  double vertx_max = max(vertx) + 0.5;
  double verty_min = min(verty) - 0.5;
  double verty_max = max(verty) + 0.5;
  
  #ifdef _OPENMP
    // Parallelize the main loop. Use 'if' to avoid overhead for tiny n.
  #pragma omp parallel for schedule(dynamic, 10) if(n > 100) num_threads(nthreads)
  #endif
  for (int i = 0; i < n; ++i) {
    // Make loop-local copies to avoid data races
    if (x[i] > vertx_min && x[i] < vertx_max) {
      if (y[i] > verty_min && y[i] < verty_max) {
        result[i] = pixelCoverPoly(x[i], y[i], vertx, verty, depth);
      } else {
        result[i] = 0.0;
      }
    } else {
      result[i] = 0.0;
    }
  }
  
  return result;
}
