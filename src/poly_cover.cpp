#include <Rcpp.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

static double in_poly_internal(double testx, double testy, const NumericVector &poly_x, const NumericVector &poly_y)
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

double in_poly(double testx, double testy, NumericVector poly_x, NumericVector poly_y)
{
  return in_poly_internal(testx, testy, poly_x, poly_y);
}

// Recursive function to determine fractional pixel coverage

double pixelCoverPoly(double x, double y, const NumericVector &poly_x, const NumericVector &poly_y, int depth) {
  if (depth == 0) {
    return in_poly_internal(x, y, poly_x, poly_y);
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

// [[Rcpp::export]]
double profoundPolyFlux(NumericMatrix image,
                        NumericVector poly_x,
                        NumericVector poly_y,
                        int depth = 4,
                        int nthreads = 1) {
  int nrow = image.nrow();
  int ncol = image.ncol();
  
  // This is different to above because the R coord system is [0.5,0.5] for the first pixel, versus [0,0] in Rcpp
  double poly_x_min = min(poly_x) - 1;
  double poly_x_max = max(poly_x);
  double poly_y_min = min(poly_y) - 1;
  double poly_y_max = max(poly_y);
  
  int start_row = std::max(0.0, floor(poly_x_min));
  int end_row = std::min(nrow - 1.0, ceil(poly_x_max));
  int start_col = std::max(0.0, floor(poly_y_min));
  int end_col = std::min(ncol - 1.0, ceil(poly_y_max));
  
  double sum = 0.0;
  
  #ifdef _OPENMP
    // Parallelize the main loop
  #pragma omp parallel for schedule(static) num_threads(nthreads)
  #endif
  for (int i = start_row; i <= end_row; ++i) {
    for (int j = start_col; j <= end_col; ++j) {
      if(!NumericMatrix::is_na(image(i, j))){
        sum += image(i, j)*pixelCoverPoly(i + 0.5, j + 0.5, poly_x, poly_y, depth);
      }
    }
  }
  
  return sum;
}
