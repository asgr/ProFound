#include <Rcpp.h>
using namespace Rcpp;

/**
 *  Represents an Akima spline
 *  Given six input (X,Z) pairs this akima spline fits a smooth curve between points 3 and 4.
 *  This akima spline that is valid for interpolations between the third and fourth points (inclusive) of the six (X,Z) pairs
 *  used in its construction
 *  Reference :  A new method of interpolation and smooth curve fitting based on local procedures, Hiroshi Akima,
 *               Journal Of The Association Of Computing Machinery, Volume 17, number 4, October, 1970, pp. 589 - 602.
 *               https://dl.acm.org/citation.cfm?id=321609
 *
 */
class Coeff
{
public:

  /**
   * Evaluate this akima spline interpolating polynomial at x
   */
  double interpValue(double x) const
  {
    double xx = x - x3;
    return y3 + xx * (w0 + xx * (w2 + xx * w3));
  }

  /**
   *   Finds the slope of the middle point in a sequence of five (x,z) pairs.
   */
  static double calcSlopeAtMiddle(const double *x, const double *z)
  {
    int i;
    double S13, S34, S12, S24, W2, W3, Z, DEN;
    double A[4], B[4];

    // Calculate differences between input points
    for (i = 0; i < 4; i++) {
      A[i] = x[i + 1] - x[i];
      B[i] = z[i + 1] - z[i];
    }

    S13 = A[0] * B[2] - A[2] * B[0];
    S34 = A[2] * B[3] - A[3] * B[2];
    S12 = A[0] * B[1] - A[1] * B[0];
    S24 = A[1] * B[3] - A[3] * B[1];
    W2  = sqrt(std::abs(S13 * S34));
    W3  = sqrt(std::abs(S12 * S24));
    Z   = W2 * B[1] + W3 * B[2];
    DEN = W2 * A[1] + W3 * A[2];

    if (DEN == 0.0) {
      return 0.0;
    }
    return Z / DEN;
  }

  // valid interpolation range
  double x3 = 0;
  double x4 = 0;

  // extra coeffs for akima
  double w0 = 0;
  double w1 = 0;
  double w2 = 0;
  double w3 = 0;
  double y3 = 0;

};

/**
 *  Represents an array of Akima splines covering a set of (at least 5) input (X,Z) pairs.
 *  Given six input (X,Z) pairs this akima spline fits a smooth curve between points 3 and 4.
 *  The array of akima splines are valid for interpolations between the first and last points (inclusive)
 *  used in its construction.
 *  The evaluation of the Akima splines for the first and last two intervals are adapted to make use of less than 5 (X,Z) pairs.
 *  Reference :  A new method of interpolation and smooth curve fitting based on local procedures, Hiroshi Akima,
 *               Journal Of The Association Of Computing Machinery, Volume 17, number 4, October, 1970, pp. 589 - 602.
 *               https://dl.acm.org/citation.cfm?id=321609
 *
 */
class adacsakima {
public:
  adacsakima();
  adacsakima(int npts, const double_t *xorig, const double_t *yorig);
  bool Initialise(int npts,const double_t *xorig, const double_t *yorig);

  bool isValid() const;
  double_t InterpValue(double_t x) const;
private:
  bool isvalid;
  int ncoeffs;
  std::vector<Coeff> coeffs;
};

//==================================================================================
/**
 * Interpolate a 2D regular grid using akima spline interpolation
 */
// [[Rcpp::export(".interpolateAkimaGrid")]]
void interpolateAkimaGrid(NumericVector xseq, NumericVector yseq,
                          NumericMatrix tempmat_sky, NumericMatrix output) {
  /*
   * An Matrix element is at (row,col)
   * The elements of a row stack vertically
   * Any row I is to the right of row I-1
   */
  int myxnpts = output.nrow();
  int myynpts = output.ncol();
  const double_t* myx=REAL(xseq);
  const double_t* myy=REAL(yseq);
  int ncol=tempmat_sky.ncol();
  int nrow=tempmat_sky.nrow();

  std::vector<double_t> xin,zin;
  std::vector<adacsakima> akimaCOL;
  xin.reserve(ncol);
  zin.reserve(nrow);
  akimaCOL.reserve(ncol);
  for (int i = 1; i <= nrow; i++)
  {
    xin.push_back(0);
    zin.push_back(0);
  }
  for (int j = 1; j <= ncol; j++) {
    int ii=(j-1)*ncol+(1-1);
    for (int i = 1; i <= nrow; i++,ii++) {
      xin[i-1] = myx[i-1];
      zin[i-1] = tempmat_sky(i-1,j-1);
    }
    adacsakima thisspline;
    thisspline.Initialise(nrow,xin.data(),zin.data());
    akimaCOL.push_back(thisspline);
  }

  // For each vertical row
  for (int i = 1; i <= myxnpts; i++) {
    // For a spline to interpolate vertically along the elements of the row
    double_t x = -0.5+i;
    for (int j = 1; j <= ncol; j++) {
      xin[j-1] = myy[j-1];
      zin[j-1] = akimaCOL[j-1].InterpValue(x);
    }
    adacsakima thisspline;
    thisspline.Initialise(ncol,xin.data(),zin.data());

    // Interpolate vertically for each element (j) in the current (i) output row
    for (int j = 1; j <= myynpts; j++) {
      double_t y = -0.5+j;
      output(i-1,j-1) = thisspline.InterpValue(y);
    }
  }
}
/**
 * Interpolate a 2D regular grid using bilinear interpolation
 */
// [[Rcpp::export(".interpolateLinearGrid")]]
void interpolateLinearGrid(NumericVector xseq, NumericVector yseq, NumericMatrix tempmat_sky,
                           NumericMatrix output) {
  /*
   * An Matrix element is at (row,col)
   * The elements of a row stack vertically
   * Any row I is to the right of row I-1
   */
  int myxnpts = output.nrow();
  int myynpts = output.ncol();
  const double_t* myx=REAL(xseq);
  const double_t* myy=REAL(yseq);
  int ncol=tempmat_sky.ncol();
  int nrow=tempmat_sky.nrow();

  // For each vertical row
  for (int i = 1; i <= myxnpts; i++) {
    // For a spline to interpolate vertically along the elements of the row
    double_t x = -0.5+i;
    // find the left and right index ibnto xseq
    int left_index = -1;
    int right_index = -1;
    for (int ii = 1; ii < nrow; ii++) {
      if (myx[ii-1] <= x && myx[ii] >= x) {
        left_index = ii-1;
        right_index = ii;
        break;
      }
    }

    //Rcpp::Rcout << "x="<<x<<" xindex="<<left_index<<" "<<right_index<<"\n";
    //Rcpp::Rcout << "x="<<x<<" xleft="<<myx[left_index]<<" "<<myx[right_index]<<"\n";
    int top_index = -1;
    int bottom_index = -1;
    for (int j = 1; j < myynpts; j++) {
      double y = -0.5+j;
      for (int jj = 1; jj < ncol; jj++) {
        if (myy[jj-1] <= y && myy[jj] >= y) {
          top_index = jj-1;
          bottom_index = jj;
          // p1...p2
          // .     .
          // .     .
          // p3...p4
          double p1 = tempmat_sky(left_index,top_index);
          double p2 = tempmat_sky(right_index,top_index);
          double p3 = tempmat_sky(left_index,bottom_index);
          double p4 = tempmat_sky(right_index,bottom_index);

          double xlambda = (x-myx[left_index])/(myx[right_index]-myx[left_index]);
          double ylambda = (y-myy[top_index])/(myy[bottom_index]-myy[top_index]);
          double ztop = p1 * (1.0-xlambda) + p2 * xlambda;
          double zbottom = p3 * (1.0-xlambda) + p4 * xlambda;
          output(i-1,j-1) = ztop * (1.0-ylambda) +zbottom * ylambda;
          //Rcpp::Rcout << "y="<<y<<" yindex="<<top_index<<" "<<bottom_index<<" "<<p1<<" "<<p2<<" "<<p3<<" "<<p4<<" result="<<output(i-1,j-1)<<"\n";
          break;
        }
      }
    }
  }
}

/**
 *  Represents an array of Akima splines covering a set of (at least 5) input (X,Z) pairs.
 *  Given six input (X,Z) pairs this akima spline fits a smooth curve between points 3 and 4.
 *  The array of akima splines are valid for interpolations between the first and last points (inclusive)
 *  used in its construction.
 *  The evaluation of the Akima splines for the first and last two intervals are adapted to make use of less than 5 (X,Z) pairs.
 *  Reference :  A new method of interpolation and smooth curve fitting based on local procedures, Hiroshi Akima,
 *               Journal Of The Association Of Computing Machinery, Volume 17, number 4, October, 1970, pp. 589 - 602.
 *               https://dl.acm.org/citation.cfm?id=321609
 *
 */
adacsakima::adacsakima() {
  ncoeffs = 0;
  isvalid = false;
}
bool adacsakima::isValid() const {
  return isvalid;
}

/**---------------------------------------------------------------
 *   PURPOSE:
 *       To determine the array of akima spline that fit a set of points (Xi,Yi) i = 1,n.  The table of coefficients stored
 * in coeffs can be used with InterpValue to perform spline interpolations.
 *
 */
adacsakima::adacsakima(int npts, const double_t *xorig, const double_t *yorig)
{
  Initialise(npts,xorig,yorig);
}
/**---------------------------------------------------------------
 *   PURPOSE:
 *       To determine the akima spline that fits a set of points (Xi,Yi) i = 1,n.  The table of coefficients stored
 * in wkarea can be used withInterp to perform spline interpolations.
 *
 *  s3 slope at start of interval s4 slope at end
 */
/**
 * Computes the array of akima splines - notes the special treatment of intervals within 2 of the start and end.
 */
bool adacsakima::Initialise(int npts, const double_t *xorig, const double_t *yorig)
{
  double_t r2 = 2, r3 = 3;
  double_t s3 = 0, s4 = 0, dx = 0., dy = 0., p2, p3, x3, x4, y3, y4;
  int i;

  if (npts < 5)	// Error status if validity test fails
  {
    return false;
  }

  // One spline for each interval (xorig[i], xorig[i+1])
  ncoeffs = npts-1;
  coeffs.reserve(ncoeffs);

  for (i = 0; i < ncoeffs; i++)
  {
    x3 = xorig[i];
    y3 = yorig[i];

    x4 = xorig[i + 1];
    y4 = yorig[i + 1];
    dx = x4 - x3;
    dy = y4 - y3;

    // check for boundary conditions
    if (i == 0)
    {
      // do first interval
      s3 = dy / dx;
      s4 = Coeff::calcSlopeAtMiddle(&(xorig[0]), &(yorig[0]));
      s4 = (s4 + s3) / 2;
    }
    else if (i == 1)
    {
      // do second interval
      s3 = dy / dx;
      s4 = Coeff::calcSlopeAtMiddle(&(xorig[0]), &(yorig[0]));
      s3 = (s4 + s3) / 2;
    }
    else if (i == ncoeffs - 2)
    {
      // to second last interval
      s3 = Coeff::calcSlopeAtMiddle(&xorig[i - 2], &yorig[i - 2]);
      s4 = dy / dx;
      s4 = (s4 + s3) / 2;
    }
    else if (i == ncoeffs - 1)
    {
      // do last interval
      s3 = Coeff::calcSlopeAtMiddle(&xorig[ncoeffs - 4], &yorig[ncoeffs - 4]);
      x3 = xorig[npts - 2];
      y3 = yorig[npts - 2];
      x4 = xorig[npts - 1];
      y4 = yorig[npts - 1];
      dx = x4 - x3;
      dy = y4 - y3;
      s4 = dy / dx;
      s3 = (s4 + s3) / 2;
    }
    else
    {
      // do the "pure" akima intervals
      // Determine slope at beginning and end of current interval.
      s3 = Coeff::calcSlopeAtMiddle(&xorig[i - 2], &yorig[i - 2]);
      s4 = Coeff::calcSlopeAtMiddle(&xorig[i - 1], &yorig[i - 1]);
    }

    //
    // Compute coefficients of cubic equation for this akima
    //
    p2 = (r3 * dy / dx - r2 * s3 - s4) / dx;
    p3 = (s3 + s4 - r2 * dy / dx) / (dx * dx);
    //
    // store slopes and interploting coefficients
    Coeff coeff;
    coeff.w0 = s3;
    coeff.w1 = s4;
    coeff.w2 = p2;
    coeff.w3 = p3;
    coeff.x3 = x3;
    coeff.x4 = x4;
    coeff.y3 = y3;
    coeffs.push_back(coeff);
  }
  return true;
}
double_t adacsakima::InterpValue(double_t x) const
{
  for (int spline = 0; spline < ncoeffs; spline++)
  {
    // Determine if this spline's interval covers x
    if (x >= coeffs[spline].x3 && x <= coeffs[spline].x4)
    {
      // Evaluate this akima spline interpolating polynomial at x
      return coeffs[spline].interpValue(x);
    }
  }
  return 0.0;
}
