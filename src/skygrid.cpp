#include <Rcpp.h>
#include "akima.h"

#define adacs_BOTH 1
#define adacs_LO 2
#define adacs_HI 3
#define adacs_SD 4
#define adacs_RBOTH 5
#define adacs_RLO 6
#define adacs_RHI 7
#define adacs_RSD 8

#define adacs_AUTO 1
#define adacs_SET 2

#define adacs_MEDIAN 1
#define adacs_MEAN 2
#define adacs_MODE 3
#define adacs_RMEDIAN 4
#define adacs_RMEAN 5
#define adacs_RMODE 6

#define adacs_CLASSIC_BILINEAR 1
#define adacs_AKIMA_BICUBIC 2

/**
 * Simple functors for testing postiveness/negativeness, and setting NaNs on
 * a vector for a given condition
 */
struct positive_or_zero {
  bool operator()(double n) { return n >= 0; }
};

struct negative_or_zero {
  bool operator()(double n) { return n <= 0; }
};

template <typename Predicate>
void set_nan_if(Rcpp::NumericVector numbers, Predicate condition)
{
  for (auto &n: numbers) {
    if (condition(n)) {
      n = NA_REAL;
    }
  }
}

/**
 * A simple histogram that first accumulates values from a vector following
 * certain conditions, and then allows quick calculation of quantiles
 */
class AdacsHistogram {

public:
  void accumulate(Rcpp::NumericVector x, int nbins=16384, double minv=NA_REAL, double maxv=NA_REAL)
  {
    accumulate(x, nbins, minv, maxv,
               [](double) { return true; });
  }

  void accumulateLO(Rcpp::NumericVector x, double offset=0, int nbins=16384, double minv=NA_REAL, double maxv=NA_REAL)
  {
    accumulate(x, nbins, minv, maxv,
               [offset](double x) { return x < offset; });
  }

  void accumulateHI(Rcpp::NumericVector x, double offset=0, int nbins=16384, double minv=NA_REAL, double maxv=NA_REAL)
  {
    accumulate(x, nbins, minv, maxv,
               [offset](double x) { return x > offset; });
  }

  double quantile(double quantile, double offset=0) const
  {
    if (_non_null_sample_count == 0) {
      return NA_REAL;
    }
    int count = 0;
    double quantileValue = _min - offset;
    double binwidth = (_max - _min) / _nbins;
    int target_count = _non_null_sample_count * quantile;
    for (auto bin_count: _histogram) {
      if (count >= target_count) {
        return quantileValue;
      }
      quantileValue += binwidth;
      count += bin_count;
    }
    return quantileValue;
  }

private:

  template <typename Predicate>
  void accumulate(Rcpp::NumericVector x, int nbins, double minv, double maxv,
                  Predicate condition)
  {
    _nbins = nbins;
    _min = std::numeric_limits<double>::max();
    _max = -_min;
    _non_null_sample_count = 0;
    _null_sample_count = 0;

    for (auto _x: x) {
      if (!std::isnan(_x) && condition(_x)) {
        _non_null_sample_count++;
        _min = std::min(_min, _x);
        _max = std::max(_max, _x);
      }
    }
    _null_sample_count = x.size() - _non_null_sample_count;

    if (_non_null_sample_count < 1)
      return;

    if (!std::isnan(minv) && !std::isnan(maxv)) {
      _min = minv;
      _max = maxv;
    }
    if (_min == _max) {
      return;
    }

    // histogram
    _toolow = 0;
    _toohigh = 0;
    _histogram.resize(_nbins);
    double value_to_bin_index = _nbins - 1;
    value_to_bin_index /= _max - _min;
    for (auto _x: x) {
      if (!std::isnan(_x) && condition(_x)) {
        int index = (_x - _min) * value_to_bin_index;
        if (index < 0) {
          _toolow++;
        } else if (index >= nbins) {
          _toohigh++;
        } else {
          _histogram[index]++;
        }
      }
    }
  }

  int _nbins;
  double _min;
  double _max;
  int _non_null_sample_count;
  int _null_sample_count;
  std::vector<int> _histogram;
  int _toolow;
  int _toohigh;
};

//==================================================================================
/*
* Search neighbourhood of (loc1, loc2) for at least skypixmin viable "sky" values.
* Expand the box by boxadd until enough found.
*/
Rcpp::NumericVector Cadacs_FindSkyCellValues(
    Rcpp::NumericMatrix image, Rcpp::Nullable<Rcpp::IntegerMatrix> objects,
    Rcpp::Nullable<Rcpp::IntegerMatrix> mask,
    const double loc1, const double loc2, const double box1, const double box2,
    const double boxadd1, const double boxadd2, const int skypixmin,
	const int boxiters)
{
  // R is 1 relative
  int iloc1 = (int)(loc1+0.5);
  int iloc2 = (int)(loc2+0.5);
  int ibox1 = (int)(box1/2);
  int ibox2 = (int)(box2/2);
  int nrow = image.nrow();
  int ncol = image.ncol();

  const double* iiimage=REAL(image);
  Rcpp::IntegerMatrix iobjects;
  const int32_t* iiobjects=NULL;
  if (objects.isNotNull()) {
    iobjects = Rcpp::as<Rcpp::IntegerMatrix>(objects);
    iiobjects=INTEGER(objects.get());
  }
  Rcpp::IntegerMatrix imask;
  const int32_t* iimask=NULL;
  if (mask.isNotNull()) {
    imask = Rcpp::as<Rcpp::IntegerMatrix>(mask);
    iimask=INTEGER(mask.get());
  }

  int iboxadd1=0;
  int iboxadd2=0;
  int skyN=0;
  int iterN=0;
  int ssrow = 1;
  int eerow = 0;
  int sscol = 1;
  int eecol = 0;

  while(skyN<skypixmin && iterN<=boxiters){
    skyN = 0;
    ibox1 += iboxadd1;
    ibox2 += iboxadd2;
    ssrow = std::max(1,iloc1-ibox1);
    eerow = std::min(nrow,iloc1+ibox1);
    sscol = std::max(1,iloc2-ibox2);
    eecol = std::min(ncol,iloc2+ibox2);

    for (int j = sscol; j <= eecol; j++) {
      int ii=(j-1)*nrow+(ssrow-1);
      for (int i = ssrow; i <= eerow; i++,ii++) {
        // Count sky cells (sky cells are those NOT masked out and NOT objects)
        if ((iiobjects!=NULL)) {
          if (iiobjects[ii]==0 && (iimask==NULL || iimask[ii]==0)) {
            skyN++;
          }
        } else if (iimask!=NULL) {
          if (iimask[ii]==0) {
            skyN++;
          }
        } else {
          skyN++;
        }
      }
    }
    iterN++;
    iboxadd1 = (int)(boxadd1/2);
    iboxadd2 = (int)(boxadd2/2);

  }
  // copy sky cell values to vec and return
  // Rcpp::Rcout << "\nCbox "<<ssrow<<" "<<eerow<<" "<<sscol<<" "<<eecol<<"\n";
  
  // Rcpp::Rcout << skyN << "\n";
  
  Rcpp::NumericVector vec(skyN);
  int k=0;
  for (int j = sscol; j <= eecol; j++) {
    int ii=(j-1)*nrow+(ssrow-1);
    for (int i = ssrow; i <= eerow; i++,ii++) {
      //Rcpp::Rcout << i << "\n";
      //Rcpp::Rcout << ii << "\n";
      //Rcpp::Rcout << k << "\n";
      if ((iiobjects!=NULL)) {
        if (iiobjects[ii]==0 && (iimask==NULL || iimask[ii]==0)) {
          vec[k++] = iiimage[ii];
        }
      } else if (iimask!=NULL) {
        if (iimask[ii]==0) {
          vec[k++] = iiimage[ii];
        }
      } else {
        vec[k++] = iiimage[ii];
      }
    }
  }
  // Rcpp::Rcout << vec[k-1] << k << " FINE!\n";
  return vec;
}

//==================================================================================
/**
 * C++ version of R quantile, low variant
 */
double Cadacs_quantileLO(Rcpp::NumericVector x, double quantile, const double offset, int nbins=16384, double minv=NA_REAL, double maxv=NA_REAL) {
  // The population we want the quantile for is x-offset where x<offset
  AdacsHistogram histogram;
  histogram.accumulateLO(x, offset, nbins, minv, maxv);
  return histogram.quantile(quantile, offset);
}

/**
 * C++ version of R quantile, high variant
 */
double Cadacs_quantileHI(Rcpp::NumericVector x, double quantile, const double offset, int nbins=16384, double minv=NA_REAL, double maxv=NA_REAL) {
  // The population we want the quantile for is x-offset where x>offset
  AdacsHistogram histogram;
  histogram.accumulateHI(x, offset, nbins, minv, maxv);
  return histogram.quantile(quantile, offset);
}

/**
 * C++ version of R stats::mean
 */
double Cadacs_mean(Rcpp::NumericVector x) {
  const double* myx=REAL(x);
  int size = x.size();
  double mean=0;
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      mean += myx[i];
    }
  }
  if (non_null_sample_count==0)
    return NA_REAL;
  return mean/non_null_sample_count;
}

double Cadacs_population_variance(Rcpp::NumericVector x, const double offset) {
  const double* myx=REAL(x);
  int size = x.size();
  double v=0;
  double sum_sq=0;
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      v = myx[i]-offset;
      v *= v;
      sum_sq += v;
    }
  }
  if (non_null_sample_count==0)
    return NA_REAL;
  double N=non_null_sample_count;
  return sum_sq/N;
}

double Cadacs_sample_variance(Rcpp::NumericVector x, const double offset) {
  const double* myx=REAL(x);
  int size = x.size();
  double v=0;
  double sum=0;
  double sum_sq=0;
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      v = myx[i]-offset;
      sum += v;
      v *= v;
      sum_sq += v;
    }
  }
  if (non_null_sample_count<=1)
    return NA_REAL;
  double N=non_null_sample_count;
  return (N*sum_sq - sum*sum)/(N * (N - 1));
  //return sqrt(sum_sq);
}

double Cadacs_median(Rcpp::NumericVector x) {
  AdacsHistogram histogram;
  histogram.accumulate(x);
  return histogram.quantile(0.5);
}

double Cadacs_mode(Rcpp::NumericVector x) {
  const double* iiix=REAL(x);
  int size = x.size();
  std::vector<double> myx (iiix, iiix+size);
  double min=std::numeric_limits<double>::max();
  double max=std::numeric_limits<double>::min();
  int non_null_sample_count=0;
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      non_null_sample_count++;
      min = std::min(min,myx[i]);
      max = std::max(max,myx[i]);
    }
  }

  // histogram
  int levels = 100;
  std::vector<int> histogram;
  histogram.resize(levels);
  for (int i=0;i<levels;i++)
  {
    histogram[i] = 0;
  }
  double value_to_bin_index = (levels-1);
  value_to_bin_index /= (max - min);
  for (int i=0;i<size;i++)
  {
    if (!std::isnan(myx[i])) {
      int index = (myx[i] - min)*value_to_bin_index;
      histogram[index]++;
    }
  }

  double current_bin_lower=min;
  double mode=current_bin_lower;

  double binwidth = (max - min)/levels;
  int max_count = 0;
  for (int i=0;i<levels;i++) {
    if (histogram[i]>max_count)
    {
      max_count = histogram[i];
      mode = current_bin_lower;
    }
    current_bin_lower += binwidth;
  }
  return mode;
}

/**
 * Sort based (rather than Histogram based) method to clip outliers (A histogram equivalent has not been evaluated)
 */
Rcpp::NumericVector Cadacs_magclip(Rcpp::NumericVector x, const int sigma, const int clipiters, const double sigmasel, const int estimate){
  const double* iiix=REAL(x);
  int nb = x.length();
  std::vector<double> myx (iiix, iiix+nb);
  int length=0;
  for (int i=0;i<nb;i++)
  {
    if (!std::isnan(myx[i])) {
      myx[length++] = myx[i];
    }
  }
  std::sort (myx.begin(), myx.begin()+length, std::less<double>()); // ascending

  int newlen = length;
  if(clipiters>0 && length>0){
    double sigcut=R::pnorm(sigmasel, 0.0, 1.0, 1, 0);

    for(int iteration=0; iteration<clipiters; iteration++){
      if(newlen<=1)
        break;
      int oldlen=newlen;
      double roughmed=myx[newlen/2-1];
      double clipsigma=sigma;
      if (sigma==1) {
        double l1=std::max(newlen,2);
        double y=1.0-2.0/l1;
        clipsigma = R::qnorm(y, 0.0, 1.0, 1, 0);
      }

      double vallims = 0;
      switch(estimate) {
      case 1:
        vallims = clipsigma*(myx[sigcut*newlen-1]-myx[(1-sigcut)*newlen-1])/2/sigmasel;
        break;
      case 2:
        vallims = clipsigma*(roughmed-myx[(1-sigcut)*newlen-1])/sigmasel;
        break;
      case 3:
        vallims = clipsigma*(myx[sigcut*newlen-1]-roughmed)/sigmasel;
        break;
      }
      newlen = 0;
      for (int i=0;i<oldlen;i++)
      {
        if(myx[i]>=(roughmed-vallims) && myx[i]<=(roughmed+vallims))
        {
          myx[newlen++] = myx[i];
        }
      }
      if(oldlen==newlen)
        break;
    }
  }
  // copy sky cell values to vec and return
  Rcpp::NumericVector vec(newlen);
  for (int i=0;i<newlen;i++)
  {
    vec[i] = myx[i];
  }
  return vec;
}

Rcpp::NumericVector Cadacs_SkyEstLoc(Rcpp::NumericMatrix image,
    Rcpp::Nullable<Rcpp::IntegerMatrix> objects, Rcpp::Nullable<Rcpp::IntegerMatrix> mask,
    const double loc1, const double loc2, const double box1, const double box2,
    const double boxadd1, const double boxadd2,
    const int skypixmin, const int boxiters, const int doclip, const int skytype,
    const int skyRMStype, const double sigmasel)
{
  Rcpp::NumericVector select = Cadacs_FindSkyCellValues(image, objects, mask, loc1, loc2, box1, box2, boxadd1, boxadd2, skypixmin, boxiters);
  Rcpp::NumericVector clip;
  Function Fquantile("quantile");
  if(doclip) {
    clip = Cadacs_magclip(select,adacs_AUTO,5,sigmasel,adacs_LO);
  } else {
    clip = select;
  }
  double skyloc=0.0;
  switch (skytype) {
  case adacs_MEDIAN:
    skyloc = Cadacs_median(clip);
    // Rcpp::Rcout << skyloc << "\n";
    break;
  case adacs_RMEDIAN:
    skyloc = Rcpp::median(clip);
    break;
  case adacs_MEAN:
    skyloc = Cadacs_mean(clip);
    break;
  case adacs_RMEAN:
    skyloc = Rcpp::mean(clip);
    break;
  case adacs_MODE:
    skyloc = Cadacs_mode(clip);
    break;
  case adacs_RMODE:
  {
    Rcpp::Environment profound = Rcpp::Environment::namespace_env("ProFound");
    Rcpp::Function mode= profound["adacs_mode"];
    skyloc = REAL(mode(clip))[0];
  }
    break;
  }

  // This is done by a couple of "case" blocks below, so let's not duplicate
  auto subtract_skyloc_from_clip = [&]() {
    for (auto &n: clip) {
      n -= skyloc;
    }
  };

  // As above, centralising calculation of pnorm
  auto get_pnorm = [&](bool is_low) {
    if (is_low) {
      return R::pnorm(-sigmasel, 0.0, 1.0, 1, 0) * 2;
    }
    return (R::pnorm(sigmasel, 0.0, 1.0, 1, 0) - 0.5)  * 2;
  };

  // As above, let's bring together what's done equally in a "case" bloocks
  auto get_rms = [&](Rcpp::NumericVector values, bool is_low, bool use_r) {
    auto pnorm = get_pnorm(is_low);
    double quantile;
    if (use_r) {
      quantile = REAL(Fquantile(values, pnorm, true))[0];
    }
    else if (is_low) {
      quantile = Cadacs_quantileLO(clip, pnorm, skyloc);
    }
    else {
      quantile = Cadacs_quantileHI(clip, pnorm, skyloc);
    }
    return std::abs(quantile) / sigmasel;
  };

  // Calculate skyRMSloc depending on the type
  double skyRMSloc = 0.0;
  double lo, hi;
  Rcpp::DoubleVector templo, temphi;
  switch (skyRMStype) {

  case adacs_LO:
    skyRMSloc = get_rms(clip, true, false);
    break;

  case adacs_RLO:
    subtract_skyloc_from_clip();
    set_nan_if(clip, positive_or_zero{});
    skyRMSloc = get_rms(clip, true, true);
    break;

  case adacs_HI:
    skyRMSloc = get_rms(clip, false, false);
    break;

  case adacs_RHI:
    subtract_skyloc_from_clip();
    set_nan_if(clip, negative_or_zero{});
    skyRMSloc = get_rms(clip, false, true);
    break;

  case adacs_BOTH:
    lo = get_rms(clip, true, false);
    hi = get_rms(clip, false, false);
    skyRMSloc = (lo + hi) / 2;
    break;

  case adacs_RBOTH:
    subtract_skyloc_from_clip();
    templo = Rcpp::clone(clip);
    temphi = Rcpp::clone(clip);
    set_nan_if(templo, positive_or_zero{});
    set_nan_if(temphi, negative_or_zero{});
    lo = get_rms(templo, true, true);
    hi = get_rms(temphi, false, true);
    skyRMSloc = (lo + hi) / 2;
    break;

  case adacs_SD:
    skyRMSloc = sqrt(Cadacs_population_variance(clip, skyloc));
    break;

  case adacs_RSD:
    // Its ok to modify clip since its a fresh object and will not be used later
    subtract_skyloc_from_clip();
    skyRMSloc = sqrt(Rcpp::var(clip));
    break;
  }

  Rcpp::NumericVector result(2);
  result[0] = skyloc;
  result[1] = skyRMSloc;
  return result;
}


// [[Rcpp::export(".Cadacs_MakeSkyGrid")]]
void Cadacs_MakeSkyGrid(
    Rcpp::NumericMatrix image, Rcpp::NumericMatrix sky, Rcpp::NumericMatrix skyRMS,
    Rcpp::Nullable<Rcpp::IntegerMatrix> objects = R_NilValue,
    Rcpp::Nullable<Rcpp::IntegerMatrix> mask = R_NilValue,
    const int box1 = 100, const int box2 = 100,
    const int grid1 =100, const int grid2 = 100,
    const int boxadd1 = 50, const int boxadd2 = 50,
    const int type = 2, const int skypixmin = 5000,
    const int boxiters = 0, const int doclip = 1,
    const int skytype = 1, const int skyRMStype = 2,
    const double sigmasel = 1)
{
  // box MUST NOT be larger than the input image
  double box[2] = {(double)box1, (double)box2};
  if(box[0]>image.nrow())
    box[0]=image.nrow();
  if(box[1]>image.ncol())
    box[1]=image.ncol();

  double grid[2] = {(double)grid1, (double)grid2};
  if(grid[0]>image.nrow())
    grid[0]=image.nrow();
  if(grid[1]>image.ncol())
    grid[1]=image.ncol();

  // tile over input image with tile size (grid) and no overlap
  // xseq,yseq give the centres of each tile
  int tile_nrows=0;
  double x_tile_centre=grid[0]/2;
  while (x_tile_centre<image.nrow()) {
    tile_nrows++;
    x_tile_centre += grid[0];
  }
  int tile_ncols=0;
  double y_tile_centre=grid[1]/2;
  while (y_tile_centre<image.ncol()) {
    tile_ncols++;
    y_tile_centre += grid[1];
  }

  // add room for linearly extrapolated padding
  tile_nrows += 2;
  tile_ncols += 2;

  // Construct the vector of tile centroids
  Rcpp::NumericVector xseq(tile_nrows);
  Rcpp::NumericVector yseq(tile_ncols);
  x_tile_centre=grid[0]/2 - grid[0];
  for (int i=0; i<tile_nrows; i++) {
    xseq[i] = x_tile_centre;
    x_tile_centre += grid[0];
  }
  y_tile_centre=grid[1]/2 - grid[1];
  for (int i=0; i<tile_ncols; i++) {
    yseq[i] = y_tile_centre;
    y_tile_centre += grid[1];
  }

  Rcpp::NumericMatrix z_sky_centre(tile_nrows, tile_ncols);
  Rcpp::NumericMatrix z_skyRMS_centre(tile_nrows, tile_ncols);

  bool hasNaNs=false;
  x_tile_centre=grid[0]/2;
  for (int i=1; i<tile_nrows-1; i++) {
    x_tile_centre = xseq[i];
    for (int j=1; j<tile_ncols-1; j++) {
      y_tile_centre = yseq[j];
      Rcpp::NumericVector z_tile_centre = Cadacs_SkyEstLoc(image, objects, mask,
                                                           x_tile_centre, y_tile_centre,
                                                           box1, box2,
                                                           boxadd1, boxadd2,
                                                           skypixmin, boxiters,
                                                           doclip, skytype, skyRMStype, sigmasel);
      if (std::isnan(z_tile_centre[0]) || std::isnan(z_tile_centre[1])) {
        hasNaNs = true;
      }

      z_sky_centre(i, j) = z_tile_centre[0];
      z_skyRMS_centre(i, j) = z_tile_centre[1];
    }
  }
  if (hasNaNs) {
    // Rcpp::Rcout << "HERE!\n";
    // Replace any NaN's with reasonable substitute
    // initialise the pad area before getting the medians
    for (int i=0; i<tile_nrows; i++) {
      // sky
      z_sky_centre(i,0) = NA_REAL;
      z_sky_centre(i,tile_ncols-1) = NA_REAL;
      // skyRMS
      z_skyRMS_centre(i,0) = NA_REAL;
      z_skyRMS_centre(i,tile_ncols-1) = NA_REAL;
    }
    for (int i=0; i<tile_ncols; i++) {
      // sky
      z_sky_centre(0, i) = NA_REAL;
      z_sky_centre(tile_nrows-1, i) = NA_REAL;
      // skyRMS
      z_skyRMS_centre(0, i) = NA_REAL;
      z_skyRMS_centre(tile_nrows-1, i) = NA_REAL;
    }
    double medianSkyCentre=Cadacs_median(z_sky_centre);
    double medianSkyRMSCentre=Cadacs_median(z_skyRMS_centre);
    // replace NaN's now
    for (int i=1; i<tile_nrows-1; i++) {
      for (int j=1; j<tile_ncols-1; j++) {
        if (std::isnan(z_sky_centre(i, j)))
          z_sky_centre(i, j) = medianSkyCentre;
        if (std::isnan(z_skyRMS_centre(i, j)))
          z_skyRMS_centre(i, j) = medianSkyRMSCentre;
      }
    }
  }

  // Padding
  //work out the second point for linear extrapolation (the first one is at 1+1 and length(seq)-1)
  int xstart=std::min(2,tile_nrows-2);
  int ystart=std::min(2,tile_ncols-2);
  int xend=std::max(tile_nrows-3,1);
  int yend=std::max(tile_ncols-3,1);
  for (int i=0; i<tile_nrows; i++) {
    // sky
    z_sky_centre(i,0) = z_sky_centre(i, 1)*2 - z_sky_centre(i, ystart);
    z_sky_centre(i,tile_ncols-1) = z_sky_centre(i, tile_ncols-2)*2 - z_sky_centre(i, yend);
    // skyRMS
    z_skyRMS_centre(i,0) = z_skyRMS_centre(i, 1)*2 - z_skyRMS_centre(i, ystart);
    z_skyRMS_centre(i,tile_ncols-1) = z_skyRMS_centre(i, tile_ncols-2)*2 - z_skyRMS_centre(i, yend);
  }
  for (int i=0; i<tile_ncols; i++) {
    // sky
    z_sky_centre(0, i) = z_sky_centre(1, i)*2 - z_sky_centre(xstart, i);
    z_sky_centre(tile_nrows-1, i) = z_sky_centre(tile_nrows-2, i)*2 - z_sky_centre(xend, i);
    // skyRMS
    z_skyRMS_centre(0, i) = z_skyRMS_centre(1, i)*2 - z_skyRMS_centre(xstart, i);
    z_skyRMS_centre(tile_nrows-1, i) = z_skyRMS_centre(tile_nrows-2, i)*2 - z_skyRMS_centre(xend, i);
  }

  // Now interpolate for each image cell
  // Rcpp::Rcout << " AKIMA!\n";
  
  switch (type) {
  case adacs_CLASSIC_BILINEAR:
    interpolateLinearGrid(xseq, yseq, z_sky_centre, sky);
    interpolateLinearGrid(xseq, yseq, z_skyRMS_centre, skyRMS);
    break;
  case adacs_AKIMA_BICUBIC:
    interpolateAkimaGrid(xseq, yseq, z_sky_centre, sky);
    interpolateAkimaGrid(xseq, yseq, z_skyRMS_centre, skyRMS);
    break;
  }

  // Rcpp::Rcout << "POST AKIMA!\n";
  
  // Apply mask
  // if (mask.isNotNull()) {
  //   Rcpp::IntegerMatrix imask = Rcpp::as<Rcpp::IntegerMatrix>(mask);
  //   int nrows=image.nrow();
  //   int ncols=image.ncol();
  //   for (int i=0; i<ncols; i++) {
  //     for (int j=0; j<nrows; j++) {
  //       if (imask(j, i)==1) {
  //         sky(j, i) = NA_REAL;
  //         skyRMS(j, i) = NA_REAL;
  //       }
  //     }
  //   }
  // }
}

