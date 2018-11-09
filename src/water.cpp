#include <Rcpp.h>
#include <cmath>
#include <unistd.h>
using namespace Rcpp;

// [[Rcpp::export(".order_cpp")]]
IntegerVector order_cpp(NumericVector x) {
  // NumericVector sorted = clone(x).sort(true);
  return match(clone(x).sort(true), x) - 1;
}
  
// [[Rcpp::export(".tabulate_cpp")]]
IntegerVector tabulate_cpp(const IntegerVector& x, const int max) {
    IntegerVector counts(max);
    int lim = x.size();
    for (int i = 0; i < lim; i++) {
        if (x[i] > 0 && x[i] <= max)
            counts[x[i] - 1]++;
    }
    return counts;
}

// [[Rcpp::export]]
IntegerVector water_cpp(const NumericVector image = 0, const int nx = 1, const int ny = 1,
                        const double abstol = 1, const double reltol = 0, const int ext = 1, 
                        const double skycut = 0, const int pixcut = 1, const bool verbose = false,
                        const int Ncheck = 1000000){
  //sanity check
  if(image.size() != nx*ny){
    stop("image size does not equate to nx.ny!");
  }
  
  Function order_int("order");
  
  // check required vector size
  int size_count = 0;
  int ilim = image.size();
  for (int i = 0; i < ilim; ++i) {
    if(image[i]>= skycut){
      ++size_count;
    }
  }
  // assign and fill up with IDs
  // checking the required size first is faster than making it too big and cutting down
  IntegerVector imvec(size_count,-1);
  size_count = 0;
  for (int i = 0; i < ilim; ++i) {
    if(image[i]>= skycut){
      imvec(size_count) = i;
      ++size_count;
    }
  }
  
  // constants
  const int Nmerge = (ext*2+1)*(ext*2+1)-1;

  // variables
  IntegerVector imord =  order_int(image[imvec], Named("decreasing",true), Named("na.last",NA_REAL));
  imvec = imvec[imord-1]; // here we order the pixels, and therefore assess brightest first
  IntegerMatrix segim(nx,ny); // the segim map we want at the end
  // IntegerVector seg_max_pos(imvec.size()/pixcut); // containing for maximum pixel flux for a segment
  IntegerVector seg_max_i(imvec.size()/pixcut); // i offset
  IntegerVector segmerge(Nmerge); // empty vector of segments that might need merging
  
  int seg_id = 0; // current segment ID
  int segmerge_count = 0; // counter for merging
  bool merge_flag = false;
  
  // loop over important pixels only
  ilim = imvec.size();
  for (int i = 0; i < ilim; ++i) {
    if((i % Ncheck)==0){
      try{
        Rcpp::checkUserInterrupt();
      }
      catch (Rcpp::internal::InterruptedException& e)
      {
        Rcout << "Caught an interrupt!" << std::endl;
        i=ilim-1;
      }
      if(verbose == true){
        Rcout << "  - Segmented pixel " << i << " out of " << ilim-1 << std::endl;
      }
    }
    segmerge_count=0;
    // reset segmerge vector
    segmerge = rep(0, Nmerge);
    // reset merge_flag
    merge_flag = false;
    int current_pos = imvec[i];
    int x_current = current_pos % nx; // current x position in image
    int y_current = current_pos / nx; // current y position in image
    int comp_pos = current_pos;
    for (int j = -ext; j <= ext; ++j) {
      for (int k = -ext; k <= ext; ++k) {
        if(j == 0 && k == 0){
          ++k;
        }
        int x_offset = x_current + j; // current x position in image
        int y_offset = y_current + k; // current y position in image
        // check we are not at the edge of the image
        if((x_offset >= 0) & (x_offset < nx) & (y_offset >= 0) & (y_offset < ny)) {
          // apply conversion to absolute pixel ref
          int offset_pos = x_offset + y_offset*nx;
          int offset_seg = segim[offset_pos];
          if(offset_seg > 0){
            // if the offset position if brighter consider doing something
            if(image[offset_pos] > image[comp_pos]){
              // if the brightest pixel in the offset segment is not brighter than the abstol flag for merging
              segmerge[segmerge_count] = offset_seg;
              segmerge_count++;
              double imref = image[imvec[seg_max_i[offset_seg-1]]];
              if(imref - image[current_pos] < abstol * pow(imref/image[current_pos],reltol)){
                merge_flag=true;
              }
              comp_pos = offset_pos;
              segim[current_pos] = offset_seg;
            }
          }
        }
      }
    }

    // is there anything to consider merging?
    if(merge_flag){
      if(is_true(any(segmerge>0)) && abstol > 0){
        segmerge=segmerge[segmerge>0];
        segmerge=sort_unique(segmerge);
        // are there at least two unique segments flagged?
        if(segmerge.size()>1){
          // looper over flagged segments that are not the brightest
          for (int m = 1; m < segmerge.size(); ++m) {
            // loop over pixels segmented to date
            int segcheck = segmerge[m];
            double imref = image[imvec[seg_max_i[segcheck-1]]];
            if(imref - image[current_pos] < abstol * pow(imref/image[current_pos],reltol)){
              for (int n = seg_max_i[segcheck-1]; n <= i; ++n) {
              // for (int n = 0; n <= i; ++n) {
                int merge_loc = imvec[n];
                // if pixel is flagged for merging, set to lowest segment value (brightest peak flux segment)
                if(segim[merge_loc] == segcheck){
                  segim[merge_loc] = segmerge[0];
                }
              }
            }
          }
        }
      }
    }

    // if nothing has a segment value in the surrounding pixels then create a new segment seed
    if(segim[current_pos] == 0){
      seg_id++;
      segim[current_pos] = seg_id;
      seg_max_i(seg_id-1)=i;
      //seg_max_pos(seg_id-1) = current_pos;
    }
  }
  
  if(pixcut>1){
    IntegerVector tabulate_seg = tabulate_cpp(segim, max(segim));
    ilim = segim.size();
    for (int i = 0; i < ilim; ++i) {
      if(segim[i] > 0){
        if(tabulate_seg[segim[i]-1] < pixcut){
          segim[i] = 0;
        }
      }
    }
  }
  return segim;
  //return imvec;
}

// Rcout << "comment=" << output << std::endl;
