#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <unistd.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix dilate_cpp(IntegerMatrix segim, IntegerMatrix kern){
  
  IntegerMatrix segim_new(segim.nrow(), segim.ncol());
  
  int srow = segim.nrow();
  int scol = segim.ncol();
  int krow = kern.nrow();
  int kcol = kern.ncol();
  int krow_off = ((krow - 1) / 2);
  int kcol_off = ((kcol - 1) / 2);
  
  for (int i = 0; i < srow; i++) {
    for (int j = 0; j < scol; j++) {
      int segID = 0;
      if(segim(i,j) == 0){
        for (int m = std::max(0,krow_off - i); m < std::min(krow, krow_off - (i - srow)); m++) {
          for (int n = std::max(0,kcol_off - j); n < std::min(kcol, krow_off - (j - scol)); n++) {
            if(kern(m,n) > 0){
              int xloc = i + m - krow_off;
              int yloc = j + n - kcol_off;
              if(segim(xloc, yloc) > 0) {
                if(segim(xloc,yloc) < segID | segID==0) {
                  segID = segim(xloc,yloc);
                }
              }
            }
          }
        }
      }
      if(segID > 0){
        segim_new(i,j) = segID;
      }else{
        segim_new(i,j) = segim(i,j);
      }
    }
  }
  
  return segim_new;
}
