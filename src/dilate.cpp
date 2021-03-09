#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".dilate_cpp")]]
IntegerMatrix dilate_cpp(IntegerMatrix segim, IntegerMatrix kern, IntegerVector expand=0){
  
  int srow = segim.nrow();
  int scol = segim.ncol();
  int krow = kern.nrow();
  int kcol = kern.ncol();
  int krow_off = ((krow - 1) / 2);
  int kcol_off = ((kcol - 1) / 2);
  bool checkseg = true;
  int max_segim = max(segim);
  IntegerMatrix segim_new(srow, scol);
  LogicalVector seglogic  (max_segim);
  
  // if(expand(0) == 0){
  //   expandall = true;
  // }else{
  //   expandall = false;
  //   
  // }
  
  if(expand(0) > 0){
    for(int k = 0; k < expand.length(); k++){
      if(expand(k) <= max_segim){
        seglogic(expand(k)-1) = true;
      }
    }
  }
  
  for (int i = 0; i < srow; i++) {
    for (int j = 0; j < scol; j++) {
      if(segim(i,j) > 0){
        if(expand(0) > 0){
          checkseg = seglogic(segim(i,j)-1);
        }
        if(checkseg){
          for (int m = std::max(0, krow_off - i); m < std::min(krow, krow_off - (i - srow)); m++) {
            for (int n = std::max(0, kcol_off - j); n < std::min(kcol, kcol_off - (j - scol)); n++) {
              if(kern(m,n) > 0){
                if(m != krow_off || n != kcol_off){
                  if(segim(i + m - krow_off,j + n - kcol_off) == 0){
                    int xloc = i + m - krow_off;
                    int yloc = j + n - kcol_off;
                    if(segim(i,j) < segim_new(xloc,yloc) || segim_new(xloc,yloc) == 0) {
                      segim_new(xloc,yloc) = segim(i,j);
                    }
                  }
                }else{
                  segim_new(i,j) = segim(i,j);
                }
              }
            }
          }
        }else{
          segim_new(i,j) = segim(i,j);
        }
      }
    }
  }
  return segim_new;
}

// IntegerMatrix dilate_cpp_old(IntegerMatrix segim, IntegerMatrix kern){
//   
//   int srow = segim.nrow();
//   int scol = segim.ncol();
//   int krow = kern.nrow();
//   int kcol = kern.ncol();
//   int krow_off = ((krow - 1) / 2);
//   int kcol_off = ((kcol - 1) / 2);
//   int maxint = std::numeric_limits<int>::max();
//   IntegerMatrix segim_new(srow, scol);
//   
//   for (int j = 0; j < scol; j++) {
//     for (int i = 0; i < srow; i++) {
//       int segID = maxint;
//       if(segim(i,j) == 0){
//         for (int n = std::max(0,kcol_off - j); n < std::min(kcol, kcol_off - (j - scol)); n++) {
//           for (int m = std::max(0,krow_off - i); m < std::min(krow, krow_off - (i - srow)); m++) {
//             if(kern(m,n) > 0){
//               int segim_segID = segim(i + m - krow_off,j + n - kcol_off);
//               if(segim_segID > 0 & segim_segID < segID) {
//                 segID = segim_segID;
//               }
//             }
//           }
//         }
//         if(segID < maxint){
//           segim_new(i,j) = segID;
//         }
//       }else{
//         segim_new(i,j) = segim(i,j);
//       }
//     }
//   }
//   
//   return segim_new;
// }
