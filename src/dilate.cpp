#include <Rcpp.h>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;

// [[Rcpp::export(".dilate_cpp")]]
IntegerMatrix dilate_cpp(IntegerMatrix segim, IntegerMatrix kern, IntegerVector expand=0, int nthreads = 1){

  const int srow = segim.nrow();
  const int scol = segim.ncol();
  const int krow = kern.nrow();
  const int kcol = kern.ncol();
  const int krow_off = ((krow - 1) / 2);
  const int kcol_off = ((kcol - 1) / 2);
  const int max_segim = max(segim);
  IntegerMatrix segim_new(srow, scol);

  const int* s = INTEGER(segim);
  const int* k = INTEGER(kern);
  int* sn = INTEGER(segim_new);

  // expand lookup table (1-based segment ids)
  const bool use_expand = expand.length() > 0 && expand(0) > 0;
  std::vector<unsigned char> in_expand(use_expand ? (size_t)std::max(max_segim, 0) + 1 : 0, 0);
  if(use_expand){
    for(int e = 0; e < expand.length(); e++){
      const int v = expand(e);
      if(v > 0 && v <= max_segim){
        in_expand[v] = 1;
      }
    }
  }

  // Precompute the non-zero, non-centre kernel offsets. R matrices are
  // column-major, so element (m,n) lives at m + n*krow. Skipping the zero
  // cells of the kernel avoids repeated tests inside the innermost loop;
  // the result is order independent (each target keeps the running minimum
  // segment id seen), so the cell order does not matter.
  const bool centre_on = k[krow_off + (size_t)kcol_off * krow] > 0;
  std::vector<int> kdm, kdn, koff;
  kdm.reserve((size_t)krow * kcol);
  kdn.reserve((size_t)krow * kcol);
  koff.reserve((size_t)krow * kcol);
  int maxdm = 0, maxdn = 0;
  for (int n = 0; n < kcol; n++) {
    for (int m = 0; m < krow; m++) {
      if(k[m + (size_t)n * krow] > 0 && (m != krow_off || n != kcol_off)){
        const int dm = m - krow_off;
        const int dn = n - kcol_off;
        kdm.push_back(dm);
        kdn.push_back(dn);
        koff.push_back(dm + dn * srow);
        if(dm > maxdm) maxdm = dm;
        if(-dm > maxdm) maxdm = -dm;
        if(dn > maxdn) maxdn = dn;
        if(-dn > maxdn) maxdn = -dn;
      }
    }
  }
  const size_t nkoff = kdm.size();
  const int* kdm_p = kdm.data();
  const int* kdn_p = kdn.data();
  const int* koff_p = koff.data();
  const long long N = (long long)srow * scol;
  const bool can_be_interior = srow > 2 * maxdm && scol > 2 * maxdn;

  // Scan the image contiguously in memory (column-major) and only touch the
  // kernel footprint for pixels that actually carry a segment. Targets only
  // ever receive values from source pixels, so the two-phase structure
  // (scan for sources, then scatter) is equivalent to the original scan.
#ifdef _OPENMP
  // Parallelize the main loop
#pragma omp parallel for schedule(static) num_threads(nthreads)
#endif
  for (long long idx = 0; idx < N; idx++) {
    const int srcval = s[idx];
    if(srcval <= 0){
      continue;
    }
    if(use_expand && in_expand[srcval] == 0){
      sn[idx] = srcval;
      continue;
    }
    const int i = (int)(idx % srow);
    const int j = (int)(idx / srow);
    if(centre_on){
      sn[idx] = srcval;
    }
    if(can_be_interior && i >= maxdm && i < srow - maxdm && j >= maxdn && j < scol - maxdn){
      for (size_t t = 0; t < nkoff; t++) {
        const long long tidx = idx + koff_p[t];
        if(s[tidx] != 0){
          continue;
        }
        const int cur = sn[tidx];
        if(cur == 0 || srcval < cur){
          sn[tidx] = srcval;
        }
      }
    } else {
      for (size_t t = 0; t < nkoff; t++) {
        const int xloc = i + kdm_p[t];
        if(xloc < 0 || xloc >= srow){
          continue;
        }
        const int yloc = j + kdn_p[t];
        if(yloc < 0 || yloc >= scol){
          continue;
        }
        const long long tidx = xloc + (long long)yloc * srow;
        if(s[tidx] != 0){
          continue;
        }
        const int cur = sn[tidx];
        if(cur == 0 || srcval < cur){
          sn[tidx] = srcval;
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
