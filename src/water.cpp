//
// Binding the watershed implementation to R using Rcpp
//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2018
// Copyright by UWA (in the framework of the ICRAR)
// All rights reserved
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston,
// MA 02111-1307  USA
//

#include <Rcpp.h>
#include "water.h"

// [[Rcpp::export]]
Rcpp::IntegerMatrix water_cpp(
    Rcpp::NumericVector image = 0, const int nx = 1, const int ny = 1,
    const double abstol = 1, const double reltol = 0, const double cliptol = 1000000,
    const int ext = 1, const double skycut = 0, const int pixcut = 1,
    const bool verbose = false, const int Ncheck = 1000000)
{
    // Use Rcpp for checking interruptions
    auto interrupt_checker = [Ncheck, verbose](const std::size_t i, const std::size_t total) {
        if ((i % Ncheck) == 0) {
            try {
                Rcpp::checkUserInterrupt();
            }
            catch (Rcpp::internal::InterruptedException& e)
            {
                if (verbose) {
                    Rcpp::Rcout << "Caught an interrupt!\n";
                }
                return true;
            }
            if (verbose == true && i > 0) {
                Rcpp::Rcout << "  - Segmented pixel " << i << " out of " << total << '\n';
            }
        }
        return false;
    };

    // Add 1 to segments so first segment is 1 and no-segments are 0
    // (since profound::watershed uses 0 and -1 respectively)
    Rcpp::IntegerMatrix segments(nx, ny);
    profound::watershed(&(image[0]), &(segments[0]), nx, ny, ext, abstol, reltol, cliptol, skycut, pixcut, interrupt_checker);
    for (R_xlen_t i = 0; i != segments.size(); i++) {
        segments[i]++;
    }
    return segments;
}
