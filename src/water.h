//
// Watershed implementation using C++11
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

#include <algorithm>
#include <cmath>
#include <functional>
#include <set>
#include <vector>

namespace profound {

/// The value signaling that no segment is present
const int NO_SEGMENT = -1;

/**
 * Given an image, filter pixels that are above the skycut, sort them from
 * brightest to faintest, and return the positions of those sorted pixels in
 * the original image
 */
static inline
std::vector<std::size_t> get_sorted_indices(const double *image, std::size_t size, double skycut)
{
    struct pix_idx {
        double pix;
        std::size_t idx;
        bool operator<(const pix_idx &rhs) const
        {
            return pix > rhs.pix;
        }
    };

    // Take positions of pixels above the skycut. We estimate that ~10%
    // of the pixels will meet the condition and try to reserve that much memory
    // at once
    std::vector<pix_idx> valid_pixels;
    valid_pixels.reserve(size / 10.);
    for (std::size_t i = 0; i < size; ++i) {
        if (image[i] >= skycut) {
            valid_pixels.push_back({image[i], i});
        }
    }

    // Stable-sort by decreasing pixel value, return only the indices
    std::stable_sort(valid_pixels.begin(), valid_pixels.end());

    std::vector<std::size_t> indices(valid_pixels.size());
    std::transform(valid_pixels.begin(), valid_pixels.end(), indices.begin(),
        [](const pix_idx &p) {
            return p.idx;
        });
    return indices;
}

static inline
std::vector<int> tabulate(const int *segments, std::size_t n, int max)
{
    std::vector<int> counts(max + 1);
    for (std::size_t i = 0; i != n; i++) {
        int segment = segments[i];
        if (segment >= 0) {
            counts[segment]++;
        }
    }
    return counts;
}

static inline
void apply_pixcut(int *segments, std::size_t size, double pixcut)
{
    if (pixcut <= 1) {
        return;
    }
    auto max = *std::max_element(segments, segments + size);
    if (max < 0) {
        return;
    }
    std::vector<int> segment_count = tabulate(segments, size, max);
    for (std::size_t i = 0; i != size; i++) {
        int segment = segments[i];
        if (segment >= 0 && segment_count[segment] < pixcut) {
            segments[i] = NO_SEGMENT;
        }
    }
}

/**
 * The watershedding problem. Inputs are an image with certain width and height
 * and some tolerance values.
 */
struct Problem {

    Problem(const double *image, int *segments, unsigned int width, unsigned int height,
        unsigned int ext, double abstol, double reltol, double cliptol, double skycut) :
        image(image), segments(segments),
        width(width), height(height), size(width * height),
        relevant_indices(get_sorted_indices(image, size, skycut)),
        abstol(abstol), reltol(reltol), cliptol(cliptol)
    {
        std::fill(segments, segments + size, NO_SEGMENT);
        merger_candidates.reserve(ext * ext * 4 - 1);
    }

    const double *image;
    int *segments;
    const std::size_t width;
    const std::size_t height;
    const std::size_t size;
    std::vector<std::size_t> relevant_indices {};
    std::vector<int> merger_candidates {};
    std::vector<int> seg_max_i {};
    const double abstol;
    const double reltol;
    const double cliptol;
    unsigned int segment_id = 0;

    bool within_merge_tolerance(int segment, double central_pixel) const
    {
        double pixel = image[relevant_indices[seg_max_i[segment]]];
        return central_pixel > cliptol ||
               pixel - central_pixel < abstol * std::pow(pixel / central_pixel, reltol);
    }
};

static inline
void merge_segments(Problem &p, int i, double central_pixel)
{
    // are there at least two unique segments flagged?
    std::set<int> mergers(p.merger_candidates.begin(), p.merger_candidates.end());
    if (mergers.size() < 2) {
        return;
    }

    // first element (brightest) will be what is merged into the rest, if they
    // pass the test
    auto it = mergers.begin();
    auto lowest_segment = *it;
    for (it++; it != mergers.end(); it++) {
        auto segment = *it;
        if (!p.within_merge_tolerance(segment, central_pixel)) {
            continue;
        }
        // loop round segments that have been allocated already
        // if pixel is flagged for merging, set to lowest segment value (brightest peak flux segment)
        for (int n = p.seg_max_i[segment]; n <= i; ++n) {
            auto merge_idx = p.relevant_indices[n];
            if (p.segments[merge_idx] == segment) {
                p.segments[merge_idx] = lowest_segment;
            }
        }
    }
}

static inline
void watershed_cetered_at(Problem &p, int i, int ext)
{
    std::size_t center_idx = p.relevant_indices[i];
    int x = center_idx % p.width;
    int y = center_idx / p.width;
    const double center_pixel = p.image[center_idx];

    // the brightest pixel we have seen so far in this area
    double brightest_pixel = center_pixel;

    bool merge = false;
    p.merger_candidates.clear();

    // Consider valid surrounding pixels; central pixel is skipped
    // Mind the looping order: x changes faster, so we should have nicer memory
    // access patterns
    for (int j = -ext; j <= ext; ++j) {
        for (int k = -ext; k <= ext; ++k) {
            if (j == 0 && k == 0) {
                ++k;
            }
            int off_x = x + k;
            int off_y = y + j;
            if (off_x < 0 || (unsigned)off_x >= p.width || off_y < 0 || (unsigned)off_y >= p.height) {
                continue;
            }

            int off_idx = off_x + off_y * p.width;

            // If segment exists, it will be considered for merging
            int segment = p.segments[off_idx];
            if (segment == NO_SEGMENT) {
                continue;
            }
            p.merger_candidates.push_back(segment);

            // do we actually need to perform merging later?
            if (!merge && p.within_merge_tolerance(segment, center_pixel)) {
                merge = true;
            }

            // Existing segment is brighter than our brightest pixel (and thus
            // our center), update both
            auto off_pixel = p.image[off_idx];
            if (off_pixel <= brightest_pixel) {
                continue;
            }
            brightest_pixel = off_pixel;
            p.segments[center_idx] = segment;

        }
    }

    if (merge && p.abstol > 0) {
        merge_segments(p, i, center_pixel);
    }

    // if nothing has a segment value in the surrounding pixels then create a new segment seed
    if (p.segments[center_idx] == NO_SEGMENT) {
        p.segments[center_idx] = p.segment_id++;
        p.seg_max_i.push_back(i);
    }
}

template <typename InterruptChecker>
void watershed(
    const double *image, int *segments, const int nx, const int ny, const int ext,
    const double abstol, const double reltol, const double cliptol,
    const double skycut, const int pixcut, InterruptChecker &&interrupt_checker)
{
    // Prepare all structures, etc
    Problem p(image, segments, nx, ny, ext, abstol, reltol, cliptol, skycut);

    // Loop over relevant indices (those above skycut, sorted by desc pixel value)
    auto n_indices = p.relevant_indices.size();
    for (std::size_t i = 0; i < n_indices; ++i) {
        if (interrupt_checker(i, n_indices)) {;
            break;
        }
        watershed_cetered_at(p, i, ext);
    }

    // Final cut by pixel count cut and add 1 to segment image (so first segment is 1 and not 0)
    apply_pixcut(p.segments, p.size, pixcut);
}

}  // namespace profound