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

    // Take positions of pixels above the skycut. Counting first lets us size
    // the vector exactly: the previous fixed size/10 guess is wrong for low
    // skycuts and the vector was left to grow (a ~2x transient allocation on
    // large, dense images). The extra scan is cheap relative to the sort.
    std::size_t count = 0;
    for (std::size_t i = 0; i < size; ++i) {
        if (image[i] > skycut) {
            ++count;
        }
    }

    std::vector<pix_idx> valid_pixels;
    valid_pixels.reserve(count);
    for (std::size_t i = 0; i < size; ++i) {
        if (image[i] > skycut) {
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
        const std::size_t nneigh = 2 * (std::size_t)ext + 1;
        merger_candidates.reserve(nneigh * nneigh - 1);
    }

    const double *image;
    int *segments;
    const std::size_t width;
    const std::size_t height;
    const std::size_t size;
    std::vector<std::size_t> relevant_indices {};
    std::vector<int> merger_candidates {};
    std::vector<int> seg_max_i {};
    // Union-find over segment ids. Merging two segments is an O(1) (amortised)
    // parent update; ids stored in `segments` may then be non-canonical and are
    // resolved with find_segment() when read. This replaces the previous
    // "relabel every pixel of the segment" loop, which made dense watersheds
    // quadratic in the number of pixels.
    std::vector<int> seg_parent {};
    const double abstol;
    const double reltol;
    const double cliptol;
    unsigned int segment_id = 0;

    int find_segment(int segment)
    {
        // path halving
        while (seg_parent[segment] != segment) {
            segment = seg_parent[segment] = seg_parent[seg_parent[segment]];
        }
        return segment;
    }

    void union_segments(int child, int parent)
    {
        seg_parent[find_segment(child)] = find_segment(parent);
    }

    bool within_merge_tolerance(int segment, double central_pixel) const
    {
        double pixel = image[relevant_indices[seg_max_i[segment]]];
        if (central_pixel > cliptol) {
            return true;
        }
        // reltol==0 is the default; pow(x, 0) is exactly 1, so skip it.
        if (reltol == 0.0) {
            return pixel - central_pixel < abstol;
        }
        return pixel - central_pixel < abstol * std::pow(pixel / central_pixel, reltol);
    }

    void apply_pixcut(double pixcut)
    {
        if (pixcut <= 1) {
            return;
        }
        // Canonicalise ids first so that the tabulation is by final segment
        int max_id = -1;
        for (std::size_t i = 0; i != size; i++) {
            int segment = segments[i];
            if (segment != NO_SEGMENT) {
                segment = segments[i] = find_segment(segment);
                if (segment > max_id) {
                    max_id = segment;
                }
            }
        }
        if (max_id < 0) {
            return;
        }
        std::vector<int> segment_count = tabulate(segments, size, max_id);
        for (std::size_t i = 0; i != size; i++) {
            int segment = segments[i];
            if (segment >= 0 && segment_count[segment] < pixcut) {
                segments[i] = NO_SEGMENT;
            }
        }
    }
};

static inline
void merge_segments(Problem &p, double central_pixel)
{
    // are there at least two unique segments flagged? Sorting a small vector
    // in place is much cheaper than building a std::set (which allocates a
    // node per element).
    auto &candidates = p.merger_candidates;
    std::sort(candidates.begin(), candidates.end());
    const auto unique_end = std::unique(candidates.begin(), candidates.end());
    const std::size_t n_unique = static_cast<std::size_t>(unique_end - candidates.begin());
    if (n_unique < 2) {
        return;
    }

    // first element (brightest) will be what is merged into the rest, if they
    // pass the test
    auto lowest_segment = candidates[0];
    for (std::size_t m = 1; m < n_unique; m++) {
        auto segment = candidates[m];
        if (!p.within_merge_tolerance(segment, central_pixel)) {
            continue;
        }
        // point the merged segment at the brightest peak flux segment
        p.union_segments(segment, lowest_segment);
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
            segment = p.find_segment(segment);
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
        merge_segments(p, center_pixel);
    }

    // if nothing has a segment value in the surrounding pixels then create a new segment seed
    if (p.segments[center_idx] == NO_SEGMENT) {
        p.segments[center_idx] = p.segment_id;
        p.seg_max_i.push_back(i);
        p.seg_parent.push_back(p.segment_id);
        p.segment_id++;
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

    // Final cut by pixel count cut; canonicalise the segment ids that were
    // lazily left pointing at merged-away segments.
    for (std::size_t i = 0; i != p.size; i++) {
        int segment = p.segments[i];
        if (segment != NO_SEGMENT) {
            p.segments[i] = p.find_segment(segment);
        }
    }
    p.apply_pixcut(pixcut);
}

}  // namespace profound