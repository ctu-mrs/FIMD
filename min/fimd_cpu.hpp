/**
 * \file fimd_cpu.hpp
 * \author Vojtech Vrba (vrba.vojtech [at] fel.cvut.cz)
 * \date January 2025
 * \brief Compile-time Bresenham's circle generation and FIMD-CPU implementation (requires C++20).
 * \copyright GNU Public License.
 */

#ifndef FIMD_CPU_HPP
#define FIMD_CPU_HPP

#include <array>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <list>
#include <type_traits>


namespace fimd {

/**
 * \brief Alias for 2D points on an integer grid.
 */
using Point2D = std::array<int, 2>;

/**
 * \brief Computes the integer square root of a number at compile-time.
 * \param x The number to compute the square root of.
 * \param lo The lower bound of the search range.
 * \param hi The upper bound of the search range.
 * \return The integer square root of x.
 */
static consteval unsigned sqrt_recursive(const unsigned x, const unsigned lo, const unsigned hi) {
    const unsigned mid = (lo + hi + 1) / 2;
    return (lo == hi) ? lo : ((x / mid < mid) ? sqrt_recursive(x, lo, mid - 1) : sqrt_recursive(x, mid, hi));
}


/**
 * \brief Compile-time estimation of the Bresenham boundary length.
 * \tparam R Radius of the circle.
 */
template<int R> requires (R > 0)
constexpr unsigned BresenhamBoundaryLengthEstimation = (R == 2) ? 4 : sqrt_recursive(2 * R * R, 0, R * R + 1) + 1;


/**
 * \brief Represents a boundary point in the first quadrant of a Bresenham's circle.
 * \tparam R Radius of the circle.
 * \tparam I Index of the point.
 * \tparam S Sequential order index.
 * \tparam X X-coordinate of the point.
 * \tparam Y Y-coordinate of the point.
 * \tparam P Decision parameter.
 * \tparam O Octant of the point.
 */
template<unsigned R, unsigned I, unsigned S, int X, int Y, int P, unsigned O>
class BresenhamBoundaryPoint : public Point2D {
private:
    static constexpr unsigned NextS = (O == 0) ? (BresenhamBoundaryLengthEstimation<R> - S - 1) : (BresenhamBoundaryLengthEstimation<R> - S);
    static constexpr Point2D XY = (O == 0) ? Point2D{X, Y} : Point2D{Y, X};

    using NextBresenhamBoundaryPoint = const std::conditional_t<(O == 1),
        std::conditional_t<(P < 0),
            BresenhamBoundaryPoint<R, I+1, NextS, X+1, Y, P+4*(X+1)+6, 0>,
            BresenhamBoundaryPoint<R, I+1, NextS, X+1, Y-1, P+4*((X+1)-(Y-1))+10, 0>
        >,
        BresenhamBoundaryPoint<R, I+1, NextS, X, Y, P, 1>
    >;
    static constexpr auto next = NextBresenhamBoundaryPoint();
public:
    constexpr BresenhamBoundaryPoint() : Point2D(XY) {};
    constexpr NextBresenhamBoundaryPoint get_next() { return next; };
};


/**
 * \brief Specialization of BresenhamBoundaryPoint for termination condition.
 * \tparam R Radius of the circle.
 * \tparam I Index of the point.
 * \tparam S Sequential order index.
 * \tparam X X-coordinate of the point.
 * \tparam Y Y-coordinate of the point.
 * \tparam P Decision parameter.
 * \tparam O Octant of the point.
 */
template<unsigned R, unsigned I, unsigned S, int X, int Y, int P, unsigned O> requires (X > Y or (X == Y and O == 1))
class BresenhamBoundaryPoint<R, I, S, X, Y, P, O> : public Point2D {
public:
    constexpr BresenhamBoundaryPoint() : Point2D({-1, -1}) {};
    constexpr BresenhamBoundaryPoint get_next() { return *this; };
};


/**
 * \brief Counts the number of boundary points in a Bresenham's circle quadrant.
 * \tparam R Radius of the circle.
 * \tparam I Index of the point.
 * \tparam S Sequential order index.
 * \tparam X X-coordinate of the point.
 * \tparam Y Y-coordinate of the point.
 * \tparam P Decision parameter.
 * \tparam O Octant of the point.
 * \param point The initial boundary point.
 * \return The number of boundary points.
 */
template<unsigned R, unsigned I, unsigned S, int X, int Y, int P, unsigned O>
consteval unsigned BresenhamBoundaryCounter(BresenhamBoundaryPoint<R, I, S, X, Y, P, O> point) {
    if (point != point.get_next()) {
        return BresenhamBoundaryCounter(point.get_next());
    }
    return I;
}


/**
 * \brief Alias for the initial boundary point in Bresenham's circle algorithm.
 * \tparam R Radius of the circle.
 */
template<int R> requires (R > 0)
using BresenhamBoundaryInit = BresenhamBoundaryPoint<R, 0, 0, 0, R, 3-2*R, 0>;


/**
 * \brief Compile-time calculation of the Bresenham boundary length.
 * \tparam R Radius of the circle.
 */
template<int R> requires (R > 0)
constexpr unsigned BresenhamBoundaryLength = BresenhamBoundaryCounter(BresenhamBoundaryInit<R>());


/**
 * \brief Represents the boundary of a circle using Bresenham's algorithm.
 * \tparam R Radius of the circle.
 */
template<int R> requires (R > 0)
class BresenhamBoundary : public std::array<Point2D, 4*(BresenhamBoundaryLength<R>-1)>
{
public:
    consteval BresenhamBoundary() {
        seq_sort_recursive(*this, BresenhamBoundaryInit<R>());
        eval_sort(*this);
    };

    static std::ostream& print_verbose(std::ostream& os) {
        os << "[Bresenham boundary for R=" << radius() << "]" << std::endl;
        os << "Quadrant length: " << length_quadrant() << " (estimated: " << length_quadrant_estimation() << ")" << std::endl;
        os << "Total boundary length: " << length() << std::endl;
        os << "Compile-time generation:" << std::endl;
        return print_recursive(os, BresenhamBoundaryInit<R>());
    }

    static std::ostream& print_sorted(std::ostream& os, const BresenhamBoundary boundary) {
        os << "Sorted points:" << std::endl << "[";
        for (Point2D point : boundary) {
            os << "(" << point[0] << "," << point[1] << "), ";
        }
        os << "\b\b]" << std::endl;
        return os;
    }

    static constexpr unsigned length_quadrant() {
        return BresenhamBoundaryLength<R>;
    }

    static constexpr unsigned length() {
        return 4 * (length_quadrant() - 1);
    }

    static constexpr unsigned length_quadrant_estimation() {
        return BresenhamBoundaryLengthEstimation<R>;
    }

    static constexpr int radius() {
        return R;
    }

private:
    template<unsigned R2, unsigned I, unsigned S, int X, int Y, int P, unsigned O>
    static constexpr std::ostream& print_recursive(std::ostream& os, BresenhamBoundaryPoint<R2, I, S, X, Y, P, O> point) {
        if (point != point.get_next()) {
            os << "-- Point #" << (I+1) << ": (" << point[0] << ", " << point[1] << ") [S=" << S << ", X=" << X << ", Y=" << Y << ", P=" << P << ", O=" << O << "]" << std::endl;
            return print_recursive(os, point.get_next());
        }
        return os;
    };

    template<unsigned R2, unsigned I, unsigned S, int X, int Y, int P, unsigned O>
    static consteval void seq_sort_recursive(std::array<Point2D, length()> &points, BresenhamBoundaryPoint<R2, I, S, X, Y, P, O> point) {
        if (point != point.get_next()) {
            constexpr unsigned index = (O == 1 and length_quadrant_estimation() > length_quadrant()) ? S - 1 : S;
            if (index < length_quadrant() - 1) {
                points[index] = {point[0], point[1]};
                points[index + (length_quadrant() - 1)] = {point[1], -point[0]};
                points[index + 2 * (length_quadrant() - 1)] = {-point[0], -point[1]};
                points[index + 3 * (length_quadrant() - 1)] = {-point[1], point[0]};
            }
            seq_sort_recursive(points, point.get_next());
        }
    };

    static consteval void eval_sort(std::array<Point2D, length()> &points) {
        auto eval = std::array<Point2D, length()>();
        auto dists = std::array<int, length_quadrant()>();
        for (int i = 0; i < length_quadrant(); i++) {
            dists[i] = std::max(R - points[i][0], points[i][1]);
        }

        int i_next = 0;
        unsigned cnt = 0;
        while (cnt < length()) {
            auto point = points[i_next];
            eval[cnt++] = {point[0], point[1]};
            eval[cnt++] = {-point[0], -point[1]};
            if (point[0] == 0) {
                eval[cnt++] = {point[1], 0};
                eval[cnt++] = {-point[1], 0};
            } else {
                eval[cnt++] = {-point[1], point[0]};
                eval[cnt++] = {point[1], -point[0]};
            }
            dists[i_next] = 0;

            int sel_d = 0;
            int sel_y = 0;
            int sel_i = 0;
            for (int i = 0; i < length_quadrant(); i++) {
                auto point2 = points[i];
                auto dist = std::max(std::abs(point2[0] - point[0]), std::abs(point2[1] - point[1]));
                if (dists[i] > dist) {
                    dists[i] = dist;
                }
                if (dists[i] > sel_d or (dists[i] == sel_d and point2[1] > sel_y)) {
                    sel_d = dists[i];
                    sel_y = point2[1];
                    sel_i = i;
                }
            }
            i_next = sel_i;
        }

        for (int i = 0; i < length(); i++) {
            points[i] = eval[i];
        }
    };
};


/**
 * \brief Represents an interior point in upper half of a Bresenham's circle.
 * \tparam I Index of the point.
 * \tparam X X-coordinate of the point.
 * \tparam Y Y-coordinate of the point.
 * \tparam P Decision parameter.
 * \tparam O Octant of the point.
 * \tparam YI Y-coordinate increment.
 */
template<unsigned I, int X, int Y, int P, unsigned O, int YI>
class BresenhamInteriorPoint : public Point2D {
    using NextBresenhamInteriorPoint = const std::conditional_t<(YI < Y-1),
        BresenhamInteriorPoint<I+1, X, Y, P, O, YI+1>,
        std::conditional_t<(O == 0),
            std::conditional_t<(X < YI) or (X == 0),
                BresenhamInteriorPoint<I+1, X, Y, P, 1, X+1>,
                BresenhamInteriorPoint<I+1, X, Y, P, 2, X>
            >,
            std::conditional_t<(O == 1) and (X > 0),
                BresenhamInteriorPoint<I+1, X, Y, P, 2, X>,
                std::conditional_t<(O == 2),
                    BresenhamInteriorPoint<I+1, X, Y, P, 3, X+1>,
                    std::conditional_t<(P < 0),
                        BresenhamInteriorPoint<I+1, X+1, Y, P+4*(X+1)+6, 0, X+1>,
                        BresenhamInteriorPoint<I+1, X+1, Y-1, P+4*((X+1)-(Y-1))+10, 0, X+1>
                    >
                >
            >
        >
    >;
    using XY = const std::conditional_t<(O == 0),
        std::integral_constant<Point2D, {X, YI}>,
        std::conditional_t<(O == 1),
            std::integral_constant<Point2D, {YI, X}>,
            std::conditional_t<(O == 2),
                std::integral_constant<Point2D, {-X, YI}>,
                std::integral_constant<Point2D, {-YI, X}>
            >
        >
    >;

    static constexpr auto next = NextBresenhamInteriorPoint();

public:
    constexpr BresenhamInteriorPoint() : Point2D(XY::value) {};
    constexpr NextBresenhamInteriorPoint get_next() { return next; };
};


/**
 * \brief Specialization of BresenhamInteriorPoint for termination condition.
 * \tparam I Index of the point.
 * \tparam X X-coordinate of the point.
 * \tparam Y Y-coordinate of the point.
 * \tparam P Decision parameter.
 * \tparam O Octant of the point.
 * \tparam YI Interior coordinate increment.
 */
template<unsigned I, int X, int Y, int P, unsigned O, int YI> requires (YI >= Y)
class BresenhamInteriorPoint<I, X, Y, P, O, YI> : public Point2D {
public:
    constexpr BresenhamInteriorPoint() : Point2D({-1, -1}) {};
    constexpr BresenhamInteriorPoint get_next() { return *this; };
};


/**
 * \brief Counts the number of interior points in a Bresenham's circle.
 * \tparam I Index of the point.
 * \tparam X X-coordinate of the point.
 * \tparam Y Y-coordinate of the point.
 * \tparam P Decision parameter.
 * \tparam O Octant of the point.
 * \tparam YI Interior coordinate increment.
 * \param point The current interior point.
 * \return The number of interior points.
 */
template<unsigned I, int X, int Y, int P, unsigned O, int YI>
consteval unsigned BresenhamInteriorCounter(BresenhamInteriorPoint<I, X, Y, P, O, YI> point) {
    if (point != point.get_next()) {
        return BresenhamInteriorCounter(point.get_next());
    }
    return I;
}


/**
 * \brief Alias for the initial interior point in Bresenham's circle algorithm.
 * \tparam R Radius of the circle.
 */
template<int R> requires (R > 0)
using BresenhamInteriorInit = BresenhamInteriorPoint<0, 0, R, 3-2*R, 0, 0>;


/**
 * \brief Compile-time calculation of the Bresenham interior length.
 * \tparam R Radius of the circle.
 */
template<int R> requires (R > 0)
constexpr unsigned BresenhamInteriorLength = BresenhamInteriorCounter(BresenhamInteriorInit<R>());


/**
 * \brief Represents the interior of a circle using Bresenham's algorithm.
 * \tparam R Radius of the circle.
 */
template<int R> requires (R > 0)
class BresenhamInterior : public std::array<Point2D, BresenhamInteriorLength<R>>
{
public:
    consteval BresenhamInterior() {
        copy_recursive(*this, BresenhamInteriorInit<R>());
    };

    static std::ostream& print_verbose(std::ostream& os) {
        os << "[Bresenham interior for R=" << radius() << "]" << std::endl;
        os << "Half area (length): " << length() << std::endl;
        os << "Compile-time generation:" << std::endl;
        return print_recursive(os, BresenhamInteriorInit<R>());
    }

    static constexpr unsigned length() {
        return BresenhamInteriorLength<R>;
    }

    static constexpr int radius() {
        return R;
    }

private:
    template<unsigned I, int X, int Y, int P, unsigned O, int YI>
    static constexpr std::ostream& print_recursive(std::ostream& os, BresenhamInteriorPoint<I, X, Y, P, O, YI> point) {
        if (point != point.get_next()) {
            os << "-- Point #" << (I+1) << ": (" << point[0] << ", " << point[1] << ") [X=" << X << ", Y=" << Y << ", P=" << P << ", O=" << O << ", YI=" << YI << "]" << std::endl;
            return print_recursive(os, point.get_next());
        }
        return os;
    };

    template<unsigned I, int X, int Y, int P, unsigned O, int YI>
    static consteval void copy_recursive(std::array<Point2D, length()> &points, BresenhamInteriorPoint<I, X, Y, P, O, YI> point) {
        if (point != point.get_next()) {
            points[I] = {point[0], point[1]};
            copy_recursive(points, point.get_next());
        }
    };
};



/**
 * \brief Converts a 2D point to a 1D coordinate.
 * \param point The 2D point.
 * \param im_width Image width.
 * \return The 1D coordinate.
 */
inline int coord2to1(const Point2D& point, const unsigned& im_width) { return (point[1] * im_width) + point[0]; };

/**
 * \brief Converts a 1D coordinate to a 2D point.
 * \param coord1d The 1D coordinate.
 * \param im_width Image width.
 * \return The 2D point.
 */
inline Point2D coord1to2(const size_t& coord1d, const unsigned& im_width) { return Point2D{static_cast<int>(coord1d % im_width), static_cast<int>(coord1d / im_width)}; };


/**
 * \brief Fast Isolated Marker Detection (FIMD) CPU implementation.
 * \tparam RADIUS Radius of the circle.
 * \tparam PIXEL Type of the pixel data (unsigned char by default).
 * \tparam TERM_SEQ Type defining the array of pixels used for the termination sequence (array of 2 PIXEL values by default).
 */
template<unsigned RADIUS, typename PIXEL=unsigned char, typename TERM_SEQ=std::array<PIXEL, 2>>
class FIMD_CPU {
public:
    /**
     * \brief FIMD-CPU constructor.
     * \param im_width Image width.
     * \param im_height Image height.
     * \param threshold_center Center (marker) threshold value.
     * \param threshold_diff Difference threshold value.
     * \param threshold_sun Sun threshold value.
     * \param termination Termination sequence (array of pixels).
     * \param max_markers_count Limit for the number of the detected markers (0 for no limit).
     * \param max_sun_points_count Limit for the number of the detected sun points (0 for no limit).
     * \param alloc_frame If true, the frame buffer will be allocated during initialization.
     */
    FIMD_CPU(unsigned im_width, unsigned im_height, PIXEL threshold_center=120, PIXEL threshold_diff=60, PIXEL threshold_sun=240, TERM_SEQ termination={0xFF, 0x00}, unsigned max_markers_count=0, unsigned max_sun_points_count=0, bool alloc_frame=true)
    : im_width_(im_width), im_height_(im_height), threshold_center_(threshold_center), threshold_diff_(threshold_diff), threshold_sun_(threshold_sun), termination_(termination)
    {
        offset_ = (im_width_ * RADIUS) + RADIUS;
        max_markers_count_ = (max_markers_count == 0) ? std::numeric_limits<unsigned>::max() : max_markers_count;
        max_sun_points_count_ = (max_sun_points_count == 0) ? std::numeric_limits<unsigned>::max() : max_sun_points_count;
        frame_ = alloc_frame ? static_cast<PIXEL*>(std::malloc(im_width_ * im_height_ * sizeof(PIXEL))) : nullptr;
    };

    ~FIMD_CPU() {
        if (frame_ != nullptr) {
            std::free(frame_);
        }
    };

    /**
     * \brief Detects markers and sun points in an image.
     * \param image The image data (array of pixels).
     * \param markers The list of detected markers (2D points).
     * \param sun_points The list of detected sun points (2D points).
     * \param make_copy If true, a copy of the input image will be used in the detection process.
     * \return The total number of processed pixels in the input image.
     */
    unsigned detect(PIXEL* image, std::list<Point2D> &markers, std::list<Point2D> &sun_points, bool make_copy=true) {
        if (im_width_ < (2*RADIUS+1) or im_height_ < (2*RADIUS+1)) {
            return 0;
        }

        PIXEL* target_image;

        if (make_copy) {
            if (frame_ == nullptr) {
                frame_ = static_cast<PIXEL*>(std::malloc(im_width_ * im_height_ * sizeof(PIXEL)));
            }
            std::copy_n(image, im_height_ * im_width_, frame_);
            target_image = frame_;
        } else {
            target_image = image;
        }

        // write the termination sequence to the image
        *reinterpret_cast<TERM_SEQ*>((target_image) + (im_width_ * im_height_) - sizeof(TERM_SEQ)) = termination_;
        PIXEL* cursor = target_image + offset_;

        while (true) {
            // check for the presence of the termination sequence
            if (*reinterpret_cast<TERM_SEQ*>((cursor) + offset_ - sizeof(TERM_SEQ)) == termination_) {
                return reinterpret_cast<size_t>(cursor) - reinterpret_cast<size_t>(target_image) - offset_;
            }

            // load new pixel value from pre-incremented address
            PIXEL pix_val = *(++cursor);
            if (pix_val <= threshold_center_) continue;

            // first boundary pixel test - decide between MARKER_TEST and SUN_TEST
            if ((pix_val - *(cursor + coord2to1(boundary[0], im_width_))) <= threshold_diff_) {
                if (pix_val >= threshold_sun_) { // SUN_TEST:
                    // testing of the boundary pixels
                    if (boundary_unroll([&](const Point2D& point) -> bool {
                        return (pix_val - *(cursor + coord2to1(point, im_width_))) > threshold_diff_;
                    })) continue;

                    // clearing interior pixels
                    interior_unroll([&](const Point2D& point) -> void {
                        *(cursor + coord2to1(point, im_width_)) = 0x00;
                    });

                    // store the center as a sun point
                    sun_points.push_back(coord1to2(reinterpret_cast<size_t>(cursor) - reinterpret_cast<size_t>(target_image), im_width_));

                    // check the current number of the detected sun points
                    if (sun_points.size() == max_sun_points_count_) {
                        *reinterpret_cast<TERM_SEQ*>((cursor) + offset_ - sizeof(TERM_SEQ)) = termination_;
                    }
                }
            } else { // MARKER_TEST:
                // testing of the boundary pixels
                if (boundary_unroll([&](const Point2D& point) -> bool {
                    return (pix_val - *(cursor + coord2to1(point, im_width_))) <= threshold_diff_;
                })) continue;

                // search for a peak in the interior
                PIXEL peak = 0;
                size_t peak_pos = 0;
                PIXEL* curr_int_ptr = nullptr;

                interior_unroll([&](const Point2D& point) -> void {
                    curr_int_ptr = cursor + coord2to1(point, im_width_);
                    if (*curr_int_ptr > peak) {
                        peak = *curr_int_ptr;
                        peak_pos = reinterpret_cast<size_t>(curr_int_ptr) - reinterpret_cast<size_t>(target_image);
                    }
                    *curr_int_ptr = 0x00;
                });

                // store the detected peak as a marker
                markers.push_back(coord1to2(peak_pos, im_width_));

                // check the current number of the detected markers
                if (markers.size() == max_markers_count_) {
                    *reinterpret_cast<TERM_SEQ*>((cursor) + offset_ - sizeof(TERM_SEQ)) = termination_;
                }
            }

            // go to the next pixel
        }
    };

    /**
     * \brief Gets the image width.
     * \return The width of the image.
     */
    unsigned get_im_width() const { return im_width_; }

    /**
     * \brief Gets the image height.
     * \return The height of the image.
     */
    unsigned get_im_height() const { return im_height_; }

    /**
     * \brief Sets the image size.
     * \param im_width The new width of the image.
     * \param im_height The new height of the image.
     */
    void set_im_size(unsigned im_width, unsigned im_height) {
        im_width_ = im_width;
        im_height_ = im_height;
        offset_ = (im_width_ * RADIUS) + RADIUS;
        if (frame_ != nullptr) {
            frame_ = static_cast<PIXEL*>(std::realloc(frame_, im_width_ * im_height_ * sizeof(PIXEL)));
        }
    }

    /**
     * \brief Gets the center threshold value.
     * \return The center threshold value.
     */
    PIXEL get_threshold_center() const { return threshold_center_; }

    /**
     * \brief Sets the center threshold value.
     * \param threshold_center The new center threshold value.
     */
    void set_threshold_center(PIXEL threshold_center) { threshold_center_ = threshold_center; }

    /**
     * \brief Gets the difference threshold value.
     * \return The difference threshold value.
     */
    PIXEL get_threshold_diff() const { return threshold_diff_; }

    /**
     * \brief Sets the difference threshold value.
     * \param threshold_diff The new difference threshold value.
     */
    void set_threshold_diff(PIXEL threshold_diff) { threshold_diff_ = threshold_diff; }

    /**
     * \brief Gets the sun threshold value.
     * \return The sun threshold value.
     */
    PIXEL get_threshold_sun() const { return threshold_sun_; }

    /**
     * \brief Sets the sun threshold value.
     * \param threshold_sun The new sun threshold value.
     */
    void set_threshold_sun(PIXEL threshold_sun) { threshold_sun_ = threshold_sun; }

    /**
     * \brief Gets the termination sequence.
     * \return The termination sequence.
     */
    TERM_SEQ get_termination() const { return termination_; }

    /**
     * \brief Sets the termination sequence.
     * \param termination The new termination sequence.
     */
    void set_termination(TERM_SEQ termination) { termination_ = termination; }

    /**
     * \brief Gets the maximum number of markers to detect.
     * \return The maximum number of markers.
     */
    unsigned get_max_markers_count() const { return max_markers_count_; }

    /**
     * \brief Sets the maximum number of markers to detect.
     * \param max_markers_count The new maximum number of markers. 0 for no limit.
     */
    void set_max_markers_count(unsigned max_markers_count) { max_markers_count_ = (max_markers_count == 0) ? std::numeric_limits<unsigned>::max() : max_markers_count; }

    /**
     * \brief Gets the maximum number of sun points to detect.
     * \return The maximum number of sun points.
     */
    unsigned get_max_sun_points_count() const { return max_sun_points_count_; }

    /**
     * \brief Sets the maximum number of sun points to detect.
     * \param max_sun_points_count The new maximum number of sun points. 0 for no limit.
     */
    void set_max_sun_points_count(unsigned max_sun_points_count) { max_sun_points_count_ = (max_sun_points_count == 0) ? std::numeric_limits<unsigned>::max() : max_sun_points_count; }

    /**
     * \brief Prints the boundary points of the Bresenham circle.
     * \param os The output stream to print to (stdout by default).
     */
    void print_boundary(std::ostream &os=std::cout) {
        boundary.print_verbose(os);
        boundary.print_sorted(os, boundary);
    }

    /**
     * \brief Prints the interior points of the Bresenham circle.
     * \param os The output stream to print to (stdout by default).
     */
    void print_interior(std::ostream &os=std::cout) {
        interior.print_verbose(os);
    }

private:
    unsigned im_width_;
    unsigned im_height_;
    unsigned offset_;
    PIXEL threshold_center_;
    PIXEL threshold_diff_;
    PIXEL threshold_sun_;
    TERM_SEQ termination_;
    unsigned max_markers_count_;
    unsigned max_sun_points_count_;
    PIXEL* frame_ = nullptr;

    static constexpr auto boundary = BresenhamBoundary<RADIUS>();
    static constexpr auto interior = BresenhamInterior<RADIUS>();

    template<unsigned... I, class F>
    static constexpr bool boundary_unroll_helper(std::integer_sequence<unsigned, I...>, F&& f) {
        return (f(boundary[I]) or ...);
    }

    template<class F>
    static constexpr bool boundary_unroll(F&& f) {
        return boundary_unroll_helper(std::make_integer_sequence<unsigned, boundary.length()>{}, std::forward<F>(f));
    }

    template<unsigned... I, class F>
    static constexpr void interior_unroll_helper(std::integer_sequence<unsigned, I...>, F&& f) {
        (f(interior[I]), ...);
    }

    template<class F>
    static constexpr void interior_unroll(F&& f) {
        interior_unroll_helper(std::make_integer_sequence<unsigned, interior.length()>{}, std::forward<F>(f));
    }
};

};

#endif // FIMD_CPU_HPP
