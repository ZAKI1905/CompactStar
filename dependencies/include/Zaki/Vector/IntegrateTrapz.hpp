// -*- lsst-c++ -*-
/*
 * Zaki's Common Library
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025 Mohammadreza Zakeri
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */
/**
 * @file IntegrateTrapz.hpp
 *
 * @brief Fast trapezoidal-rule integration utilities for Vector objects.
 *
 * This header provides:
 *  - A fast trapezoidal integrator over index ranges (i0, i1)
 *  - Convenience wrappers over physical bounds (xmin, xmax) for monotonic x
 *  - A helper to resolve (xmin, xmax) -> (i0, i1) once, then reuse indices
 *
 * The design is intentionally "Vector"-level (not Math-level) because these
 * utilities integrate discrete sampled data held in DataColumn / std::vector.
 *
 * @ingroup Vector
 */

#ifndef Zaki_Vector_IntegrateTrapz_H
#define Zaki_Vector_IntegrateTrapz_H

#include <cstddef>
#include <utility>
#include <vector>

#include "Zaki/Vector/DataColumn.hpp"

namespace Zaki::Vector
{

/**
 * @brief Resolve a physical x-range (xmin, xmax) into an index range [i0, i1]
 *        suitable for trapezoidal integration on sampled data.
 *
 * This helper assumes that `x` is **monotonic** (strictly increasing or strictly
 * decreasing). It uses binary search to find the indices of samples nearest to
 * xmin and xmax, then returns a **sorted** index pair `(i0, i1)` with `i0 <= i1`.
 *
 * Typical workflow for maximum performance:
 * @code
 * auto [i0, i1] = FindIndexRangeMonotonic(x, xmin, xmax);
 * double I1 = IntegrateTrapzIdx(x, y1, i0, i1);
 * double I2 = IntegrateTrapzIdx(x, y2, i0, i1);
 * @endcode
 *
 * @param x     Monotonic sample locations.
 * @param xmin  Lower physical bound (not necessarily ordered relative to xmax).
 * @param xmax  Upper physical bound.
 *
 * @return Pair `(i0, i1)` with `i0 <= i1`.
 *
 * @warning If `x` is not monotonic, the returned indices may be meaningless.
 *
 * @note The returned indices typically correspond to the nearest sample points
 *       to the requested bounds; there is no partial-cell correction here.
 */
std::pair<std::size_t, std::size_t>
FindIndexRangeMonotonic(const DataColumn &x, double xmin, double xmax);

/**
 * @brief Vector overload of FindIndexRangeMonotonic.
 *
 * @param x     Monotonic sample locations.
 * @param xmin  Lower physical bound (not necessarily ordered relative to xmax).
 * @param xmax  Upper physical bound.
 * @return Pair `(i0, i1)` with `i0 <= i1`.
 */
std::pair<std::size_t, std::size_t>
FindIndexRangeMonotonic(const std::vector<double> &x, double xmin, double xmax);

/**
 * @brief Fast trapezoidal integration over an index range [i0, i1] on DataColumns.
 *
 * Computes:
 * \f[
 * \int y(x)\,dx \approx \sum_{k=i0}^{i1-1} \frac{(x_{k+1}-x_k)}{2}\,(y_k+y_{k+1})
 * \f]
 *
 * This is the **fast core**: no searching, no allocation, no interpolation.
 *
 * @param x   Sample locations (same length as y).
 * @param y   Sample values.
 * @param i0  Start index (inclusive).
 * @param i1  End index (inclusive). Must satisfy `i1 >= i0`.
 *
 * @return Approximate integral using the trapezoidal rule.
 *
 * @warning Behavior is undefined if sizes mismatch or indices are out of range.
 *         In debug builds you may want to add asserts; in release this stays fast.
 */
double IntegrateTrapzIdx(const DataColumn &x, const DataColumn &y,
						 std::size_t i0, std::size_t i1) noexcept;

/**
 * @brief Fast trapezoidal integration over an index range [i0, i1] on vectors.
 *
 * @param x   Sample locations (same length as y).
 * @param y   Sample values.
 * @param i0  Start index (inclusive).
 * @param i1  End index (inclusive).
 * @return Approximate integral using the trapezoidal rule.
 */
double IntegrateTrapzIdx(const std::vector<double> &x, const std::vector<double> &y,
						 std::size_t i0, std::size_t i1) noexcept;

/**
 * @brief Fast trapezoidal integration over an index range [i0, i1] for DataColumn & vector.
 *
 * @param x   Sample locations (same length as y).
 * @param y   Sample values.
 * @param i0  Start index (inclusive).
 * @param i1  End index (inclusive).
 * @return Approximate integral using the trapezoidal rule.
 */
double IntegrateTrapzIdx(const DataColumn &x, const std::vector<double> &y,
						 std::size_t i0, std::size_t i1) noexcept;

/**
 * @brief Convenience trapezoidal integration over physical bounds (xmin, xmax) on DataColumns.
 *
 * Internally resolves indices using FindIndexRangeMonotonic(x, xmin, xmax),
 * then calls IntegrateTrapzIdx(...).
 *
 * @param x     Monotonic sample locations.
 * @param y     Sample values (same length as x).
 * @param xmin  Physical bound.
 * @param xmax  Physical bound.
 * @return Approximate integral on the sampled grid.
 */
double IntegrateTrapz(const DataColumn &x, const DataColumn &y,
					  double xmin, double xmax);

/**
 * @brief Convenience trapezoidal integration over physical bounds (xmin, xmax) on vectors.
 *
 * @param x     Monotonic sample locations.
 * @param y     Sample values (same length as x).
 * @param xmin  Physical bound.
 * @param xmax  Physical bound.
 * @return Approximate integral on the sampled grid.
 */
double IntegrateTrapz(const std::vector<double> &x, const std::vector<double> &y,
					  double xmin, double xmax);

/**
 * @brief Mixed overload: integrate DataColumn y over std::vector x.
 *
 * This is useful if your x-grid is held externally as a plain vector but
 * you want to reuse DataColumn for y values/labels elsewhere.
 *
 * @param x     Monotonic sample locations.
 * @param y     Sample values.
 * @param xmin  Physical bound.
 * @param xmax  Physical bound.
 */
double IntegrateTrapz(const std::vector<double> &x, const DataColumn &y,
					  double xmin, double xmax);

/**
 * @brief Mixed overload: integrate std::vector y over DataColumn x.
 *
 * @param x     Monotonic sample locations.
 * @param y     Sample values.
 * @param xmin  Physical bound.
 * @param xmax  Physical bound.
 */
double IntegrateTrapz(const DataColumn &x, const std::vector<double> &y,
					  double xmin, double xmax);

} // namespace Zaki::Vector

#endif // Zaki_Vector_IntegrateTrapz_H