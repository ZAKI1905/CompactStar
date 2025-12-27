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
 * @file DataColumn.hpp
 * @brief DataColumn: a labeled numeric vector with Mathematica-style indexing and algebraic operators.
 *
 * This header defines:
 *  - ::Zaki::Vector::DataColumn: labeled 1D numeric column (owns contiguous storage).
 *  - Free-function operators and common transforms (exp/log/log10) for DataColumn.
 *
 * DataColumn is designed to be:
 *  - **Data-oriented** (contiguous `double` storage),
 *  - **Label-aware** (for use in DataSet tables),
 *  - **Math-friendly** (vectorized operators and common elementwise transforms),
 *  - **Mathematica-like** (supports negative indexing via `operator[]`).
 *
 * @ingroup Vector
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@eku.edu
 */

#ifndef Zaki_Vector_DataColumn_H
#define Zaki_Vector_DataColumn_H

#include <cstddef>	   // std::size_t
#include <string>	   // std::string
#include <string_view> // std::string_view
#include <vector>	   // std::vector

//--------------------------------------------------------------
namespace Zaki::File
{
/** @brief File mode forward declaration (defined elsewhere). */
enum class FileMode;
} // namespace Zaki::File

//--------------------------------------------------------------
namespace Zaki::Vector
{

// Forward declarations
struct DataColumn;
class DataSet;

//==============================================================
// Free-function algebra on DataColumn
//==============================================================

/**
 * @name Addition (+)
 * Elementwise addition and scalar/vector broadcasting.
 * @{
 */

/// @brief Elementwise addition: `a + b`
DataColumn operator+(const DataColumn &, const DataColumn &);

/// @brief Scalar + column: `c + a`
DataColumn operator+(double, const DataColumn &);

/// @brief Column + scalar: `a + c`
DataColumn operator+(const DataColumn &, double);

/// @brief Column + std::vector: `a + v` (broadcast/elementwise, library-defined semantics)
DataColumn operator+(const DataColumn &, const std::vector<double> &);

/// @brief std::vector + column: `v + a` (broadcast/elementwise, library-defined semantics)
DataColumn operator+(const std::vector<double> &, const DataColumn &);

/** @} */

//..................................................
/**
 * @name Subtraction (-)
 * Elementwise subtraction and scalar/vector broadcasting.
 * @{
 */

/// @brief Elementwise subtraction: `a - b`
DataColumn operator-(const DataColumn &, const DataColumn &);

/// @brief Scalar - column: `c - a`
DataColumn operator-(double, const DataColumn &);

/// @brief Column - scalar: `a - c`
DataColumn operator-(const DataColumn &, double);

/// @brief Column - std::vector: `a - v` (broadcast/elementwise, library-defined semantics)
DataColumn operator-(const DataColumn &, const std::vector<double> &);

/// @brief std::vector - column: `v - a` (broadcast/elementwise, library-defined semantics)
DataColumn operator-(const std::vector<double> &, const DataColumn &);

/** @} */

//..................................................
/**
 * @name Multiplication (*)
 * Elementwise multiplication and scalar scaling.
 * @{
 */

/// @brief Column * scalar: `a * c`
DataColumn operator*(const DataColumn &, double);

/// @brief Scalar * column: `c * a`
DataColumn operator*(double, const DataColumn &);

/// @brief Column * int: `a * n`
DataColumn operator*(const DataColumn &, int);

/// @brief Int * column: `n * a`
DataColumn operator*(int, const DataColumn &);

/// @brief Elementwise multiplication: `a * b`
DataColumn operator*(const DataColumn &, const DataColumn &);

/** @} */

//..................................................
/**
 * @name Division (/)
 * Elementwise division and scalar scaling.
 * @{
 */

/// @brief Column / scalar: `a / c`
DataColumn operator/(const DataColumn &, double);

/// @brief Scalar / column: `c / a`
DataColumn operator/(double, const DataColumn &);

/// @brief Column / int: `a / n`
DataColumn operator/(const DataColumn &, int);

/// @brief Int / column: `n / a`
DataColumn operator/(int, const DataColumn &);

/// @brief Elementwise division: `a / b`
DataColumn operator/(const DataColumn &, const DataColumn &);

/** @} */

//..................................................
/**
 * @name Common transforms
 * Elementwise transforms on DataColumn.
 * @{
 */

/// @brief Elementwise exponential: `exp(x)`
DataColumn exp(const DataColumn &);

/// @brief Elementwise natural logarithm: `log(x)`
DataColumn log(const DataColumn &);

/// @brief Elementwise base-10 logarithm: `log10(x)`
DataColumn log10(const DataColumn &);

/** @} */

//==============================================================
//                        DataColumn Struct
//==============================================================

/**
 * @struct DataColumn
 * @brief Labeled contiguous numeric column with Mathematica-style indexing and basic vector operations.
 *
 * DataColumn owns a `std::vector<double>` and a label. It supports:
 *  - Elementwise arithmetic via free-function operators and assignment operators.
 *  - Transform utilities (moving average smoothing, bump removal).
 *  - Subsetting, trimming, and lookup helpers.
 *  - Negative indexing through `operator[]` for Mathematica-like access.
 *
 * ### Indexing rules (Mathematica-style)
 * - `i >= 0` accesses element `i` from the front (0-based).
 * - `i < 0` accesses from the back (`-1` is last).
 *
 * @note Range checking behavior for `operator[]` is defined in the implementation.
 *       (Your existing library often logs and falls back; keep that policy consistent.)
 */
struct DataColumn
{
  private:
	/// Human-readable label for this column (used by DataSet and exporters).
	std::string label;

	/// Contiguous storage of numeric values.
	std::vector<double> vals;

  public:
	// ----------------------------
	// Friend declarations
	// ----------------------------

	/// @name Friends: arithmetic operators and transforms
	/// @{
	friend DataColumn operator+(const DataColumn &, const DataColumn &);
	friend DataColumn operator+(double, const DataColumn &);
	friend DataColumn operator+(const DataColumn &, double);
	friend DataColumn operator+(const DataColumn &, const std::vector<double> &);
	friend DataColumn operator+(const std::vector<double> &, const DataColumn &);

	friend DataColumn operator-(const DataColumn &, const DataColumn &);
	friend DataColumn operator-(double, const DataColumn &);
	friend DataColumn operator-(const DataColumn &, double);
	friend DataColumn operator-(const DataColumn &, const std::vector<double> &);
	friend DataColumn operator-(const std::vector<double> &, const DataColumn &);

	friend DataColumn operator*(const DataColumn &, double);
	friend DataColumn operator*(double, const DataColumn &);
	friend DataColumn operator*(const DataColumn &, int);
	friend DataColumn operator*(int, const DataColumn &);
	friend DataColumn operator*(const DataColumn &, const DataColumn &);

	friend DataColumn operator/(const DataColumn &, double);
	friend DataColumn operator/(double, const DataColumn &);
	friend DataColumn operator/(const DataColumn &, int);
	friend DataColumn operator/(int, const DataColumn &);
	friend DataColumn operator/(const DataColumn &, const DataColumn &);

	friend DataColumn exp(const DataColumn &);
	friend DataColumn log(const DataColumn &);
	friend DataColumn log10(const DataColumn &);
	/// @}

	// ----------------------------
	// Constructors
	// ----------------------------

	/// @brief Default constructor. Creates an empty column with an empty label.
	DataColumn();

	/**
	 * @brief Construct a column with a known number of rows.
	 * @param rows Number of values to allocate (values are default-initialized).
	 */
	explicit DataColumn(const size_t &rows);

	/**
	 * @brief Construct a column with label and size; optionally fill with a constant.
	 * @param label Column label.
	 * @param size  Number of rows.
	 * @param val   Fill value (default = 0).
	 */
	DataColumn(const std::string &label,
			   const size_t &size,
			   double val = 0);

	/**
	 * @brief Construct a column from (label, values).
	 * @param label Column label.
	 * @param values Value storage to copy.
	 */
	DataColumn(const std::string &label, const std::vector<double> &values);

	/**
	 * @brief Replace the column contents with the given values.
	 *
	 * Capacity is reused when possible.
	 */
	void SetValues(std::vector<double> values)
	{
		vals = std::move(values);
	}

	// ----------------------------
	// Fill helpers
	// ----------------------------

	/**
	 * @brief Fill all values in the column with a constant.
	 * @param in_val Fill value.
	 */
	void Fill(double in_val);

	/**
	 * @brief Conditionally fill values based on a predicate.
	 * @param in_val Fill value.
	 * @param fill_condition Predicate applied to each existing element; if true, the element is replaced with @p in_val.
	 */
	void Fill(double in_val, bool (*fill_condition)(double));

	// ----------------------------
	// Raw data access
	// ----------------------------

	/**
	 * @brief Read-only access to the underlying values vector.
	 * @return Const reference to internal storage.
	 */
	const std::vector<double> &Values() const noexcept { return vals; }

	/**
	 * @brief Pointer to the contiguous storage (read-only).
	 * @return Pointer to first element (or nullptr if empty).
	 */
	const double *Data() const noexcept { return vals.data(); }

	/**
	 * @brief Mutable access to the underlying values vector.
	 * @return Reference to internal storage.
	 */
	std::vector<double> &ValuesMutable() noexcept { return vals; }

	/**
	 * @brief Pointer to the contiguous storage (mutable).
	 * @return Pointer to first element (or nullptr if empty).
	 */
	double *DataMutable() noexcept { return vals.data(); }

	// ----------------------------
	// Label access
	// ----------------------------

	/**
	 * @brief Set the column label.
	 * @param in_label New label.
	 */
	void SetLabel(const std::string &in_label) { label = in_label; }

	/**
	 * @brief Get the column label as a reference.
	 * @return Const reference to label.
	 */
	const std::string &Label() const noexcept { return label; }

	/**
	 * @brief Get the column label as a string_view.
	 * @return Non-owning view of label.
	 */
	std::string_view LabelView() const noexcept { return label; }

	// ----------------------------
	// Iteration
	// ----------------------------

	/// @brief Const iterator to beginning.
	auto begin() const noexcept { return vals.begin(); }

	/// @brief Const iterator to end.
	auto end() const noexcept { return vals.end(); }

	/// @brief Mutable iterator to beginning.
	auto begin() noexcept { return vals.begin(); }

	/// @brief Mutable iterator to end.
	auto end() noexcept { return vals.end(); }

	// ----------------------------
	// Append helpers
	// ----------------------------

	/**
	 * @brief Append a single value.
	 * @param v Value to append.
	 */
	void PushBack(double v) { vals.emplace_back(v); }

	/**
	 * @brief Append @p n values from a raw pointer.
	 * @param ptr Pointer to first value.
	 * @param n   Number of values to append.
	 *
	 * @warning @p ptr must point to at least @p n valid doubles.
	 */
	void Append(const double *ptr, std::size_t n)
	{
		vals.insert(vals.end(), ptr, ptr + n);
	}

	// ----------------------------
	// Smoothing / cleanup
	// ----------------------------

	/**
	 * @brief Return a smoothed copy via moving average.
	 * @param window Moving-average window size.
	 * @return Smoothed DataColumn copy.
	 */
	DataColumn GetSmooth(const short int &window) const;

	/**
	 * @brief Smooth in-place via moving average.
	 * @param window Moving-average window size.
	 */
	void MakeSmooth(const short int &window);

	/**
	 * @brief Remove bumps in-place (implementation-defined heuristic).
	 */
	void RemoveBumps();

	// ----------------------------
	// Arithmetic assignment
	// ----------------------------

	/// @brief Elementwise `+=` with another column (implementation-defined size policy).
	DataColumn &operator+=(const DataColumn &in_dc);

	/// @brief Scalar `+=`.
	DataColumn &operator+=(double in_num);

	/// @brief Elementwise `-=` with another column (implementation-defined size policy).
	DataColumn &operator-=(const DataColumn &in_dc);

	/// @brief Scalar `-=`.
	DataColumn &operator-=(double in_num);

	/// @brief Elementwise `*=` with another column (implementation-defined size policy).
	DataColumn &operator*=(const DataColumn &in_dc);

	/// @brief Scalar `*=`.
	DataColumn &operator*=(double in_num);

	/// @brief Elementwise `/=` with another column (implementation-defined size policy).
	DataColumn &operator/=(const DataColumn &in_dc);

	/// @brief Scalar `/=`.
	DataColumn &operator/=(double in_num);

	// ----------------------------
	// Unary / elementwise utilities
	// ----------------------------

	/**
	 * @brief Unary negation (returns a copy).
	 * @return `(-1) * this`.
	 */
	DataColumn operator-() const;

	/**
	 * @brief Elementwise power (returns a copy).
	 * @param in_num Power exponent.
	 */
	DataColumn pow(double in_num) const;

	/// @brief Elementwise square root (returns a copy).
	DataColumn sqrt() const;

	/// @brief Elementwise absolute value (returns a copy).
	DataColumn Abs() const;

	// ----------------------------
	// Indexing (Mathematica-style)
	// ----------------------------

	/**
	 * @brief Resolve an element index with support for negative ("Mathematica-style") indexing.
	 *
	 * Indexing rules:
	 *  - `i >= 0` → zero-based index from the front.
	 *  - `i < 0`  → index from the back (`-1` = last element).
	 *
	 * @param i   Input index (may be negative).
	 * @param out Resolved non-negative index if successful.
	 * @return true if valid; false otherwise.
	 *
	 * @note This function performs no logging and never throws.
	 */
	bool ResolveIndex(int i, std::size_t &out) const noexcept
	{
		const std::size_t n = vals.size();
		if (n == 0)
			return false;

		if (i >= 0)
		{
			const std::size_t ui = static_cast<std::size_t>(i);
			if (ui < n)
			{
				out = ui;
				return true;
			}
			return false;
		}

		const std::size_t k = static_cast<std::size_t>(-i); // -1 => 1
		if (k == 0 || k > n)
			return false;

		out = n - k;
		return true;
	}

	/**
	 * @brief Element access (supports negative indexing).
	 *
	 * @warning If index is out of range, logs an error and returns vals[0] as fallback.
	 */
	double &At(int i);

	/** @copydoc At(int) */
	double At(int i) const;

	/**
	 * @brief Element access with Mathematica-style indexing (mutable).
	 * @param i Index (supports negative indexing).
	 * @return Reference to element.
	 */
	double &operator[](int i) { return At(i); }

	/**
	 * @brief Element access with Mathematica-style indexing (const).
	 * @param i Index (supports negative indexing).
	 * @return Value copy.
	 */
	double operator[](int i) const { return At(i); }
	// ----------------------------

	// ----------------------------
	// Combine into DataSet
	// ----------------------------

	/**
	 * @brief Combine this column with another column into a DataSet (two columns).
	 * @param other Column to append as second column.
	 * @return DataSet with two columns.
	 */
	DataSet CombineColumns(const DataColumn &other) const;

	/**
	 * @brief Combine this column with a list of columns into a DataSet.
	 * @param others Additional columns to append.
	 * @return DataSet with (1 + others.size()) columns.
	 */
	DataSet CombineColumns(const std::vector<DataColumn> &others) const;

	// ----------------------------
	// Subset utilities
	// ----------------------------

	/**
	 * @brief Return subset of elements satisfying a scalar predicate.
	 * @param cond Predicate on value; elements where cond(value) is true are kept.
	 * @return Subset DataColumn.
	 */
	DataColumn GetSubSet(bool (*cond)(double)) const;

	/**
	 * @brief Return subset of elements satisfying a predicate that depends on (column, index).
	 * @param cond Predicate called as cond(*this, idx).
	 * @return Subset DataColumn.
	 */
	DataColumn GetSubSet(bool (*cond)(const DataColumn &, int idx)) const;

	/**
	 * @brief Return subset in index range [idx_1, idx_2] (implementation-defined inclusivity).
	 * @param idx_1 Start index.
	 * @param idx_2 End index.
	 * @return Subset DataColumn.
	 */
	DataColumn GetSubSet(const size_t idx_1, const size_t idx_2) const;

	/**
	 * @brief Trim the column in-place to the index range [idx_1, idx_2] (implementation-defined inclusivity).
	 * @param idx_1 Start index.
	 * @param idx_2 End index.
	 * @warning Modifies the current DataColumn.
	 * @note Exclusive of idx_2
	 */
	void Trim(const size_t idx_1, const size_t idx_2);

	// ----------------------------
	// Search helpers
	// ----------------------------

	/**
	 * @brief Get the first index meeting a monotonic threshold condition.
	 *
	 * If the column is strictly increasing:
	 *  - returns first index i with vals[i] >= x.
	 *
	 * If strictly decreasing:
	 *  - returns first index i with vals[i] <= x.
	 *
	 * @param x Target value.
	 * @return Index satisfying the condition, or implementation-defined sentinel if not found.
	 */
	int GetFirstIdx(double x) const;

	/**
	 * @brief Get index of element closest to a target value.
	 * @param x Target value.
	 * @return Index of closest element.
	 */
	int GetClosestIdx(double x) const;

	// ----------------------------
	// Capacity / size
	// ----------------------------

	/// @brief Number of stored values.
	size_t Size() const;

	/// @brief Reserve capacity for @p n values.
	void Reserve(const size_t &n);

	/// @brief Resize storage to @p n values.
	void Resize(const size_t &n);

	/**
	 * @brief Clear logical contents but keep capacity for reuse.
	 *
	 * This is the fast reset operation when you expect to refill the column.
	 */
	void Clear() noexcept { vals.clear(); }

	/**
	 * @brief Clear contents and release allocated capacity.
	 *
	 * This is the expensive reset operation used when you want to reclaim memory.
	 */
	void Release() noexcept
	{
		vals.clear();
		vals.shrink_to_fit();
	}

	/// @brief True if no values are stored.
	bool Empty() const noexcept { return vals.empty(); }

	// ----------------------------
	// Min/max helpers
	// ----------------------------

	/// @brief Minimum value (requires non-empty column; behavior defined in implementation).
	double Min() const;

	/// @brief Maximum value (requires non-empty column; behavior defined in implementation).
	double Max() const;

	/// @brief Index of minimum value.
	int MinIdx() const;

	/// @brief Index of maximum value.
	int MaxIdx() const;

	// ----------------------------
	// Convenience: Add with fill policy
	// ----------------------------

	/**
	 * @brief Add two columns with a fill policy for missing elements.
	 *
	 * Use when the sizes differ. Missing elements are replaced with @p fill.
	 *
	 * @param dc_1 First column.
	 * @param dc_2 Second column.
	 * @param fill Fill value for out-of-range elements.
	 * @return Resulting column.
	 */
	static DataColumn Add(const DataColumn &dc_1,
						  const DataColumn &dc_2,
						  double fill = 0);
};

//--------------------------------------------------------------
} // namespace Zaki::Vector
//--------------------------------------------------------------

#endif /* Zaki_Vector_DataColumn_H */