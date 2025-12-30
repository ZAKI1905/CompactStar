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
 * @file DataSet.hpp
 *
 * @brief DataColumns and DataSets generalize std::vector to Mathematica style sets.
 *
 * @ingroup Vector
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
 */

#ifndef Zaki_Vector_DataSet_H
#define Zaki_Vector_DataSet_H

// #include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include <vector>

#include <sstream>
#include <string>

// For SplineSet
#include <map>

#include <Zaki/File/VecSaver.hpp>
#include <Zaki/Math/Math_Core.hpp>
#include <Zaki/String/Directory.hpp>

#include "Zaki/Vector/DataColumn.hpp"

//--------------------------------------------------------------
namespace Zaki::Math
{

// Math forward declarations
// template <typename T>
// struct Range ;

// template<typename FuncObj, typename MemFuncPtr >
// class GSLFuncWrapper ;

} // namespace Zaki::Math

//--------------------------------------------------------------
namespace Zaki::File
{

// File forward declarations
enum class FileMode;

} // namespace Zaki::File

//--------------------------------------------------------------
namespace Zaki::Vector
{

// Vector forward declarations
// struct DataColumn;
class DataSet;

// //==============================================================
// //..................................................
// //                Addition (+)
// //..................................................
// /// Addition of a DataColumn to a DataColumn
// DataColumn operator+(const DataColumn &, const DataColumn &);

// /// Addition of a DataColumn to a double
// DataColumn operator+(const double &, const DataColumn &);

// /// Addition of a double to a DataColumn
// DataColumn operator+(const DataColumn &, const double &);

// /// Addition of a list of doubles to a DataColumn
// DataColumn operator+(const DataColumn &, const std::vector<double> &);

// /// Addition of a list of doubles to a DataColumn
// DataColumn operator+(const std::vector<double> &, const DataColumn &);

// //..................................................
// //                Subtraction (-)
// //..................................................
// /// Subtraction of a DataColumn from a DataColumn
// DataColumn operator-(const DataColumn &, const DataColumn &);

// /// Subtraction of a DataColumn from a double
// DataColumn operator-(const double &, const DataColumn &);

// /// Subtraction of a double from a DataColumn
// DataColumn operator-(const DataColumn &, const double &);

// /// Subtraction of a list of doubles from a DataColumn
// DataColumn operator-(const DataColumn &, const std::vector<double> &);

// /// Subtraction of a list of doubles from a DataColumn
// DataColumn operator-(const std::vector<double> &, const DataColumn &);

// //..................................................
// //              Multiplication (*)
// //..................................................
// /// Multiplication of a DataColumn by a double
// DataColumn operator*(const DataColumn &, const double &);

// /// Multiplication of a double by a DataColumn
// DataColumn operator*(const double &, const DataColumn &);

// /// Multiplication of a DataColumn by an int
// DataColumn operator*(const DataColumn &, const int &);

// /// Multiplication of an int by a DataColumn
// DataColumn operator*(const int &, const DataColumn &);

// /// Multiplication of a DataColumn by another DataColumn
// DataColumn operator*(const DataColumn &, const DataColumn &);
// //..................................................
// //                  Division (/)
// //..................................................

// /// Division of a DataColumn by a double
// DataColumn operator/(const DataColumn &, const double &);

// /// Division of a double by a DataColumn
// DataColumn operator/(const double &, const DataColumn &);

// /// Division of a DataColumn by an int
// DataColumn operator/(const DataColumn &, const int &);

// /// Division of an int by a DataColumn
// DataColumn operator/(const int &, const DataColumn &);

// /// Division of a DataColumn by another DataColumn
// DataColumn operator/(const DataColumn &, const DataColumn &);

// //..................................................
// //               Exponentiation (exp)
// //..................................................
// DataColumn exp(const DataColumn &);
// DataColumn log(const DataColumn &);
// DataColumn log10(const DataColumn &);

// //==============================================================

// //==============================================================
// //                        DataColumn Struct
// //==============================================================
// struct DataColumn
// {
//   private:
// 	// in case of missing data (NaN) this will be used to fill
// 	// the data:
// 	// double add_def_val = 0 ;

// 	std::string label;
// 	std::vector<double> vals;

//   public:
// 	//..................................................
// 	friend DataColumn operator+(const DataColumn &, const DataColumn &);
// 	friend DataColumn operator+(const double &, const DataColumn &);
// 	friend DataColumn operator+(const DataColumn &, const double &);
// 	friend DataColumn operator+(const DataColumn &, const std::vector<double> &);
// 	friend DataColumn operator+(const std::vector<double> &, const DataColumn &);
// 	//..................................................
// 	friend DataColumn operator-(const DataColumn &, const DataColumn &);
// 	friend DataColumn operator-(const double &, const DataColumn &);
// 	friend DataColumn operator-(const DataColumn &, const double &);
// 	friend DataColumn operator-(const DataColumn &, const std::vector<double> &);
// 	friend DataColumn operator-(const std::vector<double> &, const DataColumn &);
// 	//..................................................
// 	friend DataColumn operator*(const DataColumn &, const double &);
// 	friend DataColumn operator*(const double &, const DataColumn &);
// 	friend DataColumn operator*(const DataColumn &, const int &);
// 	friend DataColumn operator*(const int &, const DataColumn &);
// 	friend DataColumn operator*(const DataColumn &, const DataColumn &);
// 	//..................................................
// 	friend DataColumn operator/(const DataColumn &, const double &);
// 	friend DataColumn operator/(const double &, const DataColumn &);
// 	friend DataColumn operator/(const DataColumn &, const int &);
// 	friend DataColumn operator/(const int &, const DataColumn &);
// 	friend DataColumn operator/(const DataColumn &, const DataColumn &);
// 	//..................................................
// 	friend DataColumn exp(const DataColumn &);
// 	DataColumn log(const DataColumn &);
// 	DataColumn log10(const DataColumn &);
// 	//..................................................
// 	/// Default constructor
// 	DataColumn();

// 	/// Constructor from a known size for rows
// 	DataColumn(const size_t &rows);

// 	/// Constructor from a label and size and values
// 	/// by default it fills the column with zeros
// 	DataColumn(const std::string &label,
// 			   const size_t &size,
// 			   const double &val = 0);

// 	/// Constructor from a vector and label:
// 	DataColumn(const std::string &, const std::vector<double> &);

// 	/// Fills all of the elements of the data column with the input value.
// 	void Fill(const double &in_val);

// 	/// Fills specific elements of the data column with the input
// 	/// value if the input condition function returns true.
// 	void Fill(const double &in_val, bool (*fill_condition)(const double &));

// 	/**
// 	 * @brief Returns values vector
// 	 *
// 	 * @return const std::vector<double>&
// 	 */
// 	const std::vector<double> &Values() const noexcept { return vals; }

// 	/**
// 	 * @brief Returns pointer to the data
// 	 *
// 	 * @return const double*
// 	 */
// 	const double *Data() const noexcept { return vals.data(); }

// 	/**
// 	 * @brief Returns mutable values vector
// 	 *
// 	 * @return std::vector<double>&
// 	 */
// 	std::vector<double> &ValuesMutable() noexcept { return vals; }

// 	/**
// 	 * @brief Returns mutable pointer to the data
// 	 *
// 	 * @return double*
// 	 */
// 	double *DataMutable() noexcept { return vals.data(); }

// 	/**
// 	 * @brief Sets the label of the DataColumn
// 	 *
// 	 * @param in_label The new label for the DataColumn
// 	 */
// 	void SetLabel(const std::string &in_label) { label = in_label; }

// 	/**
// 	 * @brief  Returns the label of the DataColumn
// 	 * @return const std::string&
// 	 */
// 	const std::string &Label() const noexcept { return label; }

// 	/**
// 	 * @brief  Returns the label of the DataColumn
// 	 * @return std::string_view
// 	 */
// 	std::string_view LabelView() const noexcept { return label; }

// 	/**
// 	 * @brief  Returns iterator to beginning of the DataColumn
// 	 */
// 	auto begin() const noexcept { return vals.begin(); }

// 	/**
// 	 * @brief  Returns iterator to end of the DataColumn
// 	 */
// 	auto end() const noexcept { return vals.end(); }

// 	/// @brief  Returns mutable iterator to beginning of the DataColumn
// 	auto begin() noexcept { return vals.begin(); }

// 	/// @brief  Returns mutable iterator to end of the DataColumn
// 	auto end() noexcept { return vals.end(); }

// 	/// @brief  Pushes back a value to the DataColumn
// 	void PushBack(double v) { vals.emplace_back(v); } // non-throwing not guaranteed

// 	/// @brief  Appends n values from ptr to the DataColumn
// 	void Append(const double *ptr, std::size_t n)
// 	{
// 		vals.insert(vals.end(), ptr, ptr + n);
// 	}

// 	/// @brief  Generates a smoother datacolumn by moving average
// 	/// @param window window for moving average
// 	/// @return Smoother datacolumn
// 	DataColumn GetSmooth(const short int &window) const;

// 	/// @brief  Makes the datacolumn smoother by moving average
// 	/// @param window window for moving average
// 	void MakeSmooth(const short int &window);

// 	/// @brief  Makes the datacolumn smoother by removing bumps
// 	void RemoveBumps();
// 	//..................................................
// 	// Operator overloading
// 	//..................................................
// 	// Addition (+)
// 	// DataColumn operator+(const DataColumn& in_dc) const ;

// 	/// Adds the two DataColumns and fills the non-existent
// 	/// elements with '0'
// 	/// Use when the sizes of DC's are not the same
// 	static DataColumn Add(const DataColumn &dc_1,
// 						  const DataColumn &dc_2,
// 						  const double &fill = 0);

// 	/// Addition to a list of numbers (+)
// 	// DataColumn operator+(const std::vector<double>& in_vec) const ;
// 	/// Addition to a single number (+)
// 	// DataColumn operator+(const double& in_num) const ;
// 	//..................................................
// 	/// Addition assignment operator
// 	DataColumn &operator+=(const DataColumn &in_dc);
// 	DataColumn &operator+=(const double &in_num);

// 	/// Subtraction assignment operator
// 	DataColumn &operator-=(const DataColumn &in_dc);
// 	DataColumn &operator-=(const double &in_num);

// 	/// Multiplication assignment operator
// 	DataColumn &operator*=(const DataColumn &in_dc);
// 	DataColumn &operator*=(const double &in_num);

// 	/// Division assignment operator
// 	DataColumn &operator/=(const DataColumn &in_dc);
// 	DataColumn &operator/=(const double &in_num);

// 	//..................................................
// 	// // Subtraction (-)
// 	// DataColumn operator-(const DataColumn& in_dc) const ;
// 	// // Subtraction of a single number (-)
// 	// DataColumn operator-(const double& in_num) const ;
// 	// // Subtraction of a list of numbers (-)
// 	// DataColumn operator-(const std::vector<double>& in_vec) const ;
// 	//..................................................
// 	/// Unary minus (negation '-') operator
// 	DataColumn operator-() const;
// 	//..................................................
// 	/// Power
// 	DataColumn pow(const double &in_num) const;
// 	// DataColumn pow^(const int& in_num) const ;
// 	/// sqrt
// 	DataColumn sqrt() const;
// 	/// absolute value
// 	DataColumn Abs() const;
// 	// // Exponentiation
// 	// DataColumn exp(const DataColumn&) const ;
// 	//..................................................
// 	/// Overloading []
// 	double &operator[](const int);

// 	/// Overloading [] ( const )
// 	double operator[](const int) const;
// 	//..................................................

// 	DataSet CombineColumns(const DataColumn &) const;
// 	DataSet CombineColumns(const std::vector<DataColumn> &) const;

// 	/// Selecting a subset of DataColumn elements satisfying
// 	///  the input condition 'cond'
// 	DataColumn GetSubSet(bool (*cond)(const double &)) const;

// 	/// Selecting a subset of DataColumn elements satisfying
// 	///  the input condition 'cond'
// 	DataColumn GetSubSet(bool (*cond)(const DataColumn &,
// 									  const int &idx)) const;

// 	/// Selecting a subset of DataColumn ranging from
// 	///  idx_1 --> idx_2
// 	DataColumn GetSubSet(const size_t idx_1, const size_t idx_2) const;

// 	/// Trims the data column to its subset ranging from
// 	///  idx_1 --> idx_2
// 	void Trim(const size_t idx_1, const size_t idx_2);

// 	/// Returns the index to :
// 	///  If the data is:
// 	///  Ascending: the first value that is equal or greater than
// 	///  the input value.
// 	///  Descending: the first value that is equal or less than
// 	///  the input value.
// 	///  The data column must be strictly increasing or decreasing.
// 	int GetFirstIdx(const double &) const;

// 	/// Returns the index to the element that is closest.
// 	int GetClosestIdx(const double &) const;

// 	/// Returns the size of 'vals'
// 	size_t Size() const;

// 	/// Reserves space for 'vals'
// 	void Reserve(const size_t &);

// 	/// Resizes the 'vals'
// 	void Resize(const size_t &);

// 	/// Clears the 'vals'
// 	void Clear() noexcept { vals.clear(); }

// 	/// Releases the memory of 'vals'
// 	void Release() noexcept
// 	{
// 		vals.clear();
// 		vals.shrink_to_fit();
// 	}

// 	/// Checks if 'vals' is empty
// 	bool Empty() const noexcept
// 	{
// 		return vals.empty();
// 	}

// 	/// Returns the minimum element
// 	double Min() const;

// 	/// Returns the maximum element
// 	double Max() const;

// 	/// Returns the minimum element's index
// 	int MinIdx() const;

// 	/// Returns the maximum element's index
// 	int MaxIdx() const;
// };

//==============================================================
//                        Row Struct
//==============================================================
struct Row
{
  public:
	std::vector<double> vals;
	char delim = '\t';

	// Precision for printing (double) numbers
	int precision = 8;

	Row() {}

	Row(const std::vector<double> &in_vals)
	{
		vals = in_vals;
	}

	std::string Str() const
	{
		std::stringstream ss;

		// [Apr 12, 2023] : Delim issue was resolved.
		for (size_t i = 0; i < vals.size() - 1; i++)
		{
			char tmp[150];
			snprintf(tmp, sizeof(tmp), "%-*.*e%c ", 7 + precision,
					 precision, vals[i], delim);
			ss << tmp;
		}

		// Last column
		char tmp[150];
		snprintf(tmp, sizeof(tmp), "%-*.*e ", 7 + precision,
				 precision, vals[vals.size() - 1]);
		ss << tmp;

		return ss.str();
	}
};

//==============================================================
//                        SplineIdx Struct
//==============================================================
// struct SplineIdx
// {
// 	int idx = -1;
// 	gsl_spline *spline = nullptr;

// 	// const gsl_interp_type* dataset_gsl_interp_type = gsl_interp_steffen ;
// 	const gsl_interp_type *dataset_gsl_interp_type = gsl_interp_linear;

// 	SplineIdx() {};

// 	SplineIdx(const int &in_idx, const double *x,
// 			  const double *y, const size_t &in_size)
// 		: idx(in_idx)
// 	{
// 		spline = gsl_spline_alloc(dataset_gsl_interp_type, in_size);
// 		gsl_spline_init(spline, x, y, in_size);
// 	}

// 	~SplineIdx()
// 	{
// 		if (spline)
// 		{
// 			gsl_spline_free(spline);
// 			spline = nullptr;
// 		}
// 	}

// 	void Reset(const double *x, const double *y, const size_t &in_size)
// 	{
// 		if (spline)
// 		{
// 			gsl_spline_free(spline);
// 		}

// 		spline = gsl_spline_alloc(dataset_gsl_interp_type, in_size);
// 		gsl_spline_init(spline, x, y, in_size);
// 	}

// 	// Equality operator
// 	bool operator==(const int &in_idx) const
// 	{
// 		return in_idx == idx;
// 	}

// 	/// Copy Constructor
// 	SplineIdx(const SplineIdx &) = delete;

// 	/// Assignment operator
// 	SplineIdx &operator=(const SplineIdx &) = delete;

// 	/// Move constructor
// 	SplineIdx(SplineIdx &&other)
// 		: idx(other.idx),
// 		  spline(other.spline)
// 	{
// 		// mov_con_count++ ;
// 		other.spline = nullptr;
// 	}

// 	/// Move Assignment
// 	SplineIdx &operator=(SplineIdx &&other)
// 	{
// 		// mov_ass_count++ ;
// 		idx = other.idx;
// 		spline = other.spline;
// 		other.spline = nullptr;
// 		return *this;
// 	}
// };

struct SplineIdx
{
	gsl_spline *spline = nullptr;
	std::size_t n = 0; // number of samples used to initialize

	// Interpolation type
	inline static const gsl_interp_type *kInterpType = gsl_interp_linear;
	// or: inline static const gsl_interp_type* kInterpType = gsl_interp_steffen;

	SplineIdx() = default;

	SplineIdx(const double *x, const double *y, std::size_t in_size)
	{
		Reset(x, y, in_size);
	}

	~SplineIdx() { Free(); }

	SplineIdx(const SplineIdx &) = delete;
	SplineIdx &operator=(const SplineIdx &) = delete;

	SplineIdx(SplineIdx &&other) noexcept
		: spline(other.spline), n(other.n)
	{
		other.spline = nullptr;
		other.n = 0;
	}

	SplineIdx &operator=(SplineIdx &&other) noexcept
	{
		if (this == &other)
			return *this;
		Free(); // IMPORTANT: avoid leak
		spline = other.spline;
		n = other.n;
		other.spline = nullptr;
		other.n = 0;
		return *this;
	}

	bool IsReady() const noexcept { return spline != nullptr; }

	void Free() noexcept
	{
		if (spline)
		{
			gsl_spline_free(spline);
			spline = nullptr;
			n = 0;
		}
	}

	void Reset(const double *x, const double *y, std::size_t in_size);
};

//==============================================================
//                        SplineSet Struct
//==============================================================
// struct SplineSet
// {
// 	std::vector<SplineIdx> spline_set;

// 	void AssignSpline(const int &in_idx, const double *x,
// 					  const double *y, const size_t &in_size);

// 	double Evaluate(const int &in_idx, const double &in_x,
// 					gsl_interp_accel *in_accel) const;

// 	double Derivative(const int &in_idx, const double &in_x,
// 					  gsl_interp_accel *in_accel) const;

// 	double Integrate(const int &in_idx, const Zaki::Math::Range<double> &in_range,
// 					 gsl_interp_accel *in_accel) const;

// 	/// Default constructor
// 	SplineSet() {}

// 	/// Copy Constructor
// 	SplineSet(const SplineSet &) = delete;

// 	// Assignment operator
// 	SplineSet &operator=(const SplineSet &) = delete;

// 	void Resize(const size_t &in_size)
// 	{
// 		spline_set.resize(in_size);
// 	}
// };

struct SplineSet
{
	std::vector<SplineIdx> spline_set;

	SplineSet() = default;
	SplineSet(const SplineSet &) = delete;
	SplineSet &operator=(const SplineSet &) = delete;

	void Resize(std::size_t in_size)
	{
		spline_set.resize(in_size);
	}

	std::size_t Size() const noexcept { return spline_set.size(); }

	void AssignSpline(std::size_t idx, const double *x, const double *y, std::size_t n);

	double Evaluate(std::size_t idx, double x, gsl_interp_accel *accel) const;

	double Derivative(std::size_t idx, double x, gsl_interp_accel *accel) const;

	double Integrate(std::size_t idx,
					 const Zaki::Math::Range<double> &r,
					 gsl_interp_accel *accel) const;
};

//==============================================================
//                        DataSet Class
//==============================================================
class DataSet
{
  public:
	struct PlotParam
	{
	  private:
		// ----------------------------------------------
		/// Figure dimensions
		struct FigSize
		{
			/// Width of the figure
			size_t w = 800;

			/// Height of the figure
			size_t h = 600;

			/// Constructor
			FigSize(const size_t &in_w, const size_t &in_h)
				: w(in_w), h(in_h) {}
		};
		// ----------------------------------------------

		// ----------------------------------------------
		/// Axis Ticks dimensions
		struct AxisTicks
		{
			/// Location of ticks
			std::vector<double> ticks;

			/// Label of ticks
			std::vector<std::string> labels;

			/// Labels are set or no?
			bool label_flag = false;

			/// Constructor-0
			AxisTicks() {}

			/// Constructor-1 from tick locations
			AxisTicks(const std::vector<double> &in_ticks)
				: ticks(in_ticks) {}

			/// Constructor-2 from tick locations & labels
			AxisTicks(const std::vector<double> &in_ticks,
					  const std::vector<std::string> &in_labels)
				: ticks(in_ticks), labels(in_labels)
			{
				if (in_ticks.size() == in_labels.size())
				{
					label_flag = true;
				}
				else
				{
					Z_LOG_ERROR("The sizes of ticks and labels should match!");
					Z_LOG_INFO("Axis labels are ignored.");
				}
			}

			void SetTicks(const std::vector<double> &in_ticks)
			{
				ticks = in_ticks;
			}

			void SetTickLabels(const std::vector<std::string> &in_labels)
			{
				labels = in_labels;
				label_flag = true;
			}
		};
		// ----------------------------------------------

		// ----------------------------------------------
		/// @brief The legend for plots
		struct Legend
		{
			/// @brief Alignment inside the bounding box
			/// can be any of these values:
			/// {'best', 'center', 'center left',
			///  'center right', 'lower center',
			/// 'lower left', 'lower right', 'right',
			/// 'upper center', 'upper left', 'upper right'}
			std::string loc = "best";

			/// @brief  The coordinate of the bottom-left corner
			double x0 = 0, y0 = 0;

			/// @brief Width of the BBox
			double w = 0;

			/// @brief Height of the BBox
			double h = 0;

			// std::map<std::string, std::string> &keywords = {} ;

			/// @brief Default Constructor
			Legend() {}

			/// @brief Constructor-1
			Legend(const double &in_x0, const double &in_y0,
				   const double &in_w = 0, const double &in_h = 0)
				: x0(in_x0), y0(in_y0), w(in_w), h(in_h)
			{
			}

			/// @brief Constructor-2
			Legend(const std::string &in_loc,
				   const double &in_x0, const double &in_y0,
				   const double &in_w = 0, const double &in_h = 0)
				: loc(in_loc), x0(in_x0), y0(in_y0), w(in_w), h(in_h)
			{
			}
		};
		// ----------------------------------------------

		Zaki::Math::Range<double> x_ax = {0, 0};
		Zaki::Math::Range<double> y_ax = {0, 0};
		FigSize fig_size = {800, 600};
		bool grid = false;
		AxisTicks x_ticks;
		AxisTicks y_ticks;

		std::string x_ax_label;
		std::string y_ax_label;
		std::vector<
			std::pair<double,
					  std::map<std::string, std::string>>>
			ax_h_lines;

		std::vector<
			std::pair<double,
					  std::map<std::string, std::string>>>
			ax_v_lines;

		Legend legend;
		// ---------------------------
		//         Flags
		// ---------------------------
		bool x_ax_flag = false;
		bool y_ax_flag = false;

		bool fig_size_flag = false;
		bool grid_flag = false;
		bool x_ticks_flag = false;
		bool y_ticks_flag = false;
		bool x_ax_label_flag = false;
		bool y_ax_label_flag = false;
		bool legend_flag = false;
		bool ax_h_line_flag = false;
		bool ax_v_line_flag = false;

	  public:
		PlotParam() {}

		void AddAxHLine(const double &y,
						const std::map<std::string, std::string> &options = {})
		{
			ax_h_lines.emplace_back(y, options);
			ax_h_line_flag = true;
		}

		void AddAxVLine(const double &x,
						const std::map<std::string, std::string> &options = {})
		{
			ax_v_lines.emplace_back(x, options);
			ax_v_line_flag = true;
		}

		void SetXAxis(const Zaki::Math::Range<double> &in_x_ax)
		{
			x_ax = in_x_ax;
			x_ax_flag = true;
		}

		void SetYAxis(const Zaki::Math::Range<double> &in_y_ax)
		{
			y_ax = in_y_ax;
			y_ax_flag = true;
		}

		void SetFigSize(const FigSize &in_fig_Size)
		{
			fig_size = in_fig_Size;
			fig_size_flag = true;
		}

		void SetGrid(const bool &in_grid = true)
		{
			grid = in_grid;
			grid_flag = true;
		}

		void SetXTicks(const AxisTicks &in_x_ticks)
		{
			x_ticks = in_x_ticks;
			x_ticks_flag = true;
		}

		void SetYTicks(const AxisTicks &in_y_ticks)
		{
			y_ticks = in_y_ticks;
			y_ticks_flag = true;
		}

		void SetXAxisLabel(const std::string &in_x_ax_label)
		{
			x_ax_label = in_x_ax_label;
			x_ax_label_flag = true;
		}

		void SetYAxisLabel(const std::string &in_y_ax_label)
		{
			y_ax_label = in_y_ax_label;
			y_ax_label_flag = true;
		}

		void SetLegend(const Legend &in_legend)
		{
			legend = in_legend;
			legend_flag = true;
		}

		void SetParams(const PlotParam &other)
		{
			*this = other;
		}

		void Use();

		/// @brief Resets all the flags
		void Reset();
	};

  private:
	/// @brief work directory
	Zaki::String::Directory wrk_dir = "";

	///  @brief If the work directory is set
	bool set_wrk_dir_flag = false;

	size_t default_data_size = 100;
	// std::vector<Row> data_rows ;

	int interp_x_idx = 0;

	gsl_interp_accel *accel = nullptr;

	static DataColumn &NullColumn() noexcept;

	/// @brief Header text
	std::string header_text = "";

	/// @brief Footer text
	std::string footer_text = "";

	SplineSet spline_set;

	PlotParam plt_par;

	// Precision for printing (double) numbers
	int precision = 8;

	/**
	 * @brief The main data storage for the DataSet.
	 *
	 */
	std::vector<DataColumn> data_set;

  public:
	/// @brief Sets the work directory
	void SetWrkDir(const Zaki::String::Directory &);

	/// @brief Returns the work directory
	Zaki::String::Directory GetWrkDir() const;

	// -----------------------------------------------
	// Constructors & Destructors
	// -----------------------------------------------
	/// @brief Default constructor
	DataSet();

	/// @brief Constructor from a data column
	DataSet(const DataColumn &in_dc);

	/// @brief Constructor from a bunch of data columns
	DataSet(const std::vector<DataColumn> &);

	/// @brief Constructor from a file
	DataSet(const Zaki::String::Directory &w_dir,
			const std::string &file,
			const std::vector<size_t> &col_idx = {},
			const char &in_delim = '\t',
			const char &comment_symbol = '#');

	/// @brief Constructor from known sizes of columns and rows
	/// This resizes the number of columns but
	/// only reserves space in rows.
	DataSet(const size_t &column, const size_t &rows);

	/// @brief Destructor
	~DataSet();

	/// @brief Copy Constructor
	DataSet(const DataSet &);
	// -----------------------------------------------

	// -----------------------------
	// Iteration support (range-for)
	// -----------------------------

	/**
	 * @brief Iterator to the first column (const).
	 * Enables `for (const auto& c : ds)` iteration.
	 */
	auto begin() const noexcept { return data_set.begin(); }

	/**
	 * @brief Iterator past the last column (const).
	 */
	auto end() const noexcept { return data_set.end(); }

	/**
	 * @brief Iterator to the first column (mutable).
	 * Enables `for (auto& c : ds)` iteration.
	 */
	auto begin() noexcept { return data_set.begin(); }

	/**
	 * @brief Iterator past the last column (mutable).
	 */
	auto end() noexcept { return data_set.end(); }

	/**
	 * @brief Const iterator to the first column.
	 */
	auto cbegin() const noexcept { return data_set.cbegin(); }

	/**
	 * @brief Const iterator past the last column.
	 */
	auto cend() const noexcept { return data_set.cend(); }

	// -----------------------------------------------
	// Accessors
	// -----------------------------------------------

	/**
	 * @brief Returns the number of columns in the DataSet.
	 *
	 * This counts the number of DataColumn objects currently stored,
	 * independent of the number of rows in each column.
	 *
	 * @return Number of columns.
	 */
	std::size_t ColCount() const noexcept { return data_set.size(); }

	/**
	 * @brief Resolve a column index with support for negative ("Mathematica-style") indexing.
	 *
	 * Indexing rules:
	 *  - `i >= 0`   → zero-based index from the front.
	 *  - `i < 0`    → index from the back (`-1` = last column).
	 *
	 * Examples:
	 *  - `i = 0`   → first column
	 *  - `i = -1`  → last column
	 *  - `i = -N`  → first column (where N = ColCount())
	 *
	 * @param i   Input index (may be negative).
	 * @param out Resolved non-negative index if successful.
	 *
	 * @return `true` if index is valid and resolved; `false` otherwise.
	 *
	 * @note This function never throws and performs no logging.
	 *       Callers decide how to handle out-of-range indices.
	 */
	bool ResolveIndex(int i, std::size_t &out) const noexcept;

	/**
	 * @brief Access a column by index (supports negative indexing).
	 *
	 * Negative indices count from the end, e.g. `-1` refers to the last column.
	 *
	 * @param i Column index (may be negative).
	 *
	 * @return Reference to the requested DataColumn.
	 *
	 * @warning If the index is out of range:
	 *          - In debug builds, this triggers a ZAKI_ASSERT (program abort).
	 *          - In release builds, an error is logged and a safe null column
	 *            is returned as a fallback.
	 *
	 * @see ResolveIndex()
	 */
	DataColumn &Col(int i);

	/**
	 * @brief Access a column by index (const, supports negative indexing).
	 *
	 * Negative indices count from the end, e.g. `-1` refers to the last column.
	 *
	 * @param i Column index (may be negative).
	 *
	 * @return Const reference to the requested DataColumn.
	 *
	 * @warning If the index is out of range:
	 *          - In debug builds, this triggers a ZAKI_ASSERT (program abort).
	 *          - In release builds, an error is logged and a safe null column
	 *            is returned as a fallback.
	 *
	 * @see ResolveIndex()
	 */
	const DataColumn &Col(int i) const;

	/**
	 * @brief Access a column by its label.
	 *
	 * Returns the first column whose label matches the input string.
	 *
	 * @param label Column label to search for.
	 *
	 * @return Reference to the matching DataColumn.
	 *
	 * @warning If no column with the given label exists:
	 *          - In debug builds, this triggers a ZAKI_ASSERT (program abort).
	 *          - In release builds, an error is logged and a safe null column
	 *            is returned as a fallback.
	 */
	DataColumn &Col(const std::string &label);

	/**
	 * @brief Access a column by its label (const).
	 *
	 * Returns the first column whose label matches the input string.
	 *
	 * @param label Column label to search for.
	 *
	 * @return Const reference to the matching DataColumn.
	 *
	 * @warning If no column with the given label exists:
	 *          - In debug builds, this triggers a ZAKI_ASSERT (program abort).
	 *          - In release builds, an error is logged and a safe null column
	 *            is returned as a fallback.
	 */
	const DataColumn &Col(const std::string &label) const;

	/**
	 * @brief Read-only access to the underlying column container.
	 *
	 * This accessor is intended for iteration and inspection only.
	 * Direct modification of the container is not permitted.
	 *
	 * @return Const reference to the internal vector of DataColumn objects.
	 */
	const std::vector<DataColumn> &Columns() const noexcept { return data_set; }

	// -----------------------------------------------
	// Operator overloads (Mathematica-style access)
	// -----------------------------------------------

	/**
	 * @brief Column access operator (supports negative indexing).
	 *
	 * Equivalent to `Col(i)`.
	 */
	DataColumn &operator[](int i) { return Col(i); }

	/**
	 * @brief Column access operator (const, supports negative indexing).
	 *
	 * Equivalent to `Col(i) const`.
	 */
	const DataColumn &operator[](int i) const { return Col(i); }

	/**
	 * @brief Column access operator by label.
	 *
	 * Equivalent to `Col(label)`.
	 */
	DataColumn &operator[](const std::string &label) { return Col(label); }

	/**
	 * @brief Column access operator by label (const).
	 *
	 * Equivalent to `Col(label) const`.
	 */
	const DataColumn &operator[](const std::string &label) const { return Col(label); }

	/**
	 * @brief Extract a subset of columns by index list.
	 *
	 * Each index supports negative ("Mathematica-style") indexing.
	 * The returned columns are **copies** of the original DataColumn objects.
	 *
	 * @param idxs List of column indices to extract.
	 *
	 * @return Vector of selected DataColumn copies.
	 *
	 * @warning If any index is out of range, an error is logged and
	 *          an empty vector is returned.
	 */
	std::vector<DataColumn> operator[](const std::vector<int> &idxs);
	// -----------------------------------------------

	/// @brief Number of rows (max column length).
	std::size_t RowCount() const;

	/// @brief True if there are no columns or no rows.
	bool Empty() const
	{
		if (data_set.empty())
			return true;
		return RowCount() == 0;
	}

	/// @brief  Sets the precision for printing numbers
	/// @param prec number of digits
	void SetPrecision(const int &prec);

	// Returns a pointer to the vec_saver
	// Zaki::File::VecSaver* GetVecSaver() ;

	void AddHead(const std::string &in_txt);

	void AddFoot(const std::string &in_txt);

	/// Appends another dataset data to the end of
	/// the current dataset.
	DataSet Append(const DataSet &) const;

	/// Adds an empty column with the input label
	/// and fills it up with the value
	void AddColumn(const std::string &label, const double &val = 0);

	/// @brief Adds a new DataColumn to the current dataset
	/// @param dc DataColumn to be added
	void AddColumn(DataColumn dc);

	/// Appends a new row to the current dataset
	///  must have the same size as the dataset's width
	void AppendRow(const std::vector<double> &);

	/**
	 * @brief Return the dataset transposed into rows, with optional column selection.
	 *
	 * This method constructs a vector of @ref Row objects, where each Row corresponds
	 * to one logical row of the dataset and contains values drawn from selected columns.
	 *
	 * Column selection is controlled by @p col_idx:
	 *  - If @p col_idx is empty, **all columns** are used in their natural order.
	 *  - If @p col_idx is non-empty, each entry is interpreted as a column index:
	 *      - Non-negative indices select columns from the front.
	 *      - Negative indices are resolved from the end (Mathematica-style indexing).
	 *
	 * Invalid indices are skipped with a warning. If no valid columns remain after
	 * resolution, the function logs an error and returns an empty vector.
	 *
	 * ### Row count policy
	 * The number of returned rows equals the **minimum length** among the selected
	 * columns. This allows heterogeneous column sizes while guaranteeing that
	 * every returned Row is internally consistent.
	 *
	 * ### Performance notes
	 *  - Column indices are resolved once, up front.
	 *  - Raw pointers to the selected DataColumn objects are cached to avoid repeated
	 *    lookups.
	 *  - No dynamic allocation occurs inside the inner row/column loop.
	 *
	 * ### Error handling
	 *  - If the dataset has zero columns, an error is logged and `{}` is returned.
	 *  - If all requested column indices are invalid, an error is logged and `{}` is returned.
	 *  - If the minimum column size is zero, `{}` is returned silently.
	 *
	 * @param col_idx Optional list of column indices to extract. Supports negative indexing.
	 *
	 * @return Vector of Row objects representing the transposed data.
	 *
	 * @warning Returned rows are truncated to the shortest selected column.
	 * @note This function does not modify the dataset.
	 */
	std::vector<Row> GetDataRows(
		const std::vector<int> &col_idx = {}) const;

	/**
	 * @brief Return the dataset transposed into rows using pre-resolved column indices.
	 *
	 * This is a **lower-level, fast variant** of @ref GetDataRows that assumes all column
	 * indices are already resolved, non-negative, and refer directly to valid columns
	 * in the dataset.
	 *
	 * Unlike GetDataRows, this method:
	 *  - Does **not** perform negative-index resolution.
	 *  - Does **not** skip invalid indices.
	 *  - Enforces a stricter contract: all indices must be valid.
	 *
	 * This design makes it suitable for performance-critical paths where column
	 * resolution is done once (e.g., during setup), and row extraction is performed
	 * repeatedly.
	 *
	 * ### Preconditions
	 *  - `cols` must not be empty.
	 *  - Every index in `cols` must satisfy `0 <= idx < ColCount()`.
	 *
	 * If any precondition is violated, an error is logged and `{}` is returned.
	 *
	 * ### Row count policy
	 * The number of returned rows equals the **minimum length** among the selected
	 * columns, identical to @ref GetDataRows.
	 *
	 * ### Performance notes
	 *  - Column pointer lookup is performed once.
	 *  - No index resolution or validation occurs inside inner loops.
	 *  - Intended for internal use or advanced callers.
	 *
	 * @param cols Vector of resolved column indices (non-negative).
	 *
	 * @return Vector of Row objects representing the transposed data.
	 *
	 * @warning Passing invalid indices is a logic error.
	 * @note This function does not modify the dataset.
	 */
	std::vector<Row> GetDataRowsResolved(
		const std::vector<std::size_t> &cols) const;

	/// @brief This is used for reserving memory
	void SetDefaultDataSize(const size_t);

	// -----------------------------------------------
	//                Overloading []
	// -----------------------------------------------
	/// Overloading []
	// DataColumn &operator[](const int);

	/// Overloading []
	// DataColumn operator[](const int) const;

	/// Overloading []: Calling by label
	/// Returns the first column matching the input label
	// DataColumn &operator[](const std::string &in_label);

	/// Overloading [] ( const ): : Calling by label
	/// Returns the first column matching the input label
	// DataColumn operator[](const std::string &in_label) const;

	/// Overloading []
	// std::vector<DataColumn> operator[](const std::vector<int> &);
	// -----------------------------------------------

	/**
	 * @brief Erase a contiguous index range from every column in the dataset (in-place).
	 *
	 * For each column, removes elements in the half-open range [idx_1, idx_2).
	 * This operation is performed independently per column. For columns shorter than
	 * the dataset RowCount(), only the column's own Size() is relevant to validity.
	 *
	 * @param idx_1 Start index of the range to erase (inclusive).
	 * @param idx_2 End index of the range to erase (exclusive).
	 *
	 * @pre idx_1 <= idx_2
	 * @pre For every column c: idx_2 <= c.Size()
	 *      (If you maintain equal-length columns as an invariant, this reduces to
	 *      idx_2 <= RowCount().)
	 *
	 * @note If idx_1 == idx_2, this is a no-op.
	 * @warning Modifies the dataset in-place and invalidates iterators/references
	 *          to erased elements, consistent with std::vector::erase on each column.
	 *
	 * @par Debug checks
	 * In debug builds, the preconditions are enforced via ZAKI_ASSERT.
	 */
	void EraseRange(std::size_t idx_1, std::size_t idx_2);

	/// Reserve uses the resize and reserve methods in std::vector
	/// This is used when the exact number of columns,
	/// and approximate number of rows are known.
	void Reserve(const size_t &s_i, const size_t &s_j);

	// /// Use this when both the number of columns and rows
	// /// are known!
	// void Resize(const size_t& s_i, const size_t& s_j) ;

	///  Resizes the data_Set, and reserves space in spline_set
	/// This is used when only the number of columns is known
	void Resize(const size_t &columns);

	/**
	 * @brief Clears all rows in the dataset and releases allocated memory.
	 *
	 * This method removes all row data from every DataColumn in the DataSet
	 * and **explicitly releases their allocated capacity** back to the system.
	 * Any cached interpolation state (e.g., GSL accelerators and splines)
	 * is also destroyed.
	 *
	 * ### Semantics
	 * - Logical size of all columns becomes zero.
	 * - Underlying storage capacity is released (`shrink_to_fit` / swap idiom).
	 * - Any cached interpolation structures become invalid.
	 *
	 * ### Performance notes
	 * - This operation is **potentially expensive** due to memory deallocation
	 *   and reallocation on subsequent reuse.
	 * - Prefer `ClearRows()` when the DataSet will be rebuilt with a similar
	 *   size or resolution (typical in evolution and integration loops).
	 *
	 * ### Intended use cases
	 * - One-time or infrequent teardown of large datasets.
	 * - Permanent change in resolution or grid structure.
	 * - Memory reclamation before long idle periods or program shutdown.
	 *
	 * @warning All pointers, references, or views into column data become invalid.
	 *
	 * @see ClearRows()
	 */
	void ReleaseRows();

	/**
	 * @brief Clears the column data while preserving allocated capacity.
	 *
	 * This method removes all values from the column but retains the
	 * underlying memory allocation for efficient reuse.
	 *
	 * ### Semantics
	 * - Column size becomes zero.
	 * - Capacity is preserved.
	 *
	 * ### Performance notes
	 * - This is the **preferred method** for resetting columns that will be
	 *   rebuilt repeatedly (e.g., during integrations or evolution steps).
	 *
	 * @see Release()
	 */
	// void Clear() noexcept;
	void ClearRows() noexcept;

	/// @brief Clears the header and footer texts
	void ClearHeadFoot();

	/// @brief Imports a data set from a file
	DataSet *Import(const Zaki::String::Directory &f_name,
					const std::vector<size_t> &col_idx = {},
					const char &in_delim = '\t',
					const char &comment_symbol = '#');

	/// @brief Returns the dimensions of the data set:
	std::vector<size_t> Dim() const;

	/// @brief Exports the data set to a file
	void Export(const Zaki::String::Directory &f_name,
				const std::vector<int> &col_idx = {},
				const char &in_delim = '\t',
				const Zaki::File::FileMode fmode = (Zaki::File::FileMode)1)
		const;

	/// @brief Interpolating y as a function of x
	void Interpolate(const int &x_idx, const int &y_idx);

	/// @brief Interpolating y's as a function of x
	void Interpolate(const int &x_idx, const std::vector<int> &y_idx);

	/// @brief Evaluates the interpolating function (y) given x
	double Evaluate(const int &y_idx, const double &in_x) const;

	/// @brief Evaluates the interpolating function (y) given a list of x's
	std::vector<double> Evaluate(const int &y_idx,
								 const std::vector<double> &in_xs) const;

	/// @brief Evaluates the interpolating function (y) given a DataColumn
	DataColumn Evaluate(const int &y_idx,
						const DataColumn &in_xs) const;

	/// Calculates the derivative of the interpolating function
	/// at a specific point 'x'
	double Derivative(const int &y_idx, const double &in_x) const;

	/// Calculates the derivative of the interpolating function
	/// for all values of 'x', and returns a dataset with two columns
	/// the first column is x, the second is dy/dx
	DataSet Derivative(const int &y_idx) const;

	/// Integrates the interpolating function
	// Zaki::Math::Quantity Integrate(const Zaki::Math::Range<double>&) ;
	double Integrate(const int &y_idx, const Zaki::Math::Range<double> &) const;

	/// Sets integ_rel_err
	// void SetIntegrationRelErr(const double&) ;

	//..................................................
	//                Operator overloading
	//..................................................
	//                    Addition (+)
	//..................................................
	/// Adding a column
	DataSet operator+(const DataColumn &in_dc) const;
	//..................................................
	/// Multiplying by a number
	DataSet operator*(const double &in_num) const;
	//..................................................
	/// Dividing by a number
	DataSet operator/(const double &in_num) const;
	//..................................................
	/// Assignment operator
	DataSet &operator=(const DataSet &other_ds);
	// Assignment operator
	DataSet &operator=(const DataColumn &in_dc);
	// Assignment operator
	DataSet &operator=(const std::vector<DataColumn> &in_dc);
	//..................................................
	/// Addition assignment operator
	DataSet &operator+=(const DataColumn &in_dc);
	//..................................................

	/// @brief  Generates a smoother dataset by moving average
	/// @param window window for moving average
	/// @return Smoother dataset
	DataSet GetSmooth(const short int &window) const;

	/// @brief  Makes the dataset smoother by moving average
	/// @param window window for moving average
	void MakeSmooth(const short int &window);

	// /// @brief  Generates a smoother dataset by removing bumps
	// /// @param window window for moving average
	// /// @return Smoother dataset
	// DataSet RemoveBumps(const short int& window) const ;

	/// @brief  Makes the dataset smoother by removing bumps
	void RemoveBumps();
	// .....................................................
	//              Plotting Methods
	// .....................................................
	// Sets plt_par
	void SetPlotPars(const PlotParam &in_plt_par);

	/// Makes a 2D plot of y vs. x
	void Plot(const int &x_idx, const int &y_idx,
			  const Zaki::String::Directory &f_name,
			  const std::string &in_title = "");

	/// Makes a 2D plot of y vs. x with custom label
	void Plot(const int &x_idx, const std::pair<int, std::string> &y_idx,
			  const Zaki::String::Directory &f_name,
			  const std::string &in_title = "");

	/// Makes a 2D plot of y's vs. x
	void Plot(const int &x_idx,
			  const std::vector<int> &y_idx_set,
			  const Zaki::String::Directory &f_name,
			  const std::string &in_title = "");

	/// Makes a 2D plot of y's vs. x's
	// Added on Sep 21, 2023
	void Plot(const std::vector<std::pair<int, int>> &x_y_idx,
			  const Zaki::String::Directory &f_name,
			  const std::string &in_title = "");

	/// Makes a 2D plot of y's vs. x with custom labels
	void Plot(const int &x_idx,
			  const std::vector<std::pair<int, std::string>> &y_idx_set,
			  const Zaki::String::Directory &f_name,
			  const std::string &in_title = "");

	/// Makes a 2D plot of y's vs. x with style options
	void Plot(const int &x_idx,
			  const std::vector<std::pair<int,
										  std::map<std::string, std::string>>> &y_idx_set,
			  const Zaki::String::Directory &f_name,
			  const std::string &in_title = "");

	/// Makes a 2D log-log-plot of y vs. x
	void LogLogPlot(const int &x_idx, const int &y_idx,
					const Zaki::String::Directory &f_name,
					const std::string &in_title = "");

	/// Makes a 2D log-log-plot of y vs. x with custom label
	void LogLogPlot(const int &x_idx,
					const std::pair<int, std::string> &y_idx,
					const Zaki::String::Directory &f_name,
					const std::string &in_title = "");

	/// Makes a 2D log-log-plot of y's vs. x
	void LogLogPlot(const int &x_idx,
					const std::vector<int> &y_idx_set,
					const Zaki::String::Directory &f_name,
					const std::string &in_title = "");

	/// Makes a 2D log-log-plot of y's vs. x's
	// Added on Sep 21, 2023
	void LogLogPlot(const std::vector<std::pair<int, int>> &x_y_idx,
					const Zaki::String::Directory &f_name,
					const std::string &in_title = "");

	/// Makes a 2D log-log-plot of y's vs. x with custom labels
	void LogLogPlot(const int &x_idx,
					const std::vector<std::pair<int, std::string>> &y_idx_set,
					const Zaki::String::Directory &f_name,
					const std::string &in_title = "");

	/// Makes a 2D log-log-plot of y's vs. x with style options
	void LogLogPlot(const int &x_idx,
					const std::vector<std::pair<int,
												std::map<std::string, std::string>>> &y_idx_set,
					const Zaki::String::Directory &f_name,
					const std::string &in_title = "");

	/// Makes a 2D semi-log-X-plot of y vs. x
	void SemiLogXPlot(const int &x_idx, const int &y_idx,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-X-plot of y vs. x with custom label
	void SemiLogXPlot(const int &x_idx,
					  const std::pair<int, std::string> &y_idx,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-X-plot of y's vs. x
	void SemiLogXPlot(const int &x_idx,
					  const std::vector<int> &y_idx_set,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-X-plot of y's vs. x's
	// Added on Sep 21, 2023
	void SemiLogXPlot(const std::vector<std::pair<int, int>> &x_y_idx,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-X-plot of y's vs. x with custom labels
	void SemiLogXPlot(const int &x_idx,
					  const std::vector<std::pair<int, std::string>> &y_idx_set,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-X-plot of y's vs. x with style options
	void SemiLogXPlot(const int &x_idx,
					  const std::vector<std::pair<int,
												  std::map<std::string, std::string>>> &y_idx_set,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-Y-plot of y vs. x
	void SemiLogYPlot(const int &x_idx, const int &y_idx,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-Y-plot of y vs. x with custom label
	void SemiLogYPlot(const int &x_idx,
					  const std::pair<int, std::string> &y_idx,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-Y-plot of y's vs. x
	void SemiLogYPlot(const int &x_idx,
					  const std::vector<int> &y_idx_set,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-Y-plot of y's vs. x's
	// Added on Sep 21, 2023
	void SemiLogYPlot(const std::vector<std::pair<int, int>> &x_y_idx,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-Y-plot of y's vs. x with custom labels
	void SemiLogYPlot(const int &x_idx,
					  const std::vector<std::pair<int, std::string>> &y_idx_set,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// Makes a 2D semi-log-Y-plot of y's vs. x with style options
	void SemiLogYPlot(const int &x_idx,
					  const std::vector<std::pair<int,
												  std::map<std::string, std::string>>> &y_idx_set,
					  const Zaki::String::Directory &f_name,
					  const std::string &in_title = "");

	/// @brief Resets 'plt_par'
	void ResetPlotPars();
	// .....................................................

	/// Solves for x in 'y(x) = a' using Newton's method
	/// default value for a is '0'
	double Solve(const int &x_idx, const int &y_idx, const double &a = 0);
};

//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::Vector
//--------------------------------------------------------------

#endif /*Zaki_Vector_DataSet_H*/
