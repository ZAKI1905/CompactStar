// -*- lsst-c++ -*-
/*
* Zaki's Common Library
* See License file at the top of the source tree.
*
* Copyright (c) 2023 Mohammadreza Zakeri
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

#include <string>
#include <sstream>

// For SplineSet
#include <map>

#include <Zaki/String/Directory.hpp>
#include <Zaki/Math/Math_Core.hpp>
#include <Zaki/File/VecSaver.hpp>

//--------------------------------------------------------------
namespace Zaki::Math
{

// Math forward declarations
// template <typename T>
// struct Range ;

// template<typename FuncObj, typename MemFuncPtr >
// class GSLFuncWrapper ;

}

//--------------------------------------------------------------
namespace Zaki::File
{

// File forward declarations
enum class FileMode ;

}

//--------------------------------------------------------------
namespace Zaki::Vector
{

// Vector forward declarations
struct DataColumn ; 
class DataSet    ;   

//==============================================================
//..................................................
//                Addition (+)
//..................................................
/// Addition of a DataColumn to a DataColumn
DataColumn operator+(const DataColumn&, const DataColumn&) ;

/// Addition of a DataColumn to a double
DataColumn operator+(const double&, const DataColumn&) ;

/// Addition of a double to a DataColumn
DataColumn operator+(const DataColumn&, const double&) ;

/// Addition of a list of doubles to a DataColumn
DataColumn operator+(const DataColumn&, const std::vector<double>&) ;

/// Addition of a list of doubles to a DataColumn
DataColumn operator+(const std::vector<double>&, const DataColumn&) ;

//..................................................
//                Subtraction (-)
//..................................................
/// Subtraction of a DataColumn from a DataColumn
DataColumn operator-(const DataColumn&, const DataColumn&) ;

/// Subtraction of a DataColumn from a double
DataColumn operator-(const double&, const DataColumn&) ;

/// Subtraction of a double from a DataColumn
DataColumn operator-(const DataColumn&, const double&) ;

/// Subtraction of a list of doubles from a DataColumn
DataColumn operator-(const DataColumn&, const std::vector<double>&) ;

/// Subtraction of a list of doubles from a DataColumn
DataColumn operator-(const std::vector<double>&, const DataColumn&) ;

//..................................................
//              Multiplication (*)
//..................................................
/// Multiplication of a DataColumn by a double
DataColumn operator*(const DataColumn&, const double& ) ;

/// Multiplication of a double by a DataColumn
DataColumn operator*(const double& , const DataColumn&) ;

/// Multiplication of a DataColumn by an int
DataColumn operator*(const DataColumn&, const int& )  ;

/// Multiplication of an int by a DataColumn
DataColumn operator*(const int&, const DataColumn& )  ;

/// Multiplication of a DataColumn by another DataColumn
DataColumn operator*(const DataColumn&, const DataColumn& ) ;
//..................................................
//                  Division (/)
//..................................................

/// Division of a DataColumn by a double
DataColumn operator/(const DataColumn&, const double& ) ;

/// Division of a double by a DataColumn
DataColumn operator/(const double& , const DataColumn&) ;

/// Division of a DataColumn by an int
DataColumn operator/(const DataColumn&, const int& )  ;

/// Division of an int by a DataColumn
DataColumn operator/(const int&, const DataColumn& )  ;

/// Division of a DataColumn by another DataColumn
DataColumn operator/(const DataColumn&, const DataColumn& ) ;

//..................................................
//               Exponentiation (exp)
//..................................................
DataColumn exp(const DataColumn&) ;
DataColumn log(const DataColumn&) ;
DataColumn log10(const DataColumn&) ;

//==============================================================

//==============================================================
//                        DataColumn Struct
//==============================================================
struct DataColumn
{
  private:
    // in case of missing data (NaN) this will be used to fill
    // the data:
    // double add_def_val = 0 ;

  public:

    //..................................................
    friend
    DataColumn operator+(const DataColumn&, const DataColumn&) ;
    friend
    DataColumn operator+(const double&, const DataColumn&) ;
    friend
    DataColumn operator+(const DataColumn&, const double&) ;
    friend
    DataColumn operator+(const DataColumn&, const std::vector<double>&) ;
    friend
    DataColumn operator+(const std::vector<double>&, const DataColumn&) ;
    //..................................................
    friend
    DataColumn operator-(const DataColumn&, const DataColumn&) ;
    friend
    DataColumn operator-(const double&, const DataColumn&) ;
    friend
    DataColumn operator-(const DataColumn&, const double&) ;
    friend
    DataColumn operator-(const DataColumn&, const std::vector<double>&) ;
    friend
    DataColumn operator-(const std::vector<double>&, const DataColumn&) ;
    //..................................................
    friend
    DataColumn operator*(const DataColumn&, const double& ) ;
    friend
    DataColumn operator*(const double& , const DataColumn&) ;
    friend
    DataColumn operator*(const DataColumn&, const int& ) ;
    friend
    DataColumn operator*(const int& , const DataColumn&) ;
    friend
    DataColumn operator*(const DataColumn&, const DataColumn& ) ;
    //..................................................
    friend
    DataColumn operator/(const DataColumn&, const double& ) ;
    friend
    DataColumn operator/(const double& , const DataColumn&) ;
    friend
    DataColumn operator/(const DataColumn&, const int& ) ;
    friend
    DataColumn operator/(const int& , const DataColumn&) ;
    friend
    DataColumn operator/(const DataColumn&, const DataColumn& ) ;
    //..................................................
    friend
    DataColumn exp(const DataColumn&) ;
    DataColumn log(const DataColumn&) ;
    DataColumn log10(const DataColumn&) ;
    //..................................................
    /// Default constructor
    DataColumn() ;
    
    /// Constructor from a known size for rows
    DataColumn(const size_t& rows) ;

    /// Constructor from a label and size and values
    /// by default it fills the column with zeros
    DataColumn( const std::string& label, 
                const size_t& size, 
                const double& val=0) ;

    /// Constructor from a vector and label:
    DataColumn(const std::string&, const std::vector<double>&) ;

    std::string label ;
    std::vector<double> vals ;
    
    /// Fills all of the elements of the data column with the input value.
    void Fill(const double& in_val) ;

    /// Fills specific elements of the data column with the input 
    /// value if the input condition function returns true.
    void Fill(const double& in_val, bool (*fill_condition)(const double&) ) ;

    /// @brief  Generates a smoother datacolumn by moving average
    /// @param window window for moving average
    /// @return Smoother datacolumn
    DataColumn GetSmooth(const short int& window) const ;

    /// @brief  Makes the datacolumn smoother by moving average
    /// @param window window for moving average
    void MakeSmooth(const short int& window) ;

    /// @brief  Makes the datacolumn smoother by removing bumps
    void RemoveBumps() ;
    //..................................................
    // Operator overloading
    //..................................................
    // Addition (+)
    // DataColumn operator+(const DataColumn& in_dc) const ;

    /// Adds the two DataColumns and fills the non-existent 
    /// elements with '0'
    /// Use when the sizes of DC's are not the same
    static DataColumn Add(const DataColumn& dc_1, 
                          const DataColumn& dc_2, 
                          const double& fill=0) ;

    /// Addition to a list of numbers (+)
    // DataColumn operator+(const std::vector<double>& in_vec) const ;
    /// Addition to a single number (+)
    // DataColumn operator+(const double& in_num) const ;
    //..................................................
    /// Addition assignment operator
    DataColumn& operator+=(const DataColumn& in_dc) ;
    DataColumn& operator+=(const double& in_num) ;

    /// Subtraction assignment operator
    DataColumn& operator-=(const DataColumn& in_dc) ;
    DataColumn& operator-=(const double& in_num) ;

    /// Multiplication assignment operator
    DataColumn& operator*=(const DataColumn& in_dc) ;
    DataColumn& operator*=(const double& in_num) ;

    /// Division assignment operator
    DataColumn& operator/=(const DataColumn& in_dc) ;
    DataColumn& operator/=(const double& in_num) ;

    //..................................................
    // // Subtraction (-)
    // DataColumn operator-(const DataColumn& in_dc) const ;
    // // Subtraction of a single number (-)
    // DataColumn operator-(const double& in_num) const ;
    // // Subtraction of a list of numbers (-)
    // DataColumn operator-(const std::vector<double>& in_vec) const ;
    //..................................................
    /// Unary minus (negation '-') operator
    DataColumn operator-() const ;
    //..................................................

    // Multiplication (*)
    // // By a double
    // DataColumn operator*(const double& in_num) const ;
    // // By an int
    // DataColumn operator*(const int& in_num) const ;
    // // By another column
    // DataColumn operator*(const DataColumn& in_dc) const ;
    //..................................................
    // Division (/)
    // // by a number
    // DataColumn operator/(const double& in_nm) const ;
    // // By an int
    // DataColumn operator/(const int& in_num) const ;
    // // By another column
    // DataColumn operator/(const DataColumn& in_dc) const ;
    //..................................................
    /// Power 
    DataColumn pow(const double& in_num) const ;
    // DataColumn pow^(const int& in_num) const ;
    /// sqrt 
    DataColumn sqrt() const  ;
    /// absolute value 
    DataColumn Abs() const  ;
    // // Exponentiation 
    // DataColumn exp(const DataColumn&) const ;
    //..................................................
    /// Overloading []
    double& operator[] (const int) ; 

    /// Overloading [] ( const )
    double operator[] (const int) const ;
    //..................................................

    DataSet CombineColumns(const DataColumn&) const;
    DataSet CombineColumns(const std::vector<DataColumn>&) const;

    /// Selecting a subset of DataColumn elements satisfying 
    ///  the input condition 'cond'
    DataColumn GetSubSet(bool (*cond)(const double&)) const ; 

    /// Selecting a subset of DataColumn elements satisfying 
    ///  the input condition 'cond'
    DataColumn GetSubSet(bool (*cond)(const DataColumn&, 
                                      const int& idx)
                          ) const ; 

    /// Selecting a subset of DataColumn ranging from
    ///  idx_1 --> idx_2
    DataColumn GetSubSet(const size_t idx_1, const size_t idx_2) const ;  

    /// Trims the data column to its subset ranging from
    ///  idx_1 --> idx_2
    void Trim(const size_t idx_1, const size_t idx_2)  ;

    /// Returns the index to :
    ///  If the data is:
    ///  Ascending: the first value that is equal or greater than
    ///  the input value. 
    ///  Descending: the first value that is equal or less than
    ///  the input value. 
    ///  The data column must be strictly increasing or decreasing. 
    int GetFirstIdx(const double&) const ;

    /// Returns the index to the element that is closest. 
    int GetClosestIdx(const double&) const ;

    /// Returns the size of 'vals'
    size_t Size() const;

    /// Reserves space for 'vals'
    void Reserve(const size_t&) ;

    /// Resizes the 'vals'
    void Resize(const size_t&) ;

    /// Returns the minimum element
    double Min() const ;

    /// Returns the maximum element
    double Max() const ;

    /// Returns the minimum element's index
    int MinIdx() const ;

    /// Returns the maximum element's index
    int MaxIdx() const ;
}; 

//==============================================================
//                        PlotData2D Struct
//==============================================================
// struct PlotData2D
// {
//   public:
//     DataColumn x ;
//     DataColumn y ;
// };

//==============================================================
//                        Row Struct
//==============================================================
struct Row
{
  public:
    std::vector<double> vals ;
    char delim = '\t' ;
    
    // Precision for printing (double) numbers
    int precision = 8 ; 

    Row() {} 

    Row(const std::vector<double>& in_vals ) 
    {
      vals = in_vals ;
    }

  std::string Str() const
  {
    std::stringstream ss;

    // [Apr 12, 2023] : Delim issue was resolved.
    for (size_t i = 0; i < vals.size() - 1 ; i++)
    {
      char tmp[150] ;
      snprintf(tmp, sizeof(tmp), "%-*.*e%c ", 7 + precision,
                    precision, vals[i], delim) ;
      ss << tmp ;
    }



    // Last column
    char tmp[150] ;
    snprintf(tmp, sizeof(tmp), "%-*.*e ", 7 + precision, 
              precision, vals[vals.size() - 1]) ;
    ss << tmp ;

    return ss.str() ;
  }
};

//==============================================================
//                        SplineIdx Struct
//==============================================================
struct SplineIdx
{
  int  idx  = -1 ;
  gsl_spline *spline = nullptr ;

  // static inline std::atomic<size_t> def_con_count = 0;
  // static inline std::atomic<size_t> sec_con_count = 0;
  // static inline std::atomic<size_t> mov_con_count = 0;
  // static inline std::atomic<size_t> mov_ass_count = 0;
  // static inline std::atomic<size_t> des_count = 0;
  
  /// The number of times destructor is called
  // size_t des_counter = 0 ;
  // const gsl_interp_type* dataset_gsl_interp_type = gsl_interp_steffen ;
  const gsl_interp_type* dataset_gsl_interp_type = gsl_interp_linear ;

  SplineIdx() 
  {
    // def_con_count++ ;
  } ;

  SplineIdx(const int& in_idx, const double *x,
                  const double *y, const size_t& in_size)
    : idx(in_idx)
    {
      spline = gsl_spline_alloc (dataset_gsl_interp_type, in_size) ;
      gsl_spline_init (spline, x, y, in_size) ;
      // sec_con_count++ ;
    }

  ~SplineIdx()
  {
    // des_counter++ ;
    // std::cout << "SplineIdx destructor called for index = "
    //           << idx <<  " for the'"
    //           << des_counter
    //           << "' time!\n" ;
    // des_count++ ;

    // if( 1 == def_con_count + sec_con_count - des_count )
    // {
    //   std::cout << "\n\n\n Default Constructors = " << def_con_count
    //             << ", Second Constructor = "  << sec_con_count
    //             << ", Move Constructor = " << mov_con_count 
    //             << ", Move Assignment = " << mov_ass_count
    //             << ", Destructors = " << des_count << std::endl ;
    // }

    if (spline)
    {
      gsl_spline_free (spline);
      spline = nullptr ;
    }
  }

  void Reset(const double *x, const double *y, const size_t& in_size)
  {
    if (spline)
    {
      gsl_spline_free (spline);
    }

    spline = gsl_spline_alloc (dataset_gsl_interp_type, in_size) ;
    gsl_spline_init (spline, x, y, in_size) ;
  }

  // Equality operator
  bool operator==(const int& in_idx) const
  {
    return in_idx == idx ;
  }

  /// Copy Constructor
  SplineIdx(const SplineIdx &) = delete ;

  /// Assignment operator
  SplineIdx& operator=(const SplineIdx&) = delete ;

  /// Move constructor
  SplineIdx(SplineIdx&& other) 
    : idx(other.idx),
      spline( other.spline )  
  {
    // mov_con_count++ ;
    other.spline = nullptr ;
  }

  /// Move Assignment
  SplineIdx& operator=(SplineIdx&& other) 
  {
    // mov_ass_count++ ;
    idx = other.idx ;
    spline = other.spline ;
    other.spline = nullptr ;
    return *this ;
  }

};

//==============================================================
//                        SplineSet Struct
//==============================================================
struct SplineSet
{
  std::vector<SplineIdx> spline_set ;
  // std::map<int, SplineIdx>  spline_map ;

  // void AddSpline( const int& in_idx, const double *x,
  //                 const double *y, const size_t& in_size) ;

  void AssignSpline( const int& in_idx, const double *x,
                  const double *y, const size_t& in_size) ;

  double Evaluate(const int& in_idx, const double& in_x, 
                  gsl_interp_accel* in_accel) const ;

  double Derivative(const int& in_idx, const double& in_x, 
                    gsl_interp_accel* in_accel)  const;

  double Integrate(const int& in_idx, const Zaki::Math::Range<double>& in_range, 
                    gsl_interp_accel* in_accel) const ;

  /// Default constructor
  SplineSet() {}

  /// Copy Constructor
  SplineSet(const SplineSet &) = delete ;

  // Assignment operator
  SplineSet& operator=(const SplineSet&) = delete ;

  void Resize(const size_t& in_size)
  {
    spline_set.resize(in_size) ;
  }
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
        size_t w = 800 ;

        /// Height of the figure
        size_t h = 600 ;

        /// Constructor
        FigSize(const size_t& in_w, const size_t& in_h)
        : w(in_w), h(in_h) {}
      };
      // ----------------------------------------------
      
      // ----------------------------------------------
      /// Axis Ticks dimensions
      struct AxisTicks
      {
        /// Location of ticks
        std::vector<double> ticks ;

        /// Label of ticks
        std::vector<std::string> labels ;

        /// Labels are set or no?
        bool label_flag = false ;

        /// Constructor-0
        AxisTicks() {}

        /// Constructor-1 from tick locations
        AxisTicks(const std::vector<double>& in_ticks)
        : ticks(in_ticks) {}
        
        /// Constructor-2 from tick locations & labels
        AxisTicks(const std::vector<double>& in_ticks,
                  const std::vector<std::string>& in_labels)
        : ticks(in_ticks), labels(in_labels) 
        {
          if (in_ticks.size() == in_labels.size())
          {
            label_flag = true ;
          }
          else
          {
            Z_LOG_ERROR("The sizes of ticks and labels should match!") ;
            Z_LOG_INFO("Axis labels are ignored.") ;
          }
        }

        void SetTicks(const std::vector<double>& in_ticks)
        {
          ticks = in_ticks ;
        }

        void SetTickLabels(const std::vector<std::string>& in_labels)
        {
          labels = in_labels ;
          label_flag = true ;
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
        std::string loc = "best" ;

        /// @brief  The coordinate of the bottom-left corner
        double x0 = 0, y0 = 0 ;

        /// @brief Width of the BBox
        double w = 0 ;

        /// @brief Height of the BBox
        double h = 0 ;

        // std::map<std::string, std::string> &keywords = {} ;

        /// @brief Default Constructor
        Legend() {}

        /// @brief Constructor-1
        Legend(const double& in_x0, const double& in_y0, 
               const double& in_w=0, const double& in_h=0)
        : x0(in_x0), y0(in_y0), w(in_w), h(in_h) 
        {}

        /// @brief Constructor-2
        Legend(const std::string& in_loc, 
               const double& in_x0, const double& in_y0, 
               const double& in_w=0, const double& in_h=0)
        : loc(in_loc), x0(in_x0), y0(in_y0), w(in_w), h(in_h) 
        {}
      };
      // ----------------------------------------------

      Zaki::Math::Range<double> x_ax = {0,0} ;
      Zaki::Math::Range<double> y_ax = {0,0} ;
      FigSize fig_size = {800, 600} ;
      bool grid = false ;
      AxisTicks x_ticks ;
      AxisTicks y_ticks ;

      std::string x_ax_label ;
      std::string y_ax_label ;
      std::vector<
        std::pair<double, 
          std::map<std::string, std::string>
                  >> ax_h_lines ;

      std::vector<
        std::pair<double, 
          std::map<std::string, std::string>
                  >> ax_v_lines ;

      Legend legend ;
      // ---------------------------
      //         Flags
      // ---------------------------
      bool x_ax_flag = false ;
      bool y_ax_flag = false ;

      bool fig_size_flag = false ;
      bool grid_flag = false ;
      bool x_ticks_flag = false ;
      bool y_ticks_flag = false ;
      bool x_ax_label_flag = false ;
      bool y_ax_label_flag = false ;
      bool legend_flag = false ;
      bool ax_h_line_flag = false ;
      bool ax_v_line_flag = false ;

      public :
      PlotParam() {}

      void AddAxHLine(const double& y, 
        const std::map<std::string, std::string>& options = {})
      {
        ax_h_lines.emplace_back(y, options) ;
        ax_h_line_flag = true ;
      }

      void AddAxVLine(const double& x, 
        const std::map<std::string, std::string>& options = {})
      {
        ax_v_lines.emplace_back(x, options) ;
        ax_v_line_flag = true ;
      }

      void SetXAxis(const Zaki::Math::Range<double>& in_x_ax)
      {
        x_ax = in_x_ax ;
        x_ax_flag = true ;
      }

      void SetYAxis(const Zaki::Math::Range<double>& in_y_ax)
      {
        y_ax = in_y_ax ;
        y_ax_flag = true ;
      }

      void SetFigSize(const FigSize& in_fig_Size)
      {
        fig_size = in_fig_Size ;
        fig_size_flag = true ;
      }

      void SetGrid(const bool& in_grid=true)
      {
        grid = in_grid ;
        grid_flag = true ;
      }

      void SetXTicks(const AxisTicks& in_x_ticks)
      {
        x_ticks = in_x_ticks ;
        x_ticks_flag = true ;
      }

      void SetYTicks(const AxisTicks& in_y_ticks)
      {
        y_ticks = in_y_ticks ;
        y_ticks_flag = true ;
      }

      void SetXAxisLabel(const std::string& in_x_ax_label)
      {
        x_ax_label = in_x_ax_label ;
        x_ax_label_flag = true ;
      }

      void SetYAxisLabel(const std::string& in_y_ax_label)
      {
        y_ax_label = in_y_ax_label ;
        y_ax_label_flag = true ;
      }


      void SetLegend(const Legend& in_legend)
      {
        legend = in_legend ;
        legend_flag = true ;
      }

      void SetParams(const PlotParam& other )
      {
        *this = other ;
      }

      void Use() ;

      /// @brief Resets all the flags
      void Reset() ;
    };

  private:

    /// @brief work directory
    Zaki::String::Directory wrk_dir = "" ;
    
    ///  @brief If the work directory is set
    bool set_wrk_dir_flag = false ;

    size_t default_data_size = 100 ;
    // std::vector<Row> data_rows ;

    int  interp_x_idx = 0 ;

    gsl_interp_accel *accel = nullptr ;

    // Used for saving the dataset
    // Zaki::File::VecSaver vec_saver ;
    
    /// @brief Header text
    std::string header_text = "" ;

    /// @brief Footer text
    std::string footer_text = "" ;

    SplineSet spline_set ;

    // gsl_integration_workspace *integ_workspace = nullptr ;
    // Zaki::Math::GSLFuncWrapper<DataSet, double (DataSet::*)(const double& )>
    //          *wrapper_GSL = nullptr ;
    // gsl_function* func_GSL  = nullptr ;

    // double integ_rel_err = 1e-6 ;
    
    PlotParam plt_par ;

    // Precision for printing (double) numbers
    int precision = 8 ;

  public:

    /// @brief Sets the work directory
    void SetWrkDir(const Zaki::String::Directory&) ;

    /// @brief Returns the work directory
    Zaki::String::Directory GetWrkDir() const ;

    std::vector<DataColumn> data_set ;

    /// @brief Default constructor
    DataSet() ;

    /// @brief Constructor from a data column
    DataSet(const DataColumn& in_dc) ;

    /// @brief Constructor from a bunch of data columns
    DataSet(const std::vector<DataColumn>&) ;

    /// @brief Constructor from a file
    DataSet(const Zaki::String::Directory& w_dir, 
            const std::string& file, 
            const std::vector<size_t>& col_idx = {}, 
            const char& in_delim='\t',
            const char& comment_symbol='#') ;

    /// @brief Constructor from known sizes of columns and rows
    /// This resizes the number of columns but
    /// only reserves space in rows.
    DataSet(const size_t& column, const size_t& rows) ;

    /// @brief Destructor
    ~DataSet() ;

    /// @brief Copy Constructor
    DataSet(const DataSet &) ;

    /// @brief  Sets the precision for printing numbers
    /// @param prec number of digits
    void SetPrecision(const int& prec) ;

    // Returns a pointer to the vec_saver
    // Zaki::File::VecSaver* GetVecSaver() ;

    void AddHead(const std::string& in_txt) ;

    void AddFoot(const std::string& in_txt) ;

    /// Appends another dataset data to the end of
    /// the current dataset. 
    DataSet Append(const DataSet&) const ;

    /// Adds an empty column with the input label
    /// and fills it up with the value
    void AddColumn(const std::string& label, const double& val=0) ;

    /// Appends a new row to the current dataset
    ///  must have the same size as the dataset's width
    void AppendRow(const std::vector<double>&) ;


    /// @brief Gets the transpose of data (only selected columns in col_idx)
    std::vector<Row> GetDataRows(
      const std::vector<int>& col_idx = {}) const ;

    /// @brief This is used for reserving memory
    void SetDefaultDataSize(const size_t) ;

    // -----------------------------------------------
    //                Overloading []
    // -----------------------------------------------
    /// Overloading []
    DataColumn& operator[] (const int) ; 

    /// Overloading []
    DataColumn operator[] (const int) const ; 

    /// Overloading []: Calling by label
    /// Returns the first column matching the input label
    DataColumn& operator[] (const std::string& in_label) ; 

    /// Overloading [] ( const ): : Calling by label
    /// Returns the first column matching the input label
    DataColumn operator[] (const std::string& in_label) const ;

    /// Overloading []
    std::vector<DataColumn> operator[] (const std::vector<int>&) ; 
    // -----------------------------------------------

    /// Trims the data set to its subset ranging from
    ///  idx_1 --> idx_2
    void Trim(const size_t idx_1, const size_t idx_2) ;  

    /// Reserve uses the resize and reserve methods in std::vector
    /// This is used when the exact number of columns,
    /// and approximate number of rows are known.
    void Reserve(const size_t& s_i, const size_t& s_j) ;
    
    // /// Use this when both the number of columns and rows
    // /// are known!
    // void Resize(const size_t& s_i, const size_t& s_j) ;

    ///  Resizes the data_Set, and reserves space in spline_set
    /// This is used when only the number of columns is known
    void Resize(const size_t& columns) ;

    /// It clears the rows in data_set columns,
    ///  but keeps the columns! 
    void ClearRows() ;

    /// @brief Clears the header and footer texts
    void ClearHeadFoot() ;
    
    /// @brief Imports a data set from a file
    DataSet* Import(const Zaki::String::Directory& f_name, 
                    const std::vector<size_t>& col_idx = {}, 
                    const char& in_delim='\t',
                    const char& comment_symbol='#') ;

    /// @brief Returns the dimensions of the data set:
    std::vector<size_t> Dim() const ;   

    /// @brief Exports the data set to a file
    void Export(const Zaki::String::Directory& f_name, 
                const std::vector<int>& col_idx = {}, 
                const char& in_delim='\t',
                const Zaki::File::FileMode fmode=(Zaki::File::FileMode)1) 
                const ;

    /// @brief Interpolating y as a function of x
    void Interpolate(const int& x_idx, const int& y_idx) ;

    /// @brief Interpolating y's as a function of x
    void Interpolate(const int& x_idx, const std::vector<int>& y_idx) ;

    /// @brief Evaluates the interpolating function (y) given x
    double Evaluate(const int& y_idx, const double& in_x) const ;

    /// @brief Evaluates the interpolating function (y) given a list of x's
    std::vector<double> Evaluate(const int& y_idx, 
                        const std::vector<double>& in_xs) const ;

    /// @brief Evaluates the interpolating function (y) given a DataColumn
    DataColumn Evaluate(const int& y_idx, 
                        const DataColumn& in_xs) const ;

    /// Calculates the derivative of the interpolating function
    /// at a specific point 'x'
    double Derivative(const int& y_idx, const double& in_x) const ;

    /// Calculates the derivative of the interpolating function
    /// for all values of 'x', and returns a dataset with two columns
    /// the first column is x, the second is dy/dx
    DataSet Derivative(const int& y_idx) const ;

    /// Integrates the interpolating function
    // Zaki::Math::Quantity Integrate(const Zaki::Math::Range<double>&) ;
    double Integrate(const int& y_idx, const Zaki::Math::Range<double>&) const ;

    /// Sets integ_rel_err
    // void SetIntegrationRelErr(const double&) ;

    //..................................................
    //                Operator overloading
    //..................................................
    //                    Addition (+)
    //..................................................
    /// Adding a column
    DataSet operator+(const DataColumn& in_dc) const ;
    //..................................................
    /// Multiplying by a number
    DataSet operator*(const double& in_num) const ;
    //..................................................
    /// Dividing by a number
    DataSet operator/(const double& in_num) const ;
    //..................................................
    /// Assignment operator
    DataSet& operator=(const DataSet& other_ds) ;
    // Assignment operator
    DataSet& operator=(const DataColumn& in_dc) ;
    // Assignment operator
    DataSet& operator=(const std::vector<DataColumn>& in_dc) ;
    //..................................................
    /// Addition assignment operator
    DataSet& operator+=(const DataColumn& in_dc) ;
    //..................................................

    /// @brief  Generates a smoother dataset by moving average
    /// @param window window for moving average
    /// @return Smoother dataset
    DataSet GetSmooth(const short int& window) const ;

    /// @brief  Makes the dataset smoother by moving average
    /// @param window window for moving average
    void MakeSmooth(const short int& window) ;

    // /// @brief  Generates a smoother dataset by removing bumps
    // /// @param window window for moving average
    // /// @return Smoother dataset
    // DataSet RemoveBumps(const short int& window) const ;

    /// @brief  Makes the dataset smoother by removing bumps
    void RemoveBumps() ;
    // .....................................................
    //              Plotting Methods
    // .....................................................
    // Sets plt_par
    void SetPlotPars(const PlotParam& in_plt_par) ;

    /// Makes a 2D plot of y vs. x
    void Plot(const int& x_idx, const int& y_idx,
              const Zaki::String::Directory& f_name, 
              const std::string& in_title = "") ;

    /// Makes a 2D plot of y vs. x with custom label
    void Plot(const int& x_idx, const std::pair<int, std::string>& y_idx,
              const Zaki::String::Directory& f_name, 
              const std::string& in_title = "") ;

    /// Makes a 2D plot of y's vs. x
    void Plot(const int& x_idx, 
        const std::vector<int>& y_idx_set,
        const Zaki::String::Directory& f_name, 
        const std::string& in_title="") ;
    
    /// Makes a 2D plot of y's vs. x's
    // Added on Sep 21, 2023
    void Plot(const std::vector< std::pair<int, int>> & x_y_idx, 
        const Zaki::String::Directory& f_name, 
        const std::string& in_title="") ;

    /// Makes a 2D plot of y's vs. x with custom labels
    void Plot(const int& x_idx, 
        const std::vector<std::pair<int, std::string>>& y_idx_set,
        const Zaki::String::Directory& f_name, 
        const std::string& in_title="") ;

    /// Makes a 2D plot of y's vs. x with style options
    void Plot(const int& x_idx, 
            const std::vector< std::pair<int, 
              std::map<std::string, std::string>>>& y_idx_set,
            const Zaki::String::Directory& f_name, 
            const std::string& in_title="") ;

    /// Makes a 2D log-log-plot of y vs. x
    void LogLogPlot(const int& x_idx, const int& y_idx,
              const Zaki::String::Directory& f_name, 
              const std::string& in_title = "") ;

    /// Makes a 2D log-log-plot of y vs. x with custom label
    void LogLogPlot(const int& x_idx, 
              const std::pair<int, std::string>& y_idx,
              const Zaki::String::Directory& f_name, 
              const std::string& in_title = "") ;

    /// Makes a 2D log-log-plot of y's vs. x
    void LogLogPlot(const int& x_idx, 
        const std::vector<int>& y_idx_set,
        const Zaki::String::Directory& f_name, 
        const std::string& in_title= "") ;
    
    /// Makes a 2D log-log-plot of y's vs. x's
    // Added on Sep 21, 2023
    void LogLogPlot(const std::vector< std::pair<int, int>> & x_y_idx, 
        const Zaki::String::Directory& f_name, 
        const std::string& in_title="") ;

    /// Makes a 2D log-log-plot of y's vs. x with custom labels
    void LogLogPlot(const int& x_idx, 
        const std::vector<std::pair<int, std::string>>& y_idx_set,
        const Zaki::String::Directory& f_name, 
        const std::string& in_title= "") ;

    /// Makes a 2D log-log-plot of y's vs. x with style options
    void LogLogPlot(const int& x_idx, 
            const std::vector< std::pair<int, 
              std::map<std::string, std::string>>>& y_idx_set,
            const Zaki::String::Directory& f_name, 
            const std::string& in_title="") ;

    /// Makes a 2D semi-log-X-plot of y vs. x
    void SemiLogXPlot(const int& x_idx, const int& y_idx,
              const Zaki::String::Directory& f_name, 
              const std::string& in_title = "") ;

    /// Makes a 2D semi-log-X-plot of y vs. x with custom label
    void SemiLogXPlot(const int& x_idx, 
              const std::pair<int, std::string>& y_idx,
              const Zaki::String::Directory& f_name, 
              const std::string& in_title = "") ;

    /// Makes a 2D semi-log-X-plot of y's vs. x
    void SemiLogXPlot(const int& x_idx, 
        const std::vector<int>& y_idx_set,
        const Zaki::String::Directory& f_name, 
        const std::string& in_title= "") ;

    /// Makes a 2D semi-log-X-plot of y's vs. x's
    // Added on Sep 21, 2023
    void SemiLogXPlot(const std::vector< std::pair<int, int>> & x_y_idx, 
        const Zaki::String::Directory& f_name, 
        const std::string& in_title="") ;

    /// Makes a 2D semi-log-X-plot of y's vs. x with custom labels
    void SemiLogXPlot(const int& x_idx, 
        const std::vector<std::pair<int, std::string>>& y_idx_set,
        const Zaki::String::Directory& f_name, 
        const std::string& in_title= "") ;

    /// Makes a 2D semi-log-X-plot of y's vs. x with style options
    void SemiLogXPlot(const int& x_idx, 
            const std::vector< std::pair<int, 
              std::map<std::string, std::string>>>& y_idx_set,
            const Zaki::String::Directory& f_name, 
            const std::string& in_title="") ;

    /// Makes a 2D semi-log-Y-plot of y vs. x
    void SemiLogYPlot(const int& x_idx, const int& y_idx,
              const Zaki::String::Directory& f_name, 
              const std::string& in_title = "") ;

    /// Makes a 2D semi-log-Y-plot of y vs. x with custom label
    void SemiLogYPlot(const int& x_idx, 
              const std::pair<int, std::string>& y_idx,
              const Zaki::String::Directory& f_name, 
              const std::string& in_title = "") ;

    /// Makes a 2D semi-log-Y-plot of y's vs. x
    void SemiLogYPlot(const int& x_idx, 
        const std::vector<int>& y_idx_set,
        const Zaki::String::Directory& f_name, 
        const std::string& in_title= "") ;

    /// Makes a 2D semi-log-Y-plot of y's vs. x's
    // Added on Sep 21, 2023
    void SemiLogYPlot(const std::vector< std::pair<int, int>> & x_y_idx, 
        const Zaki::String::Directory& f_name, 
        const std::string& in_title="") ;

    /// Makes a 2D semi-log-Y-plot of y's vs. x with custom labels
    void SemiLogYPlot(const int& x_idx, 
        const std::vector<std::pair<int, std::string>>& y_idx_set,
        const Zaki::String::Directory& f_name, 
        const std::string& in_title= "") ;
    
    /// Makes a 2D semi-log-Y-plot of y's vs. x with style options
    void SemiLogYPlot(const int& x_idx, 
            const std::vector< std::pair<int, 
              std::map<std::string, std::string>>>& y_idx_set,
            const Zaki::String::Directory& f_name, 
            const std::string& in_title="") ;

    /// @brief Resets 'plt_par'
    void ResetPlotPars() ;
    // .....................................................

    /// Solves for x in 'y(x) = a' using Newton's method
    /// default value for a is '0'
    double Solve(const int& x_idx, const int& y_idx, const double& a=0) ;
};

//==============================================================  

//--------------------------------------------------------------
} // End of namespace Zaki::Vector
//--------------------------------------------------------------

#endif /*Zaki_Vector_DataSet_H*/
