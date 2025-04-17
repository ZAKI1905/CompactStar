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
 * @file Math_Core.hpp
 *
 * @brief Core math functions.
 *
 * @ingroup Math
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_Math_Basic_H
#define Zaki_Math_Basic_H

#include <iostream>
#include <vector>

#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

// #include "Zaki/Util/Logger.hpp"
// #include "Zaki/String/Directory.hpp"
// #include "Zaki/File/VecSaver.hpp"
// #include <Zaki/File/CSVIterator.hpp> 

// #include <matplotlibcpp.hpp>
//--------------------------------------------------------------
namespace Zaki::Vector
{
  class DataSet ;
}

//--------------------------------------------------------------
namespace Zaki::Physics
{
  class Coord3D ;
}

//--------------------------------------------------------------
namespace Zaki::String
{
  class Directory ;
}

//--------------------------------------------------------------
namespace Zaki::Math
{

//==============================================================
/// Returns the double factorial for integars
int DoubleFac(int n) ;

//==============================================================
//        Forward Declarations
//==============================================================

class Quantity ; 

//..................................................
//                Addition (+)
//..................................................
/// Addition of a Quantity to a Quantity
Quantity operator+(const Quantity& Q_L, const Quantity& Q_R) ;

/// Addition of a Quantity to a double
Quantity operator+(const double&, const Quantity& Q_R) ;

/// Addition of a double to a Quantity
Quantity operator+(const Quantity& Q_L, const double&) ;

/// Addition of a Quantity to a int
// Quantity operator+(const int&, const Quantity& Q_R) ;

/// Addition of a int to a Quantity
// Quantity operator+(const Quantity& Q_R, const int&) ;
//..................................................
//                Subtraction (-)
//..................................................
/// Subtraction of a Quantity from a Quantity
Quantity operator-(const Quantity& Q_L, const Quantity& Q_R) ;

/// Subtraction of a Quantity from a double
Quantity operator-(const double&, const Quantity& Q_R) ;

/// Subtraction of a double from a Quantity
Quantity operator-(const Quantity& Q_L, const double&) ;

/// Subtraction of a Quantity from a int
// Quantity operator-(const int&, const Quantity& Q_R) ;

/// Subtraction of a int from a Quantity
// Quantity operator-(const Quantity& Q_R, const int&) ;

//..................................................
//              Multiplication (*)
//..................................................
/// Multiplication of a Quantity by a Quantity
Quantity operator*(const Quantity& Q_L, const Quantity& Q_R) ;

/// Multiplication of a Quantity by a double
Quantity operator*(const double&, const Quantity& Q_R) ;

/// Multiplication of a double by a Quantity
Quantity operator*(const Quantity& Q_L, const double&) ;

/// Multiplication of a Quantity by a int
// Quantity operator*(const int&, const Quantity& Q_R) ;

/// Multiplication of a int by a Quantity
// Quantity operator*(const Quantity& Q_R, const int&) ;
//..................................................
//                  Division (/)
//..................................................
/// Division of a Quantity by a Quantity
Quantity operator/(const Quantity& Q_L, const Quantity& Q_R) ;

/// Division of a Quantity by a double
Quantity operator/(const double&, const Quantity& Q_R) ;

/// Division of a double by a Quantity
Quantity operator/(const Quantity& Q_L, const double&) ;

/// Division of a Quantity by a int
// Quantity operator/(const int&, const Quantity& Q_R) ;

/// Division of a int by a Quantity
// Quantity operator/(const Quantity& Q_R, const int&) ;

//..................................................
//               Exponentiation (exp)
//..................................................
// Quantity exp(const Quantity&) ;

//==============================================================
//==============================================================
class Quantity
{
  public:
  //==============================================================
  //        Friend functions
  //..................................................
  friend
  Quantity operator-(const Quantity& Q_L, const Quantity& Q_R) ;

  friend
  Quantity operator-(const double&, const Quantity& Q_R) ;

  friend
  Quantity operator-(const Quantity& Q_R, const double&) ;
  //..................................................
  friend
  Quantity operator+(const Quantity& Q_L, const Quantity& Q_R) ;

  friend
  Quantity operator+(const double&, const Quantity& Q_R) ;

  friend
  Quantity operator+(const Quantity& Q_R, const double&) ;

  // friend
  // Quantity operator+(const int&, const Quantity& Q_R) ;

  // friend
  // Quantity operator+(const Quantity& Q_R, const int&) ;

  //..................................................
  friend
  Quantity operator*(const Quantity& Q_L, const Quantity& Q_R) ;

  friend
  Quantity operator*(const double&, const Quantity& Q_R) ;

  friend
  Quantity operator*(const Quantity& Q_R, const double&) ;

  // friend
  // Quantity operator*(const int&, const Quantity& Q_R) ;

  // friend
  // Quantity operator*(const Quantity& Q_R, const int&) ;
  //..................................................
  
  friend
  Quantity operator/(const Quantity& Q_L, const Quantity& Q_R) ;

  friend
  Quantity operator/(const double&, const Quantity& Q_R) ;

  friend
  Quantity operator/(const Quantity& Q_R, const double&) ;

  // friend
  // Quantity operator/(const int&, const Quantity& Q_R) ;

  // friend
  // Quantity operator/(const Quantity& Q_R, const int&) ;

  //..................................................
  // friend
  // Quantity exp(const Quantity&) ;

  //==============================================================

  long double val, err ;
  int sigFig = 1 ;

  /// Addition
  // Quantity operator+(const Quantity& in_Q) const ;
  // Quantity operator+(const double&) const ;
  Quantity& operator+=(const Quantity& in_Q) ;
  Quantity& operator+=(const double&) ;

  /// Subtraction
  // Quantity operator-(const Quantity& in_Q) const ;
  // Quantity operator-(const double&) const ;
  Quantity& operator-=(const Quantity& in_Q) ;
  Quantity& operator-=(const double&) ;

  /// Multiplication
  // Quantity operator*(const Quantity& in_Q) const ;
  // Quantity operator*(const double&) const ;
  Quantity& operator*=(const Quantity& in_Q) ;
  Quantity& operator*=(const double&) ;

  /// Division
  // Quantity operator/(const Quantity& in_Q) const ;
  // Quantity operator/(const double&) const ;
  Quantity& operator/=(const Quantity& in_Q) ;
  Quantity& operator/=(const double&) ;

  /// power
  Quantity Pow(const double&) const ;


  /// Constructor 0
  Quantity() ;

  /// Constructor 1
  Quantity(const double& in_val, const double& in_err, const int& sigFig=1) ;

  /// Constructor 2
  Quantity(const int& in_val, const double& in_err, const int& sigFig=1) ;

  /// Constructor 3
  Quantity(const double& in_val, const int& in_err, const int& sigFig=1) ;

  /// Constructor 4
  Quantity(const int& in_val, const int& in_err, const int& sigFig=1) ;

  /// Constructor 5
  Quantity(const long double& in_val, const long double& in_err, const int& sigFig=1) ;

  // /// Constructor 5
  // Quantity(const std::vector<double>& in_val_err, const int& sigFig=1) ;
};

//==============================================================
template <typename T>
struct Range
{
  T min, max ;

  Range() : min(0), max(0) {}
  Range(T in_min, T in_max) : min(in_min), max(in_max) {}

  T Len() { return max - min ; }
  T LenAbs() {return abs(max - min) ; }
};

//==============================================================
struct Axis
{
  Range<double> range ; 
  long unsigned int res ;
  std::string scale ;
  double Min() const { return range.min; }
  double Max() const { return range.max; }

  double operator[](const size_t& i) const ;

};

//==============================================================
struct Grid2D
{
  Axis xAxis ;
  Axis yAxis ;
};

//==============================================================


//==============================================================
struct Coord2D
{
  double x, y ;

  Coord2D(const double& in_x, const double& in_y)
        : x(in_x), y(in_y) {}

  /// Constructor from a 3D point
  Coord2D(const Zaki::Physics::Coord3D& in_p) ;


  std::string Str() const ;

  /// Unary minus (negation '-') operator
  Coord2D operator-() const ;

  /// Subtraction 
  Coord2D operator-(const Coord2D&) ;
  
  /// Addition 
  Coord2D operator+(const Coord2D&) ;

  /// Addition assignment operator
  Coord2D& operator+=(const Coord2D&) ;

  /// Subtraction assignment operator
  Coord2D& operator-=(const Coord2D&) ;

  /// Multiplication assignment operator
  Coord2D& operator*=(const double&) ;

  /// Division assignment operator
  Coord2D& operator/=(const double&) ;
};
//==============================================================
//              Coord2D Operator Overloading
//==============================================================
/// Subtraction of a Coord2D from a double
Coord2D operator-(const double&, const Coord2D&) ;

/// Subtraction of a double from a Coord2D
Coord2D operator-(const Coord2D&, const double&) ;

/// Addition of a Coord2D to a double
Coord2D operator+(const double&, const Coord2D&) ;

/// Addition of a double to a Coord2D
Coord2D operator+(const Coord2D&, const double&) ;

/// Multiplication of a Coord2D by a double
Coord2D operator*(const double&, const Coord2D&) ;

/// Multiplication of a double by a Coord2D
Coord2D operator*(const Coord2D&, const double&) ;

/// Division of a Coord2D by a double
Coord2D operator/(const double&, const Coord2D&) ;

/// Division of a double by a Coord2D
Coord2D operator/(const Coord2D&, const double&) ;

//==============================================================
class Segment
{
  public:
    Coord2D p1, p2 ;

    /// Constructor from two points
    Segment(const Coord2D& p1, const Coord2D& p2) ;

    /// Constructor from two 3D points 
    /// (it projects on the first two components)
    Segment(const Zaki::Physics::Coord3D& p1, const Zaki::Physics::Coord3D& p2) ;

    /// Checks if 'x' is in its domain
    bool InDomain(const double& x) const ;

    /// Checks if 'p' is in its domain
    bool InDomain(const Coord2D& p) const ;

    /// Checks if 'y' is in its image
    bool InImage(const double& y) const ;

    /// Checks if 'p' is in its image
    bool InImage(const Coord2D& p) const ;

    /// Returns y(x)
    double y(const double& x) const ;

    /// Returns x(y)
    double x(const double& y) const ;

    /// Returns the slope
    double m() const ;

    /// Returns the coordinate for parameter
    ///  0 <= lambda <= 1
    /// P(0) = P_1, P(1) = P_2
    Coord2D P(const double& lambda) const ;

    /// Is point 'p' above this segment?
    bool IsExclusivelyAboveSeg(const Coord2D& p) const ;

    /// Is point 'p' below this segment?
    bool IsExclusivelyBelowSeg(const Coord2D& p) const ;

    /// Is point 'p' to the left of this segment?
    bool IsExclusivelyLeftSeg(const Coord2D& p) const ;

    /// Is point 'p' to the right of this segment?
    bool IsExclusivelyRightSeg(const Coord2D& p) const ;

    /// Is point 'p' on this segment?
    bool IsOnSeg(const Coord2D& p) const ;

    /// These don't check if intersection is within the domain!
    Coord2D GetIntersection(const Segment& other) const ;
    Coord2D GetIntersection(const Coord2D& p3, const Coord2D& p4) const ;
    Coord2D GetIntersection(const double& x3, const double& y3, 
                    const double& x4, const double& y4) const ;

    bool Intersects(const Segment& other) const ;
    bool Intersects(const Coord2D& p3, const Coord2D& p4) const ;
    bool Intersects(const double& x3, const double& y3, 
                    const double& x4, const double& y4) const ;
} ;

//==============================================================
//                        Curve2D
//==============================================================
class Curve2D
{
  // static_assert(std::is_base_of<Coord2D, T>::value, 
  //               "Curve2D: Class must derive from Coord2D.");

  public:
  std::vector<Coord2D> pts ;
  std::string label = "" ;

  Curve2D(const std::string& in_label="") 
    : label(in_label) 
    {}

  Curve2D(const std::vector<Coord2D>& in_pts, 
          const std::string& in_label="") 
          : pts(in_pts), label(in_label) 
    {}

  /// Sets the label
  void SetLabel(const std::string&) ;

  /// Returns the label
  std::string GetLabel() const ;

  /// Imports the curve points
  void Import(const Zaki::String::Directory& in_file) ;

  void Append(const Coord2D& in_p)
  {
    pts.emplace_back(in_p) ;
  }

  /// Makes a 2D plot of pts
  void Plot(const Zaki::String::Directory& f_name, 
            const std::string& in_title = "") ;

  /// Makes a 2D plot of multiple curves
  static void Plot(const std::vector<Curve2D>& , 
                  const Zaki::String::Directory& f_name, 
                  const std::string& in_title = "") ;

  size_t Size() const ;
  void Reserve(const size_t&) ;

  /// Overloading []
  Coord2D& operator[] (const int) ; 

  /// Overloading [] ( const )
  Coord2D operator[] (const int) const ;

  std::vector<double> GetXVals() const ;
  std::vector<double> GetYVals() const ;

  /// Uses moving average to smooth the curve
  Curve2D GetSmooth(const short int& window) const ;

  /// Uses moving average to smooth the curve
  ///  [ in-place ]
  void MakeSmooth(const short int& window) ;

  /// Transforms the curve into a set of segments
  std::vector<Segment> Segmentize() const ;

  /// Finds all the intersection points
  std::vector<Coord2D> Intersection(const Curve2D& other) const ;

  /// Exports the curve points
  void Export(const Zaki::String::Directory&) const ;
  
  /// Returns the closes index to the input point
  size_t GetIdx(const Coord2D& in_p) const ;

  /// Bisects the curve into two smaller curves
  std::pair<Curve2D, Curve2D> Bisect(const Coord2D& in_p) const ;
};

//--------------------------------------------------------------
// // Deduction guide for deduction of template arguments
// Curve2D() -> Curve2D<Coord2D>;
// Curve2D(const std::string& in_label) -> Curve2D<Coord2D>;
// Curve2D(const std::vector<Coord2D>& in_pts) -> Curve2D<Coord2D>;
// Curve2D(const std::vector<Coord2D>& in_pts, 
//           const std::string& in_label) -> Curve2D<Coord2D>;
// //--------------------------------------------------------------
// // Transforms the curve into a set of segments
// template <class T>
// std::vector<Segment> Curve2D<T>::Segmentize() const 
// {
//   std::vector<Segment> tmp_seg ;
//   tmp_seg.reserve(pts.size()-1) ;

//   for (size_t i = 0; i < pts.size()-1; i++)
//   {
//     tmp_seg.emplace_back(pts[i], pts[i+1]) ;
//   }
  
//   return tmp_seg ;
// }

// //--------------------------------------------------------------
// // Finds all the intersection points
// template <class T>
// std::vector<T> Curve2D<T>::Intersection(const Curve2D<T>& other) const 
// {
//   std::vector<T> inter_set ;

//   std::vector<Segment> this_seg = Segmentize() ;
//   std::vector<Segment> other_seg = other.Segmentize() ;

//   for (size_t i = 0; i < this_seg.size(); i++)
//   {
//     for (size_t j = 0; j < other_seg.size(); j++)
//     {
//       if ( this_seg[i].Intersects(other_seg[j]) )
//         inter_set.emplace_back(this_seg[i].GetIntersection(other_seg[j])) ;
//     }
//   }
  
//   return inter_set ;
// }
// //--------------------------------------------------------------
// template <class T>
// std::pair<Curve2D<T>, Curve2D<T>> Curve2D<T>::Bisect(const T& in_p) const
// {
//   std::size_t const cut_idx = GetIdx(in_p) ;

//   Curve2D split_lo({pts.begin(), pts.begin() + cut_idx},
//                     label+"_1") ;
//   Curve2D split_hi({pts.begin() + cut_idx, pts.end()}, 
//                     label+"_2") ;

//  return {split_lo, split_hi} ;
// }

// //--------------------------------------------------------------
// /// Sets the label
// template <class T>
// void Curve2D<T>::SetLabel(const std::string& in_label) 
// {
//   label = in_label ;
// }

// //--------------------------------------------------------------
// /// Returns the label
// template <class T>
// std::string Curve2D<T>::GetLabel() const 
// {
//   return label ;
// }

// //--------------------------------------------------------------
// // Imports the curve points
// template <class T>
// void Curve2D<T>::Import(const Zaki::String::Directory& in_file)
// {
//   Z_LOG_INFO("Importing file...") ;
//   // .........................................................
//   //                  Opening the file
//   // .........................................................
//   std::ifstream     file(in_file.Str());
//   // Error opening the file
//   if (file.fail()) 
//   {
//     Z_LOG_ERROR("File '"+(in_file).Str() +"' cannot be opened!") ;
//     Z_LOG_ERROR("Importing file failed!") ;
//     exit(EXIT_FAILURE) ;
//   }
//   // .........................................................
//   // Clearing the data_set:
//   pts.clear() ;
//   // .........................................................
//   //                  Reading the file
//   // .........................................................
//   for(Zaki::File::CSVIterator loop(file, '\t'); loop != Zaki::File::CSVIterator(); ++loop)
//   {
//     // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     // Ignore empty lines (spaces and tabs are not empty!)
//     if ((*loop).size() == 0)
//     {
//       // line_num++ ;
//       continue ;
//     }

//     // // Ignore comment (starts with '#') lines
//     // if (Zaki::String::Strip((*loop)[0], ' ')[0] == '#')
//     // {
//     //   // line_num++ ;
//     //   continue ;
//     // }

//     // Read the label (starts with '#') line
//     if (Zaki::String::Strip((*loop)[0], ' ')[0] == '#')
//     {
//       // line_num++ ;
//       label = Zaki::String::Strip((*loop)[1], ' ')[0] ;
//       continue ;
//     }
//     // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     // importing all of the columns:
//     if ((*loop).size() != 2)
//     {
//       Z_LOG_ERROR("Importing file failed: Two columns are needed!") ;
//       exit(EXIT_FAILURE) ;
//     }

//     pts.emplace_back( std::atof((*loop)[0].c_str()), 
//                       std::atof((*loop)[1].c_str()) ) ;
//   }
// }

// //--------------------------------------------------------------
// template <class T>
// size_t Curve2D<T>::Size() const 
// {
//   return pts.size() ;
// }

// //--------------------------------------------------------------
// template <class T>
// void Curve2D<T>::Reserve(const size_t& in_size) 
// {
//   pts.reserve(in_size) ;
// }

// //--------------------------------------------------------------
// /// Overloading []
// template <class T>
// T& Curve2D<T>::operator[] (const int in_i)  
// {
// if (in_i >= 0)
//   {
//     if ( in_i < pts.size() )
//       return pts[(size_t)in_i] ;
//     else
//     {
//       Z_LOG_ERROR("Index out of range!") ;
//       return pts[0] ;
//     }
//   }
//   // Negative index (similar to Mathematica)
//   else
//   {
//     if ( in_i >= -pts.size() )
//       return pts[pts.size()+ in_i] ;
//     else
//     {
//       Z_LOG_ERROR("Index out of range!") ;
//       return pts[0] ;
//     }
//   }
// }

// //--------------------------------------------------------------
// /// Overloading [] ( const )
// template <class T>
// T Curve2D<T>::operator[] (const int in_i) const 
// {
//  if (in_i >= 0)
//   {
//     if ( in_i < pts.size() )
//       return pts[(size_t)in_i] ;
//     else
//     {
//       Z_LOG_ERROR("Index out of range!") ;
//       return pts[0] ;
//     }
//   }
//   // Negative index (similar to Mathematica)
//   else
//   {
//     if ( in_i >= -pts.size() )
//       return pts[pts.size()+ in_i] ;
//     else
//     {
//       Z_LOG_ERROR("Index out of range!") ;
//       return pts[0] ;
//     }
//   }
// }

// //--------------------------------------------------------------
// template <class T>
// std::vector<double> Curve2D<T>::GetXVals() const 
// {
//   std::vector<double> tmp(pts.size()) ;

//   for (size_t i = 0; i < pts.size(); i++)
//   {
//     tmp[i] = pts[i].x ;
//   }
  
//   return tmp ;
// }

// //--------------------------------------------------------------
// template <class T>
// std::vector<double> Curve2D<T>::GetYVals() const 
// {
//   std::vector<double> tmp(pts.size()) ;

//   for (size_t i = 0; i < pts.size(); i++)
//   {
//     tmp[i] = pts[i].y ;
//   }
  
//   return tmp ;
// }

// //--------------------------------------------------------------
// template <class T>
// Zaki::Math::Curve2D<T> Curve2D<T>::GetSmooth(const short int& window) const
// {
//   Zaki::Math::Curve2D out(label+"_smooth") ;
//   out.Reserve(Size()) ;

//   // Zaki::Math::Coord2D moving_ave_val(0,0) ;
//   T moving_ave_val(0,0) ;

//   for (size_t i = 0; i < Size() - window + 1 ; i++)
//   {
//     for (size_t j = 0; j < window; j++)
//     {
//       moving_ave_val += operator[](i + j) ;
//     }
//     moving_ave_val /= window ;

//     out.Append(moving_ave_val) ;

//     moving_ave_val = { 0, 0} ;
//   }
  
//   return out ;
// }

// //--------------------------------------------------------------
// template <class T>
// void Curve2D<T>::MakeSmooth(const short int& window)
// {
//   Zaki::Math::Curve2D<T> out ;
//   out.Reserve(Size()) ;

//   // Zaki::Math::Coord2D moving_ave_val(0,0) ;
//   T moving_ave_val(0,0) ;


//   for (size_t i = 0; i < Size() - window + 1 ; i++)
//   {
//     for (size_t j = 0; j < window; j++)
//     {
//       moving_ave_val += operator[](i + j) ;
//     }
//     moving_ave_val /= window ;

//     out.Append(moving_ave_val) ;

//     moving_ave_val = { 0, 0} ;
//   }
  
//   *this = out ;
// }

// //--------------------------------------------------------------
// // Makes a 2D plot of multiple curves
// template <class T>
// void Curve2D<T>::Plot(const std::vector<Curve2D<T>>& in_curves, 
//           const Zaki::String::Directory& f_name, 
//           const std::string& in_title) 
// {
//   for (auto &&c : in_curves)
//   {
//     matplotlibcpp::plot(c.GetXVals(), c.GetYVals()) ; 
//   }

//   // matplotlibcpp::xlabel(label) ;
//   // matplotlibcpp::ylabel(label) ;

//   if (in_title != "")
//     matplotlibcpp::title(in_title) ;

//   matplotlibcpp::savefig((f_name).Str());
//   matplotlibcpp::clf() ;
// }

// //--------------------------------------------------------------
// // Makes a 2D plot of pts
// template <class T>
// void Curve2D<T>::Plot(const Zaki::String::Directory& f_name, 
//                     const std::string& in_title) 
// {
//   matplotlibcpp::plot(GetXVals(), GetYVals());

//   // matplotlibcpp::xlabel(label) ;
//   // matplotlibcpp::ylabel(label) ;

//   if (in_title != "")
//     matplotlibcpp::title(in_title) ;

//   matplotlibcpp::savefig((f_name).Str());
//   matplotlibcpp::clf() ;
// }

// //--------------------------------------------------------------
// // Exports the curve points
// template <class T>
// void Curve2D<T>::Export(const Zaki::String::Directory& in_file) const
// {
//   Zaki::File::VecSaver v_saver(in_file) ;
//   v_saver.SetHeader("\t" + label) ;
//   v_saver.Export1D(pts) ;
// }

// //--------------------------------------------------------------
// // Returns the index to the closest point to in_p
// template <class T>
// size_t Curve2D<T>::GetIdx(const T& in_p) const
// {
//   double min_dist = (pts[0].x - in_p.x)*(pts[0].x - in_p.x) 
//                     + (pts[0].y - in_p.y)*(pts[0].y - in_p.y) ;

//   size_t tmp_idx = 0 ;
//   for (size_t i = 1; i < pts.size(); i++)
//   {
//     double current_dist = (pts[i].x - in_p.x)*(pts[i].x - in_p.x) 
//                         + (pts[i].y - in_p.y)*(pts[i].y - in_p.y) ;
//     if ( current_dist < min_dist)
//     {
//       tmp_idx = i ;
//       min_dist = current_dist ;
//     }
//   }
  
//   return tmp_idx ;
// }
// //--------------------------------------------------------------

//==============================================================
enum class Excl_Cond
{
  ABOVE = 0, BELOW = 1, ON = 2, ON_ABOVE = 3, ON_BELOW = 4
};

//==============================================================
class Cond_Segment : public Segment
{
  public:
    // Segment seg ;
    Excl_Cond ex_con ; 

    Cond_Segment(const Coord2D& in_p1, 
                const Coord2D& in_p2,
                const Excl_Cond& in_ex_con) ;

    Cond_Segment(const Segment& in_seg, 
    const Excl_Cond& in_ex_con) ;

    bool IsExcluded(const Coord2D& in_p) ;
};

//==============================================================
// Doesn't have to be closed, but it has to be 
//  a convex polygon!
class Cond_Polygon
{
  private:
    std::vector<Cond_Segment> c_segs ;

  public: 
    Cond_Polygon() ;
    Cond_Polygon(const std::vector<Cond_Segment>& in_c_segs) ;

    void AddSegment(const Coord2D& in_p1, 
                    const Coord2D& in_p2,
                    const Excl_Cond& in_ex_con) ;

    void AddSegment(const Segment& in_seg, 
                    const Excl_Cond& in_ex_con) ;

    void AddSegment(const Cond_Segment& in_c_seg) ;

    bool IsExcluded(const Coord2D& in_p) ;

};
//==============================================================
//                GridVals_2D
//==============================================================
struct GridVals_2D
{
  private:
    gsl_spline2d *spline = nullptr ;
    gsl_interp_accel *x_acc = nullptr ;
    gsl_interp_accel *y_acc = nullptr ;
    std::vector<double> x_vals ;
    std::vector<double> y_vals ;
    const gsl_interp2d_type *m_GSL_interp2d_type = gsl_interp2d_bilinear ;

  public:
  double* m_GridValArr = nullptr ;

  
  /// n_x and n_y are the same as 'res' = # of points - 1
  size_t n_x, n_y ;


  /// Default constructor
  /// If there are missing values, the code will abort
  GridVals_2D(const Zaki::Vector::DataSet& in_ds, 
            const int& x_idx, const int& y_idx, 
            const int& z_val_idx) ;

  /// Second constructor
  /// Missing values will be replaced with 'def_val'
  GridVals_2D(const Zaki::Vector::DataSet& in_ds, 
              const int& x_idx, const size_t& x_res, 
              const int& y_idx, const size_t& y_res, 
              const int& z_val_idx, 
              const double& def_val=-1) ;

  ~GridVals_2D() ;

  void Interpolate(const Zaki::Math::Grid2D& in_g2d) ;

  double Evaluate(const double& in_x, const double& in_y) ;
};

//==============================================================
//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------

std::ostream& operator << ( std::ostream &, const Zaki::Math::Quantity&);
std::ostream& operator << (std::ostream &, const Zaki::Math::Range<double>&) ;
std::ostream& operator << (std::ostream &, const Zaki::Math::Range<double>&) ;
std::ostream& operator << (std::ostream &, const Zaki::Math::Range<int>&) ;
std::ostream& operator << ( std::ostream &, const Zaki::Math::Coord2D&);
std::ostream& operator << ( std::ostream &, const Zaki::Math::Segment&);
//==============================================================
#endif /*Zaki_Math_Basic_H*/
