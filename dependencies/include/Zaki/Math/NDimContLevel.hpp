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
 * @file NDimContLevel.hpp
 *
 * @brief Finds the N-dimensional Gaussian contour levels.
 *
 * @ingroup Math
 *
 * @author Mohammadreza Zakeri 
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_Math_NDimContLevel_H
#define Zaki_Math_NDimContLevel_H

// Zaki::Math
#include "Zaki/Math/Math_Core.hpp"

//==============================================================
namespace Zaki::Math
{

//==============================================================
class NDimContLevel
{
  private:

    double PSeries(const double n) const;

    double P(const double n) const ;

    double tolerance = 0.00001 ;
    int dim = 1 ;
    double conf_level = 1 ; // 1 sigma


  public:

    /// Default constructor
    NDimContLevel() ;

    /// First constructor
    NDimContLevel(const int in_Dim, const double in_CL) 
     : dim(in_Dim), conf_level(in_CL) {} ;

    /// Sets the tolerance for the solver
    void SetTolerance(const double) ;
    
    /// Sets the dimension
    void SetDim(const int) ;
    
    /// Sets the Confidence Level C.L.
    void SetCL(const double) ;

    struct Solution
    {
      int Dim ; double CL, P, Err;
      double Val, Up;
    };

    /// Solves for Delta(L) corresponding to the confidence level
    Solution Solve() ;

};

//==============================================================
//                  ostream << overloading
std::ostream& operator << (std::ostream &output,
const Zaki::Math::NDimContLevel::Solution&) ;
//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------

#endif /*Zaki_Math_NDimContLevel_H*/