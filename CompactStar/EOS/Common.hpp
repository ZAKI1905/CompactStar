// -*- lsst-c++ -*-
/*
* CompactStar
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
 * @file Common.hpp
 *
 * @brief Commonly used functions for EOS.
 *
 * @ingroup EOS
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit Aug 25, 2021
#ifndef CompactStar_Common_H
#define CompactStar_Common_H

#include <math.h>
#include <vector>
// #include <algorithm>
#include <string>

#include <Zaki/String/Directory.hpp>
//--------------------------------------------------------------
namespace CompactStar
{

//==============================================================
/// Input Fermi momentum in fm^-3
double FermiEDens(const double& kF) ;

/// Input Fermi momentum in fm^-3
double FermiPress(const double& kF) ;
//==============================================================
struct NuclMatterPar
{
  /// Binding Energy
  const double bind_e = -16.3 ; // MeV

  /// Saturation Density
  const double sat_den =  0.153 ; // fm^-3

  /// Symmetry Energy Coefficient
  const double a_sym = 32.5  ; // MeV

  /// Compression 
  double compr ;

  /// Effective Mass
  double eff_m ; // no dimension -> divided by MN =  938.9 MeV

  double kf() const
  {
    return pow( 3.*M_PI*M_PI*sat_den/2., 1./3.) ;
  }
};

//==============================================================
// DataSet ColExtractor(const Zaki::String::Directory& f_name,
//   const std::vector<size_t> _idx) ;
//==============================================================  
extern const double MN ; // 1/fm 

//--------------------------------------------------------------
} // End of namespace CompactStar
//--------------------------------------------------------------

#endif /*CompactStar_Common_H*/
