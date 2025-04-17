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
 * @file Model.hpp
 *
 * @brief Model base class for writing EoSs
 *
 * @ingroup EOS
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit BEFORE Aug 6
#ifndef CompactStar_Model_H
#define CompactStar_Model_H

#include <vector>
// #include <gsl/gsl_spline.h>

#include <Zaki/Math/Math_Core.hpp>

#include "CompactStar/Core/Prog.hpp"

//==============================================================
namespace CompactStar
{

//==============================================================
class Model : public Prog
{
  //--------------------------------------------------------------
  protected:
    std::vector<std::vector<double>> eos_results ;
    Zaki::Math::Range<double> valid_rho ;
  //--------------------------------------------------------------
  public:
    
    /// Input valid ranges of rho
    Model(const Zaki::Math::Range<double>& in_valid_rho);

    ~Model() ;
    
    void SetRhoRange(const Zaki::Math::Range<double>&) ;

    virtual double EDens(const double& in_rho) = 0 ;
    virtual double Press(const double& in_rho) = 0 ;

    /// This will generate a row of values for EOS
    virtual std::vector<double> EOSRow(const double&) ;

    /// This will generate the header for EOS
    virtual std::string EOSHeader() const ;
    
    Zaki::Math::Range<double> GetValidEDensRange();

    void ExportEOS(const Zaki::String::Directory&) const ;

    /// Input number of EOS points
    void FindEOS(const size_t& eos_pts) ;

    std::vector<std::vector<double>> GetEOS() const ;

};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Model_H*/
