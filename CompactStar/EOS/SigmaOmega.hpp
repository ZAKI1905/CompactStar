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
 * @file SigmaOmega.hpp
 *
 * @brief Sigma-Omega Model
 *
 * @ingroup EOS
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit BEFORE Aug 6
#ifndef CompactStar_SigmaOmega_H
#define CompactStar_SigmaOmega_H

// #include <Zaki/Physics/Constants.hpp>

#include "CompactStar/EOS/Model.hpp"
#include "CompactStar/EOS/SigmaOmegaPar.hpp"

//==============================================================
namespace CompactStar
{

//==============================================================
class SigmaOmega : public Model
{
  //--------------------------------------------------------------
  private:
    SigmaOmegaPar params ;
    
  //--------------------------------------------------------------
  public:
    double rho = 0.1 ; // Baryon number density (not sure yet) 

    SigmaOmega() ;

    ~SigmaOmega() ;
    double EDens(const double& in_rho) override ;
    double Press(const double& in_rho) override ;

    void SetPars(const SigmaOmegaPar& in_pars) ;
    double SolvegsigEquation() ;
    double I1(const double& in_rho, const double& x) ;
    double I2(const double& in_rho, const double& x) ;
    double I3(const double& in_rho, const double& x) ;

    double gsigEquation(double x) ;

};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_SigmaOmega_H*/
