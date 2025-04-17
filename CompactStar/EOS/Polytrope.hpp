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
 * @file Polytrope.hpp
 *
 * @brief Polytrope Model
 *
 * @ingroup EOS
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit BEFORE Aug 6
#ifndef CompactStar_Polytrope_H
#define CompactStar_Polytrope_H

#include "CompactStar/EOS/Model.hpp"

//==============================================================
namespace CompactStar
{

//==============================================================
class Polytrope : public Model
{
  //--------------------------------------------------------------
  private:
    double K ;
    double gamma ;
    
  //--------------------------------------------------------------
  public:

    Polytrope() ;

    ~Polytrope() ;

    double EDens(const double&) override ;
    double Press(const double&) override ;
    
    double GetRho(const double& in_e) const ;

    void SetParams(const double& in_K, const double& in_gamma) ;

};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Polytrope_H*/
