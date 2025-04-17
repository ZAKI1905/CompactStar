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
 * @file Fermi_Gas.hpp
 *
 * @brief Fermi gas Model
 *
 * @ingroup EOS
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Created on Apr 26, 2022
#ifndef CompactStar_Fermi_Gas_H
#define CompactStar_Fermi_Gas_H

#include "CompactStar/EOS/Model.hpp"

//==============================================================
namespace CompactStar
{

//==============================================================
class Fermi_Gas : public Model
{
  //--------------------------------------------------------------
  private:

    /// Mass of the fermion (in MeV)
    double m ;

    // We may implement temperature later.
  //--------------------------------------------------------------
  public:

    /// Constructor: Mass must be in MeV
    Fermi_Gas(const double& in_m) ;

    ~Fermi_Gas() ;

    double EDens(const double&) override ;
    double Press(const double&) override ;
    
    /// Mass must be in MeV
    void SetParams(const double& in_m) ;

};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Fermi_Gas_H*/
