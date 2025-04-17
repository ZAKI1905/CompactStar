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
 * @file Analysis.hpp
 *
 * @brief Virtual class for analyzing the NS in TOVSolver.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef CompactStar_Analysis_H
#define CompactStar_Analysis_H

// #include <vector>
// #include <gsl/gsl_spline.h>

// #include <Zaki/Math/Math_Core.hpp>

#include "CompactStar/Core/Prog.hpp"

//==============================================================
namespace CompactStar
{

class MixedStar ;
class NStar ;
//==============================================================
class Analysis : public Prog
{
  protected:
    std::string label ;
  //--------------------------------------------------------------
  public:
    
    Analysis();

    ~Analysis() ;
    
    /// Analysis during the sequence loop
    virtual void Analyze(MixedStar* in_star) = 0 ;
    virtual void Analyze(NStar* in_star) = 0 ;

    /// Saves the results
    virtual void Export(const Zaki::String::Directory&) = 0 ;

    // Sets the label
    virtual void SetLabel(const std::string& in_label) ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Analysis_H*/
