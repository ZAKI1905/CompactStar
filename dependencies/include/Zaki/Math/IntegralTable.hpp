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
 * @file IntegralTable.hpp
 *
 * @brief Table of integrals.
 *
 * @ingroup Math
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_Math_IntegralTable_H
#define Zaki_Math_IntegralTable_H

#include "Zaki/Math/Math_Core.hpp"

//==============================================================
namespace Zaki::Math
{
//==============================================================
/// I_1 integral with integrand :
///     x^2 / sqrt[ x^2 + a^2 ]
double I_1(const Range<double>& limits, const std::vector<double>& pars) ;

/// I_2 integral with integrand :
///     k^2 * sqrt[ k^2 + a^2 ]
double I_2(const Range<double>& limits, const std::vector<double>& pars) ;

/// I_3 integral with integrand :
///     k^4 / sqrt[ k^2 + a^2 ]
double I_3(const Range<double>& limits, const std::vector<double>& pars) ;

//==============================================================
//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------

#endif /*Zaki_Math_IntegralTable_H*/