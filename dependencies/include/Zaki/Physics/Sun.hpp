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
 * @file Sun.hpp
 *
 * @brief Clocks and Date-time definitions and conversions.
 *
 * @ingroup Physics
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 * 
 * Created on: May 2 2020
*/

#ifndef Zaki_Physics_Sun_hpp
#define Zaki_Physics_Sun_hpp

#include "Zaki/Physics/Coordinate.hpp"

//--------------------------------------------------------------
namespace Zaki::Physics
{
    
//==============================================================
//class Sun
//{
    
//public:
    /// Sun position in GEI coordinates
    Zaki::Physics::GEICoord GetSunPos(const double&) ;

    /// Earth position in Heliocentric Ecliptic Coordinate
    Zaki::Physics::HEclipticCoord GetEarthPos(const double&) ;
    
    ///  The ecliptic longitude of the Sun
    double GetSunEclipticLong(const double&) ;
//};
//==============================================================

} // End of namespace Zaki::Physics
//--------------------------------------------------------------
#endif /* Zaki_Physics_Sun_hpp */
