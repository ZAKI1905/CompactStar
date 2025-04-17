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
 * @file Func2D.hpp
 *
 * @brief Virtual class for a function of two variables.
 *
 * @ingroup Math
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_Math_Func2D_H
#define Zaki_Math_Func2D_H

#include <memory>

//==============================================================
namespace Zaki::Math
{

//==============================================================
class GeneralObject
{
  public:
    GeneralObject() ;
    virtual GeneralObject* Clone() = 0 ;
    virtual ~GeneralObject() ;
};
//==============================================================
class Func2D
{
  protected:
    virtual Func2D* IClone() const = 0 ;
    // {
    //   return nullptr ;
    // };
    // virtual Func2D* IThis() = 0 ;


  public:
    Func2D() ;

    virtual ~Func2D() ;
    std::unique_ptr<Func2D> Clone() const 
    { return std::unique_ptr<Func2D>(IClone()); }

    // std::unique_ptr<Func2D> This() 
    // { return std::unique_ptr<Func2D>(IThis()); }

    virtual double Eval(const double x, const double y) = 0 ;
    std::string PtrStr() const ;
    // virtual void* GetObj() = 0 ;
    // virtual void UpdateObj(void*) = 0 ;

    // virtual GeneralObject* GetGenObj() = 0;

};
//==============================================================

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------

#endif /*Zaki_Math_Func2D_H*/