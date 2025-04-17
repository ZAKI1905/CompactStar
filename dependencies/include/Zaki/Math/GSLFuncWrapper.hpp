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
 * @file GSLFuncWrapper.hpp
 *
 * @brief Function wrapper for gsl functions.
 *
 * @ingroup Math
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_Math_GSLFuncWrapper_H
#define Zaki_Math_GSLFuncWrapper_H

#include <gsl/gsl_math.h>

// Zaki::Util
#include "Zaki/Util/Logger.hpp"

//--------------------------------------------------------------
namespace Zaki::Math
{

//==============================================================
template<typename FuncObj, typename MemFuncPtr >
class GSLFuncWrapper : public gsl_function 
{
  public:
    /// Default Constructor
    GSLFuncWrapper()
    {
      function = &GSLFuncWrapper::invoke;
      params=this;
    }

    /// Constructor from Pointer-2-Obj and the member-function
    GSLFuncWrapper(FuncObj* objPtr, const MemFuncPtr& memFn)
      : fObjPtr(objPtr), fMemFunc( memFn )
    {
      function = &GSLFuncWrapper::invoke;
      params=this;
    }

    void SetMemberFunc(FuncObj* objPtr, const MemFuncPtr& memFn) 
    {
      fObjPtr   = objPtr ;
      fMemFunc  = memFn  ;
    }

  private:
    FuncObj* fObjPtr = nullptr ;
    MemFuncPtr fMemFunc = nullptr ;
    double Eval(const double x) { return (fObjPtr->*fMemFunc)(x) ;}

    static double invoke(double x, void *params) {
      return static_cast<GSLFuncWrapper*>(params)->Eval(x);
  }
};

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------
//==============================================================
#endif /*Zaki_Math_GSLFuncWrapper_H*/
