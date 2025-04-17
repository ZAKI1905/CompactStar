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
 * @file GSLMultiFWrapper.hpp
 *
 * @brief Multi-variable function wrapper for gsl functions.
 *
 * @ingroup Math
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_Math_GSLMultiFWrapper_H
#define Zaki_Math_GSLMultiFWrapper_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

// Zaki::Util
#include "Zaki/Util/Logger.hpp"

//--------------------------------------------------------------
namespace Zaki::Math
{

//==============================================================
template<typename FuncObj, typename MemFuncPtr >
class GSLMultiFWrapper : public gsl_multiroot_function 
{
  public:
    /// Default Constructor
    GSLMultiFWrapper(const size_t& in_n)
    {
      f = &GSLMultiFWrapper::invoke;
      n = in_n ;
      params   = this;
    }

    /// Constructor from Pointer-2-Obj and the member-function
    GSLMultiFWrapper(FuncObj* objPtr, const MemFuncPtr& memFn, const size_t& in_n)
      : fObjPtr(objPtr), fMemFunc( memFn )
    {
      f = &GSLMultiFWrapper::invoke;
      n = in_n ;
      params = this;
    }

    void SetMemberFunc(FuncObj* objPtr, const MemFuncPtr& memFn) 
    {
      fObjPtr   = objPtr ;
      fMemFunc  = memFn  ;
    }

  private:
    FuncObj* fObjPtr = nullptr ;
    MemFuncPtr fMemFunc = nullptr ;

    int Eval(const gsl_vector* x, gsl_vector* f) 
    { 
      (fObjPtr->*fMemFunc)(x, f) ;

      return GSL_SUCCESS ;
    }

    static int invoke(const gsl_vector * x, void * p, gsl_vector * f) 
    {
      return static_cast<GSLMultiFWrapper*>(p)->Eval(x, f);
    }
};

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------
//==============================================================
#endif /*Zaki_Math_GSLMultiFWrapper_H*/
