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
 * @file GSLMultiFdfWrapper.hpp
 *
 * @brief Multi-variable function & derivative wrapper for gsl functions.
 *
 * @ingroup Math
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_Math_GSLMultiFdfWrapper_H
#define Zaki_Math_GSLMultiFdfWrapper_H

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>

// Zaki::Util
#include "Zaki/Util/Logger.hpp"

//--------------------------------------------------------------
namespace Zaki::Math
{

//==============================================================
template< typename FuncObj, 
          typename MemFuncPtr, 
          typename MemDfPtr, 
          typename MemFdfPtr
        >
class GSLMultiFdfWrapper : public gsl_multiroot_function_fdf 
{
  public:
    /// Default Constructor
    GSLMultiFdfWrapper(const size_t& in_n)
    {
      f   = &GSLMultiFdfWrapper::invoke_f ;
      df  = &GSLMultiFdfWrapper::invoke_df ;
      fdf = &GSLMultiFdfWrapper::invoke_fdf ;

      n = in_n ;
      params   = this;

    }

    /// Constructor from Pointer-2-Obj and the member-function
    GSLMultiFdfWrapper(FuncObj* objPtr, const MemFuncPtr& memFn, 
                       const MemDfPtr& memDf, const MemFdfPtr& memFdf,
                       const size_t& in_n)
      : fObjPtr(objPtr), fMemFunc( memFn ), 
        dfMemFunc(memDf), fdfMemFunc(memFdf)
    {
      f   = &GSLMultiFdfWrapper::invoke_f ;
      df  = &GSLMultiFdfWrapper::invoke_df ;
      fdf = &GSLMultiFdfWrapper::invoke_fdf ;

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

    /// Function (f) method
    MemFuncPtr fMemFunc   = nullptr ;

    /// Derivative Function (df) method
    MemDfPtr dfMemFunc  = nullptr ;

    /// Function & Derivative Function (fdf) method
    MemFdfPtr fdfMemFunc = nullptr ;

    int Eval_f(const gsl_vector* x, gsl_vector* f) 
    { 
      (fObjPtr->*fMemFunc)(x, f) ;

      return GSL_SUCCESS ;
    }

    static int invoke_f(const gsl_vector * x, void * p, gsl_vector * f) 
    {
      return static_cast<GSLMultiFdfWrapper*>(p)->Eval_f(x, f);
    }

    int Eval_df(const gsl_vector* x, gsl_matrix* J) 
    { 
      (fObjPtr->*dfMemFunc)(x, J) ;

      return GSL_SUCCESS ;
    }

    static int invoke_df(const gsl_vector * x, void * p, gsl_matrix * J) 
    {
      return static_cast<GSLMultiFdfWrapper*>(p)->Eval_df(x, J);
    }

    int Eval_fdf(const gsl_vector* x, gsl_vector*  f, gsl_matrix* J) 
    { 
      (fObjPtr->*fdfMemFunc)(x, f, J) ;

      return GSL_SUCCESS ;
    }

    static int invoke_fdf(const gsl_vector * x, void * p, gsl_vector* f, gsl_matrix * J) 
    {
      return static_cast<GSLMultiFdfWrapper*>(p)->Eval_fdf(x, f, J);
    }
};

//--------------------------------------------------------------
} // End of namespace Zaki::Math
//--------------------------------------------------------------
//==============================================================
#endif /*Zaki_Math_GSLMultiFdfWrapper_H*/
