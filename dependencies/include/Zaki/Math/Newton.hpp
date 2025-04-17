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
 * @file Newton.hpp
 *
 * @brief Solver using Newton's method.
 *
 * @ingroup Math
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

#ifndef Zaki_Math_Solver_Newton_H
#define Zaki_Math_Solver_Newton_H

#include <iostream>
#include <vector>


// #include <matplotlibcpp.hpp>
//--------------------------------------------------------------
namespace Zaki::Math
{

//=======================================================
class Newton
{
private:
	// The function to solve
	double (*func)(const double&) ;

	// The derivative of function
	double (*der)(const double&) ;

	// Number of iterations
	size_t iter_num = 100 ;

	// Relative precision
	double prec = 1e-10 ;

	// Initial guess
	double x_0_guess ;

public:
	Newton() ;
	~Newton() ;

	void SetFunction(double (*f)(const double&) ) ;
	void SetDerivative(double (*f)(const double&) ) ;
	void SetGuessValue(const double& in_x) ;

	double Solve() ;
};

//=======================================================
//--------------------------------------------------------------
} // End of namespace Math::Solver
//--------------------------------------------------------------

//==============================================================
#endif /*Zaki_Math_Solver_Newton_H*/
