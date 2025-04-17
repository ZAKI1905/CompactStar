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

#include <Zaki/Physics/Constants.hpp>

#include <CompactStar/EOS/Baryon.hpp>

using namespace Zaki::Physics ;

// *******************************************************
int main()
{
  Zaki::String::Directory dir(__FILE__) ;

  // CompactStar::Baryon p(+1, 938.27208816*MEV_2_FM, +1./2., 1, 1, 1) ;

  // // Using the results in Table 5.5 for K = 240 MeV, and m*/m = 0.78
  // p.SetPars({9.927, 4.820, 4.791, 0.8659e-2, -0.2421e-2}) ;
  // p.SetRho(0.1) ;

  // std::cout << " p.k()= "<< p.k() << "\n" ;

  return 0;
}
// *******************************************************
