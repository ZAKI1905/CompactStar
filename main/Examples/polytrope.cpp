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

// Local Headers
#include <CompactStar/EOS/Common.hpp>
#include <CompactStar/Core/TOVSolver.hpp>
#include <CompactStar/EOS/Polytrope.hpp>

using namespace Zaki::Physics ;

// *******************************************************
int main()
{
  Zaki::String::Directory dir(__FILE__) ;

  CompactStar::Polytrope pol ;

  // Taken from page 91 ("NONRELATIVISTIC ELECTRON REGION")
  // In this limit gamma = 5/3
  // and we are below the neutron threshold:
  // rho < rho_n = 7.4919 e-9  [ fm^-3 ]
  pol.SetRhoRange({1e-13, 7.4919e-9}) ;

  double mn = 938.91897*MEV_2_INV_FM ;
  double g = 5./3.;
  double k = pow(3*M_PI*M_PI/mn, g)/15./M_PI/M_PI/ELECTRON_M_FM;

  pol.SetParams(k, g) ;

  pol.FindEOS(200) ;
  pol.SetWrkDir(dir.ParentDir() +"/results") ;

  pol.ExportEOS("EOS/Polytrope_5_3rd.eos") ;


  CompactStar::TOVSolver solver ;

  solver.SetWrkDir(dir.ParentDir() +"/results") ;
  solver.ImportEOS("EOS/Polytrope_5_3rd.eos") ;
  
  solver.Solve({{pol.GetValidEDensRange().min*100, pol.GetValidEDensRange().max},
                30, "Log"}, "Polytrope", "Pol_5_3rd") ;

  // solver.ExportSequence("Polytrope/Pol_5_3rd_sequence.tsv") ;

  return 0;
}
// *******************************************************