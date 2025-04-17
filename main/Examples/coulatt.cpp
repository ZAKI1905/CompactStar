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

// Local Headers
#include <CompactStar/EOS/Common.hpp>
#include <CompactStar/Core/TOVSolver.hpp>
#include <CompactStar/EOS/CoulombLattice.hpp>

// *******************************************************
int main()
{
  Zaki::String::Directory dir(__FILE__) ;

  CompactStar::CoulombLattice wd ;

  // wd.FindEOS(200) ;
  // wd.SetWrkDir(dir.ParentDir() +"/results") ;

  // wd.ExportEOS("EOS/WD_Carbon.eos") ;


  CompactStar::TOVSolver solver ;

  solver.SetWrkDir(dir.ParentDir() +"/results") ;
  solver.ImportEOS("EOS/WD_Carbon.eos") ;
  
  solver.Solve({{wd.GetValidEDensRange().min*100, wd.GetValidEDensRange().max},
                30, "Log"}, "WhiteDwarf", "WD_Carbon") ;

  // solver.ExportSequence("WhiteDwarf/WD_Carbon_sequence.tsv") ;

  return 0;
}
// *******************************************************