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

#include <CompactStar/EOS/Common.hpp>
#include <CompactStar/Core/TOVSolver.hpp>
#include <CompactStar/EOS/SigmaOmega.hpp>

// *******************************************************
int main()
{
  Zaki::String::Directory dir(__FILE__) ;


  //      Nuclear Matter Parameters
  CompactStar::NuclMatterPar nuc_matt_par ;
  nuc_matt_par.compr = 240 ; // MeV
  nuc_matt_par.eff_m = 0.78 ; 

  //      Sigma - Omega Parameters
  CompactStar::SigmaOmegaPar sig_par ;

  //      Sigma - Omega - Rho Model
  CompactStar::SigmaOmega sig_ome ;

  // sig_ome.SetPars(sig_par.Derive(nuc_matt_par)) ;

  // Using the results in Table 5.5 for K = 240 MeV, and m*/m = 0.78
  // sig_par.Derive(nuc_matt_par).Print() ;
  sig_ome.SetPars({9.927, 4.820, 4.791, 0.8659e-2, -0.2421e-2}) ;

  sig_ome.FindEOS(200) ;
  sig_ome.SetWrkDir(dir.ParentDir() +"/results") ;
  
  sig_ome.ExportEOS("EOS/SigOmeg_240_078.eos") ;


  //      TOV Solver
  CompactStar::TOVSolver solver ;

  solver.SetWrkDir(dir.ParentDir() +"/results") ;
  solver.ImportEOS("EOS/SigOmeg_240_078.eos") ;
  
  solver.Solve({{ sig_ome.GetValidEDensRange().min*2,
                  sig_ome.GetValidEDensRange().max},
                30, "Log"}, "NStar/SigOmeg_240_078", "SigOmeg_240_078") ;

  solver.ExportSequence("NStar/SigOmeg_240_078/SigOmeg_240_078_sequence.tsv") ;


  return 0;
}
// *******************************************************
