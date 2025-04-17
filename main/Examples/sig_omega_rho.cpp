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

// Last edited on Dec 15, 2020
#include <CompactStar/EOS/Common.hpp>
#include <CompactStar/Core/TOVSolver.hpp>
#include <CompactStar/EOS/SigmaOmegaRho.hpp>

// *******************************************************
int main()
{
  Zaki::String::Directory dir(__FILE__) ;


  //      Nuclear Matter Parameters
  // CompactStar::NuclMatterPar nuc_matt_par ;
  // nuc_matt_par.compr = 240 ; // MeV
  // nuc_matt_par.eff_m = 0.78 ; 

  //      Sigma - Omega Parameters
  // CompactStar::SigmaOmegaPar sig_par ;
  // sig_par.Derive(nuc_matt_par).Print() ;

  //      Sigma - Omega - Rho Model
  CompactStar::SigmaOmegaRho sig_ome_rho ;

  // sig_ome_rho.SetPars(sig_par.Derive(nuc_matt_par)) ;

  // Using the results in Table 5.5 for K = 240 MeV, and m*/m = 0.78
  sig_ome_rho.SetPars({9.927, 4.820, 4.791, 0.8659e-2, -0.2421e-2}) ;

  //........... Checking the equation solvers .........
//  sig_ome_rho.SetRho(1.5) ;
//  sig_ome_rho.UpdatePool(1.5) ;
//  sig_ome_rho.SolveEq() ;
  //.....................................................
  // std::cout << sig_ome_rho.Press(0.5) << "\n" ;

  sig_ome_rho.FindEOS(200) ;
  sig_ome_rho.SetWrkDir(dir.ParentDir() +"/results") ;

  sig_ome_rho.ExportEOS("EOS/SigOmegRho_240_078_Hyp_dmu.eos") ;

  //      TOV Solver
  //  CompactStar::TOVSolver solver ;

  //  solver.SetWrkDir(dir.ParentDir() +"/results") ;
  //  solver.ImportEOS("EOS/SigOmegRho_240_078_Hyp.eos") ;
  
  //  solver.Solve({{ sig_ome_rho.GetValidEDensRange().min*2,
  //                  sig_ome_rho.GetValidEDensRange().max},
  //                200, "Log"}, "NStar/SigOmegRho_240_078_Hyp_nu'/SigOmegRho_240_078_Hyp_*.tsv") ;

  //  solver.ExportSequence("NStar/SigOmegRho_240_078_Hyp_nu'/SigOmegRho_240_078_Hyp_sequence.tsv") ;

  //  solver.ExportNu("NStar/SigOmegRho_240_078_Hyp_nu'/SigOmegRho_240_078_Hyp_10.tsv") ;

  return 0;
}
// *******************************************************
