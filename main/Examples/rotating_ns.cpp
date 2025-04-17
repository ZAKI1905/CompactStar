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

// Last edited on Sep 8, 2021
// Local Headers
#include <CompactStar/EOS/Common.hpp>
#include <CompactStar/Core/RotationSolver.hpp>

#include <matplotlibcpp.hpp>

// *******************************************************
int main()
{
  Zaki::String::Directory dir(__FILE__) ;

#if 1
  // ----------------------------------------------------------------------------------------
  CompactStar::RotationSolver r_solver ;

  r_solver.SetWrkDir(dir.ParentDir() +"/results") ;
  // r_solver.ImportTOVSolution("NStar/SigOmegRho_240_078_Hyp/SigOmegRho_240_078_Hyp_29.tsv") ;
  // r_solver.ImportTOVSolution("NStar/SigOmegRho_240_078_Hyp_nu'/SigOmegRho_240_078_Hyp_29.tsv") ;

  // r_solver.ImportTOVSolution("NStar/Glendenning/Glendenning_Table_5-8_29.tsv") ;


  // r_solver.GetPressDer(0.01) ;

  // r_solver.GetNu(1) ;
  // r_solver.GetNu(15) ;

  // ................................................................
  // Fastest known pulsar: "PSR J1748 - 244ad"
  // P = 1.39595482 ms ---> Omega = 4.501 * 10^3 (s^-1)
  // Source: 
  // " Hessels, J. W. T. (2006).
  //   A Radio Pulsar Spinning at 716 Hz. 
  //   Science, 311(5769), 1901–1904. doi:10.1126/science.1123430 "
  // ................................................................
  r_solver.Solve({{ 5.66379367363 * 1e-4, 5.66379367363 * 1e-3},
                3, "Log"}, "RNStar/Glendenning_Table_5-8_29/Glen_5-8_29_Omega_*.tsv") ;

  r_solver.ExportResults("RNStar/Glendenning_Table_5-8_29/Glen_5-8_29_Omega_sequence.tsv") ;

  // r_solver.Solve({{ 0.01, 0.1},
  //               30, "Log"}, "RNStar/SigOmegRho_240_078_Hyp_29/SOR_240_78_Hy_29_Omega_*.tsv") ;

  // r_solver.ExportResults("RNStar/SigOmegRho_240_078_Hyp_29/SOR_240_78_Hyp_29_Omega_sequence.tsv") ;
  // ----------------------------------------------------------------------------------------
#else
  CompactStar::PlotData out = CompactStar::ColExtractor(dir.ParentDir() + "/results/RNStar/Glendenning_Table_5-8_29/Glen_5-8_29_Omega_1.tsv", 1, 2) ;


  // *********************************************************************
  // Checking the universal function omega(r)/Omega and comparing it to
  // Fig. 6.2 of Glendenning (page 286)
  // y-axis : omega(r)/Omega
  // x-axis: r/R_star
  // for (size_t i = 0; i < out.x_data.size(); i++)
  // {
  //   out.x_data[i] = out.x_data[i] / 1.024392e+01 ;
  //   out.y_data[i] = out.y_data[i] / 7.946961e+02 ;
  // }
  // out.x_label = "r / R" ;
  // out.y_label = "\u03C9(r) / \u03A9" ;
  // *********************************************************************

  // Set the size of output image = 1200x900 pixels
  matplotlibcpp::figure_size(1200, 900);

  // Plot line from given x and y data. Color is selected automatically.
  matplotlibcpp::plot(out.x_data, out.y_data, "darkviolet");

  // Plot a line whose name will show up as "log(x)" in the legend.
  // matplotlibcpp::named_plot("log(x)", x, z);

  // Set x-axis to interval [0,1.5] to match Fig. 6.2 of Glendenning (page 286)
  // matplotlibcpp::xlim(0., 1.5);

  // Add graph title
  // matplotlibcpp::title("Rotating Neutron Star Universal Function for an $~8$ ms Pulsar\n$\u03A9 = 7.95\u00D710^{2} s^{-1}$, R = 10.24 km, M = 2.29 $M_\u2609$, $\u03B5_c$ = $2.23\u00D7 10^{-3}$ [$g/cm^3$]\nCompare with Fig. 6.2 of Glendenning");
  matplotlibcpp::title("Frame Dragging Frequency for an $~8$ ms Pulsar\n$\u03A9 = 7.95\u00D710^{2} s^{-1}$, R = 10.24 km, M = 2.29 $M_\u2609$, $\u03B5_c$ = $2.23\u00D7 10^{-3}$ [$g/cm^3$]");

  out.x_label = "r [km]" ;
  out.y_label = "\u03C9(r) $[s^{-1}]$" ;

  matplotlibcpp::xlabel(out.x_label) ;
  matplotlibcpp::ylabel(out.y_label) ;

  // Enable legend.
  // matplotlibcpp::legend();

  // save figure
  std::string filename = "../main/results/RNStar/Glendenning_Table_5-8_29/Glen_5-8_29_Omega_1.png";
  Z_LOG_INFO("Saving result to " + filename) ;
  matplotlibcpp::save(filename);
#endif
  return 0;
}
// *******************************************************