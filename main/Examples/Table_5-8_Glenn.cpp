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

// Last edited on Aug 25, 2021
#include <CompactStar/EOS/Common.hpp>
#include <CompactStar/Core/TOVSolver.hpp>

#include <matplotlibcpp.hpp>

// *******************************************************
// We use the EOS from Table 5.8 on page 274 of Ch. 5
// of Compact Stars book by Glendenning
// K = 240 MeV  &  m*/m = 0.78  &  x_sigma = 0.6
// *******************************************************
int main()
{
  Zaki::String::Directory dir(__FILE__) ;

#if 1
  // -----------------------------------------------------------------------------
  //      TOV Solver
   CompactStar::TOVSolver solver ;

   solver.SetWrkDir(dir.ParentDir() +"/results") ;
   solver.ImportEOS("EOS/Glendenning_Table_5_8.eos") ;
  
   solver.Solve({{ 1e+06,
                   3e+15},
                 200, "Log"}, "NStar/Glendenning_Extended", 
                 "Glendenning_Table_5-8") ;

   solver.ExportSequence("NStar/Glendenning_Extended/Glendenning_Table_5-8_sequence.tsv") ;

  // -----------------------------------------------------------------------------
#else
  // -----------------------------------------------------------------------------
  CompactStar::PlotData out = CompactStar::ColExtractor(dir.ParentDir() +
     "/results/NStar/Glendenning/Glendenning_Table_5-8_sequence.tsv", 1, 2) ;

  // Set the size of output image = 1200x780 pixels
  matplotlibcpp::figure_size(1200, 780);

  // Plot line from given x and y data. Color is selected automatically.
  matplotlibcpp::semilogx(out.x_data, out.y_data, "b");

  // Add graph title
  matplotlibcpp::title("Glendenning Table 5.8 Sequence");

  matplotlibcpp::xlabel(out.x_label) ;
  matplotlibcpp::ylabel(out.y_label) ;

  // Enable legend.
  // matplotlibcpp::legend();

  // save figure
  std::string filename = "../main/results/NStar/Glendenning/Glendenning_Table_5-8_sequence.png";
  Z_LOG_INFO("Saving result to " + filename) ;
  matplotlibcpp::save(filename);
  // -----------------------------------------------------------------------------
#endif
  return 0;
}
// *******************************************************
