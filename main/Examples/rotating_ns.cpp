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
#include <CompactStar/Core/RotationSolver.hpp>
#include <CompactStar/EOS/Common.hpp>


// *******************************************************
int main()
{
	Zaki::String::Directory dir(__FILE__);

#if 1
	// ----------------------------------------------------------------------------------------
	CompactStar::Core::RotationSolver r_solver;

	r_solver.SetWrkDir(dir.ParentDir() + "/results");
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
	r_solver.Solve({{5.66379367363 * 1e-4, 5.66379367363 * 1e-3},
					3,
					"Log"},
				   "RNStar/Glendenning_Table_5-8_29/Glen_5-8_29_Omega_*.tsv");

	r_solver.ExportResults("RNStar/Glendenning_Table_5-8_29/Glen_5-8_29_Omega_sequence.tsv");

	// r_solver.Solve({{ 0.01, 0.1},
	//               30, "Log"}, "RNStar/SigOmegRho_240_078_Hyp_29/SOR_240_78_Hy_29_Omega_*.tsv") ;

	// r_solver.ExportResults("RNStar/SigOmegRho_240_078_Hyp_29/SOR_240_78_Hyp_29_Omega_sequence.tsv") ;
	// ----------------------------------------------------------------------------------------
#endif
	return 0;
}
// *******************************************************