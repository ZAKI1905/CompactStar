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

/**
 * @file SigmaOmegaRho_npemu.hpp
 *
 * @brief Analyzes BNV chemical heating [incomplete]
 *
 * @ingroup ChemicalHeating
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
#ifndef CompactStar_SigmaOmegaRho_npemu_H
#define CompactStar_SigmaOmegaRho_npemu_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "CompactStar/EOS/Model.hpp"
#include "CompactStar/EOS/SigmaOmegaPar.hpp"
#include "CompactStar/EOS/Baryon.hpp"
#include "CompactStar/EOS/Lepton.hpp"

//==============================================================
namespace CompactStar
{

//==============================================================
class SigmaOmegaRho_npemu : public Model
{
  //--------------------------------------------------------------
  private:
    SigmaOmegaPar params ;

    double rho   = -1 ; // Total Baryon number density

    double gsig  = -1 ; // g_sigma*sigma value

    Baryon proton, neutron;
    // Baryon sigma_minus, lambda, sigma_zero, sigma_plus, xi_minus, xi_zero ;
    Lepton electron ;
    Lepton muon ;

    size_t eq_size = 2 ;

    /// Solution guess for  
    // [ g_sigma * sigma, 
    //   Rho(proton)/RhoB, Rho(sigma-)/RhoB, Rho(muon)/RhoB]
    // double root_guess[4] = {0.2, 0.03, 0.01, 9.927*0.12} ;
    // double root_guess[4] = {0.21*MN, 0.325, 0.002, 0.01} ;
    std::vector<double> root_guess = {0.19*MN, 0.2, 9e-2} ;
    // double root_guess[4] = {0.086184, 0.0127884, 0.00123263, 0.825926} ;

    /// Contains the information for a good guess depending on rho
    void FindGuess(const double& in_rho) ;

    /// Sets the individual density ratio from the input x
    /// for electron & neutron (particles with index greater than eq_size)
    /// it will use the charge & baryon number conservation respectively
    void SetXRhos(const gsl_vector* x) ;

    /// Sets all the equations given a specific density 'in_rho' 
    /// It sets g_sig*sig equation & chemical equilibrium equations
    void SetEquilibriumEqs(const double& in_rho) ;

    /// vector contating the values of the equations
    std::vector<double> equilibrium_eqs ;

  //--------------------------------------------------------------
  public:

    SigmaOmegaRho_npemu() ;

    ~SigmaOmegaRho_npemu() ;

    double EDens(const double& in_rho) override ;
    double Press(const double& in_rho) override ;

    /// This will generate a row of values for EOS
    std::vector<double> EOSRow(const double&) override ;

    /// This will generate the header for EOS
    std::string EOSHeader() const override ;

    void SetPars(const SigmaOmegaPar& in_pars) ;
    
    int Eq(const gsl_vector* x, gsl_vector* f) ;
    int SolveEq() ;
    void PrintState(const size_t& iter, gsl_multiroot_fsolver * s) ;

    /// Jacobian calculator
    int Jacob(const gsl_vector* x, gsl_matrix* J) ;
    int EqFdf(const gsl_vector* x, gsl_vector* f, gsl_matrix* J) ;
    int SolveEqJacob() ;
    void PrintState(const size_t& iter, gsl_multiroot_fdfsolver * s) ;
    void UpdatePool(const double& in_rho) ;
    

    void SetRho(const double& in_rho)
    {
      rho = in_rho ;
    }

};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_SigmaOmegaRho_npemu_H*/
