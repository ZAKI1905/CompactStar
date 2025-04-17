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
 * @file Baryon.hpp
 *
 * @brief Baryon class.
 *
 * @ingroup EOS
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit BEFORE Aug 6
#ifndef CompactStar_Baryon_H
#define CompactStar_Baryon_H

#include <mutex>

#include "CompactStar/EOS/Particle.hpp"
#include "CompactStar/EOS/SigmaOmegaPar.hpp"
//==============================================================
namespace CompactStar
{

//==============================================================
class Baryon : public Particle
{
  friend class SigmaOmegaRho ;
  friend class SigmaOmegaRho_npemu ;
  friend class SigmaOmegaRho_nstar ;
  friend class Fermi_Gas_Many ;
  //--------------------------------------------------------------
  private:

    /// B-number density
    // double RhoB ;

    /// List of pointers to all created Baryons
    static inline std::vector<Baryon*> B_List ;

    /// Lock in case of multi-threading
    // static inline std::mutex mtx ;

    /// The overall parameters that are fixed by the nuclear properties
    static inline SigmaOmegaPar pars ;

    /// Flag to keep track of rhos & whether they are all set or no
    // bool rho_flag = false ;

    /// x_omega = g_omega(baryon) / g_omega(proton) 
    const double x_ome = 1 ;

    /// x_rho = g_rho(baryon) / g_rho(proton)
    const double x_rho = 1 ;

    /// x_sig = g_sig(baryon) / g_sig(proton)
    const double x_sig = 1 ;


    static double D_I_1_a(const double& k, const double& a) ; 

    Baryon( const int& idx, const double& charge, const double& mass, 
            const double& isospin) ;

    Baryon( const int& idx, const double& charge, const double& mass, 
            const double& isospin, const double& in_x_ome,
            const double& in_x_rho, const double& in_x_sig) ;

    ~Baryon() ;
  //--------------------------------------------------------------
  public:

    /// Copy constructor
    Baryon(const Baryon& ) = delete ;

    /// Assignment operator
    Baryon& operator=(const Baryon& ) = delete ;

    /// Setting the overall parameters of the model
    static void SetPars(const SigmaOmegaPar& in_pars) ;

    /// This sets the B-number density
    // void SetRho(const double& rho_B) ;

    /// Returns the value of [ g_omega(proton) * omega_0 ]
    static double gpOmega() ;

    /// Returns the value of [ g_rho(proton) * rho_03 ]
    static double gpRho03() ;

    /// Returns the total baryon density
    static double GetBaryonRho() ;

    /// Chemical potential
    double Mu(const double& gsig) const ;

    /// Equation describing g_sigma * sigma
    static double gsigEq(const double& gsig) ;

    /// Contribution to the energy
    double ETerm(const double& gsig) const ;

    /// Contribution to the pressure
    double PTerm(const double& gsig) const ;

    /// Shows the number of baryons
    static void PrintAll() ;

    /// Derivative of the chemical potential for baryons
    double Dmu(const int& in_idx, const double& gsig) ;

    /// Derivative of the g_sig * sigma equation
    static double Dgsig(const int& in_idx, const double& gsig) ;


    /// Adds the baryon to the pool of baryons under consideration
    /// Useful for dynamically changing the pool without needing
    /// heap allocation of particles. 
    /// Returns -1 if not successful and 0 if successful
    int PoolAdd() override ;

    /// Removes the baryon from the pool of baryons under consideration
    /// Returns -1 if not successful and 0 if successful
    int PoolRemove() override ;

    /// Empties the baryon's pool
    static void EmptyBPool() ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Baryon_H*/
