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
 * @file Particle.hpp
 *
 * @brief Particle base class
 *
 * @ingroup EOS
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit BEFORE Aug 6
#ifndef CompactStar_Particle_H
#define CompactStar_Particle_H

#include <mutex>

#include <Zaki/Math/Math_Core.hpp>

#include "CompactStar/Core/Prog.hpp"
// #include "CompactStar/EOS/SigmaOmegaPar.hpp"
//==============================================================
namespace CompactStar
{

//==============================================================
class Particle : public Prog
{
  friend class SigmaOmegaRho ;
  friend class SigmaOmegaRho_npemu ;
  friend class SigmaOmegaRho_nstar ;
  friend class Fermi_Gas_Many ;
  //--------------------------------------------------------------
  protected:
    /// Chemical potential
    double mu; // Added to store the chemical potential

    /// number density
    double Rho ;

    /// List of pointers to all created baryons
    Zaki::Math::Range<double> inclusion = { 0, 0 } ;

    /// List of pointers to all created particles
    static inline std::vector<Particle*> P_List ;

    /// Lock in case of multi-threading
    static inline std::mutex mtx ;

    /// Unique particle index
    int idx ;

    /// The overall parameters that are fixed by the nuclear properties
    // static inline SigmaOmegaPar pars ;

    /// Flag to keep track of number_densities & whether they are all set or no
    bool number_density_flag = false ;


  //--------------------------------------------------------------
  public:

    /// Electric charge
    const double Q ;

    /// Mass
    const double M ;

    /// Third component of isospin
    const double I3 ;

    Particle( const int& idx, const double& charge, 
              const double& mass,
              const double& isospin ) ;

    ~Particle() ;

    /// Copy constructor
    Particle(const Particle& ) = delete ;

    /// Assignment operator
    Particle& operator=(const Particle& ) = delete ;

    /// Derivative of the chemical potential
    // is this correct??
    virtual double Dmu(const int& idx) ;

    /// Derivative of the chemical potential wrt Rho
    double DMu_DRho() ;

    /// This sets the number density
    void SetRho(const double& rho) ;

    /// Fermi momentum
    double k() const ;

    /// Fermi momentum squared
    double k2() const ;

    // /// Chemical potential
    // virtual double Mu() ;

    /// Contribution to the energy
    virtual double ETerm() const ;

    /// Contribution to the pressure
    virtual double PTerm() const ;

    /// Prints particle info
    void Print() const override ;

    /// Returns the total baryon density
    // static double GetBaryonRho() ;

    /// Returns a pointer to the particle with index = idx
    /// returns null if idx is not availble
    static CompactStar::Particle* Find(const int& idx) ; 

    /// Adds the particle to the pool of particles under consideration
    /// Useful for dynamically changing the pool without needing
    /// heap allocation of particles. 
    /// Returns -1 if not successful and 0 if successful
    virtual int PoolAdd() = 0 ;

    /// Removes the particle from the pool of particles under consideration
    /// Returns -1 if not successful and 0 if successful
    virtual int PoolRemove() = 0 ;

    /// Empties the total particle pool (leptons & baryons)
    static void EmptyPool() ;

        /// Set the chemical potential
    void Set_Mu(const double& mu_in);

    /// Get the chemical potential
    virtual double Mu();

    /// Update number density based on chemical potential
    void Set_Rho_From_Mu();

    /// Get number density
    double GetRho() const;

    // /// Get particle name
    // std::string Name() const;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Particle_H*/
