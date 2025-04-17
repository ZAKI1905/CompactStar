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
 * @file Lepton.hpp
 *
 * @brief Lepton class
 *
 * @ingroup EOS
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit BEFORE Aug 6
#ifndef CompactStar_Lepton_H
#define CompactStar_Lepton_H

#include <mutex>

#include "CompactStar/EOS/Particle.hpp"
// #include "CompactStar/SigmaOmegaPar.hpp"
//==============================================================
namespace CompactStar
{

//==============================================================
class Lepton : public Particle
{
  friend class SigmaOmegaRho ;
  friend class SigmaOmegaRho_npemu ;
  friend class SigmaOmegaRho_nstar ;
  //--------------------------------------------------------------
  private:

    /// L-number density
    // double RhoL ;

    /// List of pointers to all created Leptons
    static inline std::vector<Lepton*> L_List ;


  //--------------------------------------------------------------
  public:

    Lepton( const int& idx, const double& charge, const double& mass ) ;

    ~Lepton() ;

    /// Copy constructor
    Lepton(const Lepton& ) = delete ;

    /// Assignment operator
    Lepton& operator=(const Lepton& ) = delete ;

    /// Adds the lepton to the pool of leptons under consideration
    /// Useful for dynamically changing the pool without needing
    /// heap allocation of particles. 
    /// Returns -1 if not successful and 0 if successful
    int PoolAdd() override ;

    /// Removes the lepton from the pool of leptons under consideration
    /// Returns -1 if not successful and 0 if successful
    int PoolRemove() override ;

    /// Shows the number of leptons
    static void PrintAll() ;

    /// Empties the lepton's pool
    static void EmptyLPool() ;

};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Lepton_H*/
