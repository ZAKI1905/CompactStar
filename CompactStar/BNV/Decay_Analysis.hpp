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
 * @file Decay_Analysis.hpp
 *
 * @brief Analyzes BNV decays in NS [old and obsolete]
 *
 * @ingroup BNV
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Created on June 19, 2022
#ifndef CompactStar_Decay_Analysis_H
#define CompactStar_Decay_Analysis_H

#include <Zaki/Vector/DataSet.hpp>
#include "CompactStar/Core/Analysis.hpp"


//==============================================================
namespace CompactStar
{

//==============================================================

class Pulsar ;

//==============================================================
class Decay_Analysis : public Analysis
{
  //--------------------------------------------------------------
  private:

  struct Baryon 
  {
    std::string name  = "" ;
    std::string label = "" ; 
    double mass = 0 ;
    double g = 0 ;

    Baryon() {}
    Baryon(const std::string& in_name, 
            const std::string& in_label, 
            const double& in_mass,
            const double& in_g) 
            : name(in_name), 
              label(in_label),
              mass(in_mass),
              g(in_g)
            {}
  } ;

  Baryon neutron ;
  Baryon lambda  ;

  //--------------------------------------------------------------
  // double g_lam = -1.22 ;

  /// Mixing parameter in fm^{-1}
  double eps_lam = 1 ;

  /// Chi's mass in fm^{-1} 
  double m_chi ;

  double AmplitudeSqrd(const double& m_lam_bar) const ;
  
  double PhaseSpace(const double& m_lam_bar, 
                    const double& n_lam) const ;

  /// The Lambda decay rate per volume is in units of 
  ///     fm^-3 s^-1 
  /// m_lm is in units of fm^{-1}
  /// n_lam is in units of fm^{-3}
  double n_dot( const Baryon& B, const double& m_lam, 
                const double& V_lam, const double& n_lam) const ;

  double G_func(const double& x) const ;

  /// Outputs the decay rate in [s^-1]
  ///   eps_lam input is in fm^-1
  double Baryon_Vacuum_Decay( const Baryon& B, 
                              const double& in_eps_lam) const ;

  // Masses must be converted from MeV to fm^-1
  Zaki::Vector::DataSet m_eff_ds_2 ;

  // Masses must be converted from MeV to fm^-1
  Zaki::Vector::DataSet m_eff_ds ;

  // Masses must be converted from MeV to fm^-1
  Zaki::Vector::DataSet V_self_E_ds ;

  Pulsar* pulsar ;

  //--------------------------------------------------------------
  public:

    Decay_Analysis() ;

    ~Decay_Analysis() ;

    /// Analysis during the sequence loop
    void Analyze(NStar* in_star) override ;
    void Analyze(MixedStar* in_star) override {}

    /// Saves the results
    void Export(const Zaki::String::Directory&) override ;

    /// Imports the effective mass set
    void ImportEffMass(const std::string& f_name, 
                       const Zaki::String::Directory& dir="") ;

    // Imports the vector self-energy set
    void ImportVSelfEnergy(const std::string& f_name, 
                          const Zaki::String::Directory& in_dir=""); 

    /// Plots [ V_self_E_ds + m_eff_ds ]
    void PlotRestEnergy() ;

    /// Attaches the pulsar
    void AttachPulsar(Pulsar* puls) ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_Decay_Analysis_H*/
