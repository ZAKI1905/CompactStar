// -*- lsst-c++ -*-
/*
* CompactStar
* See License file at the top of the source tree.
*
* Copyright (c) 2024 Mohammadreza Zakeri
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
 * @file LightDM_Scalar_Density.hpp
 *
 * @brief Class for evaluating baryonic background in neutron stars
 *
 * @ingroup LightDM
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// ---------------------------------------------------- 
//               Created on Mar 19, 2024
// This class is made for analyzing ...
// ---------------------------------------------------- 
#ifndef CompactStar_LightDM_Scalar_Density_H
#define CompactStar_LightDM_Scalar_Density_H

#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/Pulsar.hpp"
#include <CompactStar/EOS/CompOSE_EOS.hpp>


//==============================================================
namespace CompactStar
{

//==============================================================
class LightDM_Scalar_Density : public Prog
{
  private:

    /// Flag to indicate if the pulsar is set.
    bool pulsar_is_set = false ; 

  public:
  //--------------------------------------------------------------
   struct Baryon 
    {
      /// Full name
      std::string name  = "" ;

      /// Short name
      std::string short_name = "" ;

      /// Compose label
      std::string label = "" ; 
      
      /// LaTeX name
      std::string TeX_name = "" ;

      /// Mass in vacuum
      double mass = 0 ;

      /// The magnetic moment
      double g = 0 ;

      /// Vacuum life-time (in seconds)
      double tau = 0 ;   

      Baryon() {}
      Baryon(const std::string& in_name, 
              const std::string& in_TeX_name,
              const std::string& in_short_name,
              const std::string& in_label, 
              const double& in_mass,
              const double& in_g,
              const double& in_tau) 
              : name(in_name), 
                TeX_name(in_TeX_name),
                short_name(in_short_name),
                label(in_label),
                mass(in_mass),
                g(in_g), tau(in_tau)
              {}
    } ;

    //--------------------------------------------------------------
  // protected:
  //   //--------------------------------------------------------------
  //   struct Rate_Eps
  //   {
  //     double rate ;
  //     double eps ;
  //   };
  //   //--------------------------------------------------------------
  //   struct Process
  //   {
  //     const std::string name ;
  //     const std::string TeX ;

  //     Process(const std::string& name, const std::string& TeX)
  //     : name(name), TeX(TeX) {}

  //     std::string GetName() const
  //     {
  //       return name ;
  //     }

  //     std::string GetTeX() const
  //     {
  //       return TeX ;
  //     }
  //   };
  //   //--------------------------------------------------------------
  //   virtual Process GetSpecificProcess(const Baryon& B) const
  //   {
  //     return process ;
  //   }
    //--------------------------------------------------------------

    std::string model = "" ;

    // // double PhaseSpace_Integrand(const double& x, const double& mu_chi, 
    // //                             const double& sigma_0) ;
    // double phase_int_mu_chi = 0 ; 
    // double phase_int_sigma_0 = 0 ;

    /// Masses must be converted from MeV to fm^-1 ?
    Zaki::Vector::DataSet m_B_ds ;

    /// Masses must be converted from MeV to fm^-1 ?
    Zaki::Vector::DataSet sig0_ds ;

    /// micro parameters as a function of radius
    /// This is set via "FindPulsar"
    /// [0]: r[km],  [1]: n_n, [2]: n_lam
    /// [3]: m_n, [4]: m_lam, [5]: sig_n, [6]: sig_lam
    Zaki::Vector::DataSet micro_r ;

    /// The EoS
    CompactStar::CompOSE_EOS eos ;

    /// The individual baryon density
    /// [0] : "n_tot"
    /// [1]: "10"
    /// [2] : "100"
    Zaki::Vector::DataSet n_B ;
    
    Pulsar pulsar ;
    
    /// The 2sigma bound on dot{B}/B from pulsar binary
    /// orbital period decay rate
    // double gamma_bnv_bin_lim = 7.83129*1e-10 ;
    // double gamma_bnv_bin_lim = 1.75e-10 ;

    /// Default value for the mixing parameter in MeV
    // double default_eps = 1e-16 ;

    /// The resolution of plot vs m_chi
    // int m_chi_res = 2*1380 ;

    /// The m_chi range
    // Zaki::Math::Axis m_chi_vals ;
  
    /// Plots the rate as a function of radius (ds) in units of s^-1/fm^3
    // void hidden_Plot_Rate_vs_R(const Baryon& B, 
    //                           const double& m_chi,
    //                           Zaki::Vector::DataSet* ds) const ;

    // /// Plots the rate as a function of density (ds) in units of s^-1/fm^3
    // void hidden_Plot_Rate_vs_Density(const Baryon& B, 
    //                                 const double& m_chi,
    //                                 Zaki::Vector::DataSet* ds) const ;

    // /// Plots rates for a set of chi's in the same figure 
    // /// as a function of radius in units of s^-1/fm^3
    // /// for multiple pulsars!
    // /// [ Incomplete ]
    // void Rate_vs_R(const std::vector<double>& m_chi, 
    //                const Baryon& B, 
    //                const std::vector<CompactStar::Pulsar_Limit>& ) ;

    // /// Plots rates for a set of chi's in the same figure 
    // /// as a function of radius in units of s^-1/fm^3
    // /// for multiple pulsars!
    // /// [ Incomplete ]
    // void Rate_vs_R(const std::vector<double>& m_chi, 
    //     const std::vector<CompactStar::Pulsar_Limit>& ) ;
  public:
    Baryon neutron ;
    Baryon lambda  ;
    Baryon proton ;
    Baryon sigma_m ;
    
    // Process process ;

    LightDM_Scalar_Density();
    ~LightDM_Scalar_Density();

    double Scalar_Density_Integrand(const double& k) ;


    /// Sets the EoS model
    void SetModel(const std::string& in_eos_model) ;

    /// Sets the pulsar
    void SetPulsar(const CompactStar::Pulsar&) ;

    // / Imports the EOS
    void ImportEOS(const Zaki::String::Directory& eos_dir) ;

    double Eval_Scalar_Density_vs_Baryon_Density(const double& in_baryon_density, const Baryon& B) ;

    void Export_Scalar_Density_vs_Baryon_Density(const Baryon& B) ;

    void Export_Scalar_Density_vs_Radius(const Baryon& B) ;

    /// Generates the TOV sequence
    void GenSequence() const ;

    /// Finds the pulsar profile
    void FindPulsar(const bool& gen_plots=false) ;

    /// Exports a table with numerical values needed for calculating the escape
    ///  conditions of a particle produced from the decays of a baryon B
    void Export_Escape_Params(const Baryon& B) ;    
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_LightDM_Scalar_Density_H*/
