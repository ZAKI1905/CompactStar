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
 * @file BNV_B_Chi.hpp
 *
 * @brief Abstract class for evaluating baryon decay rates in neutron stars
 *
 * @ingroup BNV
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// ---------------------------------------------------- 
//               Created on Mar 15, 2023
// This abstract class is made for analyzing BNV
// processes that produce Chi from mixing, e.g,
//  B --> Chi, B --> Chi + photon, 
//  B + P --> Chi + P, etc.
// ---------------------------------------------------- 
#ifndef CompactStar_BNV_Chi_H
#define CompactStar_BNV_Chi_H

#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/Pulsar.hpp"
#include <CompactStar/EOS/CompOSE_EOS.hpp>


//==============================================================
namespace CompactStar
{

struct Pulsar_Limit
{
  CompactStar::Pulsar p ;
  std::vector<double> lim ;

  Pulsar_Limit(const CompactStar::Pulsar& in_p,
               const std::vector<double>& in_lim )
               : p(in_p), lim(in_lim)
  {}

  // Pulsar_Limit(const CompactStar::Pulsar& in_p,
  //              std::vector<double> in_lim )
  //              : p(in_p), lim(in_lim)
  // {}
};

//==============================================================
class BNV_Chi : public Prog
{
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
  protected:
    //--------------------------------------------------------------
    struct Rate_Eps
    {
      double rate ;
      double eps ;
    };
    //--------------------------------------------------------------
    struct Process
    {
      const std::string name ;
      const std::string TeX ;

      Process(const std::string& name, const std::string& TeX)
      : name(name), TeX(TeX) {}

      std::string GetName() const
      {
        return name ;
      }

      std::string GetTeX() const
      {
        return TeX ;
      }
    };
    //--------------------------------------------------------------
    virtual Process GetSpecificProcess(const Baryon& B) const
    {
      return process ;
    }
    //--------------------------------------------------------------

    // double g_lam = -1.22 ;

    /// Mixing parameter in fm^{-1}
    double eps_lam = 1 ;

    std::string model = "" ;

    // double PhaseSpace_Integrand(const double& x, const double& mu_chi, 
    //                             const double& sigma_0) ;
    double phase_int_mu_chi = 0 ; 
    double phase_int_sigma_0 = 0 ;

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
    double gamma_bnv_bin_lim = 1.75e-10 ;

    /// Default value for the mixing parameter in MeV
    double default_eps = 1e-16 ;

    /// The resolution of plot vs m_chi
    // int m_chi_res = 2*1380 ;

    /// The m_chi range
    Zaki::Math::Axis m_chi_vals ;
  
    /// Plots the rate as a function of radius (ds) in units of s^-1/fm^3
    void hidden_Plot_Rate_vs_R(const Baryon& B, 
                              const double& m_chi,
                              Zaki::Vector::DataSet* ds) const ;

    /// Plots the rate as a function of density (ds) in units of s^-1/fm^3
    void hidden_Plot_Rate_vs_Density(const Baryon& B, 
                                    const double& m_chi,
                                    Zaki::Vector::DataSet* ds) const ;

    /// Plots rates for a set of chi's in the same figure 
    /// as a function of radius in units of s^-1/fm^3
    /// for multiple pulsars!
    /// [ Incomplete ]
    void Rate_vs_R(const std::vector<double>& m_chi, 
                   const Baryon& B, 
                   const std::vector<CompactStar::Pulsar_Limit>& ) ;

    /// Plots rates for a set of chi's in the same figure 
    /// as a function of radius in units of s^-1/fm^3
    /// for multiple pulsars!
    /// [ Incomplete ]
    void Rate_vs_R(const std::vector<double>& m_chi, 
        const std::vector<CompactStar::Pulsar_Limit>& ) ;
  public:
    Baryon neutron ;
    Baryon lambda  ;
    Baryon proton ;
    Baryon sigma_m ;
    
    Process process ;

    BNV_Chi(const Process& in_process);
    ~BNV_Chi();

    /// Vacuum rate for B -> chi + photon
    /// Outputs the decay rate in [s^-1]
    ///   eps_lam input is in fm^-1
    double Vacuum_Decay_Br( const Baryon& B, 
                            const double& in_m_chi,
                            const double& in_eps_lam) const ;


    /// Sets the EoS model
    virtual void SetModel(const std::string& in_eos_model) ;

    /// Sets the pulsar
    void SetPulsar(const CompactStar::Pulsar&) ;

    /// Sets m_chi_vals
    void SetChiMassRange(const Zaki::Math::Axis&) ;

    /// Eventually will be replaced when I switch to using 
    /// binary pulsars
    /// sets 'gamma_bnv_bin_lim' in yr^-1
    void SetBNVLimit(const double& in_gamma) ;

    /// Analysis during the sequence loop
    // void Analyze(NStar* in_star) override;
    // void Analyze(MixedStar* in_star) override {}

    /// Imports the EOS
    /// if micro_dir is not input, it will look for it 
    /// in the same directory as eos_dir
    virtual void ImportEOS(const Zaki::String::Directory& eos_dir, 
                           const std::string& micro_dir="") ;

    /// Generates the TOV sequence
    void GenSequence() const ;

    /// Finds the pulsar profile
    virtual void FindPulsar(const bool& gen_plots=false) ;

    /// Plots the effective mass as a function of radius
    void Plot_Meff_Radius() ;

    /// Plots the baryon's rest-energy as a function of radius
    void Plot_RestEnergy_Radius() ;

    /// Plots the Fermi energy as a function of radius
    void Plot_EF_Radius() ;

    /// @brief  Plots both the rest energy and the Fermi energy
    void Plot_RestE_EF_Radius(const Baryon& B) ;

    /// Plots both the CM frame energy as a function of radius
    void Plot_CM_E_Radius(const Baryon& B) ;

    // Plots both E* band in the nm frame as a function of radius
    void Plot_Estar_Radius(const Baryon& B) ;

    virtual double PhaseSpace_Integrand(const double& x) ;
    virtual double PhaseSpace_Integral( const Baryon& B,
                                const double& m, 
                                const double& m_chi, 
                                const double& Sigma_0, 
                                const double& x_F) ;

    /// Finds and plots the vacuum branching ratio limits
    void PlotVacuumBrLim(const Baryon& B) ;
    
    /// Finds and plots the vacuum branching ratio limits
    void PlotVacuumBrLim() ;

    /// Returns rate 
    ///  as a function of density in units of s^-1/fm^3
    virtual Zaki::Vector::DataSet Rate_vs_Density(
                              const double& m_chi, 
                              const Baryon& B, 
                              const bool& gen_plots=false) = 0 ;

    /// Returns rate 
    ///  as a function of radius in units of s^-1/fm^3
    virtual Zaki::Vector::DataSet Rate_vs_R(const double& m_chi, 
                                    const Baryon& B, 
                                    const bool& gen_plots=false) = 0 ;

    // Plots rate 
    // as a function of radius in units of s^-1/fm^3
    void Rate_vs_R(const double& m_chi) ;

    // Plots rates for a set of chi's in the same figure 
    // as a function of radius in units of s^-1/fm^3
    void Rate_vs_R(const std::vector<double>& m_chi) ;

    // Plots rates for a set of chi's in the same figure 
    // as a function of radius in units of s^-1/fm^3
    void Rate_vs_R(const std::vector<double>& m_chi, 
                   const Baryon& B) ;

    // This function returns the (per baryon) rate for
    //  B -> chi decay integrated over the radius 
    //  and in units of [ yr^-1 ]
    double GetRate(const double& m_chi, 
                      const Baryon& B) ;

    /// Returns the limit on eps from 
    ///  B -> chi decay rate integrated over the radius 
    ///  and in units of [ MeV ]
    double GetEpsLim(const double& m_chi, 
                    const Baryon& B) ;       

    /// Returns the rate and limit on eps from 
    ///  B -> chi decay rate integrated over the radius 
    ///  and in units of [ s^-1 ] & MeV
    Rate_Eps GetRate_Eps(const double& m_chi, 
                          const Baryon& B) ;

    /// Plots the total rate and the limit on eps  
    ///  as a function of m_chi
    void PlotRate_Eps(const Baryon& B) ;

    /// Plots the total rate and the limit on eps 
    ///  as a function of m_chi
    void PlotRate_Eps() ; 

    /// Saves the results
    void Export(const Zaki::String::Directory&) const ; 
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_BNV_Chi_H*/
