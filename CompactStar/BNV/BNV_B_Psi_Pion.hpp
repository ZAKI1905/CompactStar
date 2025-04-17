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
 * @file BNV_B_Psi_Pion.hpp
 *
 * @brief Evaluates B -> psi pion decay rates in neutron stars
 *
 * @ingroup BNV
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// ---------------------------------------------------- 
//               Created on Oct 4, 2022
//        Updated on Mar 15, 2023
// This class is made specifically for analyzing the
//  B --> Psi Pion
//  decays in neutron stars. 
// ---------------------------------------------------- 
#ifndef CompactStar_BNV_B_Psi_Pion_H
#define CompactStar_BNV_B_Psi_Pion_H

#include <Zaki/Vector/DataSet.hpp>

#include <CompactStar/BNV/BNV_Chi.hpp>
#include <CompactStar/EOS/CompOSE_EOS.hpp>
#include "CompactStar/Core/Pulsar.hpp"


//==============================================================
namespace CompactStar
{

//==============================================================
class BNV_B_Psi_Pion : public BNV_Chi
{
  //--------------------------------------------------------------
  private:
    int gen_i, gen_j ;
    // struct Baryon 
    // {
    //   /// Full name
    //   std::string name  = "" ;

    //   /// Short name
    //   std::string short_name = "" ;

    //   /// Compose label
    //   std::string label = "" ; 
      
    //   /// LaTeX name
    //   std::string TeX_name = "" ;

    //   /// Mass in vacuum
    //   double mass = 0 ;

    //   /// The magnetic moment
    //   double g = 0 ;

    //   /// Vacuum life-time (in seconds)
    //   double tau = 0 ;   

    //   Baryon() {}
    //   Baryon(const std::string& in_name, 
    //           const std::string& in_TeX_name,
    //           const std::string& in_short_name,
    //           const std::string& in_label, 
    //           const double& in_mass,
    //           const double& in_g,
    //           const double& in_tau) 
    //           : name(in_name), 
    //             TeX_name(in_TeX_name),
    //             short_name(in_short_name),
    //             label(in_label),
    //             mass(in_mass),
    //             g(in_g), tau(in_tau)
    //           {}
    // } ;

    // ----------------------------------------------
    //                  Constants
    // ----------------------------------------------

    //------------------
    double GF = 1.1663788e-11 ; // MeV^-2
    double mW = 80.377e3 ; // MeV
    double m_X = 1e+6 ; // MeV
    double m_pion = 140 ; // MeV

    /// From lattice QCD
    double alpha = -0.014e9 ; // MeV^3

    /// Non-perturbative baryon-meson couplings
    double D = 0.80 ; 
    double F = 0.46 ;

    /// Pion decay constant (in MeV)
    double f_pi = 92.4  ;

    /// Pole term for neutron to pion (via neutron propagator)
    double b_B = - (D + F) / 2.0 ;

    /// Contact term for neutron to pion
    double c_B = 0.5 ;

    //------------------
    //   Quarks
    //------------------
    double mu = 2.16 ; // MeV 
    double mc = 1.27e3 ; // MeV
    double mt = 172.69e3 ; // MeV
    double md = 4.67 ; // MeV
    double ms = 93.4 ; // MeV
    double mb = 4.18e3 ; // MeV
    std::vector<double> m_D = {md, ms, mb} ;
    std::vector<double> m_U = {mu, mc, mt} ;
    //------------------
    //   CKM
    //------------------
    // 'ud'    'us'    'ub'
    // 'cd'    'cs'    'cb'
    // 'td'    'ts'    'tb'
    std::vector< std::vector<double> > V_CKM 
          = {   {0.97373, 0.2243, 3.82e-3},
                {0.221, 0.975, 40.8e-3},
                {8.6e-3, 41.5e-3, 1.014} } ;
    //------------------

    // ----------------------------------------------
    // Parameters needed for evaluating the phase-space integrals
    struct Params
    {
      // Effective mass of the baryons (m_B^*)
      double m_B ;

      // x = E_B^* (nm-frame) / m_B^*
      double x ;

      // mu_pi = m_pi / m_B^*
      double mu_pi ;

      // mu_psi = m_psi / m_B^*
      double mu_psi ;

      // sig = Sigma^0 (nm-frame) / m_B^*
      double sig ;

      void SetPars(const double& in_m_B,
                   const double& in_mu_pi,
                   const double& in_mu_psi,
                   const double& in_sig )
      {
        m_B = in_m_B ;
        mu_pi = in_mu_pi ;
        mu_psi = in_mu_psi ;
        sig = in_sig ;
      }
    };
    //--------------------------------------------------------------

    Params selected_pars ;
    // double g_lam = -1.22 ;

    /// Mixing parameter in fm^{-1}
    // double eps_lam = 1 ;

    // std::string model = "" ;

    // double PhaseSpace_Integrand(const double& x, const double& mu_chi, 
    //                             const double& sigma_0) ;
    // double phase_int_mu_chi = 0 ; 
    // double phase_int_sigma_0 = 0 ;

    /// Outputs the decay rate in [s^-1]
    ///   eps_lam input is in fm^-1
    // double Vacuum_Decay_Br( const Baryon& B, 
                            // const double& in_m_chi,
                            // const double& in_eps_lam) const ;

    /// Masses must be converted from MeV to fm^-1
    // Zaki::Vector::DataSet m_B_ds ;

    /// Masses must be converted from MeV to fm^-1
    // Zaki::Vector::DataSet sig0_ds ;

    // /// The EoS
    // CompactStar::CompOSE_EOS eos ;

    /// The individual baryon density
    /// [0] : "n_tot"
    /// [1]: "10"
    /// [2] : "100"
    // Zaki::Vector::DataSet n_B ;
    
    // Pulsar pulsar ;
    // bool cgs_units = false ;

    // /// The resolution of plot vs m_chi
    // int m_chi_res = 750 ;

    // Process GetSpecificProcess(const Baryon& B) const override ;
  //--------------------------------------------------------------
  public:
    // Baryon neutron ;
    // Baryon lambda  ;
    // Baryon proton ;

    // Input the generation index in the loop (i, j)
    BNV_B_Psi_Pion(const int& i, const int& j);

    ~BNV_B_Psi_Pion();
  
    /// Returns the time component of Sigma^0 in the baryon CM frame
    double get_Sig_0_CM(const BNV_B_Psi_Pion::Params& pars) ;

    /// Returns the length of 3-vector vec{Sigma} in the baryon CM frame
    double get_Sig_V_CM(const BNV_B_Psi_Pion::Params& pars) ;

    /// Returns the length of 3-vector vec{q} in the baryon CM frame
    ///    q is the momentum of psi
    double get_q_V_CM(const BNV_B_Psi_Pion::Params& pars) ;

    /// Returns the energy of psi in the baryon CM frame
    double get_E_Psi_CM(const BNV_B_Psi_Pion::Params& pars) ;

    /// Returns the energy* of baryon (E_B^*) in the nm frame
    double get_E_B_s(const BNV_B_Psi_Pion::Params& pars) ;

    /// Returns the energy of baryon (E_B) in the CM frame
    double get_E_B_CM(const BNV_B_Psi_Pion::Params& pars) ;

    /// Returns the energy* of baryon (E_B^*) in the CM frame
    double get_E_B_s_CM(const BNV_B_Psi_Pion::Params& pars) ;
    
    // ------------------------------------------
    //  Anglar momentum dot prdocuts
    // ------------------------------------------
    // ------------------------
    ///   The angular indep part of the 
    ///     4-vector dot-product of (q.p*) 
    double q_pstar_zero(const BNV_B_Psi_Pion::Params& pars) ;
    
    ///   The angular dep part (without cos(theta) ) of the 
    ///     4-vector dot-product of (q.p*) 
    double q_pstar_first(const BNV_B_Psi_Pion::Params& pars) ;
    // ------------------------

    // ------------------------
    ///   The angular indep part of the 
    ///     4-vector dot-product of (q*.p*) 
    double qstar_pstar_zero(const BNV_B_Psi_Pion::Params& pars) ;

    ///   The angular indep part of the 
    ///     4-vector dot-product of (q*.p*) 
    double qstar_pstar_first(const BNV_B_Psi_Pion::Params& pars) ;
    // ------------------------

    // ------------------------
    ///   The angular indep part of the 
    ///     4-vector dot-product of (q.q*) 
    double q_qstar_zero(const BNV_B_Psi_Pion::Params& pars) ;

    ///   The angular indep part of the 
    ///     4-vector dot-product of (q.q*) 
    double q_qstar_first(const BNV_B_Psi_Pion::Params& pars) ;
    // ------------------------

    // ------------------------
    ///   The angular indep part of the 
    ///     4-vector dot-product of (q*.q*) 
    double qstar_qstar_zero(const BNV_B_Psi_Pion::Params& pars) ;

    ///   The angular indep part of the 
    ///     4-vector dot-product of (q*.q*) 
    double qstar_qstar_first(const BNV_B_Psi_Pion::Params& pars) ;
    // ------------------------
    // ---------------------------------

    /// The W0 term with no extra cos(theta) factors
    ///   integrated over cos(theta)
    double phase_W0_zero(const BNV_B_Psi_Pion::Params& pars) ;

    /// The W0 term with an extra cos(theta) factor
    ///   integrated over cos(theta)
    double phase_W0_first(const BNV_B_Psi_Pion::Params& pars) ;

    /// The W1 term with no extra cos(theta) factors
    ///   integrated over cos(theta)
    double phase_W1_zero(const BNV_B_Psi_Pion::Params& pars) ;

    /// The W1 term with an extra cos(theta) factor
    ///   integrated over cos(theta)
    double phase_W1_first(const BNV_B_Psi_Pion::Params& pars) ;

    /// The W1 term with two extra cos(theta) factor
    ///   integrated over cos(theta)
    double phase_W1_second(const BNV_B_Psi_Pion::Params& pars) ;

    /// The W0*W1 term with no extra cos(theta) factors
    ///   integrated over cos(theta)
    double phase_W0W1_zero(const BNV_B_Psi_Pion::Params& pars) ;

    /// The W0*W1 term with one extra cos(theta) factors
    ///   integrated over cos(theta)
    double phase_W0W1_first(const BNV_B_Psi_Pion::Params& pars) ;

    /// The A-term in the amplitude sqrd
    double Amp_A(const BNV_B_Psi_Pion::Params& pars) ;

    ///  Loop integral in dimension-less form
    double LoopIntegral(const double& xD,const double& xU, const double& xY) ;
    
    /// The averaged amplitude sqrd (integrated over Cos[theta])
    double Amplitude(const BNV_B_Psi_Pion::Params& pars) ;

    /// The baryon decay rate in the CM frame
    double CM_Decay_Rate(const BNV_B_Psi_Pion::Params& pars) ;

    /// The baryon decay rate in the nm frame
    double NM_Decay_Rate(const BNV_B_Psi_Pion::Params& pars) ;

    /// Baryon density loss rate (the output is positive)
    double PhaseSpace_Integral(const double& x_F) ;

    /// Analysis during the sequence loop
    // void Analyze(NStar* in_star) override;
    // void Analyze(MixedStar* in_star) override {}

    /// Generates the TOV sequence
    // void GenSequence() const ;
  
    double PhaseSpace_Integrand(const double& x) override ;
    // double PhaseSpace_Integral( const Baryon& B,
    //                             const double& m, 
    //                             const double& m_chi, 
    //                             const double& Sigma_0, 
    //                             const double& x_F) override ;

    /// Returns rate per unit volume
    ///  as a function of density in units of s^-1/fm^3
    Zaki::Vector::DataSet Rate_vs_Density(
                              const double& m_psi, 
                              const Baryon& B, 
                              const bool& gen_plots=false) override ;

    /// Returns rate per unit volume
    ///  as a function of radius in units of s^-1/fm^3
    Zaki::Vector::DataSet Rate_vs_R(const double& m_psi, 
                                    const Baryon& B, 
                                    const bool& gen_plots=false) 
                                    override ;
    // /// Returns rate per unit volume
    // ///  as a function of radius in units of s^-1/fm^3
    // Zaki::Vector::DataSet Rate_vs_R_Slow(const double& m_chi, 
    //                                 const Baryon& B, 
    //                                 const bool& gen_plots=false) 
    //                                 ;                                
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_BNV_B_Psi_Pion_H*/
