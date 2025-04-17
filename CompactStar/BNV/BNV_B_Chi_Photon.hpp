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
 * @file BNV_B_Chi_Photon.hpp
 *
 * @brief Evaluates B -> chi gamma decay rates in neutron stars
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
//  B --> Chi Gamma
//  decays in neutron stars. It is very similar in 
//  structure to "Decay_Analysis" class which was made
//  earlier to evaluate the same decay rate.
//  "Decay_Analysis" decay formula in NS is wrong.
//  "BNV_B_Chi_Photon" fixes the mistakes in 
//  "Decay_Analysis". We are keeping the 
//  "Decay_Analysis" files for now, but might delete it
//   in the future.
// ---------------------------------------------------- 
#ifndef CompactStar_BNV_B_Chi_Photon_H
#define CompactStar_BNV_B_Chi_Photon_H

#include <Zaki/Vector/DataSet.hpp>

#include <CompactStar/BNV/BNV_Chi.hpp>
#include <CompactStar/EOS/CompOSE_EOS.hpp>
#include "CompactStar/Core/Pulsar.hpp"


//==============================================================
namespace CompactStar
{

//==============================================================
class BNV_B_Chi_Photon : public BNV_Chi
{
  //--------------------------------------------------------------
  private:
    
    // Pulsar pulsar ;
    bool cgs_units = false ;

    // /// The resolution of plot vs m_chi
    // int m_chi_res = 750 ;

    Process GetSpecificProcess(const Baryon& B) const override ;
  //--------------------------------------------------------------
  public:

    BNV_B_Chi_Photon();
    ~BNV_B_Chi_Photon();

    // -----------------------------------------------
    //            Escape Rates
    // -----------------------------------------------
    /// Returns the cosine of the angle between p_chi^{CM} and p_B^{nm}
    /// Added on Apr 23, 2024
    double Cos_Theta_CMpChi_NMpB(const double& E_chi_NM_over_mB, const double& mu_chi, 
                                  const double& x, const double& sigma_0)  ;

    // Returns the cosine of the angle between p_chi^{CM} and p_B^{nm}
    // Assuming that the escape condition is satisfied
    // This returns the minimum of cos(theta) for such condition.
    // If it is larger than 1, then it doesn't exist.
    // Added on Apr 23, 2024
    // double Cos_Theta_CMpChi_NMpB_Escape(const double& mu_chi, const double& x, 
    //                                     const double& sigma_0)  ;

    // Exp[-2nu(r)]
    // Used in ...
    double phase_int_escape_exp_2nu = -1 ;

    // Phase space integrand function (dim-less part) for 
    //    escaping chi particles
    // The total phase space integrand is the output of this function
    // multiplied by:
    //                  m^4 C / (pi^2)
    //
    double Escape_PhaseSpace_Integrand(const double& x) ;

    // The full phase-space integral for escaping particles
    // The output is in MeV^4
    double Escape_PhaseSpace_Integral(const Baryon& B,
                                              const double& m, 
                                              const double& m_chi, 
                                              const double& Sigma_0, 
                                              const double& x_F) ;

    // Evaluates the escape rate as a function of radius                                          
    Zaki::Vector::DataSet Escape_Rate_vs_R(const Baryon& B, const double& m_chi,
                                          const bool& export_flag=true) ;

    /// Evaluates the escape ratio integrated over radius                                          
    double Escape_Ratio_Total(const Baryon& B,  const double& m_chi, 
                                const bool& export_flag=false) ; 

    /// Evaluates the escape ratio integrated over radius   
    /// for a range of m_chi values                                       
    Zaki::Vector::DataSet Escape_Ratio_Total(const Baryon& B, 
                                              const Zaki::Math::Axis& m_chi_range,
                                              const bool& export_flag=true) ;
             
    // -----------------------------------------------
    //        Thermal Energy Production Rates
    // -----------------------------------------------
    
    // -----------------------------------------------
    //        Fermi-Hole Heating
    // -----------------------------------------------
    /// The phase-space integrand to find thermal energy 
    /// dumped into the NS due to Fermi-hole heating
    double Thermal_Hole_PhaseSpace_Integrand(const double& x) ;
    
    /// mu_chi for 'Thermal_Hole_PhaseSpace_Integrand'
    double thermal_hole_phase_int_mu_chi = 0 ;

    /// sigma_0 for 'Thermal_Hole_PhaseSpace_Integrand'
    double thermal_hole_phase_int_sigma_0 = 0 ;

    /// x_F for 'Thermal_Hole_PhaseSpace_Integrand'
    double thermal_hole_phase_int_x_F = 0 ;

    /// The phase-space integral to find thermal energy 
    /// dumped into the NS to Fermi-hole heating
    double Thermal_Hole_PhaseSpace_Integral( const Baryon& B,
                                const double& m, 
                                const double& m_chi, 
                                const double& Sigma_0, 
                                const double& x_F) ;

    /// Returns the rate of thermal energy dumped into the NS
    /// due to Fermi-hole heating
    ///  per unit volume as a function of density in units 
    ///  of [ MeV fm^-3 s^-1 ], or [ erg cm^-3 s^-1 ].
    Zaki::Vector::DataSet Thermal_Hole_E_Rate_vs_Density(
                              const double& m_chi, 
                              const Baryon& B, 
                              const bool& gen_plots=false) ;

    /// Plots the rate of thermal energy dumped into the NS
    /// due to Fermi-hole heating
    /// as a function of density (ds) in units of s^-1/fm^3
    void hidden_Plot_Thermal_Hole_E_Rate_vs_Density(const Baryon& B, 
                                    const double& m_chi,
                                    Zaki::Vector::DataSet* ds) ;

    /// Returns the rate of thermal energy dumped into the NS
    /// due to Fermi-hole heating
    ///  per unit volume as a function of radius in units of MeV s^-1/fm^3
    Zaki::Vector::DataSet Thermal_Hole_E_Rate_vs_R(const double& m_chi, 
                                    const Baryon& B, 
                                    const bool& gen_plots=false) ;

    /// Plots rate of thermal energy dumped into the NS
    /// due to Fermi-hole heating 
    /// as a function of radius in units of s^-1/fm^3
    void Thermal_Hole_E_Rate_vs_R(const double& m_chi) ;

    /// Plots rates of thermal energy dumped into the NS
    /// due to Fermi-hole heating 
    /// for a set of chi's in the same figure 
    /// as a function of radius in units of s^-1/fm^3
    void Thermal_Hole_E_Rate_vs_R(const std::vector<double>& m_chi) ;

    /// Plots rates of thermal energy dumped into the NS
    /// due to Fermi-hole heating 
    /// for a set of chi's in the same figure 
    /// as a function of radius in units of s^-1/fm^3
    void Thermal_Hole_E_Rate_vs_R(const std::vector<double>& m_chi, const Baryon& B) ;

    /// Plots the rate of thermal energy dumped into the NS
    /// due to Fermi-hole heating 
    /// as a function of radius (ds) in units of s^-1/fm^3
    void hidden_Plot_Thermal_Hole_E_Rate_vs_R(const Baryon& B, 
                              const double& m_chi,
                              Zaki::Vector::DataSet* ds) const ;

    /// This function returns the rate of thermal energy dumped into the NS
    /// due to Fermi-hole heating
    /// integrated over the radius and in units of [ MeV s^-1 ]
    double GetThermal_Hole_E_Rate( const double& m_chi, const Baryon& B) ;

    // -----------------------------------------------
    //                Photon Heating
    // -----------------------------------------------
    /// The phase-space integrand to find thermal energy 
    /// dumped into the NS due to photons
    double Thermal_Photon_PhaseSpace_Integrand(const double& x) ;

    /// mu_chi for 'Thermal_Photon_PhaseSpace_Integrand'
    double thermal_photon_phase_int_mu_chi = 0 ;

    /// sigma_0 for 'Thermal_Photon_PhaseSpace_Integrand'
    double thermal_photon_phase_int_sigma_0 = 0 ;

    /// Sets the units to cgs
    void UseCGSUnits() ;

    /// The phase-space integral to find thermal energy 
    /// due to photons
    /// dumped into the NS due to photon
    double Thermal_Photon_PhaseSpace_Integral( const Baryon& B,
                                const double& m, 
                                const double& m_chi, 
                                const double& Sigma_0, 
                                const double& x_F) ;

    /// Plots the rate of thermal energy dumped into the NS 
    /// due to photons
    /// as a function of radius (ds) in units of s^-1/fm^3
    void hidden_Plot_Thermal_Photon_E_Rate_vs_R(const Baryon& B, 
                              const double& m_chi,
                              Zaki::Vector::DataSet* ds) const ;

    /// Plots the rate of thermal energy dumped into the NS
    /// due to photons
    /// as a function of density (ds) in units of s^-1/fm^3
    void hidden_Plot_Thermal_Photon_E_Rate_vs_Density(const Baryon& B, 
                                    const double& m_chi,
                                    Zaki::Vector::DataSet* ds) ;

    /// Returns the rate of thermal energy dumped into the NS
    /// due to photons
    ///  per unit volume as a function of density in units of MeV s^-1/fm^3
    Zaki::Vector::DataSet Thermal_Photon_E_Rate_vs_Density(
                              const double& m_chi, 
                              const Baryon& B, 
                              const bool& gen_plots=false) ;

    /// Returns the rate of thermal energy dumped into the NS
    /// due to photons
    ///  per unit volume as a function of radius in units of MeV s^-1/fm^3
    Zaki::Vector::DataSet Thermal_Photon_E_Rate_vs_R(const double& m_chi, 
                                    const Baryon& B, 
                                    const bool& gen_plots=false) ;


    /// Plots rate of thermal energy dumped into the NS
    /// due to photons
    /// as a function of radius in units of s^-1/fm^3
    void Thermal_Photon_E_Rate_vs_R(const double& m_chi) ;

    /// Plots rates of thermal energy dumped into the NS
    /// due to photons
    /// for a set of chi's in the same figure 
    /// as a function of radius in units of s^-1/fm^3
    void Thermal_Photon_E_Rate_vs_R(const std::vector<double>& m_chi) ;

    /// Plots rates of thermal energy dumped into the NS
    /// due to photons
    /// for a set of chi's in the same figure 
    /// as a function of radius in units of s^-1/fm^3
    void Thermal_Photon_E_Rate_vs_R(const std::vector<double>& m_chi, const Baryon& B) ;

    /// This function returns the rate of thermal energy dumped into the NS
    /// due to photons
    /// integrated over the radius and in units of [ MeV s^-1 ]
    double GetThermal_Photon_E_Rate( const double& m_chi, const Baryon& B) ;

    /// Plots the total rate of thermal energy dumped into the NS
    /// due to photons
    /// as a function of m_chi assuming the default eps value
    void Plot_Thermal_Photon_E_Rate() ;

    // -----------------------------------------------
    //                Total Heating
    // -----------------------------------------------
    /// Plots total rates of thermal energy dumped into the NS
    /// due to photons & Fermi-hole heating 
    /// for a set of chi's in the same figure 
    /// as a function of radius in units of s^-1/fm^3
    void Thermal_Total_E_Rate_vs_R(const std::vector<double>& m_chi) ;

    /// Plots rates of thermal energy dumped into the NS
    /// due to photons & Fermi-hole heating 
    /// for a set of chi's in the same figure 
    /// as a function of radius in units of s^-1/fm^3
    void Thermal_Total_E_Rate_vs_R(const std::vector<double>& m_chi, const Baryon& B) ;

    /// Plots the limit on total rate of thermal energy dumped into the NS
    /// as a function of m_chi assuming the limit on eps value
    void Plot_Limited_Thermal_E_Rate() const ;

    /// Plots the total rate of thermal energy dumped into the NS
    /// due to photons & Fermi-hole heating 
    /// as a function of m_chi assuming the default eps value
    void Plot_Thermal_Total_E_Rate() ;

    /// This function finds and plots (eps_max, m_chi)
    /// values that lead to a steady-state thermal solution
    ///  assuming a core temperature T_core (in kelvin)
    void Limit_eps_From_Heating(const double& T_core) ;

    // This function finds and plots (eps_max, m_chi)
    // values that lead to a BNV rate faster than nu-rates
    //  at the transition point r_durca_thresh,
    //  assuming a core temperature T_core (in kelvin)
    void Limit_eps_From_Slowness(const double& T_core) ;

    /// @brief  Solves for the mixing parameter eps given the beta-imbalance
    /// @param x beta-imbalance = eta / T
    void Eps_From_Beta_Imbalance(const double& x) ;
    // double Eps_From_Beta_Imbalance_Eq(const double& x) ;
    // -----------------------------------------------

    double PhaseSpace_Integrand(const double& x) override ;
    double PhaseSpace_Integral( const Baryon& B,
                                const double& m, 
                                const double& m_chi, 
                                const double& Sigma_0, 
                                const double& x_F) override ;

    /// Returns rate per unit volume
    ///  as a function of density in units of s^-1/fm^3
    Zaki::Vector::DataSet Rate_vs_Density(
                              const double& m_chi, 
                              const Baryon& B, 
                              const bool& gen_plots=false) 
                              override ;

    /// Returns rate per unit volume
    ///  as a function of radius in units of s^-1/fm^3
    Zaki::Vector::DataSet Rate_vs_R(const double& m_chi, 
                                    const Baryon& B, 
                                    const bool& gen_plots=false) 
                                    override ;
    /// Returns rate per unit volume
    ///  as a function of radius in units of s^-1/fm^3
    Zaki::Vector::DataSet Rate_vs_R_Slow(const double& m_chi, 
                                    const Baryon& B, 
                                    const bool& gen_plots=false) 
                                    ;     
                  
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_BNV_B_Chi_Photon_H*/
