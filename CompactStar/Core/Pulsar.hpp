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
 * @file Pulsar.hpp
 *
 * @brief Typical pulsar.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// Last edit Nov 9, 2021
#ifndef CompactStar_Pulsar_H
#define CompactStar_Pulsar_H

#include <map>
// #include <vector>
// #include <gsl/gsl_spline.h>
// #include <string>

#include <Zaki/Math/Math_Core.hpp>
#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/Prog.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

//--------------------------------------------------------------
// namespace Zaki::Math
// {
//   struct Quantity ;
// }

//--------------------------------------------------------------
// namespace Zaki::Vector
// {
//   class DataSet ;
// }

//--------------------------------------------------------------
namespace CompactStar
{

//==============================================================
//                      Baryon_Lim Struct
//==============================================================
struct Baryon_Lim
{
  // Baryon's name
  std::string name ;

  // Fraction inside the NS
  double fraction ;

  // Similar to fraction, but has the exp(nu) time-dilation factor
  double b_dot  ;

  // BNV Limit
  double limit = 0 ;
  
  Baryon_Lim( const std::string& in_name, 
              const double& in_fraction,
              const double& in_b_dot) 
  : name(in_name), fraction(in_fraction), b_dot(in_b_dot)
  {}
};
//==============================================================

//==============================================================
//                        Pulsar Class
//==============================================================
class Pulsar : public Prog
{

  friend class Decay_Analysis ;
  
  //--------------------------------------------------------------
  private:

    // Pulsar mass
    Zaki::Math::Quantity mp ;

    // Observed pulsar spin period (ms)
    Zaki::Math::Quantity spin_p ; 

    // Observed pulsar spin period derivative
    Zaki::Math::Quantity spin_p_dot ; 

    // Distance between the Solar system Barycentre (SSB)
    //  and the position of the isolated pulsar or
    //  the binary barycentre in the case of a binary pulsar.
    // Unit: kpc
    Zaki::Math::Quantity d = {0, 0} ;

    // Proper motion of the pulsar
    // Unit: milli-arcsec per year
    Zaki::Math::Quantity mu = {0, 0} ;

    // Limit on BNV form spin-down
    double bnv_spin_down_limit ;

    /// The structure profile for this pulsar
    /// Index guide:
    /// [0] : Radius (km)
    /// [1] : Mass (solar)
    /// [4] : energy density 'epsilon' (1/km^2) 
    /// [5] : Number density (fm^-3)
    /// [6] : g_tt metric function exponent (nu)
    Zaki::Vector::DataSet profile ;

    // The sequence profile that this pulsar belongs to
    Zaki::Vector::DataSet seq_profile ;

    // The point along the sequence that this pulsar belongs to
    SeqPoint seq_point ;

    // eta is defined above Eq. 33 of our paper:
    // https://arxiv.org/pdf/2201.02637.pdf
    double eta_I ;

    std::map<std::string, std::string> baryons = {
      {"10", "Neutron"},
      {"100", "Lambda"}, {"110", "Sigma-"}, {"111", "Sigma0"},
      {"112", "Sigma+"}, {"120","Xi-"}, {"121", "Xi0"}
    } ;

    bool IsBaryon(const std::string& in_label) ;

    // ----------------------------------------------
    //         Temperature (private)
    // ----------------------------------------------
    /// Pulsar's core temperature in kelvin
    /// Not red-shifted: as measured in a local frame inside the star
    double T_core = 0 ;

    /// Pulsar's blanket (outer crust) temperature in kelvin
    /// Not red-shifted: as measured in a local frame inside the star
    double T_blanket = 0 ;

    /// Pulsar's surface temperature in kelvin
    /// Not red-shifted: as measured in a local frame inside the star
    double T_surf = 0 ;

    /// @brief  The starting location of thermal blanket
    ///         defined by eps ~ 10^{10} [ g / cm^3 ]
    /// This value is set in 'FindProfile'.
    double r_blanket = 0 ;
    int r_blanket_idx = 0 ;

    /// The threshold radius for direct Urca: p_f(n) = p_F(p) + p_F(e)
    /// This value is set in 'FindProfile'.
    double r_durca_thresh = 0 ;

    /// @brief  A datacolumn filled with zeros and one to impose the
    /// r_durca_thresh condition
    /// This value is set in 'FindProfile'.
    Zaki::Vector::DataColumn r_durca_cond ;
    
    /// @brief  The prefactor in Q_nu = (factor) T_9^6 for direct Urca
    ///   unit : [ erg cm^-3 s^-1 ]
    double DUrca_emissivity_prefactor = 3e27 ;

    /// @brief  The prefactor in Q_nu = (factor) T_9^8 for direct Urca
    ///   unit : [ erg cm^-3 s^-1 ]
    double MUrca_emissivity_prefactor = 3e21 ;
    // ----------------------------------------------

  //--------------------------------------------------------------
  public:
    
    Pulsar() ;

    Pulsar( const std::string& name, 
            const Zaki::Math::Quantity& m,
            const Zaki::Math::Quantity& spin_p,
            const Zaki::Math::Quantity& spin_pdot)  ;
    ~Pulsar() ;

    // /// Copy Constructor 
    // Pulsar(const Pulsar &) = delete ;

    // /// Assignment operator
    // Pulsar& operator= (const Pulsar&) = delete ;

    /// Sets pulsar mass
    void SetMass(const Zaki::Math::Quantity& in_mass) ;
    
    /// Sets pulsar spin period in s
    void SetSpinP(const Zaki::Math::Quantity& in_spin_p) ;
    
    /// Sets pulsar spin period derivative
    void SetSpinPDot(const Zaki::Math::Quantity& in_spin_p_dot) ;

    /// Sets distance in kpc :
    /// between the Solar system Barycentre (SSB)
    ///  and the position of the isolated pulsar or
    ///  the binary barycentre in the case of a binary pulsar
    void SetDistance(const Zaki::Math::Quantity& in_d) ;

    /// Sets proper motion of the pulsar
    /// Unit: milli-arcsec per year
    void SetProperMotion(const Zaki::Math::Quantity& in_mu) ;

    int FindProfile(const std::string& model_name, 
                    const Zaki::String::Directory& in_dir="") ;

    void ImportProfile(const std::string& model_name, 
                       const Zaki::String::Directory& in_dir="") ;

    /// Sets pulsar mass
    Zaki::Math::Quantity GetMass() const ;

    /// Returns the profile of this pulsar
    Zaki::Vector::DataSet* GetProfile() ;

    /// Returns the metric exponent (nu) in g_tt = exp(2*nu)
    Zaki::Vector::DataColumn GetMetricNu() const ;

    /// Returns the metric exponent (nu) in g_tt = exp(2*nu)
    /// as a function of radius in km
    double GetMetricNu(const double& in_r) const ;

    /// Returns the sequence point that this pulsar belongs to
    SeqPoint GetSeqPoint() const ;

    /// Returns the sequence profile that this pulsar belongs to
    const Zaki::Vector::DataSet* GetSeqProfile() const ;

    void FindBNVGammaLimits() ;
    double GetBNVSpinDownLimit() ;

    /// Returns pulsar spin-period in seconds
    Zaki::Math::Quantity GetSpinP() const ;

    /// Returns pulsar spin period derivative
    Zaki::Math::Quantity GetSpinPDot() const ;

    /// @brief Pulsar spin period derivative divided by period
    /// @return (P_dot / P) [1/s]
    Zaki::Math::Quantity GetSpinPDot_over_P() const ;

    /// Returns distance in kpc :
    /// between the Solar system Barycentre (SSB)
    ///  and the position of the isolated pulsar or
    ///  the binary barycentre in the case of a binary pulsar
    Zaki::Math::Quantity GetDistance() const ;

    /// returns proper motion of the pulsar
    /// Unit: milli-arcsec per year
    Zaki::Math::Quantity GetProperMotion() const ;

    /// Evaluates the baryon fractions
    std::vector<Baryon_Lim> EvalBaryonNumber() ;

    /// Plots the pulsar's relative composition vs radius
    void PlotRelativeComposition(const std::string& file, 
                       const Zaki::String::Directory& in_dir="") ;

    /// Plots the pulsar's composition vs radius
    void PlotAbsoluteComposition(const std::string& file, 
                       const Zaki::String::Directory& in_dir="") const ;

    // /// Plots the chemical potentials vs radius
    // /// For now only for electron!
    // void PlotChemPotential(const std::string& file, 
    //                         const Zaki::String::Directory& in_dir="") const ;
    
    /// Plots the Fermi energy of particles as a function of density
    /// The input directory is relative to the wrk_dir.
    /// Only plots electron and muon because pulsar class
    /// doesn't have access to m_eff and V_eff.
    void PlotFermiE(const std::string& file, 
                    const Zaki::String::Directory& in_dir="") const ;

    /// Plots the metric function exp[ nu(r) ] vs radius
    void PlotExpNu(const std::string& file, 
                const Zaki::String::Directory& in_dir="") const ;

    /// The surface gravity g = GM e^{−nu(R)} / R^2 
    ///  in units of 10^{14} cm s^{−2}
    double Get_Surface_Gravity() const ;

    // ----------------------------------------------
    //         Temperature Methods (public)
    // ----------------------------------------------
    /// Sets pulsar's core temperature in kelvin
    /// Not red-shifted: as measured in a local frame inside the star
    void Set_T_Core(const double&) ;

    /// Sets pulsar's blanket (outer crust) temperature in kelvin
    /// Not red-shifted: as measured in a local frame inside the star
    void Set_T_Blanket(const double&) ;

    /// Sets pulsar's surface temperature (T_s) in kelvin
    /// Not red-shifted: as measured in a local frame inside the star
    void Set_T_Surf(const double&) ;

    /// Sets pulsar's apparent (red-shifted) surface temperature (T_s) in kelvin
    /// as detected by a distant observer
    void Set_T_Surf_Apparent(const double&) ;

    /// Reurns the blanket radius
    /// defined by eps ~ 10^{10} [ g / cm^3 ]
    /// This value is set in 'FindProfile'.
    double Get_R_Blanket() const ;

    /// Reurns the threshold radius for direct Urca
    /// defined by p_f(n) = p_F(p) + p_F(e)
    /// This value is set in 'FindProfile'.
    double Get_R_Durca() const ;

    /// Returns pulsar's core temperature in kelvin
    /// Not red-shifted: as measured in a local frame inside the star
    double Get_T_Core() const ;

    /// Returns pulsar's blanket (outer crust) temperature in kelvin
    /// Not red-shifted: as measured in a local frame inside the star
    double Get_T_Blanket() const ;

    /// Returns pulsar's surface temperature in kelvin
    /// Not red-shifted: as measured in a local frame inside the star
    double Get_T_Surf() const ;

    /// Returns the apparent (red-shifted) surface temperature (T-infty)
    /// as detected by a distant observer.
    double Get_T_Surf_Apparent() const ;

    /// Returns the red-shifted temperature (T-infty)
    /// as detected by a distant observer.
    /// Ignoring the blanket thickness: R_b ~ R_NS
    /// Constant: assuming that the interior of NS is thermalized
    double Get_Redshifted_T() const ;

    /// Returns the temperature as a function of r(km)
    /// as measured in a local frame inside the star.
    /// Assuming that the interior of NS is thermalized
    double Get_T(const double& in_r) const ;

    // Returns the temperature DataColumn as a function of r(km)
    // as measured in a local frame inside the star.
    // Assuming that the interior of NS is thermalized
    Zaki::Vector::DataColumn Get_T() const ;

    /// Reurns the red-shifted surface emission luminosity L^{infty}
    /// due to black body radiation
    /// Units: [erg s^-1]
    double Get_Surface_Photon_Lumin() const ;

    /// Reurns the red-shifted surface emission luminosity L^{infty}
    /// due to black body radiation in units of the red-shifted T^inf: 
    /// it has to be multiplied by "[ (T^inf_9)^2.42 ]"
    /// Units: [erg s^-1]
    double Get_Surface_Photon_Lumin_T9242() const ;

    /// Reurns the neutrino emissivity
    /// due to DUrca
    /// Units: [erg cm^-3 s^-1]
    double Get_DUrca_Neutrino_Emissivity(const double& in_radius) const ;
    Zaki::Vector::DataColumn Get_DUrca_Neutrino_Emissivity() const ;

    /// due to DUrca, in units of the red-shifted T^inf: 
    /// it has to be multiplied by "[ (T^inf_9)^6 ]"
    /// Units: [erg cm^-3 s^-1]
    Zaki::Vector::DataColumn Get_DUrca_Neutrino_Emissivity_T96() const ;

    /// Reurns the neutrino Durca rate defined by
    /// emissivity divided by [ k.T ]
    /// Units: [fm^-3 s^-1]
    double Get_DUrca_Neutrino_Rate(const double& in_radius) const ;
    Zaki::Vector::DataColumn Get_DUrca_Neutrino_Rate() const ;

    /// Reurns the neutrino emissivity
    /// due to MUrca
    /// Units: [erg cm^-3 s^-1]
    double Get_MUrca_Neutrino_Emissivity(const double& in_radius) const ;
    Zaki::Vector::DataColumn Get_MUrca_Neutrino_Emissivity() const ;

    /// due to MUrca, in units of the red-shifted T^inf: 
    /// it has to be multiplied by "[ (T^inf_9)^8 ]"
    /// Units: [erg cm^-3 s^-1]
    Zaki::Vector::DataColumn Get_MUrca_Neutrino_Emissivity_T98() const ;

    /// Reurns the neutrino Murca rate defined by
    /// emissivity divided by [ k.T ]
    /// Units: [fm^-3 s^-1]
    double Get_MUrca_Neutrino_Rate(const double& in_radius) const ;
    Zaki::Vector::DataColumn Get_MUrca_Neutrino_Rate() const ;

    /// Reurns the neutrino emissivity
    /// due to DUrca and MUrca
    /// Units: [erg cm^-3 s^-1]
    double Get_Neutrino_Emissivity(const double& in_radius) const ;
    Zaki::Vector::DataColumn Get_Neutrino_Emissivity() const ;

    /// Reurns the neutrino total rate (DUrca and MUrca) defined by
    /// emissivity divided by [ k.T ]
    /// Units: [fm^-3 s^-1]
    double Get_Neutrino_Rate(const double& in_radius) const ;
    Zaki::Vector::DataColumn Get_Neutrino_Rate() const ;

    /// Plots the neutrino emissivity
    /// due to DUrca and MUrca
    /// Units: [erg cm^-3 s^-1]
    void Plot_Neutrino_Emissivity(const std::string& file, 
                          const Zaki::String::Directory& in_dir="") const ;

    /// Plots the neutrino rate
    /// due to DUrca and MUrca
    /// Units: [ cm^-3 s^-1]
    void Plot_Neutrino_Rate(const std::string& file, 
                          const Zaki::String::Directory& in_dir="") const ;

    /// Reurns the neutrino luminosity
    /// due to DUrca and MUrca
    /// Units: [erg s^-1]
    double Get_Neutrino_Lumin() const ;

    /// Reurns the integrated neutrino rate
    /// due to DUrca in units of the red-shifted T^inf: 
    ///  it has to be multiplied by "[ (T^inf_9)^5 ]"
    /// Units: [s^-1]
    double Get_DUrca_Neutrino_Rate_T95() const ;

    /// Reurns the integrated neutrino rate
    /// due to MUrca in units of the red-shifted T^inf: 
    ///  it has to be multiplied by "[ (T^inf_9)^7 ]"
    /// Units: [s^-1]
    double Get_MUrca_Neutrino_Rate_T97() const ;
    // ----------------------------------------------
};

//==============================================================  
//--------------------------------------------------------------
} // End of namespace CompactStar
//--------------------------------------------------------------
std::ostream& operator << (std::ostream &, const CompactStar::Baryon_Lim&) ;
//--------------------------------------------------------------

#endif /*CompactStar_Pulsar_H*/
