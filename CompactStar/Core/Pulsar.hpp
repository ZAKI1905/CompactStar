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
 * @brief Detailed pulsar model including rotational, structural, and thermal physics.
 *
 * Manages spin evolution, baryon violation constraints, structural profiles from TOV solutions,
 * and cooling via direct and modified Urca neutrino processes as well as photon surface emission.
 * Provides utilities for profiling, plotting, and computing astrophysical observables.
 *
 * Physical context:
 *  - Direct Urca (DUrca): rapid neutrino cooling when proton fraction exceeds threshold, enabling
 *    beta decay and inverse beta processes without spectator nucleons; emissivity ∝ T^6.
 *  - Modified Urca (MUrca): slower neutrino cooling involving an additional nucleon spectator,
 *    requiring energy-momentum conservation with reduced phase space; emissivity ∝ T^8.
 *  - Thermal blanket: the outer crust region where the heat blanket regulates surface temperature,
 *    defined by a density ~1e10 g/cm^3.
 *
 * @ingroup Core
 * @author Mohammadreza Zakeri
 * @contact M.Zakeri@eku.edu
 */
// Last edit Nov 9, 2021
#ifndef CompactStar_Pulsar_H
#define CompactStar_Pulsar_H

#include <map>

#include <Zaki/Math/Math_Core.hpp>
#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/Prog.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "CompactStar/EOS/CompOSE_EOS.hpp"

//==============================================================
namespace CompactStar
{

//==============================================================
//                      Baryon_Lim Struct
//==============================================================
/**
 * @struct Baryon_Lim
 * @brief Tracks baryon number violation (BNV) constraints for a species.
 *
 * Models the fractional abundance and BNV rate of a baryon species within the star,
 * including gravitational redshift of rates and experimental/theoretical limits.
 *
 * Physical definitions:
 *  - fraction: number density fraction of species relative to total baryon density.
 *  - b_dot: proper-time baryon number change rate including redshift factor e^{nu(r)},
 *    accounting for slower clocks deeper in the gravitational well.
 *  - limit: upper bound on BNV rate from particle physics experiments (e.g., di-nucleon decay).
 */
struct Baryon_Lim
{
  // Baryon's name
  std::string name ; /**< Baryon species name (e.g., "Neutron", "Lambda"). */

  // Fraction inside the NS
  double fraction ;  /**< Fractional abundance in neutron star interior. */

  // Similar to fraction, but has the exp(nu) time-dilation factor
  double b_dot  ; /**< Redshifted baryon number violation rate [s^-1]. */

  // BNV Limit
  double limit = 0 ; /**< Experimental/theoretical upper limit on BNV rate. */
  
  /**
   * @brief Construct a Baryon_Lim object.
   * @param in_name     Name of baryon species.
   * @param in_fraction Fractional abundance.
   * @param in_b_dot    Proper-time BNV rate including redshift.
   */
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

/**
 * @class Pulsar
 * @brief Comprehensive pulsar representation with rotation, structure, composition, and thermal evolution.
 *
 * Extends Prog for logging and directory management. Handles:
 *  - Spin and spin-down observables.
 *  - Baryon number violation constraints via spin-down energy budget.
 *  - Structural profiles from TOV solver: mass, radius, metric functions.
 *  - Composition profiles: baryon species fractions and Fermi energies.
 *  - Thermal physics: core, blanket, and surface temperatures, photon luminosity.
 *  - Neutrino cooling: DUrca and MUrca processes with physical emissivity prefactors.
 *  - Plotting utilities for composition, metric, and cooling diagnostics.
 *
 * See arXiv:2201.02637 for definitions of η_I and spin-down BNV limits.
 */
class Pulsar : public Prog
{

  friend class Decay_Analysis ;
  
  //--------------------------------------------------------------
  private:
    // Spin and mass properties

    // Pulsar mass
    Zaki::Math::Quantity mp ; /**< Pulsar mass (M_sun). */

    // Observed pulsar spin period (ms)
    Zaki::Math::Quantity spin_p ; /**< Spin period P (s). */

    // Observed pulsar spin period derivative
    Zaki::Math::Quantity spin_p_dot ;  /**< Spin period derivative \dot{P}. */

    // Distance between the Solar system Barycentre (SSB)
    //  and the position of the isolated pulsar or
    //  the binary barycentre in the case of a binary pulsar.
    // Unit: kpc
    Zaki::Math::Quantity d = {0, 0} ; /**< Distance from SSB (kpc). */

    // Proper motion of the pulsar
    // Unit: milli-arcsec per year
    Zaki::Math::Quantity mu = {0, 0} ; /**< Proper motion (mas/yr). */

    // Limit on BNV form spin-down
    double bnv_spin_down_limit ;  /**< Upper bound on BNV rate from energy loss. */

    /**
     * @brief Structure profile of the pulsar.
     * @details Contains radius, mass, energy density, baryon density,
     *          and metric function values.
     * @note The profile is indexed as follows:
     *  [0] : Radius (km)
     *  [1] : Mass (solar mass)
     *  [4] : Energy density (1/km^2)
     *  [5] : Baryon density (fm^-3)
     *  [6] : Metric function g_tt exponent (nu)
     */
    Zaki::Vector::DataSet profile ;

    /**
     * @brief The sequence profile that this pulsar belongs to.
     * 
     */
    Zaki::Vector::DataSet seq_profile ;

    // The point along the sequence that this pulsar belongs to
    SeqPoint seq_point ;

    CompOSE_EOS* eos = nullptr ; /**< Pointer to the EOS object. */

    /**
     * @brief The eta parameter for this pulsar.
     * @details eta is defined above Eq. 33 of our paper:
     *  https://arxiv.org/pdf/2201.02637.pdf
     */
    double eta_I ;

    /**
     * @brief Baryon species mapping
     * @details Maps baryon labels to their names.
     * 
     */
    std::map<std::string, std::string> baryons = {
      {"10", "Neutron"},
      {"100", "Lambda"}, {"110", "Sigma-"}, {"111", "Sigma0"},
      {"112", "Sigma+"}, {"120","Xi-"}, {"121", "Xi0"}
    } ;

    /**
     * @brief Check if a label corresponds to a baryon species.
     * @param in_label Label string from profile.
     * @return True if label is a baryon code.
     */
    bool IsBaryon(const std::string& in_label) ;

    // ----------------------------------------------
    //         Temperature (private)
    // ----------------------------------------------
    /**
     * @brief Core temperature of the pulsar in K.
     * @details Not red-shifted: as measured in a local frame inside the star.
     */
    double T_core = 0 ;

    /// Pulsar's blanket (outer crust) temperature in kelvin
    /// Not red-shifted: as measured in a local frame inside the star
    /**
     * @brief Temperature of the thermal blanket.
     * @details Not red-shifted: as measured in a local frame inside the star.
     * @note The blanket is the outer crust region where the heat blanket regulates surface temperature.
     *       Defined by a density ~1e10 g/cm^3.
     */
    double T_blanket = 0 ;

    /**
     * @brief Surface temperature of the pulsar.
     * @details Not red-shifted: as measured in a local frame inside the star.
     */
    double T_surf = 0 ;

    /**
     * @brief Radius of the thermal blanket region.
     * @details Defined by energy density ~1e10 g/cm^3.
     *         This value is set in 'FindProfile'.
     */
    double r_blanket = 0 ;

    /**
     * @brief Index of the thermal blanket region.
     * @details This index is used to access the thermal blanket region in the profile.
     * 
     */
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
    /**
     * @brief Default constructor for a blank Pulsar.
     */
    Pulsar() ;

    /**
     * @brief Parameterized constructor.
     * @param name       Identifier label for logging.
     * @param m          Mass quantity.
     * @param spin_p     Spin period quantity.
     * @param spin_pdot  Spin period derivative quantity.
     */
    Pulsar( const std::string& name, 
            const Zaki::Math::Quantity& m,
            const Zaki::Math::Quantity& spin_p,
            const Zaki::Math::Quantity& spin_pdot)  ;


    /**
     * @brief Destructor cleans up internal resources.
     */
    ~Pulsar() ;

    // /// Copy Constructor 
    // Pulsar(const Pulsar &) = delete ;

    // /// Assignment operator
    // Pulsar& operator= (const Pulsar&) = delete ;

    /**
     * @brief Set the pulsar mass.
     * @param in_mass Mass quantity.
     */
    void SetMass(const Zaki::Math::Quantity& in_mass) ;
    
    /**
     * @brief Set the spin period P.
     * @param in_spin_p Spin period quantity [s].
     */
    void SetSpinP(const Zaki::Math::Quantity& in_spin_p) ;
    
    /**
     * @brief Set the spin period derivative \dot{P}.
     * @param in_spin_p_dot Spin period derivative quantity.
     */
    void SetSpinPDot(const Zaki::Math::Quantity& in_spin_p_dot) ;

    /**
     * @brief Set distance from Solar System Barycentre in kpc. 
     * @param in_d Distance quantity [kpc].
     * @details Distance is between the Solar system Barycentre (SSB)
     *          and the position of the isolated pulsar or
     *          the binary barycentre in the case of a binary pulsar.
     */
    void SetDistance(const Zaki::Math::Quantity& in_d) ;

    /**
     * @brief Set proper motion \mu of the pulsar.
     * @param in_mu Proper motion quantity [mas/yr].
     */
    void SetProperMotion(const Zaki::Math::Quantity& in_mu) ;

    // --------------------------
    // Profile loading
    // --------------------------

    /**
     * @brief Compute structural profile via TOV solver.
     * @param model_name Name of EOS/sequence model.
     * @param in_dir     Directory for files (optional).
     * @return Number of radial points in profile.
     */
    int FindProfile(const std::string& model_name, 
                    const Zaki::String::Directory& in_dir="") ;

    /**
     * @brief Import a precomputed TOV profile from disk.
     * @param model_name Name of EOS/sequence model.
     * @param in_dir     Directory for input files.
     */       
    void ImportProfile(const std::string& model_name, 
                       const Zaki::String::Directory& in_dir="") ;

    /**
     * @brief Attach an EOS to the pulsar.
     * @param eos        Pointer to CompOSE_EOS object. 
     */
    void AttachEOS(CompOSE_EOS* eos) ;
    // --------------------------
    // Getters
    // --------------------------

    /**
     * @brief Get pulsar mass.
     * @return Mass quantity.
     */
    Zaki::Math::Quantity GetMass() const ;

    /**
     * @brief Access internal structural profile.
     * @return Pointer to DataSet of radius, mass, etc.
     */
    Zaki::Vector::DataSet* GetProfile() ;
    
    /**
     * @brief Returns the EOS object associated with this pulsar.
     * 
     * @return CompOSE_EOS* 
     */
    CompOSE_EOS* GetEOS() const ; 

    /**
     * @brief Retrieve metric exponent ν for g_tt = e^{2ν}.
     * @return DataColumn of ν values vs radius.
     */
    Zaki::Vector::DataColumn GetMetricNu() const ;

    /**
     * @brief Get metric exponent ν at a specific radius.
     * @param in_r Radius [km].
     * @return ν value.
     */
    double GetMetricNu(const double& in_r) const ;

    /// Returns the sequence point that this pulsar belongs to
    SeqPoint GetSeqPoint() const ;

    /**
     * @brief Get the sequence profile.
     * @return Pointer to DataSet of sequence profile.
     * @details Contains mass, radius, and other properties along the sequence.
     */
    const Zaki::Vector::DataSet* GetSeqProfile() const ;

    // --------------------------
    // BNV analysis
    // --------------------------

    /**
     * @brief Find baryon number violation limits.
     * 
     */
    void FindBNVGammaLimits() ;

    /**
     * @brief Get computed BNV spin-down limit.
     * @return BNV rate limit [yr^-1].
     */
    double GetBNVSpinDownLimit() ;

    // --------------------------
    // Spin observables
    // --------------------------

    /**
     * @brief Get the spin period P.
     * @return Spin period quantity [s].
     */
    Zaki::Math::Quantity GetSpinP() const ;

    /**
     * @brief Get the spin period derivative \dot{P}.
     * @return Spin period derivative quantity [s/s].
     */
    Zaki::Math::Quantity GetSpinPDot() const ;

    /**
     * @brief Get the spin period derivative divided by the period.
     * @return Spin period derivative divided by period quantity [1/s].
     * @details This is the fractional change in spin period per unit time.
     */
    Zaki::Math::Quantity GetSpinPDot_over_P() const ;
    
    // --------------------------
    // Distance & motion
    // --------------------------

    /**
     * @brief Get distance from Solar System Barycentre.
     * @return Distance quantity [kpc].
     * @details Distance is between the Solar system Barycentre (SSB)
     *          and the position of the isolated pulsar or
     *          the binary barycentre in the case of a binary pulsar.
     */
    Zaki::Math::Quantity GetDistance() const ;

    /**
     * @brief Get the proper motion of the pulsar.
     * @return Proper motion quantity [mas/yr].
     * @details Proper motion is the angular motion of the pulsar across the sky.
     */
    Zaki::Math::Quantity GetProperMotion() const ;

    /**
     * @brief Evaluate baryon number fractions and limits.
     * @return Vector of baryon limits with species name, fraction, and BNV rate.
     */
    std::vector<Baryon_Lim> EvalBaryonNumber() ;

    // --------------------------
    // Composition & plotting
    // --------------------------

    /**
     * @brief Plot relative composition of baryon species vs radius.
     * @param file Output filename for the plot.
     * @param in_dir Directory for input files (optional).
     * @details Only includes baryon species with significant contributions.
     */
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

/**
 * @brief Stream output for Baryon_Lim struct.
 * @param os Output stream.
 * @param bl Baryon_Lim object.
 * @return Reference to output stream.
 *
 * Formats: name, fraction, b_dot, limit.
 */
std::ostream& operator << (std::ostream &, const CompactStar::Baryon_Lim&) ;
//--------------------------------------------------------------

#endif /*CompactStar_Pulsar_H*/
