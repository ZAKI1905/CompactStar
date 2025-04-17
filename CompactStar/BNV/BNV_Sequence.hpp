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
 * @file BNV_Sequence.hpp
 *
 * @brief Analyzes BNV evolution along the one parameter sequence
 *
 * @ingroup BNV
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// ---------------------------------------------------- 
//               Created on Apr 21, 2023
// This class is made for analyzing BNV evolution
// along the one parameter sequence
// ---------------------------------------------------- 
#ifndef CompactStar_BNV_Sequence_H
#define CompactStar_BNV_Sequence_H

#include <Zaki/Vector/DataSet.hpp>
#include <Zaki/Physics/Constants.hpp>

#include <CompactStar/Core/Prog.hpp>

// #include "CompactStar/Core/Pulsar.hpp"
// #include <CompactStar/EOS/CompOSE_EOS.hpp>

//==============================================================
namespace CompactStar
{

//==============================================================
class BNV_Sequence : public Prog
{
private:
  Zaki::Vector::DataSet seq ;
  Zaki::Vector::DataSet beta ;

  // Index guide
  struct Beta_Idx
  {
    // 0-eps
    int eps = 0 ;
    // 1-M
    int M = 1 ;
    // 2-R
    int R = 2 ;
    // 3-B
    int B = 3 ;
    // 4-I
    int I = 4 ;
    // // 5-beta(M)
    // int beta_M = 5 ;
    // // 6-beta(R)
    // int beta_R = 6 ;
    // // 7-beta(I)
    // int beta_I = 7 ;
    // 5-b(R)
    int b_R = 5 ;
    // 6-b(I)
    int b_I = 6 ;
    // 7-b(beta(I))
    int b_beta_I = 7 ;

    // 8-b(M)
    int b_M = 8 ;
  };
  
  Beta_Idx beta_idx ;

  //  The initial index based on beta dataset
  int init_idx = 0 ;

  // -----------------------------------------
  //  The initial index based on beta dataset
  // -----------------------------------------
  //   It is initialized in the constructor
  std::vector<std::pair<double, std::string>> init_omega_set ;
  bool init_omega_is_set = false ;

  // //      H = 1e8
  // // -------------------
  // // std::vector<std::pair<double, std::string>> init_omega_set 
  // //   = {{2*M_PI/5e-3, "5 (ms)"}, {2*M_PI/10e-3, "10 (ms)"}, {2*M_PI/0.1, "0.1 (s)"}} ;
  // // -------------------
  // //      H = 1e9
  // // -------------------
  // std::vector<std::pair<double, std::string>> init_omega_set 
  //   = {{2*M_PI/10e-3, "10 (ms)"}, {2*M_PI/30e-3, "30 (ms)"}, {2*M_PI/0.1, "0.1 (s)"}} ;
  // // -------------------
  // //      H = 1e10
  // // -------------------
  // // std::vector<std::pair<double, std::string>> init_omega_set 
  // //   = {{2*M_PI/1e-1, "0.1 (s)"}, {2*M_PI/5e-1, "0.5 (s)"}, {2*M_PI/1, "1 (s)"}} ;
  // // -------------------
  // //      H = 1e11
  // // -------------------
  // // std::vector<std::pair<double, std::string>> init_omega_set 
  // //   = {{2*M_PI/1, "1 (s)"}, {2*M_PI/5, "5 (s)"}} ;
  // // -----------------------------------------

  // The strength of magnetic field in gauss
  double mag_field = 1e8 ;

  // int Pdot_sign = -1 ;
  // std::string Pdot_sign_str = "+" ;

  // The time resolution for ODE
  double time_res = 1e3 ;

  // The BNV rate per seconds 
  double gamma_bnv = 1e-10 / Zaki::Physics::YR_2_SEC ;
  
  // The time at which BNV turns off
  double bnv_cutoff_time = 1e25 ;

  // The mass at which BNV shuts down
  double bnv_cuttoff_mass = 0.1 ;

  // Records the time when BNV turns off
  void SetBNVCuttoff(const double& time) ;

  bool bnv_active_flag = true ;

  // -----------------------------------------------------
  //             File names of evolution plots
  // -----------------------------------------------------
  /// File names for braking index vs time evolution 
  /// (see "void Solve(t_i, t_f)" )
  std::string br_idx_evol_f_name = "BrIdx_Evolution/n_vs_t" ;

  /// File names for omega vs time evolution (see "void Solve(t_i, t_f)" )
  std::string omega_t_evol_f_name = "Omega_vs_t/Omega_vs_t" ;

  /// File names for age vs time evolution (see "void Solve(t_i, t_f)" )
  std::string age_t_evol_f_name = "Fig_9_age_vs_t/Age_vs_t" ;

  /// File names for omega vs mass evolution (see "void Solve(t_i, t_f)" )
  std::string omega_m_evol_f_name = "Fig_5_omega_vs_m/Omega_vs_m" ;
  // -----------------------------------------------------


public:

  /// Defult Constructor
  BNV_Sequence();

  /// Constructor from a sequence file
  BNV_Sequence(const Zaki::String::Directory& in_dir, const std::string& in_f_name) ;

  /// Destructor
  ~BNV_Sequence();

  // /// Copy constructor
  // BNV_Sequence(const BNV_Sequence&) ;

  // /// Assignment operator
  // BNV_Sequence& operator=(const BNV_Sequence&) ;

  /** Generates a sequence of neutron stars.
  *
  *
  * @param in_dir directory
  * @param in_model model number
  */
  void GenSequence( const Zaki::String::Directory& in_dir, 
                    const std::string& in_model) ;

  /// Imports the sequence
  void Import(const Zaki::String::Directory& in_dir, 
              const std::string& in_f_name) ;


  /// Evaluates b_O factor using the central difference method
  ///  and given the input NS mass
  void Find_b_O(const double& in_mass, const int& o_idx) ;

  /// Evaluates b_O factor using the central difference method
  void Find_b_factors() ;

  ///
  void Plot_Dimless_O() const ;

  /// Evaluates the beta factors
  void EvalBeta() ;

  /// Input in yr^-1
  void SetBNVRate(const double& in_gamma) ;

  /// Returns the observable with the given index
  /// Zaki::Vector::DataColumn O(const int&) const ;

  /// Sets the mass at which BNV shuts down
  void SetBNVCuttOffMass(const double& in_M) ;

  /// Returns the central energy density
  Zaki::Vector::DataColumn Eps_C() const ;

  /// Returns Mass
  Zaki::Vector::DataColumn M() const ;

  /// Returns quantity O as a function of time
  double O(const double& time, const int& o_idx) const  ;

  /// Returns mass as a function of time
  double M(const double& time) const ;

  /// Returns MomI
  Zaki::Vector::DataColumn I() const ;

  /// Returns I as a function of time
  double I(const double& time) const ;

  /// Returns baryon number
  Zaki::Vector::DataColumn B() const ;

  /// Returns baryon number as a function of time
  double B(const double& time) const ;

  /// Returns the closest index corrsponding to time
  int Time_to_Idx(const double& time) const ;

  struct weighted_idx
  {
    // This is the index
    int idx ;

    // weight to the left of idx
    double w_L ;

    // weight to the right of idx
    double w_R ;
  } ;

  /// Returns the weighted index corrsponding to time
  weighted_idx Time_to_Weighted_Idx(const double& time) const ; 

  /// Returns Radius
  Zaki::Vector::DataColumn R() const ;

  /// Returns R as a function of time
  double R(const double& time) const ;

  /// Returns Beta_I
  // Zaki::Vector::DataColumn Beta_I() const ;

  /// Returns Beta_R
  // Zaki::Vector::DataColumn Beta_R() const ;

  /// Returns Beta_M
  // Zaki::Vector::DataColumn Beta_M() const ;

  /// Returns b_I
  Zaki::Vector::DataColumn b_I() const ;

  /// Returns b_I as a function of time
  double b_I(const double& time) const ;

  /// Returns b(beta_I)
  Zaki::Vector::DataColumn b_beta_I() const ;

  /// Returns b(beta_I) as a function of time
  double b_beta_I(const double& time) const ;

  /// Returns b_R
  Zaki::Vector::DataColumn b_R() const ;

  /// Returns b_R as a function of time
  double b_R(const double& time) const ;

  /// Returns b_M
  Zaki::Vector::DataColumn b_M() const ;

  /// Returns Beta(O)
  Zaki::Vector::DataColumn Beta(const Zaki::Vector::DataColumn& in_dc) const ;

  /// Returns b(O)
  Zaki::Vector::DataColumn b(const Zaki::Vector::DataColumn& in_dc) const ;

  /// Returns b(O) as a function of time
  double b(const Zaki::Vector::DataColumn& in_dc, const double& time) const ;

  /// Returns the critical (death) omega assuming a central dipole H-field (per sec)
  double Omega_Death_Central(const double& time) ;

  /// Returns the critical (death) omega assuming a twisted dipole H-field (per sec)
  double Omega_Death_Twisted(const double& time) ;

  /// Returns the breaking index
  Zaki::Vector::DataColumn BrIdx(const double& in_delta) const ;

  /// Returns the breaking index as a function of time and omega
  double BrIdx(const double& time, const double& omega) const ;

  /// Returns the breaking index as a function of mass and omega
  double BrIdx_M_Omega(double mass, double omega) ;

  int FindInitIdx(const double& in_mass) ;

  /// Returns the strength of magnetic field in gauss
  double MagField() const ;

  /// Sets the strength of magnetic field in gauss
  void SetMagField(const double& in_B) ;

  /// Set the initial omega set
  void SetInitOmega() ;
  
  /// Returns (H^2 R^6 / I) in seconds
  ///  as a function of time 
  double ODE_Coeff(const double&) const ;

  /// Returns (H^2 R^6 / I) in seconds
  ///  as a DataColumn 
  Zaki::Vector::DataColumn H2R6_I() const ;

  /// Returns the BNV rate per seconds
  double Gamma_BNV(const double& time) const ;

  /// Returns the extremum (max/min) angular velocity
  double ExtremumOmega(const double& time) const ;

  /// Returns a dataset with two columns:
  /// [0]: Mass in solar mass
  /// [1]: Extremum value for angular velocity (per seconds)
  Zaki::Vector::DataSet ExtremumOmega() const  ;

  /// Returns the threshold magnetic field at which
  /// the SD-SU transition occurs in gauss
  double ThresholdMagField(const double& bnv_rate, const double& P, const double& mass=1.4) const ;

  // ......................
  // Dictionary :
  // y[0] = Omega(t)
  // f[0] = Omega'(t)
  // ......................
  static int ODE(double t, const double y[], double f[], void *params) ;

  void Solve(const double t_0, const double t_f) ;


  void SetInitOmegaSet(const std::vector<std::pair<double, std::string>>& init_omegas) ;
 
  // -----------------------------------------------------
  //     Setters for file names of evolution plots
  // -----------------------------------------------------
  /// Sets file names for braking index vs time evolution 
  /// (see "void Solve(t_i, t_f)" and br_idx_evol_f_name)
  void SetBrIdxvsTimeFileName(const std::string& f_name) ;
  /// File names for omega vs time evolution
  /// (see "void Solve(t_i, t_f)" and omega_t_evol_f_name)
  void SetOmegavsTimeFileName(const std::string& f_name) ;
  /// File names for age vs time evolution
  /// (see "void Solve(t_i, t_f)" and age_t_evol_f_name)
  void SetAgevsTimeFileName(const std::string& f_name) ;
  /// File names for omega vs mass evolution
  /// (see "void Solve(t_i, t_f)" and omega_m_evol_f_name)
  void SetOmegavsMFileName(const std::string& f_name) ;
  // -----------------------------------------------------

  /// Returns the second derivative of omega as a function of 
  ///  time and omega 
  double Omega_2ndDer(const double& time, const double& omega) const ;

  // Returns the second derivative of omega (DataColumn) 
  // as a function of omega 
  Zaki::Vector::DataColumn Omega_2ndDer(const double& omega) const ;

  /// Returns the first derivative of omega as a function of 
  ///  time and omega
  double Omega_1stDer(const double& time, const double& omega) const ;

  // Returns the first derivative of omega (DataColumn) 
  // as a function of omega 
  Zaki::Vector::DataColumn Omega_1stDer(const double& omega) const ;

  /// Returns the critical period at which the pulsar turns off
  ///  Taken from Eq. (39) of 
  ///  "Theory of Pulsars: Polar Gaps, Sparks, and Coherent Microwave Radiation"
  ///  by Ruderman & Sutherland (1974)
  double DeathPeriod(const double& t) const ;

  /// Returns the critical omgea at which the pulsar turns off
  ///  Taken from Eq. (39) of 
  ///  "Theory of Pulsars: Polar Gaps, Sparks, and Coherent Microwave Radiation"
  ///  by Ruderman & Sutherland (1974)
  double DeathOmega(const double& time) const ;

  /// Returns the second derivative of frequency (at t=0)
  /// as a function of the nu and nu_dot
  double Freq_ddot(const double& nu, const double& nu_dot) const ;

  double Gamma_Max() const ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_BNV_Sequence_H*/
