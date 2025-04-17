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
 * @file CompOSE_EOS.hpp
 *
 * @brief For working with CompOSE EoSs.
 *
 * @ingroup EOS
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// -------------------------------------------------------------
//                        CompOSE_EOS Class
// -------------------------------------------------------------
// Created on October 29, 2022
//
// Purpose: 
//  To simplify working with the files for the equation of 
//  states from the CompOSE website:
//        compose.obspm.fr
// -------------------------------------------------------------

#ifndef CompactStar_CompOSE_EOS_H
#define CompactStar_CompOSE_EOS_H

#include <map>

#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/Prog.hpp"

//==============================================================

namespace CompactStar
{

//==============================================================
/// Is the directory path relative to another path (wrk_dir) or
///  absolute ?
enum class Dir_Type
{
  relative = 0, absolute = 1
};

//==============================================================
class CompOSE_EOS : public Prog
{
  //--------------------------------------------------------------
  private:

  /// Indicates if the EoS has a crust
  bool has_crust = false ;

  /// crust core transition index
  int crust_core_x_idx = -1 ;
  
  /// crust core transition density
  double crust_core_x_den = -1 ;

  /// Constant neutron and proton masses
  double m_n, m_p ;

  /// The eos dataset
  /// eos[0] = "e(g/cm^3)" ;
  /// eos[1] = "p(dyne/cm^2)" ;
  /// eos[2] = "rho(1/fm^3)" ;
  Zaki::Vector::DataSet eos ;

  /// The m_eff dataset
  Zaki::Vector::DataSet m_eff ;

  /// The V_eff dataset
  Zaki::Vector::DataSet V_eff ;

  /// The single-particle potential (U_i) dataset
  Zaki::Vector::DataSet U ;

  /// Size of the density grid
  size_t grid_size = 0 ;

  /// The number of thermo columns
  size_t therm_size = 3 ;

  /// The number of composition columns
  size_t comp_size = 0 ;

  /// The index dictionary 
  enum class EOS_Idx
  {
    e = 0, p = 1, n = 2
  };
  
  struct CompOSE_Particle
  {
    /// Particle's name
    std::string name ;

    /// Mass in MeV
    double m ;

    /// Constructor
    CompOSE_Particle(const std::string& in_name, 
                      const double& in_m)
                      : name(in_name), m(in_m)
    {}
  };
  
  //--------------------------------------------------------------
  public:

  /// @brief The mapping between CompOSE codes with labels
  /// Taken from Table 3.3 in CompOSE manual (v3.00) from
  ///   https://compose.obspm.fr/manual
  const static std::map<int, CompOSE_Particle> Compose_Dict ;

    /// Constructor-0
    CompOSE_EOS() ;

    /// Destructor
    ~CompOSE_EOS() ;

    /// Copy constructor
    CompOSE_EOS(const CompOSE_EOS& ) = delete ;

    /// Assignment operator
    CompOSE_EOS& operator=(const CompOSE_EOS& ) = delete ;

    /// Sets the crust-core transition density
    void SetCrustCoreXDensity(const double&) ;

    /// Imports the EoS grid data
    /// @param in_dir Input director
    /// @param in_dir_type Relative to wrk_dir or absolute?
    void ImportGrid(const Zaki::String::Directory& in_dir, 
                          // const bool& gen_plots=false, 
                const Dir_Type in_dir_type = Dir_Type::relative ) ;

    /// Imports the ".thermo" file
    /// @param in_dir Input director
    /// @param in_dir_type Relative to wrk_dir or absolute?
    void ImportThermo(const Zaki::String::Directory& in_dir, 
                          const bool& gen_plots=false, 
                const Dir_Type in_dir_type = Dir_Type::relative ) ;

    /// Imports the ".compo" file
    /// @param in_dir Input director
    /// @param in_dir_type Relative to wrk_dir or absolute?
    void ImportCompo(const Zaki::String::Directory& in_dir, 
                          const bool& gen_plots=false, 
                const Dir_Type in_dir_type = Dir_Type::relative ) ;
    
    /// Imports the ".micro" file
    /// @param in_dir Input director
    /// @param in_dir_type Relative to wrk_dir or absolute?
    void ImportMicro(const Zaki::String::Directory& in_dir, 
                          const bool& gen_plots=false, 
                const Dir_Type in_dir_type = Dir_Type::relative ) ;

    /// Extends the M_eff file to the crust region by
    /// linearly extrapolating
    /// The path is by default absolute (not relative to wrk_dir) !  
    void ExtendMeffToCrust(const Zaki::String::Directory& in_dir, 
                          const std::string& f_name,
                          const bool& gen_plots=false, 
                const Dir_Type in_dir_type = Dir_Type::absolute) ;

    /// Extends the V_eff file to the crust region by
    /// linearly extrapolating
    /// The path is by default absolute (not relative to wrk_dir) !  
    void ExtendVeffToCrust(const Zaki::String::Directory& in_dir, 
                          const std::string& f_name,
                          const bool& gen_plots=false, 
                const Dir_Type in_dir_type = Dir_Type::absolute) ;

    /// Extends the U_eff file to the crust region by
    /// linearly extrapolating
    /// The path is by default absolute (not relative to wrk_dir) !  
    void ExtendUeffToCrust(const Zaki::String::Directory& in_dir, 
                          const std::string& f_name,
                          const bool& gen_plots=false, 
                const Dir_Type in_dir_type = Dir_Type::absolute) ;

    /// @brief Import the EOS in the standard format [ not Compose! ]
    /// imports everything: e, p, n, composition, m_eff, V_eff, U_eff.
    /// @param in_file Input file
    void ImportEOS(const Zaki::String::Directory& in_file, 
                          const bool& gen_plots=false) ;

    Zaki::Vector::DataSet GetEOS() const ;
    Zaki::Vector::DataColumn GetEOS(const int& idx) const ;
    Zaki::Vector::DataColumn GetEOS(const std::string& label) const ;

    /// Returns the effective mass dataset
    Zaki::Vector::DataSet GetMeff() const ;

    /// Returns the effective mass from index
    Zaki::Vector::DataColumn GetMeff(const int& idx) const ;

    /// Returns the effective mass from label
    Zaki::Vector::DataColumn GetMeff(const std::string& label) const ;

    /// Returns the self-energy dataset
    Zaki::Vector::DataSet GetVeff() const ;

    /// Returns the self-energy from index
    Zaki::Vector::DataColumn GetVeff(const int& idx) const ;

    /// Returns the self-energy from label
    Zaki::Vector::DataColumn GetVeff(const std::string& label) const ;


    /// Returns the single-particle potential dataset
    Zaki::Vector::DataSet GetU() const ;

    /// Returns the single-particle potential from index
    Zaki::Vector::DataColumn GetU(const int& idx) const ;

    /// Returns the single-particle potential from label
    Zaki::Vector::DataColumn GetU(const std::string& label) const ;

    /// Returns the neutron mass constant
    double GetNeutronMass() const ;

    /// Returns the proton mass constant
    double GetProtonMass() const ;

    /// Plots the Fermi energy of particles as a function of density
    /// The input directory should by default an absolute path.
    void PlotFermiE(const Zaki::String::Directory& in_dir,
                    const Dir_Type in_dir_type = Dir_Type::absolute) const ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_CompOSE_EOS_H*/
