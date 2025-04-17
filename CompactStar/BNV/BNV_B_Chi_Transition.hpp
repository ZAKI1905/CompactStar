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
 * @file BNV_B_Chi_Transition.hpp
 *
 * @brief Evaluates B -> chi transition rates in neutron stars [wrong]
 *
 * @ingroup BNV
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// ---------------------------------------------------- 
//               Created on Mar 13, 2023
// This class is made specifically for analyzing the
//  B --> Chi
//  transitions in neutron stars.
// ---------------------------------------------------- 
#ifndef CompactStar_BNV_B_Chi_Transition_H
#define CompactStar_BNV_B_Chi_Transition_H

#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/Pulsar.hpp"
#include <CompactStar/EOS/CompOSE_EOS.hpp>
#include <CompactStar/BNV/BNV_Chi.hpp>


//==============================================================
namespace CompactStar
{

//==============================================================
class BNV_B_Chi_Transition : public BNV_Chi
{
  //--------------------------------------------------------------
  private:

    Process GetSpecificProcess(const Baryon& B) const override ;

  public:

    BNV_B_Chi_Transition();
    ~BNV_B_Chi_Transition();

    /// Returns the lower limit on the vector 
    /// self-energy, from condition (a).
    /// The arguments are (m_chi, m_B^*, n [fm^-3] ).
    Zaki::Vector::DataColumn SigmaMinus(const double& m_chi, 
                                        const Baryon& B) ;

    // Returns the lower limit on the vector 
    // self-energy, from condition (a).
    // The arguments are (m_chi, m_B^*, n [fm^-3] ).
    double SigmaMinus(const double& m_chi, 
                      const double& m_B, 
                      const double& n) ;

    // /// Returns the upper limit on the vector 
    // /// self-energy, from condition (b). 
    // /// The arguments are (m_chi, m_B^*, n [fm^-3] ).
    // Zaki::Vector::DataColumn SigmaPlus(const double& m_chi, 
    //                                   const Baryon& B) ;

    // /// Returns the upper limit on the vector 
    // /// self-energy, from condition (b). 
    // /// The arguments are (m_chi, m_B^*, n [fm^-3] ).
    // double SigmaPlus(const double& m_chi, 
    //                   const double& m_B, 
    //                   const double& n) ;
  
  /// Returns the B -> chi decay rate per unit volume
  ///  as a function of density in units of s^-1/fm^3
  Zaki::Vector::DataSet Rate_vs_Density(
                            const double& m_chi, 
                            const Baryon& B, 
                            const bool& gen_plots=false)
                            override ;

  /// Returns the B -> chi decay rate per unit volume
  ///  as a function of radius in units of s^-1/fm^3
  Zaki::Vector::DataSet Rate_vs_R(const double& m_chi, 
                                  const Baryon& B, 
                                  const bool& gen_plots=false) 
                                  override ;

    /// Plots the transition conditions 
    /// for a given m_chi as a function of density                  
    void PlotTransCond(const double& m_chi,
                          const Baryon& B) ;

    /// Plots the transition conditions 
    /// for a given m_chi as a function of density                  
    void PlotTransCond(const double& m_chi) ;

    /// Plots the transition band 
    ///  as a function of m_chi                  
    void PlotTransBand(const Baryon& B) ;

    // Plots the transition band as a function of m_chi                  
    void PlotTransBand() ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_BNV_B_Chi_Transition_H*/
