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
 * @file BNV_B_Chi_Combo.hpp
 *
 * @brief Combined model of B -> chi gamma and B -> chi [wrong]
 *
 * @ingroup BNV
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/
// ---------------------------------------------------- 
//               Created on Mar 15, 2023
// ---------------------------------------------------- 
#ifndef CompactStar_BNV_B_Chi_Combo_H
#define CompactStar_BNV_B_Chi_Combo_H

#include <Zaki/Vector/DataSet.hpp>

#include <CompactStar/BNV/BNV_Chi.hpp>
#include <CompactStar/EOS/CompOSE_EOS.hpp>
#include "CompactStar/Core/Pulsar.hpp"


//==============================================================
namespace CompactStar
{

//==============================================================
class BNV_B_Chi_Combo : public BNV_Chi
{
  //--------------------------------------------------------------
  private:
    std::vector<BNV_Chi*> chi_reactions ;
  //--------------------------------------------------------------
  public:

    BNV_B_Chi_Combo();
    ~BNV_B_Chi_Combo();

    Process GetSpecificProcess(const Baryon& B) const override ;

    /// Sets the work directory for the member objects
    Prog* SetMemWrkDir(const Zaki::String::Directory&) override ;

    /// Sets the EoS model
    void SetModel(const std::string& in_eos_model) override ;

    /// Imports the EOS
    void ImportEOS(const Zaki::String::Directory& eos_dir, 
                   const std::string& micro_model="") override ;

    /// Finds the pulsar profile
    void FindPulsar(const bool& gen_plots=false) override ;

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

    void AddChiReaction(BNV_Chi*) ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_BNV_B_Chi_Combo_H*/
