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
 * @file BNV_Analysis.hpp
 *
 * @brief Class for analyzing baryon number violation in neutron stars
 *
 * @ingroup BNV
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@uky.edu
 *
*/

// Created on June 4, 2022
#ifndef CompactStar_BNV_Analysis_H
#define CompactStar_BNV_Analysis_H

// #include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/Analysis.hpp"


//==============================================================
namespace CompactStar
{

//==============================================================
struct BNV_Rate
{
  struct Baryon 
  {
    std::string name  = "";
    std::string label = "" ; 
    double mass = 0 ;

    Baryon() {}
    Baryon(const std::string& in_name, 
           const std::string& in_label, 
           const double& in_mass) 
           : name(in_name), 
             label(in_label),
             mass(in_mass)
            {}
  } ;

  const Baryon bar ;
   
  double fr, b_dot ;

  BNV_Rate( const Baryon& in_bar,
            const double& in_fr, 
            const double& in_b_dot)
        : bar(in_bar), fr(in_fr), b_dot(in_b_dot)
          {}
  
  std::string Str() const
  {
    std::stringstream ss;
    char tmp[200] ;
    snprintf(tmp, sizeof(tmp), "%-16.8e\t %-16.8e", 
            fr, b_dot) ;
    ss << tmp ;
  
    return ss.str() ;
  }
};
//==============================================================
// struct BNV_B
// {
//   std::vector<BNV_Rate> bnv_rates ;

//   double b_tot, b_vis, b_dar ;
// };

//==============================================================
struct BNV_Rate_Seq
{
  // struct Baryon 
  // {
  //   std::string name  = "";
  //   std::string label = "" ; 
  //   double mass = 0 ;
  // } ;

  // Baryon bar ;
  double b_tot, m_tot, ec ;
  std::vector<BNV_Rate> bnv_rates ;

  // std::vector<BNV_B> bnv_b ;

  // double crit_chi_den = -1 ;

  BNV_Rate_Seq(const double& in_b_tot,
               const double& in_m_tot, 
               const double& in_ec)
               : b_tot(in_b_tot), m_tot(in_m_tot), 
                 ec(in_ec)
  {}

  // BNV_Rate_Seq(const Baryon& in_bar)
  // {
  //   bar = in_bar ;
  // }

  void Append(const BNV_Rate& in_bnv_rate)
  {
    bnv_rates.emplace_back(in_bnv_rate) ;
  }

  // /// Sets the critical chi density that fully blocks the decay
  // void SetCriticalChiDen(const double&) ;

  std::string Str() const
  {
    std::stringstream ss;

    for (size_t i = 0; i < bnv_rates.size(); i++)
    {
      ss << bnv_rates[i].Str() << "\t ";
    }
    
    char tmp[200] ;
    snprintf(tmp, sizeof(tmp), "%-16.8e\t %-16.8e\t %-16.8e", 
                  b_tot, m_tot, ec) ;
    ss << tmp ;
  
    return ss.str() ;
  }
};

//==============================================================
struct BNV_Time
{
  double t = 0 ;

  BNV_Time() {}

  BNV_Time(const double& in_t)
    : t(in_t) {}

  BNV_Time operator/(const BNV_Time& rhs) const
  {
    return t/rhs.t ;
  }

  BNV_Time operator+(const BNV_Time& rhs) const
  {
    return t+rhs.t ;
  }

  BNV_Time& operator+=(const BNV_Time& rhs)
  {
    t    += rhs.t ;

    return *this ;
  }

  BNV_Time operator-(const BNV_Time& rhs) const
  {
    return t-rhs.t ;
  }

  std::string Str() const
  {
    char tmp[100] ;
    snprintf(tmp, sizeof(tmp), "%.6e", t) ;
    return std::string(tmp) ;
  }
};

//==============================================================
struct BNV_Output
{
  /// The decaying baryon
  BNV_Rate::Baryon bar ;

  /// Mass at t = t_0
  double m_0 ;

  /// Default constructor
  BNV_Output() {}

  // std::string Str() const
  // {
  //   char tmp[200] ;
  //   sprintf(tmp, "%-8s\t %.4f\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e\t %.6e", 
  //           bar.label.c_str(), m_0, t_0.t_ch, t_0.t,
  //           age.t_ch, age.t, death.t_ch, death.t, Age_over_Death().t_ch, Age_over_Death().t) ;
  //   return std::string(tmp) ;
  // }

  // std::string header() const
  // {
  //   char tmp[200] ;
  //   sprintf(tmp, "%-8s\t %-6s\t %-13s\t %-13s\t %-13s\t %-13s\t %-13s\t %-13s\t %-13s\t %-13s", 
  //           "Channel", "M(0)", "t_0", "t_0 [ch]", "age", 
  //           "age [ch]", "death", "death [ch]", "age/death", "age/death [ch]") ;
  //   return std::string(tmp) ;
  // }

  // std::vector<double> GetVals() const 
  // {
  //   return {m_0, t_0.t_ch, t_0.t,
  //           age.t_ch, age.t, death.t_ch, death.t,
  //           Age_over_Death().t_ch, Age_over_Death().t} ;
  // }
};

//==============================================================
class BNV_Analysis : public Analysis
{
  //--------------------------------------------------------------
  private:

  // std::vector<BNV_Rate> bnv_rates ;
  std::vector<BNV_Rate::Baryon> bar_list ; 
  std::vector<BNV_Rate_Seq> bnv_rate_seq ;

  // Zaki::Vector::DataSet neutron_out ;
  // Zaki::Vector::DataSet lambda_out ;
  // Zaki::Vector::DataSet sigmam_out ;

  //--------------------------------------------------------------
  public:

    BNV_Analysis() ;

    ~BNV_Analysis() ;

    /// Analysis during the sequence loop
    void Analyze(NStar* in_star) override ;
    void Analyze(MixedStar* in_star) override {}

    /// Saves the results
    void Export(const Zaki::String::Directory&) override ;

    void Evolve(const Zaki::String::Directory&) ;

    // Saves the BNV results
    // void ExportBNV(const Zaki::String::Directory&) ;
};

//==============================================================
} // CompactStar namespace
//==============================================================
#endif /*CompactStar_BNV_Analysis_H*/
