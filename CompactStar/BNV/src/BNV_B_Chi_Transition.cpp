/*
  BNV_B_Chi_Transition class
*/

// // Creating directories
// #include <sys/stat.h>
// #include <filesystem>

#include <gsl/gsl_integration.h>

#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Vector/DataSet.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/BNV/BNV_B_Chi_Transition.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
// #include "CompactStar/Core/Pulsar.hpp"

using namespace CompactStar ;
//==============================================================
//                        BNV_B_Chi_Transition class
//==============================================================
// Constructor
BNV_B_Chi_Transition::BNV_B_Chi_Transition() 
  : 
  BNV_Chi({"B_Chi_Transition", "B \\to \\chi"})
{ 
  // Z_LOG_WARNING("Constructor called for BNV_B_Chi_Transition.") ;
  // Z_LOG_WARNING("Process '" + process.name + "' created.") ;
}

//--------------------------------------------------------------
// Destructor
BNV_B_Chi_Transition::~BNV_B_Chi_Transition() { }
//--------------------------------------------------------------
BNV_Chi::Process BNV_B_Chi_Transition::GetSpecificProcess(const Baryon& B) const
{
  return {B.short_name + "_Chi_Transition", B.TeX_name + " \\to \\chi"} ;
}
//--------------------------------------------------------------
// This function returns the upper limit on the vector 
// self-energy, from condition (b). 
// The arguments are (m_chi, m_B^*, n [fm^-3] ).
// Zaki::Vector::DataColumn BNV_B_Chi_Transition::SigmaPlus(
//                         const double& m_chi, 
//                         const Baryon& B)
// { 

//   Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI*n_B[B.label]).pow(2.0/3.0) ;

//   kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

//   Zaki::Vector::DataColumn  out  = m_B_ds[B.label]*m_B_ds[B.label] 
//                                   + m_chi*m_chi + 2*kF_2 ;
  
//   out += 2*((m_B_ds[B.label]*m_B_ds[B.label] + kF_2)*(kF_2 + m_chi*m_chi )).sqrt() ;
//   out = out.sqrt() ;
//   out.label = "$\\Sigma^+_0 (MeV)$" ;
//   return  out ;
// }

//--------------------------------------------------------------
// This function returns the lower limit on the vector 
// self-energy, from condition (a).
// The arguments are (m_chi, m_B^*, n [fm^-3] ).
Zaki::Vector::DataColumn BNV_B_Chi_Transition::SigmaMinus(
                        const double& in_m_chi, 
                        const Baryon& B)
{ 

  Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI*n_B[B.label]).pow(2.0/3.0) ;

  kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

  Zaki::Vector::DataColumn  out  = m_B_ds[B.label]*m_B_ds[B.label] 
                                  + in_m_chi*in_m_chi + 2*kF_2 ;
  
  out -= 2*((m_B_ds[B.label]*m_B_ds[B.label] + kF_2)*(kF_2 + in_m_chi*in_m_chi )).sqrt() ;
  out = out.sqrt() ;
  out.label = "$\\Sigma^-_0 (MeV)$" ;
  return  out ;
}

//--------------------------------------------------------------
// This function returns the lower limit on the vector 
// self-energy, from condition (a).
// The arguments are (m_chi, m_B^*, n [fm^-3] ).
double BNV_B_Chi_Transition::SigmaMinus(
                        const double& m_chi, 
                        const double& m_B, 
                        const double& n_B)
{ 

  double kF_2 = pow(3*M_PI*M_PI*n_B, 2.0/3.0) ;

  kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

  double  out  = m_B*m_B + m_chi*m_chi + 2*kF_2 ;
  
  out -= 2*sqrt((m_B*m_B + kF_2)*(kF_2 + m_chi*m_chi )) ;
  out = sqrt(out) ;

  return  out ;
}

//--------------------------------------------------------------
// This function returns the upper limit on the vector 
// self-energy, from condition (b). 
// The arguments are (m_chi, m_B^*, n [fm^-3] ).
// double BNV_B_Chi_Transition::SigmaPlus(
//                   const double& m_chi, 
//                   const double& m_B, 
//                   const double& n_B)
// { 

//   double kF_2 = pow(3*M_PI*M_PI*n_B, 2.0/3.0) ;

//   kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

//   double  out  = m_B*m_B + m_chi*m_chi + 2*kF_2 ;
  
//   out += 2*sqrt((m_B*m_B + kF_2)*(kF_2 + m_chi*m_chi )) ;
//   out = sqrt(out) ;

//   return  out ;
// }

//--------------------------------------------------------------
// This function returns the B -> chi decay rate 
//  as a function of radius in units of s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Chi_Transition::Rate_vs_R(
                                const double& m_chi, 
                                const Baryon& B, 
                                const bool& gen_plots)
{ 
  // double eps = 1e-10 ;

  Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;
  Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
  Zaki::Vector::DataColumn n_B_r = n_r * pulsar.GetProfile()->operator[](B.label) ;

  Zaki::Vector::DataSet ds_m_B_n({n_B["n_tot"], eos.GetMeff(B.label)}) ;
  Zaki::Vector::DataSet sig0_n({n_B["n_tot"], eos.GetVeff(B.label) }) ;
  ds_m_B_n.Interpolate(0, 1) ;
  sig0_n.Interpolate(0, 1) ;
  Zaki::Vector::DataColumn m_B_r = ds_m_B_n.Evaluate(1, n_r) ;
  Zaki::Vector::DataColumn sig0_r = sig0_n.Evaluate(1, n_r) ;

  Zaki::Vector::DataSet  rate ;
  rate.Reserve(2, r.Size() ) ;
  for (size_t i = 0; i < r.Size() ; i++)
  {
    
    double m_B = m_B_r[i] ;
    double n_B = n_B_r[i] ;
    double sig0 = sig0_r[i] ;

    // Old condition
    // if ( abs(sig0) < SigmaMinus(m_chi, m_B, n_B) 
    //       || 
    //      abs(sig0) > SigmaPlus(m_chi, m_B, n_B)  )
    // New Condition
    if ( sig0 < SigmaMinus(m_chi, m_B, n_B))
    {
      rate.AppendRow({r[i], 0}) ;
      continue ;
    }
    
    // Old condition
    // if ( abs(sig0) > abs(m_B - m_chi) 
    //     &&
    //     abs(sig0) < m_B + m_chi )
    // New Condition
    if ( sig0 > m_chi - m_B )
    {
      rate.AppendRow({r[i], 0}) ;
      continue ;
    }

    double p = pow(sig0 - m_B, 2) - m_chi * m_chi ;
          p *= pow(sig0 + m_B, 2) - m_chi * m_chi  ;
          p  = sqrt(p) ;
          p /= 2 * abs(sig0) ;
  
    // double amp_sqrd = abs(pow(sig0, 2) - pow(m_B, 2) + m_chi * m_chi ) ;
    // amp_sqrd *= sig0 / ( 2* abs(sig0) ) ;
    // amp_sqrd += m_chi*m_chi + m_B * m_chi ;
    // amp_sqrd *= 2 * eps * eps ;

    double amp_sqrd = pow( m_chi + m_B, 2 ) - pow(sig0, 2) ;
          amp_sqrd *= default_eps * default_eps ;

    double rate_val = p * amp_sqrd / (2*M_PI*abs(sig0)) ;
    rate.AppendRow({r[i], rate_val}) ;
  }

  // Converting "MeV^4" into "MeV/fm^3"
  rate[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "MeV/fm^3" into "s^-1/fm^3"
  rate[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;

  if(gen_plots)
  {
    hidden_Plot_Rate_vs_R(B, m_chi, &rate) ;
  }

  return  rate ;
}

//--------------------------------------------------------------
// // This function returns the B -> chi decay rate 
// //  as a function of radius in units of s^-1/fm^3
// void BNV_B_Chi_Transition::Rate_vs_R(const double& m_chi)         
// {
//   Rate_vs_R(m_chi, neutron, true) ;
//   Rate_vs_R(m_chi, lambda, true) ;
// }

//--------------------------------------------------------------
/// Returns the B -> chi decay rate 
///  as a function of density in units of s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Chi_Transition::Rate_vs_Density(
                                      const double& m_chi, 
                                      const Baryon& B,
                                const bool& gen_plots)
{
//  double eps = 1e-10 ;

  // Finding the total and individual baryon density
  // Zaki::Vector::DataColumn n_tot = eos.GetEOS(2) ;
  // Zaki::Vector::DataColumn n_B = n_tot * eos.GetEOS(B.label) ;
  // Zaki::Vector::DataColumn m_B = eos.GetMeff(B.label) ;
  // Zaki::Vector::DataColumn sig0 = eos.GetVeff(B.label) ;

  // Zaki::Vector::DataColumn dc_Sigma_Plus = SigmaPlus(m_chi, B) ;
  Zaki::Vector::DataColumn dc_Sigma_Minus = SigmaMinus(m_chi, B) ; 

  Zaki::Vector::DataSet  rate ;
  rate.Reserve(2, n_B["n_tot"].Size() ) ;
  for (size_t i = 0; i < n_B["n_tot"].Size() ; i++)
  {
    // if ( abs(sig0_ds[B.label][i]) < dc_Sigma_Minus[i] 
    //       || 
    //      abs(sig0_ds[B.label][i]) > dc_Sigma_Plus[i]  )
    if ( sig0_ds[B.label][i] < dc_Sigma_Minus[i] )
    {
      rate.AppendRow({n_B["n_tot"][i], 0}) ;
      continue ;
    }
    
    // if ( abs(sig0_ds[B.label][i]) > abs(m_B_ds[B.label][i] - m_chi) 
    //     &&
    //     abs(sig0_ds[B.label][i]) < m_B_ds[B.label][i] + m_chi )
    if ( sig0_ds[B.label][i] > m_chi - m_B_ds[B.label][i] )
    {
      rate.AppendRow({n_B["n_tot"][i], 0}) ;
      continue ;
    }

    double p = pow(sig0_ds[B.label][i] - m_B_ds[B.label][i], 2) - m_chi * m_chi ;
          p *= pow(sig0_ds[B.label][i] + m_B_ds[B.label][i], 2) - m_chi * m_chi  ;
          p  = sqrt(p) ;
          p /= 2 * abs(sig0_ds[B.label][i]) ;
  
    // double amp_sqrd = abs(pow(sig0_ds[B.label][i], 2) - pow(m_B_ds[B.label][i], 2) + m_chi * m_chi ) ;
    // amp_sqrd *= sig0_ds[B.label][i] / ( 2* abs(sig0_ds[B.label][i]) ) ;
    // amp_sqrd += m_chi*m_chi + m_B_ds[B.label][i] * m_chi ;
    // amp_sqrd *= 2 * eps * eps ;
    double amp_sqrd = pow( m_chi + m_B_ds[B.label][i], 2 ) - pow(sig0_ds[B.label][i], 2) ;
          amp_sqrd *= default_eps * default_eps ;

    double rate_val = p * amp_sqrd / (2*M_PI*abs(sig0_ds[B.label][i])) ;
    rate.AppendRow({n_B["n_tot"][i], rate_val}) ;
  }

  // Converting "MeV^4" into "MeV/fm^3"
  rate[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "MeV/fm^3" into "s^-1/fm^3"
  rate[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;

  if(gen_plots)
  {
    hidden_Plot_Rate_vs_Density(B, m_chi, &rate) ;
  }

  return  rate ;
}

// //--------------------------------------------------------------
// // This function returns the rate for
// //  B -> chi decay integrated over the radius 
// //  and in units of [ s^-1 ]
// double BNV_B_Chi_Transition::GetRateLim(const double& m_chi, 
//                           const Baryon& B)
// {
//   Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
//   Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
//   Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;
 
//  Zaki::Vector::DataColumn dc_B_chi_rate_r = Rate_vs_R(m_chi, B)[1] ;

//   // Converting "MeV^3" into "1/fm^3"
//   dc_B_chi_rate_r *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

//   // Converting "1/fm^3" into "1/km^3"
//   dc_B_chi_rate_r *= 1e+54 ; 

//   // Converting "MeV" into "s^-1"
//   dc_B_chi_rate_r *= 1e-3 * Zaki::Physics::GEV_2_S ; 

//   Zaki::Vector::DataColumn decay_integrand_dc = 4*M_PI*r*r*dc_B_chi_rate_r ;  
//                           decay_integrand_dc /= (1 - 2*M_r / r).sqrt() ; 
//                           decay_integrand_dc *= exp(nu_r) ; 
//     Zaki::Vector::DataSet decay_integral_ds({r, decay_integrand_dc}) ;

//   decay_integral_ds.Interpolate(0, 1) ;    

//   double result = decay_integral_ds.Integrate(1, {r[0], r[-1]}) ;
  
//   // Convert s^-1 into yr^-1
//   result                *= 3600*24*365 ;

//   std::cout << "\n ---------------------------------------- " ;
//   std::cout << "\n\t Full decay rate (B_dot/B) : " 
//             << result / pulsar.GetSeqPoint().b << " per year." ;
//   // std::cout << "\n ---------------------------------------- " ;

//   return  result / pulsar.GetSeqPoint().b ;
// }

// //--------------------------------------------------------------
// // This function returns the limit on eps from 
// //  B -> chi decay rate integrated over the radius 
// //  and in units of [ MeV ]
// double BNV_B_Chi_Transition::GetEpsLim(const double& m_chi, 
//                                         const Baryon& B)
// {
//   // Zaki::Vector::DataColumn n_r = puls.GetProfile()->operator[](5) ;
//   Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
//   Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
//   Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

//   // Zaki::Vector::DataSet ds_B_chi_rate_n = 
//   //   B_to_Chi_Rate_vs_Density(m_chi, B, eos, dir) ;
//   // ds_B_chi_rate_n.Interpolate(0, 1) ;
  
//   // Zaki::Vector::DataColumn dc_B_chi_rate_r = 
//   //   ds_B_chi_rate_n.Evaluate(1, n_r) ;
 
//  Zaki::Vector::DataColumn dc_B_chi_rate_r = Rate_vs_R(m_chi, B)[1] ;

//   // Converting "MeV^3" into "1/fm^3"
//   dc_B_chi_rate_r *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

//   // Converting "1/fm^3" into "1/km^3"
//   dc_B_chi_rate_r *= 1e+54 ; 

//   // Converting "MeV" into "s^-1"
//   dc_B_chi_rate_r *= 1e-3 * Zaki::Physics::GEV_2_S ; 

//   Zaki::Vector::DataColumn decay_integrand_dc = 4*M_PI*r*r*dc_B_chi_rate_r ;  
//                           decay_integrand_dc /= (1 - 2*M_r / r).sqrt() ; 
//                           decay_integrand_dc *= exp(nu_r) ; 
//     Zaki::Vector::DataSet decay_integral_ds({r, decay_integrand_dc}) ;

//   decay_integral_ds.Interpolate(0, 1) ;    

//   double result = decay_integral_ds.Integrate(1, {r[0], r[-1]}) ;
  
//   // Convert s^-1 into yr^-1
//   result *= 3600*24*365 ;

//   std::cout << "\n ---------------------------------------- " ;
//   std::cout << "\n\t Full decay rate (B_dot/B) : " 
//             << result / pulsar.GetSeqPoint().b << " per year." ;
//   // std::cout << "\n ---------------------------------------- " ;

//   return 1e-10*sqrt(7.83129*1e-10 / (result / pulsar.GetSeqPoint().b) ) ;
// }

//--------------------------------------------------------------
// Plots the transition conditions 
// for a given m_chi as a function of density
void BNV_B_Chi_Transition::PlotTransCond(const double& m_chi, 
                                            const Baryon& B)
{
  Zaki::Vector::DataColumn m_eff = m_B_ds[B.label] ;
  Zaki::Vector::DataColumn V_eff = sig0_ds[B.label] ;

  // Zaki::Vector::DataSet ds_Cond_a({n_B[0], V_eff.Abs(), (m_eff - m_chi).Abs(),
                                  // SigmaMinus(m_chi, B)}) ;

  Zaki::Vector::DataSet ds_Cond_a({n_B[0], V_eff, m_chi - m_eff,
                                  SigmaMinus(m_chi, B)}) ;
  ds_Cond_a.SetWrkDir(wrk_dir) ;

  char m_chi_ch[50] ;
  snprintf(m_chi_ch, 50, "%.0f", m_chi) ;
  std::string m_chi_str(m_chi_ch) ;

  Zaki::Vector::DataSet::PlotParam plt_par_0 ;
  plt_par_0.SetLegend({"upper left", 0.0, 1.0}) ;
  plt_par_0.SetXAxis({0, 3}) ;

  ds_Cond_a.SetPlotPars(plt_par_0) ;
  ds_Cond_a.Plot(0, {{1, "$\\Sigma^0$"}, {2, "$m_{\\chi} - m_B^*$"}, {3, "$\\Sigma^-$"}},
                        model + "/" + process.name + "/Transit_Lim_a_"
                        + model + "_" + B.short_name + "_" + m_chi_str + ".pdf", 
                      "Transition limits (a) \n for $" + B.TeX_name
                      + " \\to \\chi$ in "+ model + " EoS. ($m_{\\chi} = "
                      + m_chi_str + "$ MeV)") ;
  
  // Zaki::Vector::DataSet ds_Cond_b({n_B[0], V_eff.Abs(), m_eff + m_chi,
  //                                 SigmaPlus(m_chi, B)}) ;
  // ds_Cond_b.SetWrkDir(wrk_dir) ;
  // // plt_par_0.SetLegend({"upper left", 0.0, 1.0}) ;
  // ds_Cond_b.SetPlotPars(plt_par_0) ;
  // ds_Cond_b.Plot(0, {{1, "$\\Sigma^0$"}, {2, "$m_B + m_{\\chi}$"}, {3, "$\\Sigma^+$"}},
  //                       model + "/" + process.name + "/Transit_Lim_b_"
  //                       + model+"_" + B.short_name + "_" + m_chi_str + ".pdf", 
  //                     "Transition limits (b)\n for $" + B.TeX_name
  //                     + " \\to \\chi$ in "+ model + " EoS. ($m_{\\chi} = "
  //                     + m_chi_str + "$ MeV)") ;
}

//--------------------------------------------------------------
// Plots the transition conditions 
// for a given m_chi as a function of density
void BNV_B_Chi_Transition::PlotTransCond(const double& m_chi)
{
  PlotTransCond(m_chi, neutron) ;
  PlotTransCond(m_chi, lambda) ;
}

//--------------------------------------------------------------
// Plots the transition band as a function of m_chi                  
void BNV_B_Chi_Transition::PlotTransBand(const Baryon& B) 
{
  Zaki::Vector::DataSet ds_trans_range ;
  ds_trans_range.SetWrkDir(wrk_dir) ;
  ds_trans_range.Reserve(3, m_chi_vals.res) ;
  ds_trans_range[0].label = "$m_{\\chi}$" ;
  ds_trans_range[1].label = "$n_{min} (fm^{-3})$" ;
  ds_trans_range[2].label = "$n_{max} (fm^{-3})$" ;
  // ds_trans_range[3].label = "$n_{min} (b) (fm^{-3})$" ;
  // ds_trans_range[4].label = "$n_{max} (b) (fm^{-3})$" ;

  // ds_trans_range[5].label = "$M_{max} = 2.07$" ; // {1, 0.92, 2.07}, {3, 0.96, 2.004}
  // ds_trans_range[6].label = "$M = 2.05$" ; // at {1, 0.76 fm^{-3}}, {3, NaN}
  // ds_trans_range[7].label = "$M = 1.97$" ; // at {1, 0.61 fm^{-3}}, {3, 0.76 fm^{-3}}}

  Zaki::Vector::DataColumn m_eff = m_B_ds[B.label] ;
  Zaki::Vector::DataColumn V_eff = sig0_ds[B.label] ;

  for (size_t i = 0; i < m_chi_vals.res; i++)
  {
    double m_chi = m_chi_vals[i] ;

    // Zaki::Vector::DataColumn dc_Sigma_Plus = SigmaPlus(m_chi, B) ;
    Zaki::Vector::DataColumn dc_Sigma_Minus = SigmaMinus(m_chi, B) ;

    Zaki::Vector::DataSet ds_n_min({n_B[0], dc_Sigma_Minus - V_eff.Abs()}) ;
    double n_min = ds_n_min.Solve(0, 1) ;

    Zaki::Vector::DataSet ds_n_max({n_B[0], V_eff.Abs() - (m_chi-m_eff)}) ;
    double n_max = ds_n_max.Solve(0, 1) ;

    // Zaki::Vector::DataSet ds_n_b_max({n_B[0], dc_Sigma_Plus - V_eff.Abs()}) ;
    // double n_b_max = ds_n_b_max.Solve(0, 1) ;

    // Zaki::Vector::DataSet ds_n_b_min({n_B[0], m_eff + m_chi - V_eff.Abs()}) ;
    // double n_b_min = ds_n_b_min.Solve(0, 1) ;

    // ds_trans_range.AppendRow({m_chi, dc_n_tot[n_min_a_idx], dc_n_tot[n_max_a_idx], 
    //                           dc_n_tot[n_min_b_idx], dc_n_tot[n_max_b_idx], 0.92, 0.76, 0.61}) ;
    ds_trans_range.AppendRow({m_chi, n_min, n_max}) ;
  }
  // ds_trans_range.MakeSmooth(3) ;

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetXAxis({0, 2500}) ;
  plt_par.SetYAxis({0, 3}) ;
  plt_par.SetLegend({"lower right", 1.0, 0.0}) ;
  ds_trans_range.SetPlotPars(plt_par) ;

  ds_trans_range.Plot(0, {1, 2}, 
                        model + "/" + process.name + "/Transit_Den_Ranges_"
                        + model+"_" + B.short_name + ".pdf", 
                      "Transition Density Ranges\n for $" + B.TeX_name
                      + " \\to \\chi$ in "+ model + " EoS.") ;
}
//--------------------------------------------------------------
// Plots the transition band as a function of m_chi                  
void BNV_B_Chi_Transition::PlotTransBand() 
{
  PlotTransBand(neutron) ;
  PlotTransBand(lambda) ;
}
//--------------------------------------------------------------

//==============================================================
