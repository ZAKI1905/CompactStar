/*
  Decay_Analysis class
*/

#include <Zaki/Vector/DataSet.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/BNV/Decay_Analysis.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "CompactStar/Core/Pulsar.hpp"

using namespace CompactStar ;
//==============================================================
//                        Decay_Analysis class
//==============================================================
// Constructor
Decay_Analysis::Decay_Analysis() 
  : 
  m_chi(0.8233*Zaki::Physics::NEUTRON_M_FM),
  neutron("neutron", "10", Zaki::Physics::NEUTRON_M_FM, 3.826),
  lambda("lambda", "100", Zaki::Physics::LAMBDA_ZERO_M_FM, -1.22) 
{ 
  // bar_list = { {"neutron", "10", Zaki::Physics::NEUTRON_M_FM}, 
  //              {"lambda", "100", 1115*Zaki::Physics::MEV_2_INV_FM},
  //              {"sigma-", "110", 1197*Zaki::Physics::MEV_2_INV_FM} } ;
}

//--------------------------------------------------------------
// Destructor
Decay_Analysis::~Decay_Analysis() { }
//--------------------------------------------------------------
double Decay_Analysis::G_func(const double& x) const
{
  double tmp = x*sqrt(x*x + 1);

  tmp -= log( x + sqrt(x*x + 1) ) ;

  return tmp ;
}

// //--------------------------------------------------------------
// // The Lambda decay rate per volume is in units of 
// //     fm^-3 s^-1 
// // m_lm is in units of fm^{-1}
// // n_lam is in units of fm^{-3}
// double Decay_Analysis::n_dot(const double& m_lam, 
//                              const double& n_lam) const
// {
//   if ( m_lam <= m_chi )
//     return 0 ;

//   double p_lam = pow(3*M_PI*M_PI*n_lam, 1.0/3.0) ;

//   double out =  g_lam*g_lam ;

//   out *= Zaki::Physics::Q_E*Zaki::Physics::Q_E ;
//   out *= eps_lam*eps_lam ;
//   out /= 512*M_PI*M_PI*M_PI ;
//   out *= (m_lam + m_chi)*(m_lam + m_chi) ;
//   out *= (1 - pow(m_chi / m_lam, 2) ) ;
//   out *= G_func(p_lam/m_lam) ;

//   // Converting from [ fm^-4 ] to [fm^-3 . s^-1]
//   out *= Zaki::Physics::LIGHT_C_M_S / 1e15 ;

//   return out ;
// }

//--------------------------------------------------------------
double Decay_Analysis::AmplitudeSqrd(const double& m_lam_bar) const
{
  return 0 ;
} 

//--------------------------------------------------------------
double Decay_Analysis::PhaseSpace(const double& m_lam_bar, 
                                  const double& n_lam) const
{
  return 0 ;
}

//--------------------------------------------------------------
// The Lambda decay rate per volume is in units of 
//     fm^-3 s^-1 
// m_lam is in units of fm^{-1}
// n_lam is in units of fm^{-3}
double Decay_Analysis::n_dot( const Baryon& B,
                              const double& m_eff,
                              const double& V,
                              const double& n) const
{
  if ( m_eff + V <= m_chi )
    return 0 ;

  // double out =  AmplitudeSqrd(m_lam_bar) * PhaseSpace(m_lam_bar, n_lam) ;

  double out =  B.g*B.g ;
  out *= Zaki::Physics::Q_E*Zaki::Physics::Q_E ;
  out *= eps_lam*eps_lam ;

  out *= m_eff + V - m_chi ;
  out *= pow(m_eff + V + m_chi, 3) ;
  out /= m_eff * pow(m_eff + V, 4) ;
  out /= 128*M_PI ;

  out *= n ;

  // Converting from [ fm^-4 ] to [fm^-3 . s^-1]
  // WARNING: Check this !!!!!!!!!!!!!!!!!
  out *= Zaki::Physics::LIGHT_C_M_S * 1e15 ;
  // out *= Zaki::Physics::LIGHT_C_M_S / 1e15 ;

  return out ;
}

//--------------------------------------------------------------
// Outputs the decay rate in [s^-1]
// eps_lam input is in fm^-1
double Decay_Analysis::Baryon_Vacuum_Decay(const Baryon& B, 
const double& in_eps_lam) 
  const
{
  // double m_B_0 = Zaki::Physics::LAMBDA_ZERO_M_FM ;
  double m_B_0 = B.mass ;

  if(m_B_0 < m_chi)
    return 0 ;

  double out =  B.g * B.g ;
  out *= Zaki::Physics::Q_E*Zaki::Physics::Q_E ;
  out /= 128*M_PI ;
  out *= m_B_0 * in_eps_lam * in_eps_lam ;
  out /= pow(m_B_0 - m_chi, 2) ;
  out *= pow(1 - pow(m_chi / m_B_0, 2), 3) ;

  // Converting the rate from [fm^-1] to [s^-1]
  // WARNING: Check this !!!!!!!!!!!!!!!!!
  // out *= Zaki::Physics::LIGHT_C_M_S / 1e15 ;
  out *= Zaki::Physics::LIGHT_C_M_S * 1e15 ;

  return out ;
}

//--------------------------------------------------------------
// Imports the effective mass set
void Decay_Analysis::ImportEffMass(const std::string& f_name, 
                          const Zaki::String::Directory& in_dir
                                  ) 
{
  // Zaki::Vector::DataSet m_eff_ds(wrk_dir + in_dir, f_name) ;

  m_eff_ds.SetWrkDir(wrk_dir) ;
  m_eff_ds.Import(in_dir + f_name) ;

  m_eff_ds[0].label = "$n (fm^{-3})$" ;
  m_eff_ds[1].label = "$m_n$" ;
  m_eff_ds[2].label = "$m_{\\Lambda}$" ;
  m_eff_ds[3].label = "$m_{\\Sigma^{-}}$" ;

  m_eff_ds.Plot(0, {1,2,3}, "M_eff.pdf") ;

  // Converting MeV to fm^-1
  for (size_t i = 1; i < m_eff_ds.Dim().size(); i++)
  {
    m_eff_ds[i] *= Zaki::Physics::MEV_2_INV_FM ;
  }

  m_eff_ds.Interpolate(0, {1, 2}) ;
}

//--------------------------------------------------------------
// Imports the vector self-energy set
void Decay_Analysis::ImportVSelfEnergy(const std::string& f_name, 
                          const Zaki::String::Directory& in_dir
                                  ) 
{
  // Zaki::Vector::DataSet m_eff_ds(wrk_dir + in_dir, f_name) ;

  V_self_E_ds.SetWrkDir(wrk_dir) ;
  V_self_E_ds.Import(in_dir + f_name) ;

  V_self_E_ds[0].label = "$n (fm^{-3})$" ;
  V_self_E_ds[1].label = "$V_n$" ;
  V_self_E_ds[2].label = "$V_{\\Lambda}$" ;
  V_self_E_ds[3].label = "$V_{\\Sigma^{-}}$" ;

  V_self_E_ds.Plot(0, {1,2,3}, "V_Self_E.pdf") ;

  // Converting MeV to fm^-1
  for (size_t i = 1; i < V_self_E_ds.Dim().size(); i++)
  {
    V_self_E_ds[i] *= Zaki::Physics::MEV_2_INV_FM ;
  }
  

  V_self_E_ds.Interpolate(0, {1, 2}) ;
}

//--------------------------------------------------------------
// Plots [ V_self_E_ds + m_eff_ds ]
void Decay_Analysis::PlotRestEnergy()
{
  Zaki::Vector::DataSet plt_ds({m_eff_ds[0], 
      (m_eff_ds[1] + V_self_E_ds[1])/Zaki::Physics::MEV_2_INV_FM, 
      (m_eff_ds[2] + V_self_E_ds[2])/Zaki::Physics::MEV_2_INV_FM, 
      (m_eff_ds[3] + V_self_E_ds[3])/Zaki::Physics::MEV_2_INV_FM}) ;

  plt_ds[1].label = "$E_n$" ;
  plt_ds[2].label = "$E_{\\Lambda}$" ;
  plt_ds[3].label = "$E_{\\Sigma^{-}}$" ;

  std::cout << " * --------------------------------------------------- * \n" ;
  std::cout << "\t Rest Energy at n_min = " << m_eff_ds[0][0] << " : \n" ;
  std::cout << "\t E(Neutron) :  " << plt_ds[1][0] << " MeV\n" ;
  std::cout << "\t E(Lambda)  :  " << plt_ds[2][0] << " MeV\n" ;
  std::cout << "\t E(Sigma-)  :  " << plt_ds[3][0] << " MeV\n\n" ;

  std::cout << " * --------------------------------------------------- * \n" ;
  std::cout << "\t Effective Masses at n_min = " << m_eff_ds[0][0] << " : \n" ;
  std::cout << "\t m*(Neutron) :  " << m_eff_ds[1][0]/Zaki::Physics::MEV_2_INV_FM << " MeV\n" ;
  std::cout << "\t m*(Lambda)  :  " << m_eff_ds[2][0]/Zaki::Physics::MEV_2_INV_FM << " MeV\n" ;
  std::cout << "\t m*(Sigma-)  :  " << m_eff_ds[3][0]/Zaki::Physics::MEV_2_INV_FM << " MeV\n" ;
  std::cout << " * --------------------------------------------------- * \n" ;

  plt_ds.SetWrkDir(wrk_dir) ;
  plt_ds.Plot(0, {1,2,3}, "Rest_E.pdf") ;

}

//--------------------------------------------------------------
// Attaches the pulsar
void Decay_Analysis::AttachPulsar(Pulsar* puls) 
{ 
  pulsar = puls ;
  Zaki::Vector::DataColumn r_set = pulsar->profile[0] ;
  Zaki::Vector::DataColumn M_r = pulsar->profile[1] ;
  Zaki::Vector::DataColumn nu_r = pulsar->profile[6];

  Zaki::Vector::DataColumn m_chi_set("$m_{\\chi} (m_n)$", 
                                    {0.05, 0.1, 0.2, 0.3, 0.4, 
                                    0.5, 0.6, 
                                    // 0.65, 
                                    0.7, 
                                    // 0.72, 
                                    // 0.74, 
                                    0.75, 
                                    // 0.76, 0.77, 
                                    // 0.78, 0.79, 0.795, 
                                    0.80,
                                    0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.89,
                                    0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99,
                                    1.0, 1.01, 1.02, 1.03, 1.04, 1.05, 1.07, 1.1, 1.11, 1.12, 
                                    1.13, 1.14, 1.15, 1.17, 1.1867, 1.2, 1.22, 
                                    1.25, 1.27, 1.3, 1.32, 1.35, 1.37, 1.39, 1.4, 1.42, 1.44
                                    }) ;

  Zaki::Vector::DataColumn lam_vac_br_set ;
  lam_vac_br_set.label = "Br$\\left(\\Lambda \\to \\chi + \\gamma\\right)$" ;
  
  Zaki::Vector::DataColumn eps_lam_lim_set ;
  eps_lam_lim_set.label = "$\\varepsilon_{\\Lambda} (GeV)$" ;

  Zaki::Vector::DataColumn neu_vac_br_set ;
  neu_vac_br_set.label = "Br$\\left(n \\to \\chi + \\gamma\\right)$" ;
  
  Zaki::Vector::DataColumn eps_n_lim_set ;
  eps_n_lim_set.label = "$\\varepsilon_{n} (GeV)$" ;

  for (size_t i = 0; i < m_chi_set.Size(); i++)
  {
    m_chi = m_chi_set[i]*Zaki::Physics::NEUTRON_M_FM ;

    // Unit = [fm^-3 . s^-1]
    Zaki::Vector::DataColumn n_lam_dot_dc ;
    n_lam_dot_dc.label = "n_lam_dot" ;

    Zaki::Vector::DataColumn n_neu_dot_dc ;
    n_neu_dot_dc.label = "n_lam_dot" ;

    double tmp_m_eff, tmp_V, tmp_den ;
    for (size_t i = 0; i < pulsar->profile[5].Size(); i++)
    {
      tmp_m_eff = m_eff_ds.Evaluate(2, pulsar->profile[5][i]) ;
      tmp_V = V_self_E_ds.Evaluate(2, pulsar->profile[5][i]) ;
      tmp_den = pulsar->profile[5][i] * pulsar->profile[14][i] ;

      n_lam_dot_dc.vals.emplace_back( n_dot(lambda, tmp_m_eff, tmp_V, tmp_den) ) ;

      tmp_m_eff = m_eff_ds.Evaluate(1, pulsar->profile[5][i]) ;
      tmp_V = V_self_E_ds.Evaluate(1, pulsar->profile[5][i]) ;
      tmp_den = pulsar->profile[5][i] * pulsar->profile[13][i] ;

      n_neu_dot_dc.vals.emplace_back( n_dot(neutron, tmp_m_eff, tmp_V, tmp_den) ) ;

    }

    Zaki::Vector::DataColumn integ_b_dot = 4*M_PI*r_set.pow(2) ;
    integ_b_dot *= (1. - 2*M_r / r_set ).pow(-0.5) ;
    integ_b_dot *= exp( nu_r ) ;

    n_lam_dot_dc *= integ_b_dot;
    n_neu_dot_dc *= integ_b_dot;

    Zaki::Vector::DataSet integrand({r_set, n_neu_dot_dc, n_lam_dot_dc}) ;
    integrand.Interpolate(0, {1, 2}) ;

    // Unit : [ km^3 . fm^-3 . s^-1 ]
    double neutron_dot_result = integrand.Integrate(1, {r_set[0], r_set[-1]}) ;
    double lambda_dot_result = integrand.Integrate(2, {r_set[0], r_set[-1]}) ;

    // Unit : [ s^-1 ]
    neutron_dot_result *= 1e54 ;
    lambda_dot_result *= 1e54 ;

    // Unit : [ yr^-1 ]
    neutron_dot_result *= 3600*24*365 ;
    lambda_dot_result *= 3600*24*365 ;

    double neutron_b_dot_over_b = -neutron_dot_result/pulsar->seq_point.b ;
    double lambda_b_dot_over_b = -lambda_dot_result/pulsar->seq_point.b ;

    double spin_down_limit = -abs(pulsar->GetBNVSpinDownLimit()) ;


    double scale_eps_neu = sqrt( spin_down_limit / neutron_b_dot_over_b ) ;
    double scale_eps_lam = sqrt( spin_down_limit / lambda_b_dot_over_b ) ;


    double eps_neu_limit_GeV = scale_eps_neu * eps_lam / Zaki::Physics::MEV_2_INV_FM ;
    double eps_lam_limit_GeV = scale_eps_lam * eps_lam / Zaki::Physics::MEV_2_INV_FM ;

    // Convert MeV to GeV
    eps_neu_limit_GeV *= 1e-3 ;
    eps_lam_limit_GeV *= 1e-3 ;

    double neutron_bnv_vac = Baryon_Vacuum_Decay(neutron, scale_eps_neu * eps_lam ) ;
    double neutron_bnv_vac_Br = Baryon_Vacuum_Decay(neutron, scale_eps_neu * eps_lam ) * 879.6  ;

    double lambda_bnv_vac = Baryon_Vacuum_Decay(lambda, scale_eps_lam * eps_lam ) ;
    double lambda_bnv_vac_Br = Baryon_Vacuum_Decay(lambda, scale_eps_lam * eps_lam ) * 2.632e-10  ;

    std::cout << "\n\t m_chi  = "<< m_chi/Zaki::Physics::MEV_2_INV_FM << " [MeV].\n" ;
    std::cout << "\t neutron_dot  = "<< -neutron_dot_result << " [yr^-1].\n" ;
    std::cout << "\t neutron_dot/B  = "<<  neutron_b_dot_over_b << " [yr^-1]. \n" ;
    std::cout << "\t Spin-Down (B_dot/B)  = "<<  spin_down_limit << " [yr^-1]. \n" ;
    std::cout << "\t Scale Lambda  = "<<  scale_eps_lam << ". \n" ;
    std::cout << "\t Limit on Eps_Lam  = "<<  eps_lam_limit_GeV << " [ GeV ]. \n" ;
    std::cout << "\t Lambda Exotic Vacuum Decay = "
              <<  lambda_bnv_vac << " [ s^-1 ]. \n" ;
    std::cout << "\t Lambda Exotic Vacuum Br = "
              <<  lambda_bnv_vac_Br << " . \n" ;

    neu_vac_br_set.vals.emplace_back(neutron_bnv_vac_Br) ;
    lam_vac_br_set.vals.emplace_back(lambda_bnv_vac_Br) ;

    eps_lam_lim_set.vals.emplace_back(eps_lam_limit_GeV) ;
    eps_n_lim_set.vals.emplace_back(eps_neu_limit_GeV) ;
  }

  Zaki::Vector::DataSet lam_plt_pt({m_chi_set, lam_vac_br_set, eps_lam_lim_set}) ;
  lam_plt_pt.SetWrkDir(wrk_dir) ;
  lam_plt_pt.Plot(0, 1, "Lambda_Br.pdf", 
  "Limit on Br$\\left(\\Lambda \\to \\chi + \\gamma\\right)$ \n from PSR J0348+0432 spin-down.") ;

  lam_plt_pt.SemiLogYPlot(0, 2, "Eps_Lambda.pdf", "Limit on $\\varepsilon_{\\Lambda} (GeV)$ \n from PSR J0348+0432 spin-down.") ;
  lam_plt_pt.Export("Lambda Plot Points.tsv") ;

  Zaki::Vector::DataSet neu_plt_pt({m_chi_set, neu_vac_br_set, eps_n_lim_set}) ;
  neu_plt_pt.SetWrkDir(wrk_dir) ;
  neu_plt_pt.Plot(0, 1, "n_Br.pdf", 
  "Limit on Br$\\left(n \\to \\chi + \\gamma\\right)$ \n from PSR J0348+0432 spin-down.") ;

  neu_plt_pt.SemiLogYPlot(0, 2, "Eps_n.pdf", "Limit on $\\varepsilon_{n} (GeV)$ \n from PSR J0348+0432 spin-down.") ;
  neu_plt_pt.Export("Neutron Plot Points.tsv") ;
}

//--------------------------------------------------------------
// Analysis during the sequence loop
void Decay_Analysis::Analyze(NStar* in_star) 
{
  Zaki::Vector::DataColumn r_set = *in_star->GetRadius();
  Zaki::Vector::DataColumn M_r = *in_star->GetMass();
  Zaki::Vector::DataColumn nu_r = *in_star->GetNu();
  
  // This is for making evaluations faster
  // we still need to multiply this by the 
  // individual species baryon density
  // A.K.A "b_den", in the loop below.
  Zaki::Vector::DataColumn integ_fr_0 = 4*M_PI*r_set.pow(2) ;
  integ_fr_0 *= (1. - 2*M_r / r_set ).pow(-0.5) ;

  Zaki::Vector::DataColumn integ_b_dot_0 = integ_fr_0*exp( nu_r ) ;


  // bnv_rate_seq.emplace_back(in_star->GetSequence().b, 
  //                           in_star->GetSequence().m, 
  //                           in_star->GetSequence().ec ) ;

  // for (size_t i = 0; i < bar_list.size(); i++)
  // {
  //   Zaki::Vector::DataColumn b_den = *in_star->GetRho_i(bar_list[i].label)
  //                                   * 
  //                                   *in_star->GetRho() ;

  //   Zaki::Vector::DataColumn integ_fr = integ_fr_0*b_den ;   
  //   Zaki::Vector::DataColumn integ_b_dot = integ_b_dot_0*b_den ;                             

  //   Zaki::Vector::DataSet integrand({r_set, integ_fr, integ_b_dot}) ;

  //   integrand.Interpolate(0, {1, 2}) ;
  //   double fr_result = integrand.Integrate(1, {r_set[0], r_set[-1]}) ;
  //   double b_dot_result = integrand.Integrate(2, {r_set[0], r_set[-1]}) ;

  //   bnv_rate_seq[bnv_rate_seq.size()-1].Append( 
  //                         { bar_list[i],
  //                           fr_result* 1e54,
  //                           b_dot_result* 1e54}) ;
  // }
  
} 

//--------------------------------------------------------------
// Saves the results
void Decay_Analysis::Export(const Zaki::String::Directory& in_dir)
{
  // Zaki::File::VecSaver vec_saver(in_dir + "/BNV_rates_"+label+".tsv") ;

  // std::string bnv_header = "" ;
  
  // for (size_t i = 0; i < bar_list.size(); i++)
  // {
  //   char tmp_bnv_header[100] ;
  //   sprintf(tmp_bnv_header, "%-16s\t %-16s\t ", 
  //               (bar_list[i].label + "(fr)").c_str(), 
  //               (bar_list[i].label +"(B_dot)").c_str() ) ;
  //   bnv_header += std::string(tmp_bnv_header) ;            
  // }

  // char tmp_bnv_header[100] ;
  // sprintf(tmp_bnv_header, "%-16s\t %-16s\t %-16s", 
  //               "B", "M", "ec" ) ;

  // bnv_header += std::string(tmp_bnv_header) ;            

  // vec_saver.SetHeader(bnv_header) ;
  // vec_saver.Export1D(bnv_rate_seq) ;
}

//--------------------------------------------------------------
// void Decay_Analysis::Evolve(const Zaki::String::Directory& in_dir) 
// {
//   double Gamma_BNV = 1e-10 ; // per year

//   Zaki::Vector::DataSet bnv_data(in_dir.ThisFileDir(), in_dir.ThisFile().Str()) ;
//   Zaki::Vector::DataColumn n_dot  = bnv_data[1]*Gamma_BNV ;
//   Zaki::Vector::DataColumn lam_dot = bnv_data[3]*Gamma_BNV*1e2 ;
//   Zaki::Vector::DataColumn sig_dot = bnv_data[5]*Gamma_BNV*1e2 ;

//   Zaki::Vector::DataColumn B      = bnv_data[-3] ;


//   // ---------------------------------------------
//   //              Neutron
//   // ---------------------------------------------
//   Zaki::Vector::DataColumn t_n_set ;
//   Zaki::Vector::DataColumn eps_n_set ;
//   Zaki::Vector::DataColumn M_n_set ;


//   t_n_set.Reserve(B.Size()) ;
//   t_n_set.label = "T [n]" ;
//   eps_n_set.label = "eps" ;
//   M_n_set.label = "M" ;

//   double t_n = 0 ; 

//   t_n_set.vals.emplace_back(t_n) ;
//   eps_n_set.vals.emplace_back(bnv_data[-1][-1]) ;
//   M_n_set.vals.emplace_back(bnv_data[-2][-1]) ;

//   for (size_t i = B.Size()-1; i > 0; i--)
//   // for (size_t i = 1; i < 60; i++) // Tau_age if M = 1.99
//   {
//     double delta_b = B[i] - B[i-1] ;

//     double n_dot_av = (n_dot[i] + n_dot[i-1] ) /2 ;  

//     if (!n_dot_av)
//       break ;

//     t_n += delta_b / n_dot_av ;

//     t_n_set.vals.emplace_back(t_n) ;
//     eps_n_set.vals.emplace_back(bnv_data[-1][i]) ;
//     M_n_set.vals.emplace_back(bnv_data[-2][i]) ;
//     // std::cout << t_n << "\n" ;
//   }
      
//   Zaki::Vector::DataSet n_out({t_n_set, eps_n_set, M_n_set}) ;
//   n_out.SetWrkDir(in_dir.ThisFileDir()) ;

//   Zaki::Vector::DataSet::PlotParam plt_par ;
//   plt_par.SetXAxis({1e5, 1e13 }) ;
//   plt_par.SetYAxis({1, 2.1}) ;

//   n_out.SetPlotPars(plt_par) ;

//   n_out.SemiLogXPlot(0, 2, "n_evolution.pdf") ;
//   n_out.Export("n_evolution.tsv") ;
//   // ---------------------------------------------

//   // ---------------------------------------------
//   //              Lambda
//   // ---------------------------------------------
//   Zaki::Vector::DataColumn M_lam_set ;
//   Zaki::Vector::DataColumn eps_lam_set ;
//   Zaki::Vector::DataColumn t_lam_set ;

//   t_lam_set.label = "T [Lambda]" ;
//   eps_lam_set.label = "eps" ;
//   M_lam_set.label = "M" ;
  
//   double t_lam = 0 ; 

//   t_lam_set.vals.emplace_back(t_lam) ;
//   eps_lam_set.vals.emplace_back(bnv_data[-1][-1]) ;
//   M_lam_set.vals.emplace_back(bnv_data[-2][-1]) ;

//   for (size_t i = B.Size()-1; i > 0; i--)
//   // for (size_t i = 1; i < 60; i++) // Tau_age if M = 1.99
//   {
//     double delta_b = B[i] - B[i-1] ;

//     double lam_dot_av = (lam_dot[i] + lam_dot[i-1] ) /2 ;  

//     if (!lam_dot_av)
//     {
//       t_lam_set.vals.emplace_back(t_lam*100) ;
//       eps_lam_set.vals.emplace_back(bnv_data[-1][i+1]) ;
//       M_lam_set.vals.emplace_back(bnv_data[-2][i+1]) ;
//       break ;
//     }
//     t_lam += delta_b / lam_dot_av ;
    
//     t_lam_set.vals.emplace_back(t_lam) ;
//     eps_lam_set.vals.emplace_back(bnv_data[-1][i]) ;
//     M_lam_set.vals.emplace_back(bnv_data[-2][i]) ;
//   }
      
//   Zaki::Vector::DataSet lam_out({t_lam_set, eps_lam_set, M_lam_set}) ;
//   lam_out.SetWrkDir(in_dir.ThisFileDir()) ;
//   lam_out.SetPlotPars(plt_par) ;
//   lam_out.SemiLogXPlot(0, 2, "Lambda_evolution.pdf") ;
//   lam_out.Export("Lambda_evolution.tsv") ;
//   // ---------------------------------------------

//   // ---------------------------------------------
//   //              Sigma-
//   // ---------------------------------------------
//   Zaki::Vector::DataColumn M_sig_set ;
//   Zaki::Vector::DataColumn eps_sig_set ;
//   Zaki::Vector::DataColumn t_sig_set ;

//   t_sig_set.label = "T [Sigma-]" ;
//   eps_sig_set.label = "eps" ;
//   M_sig_set.label = "M" ;

//   double t_sig = 0 ; 

//   t_sig_set.vals.emplace_back(t_sig) ;
//   eps_sig_set.vals.emplace_back(bnv_data[-1][-1]) ;
//   M_sig_set.vals.emplace_back(bnv_data[-2][-1]) ;

//   for (size_t i = B.Size()-1; i > 0; i--)
//   // for (size_t i = 1; i < 60; i++) // Tau_age if M = 1.99
//   {
//     double delta_b = B[i] - B[i-1] ;

//     double sig_dot_av = (sig_dot[i] + sig_dot[i-1] ) /2 ;  

//     if (!sig_dot_av)
//     {
//       t_sig_set.vals.emplace_back(t_sig*100) ;
//       eps_sig_set.vals.emplace_back(bnv_data[-1][i+1]) ;
//       M_sig_set.vals.emplace_back(bnv_data[-2][i+1]) ;
//       break ;
//     }
//     t_sig += delta_b / sig_dot_av ;
    
//     t_sig_set.vals.emplace_back(t_sig) ;
//     eps_sig_set.vals.emplace_back(bnv_data[-1][i]) ;
//     M_sig_set.vals.emplace_back(bnv_data[-2][i]) ;
//   }
      
//   Zaki::Vector::DataSet sig_out({t_sig_set, eps_sig_set, M_sig_set}) ;
//   sig_out.SetWrkDir(in_dir.ThisFileDir()) ;
//   sig_out.SetPlotPars(plt_par) ;
//   sig_out.SemiLogXPlot(0, 2, "Sigma_evolution.pdf") ;
//   sig_out.Export("sigma_evolution.tsv") ;
//   // ---------------------------------------------

// }
//--------------------------------------------------------------

//==============================================================
