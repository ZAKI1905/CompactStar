/*
  BNV_Analysis Abstract class
*/

#include <Zaki/Vector/DataSet.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/BNV/BNV_Analysis.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

using namespace CompactStar ;
//==============================================================
//                        BNV_Analysis class
//==============================================================
// Constructor
BNV_Analysis::BNV_Analysis() 
{ 
  // bnv_rate_seq = { {"neutron", "10", Zaki::Physics::NEUTRON_M_FM} } ;
  bar_list = { {"neutron", "10", Zaki::Physics::NEUTRON_M_FM}, 
               {"lambda", "100", 1115*Zaki::Physics::MEV_2_INV_FM},
               {"sigma-", "110", 1197*Zaki::Physics::MEV_2_INV_FM} } ;

  // neutron_out.Reserve(4, 100) ;
  // lambda_out.Reserve(4, 100) ;
  // sigmam_out.Reserve(4, 100) ;

  // neutron_out[0].label = "t" ;
  // neutron_out[1].label = "ec" ;
  // neutron_out[2].label = "M" ;
  // neutron_out[3].label = "B" ;

  // lambda_out[0].label = "t" ;
  // lambda_out[1].label = "ec" ;
  // lambda_out[2].label = "M" ;
  // lambda_out[3].label = "B" ;

  // sigmam_out[0].label = "t" ;
  // sigmam_out[1].label = "ec" ;
  // sigmam_out[2].label = "M" ;
  // sigmam_out[3].label = "B" ;
}
//--------------------------------------------------------------
// Destructor
BNV_Analysis::~BNV_Analysis() { }
//--------------------------------------------------------------
// Analysis during the sequence loop
void BNV_Analysis::Analyze(NStar* in_star) 
{

  // double m_n = Zaki::Physics::NEUTRON_M_FM ;

  // double m_chi = 0.8*Zaki::Physics::NEUTRON_M_FM ;

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
  // integ_b_dot *= exp( nu_r ) ;


  bnv_rate_seq.emplace_back(in_star->GetSequence().b, 
                            in_star->GetSequence().m, 
                            in_star->GetSequence().ec ) ;

  for (size_t i = 0; i < bar_list.size(); i++)
  {
    // std::cout << "\n -> crit. den. = " << bar_list[i].crit_chi_den << "\n";
    // std::cout << "\n -> Rho dark min = " << in_star->GetRho_Dark()->Min() << "\n";
    // std::cout << "\n -> Rho dark max = " << in_star->GetRho_Dark()->Max() << "\n";

    Zaki::Vector::DataColumn b_den = *in_star->GetRho_i(bar_list[i].label)
                                    * 
                                    *in_star->GetRho() ;

    Zaki::Vector::DataColumn integ_fr = integ_fr_0*b_den ;   
    Zaki::Vector::DataColumn integ_b_dot = integ_b_dot_0*b_den ;                             

    Zaki::Vector::DataSet integrand({r_set, integ_fr, integ_b_dot}) ;

    integrand.Interpolate(0, {1, 2}) ;
    double fr_result = integrand.Integrate(1, {r_set[0], r_set[-1]}) ;
    double b_dot_result = integrand.Integrate(2, {r_set[0], r_set[-1]}) ;

    bnv_rate_seq[bnv_rate_seq.size()-1].Append( 
                          { bar_list[i],
                            fr_result* 1e54,
                            b_dot_result* 1e54}) ;
  }
  


  // double critical_rho_val = pow(m_n - m_chi*m_chi/m_n, 3)
  //                           / ( 24 * M_PI*M_PI ) ;
                                
  // int critical_idx = in_star->ds_dar[
  //                         in_star->rho_idx 
  //                         ].GetClosestIdx(critical_rho_val) ;

  // int critical_idx = in_star->GetRho_Dark()->GetClosestIdx(critical_rho_val) ;

  // std::cout << "Critical Rho = " << critical_rho_val 
  //           << ", Index = " << critical_idx 
  //           << ", Size = " 
  //           << in_star->ds_dar.Dim()[0] 
  //           << ", rho_d [v_idx-1] = " << in_star->ds_dar[
  //                         in_star->rho_idx 
  //                         ][critical_idx-1]
  //           << ", rho_d [v_idx] = " << in_star->ds_dar[
  //                         in_star->rho_idx 
  //                         ][critical_idx]
  //           << "\n\t r [v_idx-1] = " << in_star->ds_dar[
  //                         in_star->r_idx 
  //                         ][critical_idx-1] 
  //           << ", r [v_idx] = " << in_star->ds_dar[
  //                         in_star->r_idx 
  //                         ][critical_idx] << "\n" ;

// std::map<std::string, std::string> baryons = {
//   {"10", "Neutron", 6},
//   {"100", "Lambda", 7}, {"110", "Sigma-", 10}, {"111", "Sigma0"},
//   {"112", "Sigma+"}, {"120","Xi-"}, {"121", "Xi0"}
// } ;

  // Zaki::Vector::DataColumn b_den = in_star->ds_vis[
  //                                   in_star->rho_i_v_idx[6] 
  //                                   ] * 
  //                                   *in_star->GetRho_Visible() ;
  // Zaki::Vector::DataColumn b_den = *in_star->GetRho_i_Visible("10")
  //                                   * 
  //                                   *in_star->GetRho_Visible() ;

  // Zaki::Vector::DataColumn r_set = in_star->ds_vis[
  //                                   in_star->r_idx 
  //                                   ];
  
  // Zaki::Vector::DataColumn M_r = in_star->ds_vis[
  //                                   in_star->m_idx 
  //                                   ];

  // Zaki::Vector::DataColumn nu_r = in_star->ds_vis[
  //                                   in_star->nu_idx 
  //                                   ];

  // Zaki::Vector::DataColumn integ_fr = 4*M_PI*b_den ;
  // integ_fr *= r_set.pow(2) * (1. - 2*M_r / r_set ).pow(-0.5) ;

  // Zaki::Vector::DataColumn integ_b_dot = integ_fr ;
  // integ_b_dot *= exp( nu_r ) ;

  // Zaki::Vector::DataSet integrand({r_set, integ_fr, integ_b_dot}) ;

  // std::cout << "\n integ_fr.Size() = " << integ_fr.Size() ;
  // std::cout << "\n integ_b_dot.Size() = " << integ_b_dot.Size() ;

  // integrand.Interpolate(0, {1, 2}) ;
  // double fr_result = integrand.Integrate(1, {r_set[0], r_set[-1]}) ;
  // double b_dot_result = integrand.Integrate(2, {r_set[0], r_set[-1]}) ;
  
  // double fr_result_choked = integrand.Integrate(1, {r_set[critical_idx], r_set[-1]}) ;
  // double b_dot_result_choked = integrand.Integrate(2, {r_set[critical_idx], r_set[-1]}) ;

  // std::cout << "\nFraction = " << fr_result* 1e54 / in_star->sequence.v.b 
  //           << ", [B_dot/B] = " 
  //           << b_dot_result* 1e54 / in_star->sequence.v.b  
  //           <<"\n Fraction choked = "
  //           <<  fr_result_choked* 1e54 / in_star->sequence.v.b  
  //           << ", [B_dot/B] choked = "
  //           << b_dot_result_choked* 1e54 / in_star->sequence.v.b  << "\n\n" ; 
  // std::cout.precision(9) ;
  // std::cout << "\n B_dot = " 
  //           << b_dot_result* 1e54 
  //           << ", B_dot (choked) = "
  //           << b_dot_result_choked* 1e54
  //           <<", B_vis = " << in_star->GetSequence().v.b  
  //           << ", B_tot = " << in_star->GetSequence().v.b 
  //               + in_star->GetSequence().d.b  
  //           << "\n\n" ;
  // bnv_rates.emplace_back( b_dot_result* 1e54, 
  //                         b_dot_result_choked* 1e54,
  //                         in_star->GetSequence().v.b) ;
} 

//--------------------------------------------------------------
// Saves the results
void BNV_Analysis::Export(const Zaki::String::Directory& in_dir)
{
  Zaki::File::VecSaver vec_saver(in_dir + "/BNV_rates_"+label+".tsv") ;

  std::string bnv_header = "" ;
  
  for (size_t i = 0; i < bar_list.size(); i++)
  {
    char tmp_bnv_header[100] ;
    snprintf(tmp_bnv_header, sizeof(tmp_bnv_header), "%-16s\t %-16s\t ", 
                (bar_list[i].label + "(fr)").c_str(), 
                (bar_list[i].label +"(B_dot)").c_str() ) ;
    bnv_header += std::string(tmp_bnv_header) ;            
  }

  char tmp_bnv_header[100] ;
  snprintf(tmp_bnv_header, sizeof(tmp_bnv_header), "%-16s\t %-16s\t %-16s", 
                "B", "M", "ec" ) ;

  bnv_header += std::string(tmp_bnv_header) ;            

  vec_saver.SetHeader(bnv_header) ;
  vec_saver.Export1D(bnv_rate_seq) ;
}

//--------------------------------------------------------------
void BNV_Analysis::Evolve(const Zaki::String::Directory& in_dir) 
{
  double Gamma_BNV = 1e-10 ; // per year

  Zaki::Vector::DataSet bnv_data(in_dir.ThisFileDir(), in_dir.ThisFile().Str()) ;
  Zaki::Vector::DataColumn n_dot  = bnv_data[1]*Gamma_BNV ;
  Zaki::Vector::DataColumn lam_dot = bnv_data[3]*Gamma_BNV*1e2 ;
  Zaki::Vector::DataColumn sig_dot = bnv_data[5]*Gamma_BNV*1e2 ;

  Zaki::Vector::DataColumn B      = bnv_data[-3] ;


  // ---------------------------------------------
  //              Neutron
  // ---------------------------------------------
  Zaki::Vector::DataColumn t_n_set ;
  Zaki::Vector::DataColumn eps_n_set ;
  Zaki::Vector::DataColumn M_n_set ;


  t_n_set.Reserve(B.Size()) ;
  t_n_set.label = "T [n]" ;
  eps_n_set.label = "eps" ;
  M_n_set.label = "M" ;

  double t_n = 0 ; 

  t_n_set.vals.emplace_back(t_n) ;
  eps_n_set.vals.emplace_back(bnv_data[-1][-1]) ;
  M_n_set.vals.emplace_back(bnv_data[-2][-1]) ;

  for (size_t i = B.Size()-1; i > 0; i--)
  // for (size_t i = 1; i < 60; i++) // Tau_age if M = 1.99
  {
    double delta_b = B[i] - B[i-1] ;

    double n_dot_av = (n_dot[i] + n_dot[i-1] ) /2 ;  

    if (!n_dot_av)
      break ;

    t_n += delta_b / n_dot_av ;

    t_n_set.vals.emplace_back(t_n) ;
    eps_n_set.vals.emplace_back(bnv_data[-1][i]) ;
    M_n_set.vals.emplace_back(bnv_data[-2][i]) ;
    // std::cout << t_n << "\n" ;
  }
      
  Zaki::Vector::DataSet n_out({t_n_set, eps_n_set, M_n_set}) ;
  n_out.SetWrkDir(in_dir.ThisFileDir()) ;

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetXAxis({1e5, 1e13 }) ;
  plt_par.SetYAxis({1, 2.1}) ;

  n_out.SetPlotPars(plt_par) ;

  n_out.SemiLogXPlot(0, 2, "n_evolution.pdf") ;
  n_out.Export("n_evolution.tsv") ;
  // ---------------------------------------------

  // ---------------------------------------------
  //              Lambda
  // ---------------------------------------------
  Zaki::Vector::DataColumn M_lam_set ;
  Zaki::Vector::DataColumn eps_lam_set ;
  Zaki::Vector::DataColumn t_lam_set ;

  t_lam_set.label = "T [Lambda]" ;
  eps_lam_set.label = "eps" ;
  M_lam_set.label = "M" ;
  
  double t_lam = 0 ; 

  t_lam_set.vals.emplace_back(t_lam) ;
  eps_lam_set.vals.emplace_back(bnv_data[-1][-1]) ;
  M_lam_set.vals.emplace_back(bnv_data[-2][-1]) ;

  for (size_t i = B.Size()-1; i > 0; i--)
  // for (size_t i = 1; i < 60; i++) // Tau_age if M = 1.99
  {
    double delta_b = B[i] - B[i-1] ;

    double lam_dot_av = (lam_dot[i] + lam_dot[i-1] ) /2 ;  

    if (!lam_dot_av)
    {
      t_lam_set.vals.emplace_back(t_lam*100) ;
      eps_lam_set.vals.emplace_back(bnv_data[-1][i+1]) ;
      M_lam_set.vals.emplace_back(bnv_data[-2][i+1]) ;
      break ;
    }
    t_lam += delta_b / lam_dot_av ;
    
    t_lam_set.vals.emplace_back(t_lam) ;
    eps_lam_set.vals.emplace_back(bnv_data[-1][i]) ;
    M_lam_set.vals.emplace_back(bnv_data[-2][i]) ;
  }
      
  Zaki::Vector::DataSet lam_out({t_lam_set, eps_lam_set, M_lam_set}) ;
  lam_out.SetWrkDir(in_dir.ThisFileDir()) ;
  lam_out.SetPlotPars(plt_par) ;
  lam_out.SemiLogXPlot(0, 2, "Lambda_evolution.pdf") ;
  lam_out.Export("Lambda_evolution.tsv") ;
  // ---------------------------------------------

  // ---------------------------------------------
  //              Sigma-
  // ---------------------------------------------
  Zaki::Vector::DataColumn M_sig_set ;
  Zaki::Vector::DataColumn eps_sig_set ;
  Zaki::Vector::DataColumn t_sig_set ;

  t_sig_set.label = "T [Sigma-]" ;
  eps_sig_set.label = "eps" ;
  M_sig_set.label = "M" ;

  double t_sig = 0 ; 

  t_sig_set.vals.emplace_back(t_sig) ;
  eps_sig_set.vals.emplace_back(bnv_data[-1][-1]) ;
  M_sig_set.vals.emplace_back(bnv_data[-2][-1]) ;

  for (size_t i = B.Size()-1; i > 0; i--)
  // for (size_t i = 1; i < 60; i++) // Tau_age if M = 1.99
  {
    double delta_b = B[i] - B[i-1] ;

    double sig_dot_av = (sig_dot[i] + sig_dot[i-1] ) /2 ;  

    if (!sig_dot_av)
    {
      t_sig_set.vals.emplace_back(t_sig*100) ;
      eps_sig_set.vals.emplace_back(bnv_data[-1][i+1]) ;
      M_sig_set.vals.emplace_back(bnv_data[-2][i+1]) ;
      break ;
    }
    t_sig += delta_b / sig_dot_av ;
    
    t_sig_set.vals.emplace_back(t_sig) ;
    eps_sig_set.vals.emplace_back(bnv_data[-1][i]) ;
    M_sig_set.vals.emplace_back(bnv_data[-2][i]) ;
  }
      
  Zaki::Vector::DataSet sig_out({t_sig_set, eps_sig_set, M_sig_set}) ;
  sig_out.SetWrkDir(in_dir.ThisFileDir()) ;
  sig_out.SetPlotPars(plt_par) ;
  sig_out.SemiLogXPlot(0, 2, "Sigma_evolution.pdf") ;
  sig_out.Export("sigma_evolution.tsv") ;
  // ---------------------------------------------

}

//--------------------------------------------------------------
// Saves the BNV results
// void BNV_Analysis::ExportBNV(const Zaki::String::Directory& in_dir)
// {
  // ---------------------------------------------
  //              Neutron
  // ---------------------------------------------
//   neutron_out.SetWrkDir(in_dir + "/BNV_tau") ;
//   neutron_out.Export("BNV_tau_neutron.tsv") ;

//   neutron_out.Plot(0, 7, "BNV_neutron_age_death_ratio.pdf") ;
//   neutron_out.Plot(0, {3,5}, "BNV_neutron_age_death.pdf", "Age & Death") ;

//   // Zaki::Vector::DataColumn Gamm_BNV =  1e-10*neutron_out[3]/1e7 ;

//   Zaki::Vector::DataSet n_gamma({neutron_out[0], 1e-10*neutron_out[3]/2.6e9,  1e-10*neutron_out[5]/1e10, neutron_out[-2]}) ;
//   n_gamma[1].label = "Gamma_age" ;
//   n_gamma[2].label = "Gamma_death" ;
//   n_gamma.SetWrkDir(in_dir + "/BNV_tau") ;
//   n_gamma.Plot(0, 1, "BNV_neutron_gamma_age.pdf", "Age = 2.6 Gyr") ;
//   n_gamma.Plot(0, 2, "BNV_neutron_gamma_death.pdf", "Death = 10 Gyr") ;
//   n_gamma.Plot(-1, 2, "BNV_neutron_gamma_death_ratio.pdf", "Death = 10 Gyr") ;

//   n_gamma.SemiLogYPlot(0, {1,2}, "BNV_neutron_gamma.pdf", "Age = 2.6 Gyr, Death = 10 Gyr") ;
//   // ---------------------------------------------

//   // ---------------------------------------------
//   //              Lambda
//   // ---------------------------------------------
//   lambda_out.SetWrkDir(in_dir + "/BNV_tau") ;
//   lambda_out.Export("BNV_tau_lambda.tsv") ;

//   lambda_out.Plot(0, 7, "BNV_lambda_age_death_ratio.pdf") ;
//   lambda_out.Plot(0, {3,5}, "BNV_lambda_age_death.pdf", "Age & Death") ;

//   Zaki::Vector::DataSet lambda_gamma({lambda_out[0], 1e-10*lambda_out[3]/2.6e9,  1e-10*lambda_out[5]/1e10, lambda_out[-2]}) ;
//   lambda_gamma[1].label = "Gamma_age" ;
//   lambda_gamma[2].label = "Gamma_death" ;
//   lambda_gamma.SetWrkDir(in_dir + "/BNV_tau") ;
//   lambda_gamma.Plot(0, 1, "BNV_lambda_gamma_age.pdf", "Age = 2.6 Gyr") ;
//   lambda_gamma.Plot(0, 2, "BNV_lambda_gamma_death.pdf", "Death = 10 Gyr") ;
//   lambda_gamma.Plot(-1, 2, "BNV_lambda_gamma_death_ratio.pdf", "Death = 10 Gyr") ;

//   lambda_gamma.SemiLogYPlot(0, {1,2}, "BNV_lambda_gamma.pdf", "Age = 2.6 Gyr, Death = 10 Gyr") ;
//   // ---------------------------------------------

//   // ---------------------------------------------
//   //              Sigma-
//   // ---------------------------------------------
//   sigmam_out.SetWrkDir(in_dir + "/BNV_tau") ;
//   sigmam_out.Export("BNV_tau_sigmam.tsv") ;

//   sigmam_out.Plot(0, 7, "BNV_sigmam_age_death_ratio.pdf") ;
//   sigmam_out.Plot(0, {3,5}, "BNV_sigmam_age_death.pdf", "Age & Death") ;

//   Zaki::Vector::DataSet sigmam_gamma({sigmam_out[0], 1e-10*sigmam_out[3]/2.6e9,  1e-10*sigmam_out[5]/1e10, sigmam_out[-2]}) ;
//   sigmam_gamma[1].label = "Gamma_age" ;
//   sigmam_gamma[2].label = "Gamma_death" ;
//   sigmam_gamma.SetWrkDir(in_dir + "/BNV_tau") ;
//   sigmam_gamma.Plot(0, 1, "BNV_sigma-_gamma_age.pdf", "Age = 2.6 Gyr") ;
//   sigmam_gamma.Plot(0, 2, "BNV_sigma-_gamma_death.pdf", "Death = 10 Gyr") ;
//   sigmam_gamma.Plot(-1, 2, "BNV_sigma-_gamma_death_ratio.pdf", "Death = 10 Gyr") ;

//   sigmam_gamma.SemiLogYPlot(0, {1,2}, "BNV_sigma-_gamma.pdf", "Age = 2.6 Gyr, Death = 10 Gyr") ;
//   // ---------------------------------------------
//   //                Combinations
//   // ---------------------------------------------
//   // std::cout << "neutron: [" << neutron_out[-2][0] << ", "<< neutron_out[-2][-1] << "]\n" ;
//   // std::cout << "Lambda : [" << lambda_out[-2][0] << ", "<< lambda_out[-2][-1] << "]\n" ;
//   // std::cout << "Sigma- : [" << sigmam_out[-2][0] << ", "<< sigmam_out[-2][-1] << "]\n" ;
//   // ---------------------------------------------

//   // for (size_t i = 0; i < bar_list.size(); i++)
//   // {
//   //   Zaki::File::VecSaver vec_saver(in_dir + "/BNV_tau/BNV_tau_"+bar_list[i].name+".tsv") ;

//   //   if(bar_list[i].label == "10") 
//   //   {
//   //     vec_saver.SetHeader(neutron_out[0].header()) ;
//   //     vec_saver.Export1D(neutron_out) ;
//   //   }
//   //   else if(bar_list[i].label == "100") 
//   //   {
//   //     vec_saver.SetHeader(lambda_out[0].header()) ;
//   //     vec_saver.Export1D(lambda_out) ;
//   //   }
//   //   else if(bar_list[i].label == "110")
//   //   {
//   //     vec_saver.SetHeader(sigmam_out[0].header()) ;
//   //     vec_saver.Export1D(sigmam_out) ;
//   //   }

//   //   // char tmp_bnv_header[100] ;
//   //   // sprintf(tmp_bnv_header, "%-16s\t %-16s\t %-16s\t %-16s\t ", 
//   //   //             (bar_list[i].label + "(fr)").c_str(), 
//   //   //             (bar_list[i].label +"(B_dot)").c_str(), 
//   //   //             (bar_list[i].label +"(fr_ch)").c_str(),
//   //   //             (bar_list[i].label +"(B_dot_ch)").c_str() ) ;
//   //   // bnv_header += std::string(tmp_bnv_header) ;            
//   // }

// }
//--------------------------------------------------------------

//==============================================================
