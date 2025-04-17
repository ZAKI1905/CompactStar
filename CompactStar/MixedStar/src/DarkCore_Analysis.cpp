/*
  DarkCore_Analysis Abstract class
*/

#include <Zaki/Vector/DataSet.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/MixedStar/DarkCore_Analysis.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

using namespace CompactStar ;
//==============================================================
//                        DarkCore_Analysis class
//==============================================================
// Constructor
DarkCore_Analysis::DarkCore_Analysis() 
{ 
  // bnv_rate_seq = { {"neutron", "10", Zaki::Physics::NEUTRON_M_FM} } ;
  bar_list = { {"neutron", "10", Zaki::Physics::NEUTRON_M_FM}, 
               {"lambda", "100", 1115*Zaki::Physics::MEV_2_INV_FM},
               {"sigma-", "110", 1197*Zaki::Physics::MEV_2_INV_FM} } ;

  m_chi = 0.8*Zaki::Physics::NEUTRON_M_FM ;


  for (size_t i = 0; i < bar_list.size(); i++)
  {
    // double m_b = bnv_rate_seq[i].bar.mass ;
    // double critical_rho_val = pow(m_b - m_chi*m_chi/m_b, 3)
    //                         / ( 24 * M_PI*M_PI ) ;
    // bnv_rate_seq[i].SetCriticalChiDen(pow(m_b - m_chi*m_chi/m_b, 3)
    //                         / ( 24 * M_PI*M_PI )) ;

    double m_b = bar_list[i].mass ;
    bar_list[i].SetCriticalChiDen(pow(m_b - m_chi*m_chi/m_b, 3)
                            / ( 24 * M_PI*M_PI )) ;
  }

  // std::cout << " n:      Crit_chi_den " << bar_list[0].crit_chi_den << "\n" ;
  // std::cout << " Lambda: Crit_chi_den " << bar_list[1].crit_chi_den << "\n" ;
  // std::cout << " Sigma-: Crit_chi_den " << bar_list[2].crit_chi_den << "\n" ;
  neutron_out.Reserve(9, 10) ;
  lambda_out.Reserve(9, 10) ;
  sigmam_out.Reserve(9, 10) ;

  // neutron_out[0].label = "Channel" ;
  neutron_out[0].label = "M(0)" ;
  neutron_out[1].label = "t_0 [ch]" ;
  neutron_out[2].label = "t_0" ;
  neutron_out[3].label = "age [ch]" ;
  neutron_out[4].label = "age" ;
  neutron_out[5].label = "death [ch]" ;
  neutron_out[6].label = "death" ;
  neutron_out[7].label = "age/death [ch]" ;
  neutron_out[8].label = "age/death" ;

  // lambda_out[0].label = "Channel" ;
  lambda_out[0].label = "M(0)" ;
  lambda_out[1].label = "t_0 [ch]" ;
  lambda_out[2].label = "t_0" ;
  lambda_out[3].label = "age [ch]" ;
  lambda_out[4].label = "age" ;
  lambda_out[5].label = "death [ch]" ;
  lambda_out[6].label = "death" ;
  lambda_out[7].label = "age/death [ch]" ;
  lambda_out[8].label = "age/death" ;

  // sigmam_out[0].label = "Channel" ;
  sigmam_out[0].label = "M(0)" ;
  sigmam_out[1].label = "t_0 [ch]" ;
  sigmam_out[2].label = "t_0" ;
  sigmam_out[3].label = "age [ch]" ;
  sigmam_out[4].label = "age" ;
  sigmam_out[5].label = "death [ch]" ;
  sigmam_out[6].label = "death" ;
  sigmam_out[7].label = "age/death [ch]" ;
  sigmam_out[8].label = "age/death" ;
}
//--------------------------------------------------------------
// Destructor
DarkCore_Analysis::~DarkCore_Analysis() { }
//--------------------------------------------------------------
// Analysis during the sequence loop
void DarkCore_Analysis::Analyze(MixedStar* in_star) 
{

  // double m_n = Zaki::Physics::NEUTRON_M_FM ;

  // double m_chi = 0.8*Zaki::Physics::NEUTRON_M_FM ;

  Zaki::Vector::DataColumn r_set = *in_star->GetRadius();
  Zaki::Vector::DataColumn M_r = *in_star->GetMass_Total();
  Zaki::Vector::DataColumn nu_r = *in_star->GetNu();
  
  // This is for making evaluations faster
  // we still need to multiply this by the 
  // individual species baryon density
  // A.K.A "b_den", in the loop below.
  Zaki::Vector::DataColumn integ_fr_0 = 4*M_PI*r_set.pow(2) ;
  integ_fr_0 *= (1. - 2*M_r / r_set ).pow(-0.5) ;

  Zaki::Vector::DataColumn integ_b_dot_0 = integ_fr_0*exp( nu_r ) ;
  // integ_b_dot *= exp( nu_r ) ;


  bnv_rate_seq.emplace_back(in_star->GetSequence().v.b +
                              in_star->GetSequence().d.b, 
                              in_star->GetSequence().v.b, 
                              in_star->GetSequence().d.b, 
                              in_star->GetSequence().v.m +
                              in_star->GetSequence().d.m ) ;

  for (size_t i = 0; i < bar_list.size(); i++)
  {
    // std::cout << "\n -> crit. den. = " << bar_list[i].crit_chi_den << "\n";
    // std::cout << "\n -> Rho dark min = " << in_star->GetRho_Dark()->Min() << "\n";
    // std::cout << "\n -> Rho dark max = " << in_star->GetRho_Dark()->Max() << "\n";

    int critical_idx = in_star->GetRho_Dark()->GetClosestIdx(bar_list[i].crit_chi_den) ;
    Zaki::Vector::DataColumn b_den = *in_star->GetRho_i_Visible(bar_list[i].label)
                                    * 
                                    *in_star->GetRho_Visible() ;

    Zaki::Vector::DataColumn integ_fr = integ_fr_0*b_den ;   
    Zaki::Vector::DataColumn integ_b_dot = integ_b_dot_0*b_den ;                             

    Zaki::Vector::DataSet integrand({r_set, integ_fr, integ_b_dot}) ;

    integrand.Interpolate(0, {1, 2}) ;
    double fr_result = integrand.Integrate(1, {r_set[0], r_set[-1]}) ;
    double b_dot_result = integrand.Integrate(2, {r_set[0], r_set[-1]}) ;
  
    double fr_result_choked = integrand.Integrate(1, {r_set[critical_idx], r_set[-1]}) ;
    double b_dot_result_choked = integrand.Integrate(2, {r_set[critical_idx], r_set[-1]}) ;

    bnv_rate_seq[bnv_rate_seq.size()-1].Append( 
                          { bar_list[i],
                            fr_result* 1e54,
                            b_dot_result* 1e54, 
                            fr_result_choked* 1e54,
                            b_dot_result_choked* 1e54}) ;
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
void DarkCore_Analysis::Export(const Zaki::String::Directory& in_dir)
{
  Zaki::File::VecSaver vec_saver(in_dir + "/BNV_rates_"+label+".tsv") ;

  std::string bnv_header = "" ;
  
  for (size_t i = 0; i < bar_list.size(); i++)
  {
    char tmp_bnv_header[100] ;
    snprintf(tmp_bnv_header,  sizeof(tmp_bnv_header),  
                "%-16s\t %-16s\t %-16s\t %-16s\t ", 
                (bar_list[i].label + "(fr)").c_str(), 
                (bar_list[i].label +"(B_dot)").c_str(), 
                (bar_list[i].label +"(fr_ch)").c_str(),
                (bar_list[i].label +"(B_dot_ch)").c_str() ) ;
    bnv_header += std::string(tmp_bnv_header) ;            
  }

  char tmp_bnv_header[100] ;
  snprintf(tmp_bnv_header, sizeof(tmp_bnv_header), "%-16s\t %-16s\t %-16s\t %-16s", 
                "B", "B_v", "B_d", "M" ) ;

  bnv_header += std::string(tmp_bnv_header) ;            

  vec_saver.SetHeader(bnv_header) ;
  vec_saver.Export1D(bnv_rate_seq) ;
}

//--------------------------------------------------------------
void DarkCore_Analysis::FindLimits(const double& in_mass, 
                              const Zaki::String::Directory& in_dir) 
{

  Zaki::Vector::DataSet bnv_data(in_dir.ThisFileDir(), in_dir.ThisFile().Str()) ;
  Zaki::Vector::DataColumn M          = bnv_data[-1] ;
  Zaki::Vector::DataColumn B_d        = bnv_data[-2] ;
  Zaki::Vector::DataColumn B_v        = bnv_data[-3] ;
  Zaki::Vector::DataColumn B          = bnv_data[-4] ;

  // std::vector<DarkCore_Output> bnv_out(bar_list.size()) ;

  double Gamma_BNV = 1e-10 ; // per year

  int m_idx = M.GetClosestIdx(in_mass) ;
  if(m_idx == -1)
    m_idx = M.Size()-1 ;

  std::cout << "\n\t -> m_idx = '" << m_idx 
            << "' out of '" << M.Size() << "'.\n" ;

  // double init_M = M[0] ;
  double init_B_d = B_d[0] ;

  DarkCore_Output bnv_out ;
  std::cout << bnv_out.header() << "\n" ;
  for (size_t b_i = 0; b_i < bar_list.size(); b_i++)
  {
    // std::cout << "\n -------------------------------------------------\n" ;
    // std::cout << "\n\t\t\t" << bar_list[b_i].name <<"\n" ;
    // std::cout << "\n -------------------------------------------------\n" ;
    bnv_out.bar = bar_list[b_i] ;
    bnv_out.m_0 = M[0] ;


    // double time = 0 ;
    // double time_chok = 0 ;
    // double age = 0 ;
    // double age_chok = 0 ;

    Zaki::Vector::DataColumn B_dot      = bnv_data[1 + 4*b_i]*Gamma_BNV ;
    Zaki::Vector::DataColumn B_dot_chok = bnv_data[3 + 4*b_i]*Gamma_BNV ;

    // time += init_B_d/B_dot[0]  ;
    // time_chok += init_B_d/B_dot_chok[0]  ;

    // std::cout << "delta t(0) = " << time  << " yr.\n";
    // std::cout << "delta t(0) [choked] = " << time_chok  << " yr.\n";
    bnv_out.t_0 = {init_B_d/B_dot[0], init_B_d/B_dot_chok[0]}  ;
    bnv_out.death = bnv_out.t_0 ;

    std::cout.precision(8) ;
    for (size_t i = 1; i < B_v.Size(); i++)
    // for (size_t i = 1; i < 60; i++) // Tau_age if M = 1.99
    {
      double delta_b = B_v[i] - B_v[i-1] ;
      if (delta_b >= 0)
      {
        std::cout << "\n" << i << ") " << B_v[i] << "\n" ;
        continue;
      }
      
      double b_dot_rate = -(B_dot[i] + B_dot[i-1])/2 ;
      double b_dot_chok_rate = -(B_dot_chok[i] + B_dot_chok[i-1])/2 ;

      bnv_out.death += {delta_b / b_dot_rate, delta_b / b_dot_chok_rate} ;
      // time_chok += delta_b / b_dot_chok_rate ;

      if(i == m_idx || (m_idx==0 && i==1) )
      { 
        // age = time ;
        // age_chok = time_chok ;
        bnv_out.age = bnv_out.death ;
      }

    }

    // std::cout << " Death = " << time  << " yr.\n";
    // std::cout << " Death [choked] = " << time_chok  << " yr.\n";

    // std::cout << "\n Gam_BNV  = " << Gamma_BNV / (1e10 / time) << "\n" ;
    // std::cout << " Gam_BNV [choked] = " << Gamma_BNV / (1e10 / time_chok) << "\n" ;

    // std::cout << "\n Gam_BNV * Age = " << Gamma_BNV / (1 / age) << "\n" ;
    // std::cout << " Gam_BNV [choked] * Age = " << Gamma_BNV / (1 / age_chok) << "\n" ;

    // std::cout << "\n Gam_BNV * Death = " << Gamma_BNV / (1 / time) << "\n" ;
    // std::cout << " Gam_BNV [choked] * Death = " << Gamma_BNV / (1 / time_chok) << "\n" ;
    
    // std::cout << "\n Age/Death = " << age / time << "\n" ;
    // std::cout << " Age/Death [choked] = " << age_chok / time_chok << "\n" ;
    std::cout << bnv_out.Str() << "\n" ;

    if(bnv_out.bar.label == "10") 
      neutron_out.AppendRow(bnv_out.GetVals()) ;
    else if(bnv_out.bar.label == "100") 
      lambda_out.AppendRow(bnv_out.GetVals()) ;
    else if(bnv_out.bar.label == "110")
      sigmam_out.AppendRow(bnv_out.GetVals()) ;

  }
}

//--------------------------------------------------------------
// Saves the BNV results
void DarkCore_Analysis::ExportBNV(const Zaki::String::Directory& in_dir)
{
  // ---------------------------------------------
  //              Neutron
  // ---------------------------------------------
  neutron_out.SetWrkDir(in_dir + "/BNV_tau") ;
  neutron_out.Export("BNV_tau_neutron.tsv") ;

  neutron_out.Plot(0, 7, "BNV_neutron_age_death_ratio.pdf") ;
  neutron_out.Plot(0, {3,5}, "BNV_neutron_age_death.pdf", "Age & Death") ;

  // Zaki::Vector::DataColumn Gamm_BNV =  1e-10*neutron_out[3]/1e7 ;

  Zaki::Vector::DataSet n_gamma({neutron_out[0], 1e-10*neutron_out[3]/2.6e9,  1e-10*neutron_out[5]/1e10, neutron_out[-2]}) ;
  n_gamma[1].label = "Gamma_age" ;
  n_gamma[2].label = "Gamma_death" ;
  n_gamma.SetWrkDir(in_dir + "/BNV_tau") ;
  n_gamma.Plot(0, 1, "BNV_neutron_gamma_age.pdf", "Age = 2.6 Gyr") ;
  n_gamma.Plot(0, 2, "BNV_neutron_gamma_death.pdf", "Death = 10 Gyr") ;
  n_gamma.Plot(-1, 2, "BNV_neutron_gamma_death_ratio.pdf", "Death = 10 Gyr") ;

  n_gamma.SemiLogYPlot(0, {1,2}, "BNV_neutron_gamma.pdf", "Age = 2.6 Gyr, Death = 10 Gyr") ;
  // ---------------------------------------------

  // ---------------------------------------------
  //              Lambda
  // ---------------------------------------------
  lambda_out.SetWrkDir(in_dir + "/BNV_tau") ;
  lambda_out.Export("BNV_tau_lambda.tsv") ;

  lambda_out.Plot(0, 7, "BNV_lambda_age_death_ratio.pdf") ;
  lambda_out.Plot(0, {3,5}, "BNV_lambda_age_death.pdf", "Age & Death") ;

  Zaki::Vector::DataSet lambda_gamma({lambda_out[0], 1e-10*lambda_out[3]/2.6e9,  1e-10*lambda_out[5]/1e10, lambda_out[-2]}) ;
  lambda_gamma[1].label = "Gamma_age" ;
  lambda_gamma[2].label = "Gamma_death" ;
  lambda_gamma.SetWrkDir(in_dir + "/BNV_tau") ;
  lambda_gamma.Plot(0, 1, "BNV_lambda_gamma_age.pdf", "Age = 2.6 Gyr") ;
  lambda_gamma.Plot(0, 2, "BNV_lambda_gamma_death.pdf", "Death = 10 Gyr") ;
  lambda_gamma.Plot(-1, 2, "BNV_lambda_gamma_death_ratio.pdf", "Death = 10 Gyr") ;

  lambda_gamma.SemiLogYPlot(0, {1,2}, "BNV_lambda_gamma.pdf", "Age = 2.6 Gyr, Death = 10 Gyr") ;
  // ---------------------------------------------

  // ---------------------------------------------
  //              Sigma-
  // ---------------------------------------------
  sigmam_out.SetWrkDir(in_dir + "/BNV_tau") ;
  sigmam_out.Export("BNV_tau_sigmam.tsv") ;

  sigmam_out.Plot(0, 7, "BNV_sigmam_age_death_ratio.pdf") ;
  sigmam_out.Plot(0, {3,5}, "BNV_sigmam_age_death.pdf", "Age & Death") ;

  Zaki::Vector::DataSet sigmam_gamma({sigmam_out[0], 1e-10*sigmam_out[3]/2.6e9,  1e-10*sigmam_out[5]/1e10, sigmam_out[-2]}) ;
  sigmam_gamma[1].label = "Gamma_age" ;
  sigmam_gamma[2].label = "Gamma_death" ;
  sigmam_gamma.SetWrkDir(in_dir + "/BNV_tau") ;
  sigmam_gamma.Plot(0, 1, "BNV_sigma-_gamma_age.pdf", "Age = 2.6 Gyr") ;
  sigmam_gamma.Plot(0, 2, "BNV_sigma-_gamma_death.pdf", "Death = 10 Gyr") ;
  sigmam_gamma.Plot(-1, 2, "BNV_sigma-_gamma_death_ratio.pdf", "Death = 10 Gyr") ;

  sigmam_gamma.SemiLogYPlot(0, {1,2}, "BNV_sigma-_gamma.pdf", "Age = 2.6 Gyr, Death = 10 Gyr") ;
  // ---------------------------------------------
  //                Combinations
  // ---------------------------------------------
  // std::cout << "neutron: [" << neutron_out[-2][0] << ", "<< neutron_out[-2][-1] << "]\n" ;
  // std::cout << "Lambda : [" << lambda_out[-2][0] << ", "<< lambda_out[-2][-1] << "]\n" ;
  // std::cout << "Sigma- : [" << sigmam_out[-2][0] << ", "<< sigmam_out[-2][-1] << "]\n" ;
  // ---------------------------------------------

  // for (size_t i = 0; i < bar_list.size(); i++)
  // {
  //   Zaki::File::VecSaver vec_saver(in_dir + "/BNV_tau/BNV_tau_"+bar_list[i].name+".tsv") ;

  //   if(bar_list[i].label == "10") 
  //   {
  //     vec_saver.SetHeader(neutron_out[0].header()) ;
  //     vec_saver.Export1D(neutron_out) ;
  //   }
  //   else if(bar_list[i].label == "100") 
  //   {
  //     vec_saver.SetHeader(lambda_out[0].header()) ;
  //     vec_saver.Export1D(lambda_out) ;
  //   }
  //   else if(bar_list[i].label == "110")
  //   {
  //     vec_saver.SetHeader(sigmam_out[0].header()) ;
  //     vec_saver.Export1D(sigmam_out) ;
  //   }

  //   // char tmp_bnv_header[100] ;
  //   // sprintf(tmp_bnv_header, "%-16s\t %-16s\t %-16s\t %-16s\t ", 
  //   //             (bar_list[i].label + "(fr)").c_str(), 
  //   //             (bar_list[i].label +"(B_dot)").c_str(), 
  //   //             (bar_list[i].label +"(fr_ch)").c_str(),
  //   //             (bar_list[i].label +"(B_dot_ch)").c_str() ) ;
  //   // bnv_header += std::string(tmp_bnv_header) ;            
  // }

}
//--------------------------------------------------------------

//==============================================================
