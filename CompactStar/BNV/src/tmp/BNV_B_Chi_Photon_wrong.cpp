/*
  BNV_B_Chi_Photon class
*/

#include <gsl/gsl_integration.h>

#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Vector/DataSet.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/BNV/BNV_B_Chi_Photon.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "CompactStar/Core/Pulsar.hpp"

using namespace CompactStar ;

//==============================================================
bool MassCondition(const CompactStar::NStar& in_star)
{
  // return ( 2.0099 <= in_star.GetSequence().m 
  //           && 
  //         in_star.GetSequence().m <= 2.01001 )  ;

  return ( 2.0 <= in_star.GetSequence().m 
            && 
          in_star.GetSequence().m <= 2.011 )  ;
}
//==============================================================

//==============================================================
//                        BNV_B_Chi_Photon class
//==============================================================
// Constructor
BNV_B_Chi_Photon::BNV_B_Chi_Photon() 
  : 
  m_chi_old(0.8233*Zaki::Physics::NEUTRON_M_FM),
  neutron("neutron", "$n$", "10", Zaki::Physics::NEUTRON_M_FM,
           -3.82608545, 879.6),
  lambda("lambda", "$\\Lambda$", "100", 
          Zaki::Physics::LAMBDA_ZERO_M_FM, -1.22, 2.632e-10), 
  proton("proton", "$p^+$", "11", Zaki::Physics::PROTON_M_FM,
           5.5856946893, 1e+41)
{ }

//--------------------------------------------------------------
// Destructor
BNV_B_Chi_Photon::~BNV_B_Chi_Photon() { }
//--------------------------------------------------------------
// Sets the EoS model
void BNV_B_Chi_Photon::SetModel(const std::string& in_eos_model) 
{
  model = in_eos_model ;
}

//--------------------------------------------------------------
void BNV_B_Chi_Photon::GenSequence() const
{
  CompactStar::TOVSolver solver ;

  solver.SetWrkDir( wrk_dir ) ;
  solver.ImportEOS("EOS/CompOSE/"+ model +"/"+ model +".eos") ;

  double e_c_min = EOS_ds[0].Min()*1.01 ;
  double e_c_max = EOS_ds[0].Max()*0.99 ;

  double m_best = 0 ;
  double m_goal = 2.01 ;
  double precision = 0.0005 ;
  
  solver.Solve( {{e_c_min, e_c_max}, 50, "Log"}, 
                "BNV_2022/results/B_Chi_Photon_Decay/"+model,  
                "B_Chi_Photon_Decay") ;

  Zaki::Vector::DataSet tmp_seq(wrk_dir +
                  "/BNV_2022/results/B_Chi_Photon_Decay/"+ model,
                  "B_Chi_Photon_Decay_Sequence.tsv") ;
  
  e_c_min = tmp_seq[0][tmp_seq[1].MinIdx()] ;
  e_c_max = tmp_seq[0][tmp_seq[1].MaxIdx()] ;

  solver.ClearSequence() ;

  while ( abs(m_best - m_goal) > precision ) {

  // std::cout << "\n\n\t e_c_min =  " << e_c_min 
  //         << ", \t e_c_max = " << e_c_max 
  //         << "\n\n" << std::flush ;

  solver.Solve( {{e_c_min, e_c_max}, 10, "Log"}, 
                "BNV_2022/results/B_Chi_Photon_Decay/"+model,  
                "B_Chi_Photon_Decay") ;
  
  tmp_seq.Import("B_Chi_Photon_Decay_Sequence.tsv") ;


  int m_idx = tmp_seq[1].GetClosestIdx(m_goal) ;
  m_best = tmp_seq[1][m_idx] ;
  e_c_min = tmp_seq[0][m_idx-1] ;
  e_c_max = tmp_seq[0][m_idx+1] ;

  solver.ClearSequence() ;
  }

  solver.AddNCondition(MassCondition) ;

  solver.Solve( {{e_c_min, e_c_max}, 5, "Linear"}, 
                "BNV_2022/results/B_Chi_Photon_Decay/"+model,  
                "B_Chi_Photon_Decay") ;
}

//--------------------------------------------------------------
// Phase space integrand function (dim-less part)
// The total phase space integrand is the output of this function
// multiplied by:
//                  m^6 C^2 / (2 pi^3)
// Here sigma_0 = - Sigma_0 / m.
//
double BNV_B_Chi_Photon::PhaseSpace_Integrand(const double& x) 
{
  double mu_chi = phase_int_mu_chi ; 
  double sigma_0 = phase_int_sigma_0 ;

  if ( x <= ( mu_chi*mu_chi - 1 - sigma_0*sigma_0 ) / ( 2*sigma_0 )  ) 
  {
    return 0 ;
  }

  // Old result [ WRONG ]
  // double out = x * sqrt(x*x - 1) ;
  //       out *= pow(1 + sigma_0*sigma_0 + 2*x*sigma_0 - mu_chi*mu_chi, 3);
  //       out /= (x + sigma_0) ;
  //       out /= (1 + sigma_0*sigma_0 + 2*x*sigma_0 ) ;

  double out = x * sqrt(x*x - 1) ;
        out *= 1 + sigma_0*sigma_0 + 2*x*sigma_0 - mu_chi*mu_chi ;
        out *= (1 + sigma_0*sigma_0 + 2*x*sigma_0) * (1 + x*sigma_0 + 2 * mu_chi) 
                + pow(mu_chi, 2) * (1 + x*sigma_0) ;
        out /= ;

  return out ;
}

//--------------------------------------------------------------
// The full phase-space integral
// The output is in MeV^4
// Here Sigma_0 is the vector self-energy ( Sigma_0 < 0 )
double BNV_B_Chi_Photon::PhaseSpace_Integral(const Baryon& B,
                                              const double& m, 
                                              const double& m_chi, 
                                              const double& Sigma_0, 
                                              const double& x_F) 
{ 
  double q_e = Zaki::Physics::Q_E ;
  double eps = 1e-10 ;
  
  phase_int_mu_chi = m_chi / m ; 
  phase_int_sigma_0 = -Sigma_0 / m ;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
  double err, result;

  Zaki::Math::GSLFuncWrapper<BNV_B_Chi_Photon, double (BNV_B_Chi_Photon::*)(const double&)> 
    Fp(this, &BNV_B_Chi_Photon::PhaseSpace_Integrand);     

  gsl_function F = static_cast<gsl_function> (Fp) ; 

  gsl_integration_qag(&F, 1, x_F, 1e-10, 1e-10, 2000, 1, w, &result, &err) ;
  gsl_integration_workspace_free(w);

  // if( abs(m - m_chi) < 1e-1 )
  // {
  //   std::cout << "\n\t " << "Resonance region! "<< "\n" ;
  // }

  // C^2
  double C_2 = B.g * eps * q_e / (8 * m ) ;
        C_2 *= C_2 ;

  result *= C_2 ;
  result *= pow(m, 5) / (2*M_PI*M_PI*M_PI)  ;

  return result ;
}

//--------------------------------------------------------------
// The integrated rate in the resonance region
// The output is in s^-1
double BNV_B_Chi_Photon::Resonance_B_Chi_Photon_Rate(const double& m_chi, 
                                                     const Baryon& in_B,
                                              const double& in_res_density,
                                              const double& in_res_radius
                                              ) 
{ 
  // Finding the resonance baryon density from the total density
  // First column is the total density and the secon column is the baryon's density
  Zaki::Vector::DataSet ds_EOS_find_res_density ({EOS_ds[2], EOS_ds[2] * EOS_ds[in_B.label]}) ;
  ds_EOS_find_res_density.Interpolate(0, 1) ;
  double res_baryon_density = ds_EOS_find_res_density.Evaluate(1, in_res_density) ;

  double kF_2 = pow(3*M_PI*M_PI*res_baryon_density, 2.0/3.0) ;

  kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

  // double m_B = m_eff_ds[in_B.label] ;
  double m_B = m_chi ;

  // Finding the resonance Sigma_0
  // First column is the total density and the secon column is the baryon's self energy 
  Zaki::Vector::DataSet ds_find_res_Vself ({V_self_E_ds[0], V_self_E_ds[in_B.label]}) ;
  ds_find_res_Vself.Interpolate(0, 1) ;

  double Sigma_0 = - ds_find_res_Vself.Evaluate(1, in_res_density) ;
  // double Sigma_0 = -V_self_E_ds[in_B.label][find_res_idx] ;
  // std::cout << "\n\n\t Sigma_0 = " << Sigma_0 << "\n" << std::flush ;

  double x_F = sqrt(m_B*m_B + kF_2) / m_B ;
  // std::cout << "\n\n\t x_F = " << x_F << "\n" ;

  double q_e = Zaki::Physics::Q_E ;
  double eps = 1e-10 ;
  
  phase_int_mu_chi = m_chi / m_chi ; 
  phase_int_sigma_0 = -Sigma_0 / m_chi ;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
  double err, res_integral;

  Zaki::Math::GSLFuncWrapper<BNV_B_Chi_Photon, double (BNV_B_Chi_Photon::*)(const double&)> 
    Fp(this, &BNV_B_Chi_Photon::PhaseSpace_Integrand);     

  gsl_function F = static_cast<gsl_function> (Fp) ; 

  gsl_integration_qag(&F, 1, x_F, 1e-10, 1e-10, 2000, 1, w, &res_integral, &err) ;
  gsl_integration_workspace_free(w);

  // This factor doesn't include the mixing (sin(2*theta)) anymore
  double C_2 = in_B.g * q_e / ( 16 * m_chi ) ;
        C_2 *= C_2 ;

  res_integral *= C_2 ;
  res_integral *= pow(m_chi, 6) / (2*M_PI*M_PI*M_PI)  ;

  // The radius at which the resonance occurs (in km)
  double R_res = in_res_radius ;

  
  // std::cout << "\n\n\t res_integral-1 = " << res_integral << "\n" ;
  // ..........................................................
  // Reminder: 
  //  r     -> 0
  // M(r)   -> 1
  // n(r)   -> 5
  // nu(r)  -> 6 
  Zaki::Vector::DataSet pulsar_profile = *pulsar.GetProfile() ;
  pulsar_profile.Interpolate(0, {1, 5, 6}) ;

  // ..........................................................
  // Finding he width in the mass of B for the resonace region (e.g., over R_res +- 50 m )
  // ..........................................................

  // The density at R_res + 50 m
  double R_max = R_res + R_res_width ;
  if (R_max > pulsar_profile[0][-1]) 
    R_max = pulsar_profile[0][-1] ;

  double density_at_R_max = pulsar_profile.Evaluate(5, R_max) ; 
  // std::cout << "\n\n\t density_at_R_max = " << density_at_R_max << "\n" ;

  // The density at R_res - 50 m
  double density_at_R_min = pulsar_profile.Evaluate(5, R_res - R_res_width) ; 
  // std::cout << "\n\n\t density_at_R_min = " << density_at_R_min << "\n" ;

  Zaki::Vector::DataSet ds_find_res_meff ({m_eff_ds[0], m_eff_ds[in_B.label]}) ;
  ds_find_res_meff.Interpolate(0, 1) ;

  // [ m_eff - m_chi ] at R_res + 50 m
  double delta_m_R_max = ds_find_res_meff.Evaluate(1, density_at_R_max) - m_chi ;
  std::cout << "\n\n\t delta_m_R_max = " << delta_m_R_max << "\n" ;

  // [ m_chi - m_eff ]  at R_res - 50 m
  double delta_m_R_min = m_chi - ds_find_res_meff.Evaluate(1, density_at_R_min) ;
  std::cout << "\n\n\t delta_m_R_min = " << delta_m_R_min << "\n" ;

  // The derivative of baryon's mass as afunction of radius at the resonance,
  //  i.e.,  m_B^* = m_chi :

  Zaki::Vector::DataSet ds_meff_n ({EOS_ds[2], m_eff_ds[in_B.label]}) ;
  ds_meff_n.Interpolate(0, 1) ;

  Zaki::Vector::DataColumn dc_meff_r = ds_meff_n.Evaluate(1, pulsar_profile[5]) ;

  Zaki::Vector::DataSet ds_meff_r({ pulsar_profile[0], dc_meff_r}) ;
  ds_meff_r.Interpolate(0,1) ;
  
  // unit is MeV per km
  double deriv_m_vs_r = ds_meff_r.Derivative(1, R_res) ;
  // std::cout << "\n\n\t deriv_m_vs_r = " << deriv_m_vs_r << "\n" ;
  // ..........................................................

  // we now multiply the 'result' which is the slowyly vaying part of the rate
  //  but the result of integraion of a small radial range (in the resonance region)
  //
  //  The overal unit [ res_integral ] = [ MeV^4 . km ] 
  res_integral *= 2 * eps *  (  atan(delta_m_R_max / (2 * eps) )  
                              + atan(delta_m_R_min / (2 * eps) )
                              ) / deriv_m_vs_r  ;
  // std::cout << "\n\n\t res_integral-2 = " << res_integral << "\n" ;

  // Converting "MeV^3" into "1/fm^3"
  //  The overal unit [ res_integral ] = [ MeV . fm^-3 . km ] 
  res_integral *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 
  
  // Converting "1/fm^3" into "1/km^3"
  //  The overal unit [ res_integral ] = [ MeV . km^-2 ] 
  res_integral *= 1e+54 ; 

  // Converting "MeV" into "s^-1"
  //  The overal unit [ res_integral ] = [ s^-1 . km^-2 ]
  res_integral *= 1e-3 * Zaki::Physics::GEV_2_S ; 


  // The total mass at R_res (in km)
  double M_at_R_res = pulsar_profile.Evaluate(1, R_res) ;
  // std::cout << "\n\n\t M_at_R_res = " << M_at_R_res << "\n" ;

  // The metric nu function at R_res
  double nu_at_R_res = pulsar_profile.Evaluate(6, R_res) ;
  // std::cout << "\n\n\t nu_at_R_res = " << nu_at_R_res << "\n" ;

  // Multplying the rest of the factors in the volumn integral
  //  The overal unit [ res_integral ] = [ s^-1 ]
  res_integral *= 4*M_PI*R_res*R_res ;  

  res_integral /= sqrt(1. - 2.*M_at_R_res / R_res) ; 
  res_integral *= exp(nu_at_R_res) ; 

  return res_integral ;
}

//--------------------------------------------------------------
Zaki::Vector::DataColumn BNV_B_Chi_Photon::B_Chi_Photon_Rate(
                          const double& m_chi, 
                          const Baryon& in_B)
                          // const Zaki::Vector::DataColumn& Sigma_0, 
                          // const Zaki::Vector::DataColumn& m_B, 
                          // const Zaki::Vector::DataColumn& in_n)
{  
  Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI*n_B_ds[in_B.label]).pow(2.0/3.0) ;

  kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

  Zaki::Vector::DataColumn m_B = m_eff_ds[in_B.label] ;
  Zaki::Vector::DataColumn Sigma_0 = -V_self_E_ds[in_B.label] ;
  Zaki::Vector::DataColumn x_F = (m_B*m_B + kF_2).sqrt() / m_B ;

  Zaki::Vector::DataColumn tmp ;
  for (size_t i = 0; i < m_B.Size(); i++)
  {
    tmp.vals.emplace_back( PhaseSpace_Integral(in_B, m_B[i], m_chi, Sigma_0[i], x_F[i]) ) ;
  }
  
  if(false)
  {
    // Plotting the decay rate per volume as a function of density
    Zaki::Vector::DataSet tmp_plot({EOS_ds[2], tmp}) ;
    tmp_plot[0].label = "$n (fm^{-3})$";

    // Converting "MeV^3" into "1/fm^3"
    tmp_plot[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ;

    // dc_B_chi_photon_rate_r *= 1e+54 ; 

    // Converting "MeV" into "s^-1"
    tmp_plot[1] *= 1e-3 * Zaki::Physics::GEV_2_S ; 

    // // Scaling from eps = 1e-10 to eps = ?:
    // tmp_plot[1] *= 1e-16 ;

    tmp_plot[1].label = "$\\Gamma ({\\rm fm}^{-3} \\cdot s^{-1})$";

    char tmp_char[50] ;
    snprintf(tmp_char, sizeof(tmp_char), "%.0f", m_chi) ;

    tmp_plot.SetWrkDir(wrk_dir + "/BNV_2022/results/B_Chi_Photon_Decay/"+model) ;
    tmp_plot.Plot(0, 1, "Decay_per_V_vs_n/" + in_B.label +"_Decay_per_V_vs_n/" + 
                        in_B.label+"_decay_per_V_vs_n_"+std::string(tmp_char)+".pdf",
                        "$\\varepsilon = 10^{-10}\\, {\\rm MeV}$") ;
  }

  return tmp ;
}

//--------------------------------------------------------------
// Plots the baryon's rest-energy as a function of radius
void BNV_B_Chi_Photon::Plot_RestEnergy_Radius()
{
  // Building the rest energy dataset as a function of density
  Zaki::Vector::DataSet ds_E_n ({EOS_ds[2], 
            m_eff_ds[neutron.label] + V_self_E_ds[neutron.label], 
            m_eff_ds[lambda.label] + V_self_E_ds[lambda.label],
            m_eff_ds[proton.label] + V_self_E_ds[proton.label]}) ;

  ds_E_n.Interpolate(0, {1,2,3}) ;

  // ......................................
  // Building the rest energy data-columns as a function of radius
  Zaki::Vector::DataColumn dc_neutron_E_r = 
    ds_E_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;
  
  Zaki::Vector::DataColumn dc_lambda_E_r = 
    ds_E_n.Evaluate(2, pulsar.GetProfile()->operator[](5)) ;

  Zaki::Vector::DataColumn dc_proton_E_r = 
    ds_E_n.Evaluate(3, pulsar.GetProfile()->operator[](5)) ;
  // ......................................

  // ......................................
  //  Combining the data-columns into a dataset
  Zaki::Vector::DataSet ds_E_r({ pulsar.GetProfile()->operator[](0), 
                      dc_neutron_E_r, dc_lambda_E_r, dc_proton_E_r}) ;
  // ......................................
  // Fixing the labels
  ds_E_r[1].label = "$E_n^0$" ;
  ds_E_r[2].label = "$E_{\\Lambda}^0$" ;
  ds_E_r[3].label = "$E_{p}^0$" ;
  // ......................................

  ds_E_r.SetWrkDir(wrk_dir + "/BNV_2022/results/B_Chi_Photon_Decay/"+model) ;
  
  Zaki::Vector::DataSet::PlotParam plt_par ;

  // Setting the boundaries
  plt_par.SetXAxis({ds_E_r[0].Min(), ds_E_r[0].Max()}) ;

  // Finding the lowest and highest points on the curves
  double Y_min = std::min({ds_E_r[1].Min(), ds_E_r[2].Min(), ds_E_r[3].Min()}) ;
  double Y_max = std::max({ds_E_r[1].Max(), ds_E_r[2].Max(), ds_E_r[3].Max()}) ;

  plt_par.SetYAxis({Y_min*0.98, Y_max*1.02}) ;
  
  // This is for setting Y-axis ticks
  int min_tick = Y_min*0.98 - fmod(Y_min*0.98, 50) ;
  int max_tick = Y_max*1.02 - fmod(Y_max*1.02, 50) ;
  std::vector<double> y_tick_set ;
  for (int i = min_tick; i <= max_tick; i+= 50)
  {
    y_tick_set.emplace_back(i) ;
  }
  
  // plt_par.SetYTicks({{900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350}}) ;
  plt_par.SetYTicks({y_tick_set}) ;

  plt_par.SetXAxisLabel("$R\\, [\\, {\\rm km}\\,]$") ;
  plt_par.SetYAxisLabel("$E_B\\, [\\, {\\rm MeV}\\,]$") ;
  plt_par.SetLegend({"upper right", 1.0, 1.0}) ;
  plt_par.SetGrid() ;
  ds_E_r.SetPlotPars(plt_par) ;

  ds_E_r.Plot(0, {1,2,3}, "E_vs_radius.pdf", pulsar.GetName() + "\n" + model) ;

}
//--------------------------------------------------------------
// Plots the effective mass as a function of radius
void BNV_B_Chi_Photon::Plot_Meff_Radius()
{
  Zaki::Vector::DataSet ds_meff_n ({EOS_ds[2], m_eff_ds[neutron.label], 
                                      m_eff_ds[lambda.label], 
                                      m_eff_ds[proton.label]}) ;
  ds_meff_n.Interpolate(0, {1,2,3}) ;

  Zaki::Vector::DataColumn dc_neutron_meff_r = 
    ds_meff_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;
  
  Zaki::Vector::DataColumn dc_lambda_meff_r = 
    ds_meff_n.Evaluate(2, pulsar.GetProfile()->operator[](5)) ;

  Zaki::Vector::DataColumn dc_proton_meff_r = 
    ds_meff_n.Evaluate(3, pulsar.GetProfile()->operator[](5)) ;

  Zaki::Vector::DataSet ds_meff_r({ pulsar.GetProfile()->operator[](0), 
                  dc_neutron_meff_r, dc_lambda_meff_r, dc_proton_meff_r}) ;
  // ds_meff_r[1].label = "$m^*_n$" ;
  // ds_meff_r[2].label = "$m^*_{\\Lambda}$" ;

  ds_meff_r.SetWrkDir(wrk_dir + "/BNV_2022/results/B_Chi_Photon_Decay/"+model) ;
  
  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetXAxis({ds_meff_r[0].Min()*0.99, ds_meff_r[0].Max()*1.01}) ;
  plt_par.SetYAxis({ds_meff_r[1].Min()*0.98, ds_meff_r[2].Max()*1.02}) ;
  plt_par.SetXAxisLabel("$R\\, [\\, {\\rm km}\\,]$") ;
  plt_par.SetYAxisLabel("$m_B^*\\, [\\, {\\rm MeV}\\,]$") ;
  plt_par.SetGrid() ;
  plt_par.SetLegend({"upper left", 0.0, 1.0}) ;
  ds_meff_r.SetPlotPars(plt_par) ;

  ds_meff_r.Plot(0, {{1, "$m^*_n$"} , {2, "$m^*_{\\Lambda}$"}, {3, "$m^*_{p}$"}}, 
                "m_eff_vs_radius.pdf", pulsar.GetName() + "\n" + model) ;


  // ------------------------------------------------------------
  // Plotting the resonance region thickness as a function of radius
  double delta_m_res = 1 ; // MeV

  // ------------------------------------------------------------
  //                    Neutron
  // ------------------------------------------------------------
  // std::cout << "\n \t OK-1!\n" ;
  Zaki::Vector::DataSet ds_Neutron_R_vs_M( { 
                                  ds_meff_r[0].GetSubSet(39, ds_meff_r[0].Size()-1), 
                                  ds_meff_r[1].GetSubSet(39, ds_meff_r[1].Size()-1) } ) ;
  // std::cout << "\n \t OK-2!\n" ;
  ds_Neutron_R_vs_M.Interpolate(1, 0) ;
  // std::cout << "\n \t OK-3!\n" ;
  // return ;

  Zaki::Vector::DataColumn dc_delta_r_neutron ;
  dc_delta_r_neutron.Reserve(ds_Neutron_R_vs_M[0].Size()) ;

  std::cout.precision(10) ;
  for (size_t i = 0; i < ds_Neutron_R_vs_M[0].Size()-1; i++)
  {
    double next_r = ds_Neutron_R_vs_M.Evaluate(0,  ds_Neutron_R_vs_M[1][i] + delta_m_res ) ;

    // std::cout << "\t ds_Neutron_R_vs_M[1]["<<i<<"] = " << ds_Neutron_R_vs_M[1][i] 
    //           << ", R_i = "<<ds_Neutron_R_vs_M[0][i]<<", R_next = " << next_r <<  "\n" ;
    dc_delta_r_neutron.vals.emplace_back(next_r - ds_Neutron_R_vs_M[0][i] ) ;          
  }

  // return ;

  // Zaki::Vector::DataColumn dc_delta_r_neutron ;
  // dc_delta_r_neutron.Reserve(dc_neutron_meff_r.Size()) ;

  // for (size_t i = 0; i < dc_neutron_meff_r.Size(); i++)
  // {
  //   int next_idx = dc_neutron_meff_r.GetClosestIdx ( dc_neutron_meff_r[i] + delta_m_res ) ;

  //   std::cout << "\t dc_neutron_meff_r["<<i<<"] = " << dc_neutron_meff_r[i] << ", next_idx = " << next_idx <<  "\n" ;

  //   double delta_r = ds_meff_r[0][next_idx] - ds_meff_r[0][i] ;
  //   double mid_point = (ds_meff_r[0][next_idx] + ds_meff_r[0][i] ) / 2.0 ;
  //   // dc_delta_r_neutron.vals.emplace_back( delta_r / mid_point) ;

  //   if ( i+1 != dc_neutron_meff_r.Size() && dc_neutron_meff_r[i+1] - dc_neutron_meff_r[i] >  delta_m_res )
  //   {
  //     delta_r = (delta_m_res / ( dc_neutron_meff_r[i+1] -  dc_neutron_meff_r[i]) )
  //                * (ds_meff_r[0][i+1] - ds_meff_r[0][i]) ;
  //   }

  //   dc_delta_r_neutron.vals.emplace_back( delta_r) ;
  // }

  // ------------------------------------------------------------
  //                    Lambda
  // ------------------------------------------------------------
  // Zaki::Vector::DataColumn dc_delta_r_lambda ;
  // dc_delta_r_lambda.Reserve(dc_lambda_meff_r.Size()) ;

  // for (size_t i = 0; i < dc_lambda_meff_r.Size(); i++)
  // {
  //   int next_idx = dc_lambda_meff_r.GetClosestIdx ( dc_lambda_meff_r[i] + delta_m_res ) ;
  //   double delta_r = ds_meff_r[0][next_idx] - ds_meff_r[0][i] ;
  //   double mid_point = (ds_meff_r[0][next_idx] + ds_meff_r[0][i] ) / 2.0 ;
  //   // dc_delta_r_neutron.vals.emplace_back( delta_r / mid_point) ;
  //   dc_delta_r_lambda.vals.emplace_back( delta_r) ;
  // }
  
  // ------------------------------------------------------------

  // ds_meff_r.Interpolate(0, {1,2}) ;
  Zaki::Vector::DataSet ds_dr_res({ ds_Neutron_R_vs_M[1].GetSubSet(0, ds_Neutron_R_vs_M[1].Size()-2), 
                                    dc_delta_r_neutron
                                    //, dc_lambda_meff_r, dc_delta_r_lambda
                                    // delta_m_res / ds_meff_r.Derivative(1)[1], 
                                    // delta_m_res / ds_meff_r.Derivative(2)[1]
                                    }) ;

  ds_dr_res[0].label = "$m_{\\chi} [ MeV ]$" ;
  ds_dr_res[1].label = "$\\Delta r_n [ km ]$" ;
  // ds_dr_res[2].label = "$m_{\\chi} [ MeV ]$" ;
  // ds_dr_res[3].label = "$\\Delta r_{\\Lambda} [ km ]$" ;

  // Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.Reset() ;
  plt_par.SetXAxis({300, 800}) ;
  plt_par.SetYAxis({1e-3, 1e0}) ;

  ds_dr_res.SetPlotPars(plt_par) ;

  ds_dr_res.SetWrkDir(wrk_dir + "/BNV_2022/results/B_Chi_Photon_Decay/"+model) ;
  ds_dr_res.SemiLogYPlot(0, 1, "Res_thick_vs_radius_neutron.pdf", 
                          "Thickness of regions with $| m_n - m_{\\chi} | = 1$ MeV.") ;

  // plt_par.SetXAxis({700, 1100}) ;
  // plt_par.SetYAxis({1e-3, 1.3e0}) ;
  
  // ds_dr_res.SetPlotPars(plt_par) ;
  // ds_dr_res.SemiLogYPlot(2, 3, "Res_thick_vs_radius_lambda.pdf", "Thickness of regions with $| m_{\\Lambda} - m_{\\chi} |  = 1$ MeV.") ;

}

//--------------------------------------------------------------
/// Evaluates the total decay rate for the pulsar
BNV_B_Chi_Photon::Decay_Rate_Type 
BNV_B_Chi_Photon::EvalDecayRate(const double& m_chi, const Baryon& in_B)
{
  // double m_chi = 950 ;
  std::cout << "\n -------------------- m_chi = " << m_chi 
            <<  " MeV ---- B = " << in_B.label <<"-----------------\n" ;

  Zaki::Vector::DataColumn B_chi_photon_rate =  
        B_Chi_Photon_Rate(m_chi, in_B) ;

  Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;

  Zaki::Vector::DataSet ds_B_chi_photon_rate_n ({EOS_ds[2], B_chi_photon_rate}) ;
  ds_B_chi_photon_rate_n.Interpolate(0, 1) ;
  
  Zaki::Vector::DataColumn dc_B_chi_photon_rate_r = 
    ds_B_chi_photon_rate_n.Evaluate(1, n_r) ;


  // Converting "MeV^3" into "1/fm^3"
  dc_B_chi_photon_rate_r *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "1/fm^3" into "1/km^3"
  dc_B_chi_photon_rate_r *= 1e+54 ; 

  // Converting "MeV" into "s^-1"
  dc_B_chi_photon_rate_r *= 1e-3 * Zaki::Physics::GEV_2_S ; 

  Zaki::Vector::DataColumn r = pulsar.GetProfile()->operator[](0) ;
  
  //----------------------------------------------------------------------
  //          Plotting decay per V [km^-3.s^-1] vs R [km]
  //----------------------------------------------------------------------
  if (false)
  {
    Zaki::Vector::DataSet ds_B_chi_photon_rate_r ({r, dc_B_chi_photon_rate_r}) ;
    ds_B_chi_photon_rate_r[1].label = "$\\Gamma ({\\rm km}^{-3} \\cdot s^{-1})$";
    char tmp_char[50] ;
    snprintf(tmp_char, sizeof(tmp_char), "%.0f", m_chi) ;
    ds_B_chi_photon_rate_r.SetWrkDir(wrk_dir + "/BNV_2022/results/B_Chi_Photon_Decay/"+model) ;
    ds_B_chi_photon_rate_r.Plot(0, 1, "Decay_per_V_vs_R/" + in_B.label +"_Decay_per_V_vs_R/" + 
                        in_B.label+"_decay_per_V_vs_R_"+std::string(tmp_char)+".pdf",
                        "$\\varepsilon = 10^{-10}\\, {\\rm MeV}$") ;
  }
  //----------------------------------------------------------------------
  
  Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
  Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

  Zaki::Vector::DataColumn decay_integrand_dc = 4*M_PI*r*r*dc_B_chi_photon_rate_r ;  
                          decay_integrand_dc /= (1 - 2*M_r / r).sqrt() ; 
                          decay_integrand_dc *= exp(nu_r) ; 

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // ---------------------------Testing Begins----------------------------
  bool resonance_region_exists = false ;
  double R_max = 5 ;
  double R_min = R_max ;

  double interp_n_resonance = 0 ; 
  double interp_r_resonance = 0 ;

  Zaki::Vector::DataSet ds_meff_n ({EOS_ds[2], m_eff_ds[in_B.label]}) ;
  ds_meff_n.Interpolate(0, 1) ;
  double m_B_min = ds_meff_n.Evaluate(1, n_r[0]) ;
  double m_B_max = ds_meff_n.Evaluate(1, n_r[-1]) ;

  std::cout << " m_min = " << m_B_min << " MeV,\t\t" ;
  std::cout << " m_max = " << m_B_max << " MeV.\n" ;

  if ( m_chi <= m_B_min || m_B_max <= m_chi ) 
  {
    std::cout << "\n\t No resonance region.\n" ;
  }
  else
  {
    resonance_region_exists = true ;

    double n_resonance = m_eff_ds[0] [ m_eff_ds[in_B.label].GetClosestIdx(m_chi) ] ;
    // std::cout << "\n\n\t n_resonance ~ " << n_resonance << "\n" ;

    double r_Resonance = r [ n_r.GetClosestIdx(n_resonance) ]  ;
    // std::cout << "\n\n\t r_Resonance ~ " <<  r_Resonance << "\n" ;

    Zaki::Vector::DataSet ds_n_res_finder ({1.0/m_eff_ds[in_B.label], m_eff_ds[0]}) ;
    ds_n_res_finder.Interpolate(0, 1) ;

    interp_n_resonance = ds_n_res_finder.Evaluate(1,  1.0/ m_chi ) ;
    std::cout << "\t interp_n_resonance = " << interp_n_resonance << "\n" ;

    Zaki::Vector::DataSet ds_r_res_finder ({1.0 / n_r.GetSubSet(39, r.Size() - 1), r.GetSubSet(39, r.Size() - 1)}) ;
    ds_r_res_finder.Interpolate(0, 1) ;

    if (interp_n_resonance >= n_r.Min() )
    {
      interp_r_resonance = ds_r_res_finder.Evaluate(1,  1.0/ interp_n_resonance ) ;
      R_min = interp_r_resonance - R_res_width ;
    }
    else
    {
      interp_r_resonance = r[-1] ;
      R_min = r[-1] ;
    }

    std::cout << "\t interp_r_resonance = " << interp_r_resonance << "\n" ;

    R_max = interp_r_resonance + R_res_width ;
    if ( R_max > r[-1] )
    { 
      R_max= r[-1] ;
    }
  }
  // return 0 ;
  // Zaki::Vector::DataSet ds_meff_n ({EOS_ds[2], m_eff_ds[neutron.label], 
  //                                     m_eff_ds[lambda.label]}) ;
  // ds_meff_n.Interpolate(0, {1,2}) ;

  // Zaki::Vector::DataColumn dc_neutron_meff_r = 
  //   ds_meff_n.Evaluate(1, n_r) ;
  
  // Zaki::Vector::DataColumn dc_lambda_meff_r = 
  //   ds_meff_n.Evaluate(2, n_r) ;

  // Zaki::Vector::DataSet ds_meff_r({ pulsar.GetProfile()->operator[](0), dc_neutron_meff_r, dc_lambda_meff_r}) ;
  // ds_meff_r[1].label = "$m^*_n$" ;
  // ds_meff_r[2].label = "$m^*_{\\Lambda}$" ;

  // Zaki::Vector::DataSet m_B_vs_R ( { r.GetSubSet(39, r.Size()-1), 
  //                   ds_meff_r[1].GetSubSet(39, ds_meff_r[1].Size()-1)} ) ; // 'ds_meff_r[1]' has to change

  // m_B_vs_R.Interpolate(1,0) ;
  // double r_resonance = m_B_vs_R.Evaluate(1, m_chi) ;
  // ---------------------------Testing Ends----------------------------
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  std::cout << " R_min = " << R_min << " km,\t\t" ;
  std::cout << " R_max = " << R_max << " km.\n" ;

  Zaki::Vector::DataSet decay_integral_ds({r, decay_integrand_dc}) ;

  decay_integral_ds.Interpolate(0, 1) ;    

  double result = decay_integral_ds.Integrate(1, {r[0], R_min}) ;
        result += decay_integral_ds.Integrate(1, {R_max, r[-1]}) ;

  double result_full = decay_integral_ds.Integrate(1, {r[0], r[-1]}) ;

  double result_resonance = 0 ;
  if (resonance_region_exists)
  {
    result_resonance = Resonance_B_Chi_Photon_Rate(m_chi, in_B, interp_n_resonance,
                                                                interp_r_resonance) ;
    // result_plus_resonance += tmp_resonance_B_Chi_Photon_Rate ;

    std::cout << "\n\t Resoannce decay rate (B_dot/B) (only resonance region): " 
              << result_resonance /pulsar.GetSeqPoint().b << " per year.\n" ;
  }

  // Convert s^-1 into yr^-1
  result                *= 3600*24*365 ;
  result_full           *= 3600*24*365 ;
  result_resonance      *= 3600*24*365 ;

  // std::cout << "\n\t Full decay rate (B_dot): " << result << " per year.\n" ;
  std::cout << "\n\t Full decay rate (B_dot/B) (resonance not removed): " << result_full/pulsar.GetSeqPoint().b << " per year.\n" ;

  std::cout << "\n -------------------- \n" ;

  Decay_Rate_Type decay_result( result/pulsar.GetSeqPoint().b, 
                                result_full/pulsar.GetSeqPoint().b, 
                                result_resonance/pulsar.GetSeqPoint().b) ;

  return  decay_result ;
} 

//--------------------------------------------------------------
// Finds the eps scaling factor (when resonance is included)
double BNV_B_Chi_Photon::FindEpsScaling(const Decay_Rate_Type& in_dec_rate, const double& in_gam_lim)
{
  double a = in_dec_rate.no_res ;
  double b = in_dec_rate.res_only ;

  if ( b == 0 )
  {
    return sqrt(in_gam_lim / a) ;
  }

  // return sqrt(in_gam_lim / a) ;

  if ( a/b < 1e-6 )
  {
    return in_gam_lim / b ;
  }

  double out = sqrt(b*b + 4.*a*in_gam_lim) - b ;
        out /= 2.*a ;

  // std::cout << "\n\t a = " << a << ", b = " << b << "\n";

  return out ;
}

//--------------------------------------------------------------
/// Finds and plots the decay limits
void BNV_B_Chi_Photon::Plot_DecayLimits() 
{
  // This dataset contains the upper limit on eps in MeV
  Zaki::Vector::DataSet dec_lim_ds ;
  dec_lim_ds.Reserve(9, m_chi_res) ;

  for (size_t i = 0; i < m_chi_res; i++)
  {
    double m_chi = 200 + i*2 ;
    dec_lim_ds[0].vals.emplace_back(m_chi) ;

    // ------------------------------------------
    //                  Neutron 
    // ------------------------------------------
    Decay_Rate_Type dec_result = EvalDecayRate(m_chi, neutron) ;
    // if (dec_result.max == 0)
    //   dec_result.max = 1e-60 ;

    double eps_lim_min = 1e-10*sqrt(7.83129*1e-10 / dec_result.no_res ) ;
    double eps_lim_normal = 1e-10*sqrt(7.83129*1e-10 / dec_result.normal ) ;

    // More complicated scaling due to the resonance region:
    double eps_lim_max = 1e-10*FindEpsScaling(dec_result, 7.83129*1e-10) ;
    double tmp_scaling =  FindEpsScaling(dec_result, 7.83129*1e-10) ;
    std::cout << "\n\t Scaling = " << tmp_scaling << "\n" ;
    std::cout << "\t No-Res Scaling = " <<  sqrt(7.83129*1e-10 / dec_result.no_res ) << "\n" ;

    dec_lim_ds[1].vals.emplace_back( eps_lim_min );
    dec_lim_ds[2].vals.emplace_back( eps_lim_normal );
    dec_lim_ds[3].vals.emplace_back( eps_lim_max );

    dec_lim_ds[7].vals.emplace_back( Vacuum_Decay_Br(neutron, m_chi, eps_lim_max)
                                    );

    // ------------------------------------------
    //                  Lambda 
    // ------------------------------------------
    dec_result = EvalDecayRate(m_chi, lambda) ;
    // if (dec_result.max == 0)
    //   dec_result.max = 1e-60 ;

    eps_lim_min = 1e-10*sqrt(7.83129*1e-10 / dec_result.no_res ) ;
    eps_lim_normal = 1e-10*sqrt(7.83129*1e-10 / dec_result.normal ) ;

    // More complicated scaling due to the resonance region:
    eps_lim_max = 1e-10*FindEpsScaling(dec_result, 7.83129*1e-10) ;
    // eps_lim_max = 1e-10*(7.83129*1e-10 / dec_result.max ) ;

    dec_lim_ds[4].vals.emplace_back( eps_lim_min );
    dec_lim_ds[5].vals.emplace_back( eps_lim_normal );
    dec_lim_ds[6].vals.emplace_back( eps_lim_max );

    dec_lim_ds[8].vals.emplace_back( Vacuum_Decay_Br(lambda, m_chi, eps_lim_max)
                                    );
  }
  
  dec_lim_ds[0].label = "$m_{\\chi}\\, [ MeV ]$" ;
  dec_lim_ds[1].label = "$\\varepsilon_{n}^{min}\\, [ MeV ]$" ;
  dec_lim_ds[2].label = "$\\varepsilon_{n}^{nor}\\, [ MeV ]$" ;
  dec_lim_ds[3].label = "$\\varepsilon_{n}^{max}\\, [ MeV ]$" ;
  dec_lim_ds[4].label = "$\\varepsilon_{\\Lambda}^{min}\\, [ MeV ]$" ;
  dec_lim_ds[5].label = "$\\varepsilon_{\\Lambda}^{nor}\\, [ MeV ]$" ;
  dec_lim_ds[6].label = "$\\varepsilon_{\\Lambda}^{max}\\, [ MeV ]$" ;
  dec_lim_ds[7].label = "Br$(n \\to \\chi \\gamma)$" ;
  dec_lim_ds[8].label = "Br$(\\Lambda \\to \\chi \\gamma)$" ;

  dec_lim_ds.SetWrkDir(wrk_dir + "/BNV_2022/results/B_Chi_Photon_Decay/"+model) ;
  
  // Changing the plot parameters
  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetXAxis({200, 1400}) ;
  plt_par.SetYAxis({1e-40, 1e-12}) ;
  plt_par.SetYAxisLabel("$\\varepsilon_B\\, [\\, {\\rm MeV} \\, ]$") ;
  plt_par.SetGrid() ;
  plt_par.SetLegend({"lower right", 1.0, 0.0}) ;
  dec_lim_ds.SetPlotPars(plt_par) ;

  dec_lim_ds.SemiLogYPlot(0, {{3, "$\\varepsilon_{n}^{\\rm Best}$"}, 
          {6, "$\\varepsilon_{\\Lambda}^{\\rm Best}$"}}, "Eps_vs_m_chi.pdf",
           pulsar.GetName() + "\n" + model) ;

  plt_par.SetXAxis({200, 1120}) ;
  plt_par.SetYAxis({1e-64, 1e-13}) ;
  plt_par.SetXTicks({{200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100}}) ;
  plt_par.SetYAxisLabel("${\\rm Br} \\, \\left( B \\to \\chi \\gamma \\, \\right)$") ;
  plt_par.SetGrid() ;
  plt_par.SetLegend({"center", 0.5, 0.92}) ;

  dec_lim_ds.SetPlotPars(plt_par) ;

  dec_lim_ds.SemiLogYPlot(0, {7,8}, "Br_vs_m_chi.pdf",
                          pulsar.GetName() + "\n" + model) ;
}

//--------------------------------------------------------------
// Outputs the vacuum branching ratio of 
//      B -> chi + photon
// 'in_m_chi' & 'in_eps_lam' must be in MeV
// ____________________________________
double BNV_B_Chi_Photon::Vacuum_Decay_Br(const Baryon& B, 
                const double& in_m_chi, const double& in_eps_lam) 
  const
{
  // ...................................
  //        CHANGE OF UNITS
  // ...................................
  // Chi's mass in units of fm^-1
  double m_chi_fm = in_m_chi * Zaki::Physics::MEV_2_INV_FM ;

  // Mixing in units of fm^-1
  double eps_lam_fm = in_eps_lam * Zaki::Physics::MEV_2_INV_FM ;

  // Baryon's mass in units of fm^-1
  double m_B_0 = B.mass ;
  // ...................................

  if(m_B_0 < m_chi_fm)
    return 0 ;

  double out =  B.g * B.g ;
        out *= Zaki::Physics::Q_E*Zaki::Physics::Q_E ;
        out *= eps_lam_fm * eps_lam_fm ;
        out /= 128*M_PI ;

        out *= (m_B_0 - m_chi_fm) ;
        out *= pow(m_B_0 + m_chi_fm, 3) ;
        out /= pow(m_B_0, 5) ;

  // out *= m_B_0 ;
  // out /= pow(m_B_0 - m_chi_fm, 2) ;
  // out *= pow(1 - pow(m_chi_fm / m_B_0, 2), 3) ;

  // Converting the rate from [fm^-1] to [s^-1]
  out *= Zaki::Physics::LIGHT_C_M_S * 1e15 ;

  // Converting rate into branching ratio
  out *= B.tau ;

  // Returning the vacuum branching ratio
  return out ;
}

//--------------------------------------------------------------
/// Imports the EoS
void BNV_B_Chi_Photon::ImportEOS()
{
  // Zaki::Vector::DataSet ds_EOS(dir.ParentDir().ParentDir() +"/EOS/CompOSE/DS(CMF)-"+ model,
  //                             "DS(CMF)-" + model + ".eos") ;
  EOS_ds.SetWrkDir(wrk_dir) ;
  EOS_ds.Import("/EOS/CompOSE/"+ model + "/" + model + ".eos") ;
  // EOS_ds.SetWrkDir(wrk_dir +"/results") ;
  // ds_EOS[2].label = "$n (fm^{-3})$" ;

  // Finding the baryon's individual density
  n_B_ds = { EOS_ds[2] * EOS_ds[neutron.label], EOS_ds[2] * EOS_ds[lambda.label] } ;
  n_B_ds[0].label = "10" ;
  n_B_ds[1].label = "100" ;
}

//--------------------------------------------------------------
// Imports the effective mass set
void BNV_B_Chi_Photon::ImportEffMass() 
{
  // Zaki::Vector::DataSet ds_meff(dir.ParentDir().ParentDir() +"/EOS/CompOSE/DS(CMF)-"+ model,
  //                               "DS(CMF)-" + model + "_m_eff.micro") ;
  // ds_meff.SetWrkDir(dir.ParentDir() +"/results") ;
  m_eff_ds.SetWrkDir(wrk_dir) ;
  m_eff_ds.Import("/EOS/CompOSE/"+ model + 
                    "/" + model + "_m_eff.micro") ;

  // m_eff_ds[0].label = "$n (fm^{-3})$" ;
  // m_eff_ds[1].label = "$m_n$" ;
  // m_eff_ds[2].label = "$m_{\\Lambda}$" ;
  // m_eff_ds[3].label = "$m_{\\Sigma^{-}}$" ;

  // m_eff_ds.Plot(0, {1,2,3}, "M_eff.pdf") ;

  // Converting MeV to fm^-1
  // for (size_t i = 1; i < m_eff_ds.Dim().size(); i++)
  // {
  //   m_eff_ds[i] *= Zaki::Physics::MEV_2_INV_FM ;
  // }

  // m_eff_ds.Interpolate(0, {1, 2}) ;
}

//--------------------------------------------------------------
// Imports the vector self-energy set
void BNV_B_Chi_Photon::ImportVSelfEnergy() 
{
  // Zaki::Vector::DataSet m_eff_ds(wrk_dir + in_dir, f_name) ;

  // Importing the vector self-energy
  // Zaki::Vector::DataSet ds_V(dir.ParentDir().ParentDir() +"/EOS/CompOSE/DS(CMF)-"+ model,
  //                             "DS(CMF)-" + model + "_V.micro") ;
  // ds_V.SetWrkDir(dir.ParentDir() +"/results") ;
  // ds_V[0].label = "$n (fm^{-3})$" ;
  // ds_V[baryon.V_idx].label = "$V_" + baryon.short_name + " (MeV)$" ;

  V_self_E_ds.SetWrkDir(wrk_dir) ;
  V_self_E_ds.Import("/EOS/CompOSE/"+ model + 
                      "/" + model + "_V.micro") ;

  // V_self_E_ds[0].label = "$n (fm^{-3})$" ;
  // V_self_E_ds[1].label = "$V_n$" ;
  // V_self_E_ds[2].label = "$V_{\\Lambda}$" ;
  // V_self_E_ds[3].label = "$V_{\\Sigma^{-}}$" ;

  // V_self_E_ds.Plot(0, {1,2,3}, "V_Self_E.pdf") ;

  // Converting MeV to fm^-1
  // for (size_t i = 1; i < V_self_E_ds.Dim().size(); i++)
  // {
  //   V_self_E_ds[i] *= Zaki::Physics::MEV_2_INV_FM ;
  // }
  

  // V_self_E_ds.Interpolate(0, {1, 2}) ;
}

//--------------------------------------------------------------
// Plots [ V_self_E_ds + m_eff_ds ]
void BNV_B_Chi_Photon::Plot_RestEnergy_Density()
{
  // Zaki::Vector::DataSet plt_ds({m_eff_ds[0], 
  //     (m_eff_ds[1] + V_self_E_ds[1])/Zaki::Physics::MEV_2_INV_FM, 
  //     (m_eff_ds[2] + V_self_E_ds[2])/Zaki::Physics::MEV_2_INV_FM, 
  //     (m_eff_ds[3] + V_self_E_ds[3])/Zaki::Physics::MEV_2_INV_FM}) ;

  // plt_ds[1].label = "$E_n$" ;
  // plt_ds[2].label = "$E_{\\Lambda}$" ;
  // plt_ds[3].label = "$E_{\\Sigma^{-}}$" ;

  // // std::cout << " * --------------------------------------------------- * \n" ;
  // // std::cout << "\t Rest Energy at n_min = " << m_eff_ds[0][0] << " : \n" ;
  // // std::cout << "\t E(Neutron) :  " << plt_ds[1][0] << " MeV\n" ;
  // // std::cout << "\t E(Lambda)  :  " << plt_ds[2][0] << " MeV\n" ;
  // // std::cout << "\t E(Sigma-)  :  " << plt_ds[3][0] << " MeV\n\n" ;

  // // std::cout << " * --------------------------------------------------- * \n" ;
  // // std::cout << "\t Effective Masses at n_min = " << m_eff_ds[0][0] << " : \n" ;
  // // std::cout << "\t m*(Neutron) :  " << m_eff_ds[1][0]/Zaki::Physics::MEV_2_INV_FM << " MeV\n" ;
  // // std::cout << "\t m*(Lambda)  :  " << m_eff_ds[2][0]/Zaki::Physics::MEV_2_INV_FM << " MeV\n" ;
  // // std::cout << "\t m*(Sigma-)  :  " << m_eff_ds[3][0]/Zaki::Physics::MEV_2_INV_FM << " MeV\n" ;
  // // std::cout << " * --------------------------------------------------- * \n" ;

  // plt_ds.SetWrkDir(wrk_dir) ;
  // plt_ds.Plot(0, {1,2,3}, "Rest_E.pdf") ;


  Zaki::Vector::DataSet plt_ds({m_eff_ds[0], 
      m_eff_ds[1], V_self_E_ds[1], m_eff_ds[2], V_self_E_ds[2], m_eff_ds[3], V_self_E_ds[3]}) ;
  plt_ds.SetWrkDir(wrk_dir) ;
  plt_ds.Plot(0, {1,2,3,4,5,6}, "m_V.pdf", pulsar.GetName() + "\n" + model) ;
}

//--------------------------------------------------------------
// Attaches the pulsar
void BNV_B_Chi_Photon::FindPulsar() 
{ 
  pulsar.SetName("J0348+0432") ;
  pulsar.SetMass({2.01, 0.04}) ;
  pulsar.SetSpinP({ 39.1226569017806, 5e-13 }) ;
  pulsar.SetSpinPDot({ 0.24073e-18, 4e-5 }) ;
  pulsar.SetWrkDir( wrk_dir + "/BNV_2022/results/B_Chi_Photon_Decay") ;
  pulsar.FindProfile("B_Chi_Photon_Decay", model +"/") ;   
  pulsar.PlotRelativeComposition(model + "/" + pulsar.GetName() 
                                  + "_RelComp_vs_R.pdf") ;
  pulsar.PlotAbsoluteComposition(model + "/" + pulsar.GetName() 
                                  + "_AbsComp_vs_R.pdf") ;
  pulsar.PlotChemPotential(model + "/" + pulsar.GetName() 
                                  + "_ChemPot_vs_R.pdf") ;
}

//--------------------------------------------------------------
// Saves the results
void BNV_B_Chi_Photon::Export(const Zaki::String::Directory& in_dir)
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
