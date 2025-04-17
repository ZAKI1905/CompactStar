/*
  BNV_Chi class
*/

// // Creating directories
// #include <sys/stat.h>
// #include <filesystem>

#include <gsl/gsl_integration.h>

#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Vector/DataSet.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/BNV/BNV_Chi.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
// #include "CompactStar/Core/Pulsar.hpp"

using namespace CompactStar ;

double BNV_Chi_pulsar_mass = 0 ;
//==============================================================
bool MassCondition(const CompactStar::NStar& in_star)
{
  return ( BNV_Chi_pulsar_mass - 1e-4 <= in_star.GetSequence().m 
            && 
          in_star.GetSequence().m <= BNV_Chi_pulsar_mass + 1e-4  )  ;
}
//==============================================================

//==============================================================
//                        BNV_Chi class
//==============================================================
// Constructor
BNV_Chi::BNV_Chi(const Process& in_process) 
  : 
  neutron("neutron", "n", "n", "10", Zaki::Physics::NEUTRON_M_MEV,
           -3.82608545, 879.6),
  lambda("lambda", "\\Lambda", "Lam", "100", 
          Zaki::Physics::LAMBDA_ZERO_M_MEV, -1.22, 2.632e-10), 
  proton("proton", "p^+", "p", "11", Zaki::Physics::PROTON_M_MEV,
           5.5856946893, 1e+41),
  sigma_m("Sigma-", "\\Sigma^-", "Sig-", "110", Zaki::Physics::SIGMA_MINUS_M_MEV,
           -1.160, 1.479e-10),
  process(in_process),
  // m_chi_vals({{0, 1400}, 700, "Linear"})
  m_chi_vals({{0, 1100}, 220, "Linear"})
{ }

//--------------------------------------------------------------
// Destructor
BNV_Chi::~BNV_Chi() { }
//--------------------------------------------------------------
// Outputs the vacuum branching ratio of 
//      B -> chi + photon
// 'in_m_chi' & 'in_eps_lam' must be in MeV
// ____________________________________
double BNV_Chi::Vacuum_Decay_Br(const Baryon& B, 
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
  double m_B_0 = B.mass * Zaki::Physics::MEV_2_INV_FM ;
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
// Sets the EoS model
void BNV_Chi::SetModel(const std::string& in_eos_model) 
{
  model = in_eos_model ;
}

//--------------------------------------------------------------
// Sets the pulsar
void BNV_Chi::SetPulsar(const CompactStar::Pulsar& in_pulsar) 
{
  pulsar = in_pulsar ;

  BNV_Chi_pulsar_mass = pulsar.GetMass().val ;
}

//--------------------------------------------------------------
// Sets m_chi_vals
void BNV_Chi::SetChiMassRange(const Zaki::Math::Axis& in_m_chis) 
{
  m_chi_vals = in_m_chis ;
}

//--------------------------------------------------------------
// Eventually will be replaced when I switch to using 
// binary pulsars
// sets 'gamma_bnv_bin_lim' in yr^-1
void BNV_Chi::SetBNVLimit(const double& in_gamma) 
{
  gamma_bnv_bin_lim = in_gamma ;
}

//--------------------------------------------------------------
// Imports the EoS
void BNV_Chi::ImportEOS(const Zaki::String::Directory& eos_dir, 
                        const std::string& micro_model)
{
  // Zaki::Vector::DataSet ds_EOS(dir.ParentDir().ParentDir() +"/EOS/CompOSE/DS(CMF)-"+ model,
  //                             "DS(CMF)-" + model + ".eos") ;
  
  // eos.SetWrkDir(wrk_dir + "/eos/CompOSE/" + model) ;
  eos.SetWrkDir(eos_dir + model) ;
  // eos.SetWrkDir(eos_dir) ;
  eos.SetName(model) ;

  // Imports everything at once!
  eos.ImportEOS(model) ;

  // eos.ImportGrid("eos.nb") ;
  // eos.ImportThermo("eos.thermo") ;
  // eos.ImportCompo("eos.compo") ;

  // if (micro_model != "")
  // {
  //   eos.ExtendMeffToCrust(eos_dir + micro_model, micro_model + "_m_eff.micro") ;
  //   eos.ExtendVeffToCrust(eos_dir + micro_model, micro_model + "_V.micro") ;
  // }
  // else
  // {
  //   eos.ImportMicro("eos.micro") ;
  // }
  
  // EOS_ds.SetWrkDir(wrk_dir) ;
  // EOS_ds.Import("/EOS/CompOSE/"+ model + "/" + model + ".eos") ;
  
  // std::cout << "\n\t -> " << eos.GetMeff(1).label << "\n" ;

  sig0_ds = eos.GetVeff() ;
  m_B_ds = eos.GetMeff() ;

  // Finding the baryon's individual density
  // as a function of density
  n_B = { eos.GetEOS(2), 
             eos.GetEOS(2) * eos.GetEOS(neutron.label),
             eos.GetEOS(2) * eos.GetEOS(lambda.label) } ;
  n_B[0].label = "n_tot" ;
  n_B[1].label = "10" ;
  n_B[2].label = "100" ;
}

//--------------------------------------------------------------
void BNV_Chi::GenSequence() const
{
  CompactStar::TOVSolver solver ;
  solver.ImportEOS(eos.GetWrkDir() + "/" + model + ".eos") ;
  solver.SetWrkDir( wrk_dir ) ;
  solver.SetRadialRes(10000) ;
  solver.SetMaxRadius(15) ;
  solver.SetProfilePrecision(12) ;
  // solver.SetRadialScale("Linear") ;

  // double e_c_min = eos.GetEOS(0).Min()*1.01 ;
  double e_c_min = eos.GetEOS(0).Max()*0.04 ;
  double e_c_max = eos.GetEOS(0).Max()*0.9999 ;

  double m_best = 0 ;
  double m_goal = pulsar.GetMass().val ;
  // double precision = 0.0005 ;
  double precision = 1.4e-5 ;
  
  ((Zaki::String::Directory) wrk_dir + "/" + model).Create() ;

  solver.Solve( {{e_c_min, e_c_max}, 100, "Log"}, 
                model + "/" + pulsar.GetName(),  
                model) ;

  Zaki::Vector::DataSet tmp_seq(wrk_dir + model + "/" + pulsar.GetName(), 
                                model + "_Sequence.tsv") ;
  // tmp_seq[1].MaxIdx()
  
  int m_idx = tmp_seq[1].GetSubSet(0, tmp_seq[1].MaxIdx()).GetClosestIdx(m_goal) ;
  int upper_idx = m_idx + 1;
  if (upper_idx > tmp_seq[1].MaxIdx())
  {
    upper_idx = tmp_seq[1].MaxIdx() ;
    m_goal = tmp_seq[1].GetSubSet(0, tmp_seq[1].MaxIdx())[-1]*0.99999 ;
    BNV_Chi_pulsar_mass = m_goal ;
  }
  
  e_c_min = tmp_seq[0][m_idx-1] ;
  e_c_max = tmp_seq[0][upper_idx] ;
  std::cout << "\n\n e_c_min = " << e_c_min << "\n" ;
  std::cout << " e_c_max = " << e_c_max << "\n" ;

  solver.ClearSequence() ;

  std::cout << " M_goal = " << m_goal << "\n" ;

  while ( abs(m_best - m_goal) > precision ) 
  {
    solver.Solve( {{e_c_min, e_c_max}, 10, "Log"}, 
                  model + "/" + pulsar.GetName(), model) ;
    
    tmp_seq.Import(model + "_Sequence.tsv") ;

    std::cout << "\n M_best = " << m_best << ", M_goal = " << m_goal <<  "\n" ;

    m_idx = tmp_seq[1].GetClosestIdx(m_goal) ;
    
    std::cout << "\n m_idx = " << m_idx << "\n" ;

    m_best = tmp_seq[1][m_idx] ;
    e_c_min = tmp_seq[0][m_idx-1] ;
    e_c_max = tmp_seq[0][m_idx+1] ;

    solver.ClearSequence() ;
  }

  solver.AddNCondition(MassCondition) ;
  solver.SetRadialRes(100000) ;

  solver.Solve( {{e_c_min, e_c_max}, 5, "Linear"}, 
                model + "/" + pulsar.GetName(),  
                model) ;
}

//--------------------------------------------------------------
// Attaches the pulsar
void BNV_Chi::FindPulsar(const bool& gen_plots) 
{ 
  bool has_hyperons = true ;
  pulsar.SetWrkDir( wrk_dir + "/" + model + "/" + pulsar.GetName()) ;

  std::string file_path = (wrk_dir + "/" + model + "/" + pulsar.GetName() + ".tsv").Str() ;
  if (std::filesystem::exists(file_path)) 
  {
    pulsar.ImportProfile(model);
  } else 
  {
    pulsar.FindProfile(model);
  }


  if (gen_plots) // It doesn't work when we have no hyperons!
  {
    pulsar.PlotRelativeComposition(pulsar.GetName()  + "/" + pulsar.GetName() 
                                    + "_RelComp_vs_R.pdf") ;
    pulsar.PlotAbsoluteComposition(pulsar.GetName()  + "/" + pulsar.GetName() 
                                    + "_AbsComp_vs_R.pdf") ;
    pulsar.PlotFermiE(pulsar.GetName()  + "/" + pulsar.GetName() 
                                    + "_EF_vs_R.pdf") ;


    Plot_Meff_Radius() ;
    Plot_RestEnergy_Radius() ;
    Plot_EF_Radius() ;
    Plot_RestE_EF_Radius(neutron) ;
    Plot_RestE_EF_Radius(proton) ;
    Plot_CM_E_Radius(neutron) ;

    if (has_hyperons)
      Plot_CM_E_Radius(lambda) ;
    Plot_CM_E_Radius(proton) ;
    Plot_Estar_Radius(neutron) ;
  }

  // pulsar.Set_T_Core(1e10) ;
  // std::cout << "\n Get_Redshifted_T = " << pulsar.Get_Redshifted_T() << "\n" ;
  // std::cout << "\n Get_Surface_Photon_Lumin = " << pulsar.Get_Surface_Photon_Lumin() ;
  // std::cout << "\n Get_Surface_Photon_Lumin = " 
  //           << pulsar.Get_Surface_Photon_Lumin_T9242() * pow(pulsar.Get_Redshifted_T()/1e9, 2.42) ;

  // pulsar.Set_T_Surf_Apparent(1e7) ;
  // pulsar.Set_T_Core(1e10) ;

  // std::cout << "\n\t R_s - R_b = " << pulsar.GetProfile()->operator[](0).Max()
  //                                     - pulsar.Get_R_Blanket() << " " ;
  // std::cout << "\n\t T_s_infty = " << pulsar.Get_T_Surf_Apparent() << " " ;
  // std::cout << "\n\t T_s = " << pulsar.Get_T_Surf() << " " ;
  // std::cout << "\n\t T_core = " << pulsar.Get_T_Core() << " " ;
  // std::cout << "\n\t T_b = " << pulsar.Get_T_Blanket() << " " ;
  // std::cout << "\n\t L_gamma = " << pulsar.Get_Surface_Photon_Lumin() << " " ;
  // std::cout << "\n\t R_Durca = " << pulsar.Get_R_Durca() << " " ;
  // std::cout << "\n\t DUrca Emissivity = " << pulsar.Get_DUrca_Neutrino_Emissivity(0) << " " ;
  // std::cout << "\n\t MUrca Emissivity = " << pulsar.Get_MUrca_Neutrino_Emissivity(0) << " " ;
  // std::cout << "\n\t Neutrino Luminosity = " << pulsar.Get_Neutrino_Lumin() << " " ;

  // std::cout << "\n\t g_14 = " << pulsar.Get_Surface_Gravity() << " \n\n" ;
  
  // pulsar.Plot_Neutrino_Emissivity(pulsar.GetName()  + "/" + pulsar.GetName() 
  //                                   + "_nu_emissivity_vs_R.pdf") ;

  // pulsar.Plot_Neutrino_Rate(pulsar.GetName()  + "/" + pulsar.GetName() 
  //                                   + "_nu_rate_vs_R.pdf") ;

  Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;
  Zaki::Vector::DataColumn r   = pulsar.GetProfile()->operator[](0) ;

  if (has_hyperons) 
  {
    // micro dataset as a function of density
    Zaki::Vector::DataSet micro_n({ m_B_ds[0], 
                                    n_B[neutron.label],
                                    n_B[lambda.label],
                                    m_B_ds[neutron.label],
                                    m_B_ds[lambda.label],
                                    sig0_ds[neutron.label],
                                    sig0_ds[lambda.label] }) ;
    micro_n.Interpolate(0, {1, 2, 3, 4, 5, 6}) ;

    // [0]: r[km],  [1]: n_n, [2]: n_lam
    // [3]: m_n, [4]: m_lam, [5]: sig_n, [6]: sig_lam
    micro_r = {r, micro_n.Evaluate(1, n_r), 
                micro_n.Evaluate(2, n_r), 
                micro_n.Evaluate(3, n_r),
                micro_n.Evaluate(4, n_r),
                micro_n.Evaluate(5, n_r),
                micro_n.Evaluate(6, n_r)
                                    } ;
    micro_r[0].label = "R [km]" ;
    micro_r[1].label = "n_" + neutron.label ;
    micro_r[2].label = "n_" + lambda.label ;
    micro_r[3].label = "m_" + neutron.label ;
    micro_r[4].label = "m_" + lambda.label ;
    micro_r[5].label = "sig_" + neutron.label ;
    micro_r[6].label = "sig_" + lambda.label ;
  }
  else
  {
    // micro dataset as a function of density
    Zaki::Vector::DataSet micro_n({ m_B_ds[0], 
                                    n_B[neutron.label],
                                    // n_B[lambda.label],
                                    m_B_ds[neutron.label],
                                    // m_B_ds[lambda.label],
                                    sig0_ds[neutron.label],
                                    // sig0_ds[lambda.label] 
                                    }) ;
    micro_n.Interpolate(0, {1, 2, 3}) ;
    micro_n.Evaluate(1, n_r) ;
    std::cout << "\n\n Inside else 1\n\n" ;

    // [0]: r[km],  [1]: n_n,
    // [2]: m_n, [3]: sig_n
    micro_r = {r, micro_n.Evaluate(1, n_r), 
                micro_n.Evaluate(2, n_r), 
                micro_n.Evaluate(3, n_r)
                // micro_n.Evaluate(4, n_r),
                // micro_n.Evaluate(5, n_r),
                // micro_n.Evaluate(6, n_r)
                                    } ;
    std::cout << "\n\n Inside else 2\n\n" ;

    micro_r[0].label = "R [km]" ;
    micro_r[1].label = "n_" + neutron.label ;
    // micro_r[2].label = "n_" + lambda.label ;
    micro_r[2].label = "m_" + neutron.label ;
    // micro_r[4].label = "m_" + lambda.label ;
    micro_r[3].label = "sig_" + neutron.label ;
    // micro_r[6].label = "sig_" + lambda.label ;
    std::cout << "\n\n Inside else 3 \n\n" ;
  }
}

//--------------------------------------------------------------
double BNV_Chi::PhaseSpace_Integrand(const double& x) 
{ 
  return 0 ;
}

//--------------------------------------------------------------
double BNV_Chi::PhaseSpace_Integral( const Baryon& B,
                            const double& m, 
                            const double& m_chi, 
                            const double& Sigma_0, 
                            const double& x_F)
{
  return 0 ;
}

//--------------------------------------------------------------
// // Phase space integrand function (dim-less part)
// // The total phase space integrand is the output of this function
// // multiplied by:
// //                  m^4 C / (pi^2)
// // Here sigma_0 = - Sigma_0 / m.
// //
// double BNV_Chi::PhaseSpace_Integrand(const double& x) 
// {
//   double mu_chi = phase_int_mu_chi ; 
//   double sigma_0 = phase_int_sigma_0 ;

//   if ( x <= ( mu_chi*mu_chi - 1 - sigma_0*sigma_0 ) / ( 2*sigma_0 )  ) 
//   {
//     return 0 ;
//   }

//   // Old result [ WRONG ]
//   // double out = x * sqrt(x*x - 1) ;
//   //       out *= pow(1 + sigma_0*sigma_0 + 2*x*sigma_0 - mu_chi*mu_chi, 3);
//   //       out /= (x + sigma_0) ;
//   //       out /= (1 + sigma_0*sigma_0 + 2*x*sigma_0 ) ;

//   double out = x * sqrt(x*x - 1) ;
//         out *= 1 + sigma_0*sigma_0 + 2*x*sigma_0 - mu_chi*mu_chi ;
//         out *= (1 + sigma_0*sigma_0 + 2*x*sigma_0) * (1 + x*sigma_0 + 2 * mu_chi) 
//                 + pow(mu_chi, 2) * (1 + x*sigma_0) ;
//         out /= x + sigma_0 ;
//         out /= 1 + x*sigma_0 ;
//         out /= 1 + sigma_0*sigma_0 + 2*x*sigma_0 ;

//   return out ;
// }

// //--------------------------------------------------------------
// // The full phase-space integral
// // The output is in MeV^4
// // Here Sigma_0 is the vector self-energy ( Sigma_0 < 0 )
// double BNV_Chi::PhaseSpace_Integral(const Baryon& B,
//                                               const double& m, 
//                                               const double& m_chi, 
//                                               const double& Sigma_0, 
//                                               const double& x_F) 
// { 
//   double q_e = Zaki::Physics::Q_E ;
//   double eps = 1e-10 ;
  
//   phase_int_mu_chi = m_chi / m ; 
//   phase_int_sigma_0 = -Sigma_0 / m ;

//   gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
//   double err, result;

//   Zaki::Math::GSLFuncWrapper<BNV_Chi, double (BNV_Chi::*)(const double&)> 
//     Fp(this, &BNV_Chi::PhaseSpace_Integrand);     

//   gsl_function F = static_cast<gsl_function> (Fp) ; 

//   gsl_integration_qag(&F, 1, x_F, 1e-10, 1e-10, 2000, 1, w, &result, &err) ;
//   gsl_integration_workspace_free(w);

//   // C
//   double C = pow(B.g * eps * q_e, 2) / (128 * M_PI * m * m ) ;

//   result *= C ;
//   result *= pow(m, 4) / (M_PI*M_PI)  ;

//   return result ;
// }

//--------------------------------------------------------------
// Plots the effective mass as a function of radius
void BNV_Chi::Plot_Meff_Radius()
{
  // std::cout << "\n\n \t eos.GetEOS(2).Size() = " << eos.GetEOS(2).Size() << "\n\n" ;
  // std::cout << "\n\n \t m_B_ds[neutron.label].Size() = " << m_B_ds[neutron.label].Size() << "\n\n" ;
  Zaki::Vector::DataSet ds_meff_n ({eos.GetEOS(2), m_B_ds[neutron.label], 
                                      m_B_ds[lambda.label], 
                                      m_B_ds[proton.label],
                                      m_B_ds[sigma_m.label]}) ;
  ds_meff_n.Interpolate(0, {1,2,3,4}) ;

  Zaki::Vector::DataColumn dc_neutron_meff_r = 
    ds_meff_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;

  Zaki::Vector::DataColumn dc_lambda_meff_r = 
    ds_meff_n.Evaluate(2, pulsar.GetProfile()->operator[](5)) ;

  Zaki::Vector::DataColumn dc_proton_meff_r = 
    ds_meff_n.Evaluate(3, pulsar.GetProfile()->operator[](5)) ;

  Zaki::Vector::DataColumn dc_sigma_meff_r = 
    ds_meff_n.Evaluate(4, pulsar.GetProfile()->operator[](5)) ;

  Zaki::Vector::DataSet ds_meff_r({ pulsar.GetProfile()->operator[](0), 
                  dc_neutron_meff_r, dc_lambda_meff_r, dc_proton_meff_r, 
                  dc_sigma_meff_r}) ;

  std::cout << "\n\t m_eff[neutron, r=0] = " << ds_meff_r[1][0] << " " ;
  std::cout << "\n\t m_eff[lambda, r=0] = "  << ds_meff_r[2][0] << "\n" ;

  // -------------------------------------------------------
  // So the issue is that this will include radii where
  // there isn't any lambda or sigma-, etc. To fix this
  //  we set a condition to cut the curves where their 
  //  abundance falls below 1e-8 fm^{-3}.
  // Loop over columns
  for (size_t c = 1; c < ds_meff_r.Dim().size(); c++)
  {
    for (size_t r = 0 ; r < ds_meff_r[c].Size() ; r++)
    {
      if ( pulsar.GetProfile()->operator[](ds_meff_n[c].label)[r] < 1e-8 )
      {
        ds_meff_r[c][r] = 0 ;
      }
    }
  }
  // -------------------------------------------------------

  ds_meff_r.SetWrkDir(pulsar.GetWrkDir() + "/" + pulsar.GetName() ) ;
  
  auto non_zero_cond = [](const double& v){ return v > 0 ;} ;

  // Finding the lowest and highest points on the curves
  double Y_min = std::min({ds_meff_r[1].GetSubSet(non_zero_cond).Min(), 
                            ds_meff_r[2].GetSubSet(non_zero_cond).Min(), 
                            ds_meff_r[3].GetSubSet(non_zero_cond).Min(), 
                            ds_meff_r[4].GetSubSet(non_zero_cond).Min()}) ;
  double Y_max = std::max({ds_meff_r[1].Max(), ds_meff_r[2].Max(), 
                            ds_meff_r[3].Max(), ds_meff_r[4].Max()}) ;


  Zaki::Vector::DataSet::PlotParam plt_par ;
  // plt_par.SetXAxis({ds_meff_r[0].Min()*0.99, ds_meff_r[0].Max()*1.01}) ;
  plt_par.SetXAxis({0, 14}) ;
  
  // plt_par.SetYAxis({ds_meff_r[1].Min()*0.98, ds_meff_r[2].Max()*1.02}) ;
  plt_par.SetYAxis({Y_min*0.98, Y_max*1.02}) ;
  plt_par.SetXAxisLabel("$R\\, [\\, {\\rm km}\\,]$") ;
  plt_par.SetYAxisLabel("$m_B^*\\, [\\, {\\rm MeV}\\,]$") ;
  plt_par.SetGrid() ;
  plt_par.SetLegend({"upper left", 0.0, 1.0}) ;

  // This is for setting Y-axis ticks
  int min_tick = Y_min*0.98 - fmod(Y_min*0.98, 50) ;
  int max_tick = Y_max*1.02 - fmod(Y_max*1.02, 50) ;
  std::vector<double> y_tick_set ;
  for (int i = min_tick; i <= max_tick; i+= 100)
  {
    y_tick_set.emplace_back(i) ;
  }
  
  // plt_par.SetYTicks({{900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350}}) ;
  plt_par.SetYTicks({y_tick_set}) ;

  ds_meff_r.SetPlotPars(plt_par) ;

  ds_meff_r.Plot(0, {{1, "$m^*_n$"} , {2, "$m^*_{\\Lambda}$"}, {3, "$m^*_{p}$"}, 
                     {4, "$m^*_{\\Sigma}$"}}, 
                 pulsar.GetName() + "_m_vs_R.pdf", pulsar.GetName() + "\n" + model) ;
}

//--------------------------------------------------------------
// Plots the baryon's rest-energy (m_B + Sigma_0) 
//  as a function of radius
void BNV_Chi::Plot_RestEnergy_Radius()
{
  // Building the rest energy dataset as a function of density
  Zaki::Vector::DataSet ds_E_n ({eos.GetEOS(2), 
            m_B_ds[neutron.label] + sig0_ds[neutron.label], 
            m_B_ds[proton.label] + sig0_ds[proton.label], 
            m_B_ds[lambda.label] + sig0_ds[lambda.label],
            m_B_ds[sigma_m.label] + sig0_ds[sigma_m.label]}) ;

  ds_E_n.Interpolate(0, {1,2,3,4}) ;
  
  ds_E_n[1].label = neutron.label ;
  ds_E_n[2].label = proton.label ;
  ds_E_n[3].label = lambda.label ;
  ds_E_n[4].label = sigma_m.label ;
  // ......................................
  // Building the rest energy data-columns as a function of radius
  Zaki::Vector::DataColumn dc_neutron_E_r = 
    ds_E_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;
  
  Zaki::Vector::DataColumn dc_proton_E_r = 
    ds_E_n.Evaluate(2, pulsar.GetProfile()->operator[](5)) ;

  Zaki::Vector::DataColumn dc_lambda_E_r = 
    ds_E_n.Evaluate(3, pulsar.GetProfile()->operator[](5)) ;
  
  Zaki::Vector::DataColumn dc_sigma_m_E_r = 
    ds_E_n.Evaluate(4, pulsar.GetProfile()->operator[](5)) ;
  // ......................................

  // ......................................
  //  Combining the data-columns into a dataset
  Zaki::Vector::DataSet ds_E_r({ pulsar.GetProfile()->operator[](0), 
                      dc_neutron_E_r, dc_proton_E_r, 
                      dc_lambda_E_r, dc_sigma_m_E_r}) ;

  std::cout << "\n\t E_rest[neutron, r=0] = " << ds_E_r[1][0] << " " ;
  std::cout << "\n\t E_rest[lambda, r=0] = "  << ds_E_r[3][0] << "\n" ;
  // ......................................
  // Fixing the labels
  ds_E_r[1].label = "$E_n^0$" ;
  ds_E_r[2].label = "$E_{p}^0$" ;
  ds_E_r[3].label = "$E_{\\Lambda}^0$" ;
  ds_E_r[4].label = "$E_{\\Sigma^-}^0$" ;
  // ......................................

  // -------------------------------------------------------
  // So the issue is that this will include radii where
  // there isn't any lambda or sigma-, etc. To fix this
  //  we set a condition to cut the curves where their 
  //  abundance falls below 1e-8 fm^{-3}.
  // Loop over columns
  for (size_t c = 1; c < ds_E_r.Dim().size(); c++)
  {
    for (size_t r = 0 ; r < ds_E_r[c].Size() ; r++)
    {
      if ( pulsar.GetProfile()->operator[](ds_E_n[c].label)[r] < 1e-8 )
      {
        ds_E_r[c][r] = 0 ;
      }
    }
  }
  // -------------------------------------------------------

  ds_E_r.SetWrkDir(pulsar.GetWrkDir() + "/" + pulsar.GetName() ) ;
  
  Zaki::Vector::DataSet::PlotParam plt_par ;

  // Setting the boundaries
  // plt_par.SetXAxis({ds_E_r[0].Min(), ds_E_r[0].Max()}) ;
  plt_par.SetXAxis({0, 14}) ;

  auto non_zero_cond = [](const double& v){ return v > 0 ;} ;

  // Finding the lowest and highest points on the curves
  double Y_min = std::min({ds_E_r[1].GetSubSet(non_zero_cond).Min(), 
                            ds_E_r[2].GetSubSet(non_zero_cond).Min(), 
                            ds_E_r[3].GetSubSet(non_zero_cond).Min(), 
                            ds_E_r[4].GetSubSet(non_zero_cond).Min()}) ;
  double Y_max = std::max({ds_E_r[1].Max(), ds_E_r[2].Max(), 
                            ds_E_r[3].Max(), ds_E_r[4].Max()}) ;

  plt_par.SetYAxis({Y_min*0.98, Y_max*1.02}) ;
  
  // This is for setting Y-axis ticks
  int min_tick = Y_min*0.98 - fmod(Y_min*0.98, 50) ;
  int max_tick = Y_max*1.02 - fmod(Y_max*1.02, 50) ;
  std::vector<double> y_tick_set ;
  for (int i = min_tick; i <= max_tick; i+= 100)
  {
    y_tick_set.emplace_back(i) ;
  }
  
  // plt_par.SetYTicks({{900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350}}) ;
  plt_par.SetYTicks({y_tick_set}) ;

  plt_par.SetXAxisLabel("$R\\, (\\, {\\rm km}\\,)$") ;
  plt_par.SetYAxisLabel("$E_B^0\\, (\\, {\\rm MeV}\\,)$") ;
  plt_par.SetLegend({"upper right", 1.0, 1.0}) ;
  plt_par.SetGrid() ;
  plt_par.AddAxHLine(neutron.mass, {{"c", "blue"}, {"linestyle", "--"}}) ;
  // plt_par.AddAxHLine(proton.mass) ;
  plt_par.AddAxHLine(lambda.mass, {{"c", "orange"}, {"linestyle", "--"}}) ;
  plt_par.AddAxHLine(sigma_m.mass, {{"c", "red"}, {"linestyle", "--"}}) ;

  ds_E_r.SetPlotPars(plt_par) ;

  ds_E_r.Plot(0, {1,3,2,4}, pulsar.GetName() + "_E_vs_R.pdf", 
              pulsar.GetName() + "\n" + model) ;
}

//--------------------------------------------------------------
// Plots the Fermi energy as a function of radius
void BNV_Chi::Plot_EF_Radius() 
{

  // ---------------------------------------------------
  // Finding Fermi energy as a function of density
  // ---------------------------------------------------
  Zaki::Vector::DataColumn fermi_electron =  ( 
    pow(Zaki::Physics::ELECTRON_M_FM, 2) 
    + (3*M_PI*M_PI* eos.GetEOS("0") * eos.GetEOS(2)).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM ;

  Zaki::Vector::DataColumn fermi_muon =  ( 
    pow(Zaki::Physics::MUON_M_FM, 2) 
    + (3*M_PI*M_PI* eos.GetEOS("1") * eos.GetEOS(2)).pow(2./3.)).sqrt() / Zaki::Physics::MEV_2_INV_FM ;

  Zaki::Vector::DataColumn fermi_neutron =  ( 
    m_B_ds[neutron.label].pow(2) 
    + (3*M_PI*M_PI* eos.GetEOS("10") * eos.GetEOS(2)).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() + sig0_ds[neutron.label]  ;

  Zaki::Vector::DataColumn fermi_lambda =  ( 
    m_B_ds[lambda.label].pow(2) 
    + (3*M_PI*M_PI* eos.GetEOS("100") * eos.GetEOS(2)).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() + sig0_ds[lambda.label]  ;

  Zaki::Vector::DataColumn fermi_proton =  ( 
    m_B_ds[proton.label].pow(2) 
    + (3*M_PI*M_PI* eos.GetEOS("11") * eos.GetEOS(2)).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() + sig0_ds[proton.label]  ;

  Zaki::Vector::DataSet fermi_ds_n({eos.GetEOS(2), 
                                  fermi_electron, fermi_muon, 
                                  fermi_neutron, fermi_lambda, 
                                  fermi_proton}) ;

  fermi_ds_n.Interpolate(0, {1,2,3,4,5}) ;

  // Zaki::Vector::DataSet ds_E_n ({eos.GetEOS(2), 
  //           m_B_ds[neutron.label] + sig0_ds[neutron.label], 
  //           m_B_ds[proton.label] + sig0_ds[proton.label], 
  //           m_B_ds[lambda.label] + sig0_ds[lambda.label],
  //           m_B_ds[sigma_m.label] + sig0_ds[sigma_m.label]}) ;

  // ds_E_n.Interpolate(0, {1,2,3,4}) ;
  
  fermi_ds_n[1].label = "0" ;
  fermi_ds_n[2].label = "1" ;
  fermi_ds_n[3].label = neutron.label ;
  fermi_ds_n[4].label = lambda.label ;
  fermi_ds_n[5].label = proton.label ;
  // ......................................
  // Building the Fermi energy data-columns as a function of radius
  Zaki::Vector::DataColumn dc_electron_EF_r = 
    fermi_ds_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;
  
  Zaki::Vector::DataColumn dc_muon_EF_r = 
    fermi_ds_n.Evaluate(2, pulsar.GetProfile()->operator[](5)) ;

  Zaki::Vector::DataColumn dc_neutron_EF_r = 
    fermi_ds_n.Evaluate(3, pulsar.GetProfile()->operator[](5)) ;
  
  Zaki::Vector::DataColumn dc_lambda_EF_r = 
    fermi_ds_n.Evaluate(4, pulsar.GetProfile()->operator[](5)) ;
  
  Zaki::Vector::DataColumn dc_proton_EF_r = 
    fermi_ds_n.Evaluate(5, pulsar.GetProfile()->operator[](5)) ;
  // ......................................

  // ......................................
  //  Combining the data-columns into a dataset
  Zaki::Vector::DataSet ds_EF_r({ pulsar.GetProfile()->operator[](0), 
                      dc_electron_EF_r, dc_muon_EF_r, 
                      dc_neutron_EF_r, dc_lambda_EF_r,
                      dc_proton_EF_r}) ;

  std::cout << "\n\t EF[neutron, r=0] = " << ds_EF_r[3][0] << " " ;
  std::cout << "\n\t EF[lambda, r=0] = "  << ds_EF_r[4][0] << "\n" ;
  // ......................................
  // Fixing the labels
  ds_EF_r[1].label = "$E_F(e^-)$" ;
  ds_EF_r[2].label = "$E_F(\\mu)$" ;
  ds_EF_r[3].label = "$E_F(n)$" ;
  ds_EF_r[4].label = "$E_F(\\Lambda)$" ;
  ds_EF_r[5].label = "$E_F(p^+)$" ;
  // ......................................

  // -------------------------------------------------------
  // So the issue is that this will include radii where
  // there isn't any lambda or sigma-, etc. To fix this
  //  we set a condition to cut the curves where their 
  //  abundance falls below 1e-10 fm^{-3}.
  // Loop over columns
  for (size_t c = 1; c < ds_EF_r.Dim().size(); c++)
  {
    for (size_t r = 0 ; r < ds_EF_r[c].Size() ; r++)
    {
      if ( pulsar.GetProfile()->operator[](fermi_ds_n[c].label)[r] < 1e-10 )
      {
        ds_EF_r[c][r] = 0 ;
      }
    }
  }
  // -------------------------------------------------------

  ds_EF_r.SetWrkDir(pulsar.GetWrkDir() + "/" + pulsar.GetName() ) ;
  
  Zaki::Vector::DataSet::PlotParam plt_par ;

  // Setting the boundaries
  // plt_par.SetXAxis({ds_E_r[0].Min(), ds_E_r[0].Max()}) ;
  plt_par.SetXAxis({0, 14}) ;

  auto non_zero_cond = [](const double& v){ return v > 0 ;} ;

  // Finding the lowest and highest points on the curves
  double Y_min = std::min({ds_EF_r[1].GetSubSet(non_zero_cond).Min(), 
                            ds_EF_r[2].GetSubSet(non_zero_cond).Min(), 
                            ds_EF_r[3].GetSubSet(non_zero_cond).Min(), 
                            ds_EF_r[4].GetSubSet(non_zero_cond).Min(),
                            ds_EF_r[5].GetSubSet(non_zero_cond).Min()}) ;
  double Y_max = std::max({ds_EF_r[1].Max(), ds_EF_r[2].Max(), 
                            ds_EF_r[3].Max(), ds_EF_r[4].Max(),
                            ds_EF_r[5].Max()}) ;

  // plt_par.SetYAxis({Y_min*0.98, Y_max*1.02}) ;
  plt_par.SetYAxis({0, 1400}) ;
  
  // This is for setting Y-axis ticks
  int min_tick = Y_min*0.98 - fmod(Y_min*0.98, 50) ;
  // int max_tick = Y_max*1.02 - fmod(Y_max*1.02, 50) ;
  int max_tick = 1400*1.02 - fmod(1400*1.02, 50) ;

  std::vector<double> y_tick_set ;
  for (int i = min_tick; i <= max_tick; i+= 100)
  {
    y_tick_set.emplace_back(i) ;
  }
  
  // plt_par.SetYTicks({{900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350}}) ;
  plt_par.SetYTicks({y_tick_set}) ;

  plt_par.SetXAxisLabel("$R\\, (\\, {\\rm km}\\,)$") ;
  plt_par.SetYAxisLabel("$E_F\\, (\\, {\\rm MeV}\\,)$") ;
  plt_par.SetLegend({"upper right", 1.0, 1.0}) ;
  plt_par.SetGrid() ;
  // plt_par.AddAxHLine(neutron.mass, {{"c", "blue"}, {"linestyle", "--"}}) ;
  // plt_par.AddAxHLine(proton.mass) ;
  // plt_par.AddAxHLine(lambda.mass, {{"c", "orange"}, {"linestyle", "--"}}) ;
  // plt_par.AddAxHLine(sigma_m.mass, {{"c", "red"}, {"linestyle", "--"}}) ;

  ds_EF_r.SetPlotPars(plt_par) ;

  ds_EF_r.Plot(0, {1,2,3,4,5}, pulsar.GetName() + "_EF_vs_R.pdf", 
              pulsar.GetName() + "\n" + model) ;
}

//--------------------------------------------------------------
// Plots both E* band in the nm frame as a function of radius
void BNV_Chi::Plot_Estar_Radius(const Baryon& B)
{
  // ---------------------------------------------------
  // Finding maximum E* as a function of density
  // ---------------------------------------------------
  Zaki::Vector::DataColumn Estar_max_n =  ( 
    m_B_ds[B.label].pow(2)
    + (3*M_PI*M_PI* eos.GetEOS(B.label) * eos.GetEOS(2)).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() ;

  Zaki::Vector::DataSet ds_Estar_max_n({eos.GetEOS(2), 
                                  Estar_max_n}) ;

  ds_Estar_max_n.Interpolate(0, {1}) ;

  ds_Estar_max_n[1].label = B.label ;
  // ......................................
  // Building the maximum CM energy data-columns as a function of radius
  Zaki::Vector::DataColumn dc_B_Estar_max_r = 
    ds_Estar_max_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;
  // ......................................

  // ..............................................
  // Building E* dataset as a function of density
  Zaki::Vector::DataSet ds_Estar_min_n ({eos.GetEOS(2), m_B_ds[B.label]}) ;

  ds_Estar_min_n.Interpolate(0, {1}) ;
  
  ds_Estar_min_n[1].label = B.label ;
  // ......................................
  // Building the E* data-columns as a function of radius
  Zaki::Vector::DataColumn dc_B_Estar_min_r = 
    ds_Estar_min_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;
  // ......................................

  // ......................................
  //  Combining the data-columns into a dataset
  Zaki::Vector::DataSet ds_E0_EF_r({ pulsar.GetProfile()->operator[](0), 
                                  dc_B_Estar_min_r, dc_B_Estar_max_r}) ;
  // ......................................
  // Fixing the labels
  ds_E0_EF_r[1].label = "$E_{\\rm min}^{*}("+B.TeX_name+")$" ;
  ds_E0_EF_r[2].label = "$E_{\\rm max}^{*}("+B.TeX_name+")$" ;
  // ......................................

  ds_E0_EF_r.SetWrkDir(pulsar.GetWrkDir() + "/" + pulsar.GetName() ) ;
  
  Zaki::Vector::DataSet::PlotParam plt_par ;

  // Setting the boundaries
  // plt_par.SetXAxis({ds_E_r[0].Min(), ds_E_r[0].Max()}) ;
  plt_par.SetXAxis({0, 14}) ;

  // auto non_zero_cond = [](const double& v){ return v > 0 ;} ;

  // plt_par.SetYAxis({800, 1325}) ; // proton
  plt_par.SetYAxis({200, 1000}) ; // neutron

  plt_par.SetXAxisLabel("$R\\, (\\, {\\rm km}\\,)$") ;
  plt_par.SetYAxisLabel("$E\\, (\\, {\\rm MeV}\\,)$") ;
  plt_par.SetLegend({"upper right", 1.0, 1.0}) ;
  plt_par.SetGrid() ;
  // plt_par.AddAxHLine(neutron.mass, {{"c", "blue"}, {"linestyle", "--"}}) ;
  // plt_par.AddAxHLine(B.mass, {{"linestyle", "--"}}) ;
  // plt_par.AddAxHLine(1078, {{"linestyle", "--"}, {"color", "red"}}) ;

  // plt_par.AddAxHLine(lambda.mass, {{"c", "orange"}, {"linestyle", "--"}}) ;
  // plt_par.AddAxHLine(sigma_m.mass, {{"c", "red"}, {"linestyle", "--"}}) ;

  ds_E0_EF_r.SetPlotPars(plt_par) ;

  ds_E0_EF_r.Plot(0, {1,2}, pulsar.GetName() + "_Estar_vs_R_"+B.short_name+".pdf", 
              pulsar.GetName() + "\n" + model) ;
}

//--------------------------------------------------------------
// Plots both the CM frame energy as a function of radius
void BNV_Chi::Plot_CM_E_Radius(const Baryon& B)
{
  // ---------------------------------------------------
  // Finding maximum CM energy as a function of density
  // ---------------------------------------------------
  // Zaki::Vector::DataColumn ECM_max_n =  ( 
  //   m_B_ds[B.label].pow(2)
  //   + (3*M_PI*M_PI* eos.GetEOS(B.label) * eos.GetEOS(2)).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
  //   ).sqrt() + sig0_ds[B.label]  ;

  Zaki::Vector::DataColumn ECM_max_n =  ( 
    m_B_ds[B.label].pow(2)
    + (3*M_PI*M_PI* eos.GetEOS(B.label) * eos.GetEOS(2)).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() * 2 * sig0_ds[B.label]  ;

  ECM_max_n += sig0_ds[B.label].pow(2) ; 
  ECM_max_n += m_B_ds[B.label].pow(2) ; 

  ECM_max_n = ECM_max_n.sqrt() ;

  Zaki::Vector::DataSet ds_ECM_max_n({eos.GetEOS(2), 
                                  ECM_max_n}) ;

  ds_ECM_max_n.Interpolate(0, {1}) ;

  ds_ECM_max_n[1].label = B.label ;
  // ......................................
  // Building the maximum CM energy data-columns as a function of radius
  Zaki::Vector::DataColumn dc_B_ECM_max_r = 
    ds_ECM_max_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;
  // ......................................

  // ..............................................
  // Building the rest energy dataset as a function of density
  Zaki::Vector::DataSet ds_ECM_min_n ({eos.GetEOS(2), 
            m_B_ds[B.label] + sig0_ds[B.label]}) ;

  ds_ECM_min_n.Interpolate(0, {1}) ;
  
  ds_ECM_min_n[1].label = B.label ;
  // ......................................
  // Building the rest energy data-columns as a function of radius
  Zaki::Vector::DataColumn dc_B_ECM_min_r = 
    ds_ECM_min_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;
  // ......................................

  //  Combining the density columns into a dataset
  Zaki::Vector::DataSet ds_E0_EF_n({ ds_ECM_min_n[0], 
                                  ds_ECM_min_n[1], ds_ECM_max_n[1]}) ;
  ds_E0_EF_n.SetWrkDir(pulsar.GetWrkDir() + "/" + pulsar.GetName() ) ;

  ds_E0_EF_n.Export("ECM_vs_n_"+B.short_name+".tsv") ;

  // ......................................
  //  Combining the radius columns into a dataset
  Zaki::Vector::DataSet ds_E0_EF_r({ pulsar.GetProfile()->operator[](0), 
                                  dc_B_ECM_min_r, dc_B_ECM_max_r}) ;
  // ......................................
  // Fixing the labels
  ds_E0_EF_r[1].label = "$E_{\\rm cm}^{\\rm min}("+B.TeX_name+")$" ;
  ds_E0_EF_r[2].label = "$E_{\\rm cm}^{\\rm max}("+B.TeX_name+")$" ;
  // ......................................

  ds_E0_EF_r.SetWrkDir(pulsar.GetWrkDir() + "/" + pulsar.GetName() ) ;
  
  Zaki::Vector::DataSet::PlotParam plt_par ;

  // Setting the boundaries
  // plt_par.SetXAxis({ds_E_r[0].Min(), ds_E_r[0].Max()}) ;
  plt_par.SetXAxis({0, 14}) ;

  // auto non_zero_cond = [](const double& v){ return v > 0 ;} ;

  // plt_par.SetYAxis({800, 1325}) ; // proton
  plt_par.SetYAxis({800, 1450}) ; // neutron

  plt_par.SetXAxisLabel("$R\\, (\\, {\\rm km}\\,)$") ;
  plt_par.SetYAxisLabel("$E\\, (\\, {\\rm MeV}\\,)$") ;
  plt_par.SetLegend({"upper right", 1.0, 1.0}) ;
  plt_par.SetGrid() ;
  // plt_par.AddAxHLine(neutron.mass, {{"c", "blue"}, {"linestyle", "--"}}) ;
  plt_par.AddAxHLine(B.mass, {{"linestyle", "--"}}) ;
  // plt_par.AddAxHLine(1078, {{"linestyle", "--"}, {"color", "red"}}) ;
  plt_par.AddAxHLine(500+140, {{"linestyle", "--"}, {"color", "red"}}) ;
  plt_par.AddAxHLine(1000+140, {{"linestyle", "--"}, {"color", "red"}}) ;
  plt_par.AddAxHLine(1200+140, {{"linestyle", "--"}, {"color", "red"}}) ;

  // plt_par.AddAxHLine(lambda.mass, {{"c", "orange"}, {"linestyle", "--"}}) ;
  // plt_par.AddAxHLine(sigma_m.mass, {{"c", "red"}, {"linestyle", "--"}}) ;

  ds_E0_EF_r.SetPlotPars(plt_par) ;

  ds_E0_EF_r.Plot(0, {1,2}, pulsar.GetName() + "_ECM_vs_R_"+B.short_name+".pdf", 
              pulsar.GetName() + "\n" + model) ;

  ds_E0_EF_r.Export(pulsar.GetName() + "_ECM_vs_R_"+B.short_name+".tsv") ;
}

//--------------------------------------------------------------
// Plots both the rest energy and the Fermi energy
void BNV_Chi::Plot_RestE_EF_Radius(const Baryon& B)
{

  // ---------------------------------------------------
  // Finding Fermi energy as a function of density
  // ---------------------------------------------------
  Zaki::Vector::DataColumn fermi_B =  ( 
    m_B_ds[B.label].pow(2) 
    + (3*M_PI*M_PI* eos.GetEOS(B.label) * eos.GetEOS(2)).pow(2./3.) / pow(Zaki::Physics::MEV_2_INV_FM, 2)
    ).sqrt() + sig0_ds[B.label]  ;

  Zaki::Vector::DataSet fermi_ds_n({eos.GetEOS(2), 
                                  fermi_B}) ;

  fermi_ds_n.Interpolate(0, {1}) ;

  fermi_ds_n[1].label = B.label ;
  // ......................................
  // Building the Fermi energy data-columns as a function of radius
  Zaki::Vector::DataColumn dc_B_EF_r = 
    fermi_ds_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;
  // ......................................

  // ..............................................
  // Building the rest energy dataset as a function of density
  Zaki::Vector::DataSet ds_E_n ({eos.GetEOS(2), 
            m_B_ds[B.label] + sig0_ds[B.label]}) ;

  ds_E_n.Interpolate(0, {1}) ;
  
  ds_E_n[1].label = B.label ;
  // ......................................
  // Building the rest energy data-columns as a function of radius
  Zaki::Vector::DataColumn dc_B_E_r = 
    ds_E_n.Evaluate(1, pulsar.GetProfile()->operator[](5)) ;
  // ......................................

  // ......................................
  //  Combining the data-columns into a dataset
  Zaki::Vector::DataSet ds_E0_EF_r({ pulsar.GetProfile()->operator[](0), 
                                  dc_B_E_r, dc_B_EF_r}) ;
  // ......................................
  // Fixing the labels
  ds_E0_EF_r[1].label = "$E_0("+B.TeX_name+")$" ;
  ds_E0_EF_r[2].label = "$E_F("+B.TeX_name+")$" ;
  // ......................................

  // -------------------------------------------------------
  // So the issue is that this will include radii where
  // there isn't any lambda or sigma-, etc. To fix this
  //  we set a condition to cut the curves where their 
  //  abundance falls below 1e-10 fm^{-3}.
  // Loop over columns
  // for (size_t c = 1; c < ds_EF_r.Dim().size(); c++)
  // for (size_t c = 1; c < 5; c++)
  // {
  //   for (size_t r = 0 ; r < ds_EF_r[c].Size() ; r++)
  //   {
  //     if ( pulsar.GetProfile()->operator[](fermi_ds_n[c].label)[r] < 1e-10 )
  //     {
  //       ds_EF_r[c][r] = 0 ;
  //     }
  //   }
  // }
  // -------------------------------------------------------

  ds_E0_EF_r.SetWrkDir(pulsar.GetWrkDir() + "/" + pulsar.GetName() ) ;
  
  Zaki::Vector::DataSet::PlotParam plt_par ;

  // Setting the boundaries
  // plt_par.SetXAxis({ds_E_r[0].Min(), ds_E_r[0].Max()}) ;
  plt_par.SetXAxis({0, 14}) ;

  auto non_zero_cond = [](const double& v){ return v > 0 ;} ;

  // Finding the lowest and highest points on the curves
  // double Y_min = std::min({ds_EF_r[1].GetSubSet(non_zero_cond).Min(), 
  //                           ds_EF_r[2].GetSubSet(non_zero_cond).Min(), 
  //                           ds_EF_r[3].GetSubSet(non_zero_cond).Min(), 
  //                           ds_EF_r[4].GetSubSet(non_zero_cond).Min(),
  //                           ds_EF_r[5].GetSubSet(non_zero_cond).Min(),
  //                           ds_EF_r[6].GetSubSet(non_zero_cond).Min(),
  //                           ds_EF_r[7].GetSubSet(non_zero_cond).Min(),
  //                           ds_EF_r[8].GetSubSet(non_zero_cond).Min(),
  //                           ds_EF_r[9].GetSubSet(non_zero_cond).Min()}) ;
  // double Y_max = std::max({ds_EF_r[1].Max(), ds_EF_r[2].Max(), 
  //                           ds_EF_r[3].Max(), ds_EF_r[4].Max(),
  //                           ds_EF_r[5].Max(), ds_EF_r[6].Max(),
  //                           ds_EF_r[7].Max(), ds_EF_r[8].Max(),
  //                           ds_EF_r[9].Max()}) ;

  // plt_par.SetYAxis({Y_min*0.98, Y_max*1.02}) ;
  plt_par.SetYAxis({800, 1325}) ;
  
  // This is for setting Y-axis ticks
  // int min_tick = Y_min*0.98 - fmod(Y_min*0.98, 50) ;
  // int max_tick = Y_max*1.02 - fmod(Y_max*1.02, 50) ;
  // int max_tick = 1400*1.02 - fmod(1400*1.02, 50) ;

  // std::vector<double> y_tick_set ;
  // for (int i = min_tick; i <= max_tick; i+= 100)
  // {
  //   y_tick_set.emplace_back(i) ;
  // }
  
  // plt_par.SetYTicks({{900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350}}) ;
  // plt_par.SetYTicks({y_tick_set}) ;

  plt_par.SetXAxisLabel("$R\\, (\\, {\\rm km}\\,)$") ;
  plt_par.SetYAxisLabel("$E\\, (\\, {\\rm MeV}\\,)$") ;
  plt_par.SetLegend({"upper right", 1.0, 1.0}) ;
  plt_par.SetGrid() ;
  // plt_par.AddAxHLine(neutron.mass, {{"c", "blue"}, {"linestyle", "--"}}) ;
  plt_par.AddAxHLine(B.mass, {{"linestyle", "--"}}) ;
  // plt_par.AddAxHLine(lambda.mass, {{"c", "orange"}, {"linestyle", "--"}}) ;
  // plt_par.AddAxHLine(sigma_m.mass, {{"c", "red"}, {"linestyle", "--"}}) ;

  ds_E0_EF_r.SetPlotPars(plt_par) ;

  ds_E0_EF_r.Plot(0, {1,2}, pulsar.GetName() + "_EBand_vs_R_"+B.short_name+".pdf", 
              pulsar.GetName() + "\n" + model) ;
}

//--------------------------------------------------------------
// /// Evaluates the total decay rate for the pulsar
// double
// BNV_Chi::EvalDecayRate(const double& m_chi, const Baryon& in_B)
// {
//   std::cout << "\n -------------------- m_chi = " << m_chi 
//             <<  " MeV ---- B = " << in_B.label <<"-----------------\n" ;

//   Zaki::Vector::DataColumn B_chi_photon_rate =  
//         B_Chi_Photon_Rate(m_chi, in_B) ;

//   Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;

//   Zaki::Vector::DataSet ds_B_chi_photon_rate_n ({eos.GetEOS(2), B_chi_photon_rate}) ;
//   ds_B_chi_photon_rate_n.Interpolate(0, 1) ;
  
//   Zaki::Vector::DataColumn dc_B_chi_photon_rate_r = 
//     ds_B_chi_photon_rate_n.Evaluate(1, n_r) ;


//   // Converting "MeV^3" into "1/fm^3"
//   dc_B_chi_photon_rate_r *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

//   // Converting "1/fm^3" into "1/km^3"
//   dc_B_chi_photon_rate_r *= 1e+54 ; 

//   // Converting "MeV" into "s^-1"
//   dc_B_chi_photon_rate_r *= 1e-3 * Zaki::Physics::GEV_2_S ; 

//   Zaki::Vector::DataColumn r = pulsar.GetProfile()->operator[](0) ;
  
//   //----------------------------------------------------------------------
//   //          Plotting decay per V [km^-3.s^-1] vs R [km]
//   //----------------------------------------------------------------------
//   if (false)
//   {
//     Zaki::Vector::DataSet ds_B_chi_photon_rate_r ({r, dc_B_chi_photon_rate_r}) ;
//     ds_B_chi_photon_rate_r[1].label = "$\\Gamma ({\\rm km}^{-3} \\cdot s^{-1})$";
//     char tmp_char[50] ;
//     snprintf(tmp_char, sizeof(tmp_char), "%.0f", m_chi) ;
//     ds_B_chi_photon_rate_r.SetWrkDir(wrk_dir + "/BNV_2022/results/B_Chi_Transition/"+model) ;
//     ds_B_chi_photon_rate_r.Plot(0, 1, "Decay_per_V_vs_R/" + in_B.label +"_Decay_per_V_vs_R/" + 
//                         in_B.label+"_decay_per_V_vs_R_"+std::string(tmp_char)+".pdf",
//                         "$\\varepsilon = 10^{-10}\\, {\\rm MeV}$") ;
//   }
//   //----------------------------------------------------------------------
  
//   Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
//   Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

//   Zaki::Vector::DataColumn decay_integrand_dc = 4*M_PI*r*r*dc_B_chi_photon_rate_r ;  
//                           decay_integrand_dc /= (1 - 2*M_r / r).sqrt() ; 
//                           decay_integrand_dc *= exp(nu_r) ; 

//   Zaki::Vector::DataSet decay_integral_ds({r, decay_integrand_dc}) ;

//   decay_integral_ds.Interpolate(0, 1) ;    

//   // double result = decay_integral_ds.Integrate(1, {r[0], R_min}) ;
//   //       result += decay_integral_ds.Integrate(1, {R_max, r[-1]}) ;

//   // double result_full = decay_integral_ds.Integrate(1, {r[0], r[-1]}) ;
//   double result = decay_integral_ds.Integrate(1, {r[0], r[-1]}) ;

//   // double result_resonance = 0 ;
//   // if (resonance_region_exists)
//   // {
//   //   result_resonance = Resonance_B_Chi_Photon_Rate(m_chi, in_B, interp_n_resonance,
//   //                                                               interp_r_resonance) ;
//   //   // result_plus_resonance += tmp_resonance_B_Chi_Photon_Rate ;

//   //   std::cout << "\n\t Resoannce decay rate (B_dot/B) (only resonance region): " 
//   //             << result_resonance /pulsar.GetSeqPoint().b << " per year.\n" ;
//   // }

//   // Convert s^-1 into yr^-1
//   result                *= 3600*24*365 ;
//   // result_full           *= 3600*24*365 ;
//   // result_resonance      *= 3600*24*365 ;

//   // std::cout << "\n\t Full decay rate (B_dot): " << result << " per year.\n" ;
//   std::cout << "\n\t Full decay rate (B_dot/B) : " << result /pulsar.GetSeqPoint().b << " per year.\n" ;

//   std::cout << "\n ---------------------------------------- \n" ;

//   // Decay_Rate_Type decay_result( result/pulsar.GetSeqPoint().b, 
//   //                               result_full/pulsar.GetSeqPoint().b, 
//   //                               result_resonance/pulsar.GetSeqPoint().b) ;

//   return  result / pulsar.GetSeqPoint().b ;
// } 

//--------------------------------------------------------------
// Finds and plots the vacuum branching ratio limits
void BNV_Chi::PlotVacuumBrLim(const Baryon& B)
{
  // This dataset contains the upper limit on eps in MeV
  Zaki::Vector::DataSet dec_lim_ds ;
  dec_lim_ds.Reserve(2, m_chi_vals.res) ;

  for (size_t i = 0; i < m_chi_vals.res; i++)
  {
    double m_chi = m_chi_vals[i] ;
    double eps_lim = GetEpsLim(m_chi, B) ;
    dec_lim_ds.AppendRow({m_chi, Vacuum_Decay_Br(B, m_chi, eps_lim)}) ;
  }
  
  dec_lim_ds[0].label = "$m_{\\chi}\\, [ MeV ]$" ;
  dec_lim_ds[1].label = "Br$("+B.TeX_name+" \\to \\chi \\gamma)$" ;

  dec_lim_ds.SetWrkDir(wrk_dir + model + "/" 
                        + pulsar.GetName() + "/" + process.name) ;
  
  // Changing the plot parameters
  Zaki::Vector::DataSet::PlotParam plt_par ;

  plt_par.SetXAxis({0, 1.05 * B.mass }) ;
  // plt_par.SetYAxis({1e-16, 1e-14}) ;
  // plt_par.SetXTicks({{100, 200, 300, 400, 500, 600, 700, 800, 900}}) ;
  plt_par.SetYAxisLabel("${\\rm Br} \\, \\left( "+B.TeX_name+" \\to \\chi \\gamma \\, \\right)$") ;
  plt_par.SetGrid() ;
  // plt_par.SetLegend({"center", 0.5, 0.92}) ;

  dec_lim_ds.SetPlotPars(plt_par) ;

  dec_lim_ds.SemiLogYPlot(0, 1, "Br("+B.short_name+")_vs_m_chi.pdf",
                          "Inferred from $" + GetSpecificProcess(B).TeX + "$\n" 
                          + pulsar.GetName() + "$\\qquad$" + model) ;
}

//--------------------------------------------------------------
// Finds and plots the vacuum branching ratio limits
void BNV_Chi::PlotVacuumBrLim()
{
  PlotVacuumBrLim(neutron) ;
  PlotVacuumBrLim(lambda) ;
} 

//--------------------------------------------------------------
// // Imports the effective mass set
// void BNV_Chi::ImportEffMass() 
// {
//   m_B_ds.SetWrkDir(wrk_dir) ;
//   m_B_ds.Import("/EOS/CompOSE/"+ model + 
//                     "/" + model + "_m_eff.micro") ;
// }

//--------------------------------------------------------------
// // Imports the vector self-energy set
// void BNV_Chi::ImportVSelfEnergy() 
// {
//   V_self_E_ds.SetWrkDir(wrk_dir) ;
//   V_self_E_ds.Import("/EOS/CompOSE/"+ model + 
//                       "/" + model + "_V.micro") ;
// }

//--------------------------------------------------------------
// // Plots [ V_self_E_ds + m_B_ds ]
// void BNV_Chi::Plot_RestEnergy_Density()
// {
 
//   Zaki::Vector::DataSet plt_ds({m_B_ds[0], 
//       m_B_ds[1], V_self_E_ds[1], m_B_ds[2], V_self_E_ds[2], m_B_ds[3], V_self_E_ds[3]}) ;
//   plt_ds.SetWrkDir(wrk_dir) ;
//   plt_ds.Plot(0, {1,2,3,4,5,6}, "m_V.pdf", pulsar.GetName() + "\n" + model) ;
// }

//--------------------------------------------------------------
// This function plots the total rate and limit on eps 
//  as a function of m_chi per year
void BNV_Chi::PlotRate_Eps(const Baryon& B)
{
  // Columns are:
  // rate[0] = m_chi, rate[1] = Gamma, rate[2] = eps 
  Zaki::Vector::DataSet  rate ;
  rate.Reserve(3, m_chi_vals.res) ;

  for (size_t i = 0; i < m_chi_vals.res; i++)
  {
    double m_chi = m_chi_vals[i] ;

    Rate_Eps rate_eps =  GetRate_Eps(m_chi, B) ;
    rate.AppendRow( { m_chi, rate_eps.rate, rate_eps.eps }) ;

    if (i%20 == 0)
    {
      char m_chi_str[50] ;
      snprintf(m_chi_str, 50, "%.0f", m_chi) ;
      std::cout << "\n ---------------------------------------- " ;
      std::cout << "\n\t m = " << m_chi_str  
            << "\t Gamma (B_dot/B) = " 
            << rate_eps.rate << "\t per year." ;
    }
  }

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetLegend({"upper left", 0.0, 0.0}) ;
  plt_par.SetGrid() ;

  rate.SetWrkDir(wrk_dir + model + "/" + pulsar.GetName() 
                    + "/" + process.name) ;
  std::string f_name_suff = model + "_" + pulsar.GetName() 
                            + "_"+ B.short_name + ".pdf" ;

  std::string title_str = "\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + " ]"  ;
  plt_par.SetYAxisLabel("$\\Gamma ("+ GetSpecificProcess(B).TeX 
                          + ")\\, \\left({\\rm yr}^{-1}\\right)$") ;
  plt_par.SetXAxis({0, 1600}) ;
  plt_par.SetXAxisLabel("$m_{\\chi}\\, ( {\\rm MeV} ) $") ;
  plt_par.SetYAxis({1e-9, 1e-4}) ; 
  rate.SetPlotPars(plt_par) ;
  rate.SemiLogYPlot(0, {{1, B.TeX_name}}, "Rate_" + f_name_suff, title_str ) ;

  plt_par.SetLegend({"upper left", 0.0, 0.8}) ;
  plt_par.SetYAxisLabel("$\\varepsilon\\, ( {\\rm MeV} ) $") ;
  plt_par.SetXAxis({1, 1600}) ;
  plt_par.SetYAxis({1e-19, 1e-13}) ;
  rate.SetPlotPars(plt_par) ;
  rate.SemiLogYPlot(0, {{2, B.TeX_name}}, "Eps_" + f_name_suff, title_str) ;

  // ...........................................
  rate[0].label = "m_chi [MeV]" ;
  rate[1].label = "Gamma_"+B.short_name+" [1/yr]" ;
  rate[2].label = "eps_"+B.short_name+" [MeV]" ;
  
  char tmp_gamma[50];
  snprintf(tmp_gamma, 50, "%.2e", gamma_bnv_bin_lim) ;

  char tmp_eps[50];
  snprintf(tmp_eps, 50, "%.2e", default_eps) ;

  rate.AddHead("# Rates of '" + process.name + "' assuming eps = " 
               + std::string(tmp_eps) + " MeV.\n"
               "# Limits on eps are evaluated using Gamma_BNV < " 
               + std::string(tmp_gamma) 
               + ", which itself is inferred from the" 
               + " orbital decay rate of " + pulsar.GetName() + ".\n") ;

  rate.Export("Rate_Eps_"+ model + "_" + pulsar.GetName() + "_" + B.label + ".tsv") ;
}

//--------------------------------------------------------------
// Plots the total rate as a function of m_chi
void BNV_Chi::PlotRate_Eps()
{
  // Plotting each individually:
  // PlotRate_Eps(neutron) ;
  // PlotRate_Eps(lambda) ;

  // ...........................................
  // Combining plots:
  // Columns are:
  // rate[0] = m_chi, rate[1] = Gamma_n, 
  // rate[2] = Gamma_lambda, rate[3] = eps_n,  
  // rate[4] = eps_lambda
  Zaki::Vector::DataSet  rate ;
  rate.Reserve(5, m_chi_vals.res) ;

  for (size_t i = 0; i < m_chi_vals.res; i++)
  {
    double m_chi = m_chi_vals[i] ;

    Rate_Eps rate_eps_n =  GetRate_Eps(m_chi, neutron) ;
    Rate_Eps rate_eps_lam =  GetRate_Eps(m_chi, lambda) ;
    rate.AppendRow( { m_chi, rate_eps_n.rate, 
      rate_eps_lam.rate, rate_eps_n.eps, rate_eps_lam.eps }) ;
  
    if (i%20 == 0)
    {
      char m_chi_str[50] ;
      snprintf(m_chi_str, 50, "%.0f", m_chi) ;
      std::cout << "\n ---------------------------------------- " ;
      std::cout << "\n\t m = " << m_chi_str  
            << "\t Gamma_n (B_dot/B) = " 
            << rate_eps_n.rate
            << ",\t Gamma_Lam (B_dot/B) = " 
            << rate_eps_lam.rate << "\t per year." ;
    }
    
  }

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetLegend({"lower left", 0.0, 0.0}) ;
  plt_par.SetGrid() ;

  rate.SetWrkDir(wrk_dir + model + "/" + pulsar.GetName() 
                    + "/" + process.name) ;
  std::string f_name_suff = model + "_" + pulsar.GetName() 
                            + ".pdf" ;

  std::string title_str = "\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + " ]"  ;
  plt_par.SetYAxisLabel("$\\Gamma \\, \\left({\\rm yr}^{-1}\\right)$") ;
  plt_par.SetXAxis({0, 1600}) ;
  plt_par.SetXAxisLabel("$m_{\\chi}\\, ( {\\rm MeV} ) $") ;
  // plt_par.SetYAxis({4e-1, 3e7}) ; // For Combo
  plt_par.SetYAxis({1e-15, 1e-8}) ; 

  rate.SetPlotPars(plt_par) ;
  rate.SemiLogYPlot(0, {{1, "$n$"}, {2, "$\\Lambda$"}}, 
    "Rate_" + f_name_suff, title_str ) ;

  plt_par.SetLegend({"upper left", 0.0, 0.8}) ;
  plt_par.SetYAxisLabel("$\\varepsilon\\, ( {\\rm MeV} ) $") ;
  plt_par.SetXAxis({1, 1400}) ;
  // plt_par.SetYAxis({1e-27, 1e-4}) ; // For eps_max results
  plt_par.SetYAxis({1e-19, 1e-13}) ; // For eps_max results
  // plt_par.SetYAxis({2e-19, 3e-15}) ; // For Combo
  rate.SetPlotPars(plt_par) ;
  rate.SemiLogYPlot(0, {{3, "$n$"}, {4, "$\\Lambda$"}}, 
    "Eps_" + f_name_suff, title_str) ;
  // ...........................................
  rate[0].label = "m_chi [MeV]" ;
  rate[1].label = "Gamma_n [1/s]" ;
  rate[2].label = "Gamma_lam [1/s]" ;
  rate[3].label = "eps_n [MeV]" ;
  rate[4].label = "eps_lam [MeV]" ;
  
  char tmp_gamma[50];
  snprintf(tmp_gamma, 50, "%.2e", gamma_bnv_bin_lim) ;

  char tmp_eps[50];
  snprintf(tmp_eps, 50, "%.2e", default_eps) ;

  rate.AddHead("# Rates of '" + process.name + "' assuming eps = " 
               + std::string(tmp_eps) + " MeV.\n"
               "# Limits on eps are evaluated using Gamma_BNV < " 
               + std::string(tmp_gamma) 
               + ", which itself is inferred from the" 
               + " orbital decay rate of " + pulsar.GetName() + ".\n") ;

  rate.Export("Rate_Eps_"+ model + "_" + pulsar.GetName() + ".tsv") ;
}

//--------------------------------------------------------------
// // This function returns the B -> chi decay rate 
// //  as a function of radius in units of s^-1/fm^3
// Zaki::Vector::DataSet BNV_Chi::Rate_vs_R(
//                                 const double& m_chi, 
//                                 const Baryon& B, 
//                                 const bool& gen_plots)
// { 
//   double eps = 1e-10 ;

//   // Zaki::Vector::DataColumn n_tot = eos.GetEOS(2) ;
//   // Zaki::Vector::DataColumn n_B = n_tot * eos.GetEOS(B.label) ;
//   // Zaki::Vector::DataColumn m_B = eos.GetMeff(B.label) ;
//   // Zaki::Vector::DataColumn sig0 = eos.GetVeff(B.label) ;

//   // Zaki::Vector::DataColumn dc_Sigma_Plus = SigmaPlus(m_chi, m_B, n_B) ;
//   // Zaki::Vector::DataColumn dc_Sigma_Minus = SigmaMinus(m_chi, m_B, n_B) ; 

//   Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;
//   Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
//   Zaki::Vector::DataColumn n_B_r = n_r * pulsar.GetProfile()->operator[](B.label) ;

//   Zaki::Vector::DataSet ds_m_B_n({n_B["n_tot"], eos.GetMeff(B.label)}) ;
//   Zaki::Vector::DataSet sig0_n({n_B["n_tot"], eos.GetVeff(B.label) }) ;
//   ds_m_B_n.Interpolate(0, 1) ;
//   sig0_n.Interpolate(0, 1) ;
//   Zaki::Vector::DataColumn m_B_r = ds_m_B_n.Evaluate(1, n_r) ;
//   Zaki::Vector::DataColumn sig0_r = sig0_n.Evaluate(1, n_r) ;

//   Zaki::Vector::DataSet  rate ;
//   rate.Reserve(2, r.Size() ) ;
//   for (size_t i = 0; i < r.Size() ; i++)
//   {
    
//     double m_B = m_B_r[i] ;
//     double n_B = n_B_r[i] ;
//     double sig0 = sig0_r[i] ;

//     if ( abs(sig0) < SigmaMinus(m_chi, m_B, n_B) 
//           || 
//          abs(sig0) > SigmaPlus(m_chi, m_B, n_B)  )
//     {
//       rate.AppendRow({r[i], 0}) ;
//       continue ;
//     }
    
//     if ( abs(sig0) > abs(m_B - m_chi) 
//         &&
//         abs(sig0) < m_B + m_chi )
//     {
//       rate.AppendRow({r[i], 0}) ;
//       continue ;
//     }

//     double p = pow(sig0 - m_B, 2) - m_chi * m_chi ;
//           p *= pow(sig0 + m_B, 2) - m_chi * m_chi  ;
//           p  = sqrt(p) ;
//           p /= 2 * abs(sig0) ;
  
//     double amp_sqrd = abs(pow(sig0, 2) - pow(m_B, 2) + m_chi * m_chi ) ;
//     amp_sqrd *= sig0 / ( 2* abs(sig0) ) ;
//     amp_sqrd += m_chi*m_chi + m_B * m_chi ;
//     amp_sqrd *= 2 * eps * eps ;

//     double rate_val = p * amp_sqrd / (2*M_PI*abs(sig0)) ;
//     rate.AppendRow({r[i], rate_val}) ;
//   }

//   if(gen_plots)
//   {
//     Zaki::String::Directory out_dir = wrk_dir + model + "/"
//             + pulsar.GetName() + "/Rate_vs_R/" + B.short_name ;
//     out_dir.Create() ;
//     rate.SetWrkDir(out_dir) ;
//     char m_chi_ch[100] ;
//     sprintf(m_chi_ch, "Rate_vs_R_%s_%.0f.pdf", B.label.c_str(), m_chi) ;
//     rate.Plot(0, 1, std::string(m_chi_ch)) ;
//   }

//   return  rate ;
// }

//--------------------------------------------------------------
// This function plots the rate 
//  as a function of radius in units of s^-1 fm^-3
void BNV_Chi::Rate_vs_R(const double& m_chi)       
{
  Rate_vs_R(m_chi, neutron, true) ;
  Rate_vs_R(m_chi, lambda, true) ;
}

//--------------------------------------------------------------
// Plots rates for a set of chi's in the same figure 
// as a function of radius in units of s^-1 fm^-3
void BNV_Chi::Rate_vs_R(const std::vector<double>& m_chi, const Baryon& B)
{
  Zaki::Vector::DataSet rate_vs_r(pulsar.GetProfile()->operator[](0)) ;
  std::vector<std::pair<int, std::string>> labels ;
  labels.reserve(m_chi.size()) ;

  double rate_max = 0 ;
  double rate_min = 1 ;
  for (size_t i = 0; i < m_chi.size(); i++)
  {
    rate_vs_r.data_set.emplace_back(Rate_vs_R(m_chi[i], B)[1]) ;

    if ( rate_max < rate_vs_r[i+1].Max())
      rate_max = rate_vs_r[i+1].Max() ;

    // Finding the non-zero minimum:
    double tmp_min = rate_vs_r[i+1].GetSubSet([](const double& v){return v > 0 ;}).Min() ;
    if ( rate_min > tmp_min )
      rate_min = tmp_min ;

    char tmp_ch[50] ;
    snprintf(tmp_ch, 50, "%.0f", m_chi[i]) ;
    labels.emplace_back(i+1, tmp_ch) ;
    rate_vs_r[i+1].label = tmp_ch ;
  }

  Zaki::String::Directory out_dir = wrk_dir + model + "/"
          + pulsar.GetName() + "/" + process.name 
          + "/Rate_vs_R/" + B.short_name ;
  out_dir.Create() ;
  rate_vs_r.SetWrkDir(out_dir) ;
  char file_ch[200] ;
  snprintf(file_ch, 200, "Rate_vs_R_%s.pdf", B.label.c_str()) ;
  std::string title_str ="\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + "]" ;

  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxisLabel("$R\\, ({\\rm km})$") ;
  plt_pars.SetYAxisLabel("$\\Gamma ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left(s^{-1}/{\\rm fm}^3\\right)$") ;

  // plt_pars.SetXAxis({0, rate_vs_r[0].Max()}) ;
  // if (rate_min < 1e-4 * rate_max)
  // {
  //   plt_pars.SetYAxis({1e6 * rate_min, 2 * rate_max}) ;
  // }
  
  plt_pars.SetXAxis({0, 14}) ;

  if (B.label == "10")
  {
    plt_pars.SetYAxis({1e-21, 5e-16}) ;
  }

  if (B.label == "100")
  {
    plt_pars.SetYAxis({1e-21, 1.1e-18}) ;
  }

  plt_pars.SetLegend({"lower right", 1.0, 0.0}) ;
  rate_vs_r.SetPlotPars(plt_pars) ;
  rate_vs_r.SemiLogYPlot(0, labels, std::string(file_ch), title_str) ;
  
  // ------------------------------------
  //        Exporting the dataset
  // ------------------------------------
  char tmp_eps[50];
  snprintf(tmp_eps, 50, "%.2e", default_eps) ;

  rate_vs_r.AddHead("# Rates per volume for '" 
                    + GetSpecificProcess(B).name 
                    + "' assuming eps = " 
                    + std::string(tmp_eps) 
                    + " MeV as a function of radius in '" 
                    + pulsar.GetName() + "'.\n"
                    + "# The unit for rates is in "
                    + "[ s^-1 / fm^3 ].\n"
                    ) ;

  rate_vs_r.Export("Rate_vs_R_"+ std::string(B.label.c_str()) 
                   + "_" + model + "_" 
                   + pulsar.GetName() + ".tsv") ;
  // ------------------------------------
}

//--------------------------------------------------------------
// Plots rates for a set of chi's in the same figure 
// as a function of radius in units of s^-1 fm^-3
// for multiple pulsars
//  [ Incomplete ]
void BNV_Chi::Rate_vs_R(const std::vector<double>& m_chi, 
                        const Baryon& B, 
                        const std::vector<CompactStar::Pulsar_Limit>& psr_lims)
{

  // ----------------------------------------------
  // Determine the pulsar with the largest radius
  double max_radius = 0 ;
  for (size_t i = 0; i < psr_lims.size(); i++)
  {
    CompactStar::Pulsar tmp_p = psr_lims[i].p ;
    tmp_p.SetWrkDir( wrk_dir + "/" + model + "/" + tmp_p.GetName()) ;
    tmp_p.FindProfile(model) ; 
    if (max_radius < tmp_p.GetSeqPoint().r)
    {
      max_radius = tmp_p.GetSeqPoint().r ;
    }
  }
  // ----------------------------------------------

  // Zaki::Vector::DataSet rate_vs_r(pulsar.GetProfile()->operator[](0)) ;
  // std::vector<std::pair<int, std::string>> labels ;
  // labels.reserve(m_chi.size()) ;

  // double rate_max = 0 ;
  // double rate_min = 1 ;
  // for (size_t i = 0; i < m_chi.size(); i++)
  // {
  //   rate_vs_r.data_set.emplace_back(Rate_vs_R(m_chi[i], B)[1]) ;

  //   if ( rate_max < rate_vs_r[i+1].Max())
  //     rate_max = rate_vs_r[i+1].Max() ;

  //   // Finding the non-zero minimum:
  //   double tmp_min = rate_vs_r[i+1].GetSubSet([](const double& v){return v > 0 ;}).Min() ;
  //   if ( rate_min > tmp_min )
  //     rate_min = tmp_min ;

  //   char tmp_ch[50] ;
  //   snprintf(tmp_ch, 50, "%.0f", m_chi[i]) ;
  //   labels.emplace_back(i+1, tmp_ch) ;
  //   rate_vs_r[i+1].label = tmp_ch ;
  // }

  // Zaki::String::Directory out_dir = wrk_dir + model + "/"
  //         + pulsar.GetName() + "/" + process.name 
  //         + "/Rate_vs_R/" + B.short_name ;
  // out_dir.Create() ;
  // rate_vs_r.SetWrkDir(out_dir) ;
  // char file_ch[200] ;
  // snprintf(file_ch, 200, "Rate_vs_R_%s.pdf", B.label.c_str()) ;
  // std::string title_str ="\n PSR" + pulsar.GetName() + "$\\qquad$"
  //                         + "[ " + model + "]" ;

  // Zaki::Vector::DataSet::PlotParam plt_pars ;
  // plt_pars.SetGrid() ;
  // plt_pars.SetXAxisLabel("$R\\, ({\\rm km})$") ;
  // plt_pars.SetYAxisLabel("$\\Gamma ("+ GetSpecificProcess(B).TeX 
  //                       + ")\\,\\, \\left(s^{-1}/{\\rm fm}^3\\right)$") ;

  // // plt_pars.SetXAxis({0, rate_vs_r[0].Max()}) ;
  // // if (rate_min < 1e-4 * rate_max)
  // // {
  // //   plt_pars.SetYAxis({1e6 * rate_min, 2 * rate_max}) ;
  // // }
  
  // plt_pars.SetXAxis({0, 14}) ;

  // if (B.label == "10")
  // {
  //   plt_pars.SetYAxis({1e-21, 5e-16}) ;
  // }

  // if (B.label == "100")
  // {
  //   plt_pars.SetYAxis({1e-21, 1.1e-18}) ;
  // }

  // plt_pars.SetLegend({"lower right", 1.0, 0.0}) ;
  // rate_vs_r.SetPlotPars(plt_pars) ;
  // rate_vs_r.SemiLogYPlot(0, labels, std::string(file_ch), title_str) ;
  
  // // ------------------------------------
  // //        Exporting the dataset
  // // ------------------------------------
  // char tmp_eps[50];
  // snprintf(tmp_eps, 50, "%.2e", default_eps) ;

  // rate_vs_r.AddHead("# Rates per volume for '" 
  //                   + GetSpecificProcess(B).name 
  //                   + "' assuming eps = " 
  //                   + std::string(tmp_eps) 
  //                   + " MeV as a function of radius in '" 
  //                   + pulsar.GetName() + "'.\n"
  //                   + "# The unit for rates is in "
  //                   + "[ s^-1 / fm^3 ].\n"
  //                   ) ;

  // rate_vs_r.Export("Rate_vs_R_"+ std::string(B.label.c_str()) 
  //                  + "_" + model + "_" 
  //                  + pulsar.GetName() + ".tsv") ;
  // // ------------------------------------
}

//--------------------------------------------------------------
// Plots rates for a set of chi's in the same figure 
// as a function of radius in units of s^-1/fm^3
void BNV_Chi::Rate_vs_R(const std::vector<double>& m_chi) 
{
  Rate_vs_R(m_chi, neutron) ;
  Rate_vs_R(m_chi, lambda) ;
}

//--------------------------------------------------------------
// Plots rates for a set of chi's in the same figure 
// as a function of radius in units of s^-1/fm^3
// for multiple pulsars
void BNV_Chi::Rate_vs_R(const std::vector<double>& m_chi, 
            const std::vector<CompactStar::Pulsar_Limit>& psr_lims) 
{
  Rate_vs_R(m_chi, neutron, psr_lims) ;
  Rate_vs_R(m_chi, lambda, psr_lims) ;
}

//--------------------------------------------------------------
void BNV_Chi::hidden_Plot_Rate_vs_R(const Baryon& B, 
                                    const double& m_chi,
                                    Zaki::Vector::DataSet* ds) const
{
  Zaki::String::Directory out_dir = wrk_dir + model + "/"
          + pulsar.GetName() + "/" + process.name 
          + "/Rate_vs_R/" + B.short_name ;
  out_dir.Create() ;
  ds->SetWrkDir(out_dir) ;
  char file_ch[200] ;
  snprintf(file_ch, 200, "Rate_vs_R_%s_%.0f.pdf", B.label.c_str(), m_chi) ;
  std::string title_str ="\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + "]" ;

  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxisLabel("$R\\, ({\\rm km})$") ;
  plt_pars.SetYAxisLabel("$\\Gamma ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left(s^{-1}/{\\rm fm}^3\\right)$") ;

  // double rate_max = ds->operator[](1).Max() ;
  // double rate_min = ds->operator[](1).Min() ;

  // if (rate_min < 1e-4 * rate_max)
  // {
  //   plt_pars.SetYAxis({1e-4 * rate_max, 2 * rate_max}) ;
  // }
  
  plt_pars.SetXAxis({0, 14}) ;

  if (B.label == "10")
  {
    plt_pars.SetYAxis({1e-21, 5e-16}) ;
  }

  ds->SetPlotPars(plt_pars) ;
  ds->SemiLogYPlot(0, 1, std::string(file_ch), title_str) ;
}

//--------------------------------------------------------------
void BNV_Chi::hidden_Plot_Rate_vs_Density(const Baryon& B, 
                                    const double& m_chi,
                                    Zaki::Vector::DataSet* ds) const
{
  Zaki::String::Directory out_dir = wrk_dir + model 
                  + "/" + pulsar.GetName() + "/" + process.name
                  + "/Rate_vs_Den/" + B.short_name ;

  out_dir.Create() ;
  ds->SetWrkDir(out_dir) ;
  char file_ch[200] ;
  snprintf(file_ch, 200, "Rate_vs_n_%s_%.0f.pdf", B.label.c_str(), m_chi) ;
  std::string title_str = "\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + " ]"  ;

  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxisLabel("$n\\, ({\\rm fm}^{-3})$") ;
  plt_pars.SetYAxisLabel("$\\Gamma ("+ GetSpecificProcess(B).TeX 
                          + ")\\, \\left(s^{-1}/{\\rm fm}^3\\right)$") ;

  double rate_max = ds->operator[](1).Max() ;
  double rate_min = ds->operator[](1).Min() ;
  
  if (rate_min < 1e-4 * rate_max)
  {
    plt_pars.SetYAxis({1e-4 * rate_max, 2 * rate_max}) ;
  }

  ds->SetPlotPars(plt_pars) ;
  ds->SemiLogYPlot(0, 1, std::string(file_ch), title_str) ;
}

//--------------------------------------------------------------
// /// Returns the B -> chi decay rate 
// ///  as a function of density in units of s^-1/fm^3
// Zaki::Vector::DataSet BNV_Chi::Rate_vs_Density(
//                                       const double& m_chi, 
//                                       const Baryon& B,
//                                 const bool& gen_plots)
// {
//  double eps = 1e-10 ;

//   // Finding the total and individual baryon density
//   // Zaki::Vector::DataColumn n_tot = eos.GetEOS(2) ;
//   // Zaki::Vector::DataColumn n_B = n_tot * eos.GetEOS(B.label) ;
//   // Zaki::Vector::DataColumn m_B = eos.GetMeff(B.label) ;
//   // Zaki::Vector::DataColumn sig0 = eos.GetVeff(B.label) ;

//   Zaki::Vector::DataColumn dc_Sigma_Plus = SigmaPlus(m_chi, B) ;
//   Zaki::Vector::DataColumn dc_Sigma_Minus = SigmaMinus(m_chi, B) ; 

//   Zaki::Vector::DataSet  rate ;
//   rate.Reserve(2, n_B["n_tot"].Size() ) ;
//   for (size_t i = 0; i < n_B["n_tot"].Size() ; i++)
//   {
//     if ( abs(sig0_ds[B.label][i]) < dc_Sigma_Minus[i] 
//           || 
//          abs(sig0_ds[B.label][i]) > dc_Sigma_Plus[i]  )
//     {
//       rate.AppendRow({n_B["n_tot"][i], 0}) ;
//       continue ;
//     }
    
//     if ( abs(sig0_ds[B.label][i]) > abs(m_B_ds[B.label][i] - m_chi) 
//         &&
//         abs(sig0_ds[B.label][i]) < m_B_ds[B.label][i] + m_chi )
//     {
//       rate.AppendRow({n_B["n_tot"][i], 0}) ;
//       continue ;
//     }

//     double p = pow(sig0_ds[B.label][i] - m_B_ds[B.label][i], 2) - m_chi * m_chi ;
//           p *= pow(sig0_ds[B.label][i] + m_B_ds[B.label][i], 2) - m_chi * m_chi  ;
//           p  = sqrt(p) ;
//           p /= 2 * abs(sig0_ds[B.label][i]) ;
  
//     double amp_sqrd = abs(pow(sig0_ds[B.label][i], 2) - pow(m_B_ds[B.label][i], 2) + m_chi * m_chi ) ;
//     amp_sqrd *= sig0_ds[B.label][i] / ( 2* abs(sig0_ds[B.label][i]) ) ;
//     amp_sqrd += m_chi*m_chi + m_B_ds[B.label][i] * m_chi ;
//     amp_sqrd *= 2 * eps * eps ;

//     double rate_val = p * amp_sqrd / (2*M_PI*abs(sig0_ds[B.label][i])) ;
//     rate.AppendRow({n_B["n_tot"][i], rate_val}) ;
//   }

//   if(gen_plots)
//   {
//     Zaki::String::Directory out_dir = wrk_dir + model 
//                     + "/Rate_vs_Den/" + B.short_name ;

//     out_dir.Create() ;
//     rate.SetWrkDir(out_dir) ;
//     char m_chi_ch[100] ;
//     sprintf(m_chi_ch, "Rate_vs_n_%s_%.0f.pdf", B.label.c_str(), m_chi) ;
//     rate.Plot(0, 1, std::string(m_chi_ch)) ;
//   }

//   return  rate ;
// }

// //--------------------------------------------------------------
// // This function returns the upper limit on the vector 
// // self-energy, from condition (b). 
// // The arguments are (m_chi, m_B^*, n [fm^-3] ).
// Zaki::Vector::DataColumn BNV_Chi::SigmaPlus(
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

// //--------------------------------------------------------------
// // This function returns the lower limit on the vector 
// // self-energy, from condition (a).
// // The arguments are (m_chi, m_B^*, n [fm^-3] ).
// Zaki::Vector::DataColumn BNV_Chi::SigmaMinus(
//                         const double& in_m_chi, 
//                         const Baryon& B)
// { 

//   Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI*n_B[B.label]).pow(2.0/3.0) ;

//   kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

//   Zaki::Vector::DataColumn  out  = m_B_ds[B.label]*m_B_ds[B.label] 
//                                   + in_m_chi*in_m_chi + 2*kF_2 ;
  
//   out -= 2*((m_B_ds[B.label]*m_B_ds[B.label] + kF_2)*(kF_2 + in_m_chi*in_m_chi )).sqrt() ;
//   out = out.sqrt() ;
//   out.label = "$\\Sigma^-_0 (MeV)$" ;
//   return  out ;
// }

// //--------------------------------------------------------------
// // This function returns the lower limit on the vector 
// // self-energy, from condition (a).
// // The arguments are (m_chi, m_B^*, n [fm^-3] ).
// double BNV_Chi::SigmaMinus(
//                         const double& m_chi, 
//                         const double& m_B, 
//                         const double& n_B)
// { 

//   double kF_2 = pow(3*M_PI*M_PI*n_B, 2.0/3.0) ;

//   kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

//   double  out  = m_B*m_B + m_chi*m_chi + 2*kF_2 ;
  
//   out -= 2*sqrt((m_B*m_B + kF_2)*(kF_2 + m_chi*m_chi )) ;
//   out = sqrt(out) ;

//   return  out ;
// }

// //--------------------------------------------------------------
// // This function returns the upper limit on the vector 
// // self-energy, from condition (b). 
// // The arguments are (m_chi, m_B^*, n [fm^-3] ).
// double BNV_Chi::SigmaPlus(
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
// This function returns the (per baryon) rate for
//  B -> chi decay integrated over the radius 
//  and in units of [ yr^-1 ]
double BNV_Chi::GetRate(const double& m_chi, 
                          const Baryon& B)
{
  Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
  Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
  Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;
 
  Zaki::Vector::DataColumn dc_rate_r = Rate_vs_R(m_chi, B)[1] ;

  // Converting "MeV^3" into "1/fm^3"
  // dc_rate_r *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "1/fm^3" into "1/km^3"
  dc_rate_r *= 1e+54 ; 

  // Converting "MeV" into "s^-1"
  // dc_rate_r *= 1e-3 * Zaki::Physics::GEV_2_S ; 

  Zaki::Vector::DataColumn decay_integrand_dc = 4*M_PI*r*r*dc_rate_r ;  
                          decay_integrand_dc /= (1 - 2*M_r / r).sqrt() ; 
                          decay_integrand_dc *= exp(nu_r) ; 
    Zaki::Vector::DataSet decay_integral_ds({r, decay_integrand_dc}) ;

  decay_integral_ds.Interpolate(0, 1) ;    

  double result = decay_integral_ds.Integrate(1, {r[0], r[-1]}) ;
  
  // Convert s^-1 into yr^-1
  result                *= 3600*24*365 ;

  std::cout << "\n ---------------------------------------- " ;
  std::cout << "\n\t Full decay rate (B_dot/B) : " 
            << result / pulsar.GetSeqPoint().b << " per year." ;
  // std::cout << "\n ---------------------------------------- " ;

  return  result / pulsar.GetSeqPoint().b ;
}

//--------------------------------------------------------------
// This function returns the limit on eps from 
//  B -> chi decay rate integrated over the radius 
//  and in units of [ MeV ]
double BNV_Chi::GetEpsLim(const double& m_chi, 
                          const Baryon& B)
{
  // Zaki::Vector::DataColumn n_r = puls.GetProfile()->operator[](5) ;
  Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
  Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
  Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

  // Zaki::Vector::DataSet ds_B_chi_rate_n = 
  //   B_to_Chi_Rate_vs_Density(m_chi, B, eos, dir) ;
  // ds_B_chi_rate_n.Interpolate(0, 1) ;
  
  // Zaki::Vector::DataColumn dc_B_chi_rate_r = 
  //   ds_B_chi_rate_n.Evaluate(1, n_r) ;
 
 Zaki::Vector::DataColumn dc_rate_r = Rate_vs_R(m_chi, B)[1] ;

  // Converting "MeV^3" into "1/fm^3"
  // dc_rate_r *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "1/fm^3" into "1/km^3"
  dc_rate_r *= 1e+54 ; 

  // Converting "MeV" into "s^-1"
  // dc_rate_r *= 1e-3 * Zaki::Physics::GEV_2_S ; 

  Zaki::Vector::DataColumn decay_integrand_dc = 4*M_PI*r*r*dc_rate_r ;  
                          decay_integrand_dc /= (1 - 2*M_r / r).sqrt() ; 
                          decay_integrand_dc *= exp(nu_r) ; 
    Zaki::Vector::DataSet decay_integral_ds({r, decay_integrand_dc}) ;

  decay_integral_ds.Interpolate(0, 1) ;    

  double result = decay_integral_ds.Integrate(1, {r[0], r[-1]}) ;
  
  // Convert s^-1 into yr^-1
  result *= 3600*24*365 ;

  std::cout << "\n ---------------------------------------- " ;
  std::cout << "\n\t Full decay rate (B_dot/B) : " 
            << result / pulsar.GetSeqPoint().b << " per year." ;
  // std::cout << "\n ---------------------------------------- " ;

  return default_eps*sqrt(gamma_bnv_bin_lim / (result / pulsar.GetSeqPoint().b) ) ;
}

//--------------------------------------------------------------
// This function returns the rate and limit on eps from 
//  B -> chi decay rate integrated over the radius 
//  and in units of [ yr^-1 ]
BNV_Chi::Rate_Eps BNV_Chi::GetRate_Eps(
                          const double& m_chi, 
                          const Baryon& B
                          // const CompactStar::CompOSE_EOS& eos, 
                          // const Zaki::String::Directory& dir, 
                          // CompactStar::Pulsar& puls
                          )
{
  Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
  Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
  Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

  Zaki::Vector::DataColumn dc_rate_r = Rate_vs_R(m_chi, B)[1] ;

  // // Converting "MeV^3" into "1/fm^3"
  // dc_rate_r *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "1/fm^3" into "1/km^3"
  dc_rate_r *= 1e+54 ; 

  // Converting "MeV" into "s^-1"
  // dc_rate_r *= 1e-3 * Zaki::Physics::GEV_2_S ; 

  Zaki::Vector::DataColumn decay_integrand_dc = 4*M_PI*r*r*dc_rate_r ;  
                          decay_integrand_dc /= (1 - 2*M_r / r).sqrt() ; 
                          decay_integrand_dc *= exp(nu_r) ; 
    Zaki::Vector::DataSet decay_integral_ds({r, decay_integrand_dc}) ;

  decay_integral_ds.Interpolate(0, 1) ;    

  double result = decay_integral_ds.Integrate(1, {r[0], r[-1]}) ;
  
  // Convert s^-1 into yr^-1
  result                *= 3600*24*365 ;

  // char m_chi_str[50] ;
  // snprintf(m_chi_str, 50, "%.0f", m_chi) ;
  // std::cout << "\n ---------------------------------------- " ;
  // std::cout << "\n\t m = " << m_chi_str  
  //           << "\t Full decay rate (B_dot/B) : " 
  //           << result / pulsar.GetSeqPoint().b << " per year." ;
  
  double output_rate = result / pulsar.GetSeqPoint().b ;
  double output_eps = 0 ;

  if (output_rate == 0)
  {
    output_eps = std::numeric_limits<double>::max() ;
  }
  else
  {
    output_eps = default_eps*sqrt(gamma_bnv_bin_lim / output_rate ) ;
  }
  
  return {output_rate, output_eps} ;
}

//--------------------------------------------------------------
// Saves the results
void BNV_Chi::Export(const Zaki::String::Directory& in_dir) const
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
//==============================================================
