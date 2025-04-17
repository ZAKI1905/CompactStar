/*
  BNV_B_Chi_Photon class
*/

#include <gsl/gsl_integration.h>

#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Math/Newton.hpp>

#include <Zaki/Vector/DataSet.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/BNV/BNV_B_Chi_Photon.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "CompactStar/Core/Pulsar.hpp"

using namespace CompactStar ;

//==============================================================
// bool MassCondition(const CompactStar::NStar& in_star)
// {
//   // return ( 2.0099 <= in_star.GetSequence().m 
//   //           && 
//   //         in_star.GetSequence().m <= 2.01001 )  ;

//   return ( 2.0 <= in_star.GetSequence().m 
//             && 
//           in_star.GetSequence().m <= 2.011 )  ;
// }
//==============================================================
double x_5_coeff = 0 ;
double x_7_coeff = 0 ;
double x_242_coeff = 0 ;

//==============================================================
double Eps_From_Beta_Imbalance_Eq(const double& x)
{
  return x_5_coeff * pow(x, 5) + x_7_coeff * pow(x, 7) - x_242_coeff * pow(x, 2.42) ;
}

//==============================================================
double Eps_From_Beta_Imbalance_Eq_Der(const double& x)
{
  return 5 * x_5_coeff * pow(x, 4) + 7 * x_7_coeff * pow(x, 6) 
          - 2.42 * x_242_coeff * pow(x, 1.42) ;
}

//==============================================================
//==============================================================
//                        BNV_B_Chi_Photon class
//==============================================================
// Constructor
BNV_B_Chi_Photon::BNV_B_Chi_Photon() 
  : 
  BNV_Chi({"B_Chi_Photon", "B \\to \\chi + \\gamma"})
  // neutron("neutron", "n", "n", "10", Zaki::Physics::NEUTRON_M_FM,
  //          -3.82608545, 879.6),
  // lambda("lambda", "\\Lambda", "Lam", "100", 
  //         Zaki::Physics::LAMBDA_ZERO_M_FM, -1.22, 2.632e-10), 
  // proton("proton", "p^+", "p", "11", Zaki::Physics::PROTON_M_FM,
  //          5.5856946893, 1e+41)
{ }

//--------------------------------------------------------------
// Destructor
BNV_B_Chi_Photon::~BNV_B_Chi_Photon() { }
//--------------------------------------------------------------
BNV_Chi::Process BNV_B_Chi_Photon::GetSpecificProcess(const Baryon& B) const
{
  return {B.short_name + "_Chi_Photon", B.TeX_name + " \\to \\chi + \\gamma"} ;
}


//--------------------------------------------------------------
// Returns the cosine of the angle between p_chi^{CM} and p_B^{nm}
// Added on Apr 23, 2024
double BNV_B_Chi_Photon::Cos_Theta_CMpChi_NMpB(const double& E_chi_NM_over_mB, const double& mu_chi, 
                                               const double& x, const double& sigma_0) 
{
  // If p_B = 0:
  if (x == 1)
    return +2.123456789 ; // Not valid

  // double c_theta = 2 * E_chi_NM ;
  //       c_theta -= E_B_NM * (1 + m_chi * m_chi / (E_B_NM*E_B_NM - p_B_NM * p_B_NM) ) ;
  //       c_theta /= p_B_NM ;
  //       c_theta /= (1 - m_chi * m_chi / (E_B_NM*E_B_NM - p_B_NM * p_B_NM) ) ;

  double eB_NM2_minus_pBNM2 = 1 + sigma_0*sigma_0 + 2*x*sigma_0 ;

  // Apr 29, 2024: Made the definiton dimensionless
  double c_theta = 2 * E_chi_NM_over_mB ;
        c_theta -= (x + sigma_0) * (1 + mu_chi * mu_chi / eB_NM2_minus_pBNM2 ) ;
        c_theta /= sqrt(x*x - 1) ;
        c_theta /= (1 - mu_chi * mu_chi / eB_NM2_minus_pBNM2 ) ;

  return c_theta ;
}

//--------------------------------------------------------------
// Returns the cosine of the angle between p_chi^{CM} and p_B^{nm}
// Assuming that the escape condition is satisfied
// This returns the minimum of cos(theta) for such condition.
// If it is larger than 1, then it doesn't exist.
// Added on Apr 23, 2024
// To do: Currently idx = 0 (core of the neutron star at R=0)
//      Either input radius or density, not sure yet! 
// double BNV_B_Chi_Photon::Cos_Theta_CMpChi_NMpB_Escape(const double& mu_chi, 
//                                                const double& x, const double& sigma_0) 
// {
//   int idx = 0 ;
//   double nu = pulsar.GetProfile()->operator[](6)[idx] ;
//   double exp_2nu = exp(-2*nu) ;
//   double c_theta = Cos_Theta_CMpChi_NMpB(mu_chi * exp_2nu, mu_chi,  x, sigma_0) ;

//   if ( abs(c_theta) >= 1 )
//     return 1 ;

//   return c_theta ;
// }

//--------------------------------------------------------------
// Phase space integrand function (dim-less part) for 
//    escaping chi particles
// The total phase space integrand is the output of this function
// multiplied by:
//                  m^4 C / (pi^2)
//
double BNV_B_Chi_Photon::Escape_PhaseSpace_Integrand(const double& x) 
{
  double mu_chi = phase_int_mu_chi ; 
  double sigma_0 = phase_int_sigma_0 ;
  double exp_2nu = phase_int_escape_exp_2nu ;

  if ( x <= ( mu_chi*mu_chi - 1 - sigma_0*sigma_0 ) / ( 2*sigma_0 )  ) 
  {
    return 0 ;
  }

  if (x == 0 )
  {
    double zero_mom_rate = - sigma_0 * mu_chi * mu_chi / 4.0 ;
          zero_mom_rate *= pow(2*mu_chi * exp_2nu*mu_chi*(1+sigma_0) - (1+sigma_0)*(1+sigma_0) - mu_chi * mu_chi, 2);
          zero_mom_rate /= pow(1+sigma_0, 2) ;
          zero_mom_rate /= pow(1+sigma_0, 2) - mu_chi * mu_chi ;

    return zero_mom_rate ;
  }

  double cos_theta_escape = Cos_Theta_CMpChi_NMpB(mu_chi * exp_2nu, mu_chi, x, sigma_0) ;

  if (cos_theta_escape < -1)
  {
    cos_theta_escape = -1 ;
  }

  // std::cout << "mu_chi= " << mu_chi << ", exp_2nu= " << exp_2nu << "\n" ;

  if ( cos_theta_escape > 1 )
  {
    return 0 ;
  }

  double sin_theta_escape_sqrd = 1.0 - cos_theta_escape*cos_theta_escape ;

  // Undoing the integral over cos(theta) which led to "*2"
  double polar = 0.25 * sigma_0 * sin_theta_escape_sqrd ;

        polar *= x*x - 1  ;
        polar *= 1 + sigma_0*sigma_0 + 2*x*sigma_0 - mu_chi*mu_chi ;
  
        polar /= pow(1 + sigma_0*sigma_0 + 2*x*sigma_0, 2) ;

        polar *= pow(mu_chi, 2) - 1 - sigma_0*sigma_0 - 2*x*sigma_0 ;

  // Undoing the integral over cos(theta) which led to "*2"
  double isotropic = (1.0 - cos_theta_escape) / 2.0 ;

        isotropic *= sqrt(x*x - 1) ;
        isotropic *= 1 + sigma_0*sigma_0 + 2*x*sigma_0 - mu_chi*mu_chi ;
  
        isotropic /= pow(1 + sigma_0*sigma_0 + 2*x*sigma_0, 2) ;

        isotropic *= (1 + sigma_0*sigma_0 + 2*x*sigma_0) 
                * (1 + x*sigma_0 + 2 * mu_chi) 
                + pow(mu_chi, 2) * (1 + x*sigma_0) ;

  return isotropic + polar ;
}

//--------------------------------------------------------------
// The full phase-space integral for escaping particles
// The output is in MeV^4
double BNV_B_Chi_Photon::Escape_PhaseSpace_Integral(const Baryon& B,
                                              const double& m, 
                                              const double& m_chi, 
                                              const double& Sigma_0, 
                                              const double& x_F) 
{ 
  double q_e = Zaki::Physics::Q_E ;
  
  phase_int_mu_chi = m_chi / m ; 
  phase_int_sigma_0 = Sigma_0 / m ;

  // int idx = 0 ;
  // double nu = pulsar.GetProfile()->operator[](6)[idx] ;
  // phase_int_escape_exp_2nu = exp(-2*nu) ;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
  double err, result;

  Zaki::Math::GSLFuncWrapper<BNV_B_Chi_Photon, double (BNV_B_Chi_Photon::*)(const double&)> 
    Fp(this, &BNV_B_Chi_Photon::Escape_PhaseSpace_Integrand);     

  gsl_function F = static_cast<gsl_function> (Fp) ; 

  gsl_integration_qag(&F, 1, x_F, 1e-10, 1e-10, 2000, 1, w, &result, &err) ;
  gsl_integration_workspace_free(w);

  // C
  double C = pow(B.g * default_eps * q_e, 2) / (128 * M_PI * m * m ) ;

  result *= C ;
  result *= pow(m, 4) / (M_PI*M_PI)  ;

  return result ;
}

//--------------------------------------------------------------
// Evaluates the escape rate as a function of radius                                          
Zaki::Vector::DataSet BNV_B_Chi_Photon::Escape_Rate_vs_R(const Baryon& B, 
                                          const double& m_chi,
                                          const bool& export_flag) 
{
  Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI
                        * micro_r["n_" + B.label]).pow(2.0/3.0) ;
  kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

  Zaki::Vector::DataColumn m_B = micro_r["m_" + B.label] ;
  Zaki::Vector::DataColumn Sigma_0 = micro_r["sig_" + B.label]  ;
  Zaki::Vector::DataColumn x_F = (m_B*m_B + kF_2).sqrt() / m_B ;

  Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

  Zaki::Vector::DataSet rate_vs_r ;
  rate_vs_r.Reserve(4, m_B.Size() ) ;

  for (size_t i = 0; i < m_B.Size(); i++)
  {
    phase_int_escape_exp_2nu = exp(-2*nu_r[i]) ;
    
    double escaped_rate = Escape_PhaseSpace_Integral(B, m_B[i], m_chi, Sigma_0[i], x_F[i]) ;
    double total_rate = PhaseSpace_Integral(B, m_B[i], m_chi, Sigma_0[i], x_F[i]) ;

    rate_vs_r.AppendRow({micro_r[0][i], 
      escaped_rate,
      total_rate, 
      escaped_rate / total_rate}) ;
  }
  
  // Converting "MeV^4" into "MeV/fm^3"
  rate_vs_r[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 
  rate_vs_r[2] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "MeV/fm^3" into "s^-1/fm^3"
  rate_vs_r[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;
  rate_vs_r[2] *= 1e-3 * Zaki::Physics::GEV_2_S ;

  // ...........................................
  if (export_flag)
  {
    rate_vs_r.SetWrkDir(wrk_dir + model + "/" + pulsar.GetName() 
                                          + "/" + process.name) ;

    rate_vs_r[0].label = "r [km]" ;
    rate_vs_r[1].label = "Gamma_esc_"+B.short_name+" [s^-1/fm^3]" ;
    rate_vs_r[2].label = "Gamma_tot_"+B.short_name+" [s^-1/fm^3]" ;
    rate_vs_r[3].label = "Escape_Frac_"+B.short_name ;

    char tmp_mchi[50];
    snprintf(tmp_mchi, 50, "%.2e", m_chi) ;

    rate_vs_r.AddHead("# Rates of chi-escape due to' " + process.name + 
                  " from " + pulsar.GetName() 
                + ", and assuming m_chi = " 
                + std::string(tmp_mchi) + " MeV.\n") ;

    rate_vs_r.Export("Escape_Rate_vs_R_"+ model + "_" + pulsar.GetName() + "_" + B.label 
                      + "_m=" + tmp_mchi + ".tsv") ;
  }
  // ...........................................

  return rate_vs_r ;
}

//--------------------------------------------------------------
// Evaluates the escape ratio integrated over radius                                          
double BNV_B_Chi_Photon::Escape_Ratio_Total(const Baryon& B, 
                                          const double& m_chi,
                                          const bool& export_flag) 
{
  Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
  Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
  Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

  Zaki::Vector::DataSet ds_rate_r = Escape_Rate_vs_R(B, m_chi, export_flag) ;

  // Converting "1/fm^3" into "1/km^3"
  ds_rate_r[1] *= 1e+54 ; 
  ds_rate_r[2] *= 1e+54 ; 

  ds_rate_r[1] *= 4*M_PI*r*r ;
  ds_rate_r[1] /= (1 - 2*M_r / r).sqrt() ;
  ds_rate_r[1] *= exp(nu_r) ; 

  ds_rate_r[2] *= 4*M_PI*r*r ;
  ds_rate_r[2] /= (1 - 2*M_r / r).sqrt() ;
  ds_rate_r[2] *= exp(nu_r) ; 

  ds_rate_r.Interpolate(0, {1, 2}) ;    

  double total_escape_rate = ds_rate_r.Integrate(1, {r[0], r[-1]}) ;
  double total_rate = ds_rate_r.Integrate(2, {r[0], r[-1]}) ;
  
  // Convert s^-1 into yr^-1
  // total_escape_rate  *= 3600*24*365 ;
  // total_rate  *= 3600*24*365 ;
  
  // std::cout << "\n\t Total Rate = " << total_rate << "\n" ;
  // std::cout << "\t Escape Rate = " << total_escape_rate << "\n" ;

  return total_escape_rate / total_rate ;
}

//--------------------------------------------------------------
// Evaluates the escape ratio integrated over radius   
// for a range of m_chi values                                       
Zaki::Vector::DataSet BNV_B_Chi_Photon::Escape_Ratio_Total(const Baryon& B, 
                                          const Zaki::Math::Axis& m_chi_range,
                                          const bool& export_flag) 
{
  Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
  Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
  Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

  Zaki::Vector::DataSet ds_escape_ratio_vs_mchi(4, m_chi_range.res) ;

  Zaki::Vector::DataSet ds_rate_r(3, r.Size()) ;

  ds_rate_r[0] = r ;

  for (size_t i = 0; i < m_chi_range.res + 1; i++)
  {
    // std::cout << "\t m_chi = " << m_chi_range[i] << "\r" ;
    ds_rate_r[1] = 4*M_PI*r*r ;
    ds_rate_r[1] /= (1 - 2*M_r / r).sqrt() ;
    ds_rate_r[1] *= exp(nu_r) ; 

    ds_rate_r[2] = ds_rate_r[1] ;

    Zaki::Vector::DataSet ds_rate_r_factor = Escape_Rate_vs_R(B, m_chi_range[i], false) ;

    ds_rate_r[1] *= ds_rate_r_factor[1] ;
    ds_rate_r[2] *= ds_rate_r_factor[2] ;

    // Converting "1/fm^3" into "1/km^3"  
    // This is not needed since it cancels in a ratio.
    // ds_rate_r[1] *= 1e+54 ; 
    // ds_rate_r[2] *= 1e+54 ; 

    // std::cout << " m_chi = " << m_chi_range[i] << " MeV.\n" ;
    // std::cout << " ds_rate_r[1][0] = " << ds_rate_r[1][0] << ".\n" ;
    // std::cout << " ds_rate_r[2][0] = " << ds_rate_r[2][0] << ".\n\n" ;

    ds_rate_r.Interpolate(0, {1, 2}) ;    

    double total_escape_rate = ds_rate_r.Integrate(1, {r[0], r[-1]}) ;
    double total_rate = ds_rate_r.Integrate(2, {r[0], r[-1]}) ;
    
    // Convert s^-1 into yr^-1
    // This is not needed since it cancels in a ratio.
    // total_escape_rate  *= 3600*24*365 ;
    // total_rate  *= 3600*24*365 ;

    ds_escape_ratio_vs_mchi.AppendRow({m_chi_range[i], total_escape_rate, total_rate, total_escape_rate / total_rate}) ;
  }

  // ...........................................
  if (export_flag)
  {
    ds_escape_ratio_vs_mchi.SetWrkDir(wrk_dir + model + "/" + pulsar.GetName() 
                                          + "/" + process.name) ;

    ds_escape_ratio_vs_mchi[0].label = "m_chi [MeV]" ;
    ds_escape_ratio_vs_mchi[1].label = "Escape_Rate_" + B.short_name ;
    ds_escape_ratio_vs_mchi[2].label = "Total_Rate_" + B.short_name ;
    ds_escape_ratio_vs_mchi[3].label = "Escape_Frac_" + B.short_name ;

    char tmp_mchi_min[50];
    snprintf(tmp_mchi_min, 50, "%.2e", m_chi_range.Min()) ;

    char tmp_mchi_max[50];
    snprintf(tmp_mchi_max, 50, "%.2e", m_chi_range.Max()) ;

    ds_escape_ratio_vs_mchi.AddHead("# Rates of chi-escape due to '" + process.name + 
                  " from " + pulsar.GetName() 
                + ", as a function of m_chi in the range [ " 
                + std::string(tmp_mchi_min) + ", " + std::string(tmp_mchi_max) + "] MeV.\n") ;

    ds_escape_ratio_vs_mchi.Export("Escape_Ratio_vs_mchi_"+ model + "_" + pulsar.GetName() 
                                    + "_" + B.label + ".tsv") ;
  }
  // ...........................................

  return ds_escape_ratio_vs_mchi ;
}


//--------------------------------------------------------------
/// Sets the units to cgs
void BNV_B_Chi_Photon::UseCGSUnits() 
{
  cgs_units = true ;
}
//--------------------------------------------------------------
  
// Sets the EoS model
// void BNV_B_Chi_Photon::SetModel(const std::string& in_eos_model) 
// {
//   model = in_eos_model ;
// }

//--------------------------------------------------------------
// void BNV_B_Chi_Photon::GenSequence() const
// {
//   CompactStar::TOVSolver solver ;

//   solver.SetWrkDir( wrk_dir ) ;
//   solver.ImportEOS("EOS/CompOSE/"+ model +"/"+ model +".eos") ;

//   double e_c_min = eos.GetEOS(0).Min()*1.01 ;
//   double e_c_max = eos.GetEOS(0).Max()*0.99 ;

//   double m_best = 0 ;
//   double m_goal = 2.01 ;
//   double precision = 0.0005 ;
  
//   solver.Solve( {{e_c_min, e_c_max}, 50, "Log"}, 
//                 "BNV_2022/results/B_Chi_Photon_Decay/"+model,  
//                 "B_Chi_Photon_Decay") ;

//   Zaki::Vector::DataSet tmp_seq(wrk_dir +
//                   "/BNV_2022/results/B_Chi_Photon_Decay/"+ model,
//                   "B_Chi_Photon_Decay_Sequence.tsv") ;
  
//   e_c_min = tmp_seq[0][tmp_seq[1].MinIdx()] ;
//   e_c_max = tmp_seq[0][tmp_seq[1].MaxIdx()] ;

//   solver.ClearSequence() ;

//   while ( abs(m_best - m_goal) > precision ) {

//   // std::cout << "\n\n\t e_c_min =  " << e_c_min 
//   //         << ", \t e_c_max = " << e_c_max 
//   //         << "\n\n" << std::flush ;

//   solver.Solve( {{e_c_min, e_c_max}, 10, "Log"}, 
//                 "BNV_2022/results/B_Chi_Photon_Decay/"+model,  
//                 "B_Chi_Photon_Decay") ;
  
//   tmp_seq.Import("B_Chi_Photon_Decay_Sequence.tsv") ;


//   int m_idx = tmp_seq[1].GetClosestIdx(m_goal) ;
//   m_best = tmp_seq[1][m_idx] ;
//   e_c_min = tmp_seq[0][m_idx-1] ;
//   e_c_max = tmp_seq[0][m_idx+1] ;

//   solver.ClearSequence() ;
//   }

//   solver.AddNCondition(MassCondition) ;

//   solver.Solve( {{e_c_min, e_c_max}, 5, "Linear"}, 
//                 "BNV_2022/results/B_Chi_Photon_Decay/"+model,  
//                 "B_Chi_Photon_Decay") ;
// }

//--------------------------------------------------------------
/// The phase-space integral to find thermal energy 
/// dumped into the NS
// [Apr 15, 2023]: The old formula which was incorrect is now 
//                  replaced with the new one. 
//    Also mu_chi --> mu.
double BNV_B_Chi_Photon::Thermal_Photon_PhaseSpace_Integrand(const double& x) 
{
  double mu = thermal_photon_phase_int_mu_chi ; 
  double sig = thermal_photon_phase_int_sigma_0 ;

  if ( x <= ( mu*mu - 1 - sig*sig ) / ( 2*sig )  ) 
  {
    return 0 ;
  }

  //  -----------------------------------
  //              Wrong Answer:
  //  -----------------------------------
  // double out = sqrt(x*x - 1) ;

  //     out *= x + sigma_0 ;

  //     out *= pow(1 - mu_chi*mu_chi + sigma_0*sigma_0 + 2*x*sigma_0, 2) ;

  //     out /= pow(1 + sigma_0*sigma_0 + 2*x*sigma_0, 3) ;

  //     out *= pow(mu_chi+1, 2) 
  //           + pow(sigma_0, 2) * (2*mu_chi + 2*x*x + 1 )
  //           + (mu_chi + 1) * (mu_chi + 3) * sigma_0 * x 
  //           + pow(sigma_0, 3) * x ;
  //  -----------------------------------

  double out = 2 * ( mu + 1 ) * sig * ( 2 * mu + (mu+5)*x*x + 1 ) ;
        out += pow(sig, 2) * x * ( 3 * mu * (mu + 6) + 8*x*x + 10) ;
        out += 2 * pow(sig, 3) * ( 3 * mu + 5 * x * x + 1) ;
        out += 3 * pow(mu + 1, 2) * x +  3 * pow(sig, 4) * x ;

        out *= pow(sig*sig + 2*sig*x - mu*mu + 1, 2) ;
        out /= pow(sig*sig + 2*sig*x + 1, 3) ;
        out *= sqrt(x*x - 1) ;
  
  return out ;
}

//--------------------------------------------------------------
// Added on [Apr 15, 2023] 
// The phase-space integrand to find thermal energy 
// dumped into the NS due to Fermi-hole heating
double BNV_B_Chi_Photon::Thermal_Hole_PhaseSpace_Integrand(const double& x) 
{
  double mu = thermal_hole_phase_int_mu_chi ; 
  double sig = thermal_hole_phase_int_sigma_0 ;
  double x_f = thermal_hole_phase_int_x_F ;

  if ( x <= ( mu*mu - 1 - sig*sig ) / ( 2*sig )  ) 
  {
    return 0 ;
  }

  double out = pow(mu + 1, 2) +  pow(sig,2) * ( 2 * mu + 2*x*x + 1 ) ;
        out += (mu+1) * (mu+3) *  sig * x ;
        out += pow(sig, 3) * x ;
        
        out *= sig * (sig + 2 * x ) - mu*mu + 1 ;

        out *= x_f - x ;
        out /= x * (sig * x + 1) * (sig*sig + 2*sig*x + 1) ;
        out *= x*sqrt(x*x - 1) ;
  
  return out ;
}

//--------------------------------------------------------------
// The phase-space integral to find thermal energy 
// dumped into the NS
// [Apr 15, 2023]: The old formula which was incorrect is now 
//                  replaced with the new one. 
double BNV_B_Chi_Photon::Thermal_Photon_PhaseSpace_Integral(
                        const Baryon& B,
                        const double& m, 
                        const double& m_chi, 
                        const double& Sigma_0, 
                        const double& x_F) 
{
  double q_e = Zaki::Physics::Q_E ;
  
  thermal_photon_phase_int_mu_chi = m_chi / m ; 
  thermal_photon_phase_int_sigma_0 = Sigma_0 / m ;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
  double err, result;

  Zaki::Math::GSLFuncWrapper<BNV_B_Chi_Photon, double (BNV_B_Chi_Photon::*)(const double&)> 
    Fp(this, &BNV_B_Chi_Photon::Thermal_Photon_PhaseSpace_Integrand);     

  gsl_function F = static_cast<gsl_function> (Fp) ; 

  gsl_integration_qag(&F, 1, x_F, 1e-10, 1e-10, 2000, 1, w, &result, &err) ;
  gsl_integration_workspace_free(w);

  // Factor behind the integral
  double C = pow(B.g * default_eps * q_e, 2) * pow(m, 3) 
              / (768 * M_PI * M_PI * M_PI ) ;

  result *= C ;

  return result ;
}

//--------------------------------------------------------------
/// The phase-space integral to find thermal energy 
/// dumped into the NS to Fermi-hole heating
double BNV_B_Chi_Photon::Thermal_Hole_PhaseSpace_Integral( const Baryon& B,
                            const double& m, 
                            const double& m_chi, 
                            const double& Sigma_0, 
                            const double& x_F) 
{
  double q_e = Zaki::Physics::Q_E ;
  
  thermal_hole_phase_int_mu_chi = m_chi / m ; 
  thermal_hole_phase_int_sigma_0 = Sigma_0 / m ;
  thermal_hole_phase_int_x_F = x_F ;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
  double err, result;

  Zaki::Math::GSLFuncWrapper<BNV_B_Chi_Photon, double (BNV_B_Chi_Photon::*)(const double&)> 
    Fp(this, &BNV_B_Chi_Photon::Thermal_Hole_PhaseSpace_Integrand);     

  gsl_function F = static_cast<gsl_function> (Fp) ; 

  gsl_integration_qag(&F, 1, x_F, 1e-10, 1e-10, 2000, 1, w, &result, &err) ;
  gsl_integration_workspace_free(w);

  // Factor behind the integral
  double C = pow(B.g * default_eps * q_e, 2) * pow(m, 3) 
              / (128 * M_PI * M_PI * M_PI ) ;

  result *= C ;

  return result ;
}

//--------------------------------------------------------------
// Returns the rate of thermal energy dumped into the NS
//  due to Fermi-hole heating
//  per unit volume as a function of density in units 
//  of [ MeV fm^-3 s^-1 ], or [ erg cm^-3 s^-1 ].
Zaki::Vector::DataSet BNV_B_Chi_Photon::Thermal_Hole_E_Rate_vs_Density(
                          const double& m_chi, 
                          const Baryon& B, 
                          const bool& gen_plots) 
{
  Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI*n_B[B.label]).pow(2.0/3.0) ;

  kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

  Zaki::Vector::DataColumn m_B = m_B_ds[B.label] ;
  Zaki::Vector::DataColumn Sigma_0 = sig0_ds[B.label] ;
  Zaki::Vector::DataColumn x_F = (m_B*m_B + kF_2).sqrt() / m_B ;

  Zaki::Vector::DataSet rate ;
  rate.Reserve(2, n_B["n_tot"].Size() ) ;

  for (size_t i = 0; i < m_B.Size(); i++)
  {
    rate.AppendRow({n_B["n_tot"][i], 
      Thermal_Hole_PhaseSpace_Integral(B, m_B[i], m_chi, Sigma_0[i], x_F[i]) }) ;
  }
  
  // Converting "MeV^5" into "MeV^2/fm^3"
  rate[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "MeV^2/fm^3" into "MeV*s^-1/fm^3"
  rate[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;

  if(cgs_units)
  {
    // Converting "MeV*s^-1/fm^3" into "erg*s^-1/fm^3"
    rate[1] *= 1.6021766339999e-6 ;
    
    // Converting "erg*s^-1/fm^3" into "erg*s^-1/cm^3"
    rate[1] *= 1e39 ;
  }

  if(gen_plots)
  {
    hidden_Plot_Thermal_Hole_E_Rate_vs_Density(B, m_chi, &rate) ;
  }

  return rate ;
}

//--------------------------------------------------------------
// Plots the rate of thermal energy dumped into the NS
// due to Fermi-hole heating
// as a function of density (ds) in units of s^-1/fm^3
void BNV_B_Chi_Photon::hidden_Plot_Thermal_Hole_E_Rate_vs_Density(
                              const Baryon& B, 
                              const double& m_chi,
                              Zaki::Vector::DataSet* ds) 
{
  Zaki::String::Directory out_dir = wrk_dir + model 
                  + "/" + pulsar.GetName() + "/" + process.name
                  + "/Rate_vs_Den/" + B.short_name ;

  out_dir.Create() ;
  ds->SetWrkDir(out_dir) ;
  char file_ch[200] ;
  snprintf(file_ch, 200, "Thermal_Hole_Rate_vs_n_%s_%.0f.pdf", B.label.c_str(), m_chi) ;
  std::string title_str = "\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + " ]"  ;

  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxisLabel("$n\\, ({\\rm fm}^{-3})$") ;
  plt_pars.SetYAxisLabel("$Q_{H} ("+ GetSpecificProcess(B).TeX 
                          + ")\\, \\left({\\rm MeV}\\, s^{-1}/{\\rm fm}^3\\right)$") ;

  if(cgs_units)
  {
    plt_pars.SetYAxisLabel("$Q_{H} ("+ GetSpecificProcess(B).TeX 
                          + ")\\, \\left({\\rm erg}\\, {\\rm cm}^{-3}\\, {\rm s}^{-1}\\right)$") ;
  }

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
// Returns the rate of thermal energy dumped into the NS
// due to Fermi-hole heating
//  per unit volume as a function of radius in units of 
//  MeV s^-1/fm^3 or erg cm^-3 s^-1
Zaki::Vector::DataSet BNV_B_Chi_Photon::Thermal_Hole_E_Rate_vs_R(const double& m_chi, 
                                const Baryon& B, 
                                const bool& gen_plots) 
{
  Zaki::Vector::DataSet rate_vs_n = Thermal_Hole_E_Rate_vs_Density(m_chi, B) ;
  rate_vs_n.Interpolate(0, 1) ;

  Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;
  
  Zaki::Vector::DataSet rate_vs_r({pulsar.GetProfile()->operator[](0),
                                   rate_vs_n.Evaluate(1, n_r)}) ;

  if (gen_plots)
  {
    hidden_Plot_Thermal_Hole_E_Rate_vs_R(B, m_chi, &rate_vs_r) ;
  }
  
  return rate_vs_r ;
}

//--------------------------------------------------------------
/// Plots rate of thermal energy dumped into the NS
/// due to Fermi-hole heating 
/// as a function of radius in units of s^-1/fm^3
void BNV_B_Chi_Photon::Thermal_Hole_E_Rate_vs_R(const double& m_chi) 
{
  Thermal_Hole_E_Rate_vs_R(m_chi, neutron, true) ;
  Thermal_Hole_E_Rate_vs_R(m_chi, lambda, true) ;
}

//--------------------------------------------------------------
/// Plots rates of thermal energy dumped into the NS
/// due to Fermi-hole heating 
/// for a set of chi's in the same figure 
/// as a function of radius in units of ?
void BNV_B_Chi_Photon::Thermal_Hole_E_Rate_vs_R(const std::vector<double>& m_chi) 
{
  Thermal_Hole_E_Rate_vs_R(m_chi, neutron) ;
  Thermal_Hole_E_Rate_vs_R(m_chi, lambda) ;
}

//--------------------------------------------------------------
/// Plots rates of thermal energy dumped into the NS
/// due to Fermi-hole heating 
/// for a set of chi's in the same figure 
/// as a function of radius in units of ?
void BNV_B_Chi_Photon::Thermal_Hole_E_Rate_vs_R(const std::vector<double>& m_chi, const Baryon& B) 
{
  Zaki::Vector::DataSet rate_vs_r(pulsar.GetProfile()->operator[](0)) ;
  std::vector<std::pair<int, std::string>> labels ;
  labels.reserve(m_chi.size()) ;

  double rate_max = 0 ;
  double rate_min = 1 ;
  for (size_t i = 0; i < m_chi.size(); i++)
  {
    rate_vs_r.data_set.emplace_back(Thermal_Hole_E_Rate_vs_R(m_chi[i], B)[1]) ;

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
  snprintf(file_ch, 200, "Thermal_Hole_Rate_vs_R_%s.pdf", B.label.c_str()) ;
  std::string title_str ="\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + "]" ;

  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxisLabel("$R\\, ({\\rm km})$") ;
  plt_pars.SetYAxisLabel("$Q_H ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left({\\rm MeV}\\, s^{-1}/{\\rm fm}^3\\right)$") ;

  if(cgs_units)
  {
    plt_pars.SetYAxisLabel("$Q_H ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left({\\rm erg}\\, {\\rm cm}^{-3}\\, s^{-1}\\right)$") ;
  }

  // if (rate_min < 1e-4 * rate_max)
  // {
  //   plt_pars.SetYAxis({1e6 * rate_min, 2 * rate_max}) ;
  // }
  plt_pars.SetLegend({"lower right", 1.0, 0.0}) ;
  // plt_pars.SetXAxis({0, rate_vs_r[0].Max()}) ;
  plt_pars.SetXAxis({0, 14}) ;

  if (B.label == "10")
  {
    if(cgs_units)
    {
      plt_pars.SetYAxis({1e12, 1e21}) ;
    }
    else 
    {
      plt_pars.SetYAxis({1e-21, 3e-13}) ;
    }
  }

  if (B.label == "100")
  {
      if(cgs_units)
      {
        plt_pars.SetYAxis({1e13, 2e18}) ;
      }
      else
      {
        plt_pars.SetYAxis({1e-20, 1e-15}) ;
      }
  }

  rate_vs_r.SetPlotPars(plt_pars) ;
  rate_vs_r.SemiLogYPlot(0, labels, std::string(file_ch), title_str) ;
  
  // ------------------------------------
  //        Exporting the dataset
  // ------------------------------------
  char tmp_eps[50];
  snprintf(tmp_eps, 50, "%.2e", default_eps) ;

  std::string unit_str = "[ MeV s^-1 / fm^3 ]" ;
  if (cgs_units)
  {
    unit_str = "[ erg / cm^3 s ]" ;
  }
  
  rate_vs_r.AddHead("# Rates of thermal energy dumped into the "
                    "NS due to Fermi-hole heating per unit volume for '" 
                    + GetSpecificProcess(B).name 
                    + "' assuming eps = " 
                    + std::string(tmp_eps) 
                    + " MeV as a function of radius in '" 
                    + pulsar.GetName() + "'.\n"
                    + "# The unit for rates is in "
                    + unit_str + ".\n"
                    ) ;

  rate_vs_r.Export("Thermal_Hole_Rate_vs_R_"+ std::string(B.label.c_str()) 
                   + "_" + model + "_" 
                   + pulsar.GetName() + ".tsv") ;
  // ------------------------------------
}

//--------------------------------------------------------------
// Plots the rate of thermal energy dumped into the NS
// due to Fermi-hole heating 
// as a function of radius (ds) in units of s^-1/fm^3
void BNV_B_Chi_Photon::hidden_Plot_Thermal_Hole_E_Rate_vs_R(
                          const Baryon& B, 
                          const double& m_chi,
                          Zaki::Vector::DataSet* ds) const
{
  Zaki::String::Directory out_dir = wrk_dir + model + "/"
          + pulsar.GetName() + "/" + process.name 
          + "/Rate_vs_R/" + B.short_name ;
  out_dir.Create() ;
  ds->SetWrkDir(out_dir) ;
  char file_ch[200] ;
  snprintf(file_ch, 200, "Thermal_Hole_Rate_vs_R_%s_%.0f.pdf", B.label.c_str(), m_chi) ;
  std::string title_str ="\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + "]" ;

  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxisLabel("$R\\, ({\\rm km})$") ;
  // plt_pars.SetYAxisLabel("$\\frac{dE}{d\\tau\\, dV} ("+ GetSpecificProcess(B).TeX 
  //                       + ")\\,\\, \\left({\\rm MeV} \\, s^{-1}/{\\rm fm}^3\\right)$") ;

  plt_pars.SetYAxisLabel("$Q_H ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left({\\rm MeV} \\, s^{-1}/{\\rm fm}^3\\right)$") ;

  if(cgs_units)
  {
    plt_pars.SetYAxisLabel("$Q_H ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left({\\rm erg}\\, {\\rm cm}^-3 \\, {\\rm s}^{-1}\\right)$") ;
  }

  double rate_max = ds->operator[](1).Max() ;
  double rate_min = ds->operator[](1).Min() ;

  // if (rate_min < 1e-4 * rate_max)
  // {
  //   plt_pars.SetYAxis({1e-4 * rate_max, 2 * rate_max}) ;
  // }

  plt_pars.SetXAxis({0, 14}) ;

  if (B.label == "10")
  {
    plt_pars.SetYAxis({1e-21, 3e-13}) ;
  }

  if (B.label == "100")
  {
    plt_pars.SetYAxis({1e-21, 1.1e-18}) ;
  }

  ds->SetPlotPars(plt_pars) ;
  ds->SemiLogYPlot(0, 1, std::string(file_ch), title_str) ;
}

//--------------------------------------------------------------
void BNV_B_Chi_Photon::hidden_Plot_Thermal_Photon_E_Rate_vs_R(
                                    const Baryon& B, 
                                    const double& m_chi,
                                    Zaki::Vector::DataSet* ds) const
{
  Zaki::String::Directory out_dir = wrk_dir + model + "/"
          + pulsar.GetName() + "/" + process.name 
          + "/Rate_vs_R/" + B.short_name ;
  out_dir.Create() ;
  ds->SetWrkDir(out_dir) ;
  char file_ch[200] ;
  snprintf(file_ch, 200, "Thermal_Photon_Rate_vs_R_%s_%.0f.pdf", B.label.c_str(), m_chi) ;
  std::string title_str ="\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + "]" ;

  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxisLabel("$R\\, ({\\rm km})$") ;
  // plt_pars.SetYAxisLabel("$\\frac{dE}{d\\tau\\, dV} ("+ GetSpecificProcess(B).TeX 
  //                       + ")\\,\\, \\left({\\rm MeV} \\, s^{-1}/{\\rm fm}^3\\right)$") ;

  plt_pars.SetYAxisLabel("$Q_{\\gamma} ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left({\\rm MeV} \\, s^{-1}/{\\rm fm}^3\\right)$") ;

  if(cgs_units)
  {
    plt_pars.SetYAxisLabel("$Q_{\\gamma} ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left({\\rm erg}\\, {\\rm cm}^-3 \\, {\\rm s}^{-1}\\right)$") ;
  }

  double rate_max = ds->operator[](1).Max() ;
  double rate_min = ds->operator[](1).Min() ;

  // if (rate_min < 1e-4 * rate_max)
  // {
  //   plt_pars.SetYAxis({1e-4 * rate_max, 2 * rate_max}) ;
  // }

  plt_pars.SetXAxis({0, 14}) ;

  if (B.label == "10")
  {
    plt_pars.SetYAxis({1e-21, 3e-13}) ;
  }

  if (B.label == "100")
  {
    plt_pars.SetYAxis({1e-21, 1.1e-18}) ;
  }

  ds->SetPlotPars(plt_pars) ;
  ds->SemiLogYPlot(0, 1, std::string(file_ch), title_str) ;
}

//--------------------------------------------------------------
void BNV_B_Chi_Photon::hidden_Plot_Thermal_Photon_E_Rate_vs_Density(const Baryon& B, 
                                    const double& m_chi,
                                    Zaki::Vector::DataSet* ds)
{
  Zaki::String::Directory out_dir = wrk_dir + model 
                  + "/" + pulsar.GetName() + "/" + process.name
                  + "/Rate_vs_Den/" + B.short_name ;

  out_dir.Create() ;
  ds->SetWrkDir(out_dir) ;
  char file_ch[200] ;
  snprintf(file_ch, 200, "Thermal_Photon_Rate_vs_n_%s_%.0f.pdf", B.label.c_str(), m_chi) ;
  std::string title_str = "\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + " ]"  ;

  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxisLabel("$n\\, ({\\rm fm}^{-3})$") ;
  plt_pars.SetYAxisLabel("$Q_{\\gamma} ("+ GetSpecificProcess(B).TeX 
                          + ")\\, \\left({\\rm MeV}\\, s^{-1}/{\\rm fm}^3\\right)$") ;

  if(cgs_units)
  {
    plt_pars.SetYAxisLabel("$Q_{\\gamma} ("+ GetSpecificProcess(B).TeX 
                          + ")\\, \\left({\\rm erg}\\, {\\rm cm}^{-3}\\, {\rm s}^{-1}\\right)$") ;
  }

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
// Returns the rate of thermal energy dumped into the NS
//  per unit volume as a function of radius in units of MeV s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Chi_Photon::Thermal_Photon_E_Rate_vs_Density(
                                const double& m_chi, 
                                const Baryon& B, 
                                const bool& gen_plots)
{  
  Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI*n_B[B.label]).pow(2.0/3.0) ;

  kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

  Zaki::Vector::DataColumn m_B = m_B_ds[B.label] ;
  Zaki::Vector::DataColumn Sigma_0 = sig0_ds[B.label] ;
  Zaki::Vector::DataColumn x_F = (m_B*m_B + kF_2).sqrt() / m_B ;

  Zaki::Vector::DataSet rate ;
  rate.Reserve(2, n_B["n_tot"].Size() ) ;

  for (size_t i = 0; i < m_B.Size(); i++)
  {
    rate.AppendRow({n_B["n_tot"][i], 
      Thermal_Photon_PhaseSpace_Integral(B, m_B[i], m_chi, Sigma_0[i], x_F[i]) }) ;
  }
  
  // Converting "MeV^5" into "MeV^2/fm^3"
  rate[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "MeV^2/fm^3" into "MeV*s^-1/fm^3"
  rate[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;

  if(cgs_units)
  {
    // Converting "MeV*s^-1/fm^3" into "erg*s^-1/fm^3"
    rate[1] *= 1.6021766339999e-6 ;
    
    // Converting "erg*s^-1/fm^3" into "erg*s^-1/cm^3"
    rate[1] *= 1e39 ;
  }

  if(gen_plots)
  {
    hidden_Plot_Thermal_Photon_E_Rate_vs_Density(B, m_chi, &rate) ;
  }

  return rate ;
}

//--------------------------------------------------------------
// Returns the rate of thermal energy dumped into the NS
//  per unit volume as a function of radius in units of MeV s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Chi_Photon::Thermal_Photon_E_Rate_vs_R(const double& m_chi, 
                                const Baryon& B, 
                                const bool& gen_plots)
{
  Zaki::Vector::DataSet rate_vs_n = Thermal_Photon_E_Rate_vs_Density(m_chi, B) ;
  rate_vs_n.Interpolate(0, 1) ;

  Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;
  
  Zaki::Vector::DataSet rate_vs_r({pulsar.GetProfile()->operator[](0),
                                   rate_vs_n.Evaluate(1, n_r)}) ;

  if (gen_plots)
  {
    hidden_Plot_Thermal_Photon_E_Rate_vs_R(B, m_chi, &rate_vs_r) ;
  }
  
  return rate_vs_r ;
}
//--------------------------------------------------------------
/// Plots rate of thermal energy dumped into the NS
/// as a function of radius in units of s^-1/fm^3
void BNV_B_Chi_Photon::Thermal_Photon_E_Rate_vs_R(const double& m_chi)
{
  Thermal_Photon_E_Rate_vs_R(m_chi, neutron, true) ;
  Thermal_Photon_E_Rate_vs_R(m_chi, lambda, true) ;
}

//--------------------------------------------------------------
/// Plots rates of thermal energy dumped into the NS
/// for a set of chi's in the same figure 
/// as a function of radius in units of s^-1/fm^3
void BNV_B_Chi_Photon::Thermal_Photon_E_Rate_vs_R(const std::vector<double>& m_chi)
{
  Thermal_Photon_E_Rate_vs_R(m_chi, neutron) ;
  Thermal_Photon_E_Rate_vs_R(m_chi, lambda) ;
}

//--------------------------------------------------------------
// Plots rates of thermal energy dumped into the NS
// for a set of chi's in the same figure 
// as a function of radius in units of s^-1/fm^3
void BNV_B_Chi_Photon::Thermal_Photon_E_Rate_vs_R(const std::vector<double>& m_chi, 
                                           const Baryon& B)
{
  Zaki::Vector::DataSet rate_vs_r(pulsar.GetProfile()->operator[](0)) ;
  std::vector<std::pair<int, std::string>> labels ;
  labels.reserve(m_chi.size()) ;

  double rate_max = 0 ;
  double rate_min = 1 ;
  for (size_t i = 0; i < m_chi.size(); i++)
  {
    rate_vs_r.data_set.emplace_back(Thermal_Photon_E_Rate_vs_R(m_chi[i], B)[1]) ;

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
  snprintf(file_ch, 200, "Thermal_Photon_Rate_vs_R_%s.pdf", B.label.c_str()) ;
  std::string title_str ="\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + "]" ;

  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxisLabel("$R\\, ({\\rm km})$") ;
  plt_pars.SetYAxisLabel("$Q_{\\gamma} ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left({\\rm MeV}\\, s^{-1}/{\\rm fm}^3\\right)$") ;

  if(cgs_units)
  {
    plt_pars.SetYAxisLabel("$Q_{\\gamma} ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left({\\rm erg}\\, {\\rm cm}^{-3}\\, s^{-1}\\right)$") ;
  }

  // if (rate_min < 1e-4 * rate_max)
  // {
  //   plt_pars.SetYAxis({1e6 * rate_min, 2 * rate_max}) ;
  // }
  plt_pars.SetLegend({"lower right", 1.0, 0.0}) ;
  // plt_pars.SetXAxis({0, rate_vs_r[0].Max()}) ;
  plt_pars.SetXAxis({0, 14}) ;

  if (B.label == "10")
  {
    if(cgs_units)
    {
      plt_pars.SetYAxis({1e12, 1e21}) ;
    }
    else 
    {
      plt_pars.SetYAxis({1e-21, 3e-13}) ;
    }
  }

  if (B.label == "100")
  {
      if(cgs_units)
      {
        plt_pars.SetYAxis({1e13, 2e18}) ;
      }
      else
      {
        plt_pars.SetYAxis({1e-20, 1e-15}) ;
      }
  }

  rate_vs_r.SetPlotPars(plt_pars) ;
  rate_vs_r.SemiLogYPlot(0, labels, std::string(file_ch), title_str) ;
  
  // ------------------------------------
  //        Exporting the dataset
  // ------------------------------------
  char tmp_eps[50];
  snprintf(tmp_eps, 50, "%.2e", default_eps) ;

  std::string unit_str = "[ MeV s^-1 / fm^3 ]" ;
  if (cgs_units)
  {
    unit_str = "[ erg / cm^3 s ]" ;
  }
  
  rate_vs_r.AddHead("# Rates of thermal energy dumped into " 
                    "the NS due to photons per unit volume for '" 
                    + GetSpecificProcess(B).name 
                    + "' assuming eps = " 
                    + std::string(tmp_eps) 
                    + " MeV as a function of radius in '" 
                    + pulsar.GetName() + "'.\n"
                    + "# The unit for rates is in "
                    + unit_str + ".\n"
                    ) ;

  rate_vs_r.Export("Thermal_Photon_Rate_vs_R_"+ std::string(B.label.c_str()) 
                   + "_" + model + "_" 
                   + pulsar.GetName() + ".tsv") ;
  // ------------------------------------
}

//--------------------------------------------------------------
// This function returns the rate of thermal energy 
//  dumped into the NS (per baryon)
//  integrated over the radius and in units of [ MeV or erg per sec ]
double BNV_B_Chi_Photon::GetThermal_Photon_E_Rate(
                          const double& m_chi, 
                          const Baryon& B
                          )
{
  Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
  Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
  Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

  Zaki::Vector::DataColumn dc_rate_r = Thermal_Photon_E_Rate_vs_R(m_chi, B)[1] ;

  if (cgs_units)
  {
    // Converting "erg s^-1 cm^-3" into "erg s^-1 km^-3"
    dc_rate_r *= 1e+15 ; 
  }
  else
  {
    // Converting "MeV s^-1 fm^-3" into "MeV s^-1 km^-3"
    dc_rate_r *= 1e+54 ; 
  }

  Zaki::Vector::DataColumn decay_integrand_dc = 4*M_PI*r*r*dc_rate_r ;  
                          decay_integrand_dc /= (1 - 2*M_r / r).sqrt() ;
                          // Two redshifts: 1) Time dilation 2) Energy redshift 
                          decay_integrand_dc *= exp(2*nu_r) ; 
    Zaki::Vector::DataSet decay_integral_ds({r, decay_integrand_dc}) ;

  decay_integral_ds.Interpolate(0, 1) ;    

  double result = decay_integral_ds.Integrate(1, {r[0], r[-1]}) ;
  
  // MeV per sec per baryon
  return result / pulsar.GetSeqPoint().b ;
}

//--------------------------------------------------------------
// Plots the total rate of thermal energy dumped into the NS
// as a function of m_chi
void BNV_B_Chi_Photon::Plot_Thermal_Photon_E_Rate()
{
  // ...........................................
  // Columns are:
  // rate[0] = m_chi, rate[1] = Gamma_n, 
  // rate[2] = Gamma_lambda
  Zaki::Vector::DataSet  rate ;
  rate.Reserve(3, m_chi_vals.res) ;

  for (size_t i = 0; i < m_chi_vals.res; i++)
  {
    // double m_chi = 0 + 2.0*i ;
    double m_chi = m_chi_vals[i] ;

    double rate_n =  GetThermal_Photon_E_Rate(m_chi, neutron) ;
    double rate_lam =  GetThermal_Photon_E_Rate(m_chi, lambda) ;
    rate.AppendRow( { m_chi, rate_n, rate_lam}) ;
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
  
  std::string unit_str = "MeV/s" ;

  if (cgs_units)
  {
    unit_str = "erg/s" ;
    plt_par.SetYAxisLabel("$L_{\\gamma} \\, \\left({\\rm erg}\\, s^{-1}\\right)$") ;
  }
  else
  {
    plt_par.SetYAxisLabel("$L_{\\gamma} \\, \\left({\\rm MeV}\\, s^{-1}\\right)$") ;
  }

  plt_par.SetXAxis({0, 1400}) ;
  plt_par.SetXAxisLabel("$m_{\\chi}\\, ( {\\rm MeV} ) $") ;
  // plt_par.SetYAxis({1e-1, 3e7}) ; // For Combo
  rate.SetPlotPars(plt_par) ;
  rate.SemiLogYPlot(0, {{1, "$n$"}, {2, "$\\Lambda$"}}, 
    "Thermal_Photon_Rate_" + f_name_suff, title_str ) ;
  // ...........................................
  
  rate[0].label = "m_chi [MeV]" ;
  rate[1].label = "Gamma_E_n [" + unit_str + "]" ;
  rate[2].label = "Gamma_E_lam [" + unit_str + "]" ;
  
  char tmp_gamma[50];
  snprintf(tmp_gamma, 50, "%.2e", gamma_bnv_bin_lim) ;

  char tmp_eps[50];
  snprintf(tmp_eps, 50, "%.2e", default_eps) ;

  rate.AddHead("# Total rates of thermal energy dumped into '" 
              + pulsar.GetName() + "' in units of '" + unit_str 
              + "' per baryon for '"
              + process.name + "' assuming eps = " 
               + std::string(tmp_eps) + " MeV.\n"
               ) ;

  rate.Export("Thermal_Photon_Rate_"+ model + "_" + pulsar.GetName() + ".tsv") ;
}

//--------------------------------------------------------------
// Plots the limit on total rate of thermal energy dumped into the NS
// as a function of m_chi assuming the limit on eps value
void BNV_B_Chi_Photon::Plot_Limited_Thermal_E_Rate() const
{
  std::string in_f_name_suff = model + "_" + pulsar.GetName() + ".tsv" ;
  // ...........................................
  // The evaluated thermal E rates for 
  //  'default_eps'.
  // rate[0] = m_chi, rate[1] = Gamma_E_n, 
  // rate[2] = Gamma_E_lambda                        
  Zaki::Vector::DataSet default_rate(wrk_dir + model + "/" + pulsar.GetName() 
                                      + "/" + process.name, 
                                      "Thermal_Rate_" + in_f_name_suff) ;

  // ...........................................
  // The evaluated limits on eps 
  // eps_limits[0] = m_chi, eps_limits[1] = Gamma_n, 
  // eps_limits[2] = Gamma_lam, eps_limits[3] = eps_n,
  // eps_limits[4] = eps_lam
  Zaki::Vector::DataSet eps_limits(wrk_dir + model + "/" + pulsar.GetName() 
                                      + "/" + "Combo", 
                                      "Rate_Eps_" + in_f_name_suff) ; 

  // Scaling the thermal rates:                                    
  default_rate[1] *=  (eps_limits[3] / default_eps).pow(2) ;
  default_rate[2] *=  (eps_limits[4] / default_eps).pow(2) ;

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetLegend({"lower left", 0.0, 0.0}) ;
  plt_par.SetGrid() ;

  // rate.SetWrkDir(wrk_dir + model + "/" + pulsar.GetName() 
  //                   + "/" + process.name) ;
  std::string out_f_name_suff = model + "_" + pulsar.GetName() 
                            + ".pdf" ;

  std::string title_str = "\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + " ]"  ;
  plt_par.SetYAxisLabel("$\\frac{dE}{dt} \\, \\left({\\rm MeV}\\, s^{-1}\\right)$") ;
  plt_par.SetXAxis({0, 1400}) ;
  plt_par.SetXAxisLabel("$m_{\\chi}\\, ( {\\rm MeV} ) $") ;
  default_rate.SetPlotPars(plt_par) ;
  default_rate.SemiLogYPlot(0, {{1, "$n$"}, {2, "$\\Lambda$"}}, 
    "Limited_Thermal_Rate_" + out_f_name_suff, title_str ) ;
  // ...........................................
  default_rate[0].label = "m_chi [MeV]" ;
  default_rate[1].label = "Gamma_E_n [MeV/s]" ;
  default_rate[2].label = "Gamma_E_lam [MeV/s]" ;
  
  char tmp_gamma[50];
  snprintf(tmp_gamma, 50, "%.2e", gamma_bnv_bin_lim) ;

  char tmp_eps[50];
  snprintf(tmp_eps, 50, "%.2e", default_eps) ;

  default_rate.AddHead("# Total rates of thermal energy dumped into '" 
              + pulsar.GetName() + "' per baryon for '"
              + process.name + "' assuming the limits on eps " 
              + "from PSR binary decay rates." + ".\n"
               ) ;

  default_rate.Export("Limited_Thermal_Rate_"+ model + "_" + pulsar.GetName() + ".tsv") ;
}

//--------------------------------------------------------------
// Phase space integrand function (dim-less part)
// The total phase space integrand is the output of this function
// multiplied by:
//                  m^4 C / (pi^2)
// Here sigma_0 = - Sigma_0 / m. [Ignore this line, sign flipped!]
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

  // Old result [ WRONG ]
  // double out = x * sqrt(x*x - 1) ;
  //       out *= 1 + sigma_0*sigma_0 + 2*x*sigma_0 - mu_chi*mu_chi ;
  //       out *= (1 + sigma_0*sigma_0 + 2*x*sigma_0) * (1 + x*sigma_0 + 2 * mu_chi) 
  //               + pow(mu_chi, 2) * (1 + x*sigma_0) ;
  //       out /= x + sigma_0 ;
  //       out /= 1 + x*sigma_0 ;
  //       out /= 1 + sigma_0*sigma_0 + 2*x*sigma_0 ;

  // Result [ Updated on Mar 28, 2023 ]
  double out = sqrt(x*x - 1) ;

        out *= 1 + sigma_0*sigma_0 + 2*x*sigma_0 - mu_chi*mu_chi ;
  
        out /= pow(1 + sigma_0*sigma_0 + 2*x*sigma_0, 2) ;

        out *= (1 + sigma_0*sigma_0 + 2*x*sigma_0) 
                * (1 + x*sigma_0 + 2 * mu_chi) 
                + pow(mu_chi, 2) * (1 + x*sigma_0) ;

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
  // double eps = 1e-10 ;
  
  phase_int_mu_chi = m_chi / m ; 
  // phase_int_sigma_0 = -Sigma_0 / m ;
  phase_int_sigma_0 = Sigma_0 / m ;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
  double err, result;

  Zaki::Math::GSLFuncWrapper<BNV_B_Chi_Photon, double (BNV_B_Chi_Photon::*)(const double&)> 
    Fp(this, &BNV_B_Chi_Photon::PhaseSpace_Integrand);     

  gsl_function F = static_cast<gsl_function> (Fp) ; 

  gsl_integration_qag(&F, 1, x_F, 1e-10, 1e-10, 2000, 1, w, &result, &err) ;
  gsl_integration_workspace_free(w);

  // C
  double C = pow(B.g * default_eps * q_e, 2) / (128 * M_PI * m * m ) ;

  result *= C ;
  result *= pow(m, 4) / (M_PI*M_PI)  ;

  return result ;
}

//--------------------------------------------------------------
// // The integrated rate in the resonance region
// // The output is in s^-1
// double BNV_B_Chi_Photon::Resonance_B_Chi_Photon_Rate(const double& m_chi, 
//                                                      const Baryon& in_B,
//                                               const double& in_res_density,
//                                               const double& in_res_radius
//                                               ) 
// { 
//   // Finding the resonance baryon density from the total density
//   // First column is the total density and the secon column is the baryon's density
//   Zaki::Vector::DataSet ds_EOS_find_res_density ({EOS_ds[2], EOS_ds[2] * EOS_ds[in_B.label]}) ;
//   ds_EOS_find_res_density.Interpolate(0, 1) ;
//   double res_baryon_density = ds_EOS_find_res_density.Evaluate(1, in_res_density) ;

//   double kF_2 = pow(3*M_PI*M_PI*res_baryon_density, 2.0/3.0) ;

//   kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

//   // double m_B = m_eff_ds[in_B.label] ;
//   double m_B = m_chi ;

//   // Finding the resonance Sigma_0
//   // First column is the total density and the secon column is the baryon's self energy 
//   Zaki::Vector::DataSet ds_find_res_Vself ({V_self_E_ds[0], V_self_E_ds[in_B.label]}) ;
//   ds_find_res_Vself.Interpolate(0, 1) ;

//   double Sigma_0 = - ds_find_res_Vself.Evaluate(1, in_res_density) ;
//   // double Sigma_0 = -V_self_E_ds[in_B.label][find_res_idx] ;
//   // std::cout << "\n\n\t Sigma_0 = " << Sigma_0 << "\n" << std::flush ;

//   double x_F = sqrt(m_B*m_B + kF_2) / m_B ;
//   // std::cout << "\n\n\t x_F = " << x_F << "\n" ;

//   double q_e = Zaki::Physics::Q_E ;
//   double eps = 1e-10 ;
  
//   phase_int_mu_chi = m_chi / m_chi ; 
//   phase_int_sigma_0 = -Sigma_0 / m_chi ;

//   gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
//   double err, res_integral;

//   Zaki::Math::GSLFuncWrapper<BNV_B_Chi_Photon, double (BNV_B_Chi_Photon::*)(const double&)> 
//     Fp(this, &BNV_B_Chi_Photon::PhaseSpace_Integrand);     

//   gsl_function F = static_cast<gsl_function> (Fp) ; 

//   gsl_integration_qag(&F, 1, x_F, 1e-10, 1e-10, 2000, 1, w, &res_integral, &err) ;
//   gsl_integration_workspace_free(w);

//   // This factor doesn't include the mixing (sin(2*theta)) anymore
//   double C_2 = in_B.g * q_e / ( 16 * m_chi ) ;
//         C_2 *= C_2 ;

//   res_integral *= C_2 ;
//   res_integral *= pow(m_chi, 6) / (2*M_PI*M_PI*M_PI)  ;

//   // The radius at which the resonance occurs (in km)
//   double R_res = in_res_radius ;

  
//   // std::cout << "\n\n\t res_integral-1 = " << res_integral << "\n" ;
//   // ..........................................................
//   // Reminder: 
//   //  r     -> 0
//   // M(r)   -> 1
//   // n(r)   -> 5
//   // nu(r)  -> 6 
//   Zaki::Vector::DataSet pulsar_profile = *pulsar.GetProfile() ;
//   pulsar_profile.Interpolate(0, {1, 5, 6}) ;

//   // ..........................................................
//   // Finding he width in the mass of B for the resonace region (e.g., over R_res +- 50 m )
//   // ..........................................................

//   // The density at R_res + 50 m
//   double R_max = R_res + R_res_width ;
//   if (R_max > pulsar_profile[0][-1]) 
//     R_max = pulsar_profile[0][-1] ;

//   double density_at_R_max = pulsar_profile.Evaluate(5, R_max) ; 
//   // std::cout << "\n\n\t density_at_R_max = " << density_at_R_max << "\n" ;

//   // The density at R_res - 50 m
//   double density_at_R_min = pulsar_profile.Evaluate(5, R_res - R_res_width) ; 
//   // std::cout << "\n\n\t density_at_R_min = " << density_at_R_min << "\n" ;

//   Zaki::Vector::DataSet ds_find_res_meff ({m_eff_ds[0], m_eff_ds[in_B.label]}) ;
//   ds_find_res_meff.Interpolate(0, 1) ;

//   // [ m_eff - m_chi ] at R_res + 50 m
//   double delta_m_R_max = ds_find_res_meff.Evaluate(1, density_at_R_max) - m_chi ;
//   std::cout << "\n\n\t delta_m_R_max = " << delta_m_R_max << "\n" ;

//   // [ m_chi - m_eff ]  at R_res - 50 m
//   double delta_m_R_min = m_chi - ds_find_res_meff.Evaluate(1, density_at_R_min) ;
//   std::cout << "\n\n\t delta_m_R_min = " << delta_m_R_min << "\n" ;

//   // The derivative of baryon's mass as afunction of radius at the resonance,
//   //  i.e.,  m_B^* = m_chi :

//   Zaki::Vector::DataSet ds_meff_n ({EOS_ds[2], m_eff_ds[in_B.label]}) ;
//   ds_meff_n.Interpolate(0, 1) ;

//   Zaki::Vector::DataColumn dc_meff_r = ds_meff_n.Evaluate(1, pulsar_profile[5]) ;

//   Zaki::Vector::DataSet ds_meff_r({ pulsar_profile[0], dc_meff_r}) ;
//   ds_meff_r.Interpolate(0,1) ;
  
//   // unit is MeV per km
//   double deriv_m_vs_r = ds_meff_r.Derivative(1, R_res) ;
//   // std::cout << "\n\n\t deriv_m_vs_r = " << deriv_m_vs_r << "\n" ;
//   // ..........................................................

//   // we now multiply the 'result' which is the slowyly vaying part of the rate
//   //  but the result of integraion of a small radial range (in the resonance region)
//   //
//   //  The overal unit [ res_integral ] = [ MeV^4 . km ] 
//   res_integral *= 2 * eps *  (  atan(delta_m_R_max / (2 * eps) )  
//                               + atan(delta_m_R_min / (2 * eps) )
//                               ) / deriv_m_vs_r  ;
//   // std::cout << "\n\n\t res_integral-2 = " << res_integral << "\n" ;

//   // Converting "MeV^3" into "1/fm^3"
//   //  The overal unit [ res_integral ] = [ MeV . fm^-3 . km ] 
//   res_integral *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 
  
//   // Converting "1/fm^3" into "1/km^3"
//   //  The overal unit [ res_integral ] = [ MeV . km^-2 ] 
//   res_integral *= 1e+54 ; 

//   // Converting "MeV" into "s^-1"
//   //  The overal unit [ res_integral ] = [ s^-1 . km^-2 ]
//   res_integral *= 1e-3 * Zaki::Physics::GEV_2_S ; 


//   // The total mass at R_res (in km)
//   double M_at_R_res = pulsar_profile.Evaluate(1, R_res) ;
//   // std::cout << "\n\n\t M_at_R_res = " << M_at_R_res << "\n" ;

//   // The metric nu function at R_res
//   double nu_at_R_res = pulsar_profile.Evaluate(6, R_res) ;
//   // std::cout << "\n\n\t nu_at_R_res = " << nu_at_R_res << "\n" ;

//   // Multplying the rest of the factors in the volumn integral
//   //  The overal unit [ res_integral ] = [ s^-1 ]
//   res_integral *= 4*M_PI*R_res*R_res ;  

//   res_integral /= sqrt(1. - 2.*M_at_R_res / R_res) ; 
//   res_integral *= exp(nu_at_R_res) ; 

//   return res_integral ;
// }

//--------------------------------------------------------------
/// The decay rate per unit volume
/// Out put is in s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Chi_Photon::Rate_vs_Density(
                          const double& m_chi, 
                          const Baryon& B, 
                          const bool& gen_plots)
{  
  Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI*n_B[B.label]).pow(2.0/3.0) ;

  kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

  Zaki::Vector::DataColumn m_B = m_B_ds[B.label] ;
  Zaki::Vector::DataColumn Sigma_0 = sig0_ds[B.label] ;
  Zaki::Vector::DataColumn x_F = (m_B*m_B + kF_2).sqrt() / m_B ;

  Zaki::Vector::DataSet rate ;
  rate.Reserve(2, n_B["n_tot"].Size() ) ;

  for (size_t i = 0; i < m_B.Size(); i++)
  {
    rate.AppendRow({n_B["n_tot"][i], 
      PhaseSpace_Integral(B, m_B[i], m_chi, Sigma_0[i], x_F[i]) }) ;
  }
  
  // Converting "MeV^4" into "MeV/fm^3"
  rate[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "MeV/fm^3" into "s^-1/fm^3"
  rate[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;

  if(gen_plots)
  {
    hidden_Plot_Rate_vs_Density(B, m_chi, &rate) ;
    // Plotting the decay rate per volume as a function of density
    // Zaki::Vector::DataSet tmp_plot({eos.GetEOS(2), tmp}) ;
    // tmp_plot[0].label = "$n (fm^{-3})$";

    // tmp_plot[1].label = "$\\Gamma ({\\rm fm}^{-3} \\cdot s^{-1})$";

    // char tmp_char[50] ;
    // snprintf(tmp_char, sizeof(tmp_char), "%.0f", m_chi) ;

    // tmp_plot.SetWrkDir(wrk_dir + "/BNV_2022/results/B_Chi_Photon_Decay/"+model) ;
    // tmp_plot.Plot(0, 1, "Decay_per_V_vs_n/" + in_B.label +"_Decay_per_V_vs_n/" + 
    //                     in_B.label+"_decay_per_V_vs_n_"+std::string(tmp_char)+".pdf",
    //                     "$\\varepsilon = 10^{-10}\\, {\\rm MeV}$") ;
  }

  return rate ;
}

//--------------------------------------------------------------
/// Returns rate 
///  as a function of radius in units of s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Chi_Photon::Rate_vs_R(const double& m_chi, 
                                const Baryon& B, 
                                const bool& gen_plots) 
{
  Zaki::Vector::DataSet rate_vs_n = Rate_vs_Density(m_chi, B) ;
  rate_vs_n.Interpolate(0, 1) ;

  Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;

  // Zaki::Vector::DataSet ds_B_chi_photon_rate_n ({EOS_ds[2], B_chi_photon_rate}) ;
  
  Zaki::Vector::DataSet rate_vs_r({pulsar.GetProfile()->operator[](0),
                                   rate_vs_n.Evaluate(1, n_r)}) ;

  // // Converting "MeV^4" into "MeV/fm^3"
  // rate_vs_r[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // // Converting "MeV/fm^3" into "s^-1/fm^3"
  // rate_vs_r[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;

  if (gen_plots)
  {
    hidden_Plot_Rate_vs_R(B, m_chi, &rate_vs_r) ;
  }
  
  return rate_vs_r ;
}

//--------------------------------------------------------------
// Returns rate 
//  as a function of radius in units of s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Chi_Photon::Rate_vs_R_Slow(const double& m_chi, 
                                const Baryon& B, 
                                const bool& gen_plots) 
{
  Zaki::Vector::DataColumn kF_2 = (3*M_PI*M_PI
                        * micro_r["n_" + B.label]).pow(2.0/3.0) ;
  kF_2 /= pow(Zaki::Physics::MEV_2_INV_FM, 2) ;

  Zaki::Vector::DataColumn m_B = micro_r["m_" + B.label] ;
  Zaki::Vector::DataColumn Sigma_0 = micro_r["sig_" + B.label]  ;
  Zaki::Vector::DataColumn x_F = (m_B*m_B + kF_2).sqrt() / m_B ;

  Zaki::Vector::DataSet rate_vs_r ;
  rate_vs_r.Reserve(2, m_B.Size() ) ;

  for (size_t i = 0; i < m_B.Size(); i++)
  {
    rate_vs_r.AppendRow({micro_r[0][i], 
      PhaseSpace_Integral(B, m_B[i], m_chi, Sigma_0[i], x_F[i]) }) ;
  }
  
  // Converting "MeV^4" into "MeV/fm^3"
  rate_vs_r[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // Converting "MeV/fm^3" into "s^-1/fm^3"
  rate_vs_r[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;


  // Zaki::Vector::DataSet rate_vs_n = Rate_vs_Density(m_chi, B) ;
  // rate_vs_n.Interpolate(0, 1) ;


  // // Zaki::Vector::DataSet ds_B_chi_photon_rate_n ({EOS_ds[2], B_chi_photon_rate}) ;
  
  // Zaki::Vector::DataSet rate_vs_r({pulsar.GetProfile()->operator[](0),
  //                                  rate_vs_n.Evaluate(1, n_r)}) ;

  // // Converting "MeV^4" into "MeV/fm^3"
  // rate_vs_r[1] *= pow(Zaki::Physics::MEV_2_INV_FM, 3) ; 

  // // Converting "MeV/fm^3" into "s^-1/fm^3"
  // rate_vs_r[1] *= 1e-3 * Zaki::Physics::GEV_2_S ;

  if (gen_plots)
  {
    hidden_Plot_Rate_vs_R(B, m_chi, &rate_vs_r) ;
  }
  
  return rate_vs_r ;
}

//--------------------------------------------------------------
// Plots total rates of thermal energy dumped into the NS
// due to photons & Fermi-hole heating 
// for a set of chi's in the same figure 
// as a function of radius in units of s^-1/fm^3
void BNV_B_Chi_Photon::Thermal_Total_E_Rate_vs_R(const std::vector<double>& m_chi)
{
  Thermal_Total_E_Rate_vs_R(m_chi, neutron) ;
  Thermal_Total_E_Rate_vs_R(m_chi, lambda) ;
}

//--------------------------------------------------------------
// Plots rates of thermal energy dumped into the NS
// due to photons & Fermi-hole heating 
// for a set of chi's in the same figure 
// as a function of radius in units of ?
void BNV_B_Chi_Photon::Thermal_Total_E_Rate_vs_R(
  const std::vector<double>& m_chi, const Baryon& B)
{
  Zaki::Vector::DataSet rate_vs_r(pulsar.GetProfile()->operator[](0)) ;
  std::vector<std::pair<int, std::string>> labels ;
  labels.reserve(m_chi.size()) ;

  double rate_max = 0 ;
  double rate_min = 1 ;
  for (size_t i = 0; i < m_chi.size(); i++)
  {

    rate_vs_r.data_set.emplace_back(
      Thermal_Hole_E_Rate_vs_R(m_chi[i], B)[1]
      +
      Thermal_Photon_E_Rate_vs_R(m_chi[i], B)[1]
      ) ;

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
  snprintf(file_ch, 200, "Thermal_Total_Rate_vs_R_%s.pdf", B.label.c_str()) ;
  std::string title_str ="\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + "]" ;

  Zaki::Vector::DataSet::PlotParam plt_pars ;
  plt_pars.SetGrid() ;
  plt_pars.SetXAxisLabel("$R\\, ({\\rm km})$") ;
  plt_pars.SetYAxisLabel("$Q_T ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left({\\rm MeV}\\, s^{-1}/{\\rm fm}^3\\right)$") ;

  if(cgs_units)
  {
    plt_pars.SetYAxisLabel("$Q_T ("+ GetSpecificProcess(B).TeX 
                        + ")\\,\\, \\left({\\rm erg}\\, {\\rm cm}^{-3}\\, s^{-1}\\right)$") ;
  }

  // if (rate_min < 1e-4 * rate_max)
  // {
  //   plt_pars.SetYAxis({1e6 * rate_min, 2 * rate_max}) ;
  // }
  plt_pars.SetLegend({"lower right", 1.0, 0.0}) ;
  // plt_pars.SetXAxis({0, rate_vs_r[0].Max()}) ;
  plt_pars.SetXAxis({0, 14}) ;

  if (B.label == "10")
  {
    if(cgs_units)
    {
      plt_pars.SetYAxis({1e12, 1e21}) ;
    }
    else 
    {
      plt_pars.SetYAxis({1e-21, 3e-13}) ;
    }
  }

  if (B.label == "100")
  {
      if(cgs_units)
      {
        plt_pars.SetYAxis({1e13, 2e18}) ;
      }
      else
      {
        plt_pars.SetYAxis({1e-20, 1e-15}) ;
      }
  }

  rate_vs_r.SetPlotPars(plt_pars) ;
  rate_vs_r.SemiLogYPlot(0, labels, std::string(file_ch), title_str) ;
  
  // ------------------------------------
  //        Exporting the dataset
  // ------------------------------------
  char tmp_eps[50];
  snprintf(tmp_eps, 50, "%.2e", default_eps) ;

  std::string unit_str = "[ MeV s^-1 / fm^3 ]" ;
  if (cgs_units)
  {
    unit_str = "[ erg / cm^3 s ]" ;
  }
  
  rate_vs_r.AddHead("# Rates of thermal energy dumped into the "
                    "NS due to photons & Fermi-hole heating per unit volume for '" 
                    + GetSpecificProcess(B).name 
                    + "' assuming eps = " 
                    + std::string(tmp_eps) 
                    + " MeV as a function of radius in '" 
                    + pulsar.GetName() + "'.\n"
                    + "# The unit for rates is in "
                    + unit_str + ".\n"
                    ) ;

  rate_vs_r.Export("Thermal_Total_Rate_vs_R_"+ std::string(B.label.c_str()) 
                   + "_" + model + "_" 
                   + pulsar.GetName() + ".tsv") ;
  // ------------------------------------
}

//--------------------------------------------------------------
// This function returns the rate of thermal energy dumped into the NS
// due to Fermi-hole heating
// integrated over the radius and in units of [ MeV or erg per sec ]
double BNV_B_Chi_Photon::GetThermal_Hole_E_Rate( const double& m_chi, const Baryon& B) 
{
  Zaki::Vector::DataColumn   r = pulsar.GetProfile()->operator[](0) ;
  Zaki::Vector::DataColumn M_r = pulsar.GetProfile()->operator[](1) ;
  Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

  Zaki::Vector::DataColumn dc_rate_r = Thermal_Hole_E_Rate_vs_R(m_chi, B)[1] ;

  if (cgs_units)
  {
    // Converting "erg s^-1 cm^-3" into "erg s^-1 km^-3"
    dc_rate_r *= 1e+15 ; 
  }
  else
  {
    // Converting "MeV s^-1 fm^-3" into "MeV s^-1 km^-3"
    dc_rate_r *= 1e+54 ; 
  }

  Zaki::Vector::DataColumn decay_integrand_dc = 4*M_PI*r*r*dc_rate_r ;  
                          decay_integrand_dc /= (1 - 2*M_r / r).sqrt() ;
                          // Two redshifts: 1) Time dilation 2) Energy redshift 
                          decay_integrand_dc *= exp(2*nu_r) ; 
    Zaki::Vector::DataSet decay_integral_ds({r, decay_integrand_dc}) ;

  decay_integral_ds.Interpolate(0, 1) ;    

  double result = decay_integral_ds.Integrate(1, {r[0], r[-1]}) ;
  
  // MeV/erg per sec per baryon
  return result / pulsar.GetSeqPoint().b ;
}

//--------------------------------------------------------------
// Plots the total rate of thermal energy dumped into the NS
// as a function of m_chi
void BNV_B_Chi_Photon::Plot_Thermal_Total_E_Rate()
{
  // ...........................................
  // Columns are:
  // rate[0] = m_chi, rate[1] = Gamma_n, 
  // rate[2] = Gamma_lambda
  Zaki::Vector::DataSet  rate ;
  rate.Reserve(3, m_chi_vals.res) ;

  for (size_t i = 0; i < m_chi_vals.res; i++)
  {
    // double m_chi = 0 + 2.0*i ;
    double m_chi = m_chi_vals[i] ;

    double rate_n =  GetThermal_Photon_E_Rate(m_chi, neutron)
                    + GetThermal_Hole_E_Rate(m_chi, neutron) ;
    double rate_lam =  GetThermal_Photon_E_Rate(m_chi, lambda)
                    + GetThermal_Hole_E_Rate(m_chi, lambda) ;

    rate.AppendRow( { m_chi, rate_n, rate_lam}) ;
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
  
  std::string unit_str = "MeV/s" ;

  if (cgs_units)
  {
    unit_str = "erg/s" ;
    plt_par.SetYAxisLabel("$L_{\\rm BNV} \\, \\left({\\rm erg}\\, s^{-1}\\right)$") ;
  }
  else
  {
    plt_par.SetYAxisLabel("$L_{\\rm BNV} \\, \\left({\\rm MeV}\\, s^{-1}\\right)$") ;
  }

  plt_par.SetXAxis({0, 1400}) ;
  plt_par.SetXAxisLabel("$m_{\\chi}\\, ( {\\rm MeV} ) $") ;
  // plt_par.SetYAxis({1e-1, 3e7}) ; // For Combo
  rate.SetPlotPars(plt_par) ;
  rate.SemiLogYPlot(0, {{1, "$n$"}, {2, "$\\Lambda$"}}, 
    "Thermal_Total_Rate_" + f_name_suff, title_str ) ;
  // ...........................................
  rate[0].label = "m_chi [MeV]" ;
  rate[1].label = "Gamma_E_n [" + unit_str + "]" ;
  rate[2].label = "Gamma_E_lam [" + unit_str + "]" ;
  
  char tmp_gamma[50];
  snprintf(tmp_gamma, 50, "%.2e", gamma_bnv_bin_lim) ;

  char tmp_eps[50];
  snprintf(tmp_eps, 50, "%.2e", default_eps) ;

  rate.AddHead("# Total rates of thermal energy dumped into '" 
              + pulsar.GetName() + "' in units of '" + unit_str
              + "' per baryon for '"
              + process.name + "' assuming eps = " 
               + std::string(tmp_eps) + " MeV.\n"
               ) ;

  rate.Export("Thermal_Total_Rate_"+ model + "_" + pulsar.GetName() + ".tsv") ;
}

//--------------------------------------------------------------
// This function finds and plots (eps_max, m_chi)
// values that lead to a steady-state thermal solution
//  assuming a core temperature T_core (in kelvin)
void BNV_B_Chi_Photon::Limit_eps_From_Heating(const double& T_core) 
{
  if(!cgs_units)
  {
    Z_LOG_ERROR("Use CGS units!") ;
    return ;
  }

  pulsar.Set_T_Core(T_core) ;

  // L_nu + L_gamma [erg s^-1]
  double L_cool = pulsar.Get_Neutrino_Lumin() 
                  + pulsar.Get_Surface_Photon_Lumin() ;
  // ...........................................
  // Columns are:
  // rate[0] = m_chi, rate[1] = eps_n [MeV], 
  // rate[2] = eps_lambda [MeV]
  Zaki::Vector::DataSet  limits ;
  limits.Reserve(3, m_chi_vals.res) ;

  for (size_t i = 0; i < m_chi_vals.res; i++)
  {
    // double m_chi = 0 + 2.0*i ;
    double m_chi = m_chi_vals[i] ;

    double L_BNV_n =  GetThermal_Photon_E_Rate(m_chi, neutron)
                    + GetThermal_Hole_E_Rate(m_chi, neutron) ;
            // Remember to multiply by B
            L_BNV_n *= pulsar.GetSeqPoint().b ;

    double L_BNV_lam =  GetThermal_Photon_E_Rate(m_chi, lambda)
                    + GetThermal_Hole_E_Rate(m_chi, lambda) ;
            // Remember to multiply by B
            L_BNV_lam *= pulsar.GetSeqPoint().b ;

    // std::cout << "\n\t M_chi = " << m_chi << ", L_cool = " << L_cool ;
    // std::cout << "\t L_BNV_n = " << L_BNV_n << " " ;

    // Scale epsilon such that 'L_BNV = L_cool = L_nu + L_gamma' holds:
    double eps_n_scaled = sqrt(L_cool/L_BNV_n) * default_eps ;
    double eps_lam_scaled = sqrt(L_cool/L_BNV_lam) * default_eps ;

    limits.AppendRow( { m_chi, eps_n_scaled, eps_lam_scaled}) ;
  }

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetLegend({"lower left", 0.0, 0.0}) ;
  plt_par.SetGrid() ;

  limits.SetWrkDir(wrk_dir + model + "/" + pulsar.GetName() 
                    + "/" + process.name) ;
  std::string f_name_suff = model + "_" + pulsar.GetName() 
                            + ".pdf" ;

  std::string title_str = "\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + " ]"  ;
  
  plt_par.SetYAxisLabel("$\\varepsilon_{\\cal B} \\, \\left({\\rm MeV}\\right)$") ;

  plt_par.SetXAxis({0, 1400}) ;
  plt_par.SetXAxisLabel("$m_{\\chi}\\, ( {\\rm MeV} ) $") ;
  plt_par.SetYAxis({1e-27, 1e-4}) ; 
  limits.SetPlotPars(plt_par) ;
  limits.SemiLogYPlot(0, {{1, "$n$"}, {2, "$\\Lambda$"}}, 
    "Limi_eps_thermal_" + f_name_suff, title_str ) ;
  // ...........................................
  limits[0].label = "m_chi [MeV]" ;
  limits[1].label = "eps_n [MeV]" ;
  limits[2].label = "eps_lam [MeV]" ;

  char T_core_str[50];
  snprintf(T_core_str, 50, "%.2e", T_core) ;

  char T_surf_str[50];
  snprintf(T_surf_str, 50, "%.2e", pulsar.Get_T_Surf_Apparent()) ;

  limits.AddHead("# The mixing values corresponding to a steady-sate at T_core = '"
                  + std::string(T_core_str) +
                  + "' kelvin, i.e., apparent (red-shifted) " 
                  "surface temperature T_s = '"
                  + std::string(T_surf_str) +
                  "' kelvin for '" 
                  + pulsar.GetName() + "' in units of MeV "
                  " for different values of m_chi [MeV].\n"
               ) ;

  limits.Export("Limi_eps_thermal__"+ model + "_" + pulsar.GetName() + ".tsv") ;
}

//--------------------------------------------------------------
// This function finds and plots (eps_max, m_chi)
// values that lead to a BNV rate faster than nu-rates
//  at the transition point r_durca_thresh,
//  assuming a core temperature T_core (in kelvin)
void BNV_B_Chi_Photon::Limit_eps_From_Slowness(const double& T_core) 
{
  // if(!cgs_units)
  // {
  //   Z_LOG_ERROR("Use CGS units!") ;
  //   return ;
  // }

  pulsar.Set_T_Core(T_core) ;

  // G_nu [fm^-3 s^-1]
  Zaki::Vector::DataColumn G_nu = pulsar.Get_Neutrino_Rate()  ;

  double r_durca = pulsar.Get_R_Durca() ;
  int r_durca_idx = pulsar.GetProfile()->operator[](0).GetClosestIdx(r_durca) ;
  double G_nu_thresh = G_nu[r_durca_idx+10] ;
  double G_nu_core = G_nu[0] ;

  std::cout << "\n\t --> R_DU = " << r_durca 
            << ", G_nu[core] = "<< G_nu_core 
            << ", G_nu[R_DU] = "<< G_nu_thresh << "\t [fm^-3 s^-1].\n" ;

  // ...........................................
  // Columns are:
  // rate[0] = m_chi, rate[1] = eps_n [MeV], 
  // rate[2] = eps_lambda [MeV]
  Zaki::Vector::DataSet  limits ;
  limits.Reserve(3, m_chi_vals.res) ;

  for (size_t i = 0; i < m_chi_vals.res; i++)
  {
    // double m_chi = 0 + 2.0*i ;
    double m_chi = m_chi_vals[i] ;

    // Compare at DU threshold
    double G_BNV_n = Rate_vs_R(m_chi, neutron)[1][r_durca_idx+10] ;
    // Compare at the core
    double G_BNV_lambda = Rate_vs_R(m_chi, lambda)[1][0] ;

    // Scale epsilon such that 'L_BNV = L_cool = L_nu + L_gamma' holds:
    double eps_n_scaled = sqrt(G_nu_thresh/G_BNV_n) * default_eps ;
    double eps_lam_scaled = sqrt(G_nu_core/G_BNV_lambda) * default_eps ;

    limits.AppendRow( { m_chi, eps_n_scaled, eps_lam_scaled}) ;
  }

  Zaki::Vector::DataSet::PlotParam plt_par ;
  plt_par.SetLegend({"lower left", 0.0, 0.0}) ;
  plt_par.SetGrid() ;

  limits.SetWrkDir(wrk_dir + model + "/" + pulsar.GetName() 
                    + "/" + process.name) ;
  std::string f_name_suff = model + "_" + pulsar.GetName() 
                            + ".pdf" ;

  std::string title_str = "\n PSR" + pulsar.GetName() + "$\\qquad$"
                          + "[ " + model + " ]"  ;
  
  plt_par.SetYAxisLabel("$\\varepsilon_{\\cal B} \\, \\left({\\rm MeV}\\right)$") ;

  plt_par.SetXAxis({0, 1400}) ;
  plt_par.SetXAxisLabel("$m_{\\chi}\\, ( {\\rm MeV} ) $") ;
  plt_par.SetYAxis({1e-27, 1e-4}) ;
  limits.SetPlotPars(plt_par) ;

  char T_core_fname_str[50] ;
  snprintf(T_core_fname_str, 50, "%.1e", T_core) ;

  limits.SemiLogYPlot(0, {{1, "$n$"}, {2, "$\\Lambda$"}}, 
    "Limi_eps_slow_" + std::string(T_core_fname_str) + "_" + f_name_suff, title_str ) ;
  // ...........................................
  limits[0].label = "m_chi [MeV]" ;
  limits[1].label = "eps_n [MeV]" ;
  limits[2].label = "eps_lam [MeV]" ;

  char T_core_str[50];
  snprintf(T_core_str, 50, "%.2e", T_core) ;

  char T_surf_str[50];
  snprintf(T_surf_str, 50, "%.2e", pulsar.Get_T_Surf_Apparent()) ;

  char r_durca_str[50];
  snprintf(r_durca_str, 50, "%.3e", r_durca) ;

  limits.AddHead("# The mixing values corresponding to a BNV rate that "
                  " is competitive with MUrca at the transition radius = '"
                  + std::string(r_durca_str) + "' with a core"
                  " temperature of '"
                  + std::string(T_core_str) +
                  + "' kelvin, i.e., apparent (red-shifted) " 
                  "surface temperature T_s = '"
                  + std::string(T_surf_str) +
                  "' kelvin for '" 
                  + pulsar.GetName() + "' in units of MeV "
                  " for different values of m_chi [MeV].\n"
               ) ;

  limits.Export("Limi_eps_slow_" + std::string(T_core_fname_str)
                + "_"+ model + "_" + pulsar.GetName() + ".tsv") ;
}

//--------------------------------------------------------------
// Solves for the mixing parameter eps given the beta-imbalance
// beta-imbalance: x = eta / T
void BNV_B_Chi_Photon::Eps_From_Beta_Imbalance(const double& x) 
{
  // The dim-less function that multiples the equilibrium DUrca rates
  double f_D = 714. * x / pow(M_PI, 2) ;
        f_D += 420. * pow(x, 3) / pow(M_PI, 4) ;
        f_D += 42. * pow(x, 5) / pow(M_PI, 6) ;
        f_D /= 457. ;

  // The dim-less function that multiples the equilibrium MUrca rates
  double f_M = 14680. * x / pow(M_PI, 2) ;
        f_M += 7560. * pow(x, 3) / pow(M_PI, 4) ;
        f_M += 840. * pow(x, 5) / pow(M_PI, 6) ;
        f_M += 24. * pow(x, 7) / pow(M_PI, 6) ;
        f_M /= 11513. ;

  // The [T^{infty}_9]^5 coefficient
  double t_95_coeff = f_D * pulsar.Get_DUrca_Neutrino_Rate_T95() ;

  // The [T^{infty}_9]^7 coefficient
  double t_97_coeff = f_M * pulsar.Get_MUrca_Neutrino_Rate_T97() ;

  double m_chi = 0 ;

  // -------------------------
  //        BNV Rate 
  // -------------------------
  // The integrated BNV rate assuming eps = 10^-16 MeV
  double bnv_rate = GetRate(m_chi, neutron) ;
  
  // Remember to multiply by B
  bnv_rate *= pulsar.GetSeqPoint().b ;

  // Change the units to [ s^-1]
  bnv_rate /= 365 * 24 * 3600 ; 
  // -------------------------

  // -------------------------
  //      BNV Thermal 
  // -------------------------
  UseCGSUnits() ;

  // Unit = [ erg s^-1 ]
  double bnv_therm_n = GetThermal_Photon_E_Rate(m_chi, neutron)
                    + GetThermal_Hole_E_Rate(m_chi, neutron) ;
  // Remember to multiply by B
  bnv_therm_n *= pulsar.GetSeqPoint().b ;

  double L_phot = pulsar.Get_Surface_Photon_Lumin_T9242() ;

  // Thermal Eq: (eps / 1e-16)^2 * bnv_therm_n = L_phot * (T^infty_9)^2.42
  double eps_16_2 = L_phot / bnv_therm_n ;
  // -------------------------
  std::cout << "eps_16_2 = " << eps_16_2 << "\n" ;

  std::cout << "f_D = " << f_D << "\n" ;
  std::cout << "bnv_therm_n = " << bnv_therm_n << "\n" ;
  std::cout << "BNV Rate = " << bnv_rate << "\n" ;
  std::cout << "L_phot = " << L_phot << "\n" ;
  std::cout << "t_95_coeff = " << t_95_coeff << "\n" ;
  std::cout << "\n\t T_9 = " << pow(L_phot * bnv_rate / ( bnv_therm_n * t_95_coeff ), 0.387597) << "\n" ;

  // The [T^{infty}_9]^2.42 coefficient
  double t_9242_coeff = eps_16_2 * bnv_rate ;

  // Rate Eq. : t_95_coeff * T9^5 + t_97_coeff * T9^7 - t_9242_coeff * T9^2.42 = 0

  x_5_coeff = t_95_coeff ;
  x_7_coeff = t_97_coeff ;
  x_242_coeff = t_9242_coeff ;

  Zaki::Math::Newton solver ;

  solver.SetFunction(Eps_From_Beta_Imbalance_Eq)  ;
  solver.SetDerivative(Eps_From_Beta_Imbalance_Eq_Der)  ;

  solver.SetGuessValue(1) ;
  std::cout << "\n\n\n" ;
  std::cout << "x_5_coeff = " << x_5_coeff ;
  std::cout << "\t x_7_coeff = " << x_7_coeff ;
  std::cout << "\t x_242_coeff = " << x_242_coeff << "\n" ;

  double T_9 = solver.Solve() ;
  std::cout << "\n\n\n T_9 = " << T_9 ;
}

//--------------------------------------------------------------


//==============================================================
// 8.0e-01 	  	 9.35384360e+02 	  	  
// 8.1e-01 	  	 9.43499350e+02 	   
// 8.2e-01 	  	 9.51591970e+02 	  	 	   
// 8.3e-01 	  	 9.59662510e+02 	  	 
// 8.4e-01 	   	 9.67711230e+02 	  	  	 	  
// 8.5e-01 	   	 9.75738400e+02 	  	  	  	   
// 8.6e-01 	   	 9.83744270e+02 	  	  	  	   
// 8.7e-01 	 		 9.91729120e+02 	  	  	  	   
// 8.8e-01 	 		 9.99693170e+02   	  	   
// 8.9e-01 	 		 1.00763670e+03 	  	  	  	   
// 9.0e-01 	   	 1.01555990e+03 	  	 	  	   