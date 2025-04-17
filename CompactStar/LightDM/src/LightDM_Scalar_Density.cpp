
// -*- lsst-c++ -*-
/*
* CompactStar
* See License file at the top of the source tree.
*
* Copyright (c) 2024 Mohammadreza Zakeri
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

// -------------------------------------------------------------
//
//    This file is created for calculating the scalar density
//      due to the baryonic background
//    
//    Created on Mar 19, 2024
//
// -------------------------------------------------------------
#include <cmath>
#include <gsl/gsl_integration.h>

#include <Zaki/Util/Logger.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <Zaki/Math/GSLFuncWrapper.hpp>

#include <CompactStar/LightDM/LightDM_Scalar_Density.hpp>

using namespace CompactStar ;

double tmp_meff_val = 0 ;

double pulsar_mass = 0 ;
//==============================================================
bool MassCondition(const CompactStar::NStar& in_star)
{
  return ( pulsar_mass - 1e-4 <= in_star.GetSequence().m 
            && 
          in_star.GetSequence().m <= pulsar_mass + 1e-4  )  ;
}
//==============================================================

//==============================================================
//         LightDM_Scalar_Density class
//==============================================================
// Constructor
LightDM_Scalar_Density::LightDM_Scalar_Density() 
  : 
  neutron("neutron", "n", "n", "10", Zaki::Physics::NEUTRON_M_MEV,
           -3.82608545, 879.6),
  lambda("lambda", "\\Lambda", "Lam", "100", 
          Zaki::Physics::LAMBDA_ZERO_M_MEV, -1.22, 2.632e-10), 
  proton("proton", "p^+", "p", "11", Zaki::Physics::PROTON_M_MEV,
           5.5856946893, 1e+41),
  sigma_m("Sigma-", "\\Sigma^-", "Sig-", "110", Zaki::Physics::SIGMA_MINUS_M_MEV,
           -1.160, 1.479e-10)
  // process(in_process),
  // m_chi_vals({{0, 1400}, 700, "Linear"})
  // m_chi_vals({{0, 1100}, 220, "Linear"})
{ }

//--------------------------------------------------------------
// Destructor
LightDM_Scalar_Density::~LightDM_Scalar_Density() { }

//--------------------------------------------------------------
// Sets the EoS model
void LightDM_Scalar_Density::SetModel(const std::string& in_eos_model) 
{
  model = in_eos_model ;
}

//--------------------------------------------------------------
// Imports the EoS
void LightDM_Scalar_Density::ImportEOS(const Zaki::String::Directory& eos_dir)
{
  eos.SetWrkDir(eos_dir + model) ;
  // eos.SetWrkDir(eos_dir) ;
  eos.SetName(model) ;

  // Imports everything at once!
  eos.ImportEOS(model) ;

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
double LightDM_Scalar_Density::Scalar_Density_Integrand(const double& k)
{
  double out = 2 * tmp_meff_val / (M_PI * M_PI) ;
        out *= k*k ;
        out /= sqrt(k*k + tmp_meff_val*tmp_meff_val) ;

  return out ;
}

//--------------------------------------------------------------
// in_tot_B_dens should be input in fm^{-3}
double LightDM_Scalar_Density::Eval_Scalar_Density_vs_Baryon_Density(const double& in_tot_B_dens, 
                                                    const Baryon& B)
{

  int idx = eos.GetEOS(2).GetClosestIdx(in_tot_B_dens) ;

  tmp_meff_val =  m_B_ds[B.label][idx] ;
  Zaki::Vector::DataColumn B_dens_dc = eos.GetEOS(2) * eos.GetEOS(B.label) ;

  double B_dens = B_dens_dc[idx] ;

  double k_F = pow(3*M_PI*M_PI*B_dens, 1.0/3.0) ;
  k_F /= Zaki::Physics::MEV_2_INV_FM ;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
  double err, result;

  Zaki::Math::GSLFuncWrapper<LightDM_Scalar_Density, double (LightDM_Scalar_Density::*)(const double&)> 
    Fp(this, &LightDM_Scalar_Density::Scalar_Density_Integrand);     

  gsl_function F = static_cast<gsl_function> (Fp) ; 

  gsl_integration_qag(&F, 0, k_F, 1e-10, 1e-10, 2000, 1, w, &result, &err) ;
  gsl_integration_workspace_free(w);

  return result ;
}

//--------------------------------------------------------------
void LightDM_Scalar_Density::Export_Scalar_Density_vs_Baryon_Density(const Baryon& B)
{
  Zaki::Vector::DataSet output( 2, eos.GetEOS(2).Size()) ;
  output[0] = eos.GetEOS(2) ;
  output[1].label = "VEV(scalar)" ;

  Zaki::Vector::DataColumn B_dens_dc = eos.GetEOS(2) * eos.GetEOS(B.label) ;

  // double B_dens = B_dens_dc[idx] ;
  Zaki::Vector::DataColumn k_F = (3*M_PI*M_PI*B_dens_dc).pow(1.0/3.0) ;
  k_F /= Zaki::Physics::MEV_2_INV_FM ;

  for (size_t i = 0; i < m_B_ds[B.label].Size(); i++)
  {
    tmp_meff_val =  m_B_ds[B.label][i] ;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
    double err, result;

    Zaki::Math::GSLFuncWrapper<LightDM_Scalar_Density, double (LightDM_Scalar_Density::*)(const double&)> 
      Fp(this, &LightDM_Scalar_Density::Scalar_Density_Integrand);     

    gsl_function F = static_cast<gsl_function> (Fp) ; 

    gsl_integration_qag(&F, 0, k_F[i], 1e-10, 1e-10, 2000, 1, w, &result, &err) ;
    gsl_integration_workspace_free(w);

    output[1].vals.emplace_back(result) ;
  }

  output.SetWrkDir(wrk_dir + "/" + model) ;
  output.Export("VEV_Scalar.txt") ;

}
//--------------------------------------------------------------
// Sets the pulsar
void LightDM_Scalar_Density::SetPulsar(const CompactStar::Pulsar& in_pulsar) 
{
  pulsar = in_pulsar ;

  pulsar_mass = pulsar.GetMass().val ;
}
//--------------------------------------------------------------
void LightDM_Scalar_Density::GenSequence() const
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
    pulsar_mass = m_goal ;
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
void LightDM_Scalar_Density::FindPulsar(const bool& gen_plots) 
{ 
  bool has_hyperons = false ;
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


    // Plot_Meff_Radius() ;
    // Plot_RestEnergy_Radius() ;
    // Plot_EF_Radius() ;
    // Plot_RestE_EF_Radius(neutron) ;
    // Plot_RestE_EF_Radius(proton) ;
    // Plot_CM_E_Radius(neutron) ;

    // if (has_hyperons)
    //   Plot_CM_E_Radius(lambda) ;
    // Plot_CM_E_Radius(proton) ;
    // Plot_Estar_Radius(neutron) ;
  }

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
void LightDM_Scalar_Density::Export_Scalar_Density_vs_Radius(const Baryon& B)
{
  Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;
  Zaki::Vector::DataColumn r   = pulsar.GetProfile()->operator[](0) ;
  Zaki::Vector::DataColumn B_dens_dc = pulsar.GetProfile()->operator[](B.label) ;
                          B_dens_dc *= n_r ;
                          
  // micro dataset as a function of density
  Zaki::Vector::DataSet m_B_n({ m_B_ds[0], m_B_ds[B.label]}) ;
  m_B_n.Interpolate(0, 1) ;

  // [0]: r[km],  [1]: m_n
  Zaki::Vector::DataSet  m_B_r({r, m_B_n.Evaluate(1, n_r)}) ;

  Zaki::Vector::DataSet output( 2, r.Size()) ;
  output[0] = r ;
  output[0].label = "R[km]" ; // Changed on Jun 3, 2024!
  output[1].label = "VEV(scalar)" ;


  // double B_dens = B_dens_dc[idx] ;
  Zaki::Vector::DataColumn k_F = (3*M_PI*M_PI*B_dens_dc).pow(1.0/3.0) ;
  k_F /= Zaki::Physics::MEV_2_INV_FM ;

  for (size_t i = 0; i < m_B_r[0].Size(); i++)
  {
    tmp_meff_val =  m_B_r[1][i] ;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(2000);
    double err, result;

    Zaki::Math::GSLFuncWrapper<LightDM_Scalar_Density, double (LightDM_Scalar_Density::*)(const double&)> 
      Fp(this, &LightDM_Scalar_Density::Scalar_Density_Integrand);     

    gsl_function F = static_cast<gsl_function> (Fp) ; 

    gsl_integration_qag(&F, 0, k_F[i], 1e-10, 1e-10, 2000, 1, w, &result, &err) ;
    gsl_integration_workspace_free(w);

    output[1].vals.emplace_back(result) ;
  }

  output.SetWrkDir(wrk_dir + "/" + model) ;
  output.Export("VEV_Scalar_" + pulsar.GetName() + ".tsv") ;
}

//--------------------------------------------------------------
// Exports a table with numerical values needed for calculating the escape
//  conditions of a particle produced from the decays of a baryon B
void LightDM_Scalar_Density::Export_Escape_Params(const Baryon& B)
{
  Zaki::Vector::DataColumn r   = pulsar.GetProfile()->operator[](0) ;
  Zaki::Vector::DataColumn n_r = pulsar.GetProfile()->operator[](5) ;
  Zaki::Vector::DataColumn f_B_r = pulsar.GetProfile()->operator[](B.label) ;
  Zaki::Vector::DataColumn nu_r = pulsar.GetProfile()->operator[](6) ;

  Zaki::Vector::DataColumn kF_B = (3*M_PI*M_PI* n_r * f_B_r).pow(1.0/3.0) ;
  kF_B /= Zaki::Physics::MEV_2_INV_FM ;

  sig0_ds = eos.GetVeff() ;
  m_B_ds = eos.GetMeff() ;

  // --------------------------------------------------------
  // micro dataset as a function of density
  Zaki::Vector::DataSet m_B_n({ m_B_ds[0], m_B_ds[B.label]}) ;
  m_B_n.Interpolate(0, 1) ;
  Zaki::Vector::DataColumn  m_B_r(m_B_n.Evaluate(1, n_r)) ;
  // ........................
  Zaki::Vector::DataSet sig0_n({ sig0_ds[0], sig0_ds[B.label]}) ;
  sig0_n.Interpolate(0, 1) ;
  Zaki::Vector::DataColumn  sig0_r(sig0_n.Evaluate(1, n_r)) ;
  // --------------------------------------------------------

  Zaki::Vector::DataSet ds(5, r.Size()) ;
  ds.SetWrkDir(wrk_dir) ;

  ds[0] = r ; // Radius [km]
  ds[1] = exp(-nu_r) ; // Exp(-nu)
  ds[2] = kF_B ; // Fermi momentum [MeV]
  ds[3] = m_B_r ; // Effective mass [MeV]
  ds[4] = sig0_r ; // Self-energy [MeV]

  ds[0].label = "R[km]";
  ds[1].label = "Exp(-nu)" ;
  ds[2].label = "P_F_" + B.short_name + "[MeV]" ;
  ds[3].label = "m*_" + B.short_name + "[MeV]" ;
  ds[4].label = "Sigma0_" + B.short_name + "[MeV]" ;

  ds.Export(model + "/" + pulsar.GetName() + "/" + "Escape_Parameters.tsv") ;
}

//--------------------------------------------------------------
