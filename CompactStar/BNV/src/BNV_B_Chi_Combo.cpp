/*
  BNV_B_Chi_Combo class
*/

#include <gsl/gsl_integration.h>

#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Vector/DataSet.hpp>
#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/BNV/BNV_B_Chi_Combo.hpp"
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

//==============================================================
//                        BNV_B_Chi_Combo class
//==============================================================
// Constructor
BNV_B_Chi_Combo::BNV_B_Chi_Combo() 
  : BNV_Chi({"Combo", "Combo"})
{ }

//--------------------------------------------------------------
// Destructor
BNV_B_Chi_Combo::~BNV_B_Chi_Combo() { }
//--------------------------------------------------------------
BNV_Chi::Process BNV_B_Chi_Combo::GetSpecificProcess(const Baryon& B) const
{
  return {"Combo_" + B.short_name, "Combo [" + B.TeX_name + "]"} ;
}
//--------------------------------------------------------------
// Sets the EoS model
void BNV_B_Chi_Combo::SetModel(const std::string& in_eos_model) 
{
  model = in_eos_model ;

  for (auto &&i : chi_reactions)
  {
    i->SetModel(in_eos_model) ;
  }
}

//--------------------------------------------------------------
/// Imports the EoS
void BNV_B_Chi_Combo::ImportEOS(const Zaki::String::Directory& eos_dir, 
                                const std::string& micro_model)
{
  eos.SetWrkDir(eos_dir + model) ;
  eos.SetName(model) ;

  eos.ImportGrid("eos.nb") ;
  eos.ImportThermo("eos.thermo") ;
  eos.ImportCompo("eos.compo") ;

  if (micro_model != "")
  {
    eos.ExtendMeffToCrust(eos_dir + micro_model, micro_model + "_m_eff.micro") ;
    eos.ExtendVeffToCrust(eos_dir + micro_model, micro_model + "_V.micro") ;
  }
  else
  {
    eos.ImportMicro("eos.micro") ;
  } 

  sig0_ds = eos.GetVeff() ;
  m_B_ds = eos.GetMeff() ;

  // Finding the baryon's individual density
  n_B = { eos.GetEOS(2), 
             eos.GetEOS(2) * eos.GetEOS(neutron.label),
             eos.GetEOS(2) * eos.GetEOS(lambda.label) } ;
  n_B[0].label = "n_tot" ;
  n_B[1].label = "10" ;
  n_B[2].label = "100" ;

  for (auto &&i : chi_reactions)
  {
    i->ImportEOS(eos_dir) ;
  }
}

//--------------------------------------------------------------
// Attaches the pulsar
void BNV_B_Chi_Combo::FindPulsar(const bool& gen_plots) 
{ 
  pulsar.SetName("J0348+0432") ;
  pulsar.SetMass({2.01, 0.04}) ;
  pulsar.SetSpinP({ 39.1226569017806, 5e-13 }) ;
  pulsar.SetSpinPDot({ 0.24073e-18, 4e-5 }) ;
  pulsar.SetWrkDir( wrk_dir + "/" + model) ;
  pulsar.FindProfile(model) ;   

  if (gen_plots)
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
  }

  for (auto &&i : chi_reactions)
  {
    i->FindPulsar(false) ;
  }
}

//--------------------------------------------------------------
Prog* BNV_B_Chi_Combo::SetMemWrkDir(const Zaki::String::Directory& input) 
{
  for (auto &&i : chi_reactions)
  {
    i->SetWrkDir(input) ;
  }
  return this ;
}

//--------------------------------------------------------------
void BNV_B_Chi_Combo::AddChiReaction(BNV_Chi* bnv_chi_ptr) 
{
  chi_reactions.emplace_back(bnv_chi_ptr) ;
}

//--------------------------------------------------------------
/// The decay rate per unit volume
/// Out put is in s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Chi_Combo::Rate_vs_Density(
                          const double& m_chi, 
                          const Baryon& B, 
                          const bool& gen_plots)
{  
  Zaki::Vector::DataSet rate ;
  rate.data_set.emplace_back(n_B["n_tot"]) ;
  rate.AddColumn("rate", 0) ;
  for (auto &&i : chi_reactions)
  {
    rate["rate"] += i->Rate_vs_Density(m_chi, B, gen_plots)[1] ;
  }
  
  if(gen_plots)
  {
    hidden_Plot_Rate_vs_Density(B, m_chi, &rate) ;
  }

  return rate ;
}

//--------------------------------------------------------------
/// Returns rate 
///  as a function of radius in units of s^-1/fm^3
Zaki::Vector::DataSet BNV_B_Chi_Combo::Rate_vs_R(const double& m_chi, 
                                const Baryon& B, 
                                const bool& gen_plots) 
{
  Zaki::Vector::DataSet rate ;
  rate.data_set.emplace_back(pulsar.GetProfile()->operator[](0)) ;
  rate.AddColumn("rate", 0) ;
  for (auto &&i : chi_reactions)
  {
    rate[1] += i->Rate_vs_R(m_chi, B, gen_plots)[1] ;
  }
  
  if (gen_plots)
  {
    hidden_Plot_Rate_vs_R(B, m_chi, &rate) ;
  }
  
  return rate ;
}

//--------------------------------------------------------------

//--------------------------------------------------------------
//==============================================================
