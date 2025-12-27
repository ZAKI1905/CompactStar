/*
  BNV_B_Chi_Combo class
*/

#include <gsl/gsl_integration.h>

#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <Zaki/Vector/DataSet.hpp>

#include "CompactStar/Core/Pulsar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"
#include "CompactStar/Microphysics/BNV/Channels/BNV_B_Chi_Combo.hpp"

using namespace CompactStar;
namespace MicroBNVCh = Microphysics::BNV::Channels;
namespace MicroBNVInt = Microphysics::BNV::Internal;

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
MicroBNVCh::BNV_B_Chi_Combo::BNV_B_Chi_Combo()
	: BNV_Chi({"Combo", "Combo"})
{
}

//--------------------------------------------------------------
// Destructor
MicroBNVCh::BNV_B_Chi_Combo::~BNV_B_Chi_Combo() {}
//--------------------------------------------------------------
MicroBNVInt::BNV_Chi::Process MicroBNVCh::BNV_B_Chi_Combo::GetSpecificProcess(const Baryon &B) const
{
	return {"Combo_" + B.short_name, "Combo [" + B.TeX_name + "]"};
}
//--------------------------------------------------------------
// Sets the EoS model
void MicroBNVCh::BNV_B_Chi_Combo::SetModel(const std::string &in_eos_model)
{
	model = in_eos_model;

	for (auto &&i : chi_reactions)
	{
		i->SetModel(in_eos_model);
	}
}

//--------------------------------------------------------------
/// Imports the EoS
void MicroBNVCh::BNV_B_Chi_Combo::ImportEOS(const Zaki::String::Directory &eos_dir,
											const std::string &micro_model)
{
	eos.SetWrkDir(eos_dir + model);
	eos.SetName(model);

	eos.ImportGrid("eos.nb");
	eos.ImportThermo("eos.thermo");
	eos.ImportCompo("eos.compo");

	if (micro_model != "")
	{
		eos.ExtendMeffToCrust(eos_dir + micro_model, micro_model + "_m_eff.micro");
		eos.ExtendVeffToCrust(eos_dir + micro_model, micro_model + "_V.micro");
	}
	else
	{
		eos.ImportMicro("eos.micro");
	}

	sig0_ds = eos.GetVeff();
	m_B_ds = eos.GetMeff();

	// Finding the baryon's individual density
	n_B = {eos.GetEOS(2),
		   eos.GetEOS(2) * eos.GetEOS(neutron.label),
		   eos.GetEOS(2) * eos.GetEOS(lambda.label)};
	n_B[0].SetLabel("n_tot");
	n_B[1].SetLabel("10");
	n_B[2].SetLabel("100");

	for (auto &&i : chi_reactions)
	{
		i->ImportEOS(eos_dir);
	}
}

//--------------------------------------------------------------
// Attaches the pulsar
void MicroBNVCh::BNV_B_Chi_Combo::FindPulsar(const bool &gen_plots)
{
	pulsar.SetName("J0348+0432");
	pulsar.SetMass({2.01, 0.04});
	pulsar.SetSpinP({39.1226569017806, 5e-13});
	pulsar.SetSpinPDot({0.24073e-18, 4e-5});
	pulsar.SetWrkDir(wrk_dir_ + "/" + model);
	pulsar.FindProfile(model);

	if (gen_plots)
	{
		Z_LOG_WARNING("Plot generation is currently disabled inside MicroBNVCh::BNV_B_Chi_Combo::FindPulsar");
		// pulsar.PlotRelativeComposition(pulsar.GetName()  + "/" + pulsar.GetName()
		//                                 + "_RelComp_vs_R.pdf") ;
		// pulsar.PlotAbsoluteComposition(pulsar.GetName()  + "/" + pulsar.GetName()
		//                                 + "_AbsComp_vs_R.pdf") ;
		// pulsar.PlotFermiE(pulsar.GetName()  + "/" + pulsar.GetName()
		//                                 + "_EF_vs_R.pdf") ;

		// Plot_Meff_Radius() ;
		// Plot_RestEnergy_Radius() ;
		// Plot_EF_Radius() ;
	}

	for (auto &&i : chi_reactions)
	{
		i->FindPulsar(false);
	}
}

//--------------------------------------------------------------
void MicroBNVCh::BNV_B_Chi_Combo::OnWorkDirChanged(const Zaki::String::Directory &input)
{
	for (auto &&i : chi_reactions)
	{
		i->SetWrkDir(input);
	}
	// return this;
}

//--------------------------------------------------------------
void MicroBNVCh::BNV_B_Chi_Combo::AddChiReaction(BNV_Chi *bnv_chi_ptr)
{
	chi_reactions.emplace_back(bnv_chi_ptr);
}

//--------------------------------------------------------------
/// The decay rate per unit volume
/// Out put is in s^-1/fm^3
Zaki::Vector::DataSet MicroBNVCh::BNV_B_Chi_Combo::Rate_vs_Density(
	const double &m_chi,
	const Baryon &B,
	const bool &gen_plots)
{
	Zaki::Vector::DataSet rate;
	rate.AddColumn(n_B["n_tot"]);
	rate.AddColumn("rate", 0);
	for (auto &&i : chi_reactions)
	{
		rate["rate"] += i->Rate_vs_Density(m_chi, B, gen_plots)[1];
	}

	if (gen_plots)
	{
		hidden_Plot_Rate_vs_Density(B, m_chi, &rate);
	}

	return rate;
}

//--------------------------------------------------------------
/// Returns rate
///  as a function of radius in units of s^-1/fm^3
Zaki::Vector::DataSet MicroBNVCh::BNV_B_Chi_Combo::Rate_vs_R(const double &m_chi,
															 const Baryon &B,
															 const bool &gen_plots)
{
	Zaki::Vector::DataSet rate;
	rate.AddColumn(*pulsar.GetProfile()->GetRadius());
	rate.AddColumn("rate", 0);
	for (auto &&i : chi_reactions)
	{
		rate[1] += i->Rate_vs_R(m_chi, B, gen_plots)[1];
	}

	if (gen_plots)
	{
		hidden_Plot_Rate_vs_R(B, m_chi, &rate);
	}

	return rate;
}

//--------------------------------------------------------------

//--------------------------------------------------------------
//==============================================================
