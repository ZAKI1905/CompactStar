/*
  MixedStar class
*/

// For the timestamps
#include <chrono>
#include <ctime>

#include <Zaki/Util/Instrumentor.hpp>
// #include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Physics/Constants.hpp>

// #include <gsl/gsl_integration.h>

#include "CompactStar/Core/MixedStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

using namespace CompactStar::Core;
//==============================================================
//                        MixedStar class
//==============================================================
// Default Constructor
MixedStar::MixedStar()
	: Prog("MixedStar"),
	  core_region("Core"), mantle_region("Mantle")
// Be sure to run SurfaceIsReached if using this constructor!
{
	// integ_wrk_space = gsl_integration_workspace_alloc(1500) ;

	// FW_B_Dark.SetMemberFunc(this, &MixedStar::BaryonNumIntegrand_Dark) ;
	// F_B_Dark = static_cast<gsl_function> (FW_B_Dark) ;

	// FW_B_Visi.SetMemberFunc(this, &MixedStar::BaryonNumIntegrand) ;
	// F_B_Visi = static_cast<gsl_function> (FW_B_Visi) ;

	rot_solver.AttachMixedStar(this);
}

//--------------------------------------------------------------
// Initializes the visible dataset
void MixedStar::InitVisible(
	const TOVSolver *in_tov_solver)
{
	// Resizes the columns to '7', and
	// reserves 'radial_res' space for each column
	ds_vis.Reserve(7 + in_tov_solver->eos_tab.rho_i.size(),
				   in_tov_solver->radial_res);

	ds_vis[m_idx].SetLabel("m(km)");
	ds_vis[r_idx].SetLabel("r(km)");
	ds_vis[rho_idx].SetLabel("rho(fm^-3)");
	ds_vis[eps_idx].SetLabel("eps(km^-2)");
	ds_vis[pre_idx].SetLabel("p(km^-2)");
	ds_vis[nu_der_idx].SetLabel("nu'_v(km^-1)");
	ds_vis[nu_idx].SetLabel("nu_v");

	for (size_t i = 0; i < in_tov_solver->eos_tab.rho_i.size(); i++)
	{
		rho_i_v_idx.emplace_back(i + 7);
		ds_vis[rho_i_v_idx[i]].SetLabel(
			in_tov_solver->eos_tab.extra_labels[i]);
	}

	B_vis_integrand.Reserve(2, in_tov_solver->radial_res);

	B_vis_integrand[0].SetLabel("r(km)");
	B_vis_integrand[1].SetLabel("B_v");
}

//--------------------------------------------------------------
// Initializes the dark dataset
void MixedStar::InitDark(
	const TOVSolver *in_tov_solver)
{
	ds_dar.Reserve(7 + in_tov_solver->eos_tab_dark.rho_i.size(),
				   in_tov_solver->radial_res);

	ds_dar[m_idx].SetLabel("m_d(km)");
	ds_dar[r_idx].SetLabel("r_d(km)");
	ds_dar[rho_idx].SetLabel("rho_d(fm^-3)");
	ds_dar[eps_idx].SetLabel("eps_d(km^-2)");
	ds_dar[pre_idx].SetLabel("p_d(km^-2)");
	ds_dar[nu_der_idx].SetLabel("nu'_d(km^-1)");
	ds_dar[nu_idx].SetLabel("nu_d");

	for (size_t i = 0; i < in_tov_solver->eos_tab_dark.rho_i.size(); i++)
	{
		rho_i_d_idx.emplace_back(i + 7);
		ds_dar[rho_i_d_idx[i]].SetLabel(
			in_tov_solver->eos_tab_dark.extra_labels[i]);
	}

	B_dar_integrand.Reserve(2, in_tov_solver->radial_res);

	B_dar_integrand[0].SetLabel("r(km)");
	B_dar_integrand[1].SetLabel("B_d");
}
//--------------------------------------------------------------
// This has to be run so the class
// knows when to initialize all the splines
void MixedStar::SurfaceIsReached(
	const size_t &in_v_idx, const size_t &in_d_idx)
{
	PROFILE_FUNCTION();

	sequence.v_idx = in_v_idx;
	sequence.d_idx = in_d_idx;

	// Checking if we have a dark core or halo
	// This returns true if the dark r is smaller than visible r
	dark_core = mantle_region.SetRange(ds_dar[r_idx][-1],
									   ds_vis[r_idx][-1]);

	core_region.SetRange(ds_vis[r_idx][0],
						 mantle_region.Min());

	ds_vis.Interpolate(r_idx, {m_idx, nu_der_idx, rho_idx, eps_idx, pre_idx});

	ds_dar.Interpolate(r_idx, {m_idx, nu_der_idx, rho_idx, eps_idx, pre_idx});

	EvaluateNu();

	// ..................................................................
	//              Finding the total mass DataColumn
	// ..................................................................
	// if (dark_core)
	// {
	// 	mass_tot_dc.Resize(ds_vis[r_idx].Size());

	// 	for (size_t i = 0; i < ds_dar[r_idx].Size(); i++)
	// 	{
	// 		mass_tot_dc[i] = ds_dar[m_idx][i] + ds_vis[m_idx][i];
	// 	}
	// 	for (size_t i = ds_dar[r_idx].Size(); i < ds_vis[r_idx].Size(); i++)
	// 	{
	// 		mass_tot_dc[i] = ds_dar[m_idx][-1] + ds_vis[m_idx][i];
	// 	}
	// }
	// else
	// {
	// 	mass_tot_dc.Resize(ds_dar[r_idx].Size());

	// 	for (size_t i = 0; i < ds_vis[r_idx].Size(); i++)
	// 	{
	// 		mass_tot_dc[i] = ds_dar[m_idx][i] + ds_vis[m_idx][i];
	// 	}
	// 	for (size_t i = ds_vis[r_idx].Size(); i < ds_dar[r_idx].Size(); i++)
	// 	{
	// 		mass_tot_dc[i] = ds_dar[m_idx][i] + ds_vis[m_idx][-1];
	// 	}
	// }
	// ..................................................................
	// ..................................................................
	//      Build master-grid total profiles: r_master, M_tot, P_tot, eps_tot
	// ..................................................................
	//
	// Master grid choice:
	//   - dark_core == true  => visible extends further => master grid = ds_vis[r]
	//   - dark_core == false => dark extends further    => master grid = ds_dar[r]
	//
	// Semantics beyond a component's surface:
	//   - Mass: hold shorter component’s mass fixed at its surface value.
	//   - Pressure/Energy density: treat as zero beyond its surface (vacuum).
	//
	// These conventions match your existing baryon-integrand logic and keep the
	// physics consistent for exterior regions.

	if (dark_core)
	{
		// Master grid is the visible radius column.
		r_master_dc = ds_vis[r_idx];

		const size_t N_vis = ds_vis[r_idx].Size();
		const size_t N_dar = ds_dar[r_idx].Size();

		mass_tot_dc.Resize(N_vis);
		pre_tot_dc.Resize(N_vis);
		eps_tot_dc.Resize(N_vis);

		// Core overlap region: both sectors exist on the same index range [0, N_dar).
		for (size_t i = 0; i < N_dar; ++i)
		{
			mass_tot_dc[i] = ds_vis[m_idx][i] + ds_dar[m_idx][i];
			pre_tot_dc[i] = ds_vis[pre_idx][i] + ds_dar[pre_idx][i];
			eps_tot_dc[i] = ds_vis[eps_idx][i] + ds_dar[eps_idx][i];
		}

		// Visible-only mantle: dark sector absent -> add only visible pressure/eps,
		// while dark mass remains fixed at its surface value.
		const double M_dar_surf = ds_dar[m_idx][-1];
		for (size_t i = N_dar; i < N_vis; ++i)
		{
			mass_tot_dc[i] = ds_vis[m_idx][i] + M_dar_surf;
			pre_tot_dc[i] = ds_vis[pre_idx][i];
			eps_tot_dc[i] = ds_vis[eps_idx][i];
		}
	}
	else
	{
		// Master grid is the dark radius column.
		r_master_dc = ds_dar[r_idx];

		const size_t N_vis = ds_vis[r_idx].Size();
		const size_t N_dar = ds_dar[r_idx].Size();

		mass_tot_dc.Resize(N_dar);
		pre_tot_dc.Resize(N_dar);
		eps_tot_dc.Resize(N_dar);

		// Core overlap region: both sectors exist on the same index range [0, N_vis).
		for (size_t i = 0; i < N_vis; ++i)
		{
			mass_tot_dc[i] = ds_vis[m_idx][i] + ds_dar[m_idx][i];
			pre_tot_dc[i] = ds_vis[pre_idx][i] + ds_dar[pre_idx][i];
			eps_tot_dc[i] = ds_vis[eps_idx][i] + ds_dar[eps_idx][i];
		}

		// Dark-only mantle/halo: visible sector absent -> add only dark pressure/eps,
		// while visible mass remains fixed at its surface value.
		const double M_vis_surf = ds_vis[m_idx][-1];
		for (size_t i = N_vis; i < N_dar; ++i)
		{
			mass_tot_dc[i] = ds_dar[m_idx][i] + M_vis_surf;
			pre_tot_dc[i] = ds_dar[pre_idx][i];
			eps_tot_dc[i] = ds_dar[eps_idx][i];
		}
	}

	totals_ready_ = true;
	// ..................................................................

	// if (!totals_ready_)
	// 	Z_LOG_ERROR("MixedStar totals not initialized before building baryon integrands.");
	// ..................................................................
	//            Visible B Integrand
	// ..................................................................
	B_vis_integrand[0] = ds_vis[r_idx];
	B_vis_integrand[1] = ds_vis[r_idx].pow(2);
	B_vis_integrand[1] *= 4 * M_PI * ds_vis[rho_idx];

	if (dark_core)
		B_vis_integrand[1] /= (1. - 2. * mass_tot_dc / ds_vis[r_idx]).sqrt();
	else
		B_vis_integrand[1] /= (1. - 2. * mass_tot_dc.GetSubSet(0, ds_vis[r_idx].Size() - 1) / ds_vis[r_idx]).sqrt();

	// converting fm^{-3} to km^{-3}
	B_vis_integrand[1] *= pow(10, 54);

	B_vis_integrand.Interpolate(0, 1);
	// ..................................................................
	//            Dark B Integrand
	// ..................................................................
	B_dar_integrand[0] = ds_dar[r_idx];
	B_dar_integrand[1] = ds_dar[r_idx].pow(2);
	B_dar_integrand[1] *= 4 * M_PI * ds_dar[rho_idx];

	if (!dark_core)
		B_dar_integrand[1] /= (1. - 2. * mass_tot_dc / ds_dar[r_idx]).sqrt();
	else
		B_dar_integrand[1] /= (1. - 2. * mass_tot_dc.GetSubSet(0, ds_dar[r_idx].Size() - 1) / ds_dar[r_idx]).sqrt();

	// converting fm^{-3} to km^{-3}
	B_dar_integrand[1] *= pow(10, 54);

	B_dar_integrand.Interpolate(0, 1);
	// ..................................................................

	// Units must be converted from km
	sequence.v.ec = ds_vis[eps_idx][0] * Zaki::Physics::INV_FM4_2_G_CM3 /
					Zaki::Physics::INV_FM4_2_INV_KM2;
	sequence.v.m = ds_vis[m_idx][-1] / Zaki::Physics::SUN_M_KM;
	sequence.v.r = ds_vis[r_idx][-1];
	sequence.v.pc = ds_vis[pre_idx][0] * Zaki::Physics::INV_FM4_2_Dyn_CM2 /
					Zaki::Physics::INV_FM4_2_INV_KM2;
	// sequence.v.b  = Find_BaryonNum_Visible() ;
	sequence.v.b = B_vis_integrand.Integrate(1, {ds_vis[r_idx][0], ds_vis[r_idx][-1]});

	sequence.d.ec = ds_dar[eps_idx][0] * Zaki::Physics::INV_FM4_2_G_CM3 /
					Zaki::Physics::INV_FM4_2_INV_KM2;
	sequence.d.m = ds_dar[m_idx][-1] / Zaki::Physics::SUN_M_KM;
	sequence.d.r = ds_dar[r_idx][-1];
	sequence.d.pc = ds_dar[pre_idx][0] * Zaki::Physics::INV_FM4_2_Dyn_CM2 /
					Zaki::Physics::INV_FM4_2_INV_KM2;
	// sequence.d.b  = Find_BaryonNum_Dark() ;
	sequence.d.b = B_dar_integrand.Integrate(1, {ds_dar[r_idx][0], ds_dar[r_idx][-1]});

	sequence.d.I = 0;
	sequence.v.I = Find_MomInertia();
}
//--------------------------------------------------------------
// May 1, 2022: Not needed anymore, check the constructors!
// Reserving the needed space for all datacolumns
// void MixedStar::Reserve(const size_t& space_size)
// {
//   for (size_t c = 0; c < ds_vis.Dim().size() ; c++)
//   {
//     ds_vis[c].Reserve(space_size) ;
//   }

//   for (size_t c = 0; c < ds_dar.Dim().size() ; c++)
//   {
//     ds_dar[c].Reserve(space_size) ;
//   }

// nu_der.Reserve(space_size) ;
// nu.Reserve(space_size) ;

// nu_der_d.Reserve(space_size) ;
// nu_der_v.Reserve(space_size) ;
// mass.Reserve(space_size) ;
// mass_d.Reserve(space_size) ;
// radius.Reserve(space_size) ;
// radius_d.Reserve(space_size) ;
// rho.Reserve(space_size) ;
// rho_d.Reserve(space_size) ;
// eps.Reserve(space_size) ;
// eps_d.Reserve(space_size) ;
// press.Reserve(space_size) ;
// press_d.Reserve(space_size) ;

// }
//--------------------------------------------------------------
// Appends tov points to the MixedStar in the core region
void MixedStar::Append_Core(const TOVPoint &in_tov,
							const TOVPoint &in_dark_tov)
{
	// std::cout << radius_d.Size()<< ", r = " << in_tov.r << "\n" ;
	ds_vis[r_idx].PushBack(in_tov.r);																		  // in km
	ds_vis[m_idx].PushBack(Zaki::Physics::SUN_M_KM * in_tov.m);												  // in km
	ds_vis[rho_idx].PushBack(in_tov.rho);																	  // in fm^{-3}
	ds_vis[eps_idx].PushBack(in_tov.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3);	  // in km^{-2}
	ds_vis[pre_idx].PushBack(in_tov.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2); // in km^{-2}
	ds_vis[nu_der_idx].PushBack(in_tov.nu_der * 1e+5);														  // convert 1/cm to 1/km

	for (size_t i = 0; i < in_tov.rho_i.size(); i++)
	{
		// in fm^{-3}
		ds_vis[rho_i_v_idx[i]].PushBack(in_tov.rho_i[i]);
	}

	ds_dar[r_idx].PushBack(in_dark_tov.r);																		   // in km
	ds_dar[m_idx].PushBack(Zaki::Physics::SUN_M_KM * in_dark_tov.m);											   // in km
	ds_dar[rho_idx].PushBack(in_dark_tov.rho);																	   // in fm^{-3}
	ds_dar[eps_idx].PushBack(in_dark_tov.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3);   // in km^{-2}
	ds_dar[pre_idx].PushBack(in_dark_tov.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2); // in km^{-2}
	ds_dar[nu_der_idx].PushBack(in_dark_tov.nu_der * 1e+5);														   // convert 1/cm to 1/km

	for (size_t j = 0; j < in_dark_tov.rho_i.size(); j++)
	{
		// in fm^{-3}
		ds_dar[rho_i_d_idx[j]].PushBack(in_dark_tov.rho_i[j]);
	}
}

//--------------------------------------------------------------
// Appends tov points to the MixedStar in the in the mantle
// from the dark solutions
void MixedStar::Append_Dark_Mantle(const TOVPoint &in_dark_tov)
{
	// std::cout << radius_d.Size()<< ", r = " << in_dark_tov.r << "\n" ;

	// if ( radius_d[-1] == in_dark_tov.r )
	// std::cout << "\n\tAppend_Dark_Mantle: "
	//           << radius_d.Size()<< ", r = " << in_dark_tov.r << "\n" ;

	ds_dar[r_idx].PushBack(in_dark_tov.r);																		   // in km
	ds_dar[m_idx].PushBack(Zaki::Physics::SUN_M_KM * in_dark_tov.m);											   // in km
	ds_dar[rho_idx].PushBack(in_dark_tov.rho);																	   // in fm^{-3}
	ds_dar[eps_idx].PushBack(in_dark_tov.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3);   // in km^{-2}
	ds_dar[pre_idx].PushBack(in_dark_tov.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2); // in km^{-2}
	ds_dar[nu_der_idx].PushBack(in_dark_tov.nu_der * 1e+5);														   // convert 1/cm to 1/km

	for (size_t j = 0; j < in_dark_tov.rho_i.size(); j++)
	{
		// in fm^{-3}
		ds_dar[rho_i_d_idx[j]].PushBack(in_dark_tov.rho_i[j]);
	}
}
//--------------------------------------------------------------
// Appends tov points to the MixedStar in the mantle
// from the visible solutions
void MixedStar::Append_Visible_Mantle(const TOVPoint &in_tov)
{
	ds_vis[r_idx].PushBack(in_tov.r);																		  // in km
	ds_vis[m_idx].PushBack(Zaki::Physics::SUN_M_KM * in_tov.m);												  // in km
	ds_vis[rho_idx].PushBack(in_tov.rho);																	  // in fm^{-3}
	ds_vis[eps_idx].PushBack(in_tov.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3);	  // in km^{-2}
	ds_vis[pre_idx].PushBack(in_tov.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2); // in km^{-2}
	ds_vis[nu_der_idx].PushBack(in_tov.nu_der * 1e+5);														  // convert 1/cm to 1/km

	for (size_t i = 0; i < in_tov.rho_i.size(); i++)
	{
		// in fm^{-3}
		ds_vis[rho_i_v_idx[i]].PushBack(in_tov.rho_i[i]);
	}
}
//--------------------------------------------------------------
// Constructor from TOV Solutions
MixedStar::MixedStar(
	const size_t &in_v_idx, const size_t &in_d_idx,
	const std::vector<TOVPoint> &in_tov,
	const std::vector<TOVPoint> &in_dark_tov)
	: Prog("MixedStar"),
	  core_region("Core"), mantle_region("Mantle"),
	  ds_vis(7 + in_tov[0].rho_i.size(), in_tov.size()),
	  ds_dar(7 + in_dark_tov[0].rho_i.size(), in_dark_tov.size()),
	  B_vis_integrand(2, in_tov.size()),
	  B_dar_integrand(2, in_dark_tov.size())
{
	sequence.v_idx = in_v_idx;
	sequence.d_idx = in_d_idx;

	B_vis_integrand[0].SetLabel("r(km)");
	B_vis_integrand[1].SetLabel("B_v");
	B_dar_integrand[0].SetLabel("r(km)");
	B_dar_integrand[1].SetLabel("B_d");

	// integ_wrk_space = gsl_integration_workspace_alloc(1500) ;

	// FW_B_Dark.SetMemberFunc(this, &MixedStar::BaryonNumIntegrand_Dark) ;
	// F_B_Dark = static_cast<gsl_function> (FW_B_Dark) ;

	// FW_B_Visi.SetMemberFunc(this, &MixedStar::BaryonNumIntegrand) ;
	// F_B_Visi = static_cast<gsl_function> (FW_B_Visi) ;

	rot_solver.AttachMixedStar(this);

	// This following assumes we have a dark core
	// so we should swap the components if we don't
	// have a dark core!
	// mantle_region.r = {  in_dark_tov[in_dark_tov.size() - 1].r,
	//             in_tov[in_tov.size() - 1].r } ;

	// Checking if we have a dark core or halo
	// This returns true if the dark r is smaller than visible r
	dark_core = mantle_region.SetRange(
		in_dark_tov[in_dark_tov.size() - 1].r,
		in_tov[in_tov.size() - 1].r);

	core_region.SetRange(in_tov[0].r,
						 mantle_region.Min());

	ds_vis[m_idx].SetLabel("m");
	ds_vis[r_idx].SetLabel("r");
	ds_vis[rho_idx].SetLabel("rho");
	ds_vis[eps_idx].SetLabel("eps");
	ds_vis[pre_idx].SetLabel("p");
	ds_vis[nu_der_idx].SetLabel("nu'_v");
	ds_vis[nu_idx].SetLabel("nu_v");

	for (size_t j = 0; j < in_tov[0].rho_i.size(); j++)
	{
		rho_i_v_idx.emplace_back(6 + j);
		ds_vis[rho_i_v_idx[j]].SetLabel("X"); // !!!
	}

	ds_dar[m_idx].SetLabel("m_d");
	ds_dar[r_idx].SetLabel("r_d");
	ds_dar[rho_idx].SetLabel("rho_d");
	ds_dar[eps_idx].SetLabel("eps_d");
	ds_dar[pre_idx].SetLabel("p_d");
	ds_dar[nu_der_idx].SetLabel("nu'_d");
	ds_dar[nu_idx].SetLabel("nu_d");

	for (size_t j = 0; j < in_dark_tov[0].rho_i.size(); j++)
	{
		rho_i_d_idx.emplace_back(6 + j);
		ds_dar[rho_i_d_idx[j]].SetLabel("Y"); // !!!
	}

	// ...............................................
	for (auto &&i : in_tov)
	{
		ds_vis[r_idx].PushBack(i.r);																		 // in km
		ds_vis[m_idx].PushBack(Zaki::Physics::SUN_M_KM * i.m);												 // in km
		ds_vis[rho_idx].PushBack(i.rho);																	 // in fm^{-3}
		ds_vis[eps_idx].PushBack(i.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3);	 // in km^{-2}
		ds_vis[pre_idx].PushBack(i.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2); // in km^{-2}
		ds_vis[nu_der_idx].PushBack(i.nu_der * 1e+5);														 // convert 1/cm to 1/km

		for (size_t j = 0; j < i.rho_i.size(); j++)
		{
			ds_vis[rho_i_v_idx[j]].PushBack(i.rho_i[j]);
		}
	}

	for (auto &&i : in_dark_tov)
	{
		ds_dar[r_idx].PushBack(i.r);																		 // in km
		ds_dar[m_idx].PushBack(Zaki::Physics::SUN_M_KM * i.m);												 // in km
		ds_dar[rho_idx].PushBack(i.rho);																	 // in fm^{-3}
		ds_dar[eps_idx].PushBack(i.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3);	 // in km^{-2}
		ds_dar[pre_idx].PushBack(i.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2); // in km^{-2}
		ds_dar[nu_der_idx].PushBack(i.nu_der * 1e+5);														 // convert 1/cm to 1/km

		for (size_t j = 0; j < i.rho_i.size(); j++)
		{
			ds_dar[rho_i_d_idx[j]].PushBack(i.rho_i[j]);
		}
	}

	// for(size_t i = 0 ; i < in_tov.size() ; ++i)
	// {
	//   ds_vis[r_idx][i]    = in_tov[i].r ; // in km
	//   ds_vis[m_idx][i]    = Zaki::Physics::SUN_M_KM*in_tov[i].m ; // in km
	//   ds_vis[rho_idx][i]  = in_tov[i].rho ; // in fm^{-3}
	//   ds_vis[eps_idx][i]  = in_tov[i].e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3 ; // in km^{-2}
	//   ds_vis[pre_idx][i]  = in_tov[i].p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2 ; // in km^{-2}
	//   ds_vis[nu_der_idx][i] = in_tov[i].nu_der * 1e+5 ; // convert 1/cm to 1/km
	// }

	//   for(size_t i = 0 ; i < in_dark_tov.size() ; ++i)
	// {
	//   ds_dar[r_idx][i]    = in_dark_tov[i].r ; // in km
	//   ds_dar[m_idx][i]    = Zaki::Physics::SUN_M_KM*in_dark_tov[i].m ; // in km
	//   ds_dar[rho_idx][i]  = in_dark_tov[i].rho ; // in fm^{-3}
	//   ds_dar[eps_idx][i]  = in_dark_tov[i].e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3 ; // in km^{-2}
	//   ds_dar[pre_idx][i]  = in_dark_tov[i].p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2 ; // in km^{-2}
	//   ds_dar[nu_der_idx][i] = in_dark_tov[i].nu_der * 1e+5 ; // convert 1/cm to 1/km
	// }

	// gsl_spline_init (mass_r_spline, &radius.vals[0], &mass.vals[0], in_tov.size()) ;
	// gsl_spline_init (rho_r_spline, &radius.vals[0], &rho.vals[0], in_tov.size()) ;
	// gsl_spline_init (eps_r_spline, &radius.vals[0], &eps.vals[0], in_tov.size()) ;
	// gsl_spline_init (press_r_spline, &radius.vals[0], &press.vals[0], in_tov.size()) ;

	// gsl_spline_init (mass_r_d_spline, &radius_d.vals[0], &mass_d.vals[0], in_dark_tov.size()) ;
	// gsl_spline_init (rho_r_d_spline, &radius_d.vals[0], &rho_d.vals[0], in_dark_tov.size()) ;
	// gsl_spline_init (eps_r_d_spline, &radius_d.vals[0], &eps_d.vals[0], in_dark_tov.size()) ;
	// gsl_spline_init (press_r_d_spline, &radius_d.vals[0], &press_d.vals[0], in_dark_tov.size()) ;

	// if(dark_core)
	// {
	//   gsl_spline_init (nu_der_r_spline, &radius.vals[0], &nu_der.vals[0], in_tov.size()) ;
	// }
	// else
	// {
	//   gsl_spline_init (nu_der_r_spline, &radius_d.vals[0], &nu_der.vals[0], in_dark_tov.size()) ;
	// }

	ds_vis.Interpolate(r_idx, {m_idx, nu_der_idx, rho_idx, eps_idx, pre_idx});
	ds_dar.Interpolate(r_idx, {m_idx, nu_der_idx, rho_idx, eps_idx, pre_idx});

	EvaluateNu();

	// // ..................................................................
	// //              Finding the total mass DataColumn
	// // ..................................................................
	// if (dark_core)
	// {
	// 	mass_tot_dc.Resize(ds_vis[r_idx].Size());

	// 	for (size_t i = 0; i < ds_dar[r_idx].Size(); i++)
	// 	{
	// 		mass_tot_dc[i] = ds_dar[m_idx][i] + ds_vis[m_idx][i];
	// 	}
	// 	for (size_t i = ds_dar[r_idx].Size(); i < ds_vis[r_idx].Size(); i++)
	// 	{
	// 		mass_tot_dc[i] = ds_dar[m_idx][-1] + ds_vis[m_idx][i];
	// 	}
	// }
	// else
	// {
	// 	mass_tot_dc.Resize(ds_dar[r_idx].Size());

	// 	for (size_t i = 0; i < ds_vis[r_idx].Size(); i++)
	// 	{
	// 		mass_tot_dc[i] = ds_dar[m_idx][i] + ds_vis[m_idx][i];
	// 	}
	// 	for (size_t i = ds_vis[r_idx].Size(); i < ds_dar[r_idx].Size(); i++)
	// 	{
	// 		mass_tot_dc[i] = ds_dar[m_idx][i] + ds_vis[m_idx][-1];
	// 	}
	// }
	// // ..................................................................
	// ..................................................................
	//      Build master-grid total profiles: r_master, M_tot, P_tot, eps_tot
	// ..................................................................
	//
	// Master grid choice:
	//   - dark_core == true  => visible extends further => master grid = ds_vis[r]
	//   - dark_core == false => dark extends further    => master grid = ds_dar[r]
	//
	// Semantics beyond a component's surface:
	//   - Mass: hold shorter component’s mass fixed at its surface value.
	//   - Pressure/Energy density: treat as zero beyond its surface (vacuum).
	//
	// These conventions match your existing baryon-integrand logic and keep the
	// physics consistent for exterior regions.

	if (dark_core)
	{
		// Master grid is the visible radius column.
		r_master_dc = ds_vis[r_idx];

		const size_t N_vis = ds_vis[r_idx].Size();
		const size_t N_dar = ds_dar[r_idx].Size();

		mass_tot_dc.Resize(N_vis);
		pre_tot_dc.Resize(N_vis);
		eps_tot_dc.Resize(N_vis);

		// Core overlap region: both sectors exist on the same index range [0, N_dar).
		for (size_t i = 0; i < N_dar; ++i)
		{
			mass_tot_dc[i] = ds_vis[m_idx][i] + ds_dar[m_idx][i];
			pre_tot_dc[i] = ds_vis[pre_idx][i] + ds_dar[pre_idx][i];
			eps_tot_dc[i] = ds_vis[eps_idx][i] + ds_dar[eps_idx][i];
		}

		// Visible-only mantle: dark sector absent -> add only visible pressure/eps,
		// while dark mass remains fixed at its surface value.
		const double M_dar_surf = ds_dar[m_idx][-1];
		for (size_t i = N_dar; i < N_vis; ++i)
		{
			mass_tot_dc[i] = ds_vis[m_idx][i] + M_dar_surf;
			pre_tot_dc[i] = ds_vis[pre_idx][i];
			eps_tot_dc[i] = ds_vis[eps_idx][i];
		}
	}
	else
	{
		// Master grid is the dark radius column.
		r_master_dc = ds_dar[r_idx];

		const size_t N_vis = ds_vis[r_idx].Size();
		const size_t N_dar = ds_dar[r_idx].Size();

		mass_tot_dc.Resize(N_dar);
		pre_tot_dc.Resize(N_dar);
		eps_tot_dc.Resize(N_dar);

		// Core overlap region: both sectors exist on the same index range [0, N_vis).
		for (size_t i = 0; i < N_vis; ++i)
		{
			mass_tot_dc[i] = ds_vis[m_idx][i] + ds_dar[m_idx][i];
			pre_tot_dc[i] = ds_vis[pre_idx][i] + ds_dar[pre_idx][i];
			eps_tot_dc[i] = ds_vis[eps_idx][i] + ds_dar[eps_idx][i];
		}

		// Dark-only mantle/halo: visible sector absent -> add only dark pressure/eps,
		// while visible mass remains fixed at its surface value.
		const double M_vis_surf = ds_vis[m_idx][-1];
		for (size_t i = N_vis; i < N_dar; ++i)
		{
			mass_tot_dc[i] = ds_dar[m_idx][i] + M_vis_surf;
			pre_tot_dc[i] = ds_dar[pre_idx][i];
			eps_tot_dc[i] = ds_dar[eps_idx][i];
		}
	}

	totals_ready_ = true;
	// ..................................................................

	// ..................................................................
	//            Visible B Integrand
	// ..................................................................
	B_vis_integrand[0] = ds_vis[r_idx];
	B_vis_integrand[1] = ds_vis[r_idx].pow(2);
	B_vis_integrand[1] *= 4 * M_PI * ds_vis[rho_idx];

	if (dark_core)
		B_vis_integrand[1] /= (1. - 2. * mass_tot_dc / ds_vis[r_idx]).sqrt();
	else
		B_vis_integrand[1] /= (1. - 2. * mass_tot_dc.GetSubSet(0, ds_vis[r_idx].Size() - 1) / ds_vis[r_idx]).sqrt();

	// converting fm^{-3} to km^{-3}
	B_vis_integrand[1] *= pow(10, 54);

	B_vis_integrand.Interpolate(0, 1);
	// ..................................................................
	//            Dark B Integrand
	// ..................................................................
	B_dar_integrand[0] = ds_dar[r_idx];
	B_dar_integrand[1] = ds_dar[r_idx].pow(2);
	B_dar_integrand[1] *= 4 * M_PI * ds_dar[rho_idx];

	if (!dark_core)
		B_dar_integrand[1] /= (1. - 2. * mass_tot_dc / ds_dar[r_idx]).sqrt();
	else
		B_dar_integrand[1] /= (1. - 2. * mass_tot_dc.GetSubSet(0, ds_dar[r_idx].Size() - 1) / ds_dar[r_idx]).sqrt();

	// converting fm^{-3} to km^{-3}
	B_dar_integrand[1] *= pow(10, 54);

	B_dar_integrand.Interpolate(0, 1);
	// ..................................................................

	sequence.v.ec = ds_vis[eps_idx][0] * Zaki::Physics::INV_FM4_2_G_CM3 /
					Zaki::Physics::INV_FM4_2_INV_KM2;
	sequence.v.m = ds_vis[m_idx][-1] / Zaki::Physics::SUN_M_KM;
	sequence.v.r = ds_vis[r_idx][-1];
	sequence.v.pc = ds_vis[pre_idx][0] * Zaki::Physics::INV_FM4_2_Dyn_CM2 /
					Zaki::Physics::INV_FM4_2_INV_KM2;
	// sequence.v.b  = Find_BaryonNum_Visible() ;
	sequence.v.b = B_vis_integrand.Integrate(1, {ds_vis[r_idx][0], ds_vis[r_idx][-1]});

	sequence.v.I = Find_MomInertia();

	sequence.d.ec = ds_dar[eps_idx][0] * Zaki::Physics::INV_FM4_2_G_CM3 /
					Zaki::Physics::INV_FM4_2_INV_KM2;
	sequence.d.m = ds_dar[m_idx][-1] / Zaki::Physics::SUN_M_KM;
	sequence.d.r = ds_dar[r_idx][-1];
	sequence.d.pc = ds_dar[pre_idx][0] * Zaki::Physics::INV_FM4_2_Dyn_CM2 /
					Zaki::Physics::INV_FM4_2_INV_KM2;
	// sequence.d.b  = Find_BaryonNum_Dark() ;
	sequence.d.b = B_dar_integrand.Integrate(1, {ds_dar[r_idx][0], ds_dar[r_idx][-1]});

	sequence.d.I = 0;
}
//--------------------------------------------------------------
// Sets the work directory for the member objects
void MixedStar::OnWorkDirChanged(
	const Zaki::String::Directory &in_dir)
{
	ds_vis.SetWrkDir(in_dir);
	ds_dar.SetWrkDir(in_dir);

	// return this ;
}

//--------------------------------------------------------------
// Similar to the destructor
void MixedStar::Reset()
{
	ds_dar.ClearRows();
	ds_vis.ClearRows();

	B_vis_integrand.ClearRows();
	B_dar_integrand.ClearRows();

	sequence.Reset();

	totals_ready_ = false;
	r_master_dc.Clear();
	pre_tot_dc.Clear();
	eps_tot_dc.Clear();
	mass_tot_dc.Clear();
}

//--------------------------------------------------------------
// Destructor
MixedStar::~MixedStar()
{
	// Memory leak fixed on May 6, 2022
	// gsl_integration_workspace_free(integ_wrk_space) ;
}

//--------------------------------------------------------------
// Returns the nu_der value given the radius input
// r is in km !
// double MixedStar::GetNuDerSpline_Vis(const double& in_r)
// {
//   double tmp = gsl_spline_eval(nu_der_r_spline, in_r, mixed_r_accel);
//   return tmp ;
// }

// //--------------------------------------------------------------
// // Returns the nu_der value given the radius input
// // r is in km !
// double MixedStar::GetNuDerSpline_Dar(const double& in_r)
// {
//   double tmp = gsl_spline_eval(nu_der_r_spline, in_r, mixed_r_accel);
//   return tmp ;
// }

//--------------------------------------------------------------
// Evaluate the metric function
void MixedStar::EvaluateNu()
{
	PROFILE_FUNCTION();
	// -----------------------------------
	// Integrate to find nu(r) :
	// -----------------------------------
	ds_vis[nu_idx].Resize(ds_vis[r_idx].Size());
	ds_dar[nu_idx].Resize(ds_dar[r_idx].Size());

	// Boundary condition
	double nu_at_R = 0.5 * log(
							   1 - 2 * (ds_vis[m_idx][-1] + ds_dar[m_idx][-1]) / mantle_region.Max());

	// if (dark_core)
	// {
	// ds_vis.AddColumn("nu") ;

	// Calculate the surface term first to find out delta_nu_r
	ds_vis[nu_idx][-1] = ds_vis.Integrate(nu_der_idx, {ds_vis[r_idx][0], ds_vis[r_idx][-1]});
	double delta_nu_r = ds_vis[nu_idx][-1] - nu_at_R;

	ds_vis[nu_idx][-1] = nu_at_R;

	// We don't need to calculate the surface term anymore,
	//  therefore we have 'i < ds_vis[r_idx].Size()-1'
	for (size_t i = 0; i < ds_vis[r_idx].Size() - 1; ++i)
	{
		ds_vis[nu_idx][i] = ds_vis.Integrate(nu_der_idx, {ds_vis[r_idx][0], ds_vis[r_idx][i]}) - delta_nu_r;
	}

	ds_vis.Interpolate(0, nu_idx);
	// }
	// else
	// {
	// ds_dar.AddColumn("nu") ;

	// Calculate the surface term first to find out delta_nu_r
	ds_dar[nu_idx][-1] = ds_dar.Integrate(nu_der_idx, {ds_dar[r_idx][0], ds_dar[r_idx][-1]});
	delta_nu_r = ds_dar[nu_idx][-1] - nu_at_R;

	ds_dar[nu_idx][-1] = nu_at_R;

	// We don't need to calculate the surface term anymore,
	//  therefore we have 'i < ds_dar[r_idx].Size()-1'
	// -----------------------------------
	// [ May 25, 2022 ] We are calculating nu in the core region twice!
	// We can just use the values that are already calculated.
	// -----------------------------------
	// Core loop
	for (size_t i = 0;
		 i < std::min(ds_dar[r_idx].Size() - 1, ds_vis[r_idx].Size() - 1);
		 ++i)
	{
		ds_dar[nu_idx][i] = ds_vis[nu_idx][i];
	}

	for (size_t i = std::min(ds_dar[r_idx].Size() - 1, ds_vis[r_idx].Size() - 1);
		 i < ds_dar[r_idx].Size() - 1; ++i)
	{
		ds_dar[nu_idx][i] = ds_dar.Integrate(nu_der_idx, {ds_dar[r_idx][0], ds_dar[r_idx][i]}) - delta_nu_r;
	}

	ds_dar.Interpolate(0, nu_idx);
}

//--------------------------------------------------------------
// Metric function as a function of radius (in km)
double MixedStar::GetNu(const double &in_r) const
{
	// return gsl_spline_eval(nu_r_spline, in_r, mixed_r_accel) ;
	if (dark_core)
		return ds_vis.Evaluate(nu_idx, in_r);
	else
		return ds_dar.Evaluate(nu_idx, in_r);
}

//--------------------------------------------------------------
// Metric function as a function of radius (in km)
Zaki::Vector::DataColumn *MixedStar::GetNu()
{
	if (dark_core)
		return &ds_vis[nu_idx];
	else
		return &ds_dar[nu_idx];
}
//--------------------------------------------------------------
/// Returns the visible radius dataset
Zaki::Vector::DataColumn *MixedStar::GetRadius_Visible()
{
	return &ds_vis[r_idx];
}

//--------------------------------------------------------------
// Returns the dark radius dataset
Zaki::Vector::DataColumn *MixedStar::GetRadius_Dark()
{
	return &ds_dar[r_idx];
}

//--------------------------------------------------------------
// Returns the larger radius dataset
Zaki::Vector::DataColumn *MixedStar::GetRadius()
{
	if (dark_core)
		return &ds_vis[r_idx];
	else
		return &ds_dar[r_idx];
}

//--------------------------------------------------------------
// Visible mass (in km) as a function of radius
double MixedStar::GetMass_Visible(const double &in_r)
	const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Input radius is must be positive!");

		return -1;
	}

	if (in_r < ds_vis[r_idx][0])
		return 0;

	if (in_r > ds_vis[r_idx][-1])
		return ds_vis[m_idx][-1];

	// return gsl_spline_eval(mass_r_spline, in_r, visi_r_accel) ;
	return ds_vis.Evaluate(m_idx, in_r);
}
//--------------------------------------------------------------
// Visible mass (in km) as a function of radius
Zaki::Vector::DataColumn *
MixedStar::GetMass_Visible()
{
	return &ds_vis[m_idx];
}

//--------------------------------------------------------------
// Dark mass (in km) as a function of radius
double MixedStar::GetMass_Dark(const double &in_r)
	const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Input radius is must be positive!");

		return -1;
	}

	if (in_r < ds_dar[r_idx][0])
		return 0;

	if (in_r > ds_dar[r_idx][-1])
		return ds_dar[m_idx][-1];

	// return gsl_spline_eval(mass_r_d_spline, in_r, dark_r_accel) ;
	return ds_dar.Evaluate(m_idx, in_r);
}

//--------------------------------------------------------------
// Dark mass (in km) as a function of radius
Zaki::Vector::DataColumn *
MixedStar::GetMass_Dark()
{
	return &ds_dar[m_idx];
}

//--------------------------------------------------------------
// Dark + visible mass (in km) as a function of radius
double MixedStar::GetMass_Total(const double &in_r)
	const
{
	return GetMass_Dark(in_r) + GetMass_Visible(in_r);
}

//--------------------------------------------------------------
// Returns a data column holding total mass values
// extending from the origin to the surface of the star
Zaki::Vector::DataColumn *MixedStar::GetMass_Total()
{
	return &mass_tot_dc;
}

//--------------------------------------------------------------
// Visible baryon number density (fm^{-3}) as a function of radius
double MixedStar::GetRho_Visible(const double &in_r)
	const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Input radius is must be positive!");

		return -1;
	}

	if (in_r < ds_vis[r_idx][0] || in_r > ds_vis[r_idx][-1])
	{
		return 0;
	}

	// return gsl_spline_eval(rho_r_spline, in_r, visi_r_accel) ;
	return ds_vis.Evaluate(rho_idx, in_r);
}

//--------------------------------------------------------------
// Visible baryon number density (fm^{-3}) as a function of radius
Zaki::Vector::DataColumn *
MixedStar::GetRho_Visible()
{
	return &ds_vis[rho_idx];
}

//--------------------------------------------------------------
// Visible baryon number density (fm^{-3}) as a function of radius
// for a specific species labeled as (in_label)
Zaki::Vector::DataColumn *
MixedStar::GetRho_i_Visible(const std::string &in_label)
{

	for (auto &&c : rho_i_v_idx)
	{
		if (ds_vis[c].Label() == in_label)
		{
			return &ds_vis[c];
		}
	}

	Z_LOG_ERROR("Species density with label '" +
				in_label + "' was not found. Returning nullptr.");

	return nullptr;
}

//--------------------------------------------------------------
// Dark baryon number density (fm^{-3}) as a function of radius
double MixedStar::GetRho_Dark(const double &in_r)
	const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Input radius is must be positive!");

		return -1;
	}

	if (in_r < ds_dar[r_idx][0] || in_r > ds_dar[r_idx][-1])
	{
		return 0;
	}

	// return gsl_spline_eval(rho_r_d_spline, in_r, dark_r_accel) ;
	return ds_dar.Evaluate(rho_idx, in_r);
}

//--------------------------------------------------------------
// Dark baryon number density (fm^{-3}) as a function of radius
Zaki::Vector::DataColumn *
MixedStar::GetRho_Dark()
{
	return &ds_dar[rho_idx];
}

//--------------------------------------------------------------
// Energy density (in km^{-2}) as a function of radius
double MixedStar::GetEps_Visible(const double &in_r)
	const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Input radius is must be positive!");

		return -1;
	}

	if (in_r < ds_vis[r_idx][0] || in_r > ds_vis[r_idx][-1])
	{
		return 0;
	}

	// return gsl_spline_eval(eps_r_spline, in_r, visi_r_accel) ;
	return ds_vis.Evaluate(eps_idx, in_r);
}

//--------------------------------------------------------------
// Energy density (in km^{-2}) as a function of radius
Zaki::Vector::DataColumn *
MixedStar::GetEps_Visible()
{
	return &ds_vis[eps_idx];
}

//--------------------------------------------------------------
// Dark energy density (in km^{-2}) as a function of radius
double MixedStar::GetEps_Dark(const double &in_r)
	const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Input radius is must be positive!");

		return -1;
	}

	if (in_r < ds_dar[r_idx][0] || in_r > ds_dar[r_idx][-1])
	{
		return 0;
	}

	// return gsl_spline_eval(eps_r_d_spline, in_r, dark_r_accel) ;
	return ds_dar.Evaluate(eps_idx, in_r);
}

//--------------------------------------------------------------
// Dark energy density (in km^{-2}) as a function of radius
Zaki::Vector::DataColumn *
MixedStar::GetEps_Dark()
{
	return &ds_dar[eps_idx];
}

//--------------------------------------------------------------
// Pressure (in km^{-2}) as a function of radius
double MixedStar::GetPress_Visible(const double &in_r)
	const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Input radius is must be positive!");

		return -1;
	}

	if (in_r < ds_vis[r_idx][0] || in_r > ds_vis[r_idx][-1])
	{
		return 0;
	}

	// return gsl_spline_eval(press_r_spline, in_r, visi_r_accel) ;
	return ds_vis.Evaluate(pre_idx, in_r);
}

//--------------------------------------------------------------
// Pressure (in km^{-2}) as a function of radius
Zaki::Vector::DataColumn *
MixedStar::GetPress_Visible()
{
	return &ds_vis[pre_idx];
}

//--------------------------------------------------------------
// Dark pressure (in km^{-2}) as a function of radius
double MixedStar::GetPress_Dark(const double &in_r)
	const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Input radius is must be positive!");

		return -1;
	}

	if (in_r < ds_dar[r_idx][0] || in_r > ds_dar[r_idx][-1])
	{
		return 0;
	}

	// return gsl_spline_eval(press_r_d_spline, in_r, dark_r_accel) ;
	return ds_dar.Evaluate(pre_idx, in_r);
}

//--------------------------------------------------------------
// Dark pressure (in km^{-2}) as a function of radius
Zaki::Vector::DataColumn *
MixedStar::GetPress_Dark()
{
	return &ds_dar[pre_idx];
}

//--------------------------------------------------------------
// Returns 'sequence'
MixedSeqPoint MixedStar::GetSequence() const
{
	return sequence;
}

//--------------------------------------------------------------
// Visible baryon number inegrand
double MixedStar::BaryonNumIntegrand(double in_r)
{
	double out = in_r * in_r;
	out *= 4 * M_PI * GetRho_Visible(in_r);

	// // converting fm^{-3} to km^{-3}
	// out *= pow(10, 54) ;

	out /= sqrt(1 - 2 * GetMass_Total(in_r) / in_r);

	return out;
}

//--------------------------------------------------------------
// Dark baryon number inegrand
double MixedStar::BaryonNumIntegrand_Dark(double in_r)
{
	double out = in_r * in_r;
	out *= 4 * M_PI * GetRho_Dark(in_r);

	// // converting fm^{-3} to km^{-3}
	// out *= pow(10, 54) ;

	out /= sqrt(1 - 2 * GetMass_Total(in_r) / in_r);

	return out;
}

//--------------------------------------------------------------
// Visible baryon number as a function of radius
// [ Deactivated on May 25, 2022 ]
double MixedStar::Find_BaryonNum_Visible(const double &in_r)
{
	PROFILE_FUNCTION();
	// double result, err ;

	// gsl_integration_workspace *w =
	//  gsl_integration_workspace_alloc(1500);

	// Zaki::Math::GSLFuncWrapper<MixedStar, double (MixedStar::*)(double)>
	//     Fp(this, &MixedStar::BaryonNumIntegrand);

	// gsl_function F = static_cast<gsl_function> (Fp) ;

	// gsl_integration_qag(&F_B_Visi, ds_vis[r_idx][0], in_r, 1e-9, 1e-9,
	//                     1500, 1, integ_wrk_space, &result, &err) ;

	// // Memory leak fixed on May 6, 2022
	// gsl_integration_workspace_free(w) ;

	// converting fm^{-3} to km^{-3}
	// result *= pow(10, 54) ;

	// return result ;
	return 0;
}

//--------------------------------------------------------------
// Visible baryon number as a function of radius
// [ Deactivated on May 25, 2022 ]
// double MixedStar::Find_BaryonNum_Dark(const double& in_r)
// {
//   PROFILE_FUNCTION() ;
// double result, err ;

// gsl_integration_workspace *w =
//  gsl_integration_workspace_alloc(1500);

// Zaki::Math::GSLFuncWrapper<MixedStar, double (MixedStar::*)(double)>
//     FW_B_Dark(this, &MixedStar::BaryonNumIntegrand_Dark);

// gsl_function F_B_Dark = static_cast<gsl_function> (FW_B_Dark) ;

// gsl_integration_qag(&F_B_Dark, ds_dar[r_idx][0], in_r, 1e-9, 1e-9,
//                     1500, 1, integ_wrk_space, &result, &err);

// Memory leak fixed on May 6, 2022
// gsl_integration_workspace_free(w) ;

// converting fm^{-3} to km^{-3}
// result *= pow(10, 54) ;

// return result ;
// return 0 ;
// }

//--------------------------------------------------------------
// Visible baryon number
// double MixedStar::Find_BaryonNum_Visible()
// {
// return Find_BaryonNum_Visible(ds_vis[r_idx][-1]) ;
// }

//--------------------------------------------------------------
// Dark baryon number
// double MixedStar::Find_BaryonNum_Dark()
// {
// return Find_BaryonNum_Dark(ds_dar[r_idx][-1]) ;
// }

//--------------------------------------------------------------
//// Moment of inertia (in km^3) inegrand
// This is wrong and missing a factor of
//  [ 1 - \omega(r) / \Omega ]
// double MixedStar::MomInertiaIntegrand(double in_r)
// {
//   double out = pow(in_r, 4) ;
//         out *= 8*M_PI/3. ;
//         out *= ( GetEps(in_r) + GetPress(in_r) ) ;

//         // // converting fm^{-3} to km^{-3}
//         // out *= pow(10, 54) ;

//         out *= ::exp(-GetNu(in_r)) ;
//         out /= sqrt(1 - 2*GetMass(in_r) / in_r) ;

//   return out ;
// }

//--------------------------------------------------------------
// Total moment of inertia (in km^3)
// MomInertiaIntegrand is wrong!
double MixedStar::Find_MomInertia()
{
	PROFILE_FUNCTION();

	// double result = 0 ;

	// rot_solver.SetWrkDir(wrk_dir) ;

	// rot_solver.Solve_Dark(5e-3) ;

	// OmegaSeqPoint omega_seq_pt = rot_solver.GetOmegaSeq()[0] ;
	// result = omega_seq_pt.J / ( omega_seq_pt.Omega / Zaki::Physics::LIGHT_C_KM_S ) ;

	rot_solver.FindMixedMomInertia();

	// std::cout << "\t Mom. I = " <<  MomI << "\n" ;
	// if ( gsl_isnan(omega_seq_pt.J) )
	// {
	//   Z_LOG_ERROR("The moment of inertia is Nan!") ;
	//   exit(EXIT_FAILURE) ;
	// }

	return MomI;
}

//--------------------------------------------------------------
// Exports the mixed star profile
void MixedStar::Export(
	const Zaki::String::Directory &in_dir)
{
	char seq_header[200];
	snprintf(seq_header, sizeof(seq_header),
			 "    %-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s",
			 "ec(g/cm^3)", "M", "R(km)", "pc(dyne/cm^2)", "B", "I(km^3)");

	std::time_t end_time = std::chrono::system_clock::to_time_t(
		std::chrono::system_clock::now());

	ds_vis.AddHead("# --------------------------------------------------------"
				   "--------------------------------------------------------\n");

	ds_vis.AddHead("# Visible profile generated on ");
	ds_vis.AddHead(std::string(std::ctime(&end_time)));
	ds_vis.AddHead("# --------------------------------------------------------"
				   "--------------------------------------------------------\n");
	ds_vis.AddHead("# Sequence point info: \n");
	ds_vis.AddHead("#         " + std::string(seq_header) + "\n");
	ds_vis.AddHead("# Visible (" + std::to_string(sequence.v_idx) + ") " + sequence.v.Str() + "\n");
	ds_vis.AddHead("# Dark    (" + std::to_string(sequence.d_idx) + ") " + sequence.d.Str() + "\n");
	ds_vis.AddHead("# --------------------------------------------------------"
				   "--------------------------------------------------------\n");
	ds_vis.AddFoot("# --------------------------------------------------------"
				   "--------------------------------------------------------\n");

	ds_dar.AddHead("# --------------------------------------------------------"
				   "--------------------------------------------------------\n");
	ds_dar.AddHead("# Dark profile generated on ");
	ds_dar.AddHead(std::string(std::ctime(&end_time)));
	ds_dar.AddHead("# --------------------------------------------------------"
				   "--------------------------------------------------------\n");
	ds_dar.AddHead("# Sequence point info: \n");
	ds_dar.AddHead("#         " + std::string(seq_header) + "\n");
	ds_dar.AddHead("# Visible (" + std::to_string(sequence.v_idx) + ") " + sequence.v.Str() + "\n");
	ds_dar.AddHead("# Dark    (" + std::to_string(sequence.d_idx) + ") " + sequence.d.Str() + "\n");
	ds_dar.AddHead("# --------------------------------------------------------"
				   "--------------------------------------------------------\n");
	ds_dar.AddFoot("# --------------------------------------------------------"
				   "--------------------------------------------------------\n");

	ds_vis.Export(in_dir.ThisFileDir() + ("/V_" + in_dir.ThisFile().Str()));
	ds_dar.Export(in_dir.ThisFileDir() + ("/D_" + in_dir.ThisFile().Str()));

	ds_vis.ClearHeadFoot();
	ds_dar.ClearHeadFoot();
}
//--------------------------------------------------------------
//==============================================================
