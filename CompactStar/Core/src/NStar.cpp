/*
  NStar class
*/

// #include <gsl/gsl_math.h>

#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <Zaki/Util/Instrumentor.hpp>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

#include <gsl/gsl_integration.h>
#include <limits>

using namespace CompactStar::Core;
//==============================================================
//               Conversion Factors
//==============================================================
static constexpr double FM3_TO_KM3 = 1e54;
//==============================================================
//                        NStar class
//==============================================================
// Default Constructor
// Be sure to run SurfaceIsReached if using this constructor!
NStar::NStar() : Prog("NStar"), B_integrand{}
{
	rot_solver.AttachNStar(this);
}

//--------------------------------------------------------------
// Constructor from TOV Solutions
NStar::NStar(const std::vector<TOVPoint> &in_tov)
	: Prog("NStar"), B_integrand{}
{
	rot_solver.AttachNStar(this);
	BuildFromTOV(in_tov, /*species_labels=*/nullptr);
}

//--------------------------------------------------------------
// Constructor from TOV Solutions with species labels
NStar::NStar(const std::vector<TOVPoint> &in_tov,
			 const std::vector<std::string> &species_labels)
	: Prog("NStar"), B_integrand{}
{
	rot_solver.AttachNStar(this);
	BuildFromTOV(in_tov, &species_labels);
}

//--------------------------------------------------------------
// Builder — Populate profile from TOV data
// ---------------------------------------------------------
void NStar::BuildFromTOV(const std::vector<TOVPoint> &in_tov,
						 const std::vector<std::string> *species_labels)
{
	PROFILE_FUNCTION();

	// -------- guardrails & basic sizes --------
	if (in_tov.empty())
	{
		Z_LOG_ERROR("Empty TOV vector; leaving object uninitialized.");
		surface_ready = false;
		return;
	}

	const std::size_t n_rows = in_tov.size();

	// Be robust: infer species count as the maximum rho_i length across rows.
	std::size_t n_species = 0;
	for (const auto &tp : in_tov)
		n_species = std::max(n_species, tp.rho_i.size());

	// If labels provided, require consistency (warn if mismatch).
	if (species_labels && species_labels->size() != n_species)
	{
		Z_LOG_ERROR("Species_labels.size() = " +
					std::to_string(species_labels->size()) +
					" != inferred n_species = " + std::to_string(n_species) +
					". Proceeding with inferred n_species; extra/missing labels ignored.");
	}
	// ------------------------------------------------------------
	// 0) Fresh start
	// ------------------------------------------------------------
	Reset(); // clears ds, B_integrand, sequence, rho_i_idx, prof_

	// ============================================================
	// 1) PROFILE-FIRST BUILD
	// ============================================================
	{
		auto &radial = prof_.radial;
		radial.ClearRows();
		radial.Reserve(8 + n_species, n_rows);
		B_integrand.Reserve(2, n_rows);

		// ----------------------------------
		// REMOVE THIS BLOCK
		// ----------------------------------
		// create canonical columns in the exact order we settled on:
		// 0 r, 1 m, 2 nu', 3 p, 4 eps, 5 rho, 6 nu, 7 lambda
		// radial.AddColumn("r(km)");
		// prof_.idx_r = 0;

		// radial.AddColumn("m(km)");
		// prof_.idx_m = 1;

		// radial.AddColumn("nu_prime(km^-1)");
		// prof_.idx_nuprime = 2;

		// radial.AddColumn("p(km^-2)");
		// prof_.idx_p = 3;

		// radial.AddColumn("eps(km^-2)");
		// prof_.idx_eps = 4;

		// radial.AddColumn("nB(fm^-3)");
		// prof_.idx_nb = 5;

		// radial.AddColumn("nu");
		// prof_.idx_nu = 6;

		// radial.AddColumn("lambda");
		// prof_.idx_lambda = 7;
		// ----------------------------------

		// label the ones we just reserved
		radial[0].SetLabel("r(km)");
		prof_.idx_r = 0;
		radial[1].SetLabel("m(km)");
		prof_.idx_m = 1;
		radial[2].SetLabel("nu_prime(km^-1)");
		prof_.idx_nuprime = 2;
		radial[3].SetLabel("p(km^-2)");
		prof_.idx_p = 3;
		radial[4].SetLabel("eps(km^-2)");
		prof_.idx_eps = 4;
		radial[5].SetLabel("nB(fm^-3)");
		prof_.idx_nb = 5;
		radial[6].SetLabel("nu");
		prof_.idx_nu = 6;
		radial[7].SetLabel("lambda");
		prof_.idx_lambda = 7;

		// species (after the fixed ones)
		// species columns: indices 8..(8+n_species-1)
		prof_.species_labels.clear();
		prof_.species_idx.clear();
		prof_.species_labels.reserve(n_species);
		prof_.species_idx.reserve(n_species);

		for (std::size_t j = 0; j < n_species; ++j)
		{
			std::string lbl = (species_labels && j < species_labels->size())
								  ? (*species_labels)[j]
								  : "rho_i_" + std::to_string(j);

			// const int col_idx = static_cast<int>(radial.Dim().size()); // REMOVE THIS LINE
			const int col_idx = 8 + static_cast<int>(j);
			radial[col_idx].SetLabel(lbl);
			// radial.AddColumn(lbl); // REMOVE THIS LINE
			prof_.AddSpecies(lbl, col_idx);
		}

		// now fill rows
		for (const auto &tp : in_tov)
		{
			// r (km)
			double r_km = tp.r;
			radial[prof_.idx_r].PushBack(r_km);

			double m_km = Zaki::Physics::SUN_M_KM * tp.m;
			// m (km): tp.m is in solar masses
			radial[prof_.idx_m].PushBack(m_km);

			// nu' (1/cm -> 1/km)
			radial[prof_.idx_nuprime].PushBack(tp.nu_der * 1e5);

			// p (→ km^-2)
			radial[prof_.idx_p].PushBack(
				tp.p * Zaki::Physics::INV_FM4_2_INV_KM2 /
				Zaki::Physics::INV_FM4_2_Dyn_CM2);

			// eps (→ km^-2)
			radial[prof_.idx_eps].PushBack(
				tp.e * Zaki::Physics::INV_FM4_2_INV_KM2 /
				Zaki::Physics::INV_FM4_2_G_CM3);

			// nB (fm^-3)
			radial[prof_.idx_nb].PushBack(tp.rho);

			// nu (will be built later)
			radial[prof_.idx_nu].PushBack(0.0);

			// ------------------------------------------------------------
			// Compute and append λ to the profile
			// ------------------------------------------------------------
			// Compute the Schwarzschild-like factor:
			// const double r_km = tp.r;
			// const double m_km = Zaki::Physics::SUN_M_KM * tp.m;

			double denom = 1.0;
			if (r_km > 0.0)
			{
				denom = 1.0 - 2.0 * m_km / r_km;
				if (denom <= 0.0)
					denom = 1e-15; // protect against log of non-positive
			}

			// g_rr = 1 / (1 - 2M/r)
			// const double grr = 1.0 / denom;

			// λ = 0.5 * ln(g_rr) = -0.5 * ln(denom)
			// const double lambda_geom = 0.5 * std::log(grr);
			// or equivalently:
			const double lambda_geom = -0.5 * std::log(denom);

			radial[prof_.idx_lambda].PushBack(lambda_geom);
			// ------------------------------------------------------------

			// ------------------------------
			// species values (pad with 0.0)
			if (!tp.rho_i.empty())
			{
				for (std::size_t j = 0; j < n_species; ++j)
				{
					const double val = (j < tp.rho_i.size()) ? tp.rho_i[j] : 0.0;
					const int col_idx = prof_.species_idx[j];
					radial[col_idx].PushBack(val);
				}
			}
			else
			{
				// no species in this row → append zeros to all species cols
				for (std::size_t j = 0; j < n_species; ++j)
				{
					const int col_idx = prof_.species_idx[j];
					radial[col_idx].PushBack(0.0);
				}
			}
		}

		// interpolate the columns we actually have
		{
			std::vector<int> interp_cols;
			interp_cols.push_back(prof_.idx_m);
			interp_cols.push_back(prof_.idx_nuprime);
			interp_cols.push_back(prof_.idx_nb);
			interp_cols.push_back(prof_.idx_eps);
			interp_cols.push_back(prof_.idx_p);
			radial.Interpolate(prof_.idx_r, interp_cols);
		}

		// build ν(r) from ν'(r)
		EvaluateNu(); // profile-aware version

		// build baryon number integrand from profile
		{
			B_integrand[0] = radial[prof_.idx_r];
			B_integrand[0].SetLabel("r(km)");

			B_integrand[1] = radial[prof_.idx_r].pow(2);
			B_integrand[1].SetLabel("B_v");
			B_integrand[1] *= 4.0 * M_PI * radial[prof_.idx_nb];
			B_integrand[1] /= (1.0 - 2.0 * radial[prof_.idx_m] / radial[prof_.idx_r]).sqrt();
			B_integrand[1] *= FM3_TO_KM3;
			B_integrand.Interpolate(0, 1);
		}

		// fill profile's sequence
		prof_.seq_point.clear();
		// central energy density
		prof_.seq_point.ec =
			radial[prof_.idx_eps][0] *
			Zaki::Physics::INV_FM4_2_G_CM3 /
			Zaki::Physics::INV_FM4_2_INV_KM2; // g/cm^3

		// surface radius
		prof_.seq_point.r = radial[prof_.idx_r][-1]; // km
		prof_.R = prof_.seq_point.r;

		// surface mass (in both solar masses and km)
		const double Msurf_km = radial[prof_.idx_m][-1];			  // km
		const double Msurf_Msun = Msurf_km / Zaki::Physics::SUN_M_KM; // Msun

		prof_.M = Msurf_km;				// km
		prof_.seq_point.m = Msurf_Msun; // Msun

		// central pressure
		prof_.seq_point.pc =
			radial[prof_.idx_p][0] *
			Zaki::Physics::INV_FM4_2_Dyn_CM2 /
			Zaki::Physics::INV_FM4_2_INV_KM2; // dyne/cm^2

		// baryon number
		prof_.seq_point.b =
			B_integrand.Integrate(1, {radial[prof_.idx_r][0], radial[prof_.idx_r][-1]});

		// moment of inertia
		prof_.seq_point.I = Find_MomInertia();

		// surface redshift factor (if we have ν)
		if (prof_.HasColumn(StarProfile::Column::MetricNu))
		{
			const auto &nu_col = radial[prof_.idx_nu];
			if (nu_col.Size() > 0)
				prof_.z_surf = std::exp(nu_col[-1]);
			else
				prof_.z_surf = 0.0;
		}

		surface_ready = true;
	}
}

/* -------------------- LEGACY VERSION FOR REFERENCE --------------------

// ---------------------------------------------------------
// Builder — Populate ds/B_integrand/sequence from TOV data
// ---------------------------------------------------------
void NStar::BuildFromTOV(const std::vector<TOVPoint> &in_tov,
									  const std::vector<std::string> *species_labels)
{
	PROFILE_FUNCTION();

	// -------- guardrails & basic sizes --------
	if (in_tov.empty())
	{
		Z_LOG_ERROR("Empty TOV vector; leaving object uninitialized.");
		surface_ready = false;
		return;
	}

	const size_t n_rows = in_tov.size();

	// Be robust: infer species count as the maximum rho_i length across rows.
	size_t n_species = 0;
	for (const auto &tp : in_tov)
		n_species = std::max(n_species, tp.rho_i.size());

	// If labels provided, require consistency (warn if mismatch).
	if (species_labels && species_labels->size() != n_species)
	{
		Z_LOG_ERROR("Species_labels.size() = " +
					std::to_string(species_labels->size()) +
					" != inferred n_species = " + std::to_string(n_species) +
					". Proceeding with inferred n_species; extra/missing labels ignored.");
	}

	// -------- allocate datasets --------
	// 7 fixed columns: r, m, nu', p, eps, rho, nu  (+ per-species rho_i)
	ds.Reserve(7 + n_species, n_rows);
	B_integrand.Reserve(2, n_rows);

	// -------- label fixed columns --------
	col(Col::R).label      = "r(km)";
	col(Col::M).label      = "m(km)";
	col(Col::NuPrime).label= "nu'(km^-1)";
	col(Col::P).label      = "p(km^-2)";
	col(Col::Eps).label    = "eps(km^-2)";
	col(Col::Rho).label    = "rho(fm^-3)";
	col(Col::Nu).label     = "nu";

	// -------- set up per-species columns & labels --------
	rho_i_idx.clear();
	rho_i_idx.reserve(n_species);
	for (size_t j = 0; j < n_species; ++j)
	{
		const int idx_col = 6 + static_cast<int>(j);
		rho_i_idx.emplace_back(idx_col);
		if (species_labels && j < species_labels->size())
			ds[idx_col].label = (*species_labels)[j];
		else
			ds[idx_col].label = "X" + std::to_string(idx_col); // fallback
	}

	// -------- fill rows --------
	for (const auto &tp : in_tov)
	{
		col(Col::R).PushBack(tp.r);
		col(Col::M).PushBack(Zaki::Physics::SUN_M_KM * tp.m);
		col(Col::Rho).PushBack(tp.rho);
		col(Col::Eps).PushBack(
			tp.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3);
		col(Col::P).PushBack(
			tp.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2);
		col(Col::NuPrime).PushBack(tp.nu_der * 1e5);

		for (size_t j = 0; j < n_species; ++j)
		{
			const double val = (j < tp.rho_i.size()) ? tp.rho_i[j] : 0.0;
			ds[rho_i_idx[j]].PushBack(val);
		}
	}

	// -------- make splines --------
	ds.Interpolate(idx(Col::R),
				   { idx(Col::M), idx(Col::NuPrime), idx(Col::Rho),
					 idx(Col::Eps), idx(Col::P) });

	// -------- build ν(r) --------
	EvaluateNu();

	// -------- Baryon-number integrand dataset --------
	B_integrand[0] = col(Col::R);
	B_integrand[0].label = "r(km)";
	B_integrand[1] = col(Col::R).pow(2);
	B_integrand[1].label = "B_v";
	B_integrand[1] *= 4.0 * M_PI * col(Col::Rho);
	B_integrand[1] /= (1.0 - 2.0 * col(Col::M) / col(Col::R)).sqrt();
	B_integrand[1] *= FM3_TO_KM3;
	B_integrand.Interpolate(0, 1);

	// -------- sequence point --------
	surface_ready = true;

	sequence.ec = col(Col::Eps)[0] * Zaki::Physics::INV_FM4_2_G_CM3 / Zaki::Physics::INV_FM4_2_INV_KM2;
	sequence.m  = col(Col::M)[-1] / Zaki::Physics::SUN_M_KM;
	sequence.r  = col(Col::R)[-1];
	sequence.pc = col(Col::P)[0] * Zaki::Physics::INV_FM4_2_Dyn_CM2 / Zaki::Physics::INV_FM4_2_INV_KM2;
	sequence.b  = B_integrand.Integrate(1, { col(Col::R)[0], col(Col::R)[-1] });
	sequence.I  = Find_MomInertia();
}

------------------------------------------------------------------------- */
// //--------------------------------------------------------------
// // Constructor from TOV Solutions
// NStar::NStar(const std::vector<CompactStar::TOVPoint> &in_tov)
// 	: Prog("NStar", true),
// 	  ds(7 + in_tov[0].rho_i.size(), in_tov.size()),
// 	  B_integrand(2, in_tov.size())
// {

// 	B_integrand[0].label = "r(km)";
// 	B_integrand[1].label = "B_v";

// 	rot_solver.AttachNStar(this);

// 	col(Col::M).label = "m(km)";
// 	col(Col::R).label = "r(km)";
// 	col(Col::Rho).label = "rho(fm^-3)";
// 	col(Col::Eps).label = "eps(km^-2)";
// 	col(Col::P).label = "p(km^-2)";
// 	col(Col::NuPrime).label = "nu'(km^-1)";
// 	col(Col::Nu).label = "nu";

// 	// ds[m_idx].label = "m";
// 	// ds[r_idx].label = "r";
// 	// ds[rho_idx].label = "rho";
// 	// ds[eps_idx].label = "eps";
// 	// ds[pre_idx].label = "p";
// 	// ds[nu_der_idx].label = "nu'";
// 	// ds[nu_idx].label = "nu";

// 	for (size_t j = 0; j < in_tov[0].rho_i.size(); j++)
// 	{
// 		rho_i_idx.emplace_back(6 + j);
// 		ds[rho_i_idx[j]].label = "X" + std::to_string(6 + j); // Fix?!!!
// 	}

// 	// ...............................................
// 	for (auto &&i : in_tov)
// 	{
// 		// ds[r_idx].PushBack(i.r);																		  // in km
// 		// ds[m_idx].PushBack(Zaki::Physics::SUN_M_KM * i.m);												  // in km
// 		// ds[rho_idx].PushBack(i.rho);																	  // in fm^{-3}
// 		// ds[eps_idx].PushBack(i.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3);	  // in km^{-2}
// 		// ds[pre_idx].PushBack(i.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2); // in km^{-2}
// 		// ds[nu_der_idx].PushBack(i.nu_der * 1e+5);
// 		col(Col::R).PushBack(i.r);																		  // in km
// 		col(Col::M).PushBack(Zaki::Physics::SUN_M_KM * i.m);											  // in km
// 		col(Col::Rho).PushBack(i.rho);																	  // in fm^{-3}
// 		col(Col::Eps).PushBack(i.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3); // in km^{-2}
// 		col(Col::P).PushBack(i.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2); // in km^{-2}
// 		col(Col::NuPrime).PushBack(i.nu_der * 1e+5);													  // convert 1/cm to 1/km

// 		for (size_t j = 0; j < i.rho_i.size(); j++)
// 		{
// 			ds[rho_i_idx[j]].PushBack(i.rho_i[j]);
// 		}
// 	}
// 	// ...............................................

// 	ds.Interpolate(idx(Col::R), {idx(Col::M), idx(Col::NuPrime),
// 								 idx(Col::Rho), idx(Col::Eps), idx(Col::P)});

// 	EvaluateNu();

// 	// ..................................................................
// 	//            B Integrand
// 	// ..................................................................
// 	B_integrand[0] = col(Col::R);
// 	B_integrand[1] = col(Col::R).pow(2);
// 	B_integrand[1] *= 4 * M_PI * col(Col::Rho);

// 	B_integrand[1] /= (1. - 2. * col(Col::M) / col(Col::R)).sqrt();

// 	// converting fm^{-3} to km^{-3}
// 	// B_integrand[1] *= pow(10, 54);
// 	B_integrand[1] *= FM3_TO_KM3;

// 	B_integrand.Interpolate(0, 1);
// 	// ..................................................................

// 	sequence.ec = col(Col::Eps)[0] * Zaki::Physics::INV_FM4_2_G_CM3 /
// 				  Zaki::Physics::INV_FM4_2_INV_KM2;
// 	sequence.m = col(Col::M)[-1] / Zaki::Physics::SUN_M_KM;
// 	sequence.r = col(Col::R)[-1];
// 	sequence.pc = col(Col::P)[0] * Zaki::Physics::INV_FM4_2_Dyn_CM2 /
// 				  Zaki::Physics::INV_FM4_2_INV_KM2;
// 	// sequence.b  = Find_BaryonNum_Visible() ;
// 	sequence.b = B_integrand.Integrate(1, {col(Col::R)[0], col(Col::R)[-1]});

// 	sequence.I = Find_MomInertia();
// 	// ..................................................................
// }

//--------------------------------------------------------------
//--------------------------------------------------------------
// Initializes the dataset from a TOV solver
void NStar::InitFromTOVSolver(const TOVSolver *in_tov_solver)
{
	if (!in_tov_solver)
	{
		Z_LOG_ERROR("Null TOVSolver pointer.");
		return;
	}

	// ------------------------------------------------------------
	// 1) Fresh start
	// ------------------------------------------------------------
	Reset(); // clears ds, B_integrand, sequence, rho_i_idx, prof_

	// ------------------------------------------------------------
	// 2) Read solver metadata
	// ------------------------------------------------------------
	const std::size_t n_species = in_tov_solver->eos_tab.rho_i.size(); // number of species columns
	const std::size_t n_rows_expect = in_tov_solver->radial_res;	   // expected radial samples
	B_integrand.Reserve(2, n_rows_expect);
	// ------------------------------------------------------------
	// 3) ----- PROFILE-FIRST INITIALIZATION ----------------------
	// ------------------------------------------------------------
	{
		auto &radial = prof_.radial;
		radial.ClearRows();
		radial.Reserve(8 + n_species, n_rows_expect);

		// label the ones we just reserved
		radial[0].SetLabel("r(km)");
		prof_.idx_r = 0;
		radial[1].SetLabel("m(km)");
		prof_.idx_m = 1;
		radial[2].SetLabel("nu_prime(km^-1)");
		prof_.idx_nuprime = 2;
		radial[3].SetLabel("p(km^-2)");
		prof_.idx_p = 3;
		radial[4].SetLabel("eps(km^-2)");
		prof_.idx_eps = 4;
		radial[5].SetLabel("nB(fm^-3)");
		prof_.idx_nb = 5;
		radial[6].SetLabel("nu");
		prof_.idx_nu = 6;
		radial[7].SetLabel("lambda");
		prof_.idx_lambda = 7;

		// Species columns start right after the 8 fixed ones
		prof_.species_labels.clear();
		prof_.species_idx.clear();
		prof_.species_labels.reserve(n_species);
		prof_.species_idx.reserve(n_species);

		for (std::size_t j = 0; j < n_species; ++j)
		{
			const std::string lbl =
				(j < in_tov_solver->eos_tab.extra_labels.size())
					? in_tov_solver->eos_tab.extra_labels[j]
					: ("rho_i_" + std::to_string(j));
			const int col_idx = 8 + static_cast<int>(j);
			radial[col_idx].SetLabel(lbl);
			prof_.AddSpecies(lbl, col_idx);
		}
		// fill profile’s sequence with zeros; it will be filled in FinalizeSurface()
		prof_.seq_point.clear();

		// also set surface M, R, z_surf to 0 here
		prof_.M = 0.0; // in solar masses
		prof_.R = 0.0;
		prof_.z_surf = 0.0;
	}

	// ------------------------------------------------------------
	// 4) ----- LEGACY DS INITIALIZATION -----------
	// NOTE: profile uses 8 fixed cols; legacy used 7 → species index shifts by 1
	// ------------------------------------------------------------
	// {
	// 	// Allocate columns and reserve row capacity
	// 	ds.Reserve(7 + n_species, n_rows_expect);
	// 	B_integrand.Reserve(2, n_rows_expect);

	// 	// Label fixed columns
	// 	col(Col::M).label = "m(km)";
	// 	col(Col::R).label = "r(km)";
	// 	col(Col::Rho).label = "rho(fm^-3)";
	// 	col(Col::Eps).label = "eps(km^-2)";
	// 	col(Col::P).label = "p(km^-2)";
	// 	col(Col::NuPrime).label = "nu'(km^-1)";
	// 	col(Col::Nu).label = "nu";

	// NOTE: profile uses 8 fixed cols; legacy used 7 → species index shifts by 1
	// 	// per-species columns in legacy ds
	// 	rho_i_idx.clear();
	// 	rho_i_idx.reserve(n_species);
	// 	for (std::size_t j = 0; j < n_species; ++j)
	// 	{
	// 		const int idx_col = 7 + static_cast<int>(j); // NOTE: legacy had 7 fixed cols
	// 		rho_i_idx.emplace_back(idx_col);

	// 		if (j < in_tov_solver->eos_tab.extra_labels.size())
	// 			ds[idx_col].label = in_tov_solver->eos_tab.extra_labels[j];
	// 		else
	// 			ds[idx_col].label = "X" + std::to_string(idx_col);
	// 	}

	// 	// Simple labels for the integrand dataset
	// 	B_integrand[0].label = "r(km)";
	// 	B_integrand[1].label = "B";
	// }
}

/* -------------- LEGACY VERSION FOR REFERENCE -----------------

//--------------------------------------------------------------
// Initializes the dataset
void NStar::InitFromTOVSolver(const TOVSolver *in_tov_solver)
{
	if (!in_tov_solver)
	{
		Z_LOG_ERROR("null TOVSolver pointer.");
		return;
	}

	// Fresh start
	Reset(); // clears ds rows, B_integrand rows, sequence, flags, rho_i_idx

	// Infer species count and planned row count from solver metadata
	const size_t n_species     = in_tov_solver->eos_tab.rho_i.size(); // number of species columns
	const size_t n_rows_reserve = in_tov_solver->radial_res;          // expected radial samples

	// Allocate columns and reserve row capacity
	ds.Reserve(7 + n_species, n_rows_reserve);
	B_integrand.Reserve(2, n_rows_reserve);

	// Label fixed columns
	col(Col::M).label       = "m(km)";
	col(Col::R).label       = "r(km)";
	col(Col::Rho).label     = "rho(fm^-3)";
	col(Col::Eps).label     = "eps(km^-2)";
	col(Col::P).label       = "p(km^-2)";
	col(Col::NuPrime).label = "nu'(km^-1)";
	col(Col::Nu).label      = "nu";

	// Set up per-species columns & labels
	rho_i_idx.clear();
	rho_i_idx.reserve(n_species);
	for (size_t j = 0; j < n_species; ++j)
	{
		const int idx_col = 6 + static_cast<int>(j);
		rho_i_idx.emplace_back(idx_col);

		// Use EOS extra_labels if available, otherwise fallback
		if (j < in_tov_solver->eos_tab.extra_labels.size())
			ds[idx_col].label = in_tov_solver->eos_tab.extra_labels[j];
		else
			ds[idx_col].label = "X" + std::to_string(idx_col);
	}

	// Simple labels for the integrand dataset
	B_integrand[0].label = "r(km)";
	B_integrand[1].label = "B";
}

---------------------------------------------------------------- */

// ------------------------------------------------------------
// Init interpolants from profile
// ------------------------------------------------------------
void NStar::InitInterpolantsFromProfile_()
{
	// nothing to do if profile is empty
	if (prof_.empty())
		return;

	// We only interpolate columns that actually exist.
	// NOTE: our StarProfile throws in Get(...) if column is missing,
	// so we MUST guard with HasColumn(...) before calling Interpolate.

	// Radius column must exist for interpolation to make sense
	if (!prof_.HasColumn(StarProfile::Column::Radius))
	{
		Z_LOG_ERROR("StarProfile has no radius column; cannot finalize.");
		surface_ready = false;
		return;
	}

	// auto &radial = prof_.radial;
	const int rcol = prof_.GetColumnIndex(StarProfile::Column::Radius);

	if (prof_.radial[rcol].Size() == 0)
	{
		Z_LOG_ERROR("StarProfile radius column is empty; nothing to finalize.");
		surface_ready = false;
		return;
	}

	if (prof_.HasColumn(StarProfile::Column::Mass))
		prof_.radial.Interpolate(rcol, prof_.GetColumnIndex(StarProfile::Column::Mass));

	if (prof_.HasColumn(StarProfile::Column::Pressure))
		prof_.radial.Interpolate(rcol, prof_.GetColumnIndex(StarProfile::Column::Pressure));

	if (prof_.HasColumn(StarProfile::Column::EnergyDensity))
		prof_.radial.Interpolate(rcol, prof_.GetColumnIndex(StarProfile::Column::EnergyDensity));

	if (prof_.HasColumn(StarProfile::Column::BaryonDensity))
		prof_.radial.Interpolate(rcol, prof_.GetColumnIndex(StarProfile::Column::BaryonDensity));

	if (prof_.HasColumn(StarProfile::Column::MetricNuPrime))
	{
		const int nupcol = prof_.GetColumnIndex(StarProfile::Column::MetricNuPrime);
		if (prof_.radial[nupcol].Size() == prof_.radial[rcol].Size())
			prof_.radial.Interpolate(rcol, nupcol);
		else
			Z_LOG_ERROR("MetricNuPrime column size mismatch; skipping ν′ interpolation.");
	}

	// if (prof_.HasColumn(StarProfile::Column::MetricLambda))
	// 	prof_.radial.Interpolate(rcol, prof_.GetColumnIndex(StarProfile::Column::MetricLambda));

	if (prof_.HasColumn(StarProfile::Column::MetricLambda))
	{
		const int lcol = prof_.GetColumnIndex(StarProfile::Column::MetricLambda);
		if (prof_.radial[lcol].Size() == prof_.radial[rcol].Size())
			prof_.radial.Interpolate(rcol, lcol);
		else
			Z_LOG_ERROR("lambda column size mismatch; skipping lambda interpolation.");
	}

	// per-species
	for (std::size_t i = 0; i < prof_.species_idx.size(); ++i)
	{
		const int scol = prof_.species_idx[i];
		if (prof_.IsValidColumnIndex(scol))
			prof_.radial.Interpolate(rcol, scol);
	}
}

//--------------------------------------------------------------
// Print profile column sizes
void NStar::PrintProfileColumnSizes() const
{
	std::cout << "[NStar] profile column sizes:\n";
	for (const auto &col : prof_.radial)
	{
		std::cout << "  - " << col.Label() << " : " << col.Size() << "\n";
	}
	std::cout << "------------------------------\n";
}

//--------------------------------------------------------------
// This has to be run so the class
// knows when to initialize all the splines
// void NStar::FinalizeSurface()
// {
// 	PROFILE_FUNCTION();

// 	// Must have rows before we can interpolate
// 	const auto &r = col(Col::R);
// 	if (r.Size() == 0)
// 	{
// 		Z_LOG_ERROR("FinalizeAfterSurface: no rows appended; nothing to finalize.");
// 		surface_ready = false;
// 		return;
// 	}

// 	// Build splines for (M, ν', ρ, ε, p) as functions of r
// 	ds.Interpolate(idx(Col::R), {idx(Col::M), idx(Col::NuPrime), idx(Col::Rho), idx(Col::Eps), idx(Col::P)});

// 	// Build ν(r) from ν'(r) with surface BC (O(N) cumulative integral of segments)
// 	EvaluateNu();

// 	// ..................................................................
// 	//                         B Integrand
// 	// ..................................................................
// 	// Build baryon-number integrand as a dataset for convenience + fast integrate
// 	B_integrand[0] = col(Col::R);
// 	B_integrand[1] = col(Col::R).pow(2);
// 	B_integrand[1] *= 4.0 * M_PI * col(Col::Rho);

// 	B_integrand[1] /= (1.0 - 2.0 * col(Col::M) / col(Col::R)).sqrt();

// 	// converting fm^{-3} to km^{-3}
// 	B_integrand[1] *= FM3_TO_KM3;

// 	B_integrand.Interpolate(0, 1);
// 	// ..................................................................

// 	// Fill the sequence point
// 	// Units must be converted from km
// 	sequence.ec = col(Col::Eps)[0] * Zaki::Physics::INV_FM4_2_G_CM3 /
// 				  Zaki::Physics::INV_FM4_2_INV_KM2;			// g/cm^3
// 	sequence.m = col(Col::M)[-1] / Zaki::Physics::SUN_M_KM; // M_sun
// 	sequence.r = col(Col::R)[-1];							// km
// 	sequence.pc = col(Col::P)[0] * Zaki::Physics::INV_FM4_2_Dyn_CM2 /
// 				  Zaki::Physics::INV_FM4_2_INV_KM2; // dyne/cm^2
// 	// sequence.b  = Find_BaryonNum_Visible() ;
// 	sequence.b = B_integrand.Integrate(1, {col(Col::R)[0], col(Col::R)[-1]});

// 	// RotationSolver updates MomI; returns it
// 	sequence.I = Find_MomInertia();

// 	surface_ready = true;
// }
//--------------------------------------------------------------
// This has to be run so the class
// knows when to initialize all the splines
void NStar::FinalizeSurface()
{
	PROFILE_FUNCTION();

	// ============================================================
	// 1) PROFILE-FIRST PATH
	// ============================================================
	if (!prof_.empty())
	{
		// --------------------------------------------------------
		// 1.a) Build interpolants from profile columns
		// --------------------------------------------------------
		InitInterpolantsFromProfile_();

		const int rcol = prof_.GetColumnIndex(StarProfile::Column::Radius);

		// --------------------------------------------------------
		// 1.b) Build ν(r) from ν′(r) with surface BC
		//      (this will also register interpolation for ν)
		// --------------------------------------------------------
		EvaluateNu(); // profile-aware version we just wrote

		// --------------------------------------------------------
		// 1.c) Build baryon-number integrand
		// --------------------------------------------------------
		if (!prof_.HasColumn(StarProfile::Column::BaryonDensity) ||
			!prof_.HasColumn(StarProfile::Column::Mass))
		{
			Z_LOG_ERROR("Missing nB or M column in StarProfile; cannot build B integrand.");
		}
		else
		{
			Z_LOG_INFO("Building baryon-number integrand from StarProfile.");
			const int mcol = prof_.GetColumnIndex(StarProfile::Column::Mass);
			const int nbcol = prof_.GetColumnIndex(StarProfile::Column::BaryonDensity);

			B_integrand[0] = prof_.radial[rcol];				// r
			B_integrand[1] = prof_.radial[rcol].pow(2);			// r^2
			B_integrand[1] *= 4.0 * M_PI * prof_.radial[nbcol]; // 4π r^2 nB

			// divide by sqrt(1 - 2M/r)
			// M is already in km, r in km
			B_integrand[1] /= (1.0 - 2.0 * prof_.radial[mcol] / prof_.radial[rcol]).sqrt();

			// fm^{-3} → km^{-3}
			B_integrand[1] *= FM3_TO_KM3;

			B_integrand.Interpolate(0, 1);
		}

		// --------------------------------------------------------
		// 1.d) Fill the sequence point from the profile
		// --------------------------------------------------------
		// energy density at center
		if (prof_.HasColumn(StarProfile::Column::EnergyDensity))
		{
			const auto &eps0 = prof_.radial[prof_.GetColumnIndex(StarProfile::Column::EnergyDensity)][0];
			prof_.seq_point.ec = eps0 * Zaki::Physics::INV_FM4_2_G_CM3 /
								 Zaki::Physics::INV_FM4_2_INV_KM2; // g/cm^3
		}
		else
		{
			prof_.seq_point.ec = 0.0;
		}

		// mass and radius at surface
		prof_.seq_point.r = prof_.GetRadius()->operator[](-1); // km

		if (prof_.HasColumn(StarProfile::Column::Mass))
		{
			const auto &m_last = prof_.GetMass()->operator[](-1);
			prof_.seq_point.m = m_last / Zaki::Physics::SUN_M_KM; // M_sun
		}
		else
		{
			prof_.seq_point.m = 0.0; // M_sun
		}

		// central pressure
		if (prof_.HasColumn(StarProfile::Column::Pressure))
		{
			const auto &p0 = prof_.GetPressure()->operator[](0);
			prof_.seq_point.pc = p0 * Zaki::Physics::INV_FM4_2_Dyn_CM2 /
								 Zaki::Physics::INV_FM4_2_INV_KM2; // dyne/cm^2
		}
		else
		{
			prof_.seq_point.pc = 0.0;
		}

		// total baryon number (visible) from integrand
		if (B_integrand[0].Size() > 0)
		{
			prof_.seq_point.b = B_integrand.Integrate(
				1, {B_integrand[0][0], B_integrand[0][-1]});
		}
		else
		{
			prof_.seq_point.b = 0.0;
		}

		// moment of inertia
		prof_.seq_point.I = Find_MomInertia();

		surface_ready = true;
		return;
	}
}

//--------------------------------------------------------------
// Check if surface is initialized
//  Returns true if SurfaceIsReached has been called
//  and internal splines are ready.
//  It is set to true in SurfaceIsReached()
// bool NStar::IsSurfaceFinalized() const noexcept
// {
// 	return surface_ready;
// }

//--------------------------------------------------------------
// Appends one TOV point to the NStar (legacy ds + new StarProfile)
void NStar::Append(const TOVPoint &in_tov)
{
	// ============================================================
	// 1) ---- LEGACY DS PATH (keep as-is) -------------------------
	// ============================================================
	// col(Col::R).PushBack(in_tov.r); // in km

	// col(Col::M).PushBack(
	// 	Zaki::Physics::SUN_M_KM * in_tov.m); // solar-mass → km

	// col(Col::Rho).PushBack(
	// 	in_tov.rho); // in fm^{-3}

	// col(Col::Eps).PushBack(
	// 	in_tov.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3); // to km^{-2}

	// col(Col::P).PushBack(
	// 	in_tov.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2); // to km^{-2}

	// col(Col::NuPrime).PushBack(in_tov.nu_der * 1e5); // 1/cm → 1/km

	// // per-species for legacy ds
	// if (in_tov.rho_i.size() != rho_i_idx.size())
	// {
	// 	Z_LOG_ERROR("Append (legacy ds): mismatch in species count: rho_i.size() = " +
	// 				std::to_string(in_tov.rho_i.size()) + " vs rho_i_idx = " +
	// 				std::to_string(rho_i_idx.size()));
	// }
	// for (std::size_t i = 0; i < in_tov.rho_i.size() && i < rho_i_idx.size(); ++i)
	// {
	// 	ds[rho_i_idx[i]].PushBack(in_tov.rho_i[i]); // fm^{-3}
	// }

	// ============================================================
	// 2) ---- PROFILE PATH (preferred going forward) --------------
	// ============================================================
	// Make sure the radial dataset inside the profile is large enough
	// and has the canonical columns.
	auto &radial = prof_.radial;

	// Now actually append the values (same unit conversions as legacy)
	radial[prof_.idx_r].PushBack(in_tov.r);
	radial[prof_.idx_m].PushBack(Zaki::Physics::SUN_M_KM * in_tov.m); // solar-mass → km
	radial[prof_.idx_nuprime].PushBack(in_tov.nu_der * 1e5);		  // 1/cm→1/km
	radial[prof_.idx_p].PushBack(
		in_tov.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2);
	radial[prof_.idx_eps].PushBack(
		in_tov.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3);
	radial[prof_.idx_nb].PushBack(in_tov.rho);

	// we do not set nu here — EvaluateNu() will overwrite this column later
	radial[prof_.idx_nu].PushBack(0.0);

	// std::cout << "radial[idx_r].Size() - radial[idx_m].Size() = "
	// 		  << radial[prof_.idx_r].Size() - radial[prof_.idx_m].Size() << "\n";

	// ------------------------------------------------------------
	// Compute and append λ to the profile
	// ------------------------------------------------------------
	// Compute the Schwarzschild-like factor:
	const double r_km = in_tov.r;
	const double m_km = Zaki::Physics::SUN_M_KM * in_tov.m;

	double denom = 1.0;
	if (r_km > 0.0)
	{
		denom = 1.0 - 2.0 * m_km / r_km;
		if (denom <= 0.0)
			denom = 1e-15; // protect against log of non-positive
	}

	// g_rr = 1 / (1 - 2M/r)
	// const double grr = 1.0 / denom;

	// λ = 0.5 * ln(g_rr) = -0.5 * ln(denom)
	// const double lambda_geom = 0.5 * std::log(grr);
	// or equivalently:
	const double lambda_geom = -0.5 * std::log(denom);

	radial[prof_.idx_lambda].PushBack(lambda_geom);
	// ------------------------------------------------------------

	// ------------------------------------------------------------
	// 2.a) per-species for the profile
	// ------------------------------------------------------------
	// if (!in_tov.rho_i.empty())
	// {
	// 	// If the profile doesn’t yet know these species, we auto-register them.
	// 	// We assume order-matching: in_tov.rho_i[k] ↔ prof_.species_labels[k]
	// 	if (prof_.species_labels.size() < in_tov.rho_i.size())
	// 	{
	// 		// extend labels and columns
	// 		for (std::size_t k = prof_.species_labels.size();
	// 			 k < in_tov.rho_i.size(); ++k)
	// 		{
	// 			// make up a label if we don't have one
	// 			const std::string lbl = "rho_i_" + std::to_string(k);
	// 			radial.AddColumn(lbl);
	// 			const int col_idx = static_cast<int>(radial.Dim().size() - 1);
	// 			prof_.AddSpecies(lbl, col_idx);
	// 		}
	// 	}

	// now append the data
	for (std::size_t k = 0; k < in_tov.rho_i.size(); ++k)
	{
		const int col_idx = prof_.species_idx[k];
		radial[col_idx].PushBack(in_tov.rho_i[k]); // fm^{-3}
	}
	// }

	// ------------------------------------------------------------
	// 2.b) update surface guess (sequence) incrementally
	// ------------------------------------------------------------
	// We can’t finalize here (that’s FinalizeSurface’s job) but we can
	// at least update the running “last point” so if we inspect sequence
	// before finalization, we get something reasonable.
	prof_.seq_point.r = in_tov.r; // km
	prof_.seq_point.m = in_tov.m; // M_sun  ← in solar masses
								  // pc, ec, B, I will be filled / fixed at finalize time
}

/* ------------- LEGACY VERSION FOR REFERENCE ------------------

void NStar::Append(const TOVPoint &in_tov)
{
	col(Col::R).PushBack(in_tov.r);                                                                       // in km
	col(Col::M).PushBack(Zaki::Physics::SUN_M_KM * in_tov.m);                                             // in km
	col(Col::Rho).PushBack(in_tov.rho);                                                                   // in fm^{-3}
	col(Col::Eps).PushBack(in_tov.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3); // in km^{-2}
	col(Col::P).PushBack(in_tov.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2); // in km^{-2}
	col(Col::NuPrime).PushBack(in_tov.nu_der * 1e+5);                                                     // convert 1/cm to 1/km

	if (in_tov.rho_i.size() != rho_i_idx.size())
	{
		Z_LOG_ERROR("Append: mismatch in species count: rho_i.size() = " +
					std::to_string(in_tov.rho_i.size()) + " vs rho_i_idx = " +
					std::to_string(rho_i_idx.size()));
	}

	for (size_t i = 0; i < in_tov.rho_i.size(); i++)
	{
		// in fm^{-3}
		ds[rho_i_idx[i]].PushBack(in_tov.rho_i[i]);
	}
}

//-------------------------------------------------------------- */
// Sets the work directory for the member objects
void NStar::OnWorkDirChanged(
	const Zaki::String::Directory &in_dir)
{
	std::cout << "[NStar::OnWorkDirChanged] Setting work directory to: "
			  << in_dir << "\n";
	prof_.radial.SetWrkDir(in_dir);
	std::cout << "[NStar::OnWorkDirChanged] prof_.radial.GetWrkDir(): "
			  << prof_.radial.GetWrkDir() << "\n";
}

//--------------------------------------------------------------
/// Similar to the destructor
void NStar::Reset()
{
	// ds.ClearRows();
	prof_.Reset();
	B_integrand.ClearRows();
	// sequence.clear();
	surface_ready = false;
}

//--------------------------------------------------------------
/// Destructor
NStar::~NStar()
{
}

//--------------------------------------------------------------
// /// Evaluate the metric function
// Out of commission from October 28, 2025
// Reason: This was an O(N^2) algorithm due to the repeated integrations
// void NStar::EvaluateNu()
// {
// 	PROFILE_FUNCTION();

// 	// -----------------------------------
// 	// Integrate to find nu(r) :
// 	// -----------------------------------
// 	col(Col::Nu).Resize(col(Col::R).Size());

// 	// Boundary condition
// 	double nu_at_R = 0.5 * log(
// 							   1 - 2 * col(Col::M)[-1] / col(Col::R)[-1]);

// 	// Calculate the surface term first to find out delta_nu_r
// 	col(Col::Nu)[-1] = ds.Integrate(idx(Col::NuPrime), {col(Col::R)[0], col(Col::R)[-1]});
// 	double delta_nu_r = col(Col::Nu)[-1] - nu_at_R;

// 	col(Col::Nu)[-1] = nu_at_R;

// 	// We don't need to calculate the surface term anymore,
// 	//  therefore we have 'i < col(Col::R).Size()-1'
// 	for (size_t i = 0; i < col(Col::R).Size() - 1; ++i)
// 	{
// 		col(Col::Nu)[i] = ds.Integrate(idx(Col::NuPrime), {col(Col::R)[0], col(Col::R)[i]}) - delta_nu_r;
// 	}

// 	ds.Interpolate(0, idx(Col::Nu));
// }

//--------------------------------------------------------------
/// Evaluate the metric function
// void NStar::EvaluateNu()
// {
// 	PROFILE_FUNCTION();

// 	auto &r = col(Col::R);

// 	const std::size_t N = r.Size();
// 	if (N == 0)
// 	{
// 		return;
// 	}

// 	auto &m = col(Col::M);
// 	auto &nup = col(Col::Nu);
// 	nup.Resize(N);

// 	// Surface boundary condition at R
// 	const double R = r[-1];
// 	const double MR = m[-1];
// 	double x = 1.0 - 2.0 * MR / R;
// 	if (x <= 0.0)
// 	{
// 		Z_LOG_ERROR("Non-physical 2M/R >= 1 in EvaluateNu(); clamping.");
// 		x = 1e-15;
// 	}
// 	const double nu_R = 0.5 * std::log(x);

// 	// Build J[i] = int_{r_i}^{R} ν'(r) dr with adjacent-interval integrals
// 	std::vector<double> J(N, 0.0);
// 	// J[N-1] = 0
// 	for (std::size_t i = N - 1; i > 0; --i)
// 	{
// 		// integrate just the small interval [r[i-1], r[i]] using our GSL-backed DataSet
// 		const double seg = ds.Integrate(idx(Col::NuPrime), {r[i - 1], r[i]});
// 		J[i - 1] = J[i] + seg; // accumulate inward
// 	}

// 	// ν(r_i) = ν(R) - J[i]
// 	for (std::size_t i = 0; i < N; ++i)
// 	{
// 		nup[i] = nu_R - J[i];
// 	}

// 	ds.Interpolate(idx(Col::R), idx(Col::Nu));
// }
void NStar::EvaluateNu()
{
	PROFILE_FUNCTION();

	// prefer StarProfile path
	if (!prof_.empty() &&
		prof_.HasColumn(StarProfile::Column::Radius) &&
		prof_.HasColumn(StarProfile::Column::Mass) &&
		prof_.HasColumn(StarProfile::Column::MetricNuPrime))
	{
		auto &radial = prof_.radial;
		const int rcol = prof_.GetColumnIndex(StarProfile::Column::Radius);
		const int mcol = prof_.GetColumnIndex(StarProfile::Column::Mass);
		const int nupcol = prof_.GetColumnIndex(StarProfile::Column::MetricNuPrime);

		// ensure ν column exists (create or clear)
		int nucol = prof_.GetColumnIndex(StarProfile::Column::MetricNu);
		if (!prof_.IsValidColumnIndex(nucol))
		{
			// append a new column for ν(r)
			nucol = static_cast<int>(radial.Dim().size());
			radial.AddColumn("nu");
			prof_.SetColumnIndex(StarProfile::Column::MetricNu, nucol);
		}

		auto &nu = radial[nucol];

		const auto &r = radial[rcol];
		const auto &m = radial[mcol];

		const std::size_t N = r.Size();
		if (N == 0)
			return;
		nu.Resize(N);

		// surface boundary condition
		const double R = r[N - 1];
		const double MR = m[N - 1];
		double x = 1.0 - 2.0 * MR / R;
		if (x <= 0.0)
		{
			Z_LOG_ERROR("Non-physical 2M/R ≥ 1 in EvaluateNu(); clamping.");
			x = 1e-15;
		}
		const double nu_R = 0.5 * std::log(x);

		// accumulate ∫ ν′ dr inward
		std::vector<double> J(N, 0.0);
		for (std::size_t i = N - 1; i > 0; --i)
		{
			const double seg = radial.Integrate(nupcol, {r[i - 1], r[i]});
			J[i - 1] = J[i] + seg;
		}

		for (std::size_t i = 0; i < N; ++i)
			nu[i] = nu_R - J[i];

		// register spline
		radial.Interpolate(rcol, nucol);
		return;
	}

	// // fallback: legacy dataset
	// if (!ds.Dim().empty())
	// {
	//     auto &r = col(Col::R);
	//     const std::size_t N = r.Size();
	//     if (N == 0) return;

	//     auto &m   = col(Col::M);
	//     auto &nu  = col(Col::Nu);
	//     auto &nup = col(Col::NuPrime);
	//     nu.Resize(N);

	//     const double R  = r[N-1];
	//     const double MR = m[N-1];
	//     double x = 1.0 - 2.0 * MR / R;
	//     if (x <= 0.0)
	//     {
	//         Z_LOG_ERROR("Non-physical 2M/R ≥ 1 in EvaluateNu(); clamping.");
	//         x = 1e-15;
	//     }
	//     const double nu_R = 0.5 * std::log(x);

	//     std::vector<double> J(N, 0.0);
	//     for (std::size_t i = N - 1; i > 0; --i)
	//     {
	//         const double seg = ds.Integrate(idx(Col::NuPrime), {r[i-1], r[i]});
	//         J[i-1] = J[i] + seg;
	//     }

	//     for (std::size_t i = 0; i < N; ++i)
	//         nu[i] = nu_R - J[i];

	//     ds.Interpolate(idx(Col::R), idx(Col::Nu));
	// }
}

//--------------------------------------------------------------
/// Metric function as a function of radius (in km)
double NStar::GetMetricNu(const double &in_r) const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("radius must be non-negative.");
		return std::numeric_limits<double>::quiet_NaN();
	}

	// profile path
	if (!prof_.empty() && prof_.HasColumn(StarProfile::Column::Radius) &&
		prof_.HasColumn(StarProfile::Column::MetricNu))
	{
		const auto &rcol = prof_.GetRadius();
		const std::size_t N = rcol->Size();
		if (N == 0)
			return 0.0;
		const double r0 = rcol->operator[](0);
		const double rR = rcol->operator[](N - 1);
		if (in_r < r0 || in_r > rR)
			return 0.0;

		return prof_.radial.Evaluate(
			prof_.GetColumnIndex(StarProfile::Column::MetricNu),
			in_r);
	}

	// legacy path
	// auto dims = ds.Dim();
	// if (!dims.empty())
	// {
	// 	const auto &r = col(Col::R);
	// 	const std::size_t N = r.Size();
	// 	if (N == 0)
	// 		return 0.0;
	// 	if (in_r < r[0] || in_r > r[N - 1])
	// 		return 0.0;
	// 	return ds.Evaluate(idx(Col::Nu), in_r);
	// }

	return 0.0;
}

//--------------------------------------------------------------
// Metric function as a function of radius (in km)
// Zaki::Vector::DataColumn *NStar::GetNu() noexcept
// {
// 	return &col(Col::Nu);
// }
//--------------------------------------------------------------
// Get radius at star surface.
//  returns Radius at surface (km).
// double NStar::RadiusSurface() const noexcept // sequence.r
// {
// 	if (!prof_.empty() && prof_.R > 0.0)
// 		return prof_.R;
// 	return sequence.r;
// }

//--------------------------------------------------------------
// Get mass at star surface.
// returns Mass at surface in solar mass units.
// prof_.M is in km
// sequence.m is in solar mass units
// double NStar::MassSurface() const noexcept // sequence.m
// {
// 	if (!prof_.empty() && prof_.M > 0.0)
// 		return prof_.M / Zaki::Physics::SUN_M_KM;
// 	return sequence.m;
// }

//--------------------------------------------------------------
//  Get the number of radial grid points.
// returns Number of points.
// std::size_t NStar::Size() const noexcept
// {
// 	if (!prof_.empty())
// 		return prof_.size();

// 	// Legacy behavior:
// 	// auto ds_dims = ds.Dim();
// 	// if (!ds_dims.empty())
// 	// 	return ds[0].Size();

// 	return 0;
// }

//--------------------------------------------------------------
// Metric function as a function of radius (in km) - const version
// const Zaki::Vector::DataColumn *NStar::GetNu() const noexcept
// {
// 	return &col(Col::Nu);
// }

//--------------------------------------------------------------
/// Mass (in km) as a function of radius
double NStar::GetMass(const double &in_r) const
{
	{
		if (in_r < 0)
		{
			Z_LOG_ERROR("Radius must be non-negative.");
			return std::numeric_limits<double>::quiet_NaN();
		}

		// ----- profile path -----
		if (!prof_.empty() && prof_.HasColumn(StarProfile::Column::Radius) &&
			prof_.HasColumn(StarProfile::Column::Mass))
		{
			const auto &radial = prof_.radial;
			const auto &rcol = prof_.GetRadius();
			const std::size_t N = rcol->Size();
			if (N == 0)
				return 0.0;

			const double r0 = rcol->operator[](0);
			const double rR = rcol->operator[](N - 1);

			if (in_r < r0)
				return 0.0;
			if (in_r > rR)
				return radial[prof_.GetColumnIndex(StarProfile::Column::Mass)][N - 1];

			return radial.Evaluate(
				prof_.GetColumnIndex(StarProfile::Column::Mass),
				in_r);
		}

		// ----- legacy path -----
		// auto dims = ds.Dim();
		// if (!dims.empty())
		// {
		// 	const auto &r = col(Col::R);
		// 	const std::size_t N = r.Size();
		// 	if (N == 0)
		// 		return 0.0;

		// 	if (in_r < r[0])
		// 		return 0.0;
		// 	if (in_r > r[N - 1])
		// 		return col(Col::M)[N - 1];

		// 	return ds.Evaluate(idx(Col::M), in_r);
		// }

		return 0.0;
	}
}

//--------------------------------------------------------------
// // Mass (in km) as a function of radius
// Zaki::Vector::DataColumn *NStar::GetMass() noexcept
// {
// 	return &col(Col::M);
// }
// //--------------------------------------------------------------
// // Mass (in km) as a function of radius - const version
// const Zaki::Vector::DataColumn *NStar::GetMass() const noexcept
// {
// 	return &col(Col::M);
// }

//--------------------------------------------------------------
// /// Returns the radius dataset
// Zaki::Vector::DataColumn *NStar::GetRadius() noexcept
// {
// 	return &col(Col::R);
// }

// //--------------------------------------------------------------
// /// Returns the radius dataset - const version
// const Zaki::Vector::DataColumn *NStar::GetRadius() const noexcept
// {
// 	return &col(Col::R);
// }

//--------------------------------------------------------------
/// Baryon number density (fm^{-3}) as a function of radius
double NStar::GetBaryonDensity(const double &in_r) const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Radius must be non-negative.");
		return std::numeric_limits<double>::quiet_NaN();
	}

	// profile path
	if (!prof_.empty() && prof_.HasColumn(StarProfile::Column::Radius) &&
		prof_.HasColumn(StarProfile::Column::BaryonDensity))
	{
		const auto &rcol = prof_.GetRadius();
		const std::size_t N = rcol->Size();
		if (N == 0)
			return 0.0;

		const double r0 = rcol->operator[](0);
		const double rR = rcol->operator[](N - 1);
		if (in_r < r0 || in_r > rR)
			return 0.0;

		return prof_.radial.Evaluate(
			prof_.GetColumnIndex(StarProfile::Column::BaryonDensity),
			in_r);
	}

	// legacy path
	// auto dims = ds.Dim();
	// if (!dims.empty())
	// {
	// 	const auto &r = col(Col::R);
	// 	const std::size_t N = r.Size();
	// 	if (N == 0)
	// 		return 0.0;

	// 	if (in_r < r[0] || in_r > r[N - 1])
	// 		return 0.0;

	// 	return ds.Evaluate(idx(Col::Rho), in_r);
	// }

	return 0.0;
}

//--------------------------------------------------------------
// Baryon number density (fm^{-3}) as a function of radius
// Zaki::Vector::DataColumn *NStar::GetRho() noexcept
// {
// 	return &col(Col::Rho);
// }

//--------------------------------------------------------------
// Baryon number density (fm^{-3}) as a function of radius - const version
// const Zaki::Vector::DataColumn *NStar::GetRho() const noexcept
// {
// 	return &col(Col::Rho);
// }

//--------------------------------------------------------------
// Check if density data exists for a labeled species.
//  label  Species label string.
// Returns True if data exists, false otherwise.
// bool NStar::HasRho_i(std::string_view label) const noexcept
// {
// 	for (auto &&c : rho_i_idx)
// 	{
// 		if (ds[c].label == label)
// 		{
// 			return true;
// 		}
// 	}

// 	return false;
// }

//--------------------------------------------------------------
// Visible baryon number density (fm^{-3}) as a function of radius
// for a specific species labeled as (in_label)
// Zaki::Vector::DataColumn *NStar::GetRho_i(const std::string &in_label) noexcept
// {

// 	for (auto &&c : rho_i_idx)
// 	{
// 		if (ds[c].label == in_label)
// 		{
// 			return &ds[c];
// 		}
// 	}

// 	Z_LOG_ERROR("Species density with label '" +
// 				in_label + "' was not found. Returning nullptr.");

// 	return nullptr;
// }

//--------------------------------------------------------------
// Visible baryon number density (fm^{-3}) as a function of radius
// for a specific species labeled as (in_label) - const version
// const Zaki::Vector::DataColumn *
// NStar::GetRho_i(const std::string_view &in_label) const
// {
// 	for (auto &&c : rho_i_idx)
// 	{
// 		if (ds[c].label == in_label)
// 			return &ds[c];
// 	}
// 	Z_LOG_ERROR(std::string("Species density with label '") +
// 				std::string(in_label) + "' was not found. Returning nullptr.");
// 	return nullptr;
// }

// ------------------------------------------------------------
// profile-aware species
// ------------------------------------------------------------
// Check if density data exists for a labeled species.
//  label  Species label string.
// Returns True if data exists, false otherwise.
// bool NStar::HasRho_i(std::string_view label) const noexcept
// {
// 	if (!prof_.empty() && prof_.HasSpecies(std::string(label)))
// 		return true;
// 	return false;
// }

// ------------------------------------------------------------
// Visible baryon number density (fm^{-3}) as a function of radius
// for a specific species labeled as (in_label) - const version
const Zaki::Vector::DataColumn *
NStar::GetRho_i(const std::string_view &in_label) const
{
	if (prof_.empty())
		return nullptr;

	// use the pointer-returning helper inside StarProfile
	const auto *col = prof_.GetSpeciesPtr(std::string(in_label));

	return col; // either a valid pointer to the real column or nullptr
}

// ------------------------------------------------------------
// Visible baryon number density (fm^{-3}) as a function of radius
// for a specific species labeled as (in_label)
// Zaki::Vector::DataColumn *
// NStar::GetRho_i(const std::string &in_label)
// {
// 	if (prof_.empty())
// 		return nullptr;

// 	return prof_.GetSpeciesPtr(in_label); // calls the non-const overload
// }

//--------------------------------------------------------------
/// Energy density (in km^{-2}) as a function of radius
double NStar::GetEnergyDensity(const double &in_r) const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Radius must be non-negative.");
		return std::numeric_limits<double>::quiet_NaN();
	}

	// profile path
	if (!prof_.empty() && prof_.HasColumn(StarProfile::Column::Radius) &&
		prof_.HasColumn(StarProfile::Column::EnergyDensity))
	{
		const auto &rcol = prof_.GetRadius();
		const std::size_t N = rcol->Size();
		if (N == 0)
			return 0.0;
		const double r0 = rcol->operator[](0);
		const double rR = rcol->operator[](N - 1);
		if (in_r < r0 || in_r > rR)
			return 0.0;

		return prof_.radial.Evaluate(
			prof_.GetColumnIndex(StarProfile::Column::EnergyDensity),
			in_r);
	}

	// legacy path
	// auto dims = ds.Dim();
	// if (!dims.empty())
	// {
	// 	const auto &r = col(Col::R);
	// 	const std::size_t N = r.Size();
	// 	if (N == 0)
	// 		return 0.0;
	// 	if (in_r < r[0] || in_r > r[N - 1])
	// 		return 0.0;
	// 	return ds.Evaluate(idx(Col::Eps), in_r);
	// }

	return 0.0;
}

//--------------------------------------------------------------
// Energy density (in km^{-2}) as a function of radius
// Zaki::Vector::DataColumn *
// NStar::GetEps() noexcept
// {
// 	return &col(Col::Eps);
// }
// //--------------------------------------------------------------
// // Energy density (in km^{-2}) as a function of radius - const version
// const Zaki::Vector::DataColumn *
// NStar::GetEps() const noexcept
// {
// 	return &col(Col::Eps);
// }

//--------------------------------------------------------------
/// Pressure (in km^{-2}) as a function of radius
double NStar::GetPressure(const double &in_r) const
{
	if (in_r < 0)
	{
		Z_LOG_ERROR("Radius must be non-negative.");
		return std::numeric_limits<double>::quiet_NaN();
	}

	// profile path
	if (!prof_.empty() && prof_.HasColumn(StarProfile::Column::Radius) &&
		prof_.HasColumn(StarProfile::Column::Pressure))
	{
		const auto &rcol = prof_.GetRadius();
		const std::size_t N = rcol->Size();
		if (N == 0)
			return 0.0;
		const double r0 = rcol->operator[](0);
		const double rR = rcol->operator[](N - 1);
		if (in_r < r0 || in_r > rR)
			return 0.0;

		return prof_.radial.Evaluate(
			prof_.GetColumnIndex(StarProfile::Column::Pressure),
			in_r);
	}

	// legacy path
	// auto dims = ds.Dim();
	// if (!dims.empty())
	// {
	// 	const auto &r = col(Col::R);
	// 	const std::size_t N = r.Size();
	// 	if (N == 0)
	// 		return 0.0;
	// 	if (in_r < r[0] || in_r > r[N - 1])
	// 		return 0.0;
	// 	return ds.Evaluate(idx(Col::P), in_r);
	// }

	return 0.0;
}

//--------------------------------------------------------------
// Pressure (in km^{-2}) as a function of radius
// Zaki::Vector::DataColumn *
// NStar::GetPress() noexcept
// {
// 	return &col(Col::P);
// }

//--------------------------------------------------------------
// Pressure (in km^{-2}) as a function of radius - const version
// const Zaki::Vector::DataColumn *
// NStar::GetPress() const noexcept
// {
// 	return &col(Col::P);
// }

// ------------------------------------------------------------
// Sequence forwarding
const SeqPoint &NStar::GetSequence() const noexcept
{
	if (!prof_.empty())
		return prof_.seq_point;

	// fallback: static empty
	static SeqPoint empty{};
	return empty;
}
// ------------------------------------------------------------
// Sequence forwarding (non-const)
SeqPoint &NStar::GetSequence() noexcept
{
	if (!prof_.empty())
		return prof_.seq_point;

	// fallback: static singleton we can write to (rare path)
	static SeqPoint empty{};
	return empty;
}

//--------------------------------------------------------------
// Baryon number integrand
double NStar::BaryonNumIntegrand(double in_r) const
{
	// if (in_r <= 0)
	// 	return 0.0; // center contributes 0 anyway (r^2 factor)
	// const double nb = GetRho(in_r);
	// if (nb <= 0)
	// 	return 0.0; // outside star or invalid -> no contribution

	// const double M = GetMass(in_r);
	// const double f = 1.0 - 2.0 * M / in_r;

	// if (f <= 0.0)
	// 	return 0.0; // numerical guard; physically should not be ≤0 inside R

	// return 4.0 * M_PI * in_r * in_r * nb / std::sqrt(f);

	// double out = in_r * in_r;
	// out *= 4 * M_PI * GetRho(in_r);

	// // // converting fm^{-3} to km^{-3}
	// // out *= pow(10, 54) ;

	// out /= sqrt(1 - 2 * GetMass(in_r) / in_r);

	// return out;
	if (in_r <= 0.0)
		return 0.0;

	const double nb = GetBaryonDensity(in_r);
	if (nb <= 0.0)
		return 0.0;

	const double M = GetMass(in_r);
	const double f = 1.0 - 2.0 * M / in_r;
	if (f <= 0.0)
		return 0.0;

	return 4.0 * M_PI * in_r * in_r * nb / std::sqrt(f);
}

//--------------------------------------------------------------
/// Total moment of inertia (in km^3)
// MomInertiaIntegrand is wrong!
double NStar::Find_MomInertia()
{
	PROFILE_FUNCTION();

	rot_solver.FindNMomInertia();

	return MomI;
}

//--------------------------------------------------------------
/// Precision for printing the profile
/// by default it is set to '9' digits
void NStar::SetProfilePrecision(const int &in_prec)
{
	prof_.SetProfilePrecision(in_prec);
}

//--------------------------------------------------------------
/// Export the profile to file
// void NStar::Export(const Zaki::String::Directory &rel_path)
// {
// 	if (prof_.empty())
// 	{
// 		Z_LOG_ERROR("StarProfile is empty; nothing to export.");
// 		return;
// 	}

// 	const auto &radial = prof_.radial;
// 	const auto dims = radial.Dim();

// 	Z_LOG_INFO("radial.Dim().size() = " +
// 			   std::to_string(dims.size()));

// 	if (!dims.empty())
// 	{
// 		Z_LOG_INFO("first column size = " +
// 				   std::to_string(radial[0].Size()));
// 	}
// 	// if (!prof_.empty())
// 	// {
// 	// 	prof_.Export(in_dir);
// 	// 	return;
// 	// }
// 	// Make sure the profile's DataSet uses the same work directory as NStar.
// 	// (OnWorkDirChanged() already does this when SetWrkDir() is called,
// 	// but calling it here is cheap and makes this path bulletproof.)
// 	prof_.radial.SetWrkDir(wrk_dir_);

// 	Z_LOG_INFO("wrk_dir_ = " + wrk_dir_.Str());
// 	Z_LOG_INFO("exporting profile to relative path: " + rel_path.Str());

// 	// Delegate actual write
// 	prof_.Export(rel_path);

// 	// For sanity, print what DataSet thinks its work directory is:
// 	Z_LOG_INFO("radial.GetWrkDir() = " + prof_.radial.GetWrkDir().Str());
// }

void NStar::Export(const Zaki::String::Directory &rel_path)
{
	PROFILE_FUNCTION();

	if (prof_.empty())
	{
		Z_LOG_ERROR("StarProfile is empty; nothing to export.");
		return;
	}

	// 1) Ensure the profile's DataSet uses the same work directory as NStar.
	prof_.radial.SetWrkDir(wrk_dir_);

	// Z_LOG_INFO("wrk_dir_ = " + wrk_dir_.Str());
	// Z_LOG_INFO("exporting profile to relative path: " + rel_path.Str());

	// 2) Build the full path and force-create its parent directory
	Zaki::String::Directory full_path = wrk_dir_ + rel_path;
	Zaki::String::Directory parent_dir = full_path.ParentDir();

	// Z_LOG_INFO("full_path      = " + full_path.Str());
	// Z_LOG_INFO("parent_dir     = " + parent_dir.Str());

	// parent_dir.Create(); // should create .../results/tov_debug if missing

	// 3) Write a tiny debug file *directly* using std::ofstream
	// {
	// 	const std::string dbg_name = full_path.Str() + ".debug";

	// 	std::ofstream dbg(dbg_name);
	// 	if (!dbg)
	// 	{
	// 		Z_LOG_ERROR("FAILED to open debug file for write: " + dbg_name);
	// 	}
	// 	else
	// 	{
	// 		dbg << "# NStar export debug file\n";
	// 		dbg << "# wrk_dir_   = " << wrk_dir_.Str() << "\n";
	// 		dbg << "# rel_path   = " << rel_path.Str() << "\n";
	// 		dbg << "# full_path  = " << full_path.Str() << "\n";
	// 		dbg << "# " << prof_.column_count() << " columns in profile\n";
	// 		dbg << "# " << prof_.size() << " rows in profile\n";
	// 		dbg << "# " << prof_.GetRadius()->Min() << " km ≤ r ≤ "
	// 			<< prof_.GetRadius()->Max() << " km\n";
	// 		dbg.close();
	// 		Z_LOG_INFO("wrote debug file: " + dbg_name);
	// 	}
	// }

	// 4) Delegate actual TSV export via StarProfile
	prof_.Export(rel_path);
	// prof_.Export(full_path);

	// Z_LOG_INFO("radial.GetWrkDir() = " +
	// 		   prof_.radial.GetWrkDir().Str());
}

//--------------------------------------------------------------
//--------------------------------------------------------------
// NStar::SolveTOV_Profile — run TOV internally and build profile
//--------------------------------------------------------------
int NStar::SolveTOV_Profile(const Zaki::String::Directory &eos_file,
							double target_M_solar,
							const Zaki::String::Directory &rel_out_dir)
{
	PROFILE_FUNCTION();

	// Fresh start: clear any existing profile / sequence / B-integrand.
	Reset();

	// ------------------------------------------------------------
	// 1) Set up a local TOVSolver
	// ------------------------------------------------------------
	CompactStar::Core::TOVSolver tov;

	// Working directory for the solver:
	//   base = this NStar's wrk_dir_
	//   plus optional relative subdirectory (e.g. "tov_debug/")
	Zaki::String::Directory out_dir = wrk_dir_ + rel_out_dir; // no-op if rel_out_dir is ""
	// out_dir = out_dir + rel_out_dir; // no-op if rel_out_dir is ""

	// std::cout << "NStar::SolveTOV_Profile: setting TOVSolver work dir to: "
	// 		  << out_dir << std::endl;
	tov.SetWrkDir(out_dir);

	// Import the EOS table directly from the provided file path.
	// Caller is responsible for passing something like:
	//   eos_root + eos_name + "/" + eos_name + ".eos"
	tov.ImportEOS(eos_file, true);

	// (Optional) If we want to use the same profile precision as NStar's
	// exports, we could mirror that here by adding a getter for
	// profile_precision and calling tov.SetProfilePrecision(...).
	// For now we just use TOVSolver's internal default.

	// ------------------------------------------------------------
	// 2) Single-star solve to target mass
	// ------------------------------------------------------------
	std::vector<CompactStar::Core::TOVPoint> tov_points;
	std::vector<std::string> species_labels;

	const int n_pts = tov.SolveToProfile(target_M_solar,
										 tov_points,
										 &species_labels);
	if (n_pts <= 0 || tov_points.empty())
	{
		Z_LOG_ERROR("SolveToProfile failed for "
					"target mass = " +
					std::to_string(target_M_solar) + " Msun "
													 "with EOS file: " +
					eos_file.Str());
		surface_ready = false;
		return 0;
	}

	// ------------------------------------------------------------
	// 3) Use the existing builder to populate StarProfile
	// ------------------------------------------------------------
	// BuildFromTOV:
	//  - fills prof_.radial (r, m, nu', p, eps, rho, nu, lambda, species…)
	//  - constructs B_integrand and integrates B
	//  - fills prof_.seq_point (ec, pc, M, R, B, I)
	//  - sets prof_.M, prof_.R, prof_.z_surf
	//  - sets surface_ready = true
	BuildFromTOV(tov_points,
				 species_labels.empty() ? nullptr : &species_labels);

	// DO NOT call InitInterpolantsFromProfile_();
	// DO NOT call FinalizeSurface();
	//
	// BuildFromTOV already:
	//   - resets NStar
	//   - lays out columns
	//   - interpolates
	//   - EvaluateNu
	//   - builds B_integrand
	//   - sets seq_point, M, R, z_surf
	//   - sets surface_ready = true

	return n_pts;
}

//--------------------------------------------------------------
//==============================================================
