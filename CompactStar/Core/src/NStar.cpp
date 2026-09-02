/*
  NStar class
*/

// #include <gsl/gsl_math.h>

#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <Zaki/Util/Instrumentor.hpp>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Geometry.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

#include <cmath>
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

	// Batch all profile mutations in this build into ONE version bump.
	auto edit = prof_.Edit();
	// ============================================================
	// 1) PROFILE-FIRST BUILD
	// ============================================================
	{
		auto &radial = prof_.RadialMutable();
		radial.ClearRows();
		radial.Reserve(8 + n_species, n_rows);
		B_integrand.Reserve(2, n_rows);

		// label the ones we just reserved
		radial[0].SetLabel("r(km)");
		prof_.SetColumnIndex(StarProfile::Column::Radius, 0);

		radial[1].SetLabel("m(km)");
		prof_.SetColumnIndex(StarProfile::Column::Mass, 1);

		radial[2].SetLabel("nu_prime(km^-1)");
		prof_.SetColumnIndex(StarProfile::Column::MetricNuPrime, 2);

		radial[3].SetLabel("p(km^-2)");
		prof_.SetColumnIndex(StarProfile::Column::Pressure, 3);

		radial[4].SetLabel("eps(km^-2)");
		prof_.SetColumnIndex(StarProfile::Column::EnergyDensity, 4);

		radial[5].SetLabel("nB(fm^-3)");
		prof_.SetColumnIndex(StarProfile::Column::BaryonDensity, 5);

		radial[6].SetLabel("nu");
		prof_.SetColumnIndex(StarProfile::Column::MetricNu, 6);

		radial[7].SetLabel("lambda");
		prof_.SetColumnIndex(StarProfile::Column::MetricLambda, 7);

		// species (after the fixed ones)
		// species columns: indices 8..(8+n_species-1)
		// prof_.species_labels.clear();
		// prof_.species_idx.clear();
		// prof_.species_labels.reserve(n_species);
		// prof_.species_idx.reserve(n_species);

		// for (std::size_t j = 0; j < n_species; ++j)
		// {
		// 	std::string lbl = (species_labels && j < species_labels->size())
		// 						  ? (*species_labels)[j]
		// 						  : "rho_i_" + std::to_string(j);

		// 	// const int col_idx = static_cast<int>(radial.Dim().size()); // REMOVE THIS LINE
		// 	const int col_idx = 8 + static_cast<int>(j);
		// 	radial[col_idx].SetLabel(lbl);
		// 	// radial.AddColumn(lbl); // REMOVE THIS LINE
		// 	prof_.AddSpecies(lbl, col_idx);
		// }
		// Species registry rebuild (no public ClearSpecies() yet)
		// prof_.species_labels.clear();
		// prof_.species_idx.clear();
		// prof_.species_labels.reserve(n_species);
		// prof_.species_idx.reserve(n_species);
		prof_.ResetSpecies(n_species);

		for (std::size_t j = 0; j < n_species; ++j)
		{
			std::string lbl = (species_labels && j < species_labels->size())
								  ? (*species_labels)[j]
								  : ("rho_i_" + std::to_string(j));

			const int col_idx = 8 + static_cast<int>(j);
			radial[col_idx].SetLabel(lbl);

			// Register species mapping (creates/updates entry)
			prof_.SetSpeciesColumn(lbl, col_idx);
		}

		// ------------------------------------------------------------
		const int idx_r = prof_.GetColumnIndex(StarProfile::Column::Radius);
		const int idx_m = prof_.GetColumnIndex(StarProfile::Column::Mass);
		const int idx_nup = prof_.GetColumnIndex(StarProfile::Column::MetricNuPrime);
		const int idx_p = prof_.GetColumnIndex(StarProfile::Column::Pressure);
		const int idx_eps = prof_.GetColumnIndex(StarProfile::Column::EnergyDensity);
		const int idx_nb = prof_.GetColumnIndex(StarProfile::Column::BaryonDensity);
		const int idx_nu = prof_.GetColumnIndex(StarProfile::Column::MetricNu);
		const int idx_lam = prof_.GetColumnIndex(StarProfile::Column::MetricLambda);
		// ------------------------------------------------------------
		// now fill rows
		for (const auto &tp : in_tov)
		{
			// r (km)
			double r_km = tp.r;
			radial[idx_r].PushBack(r_km);

			double m_km = Zaki::Physics::SUN_M_KM * tp.m;
			// m (km): tp.m is in solar masses
			radial[idx_m].PushBack(m_km);

			// nu' (1/cm -> 1/km)
			radial[idx_nup].PushBack(tp.nu_der * 1e5);

			// p (→ km^-2)
			radial[idx_p].PushBack(
				tp.p * Zaki::Physics::INV_FM4_2_INV_KM2 /
				Zaki::Physics::INV_FM4_2_Dyn_CM2);

			// eps (→ km^-2)
			radial[idx_eps].PushBack(
				tp.e * Zaki::Physics::INV_FM4_2_INV_KM2 /
				Zaki::Physics::INV_FM4_2_G_CM3);

			// nB (fm^-3)
			radial[idx_nb].PushBack(tp.rho);

			// nu (will be built later)
			radial[idx_nu].PushBack(0.0);

			// ------------------------------------------------------------
			// Compute and append λ to the profile
			// ------------------------------------------------------------
			// Compute the Schwarzschild-like factor:
			// const double r_km = tp.r;
			// const double m_km = Zaki::Physics::SUN_M_KM * tp.m;

			// ADR-0004: the metric factor has ONE mathematical owner,
			// CompactStar/Geometry.hpp. This site no longer defines it.
			//
			// g_rr = 1/(1 - 2m/r);  lambda = 0.5*ln(g_rr) = -0.5*ln(1 - 2m/r).
			//
			// Bit-identical to the former inline form on every ordinary node: the primitive
			// performs the same `1.0 - 2.0 * m_km / r_km` subtraction and the same
			// `-0.5 * std::log(denom)`. The former `denom <= 0 -> 1e-15` clamp is gone; under
			// the accepted hybrid physical-domain contract (ADR-0004 §0-Q3) such input fails
			// closed rather than being regularized. It was unreachable here in any case:
			// measured max 2m/r = 0.481 across the authenticated stars (ADR-0004 §11).
			const double lambda_geom = CompactStar::Geometry::Lambda(r_km, m_km);

			radial[idx_lam].PushBack(lambda_geom);
			// ------------------------------------------------------------

			// ------------------------------
			// species values (pad with 0.0)
			if (!tp.rho_i.empty())
			{
				for (std::size_t j = 0; j < n_species; ++j)
				{
					const double val = (j < tp.rho_i.size()) ? tp.rho_i[j] : 0.0;
					// const int col_idx = prof_.species_idx[j];
					// radial[col_idx].PushBack(val);
					const int col_idx = prof_.SpeciesColumnIndex(j);
					if (col_idx >= 0)
						radial[col_idx].PushBack(val);
				}
			}
			else
			{
				// no species in this row → append zeros to all species cols
				for (std::size_t j = 0; j < n_species; ++j)
				{
					// const int col_idx = prof_.species_idx[j];
					// radial[col_idx].PushBack(0.0);
					const int col_idx = prof_.SpeciesColumnIndex(j);
					if (col_idx >= 0)
						radial[col_idx].PushBack(0.0);
				}
			}
		}

		// interpolate the columns we actually have
		{
			std::vector<int> interp_cols;
			interp_cols.push_back(idx_m);
			interp_cols.push_back(idx_nup);
			interp_cols.push_back(idx_nb);
			interp_cols.push_back(idx_eps);
			interp_cols.push_back(idx_p);
			radial.Interpolate(idx_r, interp_cols);
		}

		// build ν(r) from ν'(r)
		EvaluateNu(); // profile-aware version

		// build baryon number integrand from profile
		{
			B_integrand[0] = radial[idx_r];
			B_integrand[0].SetLabel("r(km)");

			// ADR-0004: the proper-volume measure w_V = 4*pi*r^2*exp(Lambda) comes from the
			// single mathematical owner. This site owns only its own physics factor (n_B) and
			// its own unit conversion (FM3_TO_KM3 = 1e54, which belongs to INV-14, NOT to dV).
			//
			// The former inline form divided by `(1 - 2m/r).sqrt()` through
			// Zaki::Vector::DataColumn::operator/=, which silently substitutes 1 for a zero
			// divisor -- dropping the metric factor at that node with only a log line. The
			// primitive fails closed instead.
			//
			// This reassociates the arithmetic and is therefore NOT bit-identical. The
			// tolerance |dB|/B <= 1.0e-15 was predeclared in ADR-0004 §7.1 from the operation
			// counts, before implementation; §7.2 measured at most 1.64e-16 (1 ULP).
			// Guarded by the `baryon_number_cmf` reference test.
			const std::size_t n_nodes = radial[idx_r].Size();
			Zaki::Vector::DataColumn wV;
			wV.Reserve(n_nodes);
			for (std::size_t i = 0; i < n_nodes; ++i)
				wV.PushBack(CompactStar::Geometry::ProperVolumeWeight(radial[idx_r][i],
																	  radial[idx_m][i]));

			B_integrand[1] = wV * radial[idx_nb];
			B_integrand[1].SetLabel("B_v");
			B_integrand[1] *= FM3_TO_KM3;
			B_integrand.Interpolate(0, 1);
		}

		// fill profile's sequence
		// prof_.seq_point.clear();
		auto &seq = prof_.SeqMutable();
		seq.clear();
		// central energy density
		seq.ec = radial[idx_eps][0] *
				 Zaki::Physics::INV_FM4_2_G_CM3 /
				 Zaki::Physics::INV_FM4_2_INV_KM2; // g/cm^3

		// surface radius
		// prof_.seq_point.r = radial[idx_r][-1]; // km
		// prof_.R = prof_.seq_point.r;
		seq.r = radial[idx_r][-1]; // km
		const double R_km = seq.r;

		// surface mass (in both solar masses and km)
		const double Msurf_km = radial[idx_m][-1];					  // km
		const double Msurf_Msun = Msurf_km / Zaki::Physics::SUN_M_KM; // Msun

		// prof_.M = Msurf_km;				// km
		// prof_.seq_point.m = Msurf_Msun; // Msun
		seq.m = Msurf_Msun;

		// central pressure
		seq.pc = radial[idx_p][0] *
				 Zaki::Physics::INV_FM4_2_Dyn_CM2 /
				 Zaki::Physics::INV_FM4_2_INV_KM2; // dyne/cm^2

		// baryon number
		seq.b = B_integrand.Integrate(1,
									  {radial[idx_r][0],
									   radial[idx_r][-1]});

		// moment of inertia
		seq.I = Find_MomInertia();

		// surface redshift factor (if we have ν)
		// if (prof_.HasColumn(StarProfile::Column::MetricNu))
		// {
		// 	const auto &nu_col = radial[prof_.idx_nu];
		// 	if (nu_col.Size() > 0)
		// 		prof_.z_surf = std::exp(nu_col[-1]);
		// 	else
		// 		prof_.z_surf = 0.0;
		// }

		double zsurf = 0.0;
		if (prof_.HasColumn(StarProfile::Column::MetricNu))
		{
			const auto &nu_col = radial[idx_nu];
			if (nu_col.Size() > 0)
				zsurf = std::exp(nu_col[-1]);
		}

		prof_.SetSurfaceScalars(Msurf_km, R_km, zsurf);

		surface_ready = true;
	}
}
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
	Reset();				  // clears ds, B_integrand, sequence, rho_i_idx, prof_
	auto edit = prof_.Edit(); // batches all profile mutations below into one bump
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
		auto &radial = prof_.RadialMutable();
		radial.ClearRows();
		radial.Reserve(8 + n_species, n_rows_expect);

		// label the ones we just reserved
		radial[0].SetLabel("r(km)");
		prof_.SetColumnIndex(StarProfile::Column::Radius, 0);

		radial[1].SetLabel("m(km)");
		prof_.SetColumnIndex(StarProfile::Column::Mass, 1);

		radial[2].SetLabel("nu_prime(km^-1)");
		prof_.SetColumnIndex(StarProfile::Column::MetricNuPrime, 2);

		radial[3].SetLabel("p(km^-2)");
		prof_.SetColumnIndex(StarProfile::Column::Pressure, 3);

		radial[4].SetLabel("eps(km^-2)");
		prof_.SetColumnIndex(StarProfile::Column::EnergyDensity, 4);

		radial[5].SetLabel("nB(fm^-3)");
		prof_.SetColumnIndex(StarProfile::Column::BaryonDensity, 5);

		radial[6].SetLabel("nu");
		prof_.SetColumnIndex(StarProfile::Column::MetricNu, 6);

		radial[7].SetLabel("lambda");
		prof_.SetColumnIndex(StarProfile::Column::MetricLambda, 7);

		// Species columns start right after the 8 fixed ones
		// prof_.species_labels.clear();
		// prof_.species_idx.clear();
		// prof_.species_labels.reserve(n_species);
		// prof_.species_idx.reserve(n_species);
		prof_.ResetSpecies(n_species);

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
		// prof_.seq_point.clear();
		auto &seq = prof_.SeqMutable();
		seq.clear();

		// also set surface M, R, z_surf to 0 here
		// prof_.M = 0.0; // in solar masses
		// prof_.R = 0.0;
		// prof_.z_surf = 0.0;
		prof_.SetSurfaceScalars(0.0, 0.0, 0.0);
	}
}

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

	auto &radial = prof_.RadialMutable();

	if (radial[rcol].Size() == 0)
	{
		Z_LOG_ERROR("StarProfile radius column is empty; nothing to finalize.");
		surface_ready = false;
		return;
	}

	if (prof_.HasColumn(StarProfile::Column::Mass))
		radial.Interpolate(rcol, prof_.GetColumnIndex(StarProfile::Column::Mass));

	if (prof_.HasColumn(StarProfile::Column::Pressure))
		radial.Interpolate(rcol, prof_.GetColumnIndex(StarProfile::Column::Pressure));

	if (prof_.HasColumn(StarProfile::Column::EnergyDensity))
		radial.Interpolate(rcol, prof_.GetColumnIndex(StarProfile::Column::EnergyDensity));

	if (prof_.HasColumn(StarProfile::Column::BaryonDensity))
		radial.Interpolate(rcol, prof_.GetColumnIndex(StarProfile::Column::BaryonDensity));

	if (prof_.HasColumn(StarProfile::Column::MetricNuPrime))
	{
		const int nupcol = prof_.GetColumnIndex(StarProfile::Column::MetricNuPrime);
		if (radial[nupcol].Size() == radial[rcol].Size())
			radial.Interpolate(rcol, nupcol);
		else
			Z_LOG_ERROR("MetricNuPrime column size mismatch; skipping ν′ interpolation.");
	}

	// if (prof_.HasColumn(StarProfile::Column::MetricLambda))
	// 	radial.Interpolate(rcol, prof_.GetColumnIndex(StarProfile::Column::MetricLambda));

	if (prof_.HasColumn(StarProfile::Column::MetricLambda))
	{
		const int lcol = prof_.GetColumnIndex(StarProfile::Column::MetricLambda);
		if (radial[lcol].Size() == radial[rcol].Size())
			radial.Interpolate(rcol, lcol);
		else
			Z_LOG_ERROR("lambda column size mismatch; skipping lambda interpolation.");
	}

	// per-species
	// for (std::size_t i = 0; i < prof_.species_idx.size(); ++i)
	// {
	// 	const int scol = prof_.species_idx[i];
	// 	if (prof_.IsValidColumnIndex(scol))
	// 		radial.Interpolate(rcol, scol);
	// }
	for (std::size_t i = 0; i < prof_.SpeciesCount(); ++i)
	{
		const int scol = prof_.SpeciesColumnIndex(i);
		if (prof_.IsValidColumnIndex(scol))
			radial.Interpolate(rcol, scol);
	}
}

//--------------------------------------------------------------
// Print profile column sizes
void NStar::PrintProfileColumnSizes() const
{
	std::cout << "[NStar] profile column sizes:\n";
	for (const auto &col : prof_.Radial())
	{
		std::cout << "  - " << col.Label() << " : " << col.Size() << "\n";
	}
	std::cout << "------------------------------\n";
}

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
		auto edit = prof_.Edit(); // batches all touches below into one bump
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

			// ADR-0004: the proper-volume measure w_V = 4*pi*r^2*exp(Lambda) comes from the
			// single mathematical owner. This site owns only its own physics factor (n_B) and
			// its own unit conversion (FM3_TO_KM3 = 1e54, which belongs to INV-14, NOT to dV).
			//
			// The former inline form divided by `(1 - 2m/r).sqrt()` through
			// Zaki::Vector::DataColumn::operator/=, which silently substitutes 1 for a zero
			// divisor -- dropping the metric factor at that node with only a log line. The
			// primitive fails closed instead.
			//
			// Composition deliberately mirrors the already-validated BuildFromTOV form, so the
			// two ordinary visible-sector NStar construction paths now share one mathematical
			// owner AND one arithmetic ordering. Governed by the |dB|/B <= 1.0e-15 predeclared
			// in ADR-0004 §7.1 before any implementation existed. Deferred out of Phase 3D
			// because Path 1 had no coverage; Phase 3E-0 supplied it, so this migrated in 3E-I2.
			const std::size_t n_nodes = prof_.Radial()[rcol].Size();
			Zaki::Vector::DataColumn wV;
			wV.Reserve(n_nodes);
			for (std::size_t i = 0; i < n_nodes; ++i)
				wV.PushBack(CompactStar::Geometry::ProperVolumeWeight(
					prof_.Radial()[rcol][i], prof_.Radial()[mcol][i]));

			B_integrand[0] = prof_.Radial()[rcol]; // r
			B_integrand[1] = wV * prof_.Radial()[nbcol];
			// fm^{-3} → km^{-3}
			B_integrand[1] *= FM3_TO_KM3;

			B_integrand.Interpolate(0, 1);
		}

		// --------------------------------------------------------
		// 1.d) Fill the sequence point from the profile
		// --------------------------------------------------------
		auto &seq = prof_.SeqMutable();
		// energy density at center
		if (prof_.HasColumn(StarProfile::Column::EnergyDensity))
		{
			const auto &eps0 = prof_.Radial()[prof_.GetColumnIndex(StarProfile::Column::EnergyDensity)][0];
			seq.ec = eps0 * Zaki::Physics::INV_FM4_2_G_CM3 /
					 Zaki::Physics::INV_FM4_2_INV_KM2; // g/cm^3
		}
		else
		{
			seq.ec = 0.0;
		}

		// mass and radius at surface
		seq.r = prof_.GetRadius()->operator[](-1); // km

		if (prof_.HasColumn(StarProfile::Column::Mass))
		{
			const auto &m_last = prof_.GetMass()->operator[](-1);
			seq.m = m_last / Zaki::Physics::SUN_M_KM; // M_sun
		}
		else
		{
			seq.m = 0.0; // M_sun
		}

		// central pressure
		if (prof_.HasColumn(StarProfile::Column::Pressure))
		{
			const auto &p0 = prof_.GetPressure()->operator[](0);
			seq.pc = p0 * Zaki::Physics::INV_FM4_2_Dyn_CM2 /
					 Zaki::Physics::INV_FM4_2_INV_KM2; // dyne/cm^2
		}
		else
		{
			seq.pc = 0.0;
		}

		// total baryon number (visible) from integrand
		if (B_integrand[0].Size() > 0)
		{
			seq.b = B_integrand.Integrate(
				1, {B_integrand[0][0], B_integrand[0][-1]});
		}
		else
		{
			seq.b = 0.0;
		}

		// moment of inertia
		seq.I = Find_MomInertia();

		surface_ready = true;
		return;
	}
}

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
	// auto &radial = prof_.radial;
	auto &radial = prof_.RadialMutable();

	const int idx_r = prof_.GetColumnIndex(StarProfile::Column::Radius);
	const int idx_m = prof_.GetColumnIndex(StarProfile::Column::Mass);
	const int idx_nup = prof_.GetColumnIndex(StarProfile::Column::MetricNuPrime);
	const int idx_p = prof_.GetColumnIndex(StarProfile::Column::Pressure);
	const int idx_eps = prof_.GetColumnIndex(StarProfile::Column::EnergyDensity);
	const int idx_nb = prof_.GetColumnIndex(StarProfile::Column::BaryonDensity);
	const int idx_nu = prof_.GetColumnIndex(StarProfile::Column::MetricNu);
	const int idx_lam = prof_.GetColumnIndex(StarProfile::Column::MetricLambda);

	// Now actually append the values (same unit conversions as legacy)
	radial[idx_r].PushBack(in_tov.r);
	radial[idx_m].PushBack(Zaki::Physics::SUN_M_KM * in_tov.m); // solar-mass → km
	radial[idx_nup].PushBack(in_tov.nu_der * 1e5);				// 1/cm→1/km
	radial[idx_p].PushBack(
		in_tov.p * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_Dyn_CM2);
	radial[idx_eps].PushBack(
		in_tov.e * Zaki::Physics::INV_FM4_2_INV_KM2 / Zaki::Physics::INV_FM4_2_G_CM3);
	radial[idx_nb].PushBack(in_tov.rho);
	// we do not set nu here — EvaluateNu() will overwrite this column later
	radial[idx_nu].PushBack(0.0);

	// ------------------------------------------------------------
	// Compute and append λ to the profile
	// ------------------------------------------------------------
	// Compute the Schwarzschild-like factor:
	const double r_km = in_tov.r;
	const double m_km = Zaki::Physics::SUN_M_KM * in_tov.m;

	// ADR-0004: the metric factor has ONE mathematical owner,
	// CompactStar/Geometry.hpp. This Path-1 site no longer defines it.
	//
	// g_rr = 1/(1 - 2m/r);  lambda = 0.5*ln(g_rr) = -0.5*ln(1 - 2m/r).
	//
	// Bit-identical to the former inline form on every ordinary node: the primitive
	// performs the same `1.0 - 2.0 * m_km / r_km` subtraction and the same
	// `-0.5 * std::log(denom)`. The former `denom <= 0 -> 1e-15` clamp is gone; under
	// the accepted hybrid physical-domain contract (ADR-0004 §0-Q3) such input fails
	// closed rather than being regularized. It was unreachable here in any case:
	// measured max 2m/r = 0.481 across the authenticated stars (ADR-0004 §11).
	//
	// Deferred out of Phase 3D because Path 1 had no coverage; Phase 3E-0 supplied it,
	// so this migrated in Phase 3E-I2.
	const double lambda_geom = CompactStar::Geometry::Lambda(r_km, m_km);

	radial[idx_lam].PushBack(lambda_geom);
	// ------------------------------------------------------------

	// ------------------------------------------------------------
	// 2.a) per-species for the profile
	// ------------------------------------------------------------
	// now append the data
	// for (std::size_t k = 0; k < in_tov.rho_i.size(); ++k)
	// {
	// 	const int col_idx = prof_.species_idx[k];
	// 	radial[col_idx].PushBack(in_tov.rho_i[k]); // fm^{-3}
	// }

	const std::size_t n_species = prof_.SpeciesCount();
	for (std::size_t k = 0; k < n_species; ++k)
	{
		const double val = (k < in_tov.rho_i.size()) ? in_tov.rho_i[k] : 0.0;
		radial[prof_.SpeciesColumnIndex(k)].PushBack(val);
	}

	// ------------------------------------------------------------
	// 2.b) update surface guess (sequence) incrementally
	// ------------------------------------------------------------
	// We can’t finalize here (that’s FinalizeSurface’s job) but we can
	// at least update the running “last point” so if we inspect sequence
	// before finalization, we get something reasonable.
	auto &seq = prof_.SeqMutable();
	seq.r = in_tov.r; // km
	seq.m = in_tov.m; // M_sun  ← in solar masses
					  // pc, ec, B, I will be filled / fixed at finalize time
}

//-------------------------------------------------------------- */
// Sets the work directory for the member objects
void NStar::OnWorkDirChanged(
	const Zaki::String::Directory &in_dir)
{
	// std::cout << "[NStar::OnWorkDirChanged] Setting work directory to: "
	// 		  << in_dir << "\n";
	prof_.RadialMutable().SetWrkDir(in_dir);
	// std::cout << "[NStar::OnWorkDirChanged] prof_.radial.GetWrkDir(): "
	// 		  << prof_.radial.GetWrkDir() << "\n";
}

//--------------------------------------------------------------
/// Similar to the destructor
void NStar::Reset()
{
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
/// Evaluate the metric function
void NStar::EvaluateNu()
{
	PROFILE_FUNCTION();

	// prefer StarProfile path
	if (!prof_.empty() &&
		prof_.HasColumn(StarProfile::Column::Radius) &&
		prof_.HasColumn(StarProfile::Column::Mass) &&
		prof_.HasColumn(StarProfile::Column::MetricNuPrime))
	{
		// auto &radial = prof_.radial;
		auto &radial = prof_.RadialMutable();
		const int rcol = prof_.GetColumnIndex(StarProfile::Column::Radius);
		const int mcol = prof_.GetColumnIndex(StarProfile::Column::Mass);
		const int nupcol = prof_.GetColumnIndex(StarProfile::Column::MetricNuPrime);

		// ensure ν column exists (create or clear)
		int nucol = prof_.GetColumnIndex(StarProfile::Column::MetricNu);
		if (!prof_.IsValidColumnIndex(nucol))
		{
			// append a new column for ν(r)
			// nucol = static_cast<int>(radial.Dim().size());
			// radial.AddColumn("nu");
			// prof_.SetColumnIndex(StarProfile::Column::MetricNu, nucol);
			Z_LOG_ERROR("MetricNu column is not present; cannot EvaluateNu(). "
						"Ensure profile layout includes nu during initialization.");
			return;
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

		return prof_.Radial().Evaluate(
			prof_.GetColumnIndex(StarProfile::Column::MetricNu),
			in_r);
	}

	return 0.0;
}

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
			const auto &radial = prof_.Radial();
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

		return 0.0;
	}
}

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

		return prof_.Radial().Evaluate(
			prof_.GetColumnIndex(StarProfile::Column::BaryonDensity),
			in_r);
	}

	return 0.0;
}

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

		return prof_.Radial().Evaluate(
			prof_.GetColumnIndex(StarProfile::Column::EnergyDensity),
			in_r);
	}

	return 0.0;
}

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

		return prof_.Radial().Evaluate(
			prof_.GetColumnIndex(StarProfile::Column::Pressure),
			in_r);
	}

	return 0.0;
}

// ------------------------------------------------------------
// Sequence forwarding
const SeqPoint &NStar::GetSequence() const noexcept
{
	if (!prof_.empty())
		return prof_.Seq();

	// fallback: static empty
	static SeqPoint empty{};
	return empty;
}
// ------------------------------------------------------------
// Sequence forwarding (non-const)
SeqPoint &NStar::GetSequence() noexcept
{
	if (!prof_.empty())
		return prof_.SeqMutable();

	// fallback: static singleton we can write to (rare path)
	static SeqPoint empty{};
	return empty;
}

//--------------------------------------------------------------
// Baryon number integrand
double NStar::BaryonNumIntegrand(double in_r) const
{
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
void NStar::Export(const Zaki::String::Directory &rel_path)
{
	PROFILE_FUNCTION();

	if (prof_.empty())
	{
		Z_LOG_ERROR("StarProfile is empty; nothing to export.");
		return;
	}

	// 1) Ensure the profile's DataSet uses the same work directory as NStar.
	prof_.RadialMutable().SetWrkDir(wrk_dir_);

	// 2) Build the full path and force-create its parent directory
	Zaki::String::Directory full_path = wrk_dir_ + rel_path;
	Zaki::String::Directory parent_dir = full_path.ParentDir();

	// 4) Delegate actual TSV export via StarProfile
	prof_.Export(rel_path);
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
