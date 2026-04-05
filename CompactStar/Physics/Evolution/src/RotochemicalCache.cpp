/*
 * CompactStar
 * RotochemicalCache — implementation of rotochemical heating coefficients.
 *
 * Copyright (c) 2025 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

#include <cmath>
#include <stdexcept>

#include <Zaki/Physics/Constants.hpp>

#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/RotochemicalCache.hpp"

using namespace CompactStar::Physics::Evolution;

// Conversion factor: fm^{-3} -> km^{-3}
static constexpr double FM3_TO_KM3 = 1.0e54;

//--------------------------------------------------------------
double RotochemicalCache::ComputeEnclosedNumber(
	const Zaki::Vector::DataColumn &n_i,
	const Zaki::Vector::DataColumn &r_col,
	const GeometryCache &geo)
{
	const auto &wV = geo.WV(); // 4*pi*r^2*exp(Lambda)
	const size_t N = n_i.Size();

	// Trapezoidal integration: integral of n_i * wV * dr
	double integral = 0.0;
	for (size_t i = 1; i < N; i++)
	{
		double dr = r_col[i] - r_col[i - 1];
		double f_prev = n_i[i - 1] * wV[i - 1];
		double f_curr = n_i[i] * wV[i];
		integral += 0.5 * (f_prev + f_curr) * dr;
	}

	return integral * FM3_TO_KM3;
}

//--------------------------------------------------------------
double RotochemicalCache::ComputeStructuralDerivative(
	const Zaki::Vector::DataColumn &n_i,
	const Zaki::Vector::DataColumn &r_col,
	const Zaki::Vector::DataColumn &p_col,
	const CompactStar::Core::HartleResult &hartle,
	const GeometryCache &geo)
{
	// A_i = integral [ (dn_i/dp) * p0(r) * wV(r) ] dr
	//
	// where dn_i/dp is the derivative of the species density with respect
	// to pressure, computed numerically from the profile.
	//
	// This is the leading-order term. Additional metric correction terms
	// (from the change in the volume element at O(Omega^2)) should be
	// included for full accuracy but are subdominant.
	//
	// Ref: Reisenegger (1995), Fernandez & Reisenegger (2005) Eq. (A.14)

	const auto &wV = geo.WV();
	const auto &p0 = hartle.p0;
	const size_t N = n_i.Size();

	if (N < 3)
		return 0.0;

	// Compute dn_i/dp via centered finite differences
	Zaki::Vector::DataColumn dnidp("dnidp", N, 0.0);
	for (size_t i = 1; i < N - 1; i++)
	{
		double dp = p_col[i + 1] - p_col[i - 1];
		if (std::abs(dp) > 1.e-30)
			dnidp[i] = (n_i[i + 1] - n_i[i - 1]) / dp;
	}
	// Boundary: one-sided
	{
		double dp = p_col[1] - p_col[0];
		if (std::abs(dp) > 1.e-30)
			dnidp[0] = (n_i[1] - n_i[0]) / dp;
	}
	{
		size_t last = N - 1;
		double dp = p_col[last] - p_col[last - 1];
		if (std::abs(dp) > 1.e-30)
			dnidp[last] = (n_i[last] - n_i[last - 1]) / dp;
	}

	// Trapezoidal integration: integral of (dn_i/dp) * p0 * wV * dr
	double integral = 0.0;
	for (size_t i = 1; i < N; i++)
	{
		double dr = r_col[i] - r_col[i - 1];
		double f_prev = dnidp[i - 1] * p0[i - 1] * wV[i - 1];
		double f_curr = dnidp[i] * p0[i] * wV[i];
		integral += 0.5 * (f_prev + f_curr) * dr;
	}

	return integral * FM3_TO_KM3;
}

//--------------------------------------------------------------
void RotochemicalCache::Build(
	const CompactStar::Core::NStar &star,
	const CompactStar::Core::HartleResult &hartle,
	const GeometryCache &geo,
	const CompactStar::Core::NStar &star_plus,
	const CompactStar::Core::NStar &star_minus,
	double delta_eps_c)
{
	const auto &prof = star.Profile();
	const auto *r_col = prof.GetRadius();
	const auto *p_col = prof.GetPressure();
	const auto *nb_col = prof.GetBaryonDensity();

	if (!r_col || !p_col || !nb_col)
		return;

	// --- Total baryon number: A_B ---
	A_B = ComputeStructuralDerivative(*nb_col, *r_col, *p_col, hartle, geo);

	// --- Total baryon number: B_B (sequence derivative via finite difference) ---
	const auto *r_plus = star_plus.Profile().GetRadius();
	const auto *nb_plus = star_plus.Profile().GetBaryonDensity();
	const auto *r_minus = star_minus.Profile().GetRadius();
	const auto *nb_minus = star_minus.Profile().GetBaryonDensity();

	if (!r_plus || !nb_plus || !r_minus || !nb_minus)
		return;

	double N_B_plus = ComputeEnclosedNumber(*nb_plus, *r_plus, geo);
	double N_B_minus = ComputeEnclosedNumber(*nb_minus, *r_minus, geo);
	B_B = (N_B_plus - N_B_minus) / (2.0 * delta_eps_c);

	// --- Per-species coefficients ---
	species.clear();

	// Get species labels from profile
	const auto &labels = prof.species_labels;

	for (const auto &label : labels)
	{
		const auto *n_i = prof.GetSpeciesPtr(label);
		if (!n_i)
			continue;

		RotochemicalSpeciesData sd;
		sd.label = label;
		sd.N_total = ComputeEnclosedNumber(*n_i, *r_col, geo);
		sd.A = ComputeStructuralDerivative(*n_i, *r_col, *p_col, hartle, geo);

		// B_i from perturbed stars
		const auto *n_i_plus = star_plus.Profile().GetSpeciesPtr(label);
		const auto *n_i_minus = star_minus.Profile().GetSpeciesPtr(label);

		if (n_i_plus && n_i_minus)
		{
			double N_plus = ComputeEnclosedNumber(*n_i_plus, *r_plus, geo);
			double N_minus = ComputeEnclosedNumber(*n_i_minus, *r_minus, geo);
			sd.B = (N_plus - N_minus) / (2.0 * delta_eps_c);
		}

		// Z_i = A_i - B_i * (A_B / B_B)
		if (std::abs(B_B) > 1.e-30)
			sd.Z = sd.A - sd.B * (A_B / B_B);
		else
			sd.Z = sd.A;

		species.push_back(std::move(sd));
	}

	valid = true;
}
