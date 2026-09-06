#include "tests/relativity/fixture_units.hpp"
// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file eos_table_knots.hpp
 * @brief TEST-ONLY reader of a CompactStar EOS table file (`e(g/cm^3)  p(dyne/cm^2)  rho ...`)
 *        into geometric (p, eps) knots [km^-2], for the Phase 4D-RV EOS-knot-refined measure
 *        partition of `hartle_mono_ref::SolveStieltjes`.
 *
 * The knots are INPUT to the independent oracle (the star's own EOS values), never an oracle
 * themselves. Nothing here touches production's importer or interpolant; the conversion uses
 * the same two Zaki constants the profile uses for eps and p (INV-02).
 */

#ifndef CompactStar_Tests_EosTableKnots_H
#define CompactStar_Tests_EosTableKnots_H

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <Zaki/Physics/Constants.hpp>

namespace eos_knots
{

struct Knots
{
	std::vector<double> p, eps; // km^-2, in file order (increasing pressure for a valid table)
	bool ok = false;
	std::size_t rows = 0;
};

/// Read `path`; the first line is a header, every further line has at least two numeric
/// columns e[g/cm^3], p[dyne/cm^2]. Non-numeric lines are skipped.
inline Knots Read(const std::string &path)
{
	Knots k;
	std::ifstream in(path);
	if (!in)
		return k;
	const double e_to_km2 = relativity_fixture::rho_to_eps;
	const double p_to_km2 = relativity_fixture::pressure_to_geo;
	std::string line;
	bool first = true;
	while (std::getline(in, line))
	{
		if (first)
		{
			first = false;
			continue;
		}
		std::istringstream ss(line);
		double e = 0, p = 0;
		if (!(ss >> e >> p))
			continue;
		k.eps.push_back(e * e_to_km2);
		k.p.push_back(p * p_to_km2);
		++k.rows;
	}
	k.ok = k.rows >= 2;
	return k;
}

} // namespace eos_knots

#endif /* CompactStar_Tests_EosTableKnots_H */
