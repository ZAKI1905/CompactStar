// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file hartle_thorne_1968_hw_eos.hpp
 * @brief TEST-ONLY historical fixture: the Harrison-Wheeler equation of state exactly as
 *        printed in Hartle & Thorne (1968), ApJ 153, 807, Table 1 (p. 814), plus the published
 *        model tables (Table 3, p. 816; Table 5, p. 818) used by Phase 4D Experiment K.
 *
 * PROVENANCE. Transcribed by eye from the NASA ADS journal scan of ApJ 153, 807, rendered at
 * 300 dpi and read column by column (Phase 4D, 2026-09-03). Table 1 gives P, E and eps all in
 * cm^-2 (G = c = 1); the footnote reads: "To convert to g cm^-3, divide by G c^-2 = 0.742 x
 * 10^-28 cm g^-1." Computer notation: 8.31E-41 = 8.31 x 10^-41.
 *
 * TRANSCRIPTION UNCERTAINTY. Every entry was read twice (page-level render, then a 300-dpi
 * crop of each column group). One entry is smudged in the scan: the eps value of the row
 * P = 5.19E-29 reads "?.82E-29" with an illegible leading digit; log-log interpolation between
 * its neighbours (3.03E-29 at E = 4.68E-24 and 2.27E-28 at E = 1.48E-23) gives 6.8E-29, so
 * 6.82E-29 is recorded. eps enters ONLY the baryon-density column (rho = (E - eps)/mu), which
 * neither the TOV integration nor the Hartle equations consume, so this uncertainty cannot
 * affect any compared quantity. P and E are three-significant-figure values; the models are
 * "accurate to about 1 per cent or better" (Table 3 footnote), which is why ADR-0007 §7 item 8
 * predeclares 2e-2 for this comparison.
 *
 * HT68 interpolated "logarithmically between table entries" (p. 813). CompactStar's importer
 * builds a Steffen spline in (p, eps); on a table with ~2.2 entries per decade that would be a
 * different EOS model. To reproduce THEIR model faithfully the table is densified here by
 * log-log linear interpolation between the printed rows before import; that is their
 * interpolant, not an improvement of it.
 *
 * NOT A PRODUCTION EOS. Nothing here is used outside the Phase-4D published-comparison harness.
 */

#ifndef CompactStar_Tests_HartleThorne1968HW_H
#define CompactStar_Tests_HartleThorne1968HW_H

#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace ht68
{

/// HT68 Table 1, all in cm^-2: {P, E, eps}. Left column group (36 rows) then right (36 rows).
struct Row
{
	double P, E, eps;
};

inline const std::vector<Row> &HWTable1()
{
	static const std::vector<Row> t = {
		// ---- left group, p. 814 ----
		{8.31e-41, 5.82e-28, 1.00e-45}, {4.17e-40, 5.84e-28, 7.11e-43}, {8.31e-40, 5.86e-28, 2.77e-42},
		{4.17e-39, 5.97e-28, 4.16e-41}, {8.31e-39, 6.04e-28, 1.12e-40}, {2.79e-38, 6.32e-28, 8.64e-40},
		{2.38e-37, 8.52e-28, 3.36e-38}, {1.37e-36, 1.23e-27, 3.21e-37}, {7.00e-36, 2.34e-27, 3.47e-36},
		{6.96e-35, 6.54e-27, 5.01e-35}, {4.79e-34, 1.56e-26, 3.76e-34}, {1.74e-33, 2.95e-26, 1.52e-33},
		{5.95e-33, 5.22e-26, 5.18e-33}, {1.56e-32, 8.52e-26, 1.45e-32}, {4.62e-32, 1.53e-25, 4.73e-32},
		{2.67e-31, 3.71e-25, 2.72e-31}, {9.63e-31, 7.41e-25, 1.05e-30}, {4.83e-30, 1.86e-24, 5.82e-30},
		{2.32e-29, 4.68e-24, 3.03e-29}, {5.19e-29, 7.41e-24, 6.82e-29} /* eps leading digit smudged */,
		{1.65e-28, 1.48e-23, 2.27e-28}, {8.23e-28, 3.71e-23, 1.11e-27}, {2.37e-27, 7.41e-23, 3.59e-27},
		{7.19e-27, 1.48e-22, 1.12e-26}, {2.69e-26, 3.71e-22, 4.84e-26}, {7.86e-26, 7.41e-22, 1.42e-25},
		{1.93e-25, 1.48e-21, 4.03e-25}, {6.60e-25, 3.71e-21, 1.53e-24}, {1.65e-24, 7.41e-21, 4.07e-24},
		{4.18e-24, 1.48e-20, 1.07e-23}, {1.35e-23, 3.71e-20, 3.77e-23}, {3.29e-23, 7.41e-20, 9.57e-23},
		{8.07e-23, 1.48e-19, 2.41e-22}, {2.67e-22, 3.71e-19, 8.17e-22}, {6.53e-22, 7.41e-19, 2.04e-21},
		{1.21e-21, 1.17e-18, 3.72e-21},
		// ---- right group, p. 814 ----
		{2.73e-21, 2.34e-18, 9.21e-21}, {6.49e-21, 4.68e-18, 2.25e-20}, {1.10e-20, 7.41e-18, 4.05e-20},
		{1.88e-20, 1.17e-17, 7.21e-20}, {3.05e-20, 1.86e-17, 1.28e-19}, {4.58e-20, 2.95e-17, 2.25e-19},
		{6.59e-20, 4.68e-17, 3.88e-19}, {9.55e-20, 7.41e-17, 6.60e-19}, {1.50e-19, 1.17e-16, 1.11e-18},
		{2.54e-19, 1.86e-16, 1.88e-18}, {4.49e-19, 2.95e-16, 3.17e-18}, {9.14e-19, 4.68e-16, 5.39e-18},
		{1.88e-18, 7.41e-16, 9.28e-18}, {6.09e-18, 1.48e-15, 2.18e-17}, {2.63e-17, 3.71e-15, 7.27e-17},
		{8.23e-17, 7.41e-15, 1.90e-16}, {2.60e-16, 1.48e-14, 5.16e-16}, {1.09e-15, 3.71e-14, 2.02e-15},
		{3.25e-15, 7.41e-14, 5.70e-15}, {9.71e-15, 1.48e-13, 1.61e-14}, {3.93e-14, 3.71e-13, 6.30e-14},
		{9.71e-14, 7.41e-13, 1.69e-13}, {2.42e-13, 1.48e-12, 4.34e-13}, {7.34e-13, 3.71e-12, 1.43e-12},
		{1.65e-12, 7.41e-12, 3.37e-12}, {3.60e-12, 1.48e-11, 7.71e-12}, {1.01e-11, 3.71e-11, 2.23e-11},
		{2.08e-11, 7.41e-11, 4.87e-11}, {4.28e-11, 1.48e-10, 1.04e-10}, {1.10e-10, 3.71e-10, 2.82e-10},
		{2.26e-10, 7.41e-10, 5.89e-10}, {4.64e-10, 1.48e-09, 1.22e-09}, {1.19e-09, 3.71e-09, 3.19e-09},
		{2.43e-09, 7.41e-09, 6.52e-09}, {4.91e-09, 1.48e-08, 1.33e-08}, {1.23e-08, 3.71e-08, 3.41e-08},
	};
	return t;
}

/// HT68 footnote constant: G c^-2 = 0.742e-28 cm g^-1. Used deliberately (their model), with the
/// exact defining c; the 3-digit constant contributes < 1e-3 to any compared quantity.
inline constexpr double kGoverC2_cm_per_g = 0.742e-28;
inline constexpr double kC_cm_per_s = 2.99792458e10;
/// Rest mass per baryon, the ^56Fe atom / 56 (HT68 Table 3 footnote), grams.
inline constexpr double kMuBaryon_g = 55.9349e0 * 1.66053907e-24 / 56.0;

/// HT68 Table 3 (non-rotating HWW configurations), p. 816: E_c [g/cm^3], R [km], M/Msun.
struct NonRot
{
	double ec, R_km, M_sun;
};
inline const std::vector<NonRot> &Table3()
{
	static const std::vector<NonRot> t = {
		{1.00e14, 3.60e1, 2.66e-1}, {3.00e14, 2.08e1, 4.05e-1}, {1.00e15, 1.42e1, 5.54e-1},
		{3.00e15, 1.02e1, 6.61e-1}, {6.00e15, 8.41e0, 6.84e-1}, {1.00e16, 7.48e0, 6.68e-1},
		{3.00e16, 5.96e0, 5.77e-1}, {1.00e17, 5.15e0, 4.62e-1},
	};
	return t;
}

/// HT68 Table 5 (slowly rotating HWW configurations at Omega^2 = M/R^3), p. 818:
/// E_c, R_g/R, Omega [s^-1], omega_s/Omega, omega_c/Omega, deltaR/R, deltaM/M.
struct Rot
{
	double ec, Rg_over_R, Omega_s, ws_over_W, wc_over_W, dR_over_R, dM_over_M;
};
inline const std::vector<Rot> &Table5()
{
	static const std::vector<Rot> t = {
		{1.00e14, 0.255, 8.67e2, 1.41e-3, 6.20e-2, 0.246, 0.0522}, {3.00e14, 0.355, 2.43e3, 7.19e-3, 1.18e-1, 0.201, 0.128},
		{1.00e15, 0.404, 5.07e3, 1.88e-2, 2.17e-1, 0.181, 0.163},  {3.00e15, 0.433, 9.10e3, 3.59e-2, 3.50e-1, 0.163, 0.162},
		{6.00e15, 0.444, 1.23e4, 4.73e-2, 4.43e-1, 0.154, 0.149},  {1.00e16, 0.444, 1.45e4, 5.20e-2, 5.03e-1, 0.152, 0.137},
		{3.00e16, 0.425, 1.90e4, 5.16e-2, 6.16e-1, 0.161, 0.107},  {1.00e17, 0.380, 2.12e4, 3.81e-2, 7.11e-1, 0.188, 0.0818},
	};
	return t;
}

/**
 * Write the HW EOS in CompactStar's import format (e[g/cm^3], p[dyne/cm^2], rho[1/fm^3]),
 * densified by log-log linear interpolation with @p per_decade points per decade of P — HT68's
 * own "logarithmic interpolation between table entries".
 */
inline bool WriteDensifiedEos(const std::string &path, int per_decade = 40)
{
	std::ofstream out(path);
	if (!out)
		return false;
	out << "e(g/cm^3)\tp(dyne/cm^2)\trho(1/fm^3)\n";
	out << std::scientific << std::setprecision(12);

	const auto &T = HWTable1();
	const double GoverC4 = kGoverC2_cm_per_g / (kC_cm_per_s * kC_cm_per_s); // cm^-2 per (dyne/cm^2)

	auto emit = [&](double P, double E, double eps) {
		const double e_cgs = E / kGoverC2_cm_per_g;
		const double p_cgs = P / GoverC4;
		const double n_fm3 = ((E - eps) / kGoverC2_cm_per_g) / kMuBaryon_g * 1.0e-39;
		out << e_cgs << "\t" << p_cgs << "\t" << n_fm3 << "\n";
	};

	for (std::size_t i = 0; i + 1 < T.size(); ++i)
	{
		const double lP0 = std::log10(T[i].P), lP1 = std::log10(T[i + 1].P);
		const double lE0 = std::log10(T[i].E), lE1 = std::log10(T[i + 1].E);
		const double lx0 = std::log10(T[i].eps), lx1 = std::log10(T[i + 1].eps);
		const int n = std::max(1, static_cast<int>(std::ceil((lP1 - lP0) * per_decade)));
		for (int k = 0; k < n; ++k)
		{
			const double t = static_cast<double>(k) / static_cast<double>(n);
			const double lP = lP0 + t * (lP1 - lP0);
			emit(std::pow(10.0, lP), std::pow(10.0, lE0 + t * (lE1 - lE0)),
				 std::pow(10.0, lx0 + t * (lx1 - lx0)));
		}
	}
	emit(T.back().P, T.back().E, T.back().eps);
	return true;
}

} // namespace ht68

#endif /* CompactStar_Tests_HartleThorne1968HW_H */
