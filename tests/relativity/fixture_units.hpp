// TEST SIDE ONLY: independently written physical adapters for analytic geometric stars.
// Do not include production RelativityUnits here: the exact Schwarzschild and Hartle
// oracles must retain independence. Declared decimal GSL 2.7.1 CGS literals (ADR-0012).
#pragma once
namespace relativity_fixture
{
inline constexpr double c = 29979245800.;
inline constexpr double G = 6.673e-8;
inline constexpr double solar_grams = 1.98892e33;
inline constexpr double solar_km = G * solar_grams / (c * c) / 1e5;
inline constexpr double rho_to_eps = G / (c * c) * 1e10;
inline constexpr double pressure_to_geo = G / (c * c * c * c) * 1e10;
inline constexpr double eps_to_rho = (c * c) / G / 1e10;
inline constexpr double pressure_to_cgs = (c * c * c * c) / G / 1e10;
}
