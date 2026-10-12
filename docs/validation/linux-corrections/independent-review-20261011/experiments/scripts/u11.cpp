#include <cstdio>
#include <cstddef>
#include <gsl/gsl_const_cgsm.h>
// Production conversion (CompactStar/RelativityUnits.hpp): pure double, compile-time constants.
constexpr double Gd = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT, cd = GSL_CONST_CGSM_SPEED_OF_LIGHT;
constexpr double TovC2 = cd * cd;
constexpr double EnergyKmMinus2PerGramCmMinus3 = Gd / TovC2 * 1e10;
constexpr double rho_b_km2 = 1e10 * EnergyKmMinus2PerGramCmMinus3;
// Fixture as written in tests/relativity/relativistic_unit_boundary.cpp (x86_64: 80-bit long double)
void x86() {
  constexpr long double G = 6.673e-8L, c = 29979245800.L;
  double col[3]; for (std::size_t i = 0; i < 3; ++i) col[i] = static_cast<double>((1.0001 - i * 0.0001) * 1e10L * G / c / c * 1e10L);
  std::size_t sel = 3; for (std::size_t i = 3; i-- > 0;) if (col[i] >= rho_b_km2) { sel = i; break; }
  std::printf("x86_64 80-bit: col1=%a threshold=%a  col1>=thr? %d  selected=%zu  col1 inverse=%.17g g/cm3\n", col[1], rho_b_km2, col[1] >= rho_b_km2, sel, col[1] / EnergyKmMinus2PerGramCmMinus3);
}
// Same fixture with long double == double (arm64 semantics); literal 6.673e-8L parses as double there.
void arm64() {
  constexpr double G = 6.673e-8, c = 29979245800.;
  double col[3]; for (std::size_t i = 0; i < 3; ++i) col[i] = static_cast<double>((1.0001 - i * 0.0001) * 1e10 * G / c / c * 1e10);
  std::size_t sel = 3; for (std::size_t i = 3; i-- > 0;) if (col[i] >= rho_b_km2) { sel = i; break; }
  std::printf("arm64-equiv   : col1=%a threshold=%a  col1>=thr? %d  selected=%zu  col1 inverse=%.17g g/cm3\n", col[1], rho_b_km2, col[1] >= rho_b_km2, sel, col[1] / EnergyKmMinus2PerGramCmMinus3);
}
int main() { std::printf("(1.0001-1*0.0001) as double = %a (exactly 1? %d)\n", 1.0001 - 1 * 0.0001, (1.0001 - 1 * 0.0001) == 1.0); x86(); arm64(); return 0; }
