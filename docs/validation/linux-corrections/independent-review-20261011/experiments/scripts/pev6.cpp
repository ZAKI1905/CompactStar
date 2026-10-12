#include <cmath>
#include <cstdio>
#include <limits>
// Exact Zaki double constants (sources/Zaki/Zaki/Physics/src/Constants.cpp)
static const double MN = 939.56542052, MP = 938.27208816, ME = 0.51099895000, MEV2INVFM = 5.067730716156395e-3;
// Production path: CompactStar/EOS/src/TrackRFreeGasThermodynamics.cpp:219 + LocalThermodynamics.cpp:198
template<class LD> double production() {
  const double mu_e = ((MN - MP) * (MN + MP) + ME * ME) / (2.0 * MN);
  const double hbar_c_member = 1.0 / MEV2INVFM;
  const LD chemical_potential = mu_e; const LD mass = ME; const LD hbar_c = hbar_c_member;
  const LD p_F = std::sqrt((chemical_potential - mass) * (chemical_potential + mass));
  const LD pi = std::acos((LD)-1.0);
  const LD ratio = p_F / hbar_c;
  return static_cast<double>(ratio * ratio * ratio / ((LD)3.0 * pi * pi));
}
// Independent path: tests/eos/rotochemical_trackr_pe.cpp:112
template<class LD> LD independent() {
  const LD neutron_mass = MN, proton_mass = MP, electron_mass = ME;
  const LD hbar_c = (LD)1.0 / static_cast<LD>(MEV2INVFM);
  const LD mu_e = ((neutron_mass - proton_mass) * (neutron_mass + proton_mass) + electron_mass * electron_mass) / ((LD)2.0 * neutron_mass);
  const LD momentum = std::sqrt((mu_e - electron_mass) * (mu_e + electron_mass));
  const LD high_pi = std::acos((LD)-1.0);
  return momentum * momentum * momentum / ((LD)3.0 * high_pi * high_pi * hbar_c * hbar_c * hbar_c);
}
int main() {
  const long double reviewed = 7.3567289037328326656352e-9L;
  std::printf("long double digits=%d\n", std::numeric_limits<long double>::digits);
  double pL = production<long double>(); long double iL = independent<long double>();
  double pD = production<double>();      double iD = independent<double>();
  std::printf("x86_64 (80-bit long double):\n  production=%.25e\n  independent=%.25Le\n  |prod/indep-1|=%.3Le  |indep-rev|/rev=%.3Le\n", pL, iL, fabsl((long double)pL/iL-1.0L), fabsl((iL-reviewed)/reviewed));
  std::printf("arm64-equivalent (long double == double):\n  production=%.25e\n  independent=%.25e\n  |prod/indep-1|=%.3Le  |indep-rev|/rev=%.3Le\n", pD, iD, fabsl((long double)pD/(long double)iD-1.0L), fabsl(((long double)iD-reviewed)/reviewed));
  std::printf("hex: prod_x86=%a prod_arm=%a indep_arm=%a\n", pL, pD, iD);
  return 0;
}
