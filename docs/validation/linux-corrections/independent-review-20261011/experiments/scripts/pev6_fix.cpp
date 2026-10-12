#include <cmath>
#include <cstdio>
static const double MN = 939.56542052, MP = 938.27208816, ME = 0.51099895000, MEV2INVFM = 5.067730716156395e-3;
// Variant A (current production): mu_e and hbar_c rounded to double, then long double section.
template<class LD> double prodA(bool fma_mu) {
  double mu_e = fma_mu ? std::fma(MN - MP, MN + MP, ME * ME) / (2.0 * MN) : ((MN - MP) * (MN + MP) + ME * ME) / (2.0 * MN);
  const double hbar_c_member = 1.0 / MEV2INVFM;
  const LD cp = mu_e, mass = ME, hbar_c = hbar_c_member;
  const LD p_F = std::sqrt((cp - mass) * (cp + mass)); const LD pi = std::acos((LD)-1.0); const LD ratio = p_F / hbar_c;
  return static_cast<double>(ratio * ratio * ratio / ((LD)3.0 * pi * pi));
}
// Variant B (proposed): mu_e computed in long double; hbar_c member still double.
template<class LD> double prodB() {
  const LD mn = MN, mp = MP, me = ME;
  const LD mu_e = ((mn - mp) * (mn + mp) + me * me) / ((LD)2.0 * mn);
  const double hbar_c_member = 1.0 / MEV2INVFM;
  const LD cp = mu_e, mass = ME, hbar_c = hbar_c_member;
  const LD p_F = std::sqrt((cp - mass) * (cp + mass)); const LD pi = std::acos((LD)-1.0); const LD ratio = p_F / hbar_c;
  return static_cast<double>(ratio * ratio * ratio / ((LD)3.0 * pi * pi));
}
// Variant C (proposed): mu_e and hbar_c both in long double.
template<class LD> double prodC() {
  const LD mn = MN, mp = MP, me = ME;
  const LD mu_e = ((mn - mp) * (mn + mp) + me * me) / ((LD)2.0 * mn);
  const LD hbar_c = (LD)1.0 / (LD)MEV2INVFM;
  const LD cp = mu_e, mass = ME;
  const LD p_F = std::sqrt((cp - mass) * (cp + mass)); const LD pi = std::acos((LD)-1.0); const LD ratio = p_F / hbar_c;
  return static_cast<double>(ratio * ratio * ratio / ((LD)3.0 * pi * pi));
}
int main() {
  const double truth_cr = 0x1.f98cf18cad246p-28;  // correctly rounded double of the exact value (mpmath)
  auto ulp_off = [&](double v) { return (long long)((v - truth_cr) / 0x1p-80); };
  std::printf("correctly rounded truth = %a\n", truth_cr);
  std::printf("x86 (80-bit LD):   A(no fma)=%a [%+lld ulp]  A(fma)=%a [%+lld ulp]  B=%a [%+lld ulp]  C=%a [%+lld ulp]\n",
    prodA<long double>(false), ulp_off(prodA<long double>(false)), prodA<long double>(true), ulp_off(prodA<long double>(true)),
    prodB<long double>(), ulp_off(prodB<long double>()), prodC<long double>(), ulp_off(prodC<long double>()));
  std::printf("arm64-equiv (LD=double): A(no fma)=%a [%+lld ulp]  A(fma)=%a [%+lld ulp]  B=%a [%+lld ulp]  C=%a [%+lld ulp]\n",
    prodA<double>(false), ulp_off(prodA<double>(false)), prodA<double>(true), ulp_off(prodA<double>(true)),
    prodB<double>(), ulp_off(prodB<double>()), prodC<double>(), ulp_off(prodC<double>()));
  return 0;
}
