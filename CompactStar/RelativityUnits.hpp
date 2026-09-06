// Ordinary canonical TOV unit boundary — accepted ADR-0012 A1.
// Scientific candidate: independent review and human ratification precede baseline supersession.
#ifndef CompactStar_RelativityUnits_H
#define CompactStar_RelativityUnits_H

#include <gsl/gsl_const_cgsm.h>

namespace CompactStar::RelativityUnits
{
// These aliases authenticate the numerical solve convention, not modern metrology.
inline constexpr double TovGravitationalConstantCgs = GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT;
inline constexpr double LiteralGslSolarMassGrams = GSL_CONST_CGSM_SOLAR_MASS;
namespace detail
{
// Conversion implementation only. AngularVelocity remains the physical/geometric spin owner.
inline constexpr double TovC2 = GSL_CONST_CGSM_SPEED_OF_LIGHT * GSL_CONST_CGSM_SPEED_OF_LIGHT;
inline constexpr double MassKmPerGram = TovGravitationalConstantCgs / TovC2 / 1e5;
inline constexpr double EnergyKmMinus2PerGramCmMinus3 = TovGravitationalConstantCgs / TovC2 * 1e10;
inline constexpr double PressureKmMinus2PerDynCmMinus2 = TovGravitationalConstantCgs / (TovC2 * TovC2) * 1e10;
}
/// Physical mass [g] <-> geometric mass length [km] for the SAME GSL-solved spacetime.
inline constexpr double MassGramsToKm(double grams) { return grams * detail::MassKmPerGram; }
inline constexpr double MassKmToGrams(double km) { return km / detail::MassKmPerGram; }
/// Mass-equivalent energy density [g cm^-3] <-> geometric energy density [km^-2].
inline constexpr double MassDensityGcm3ToEnergyKmMinus2(double rho) { return rho * detail::EnergyKmMinus2PerGramCmMinus3; }
inline constexpr double EnergyKmMinus2ToMassDensityGcm3(double eps) { return eps / detail::EnergyKmMinus2PerGramCmMinus3; }
/// Physical pressure [dyn cm^-2] <-> geometric pressure [km^-2]; requires c^4.
inline constexpr double PressureDynCm2ToKmMinus2(double pressure) { return pressure * detail::PressureKmMinus2PerDynCmMinus2; }
inline constexpr double PressureKmMinus2ToDynCm2(double pressure) { return pressure / detail::PressureKmMinus2PerDynCmMinus2; }
/// Radial metric derivative nu' [cm^-1] <-> [km^-1]; exact length prefix.
inline constexpr double NuPrimeCmInverseToKmInverse(double derivative) { return derivative * 1e5; }
inline constexpr double NuPrimeKmInverseToCmInverse(double derivative) { return derivative / 1e5; }
/// Public mass is the literal GSL solar-mass ratio, NOT an IAU nominal GM_sun ratio.
inline constexpr double MassGramsToLiteralSolarMass(double grams) { return grams / LiteralGslSolarMassGrams; }
inline constexpr double LiteralSolarMassToGrams(double ratio) { return ratio * LiteralGslSolarMassGrams; }
inline constexpr double LiteralSolarMassToMassKm(double ratio) { return MassGramsToKm(LiteralSolarMassToGrams(ratio)); }
inline constexpr double MassKmToLiteralSolarMass(double km) { return MassGramsToLiteralSolarMass(MassKmToGrams(km)); }
} // namespace CompactStar::RelativityUnits
#endif
