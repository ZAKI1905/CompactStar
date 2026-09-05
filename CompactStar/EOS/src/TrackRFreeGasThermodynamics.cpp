#include "CompactStar/EOS/TrackRFreeGasThermodynamics.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

#include <Zaki/Physics/Constants.hpp>

namespace CompactStar
{
namespace
{

template <typename Function>
double BisectDecreasingRoot(double lower, double upper, Function &&function,
							double residual_tolerance, const char *context)
{
	double lower_value = function(lower);
	double upper_value = function(upper);
	if (!std::isfinite(lower_value) || !std::isfinite(upper_value) ||
		lower_value < 0.0 || upper_value > 0.0)
		throw std::runtime_error(std::string(context) +
								 " does not have a finite decreasing-sign bracket.");

	for (int iteration = 0; iteration < 256; ++iteration)
	{
		const double midpoint = lower + 0.5 * (upper - lower);
		const double midpoint_value = function(midpoint);
		if (!std::isfinite(midpoint_value))
			throw std::runtime_error(std::string(context) +
									 " produced a non-finite residual.");
		if (std::abs(midpoint_value) <= residual_tolerance)
			return midpoint;
		if (midpoint_value > 0.0)
		{
			lower = midpoint;
			lower_value = midpoint_value;
		}
		else
		{
			upper = midpoint;
			upper_value = midpoint_value;
		}
		const double width_tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
			std::max({1.0, std::abs(lower), std::abs(upper)});
		if (upper - lower <= width_tolerance)
			return lower + 0.5 * (upper - lower);
	}
	throw std::runtime_error(std::string(context) +
							 " did not converge within 256 bisections.");
}

void RequireStrictlyPositive(double value, const char *name)
{
	if (!std::isfinite(value) || value <= 0.0)
		throw std::runtime_error(std::string(name) +
							 " must be finite and strictly positive in the active npe-mu domain.");
}

} // namespace

TrackRFreeGasThermodynamicProvider::TrackRFreeGasThermodynamicProvider()
	: neutron_(ColdRelativisticIdealFermion::Neutron()),
	  proton_(ColdRelativisticIdealFermion::Proton()),
	  electron_(ColdRelativisticIdealFermion::Electron()),
	  muon_(ColdRelativisticIdealFermion::Muon()),
	  muon_onset_n_B_fm3_(0.0), sigma_minus_onset_n_B_fm3_(0.0)
{
	muon_onset_n_B_fm3_ = ComputeMuonOnsetBaryonDensityFm3();
	sigma_minus_onset_n_B_fm3_ = ComputeSigmaMinusOnsetBaryonDensityFm3();

	std::ostringstream domain;
	domain << std::setprecision(17)
		   << "Evaluate: n_n,n_p,n_e,n_mu>0 and n_B<"
		   << sigma_minus_onset_n_B_fm3_
		   << " fm^-3; EquilibriumStateAt: " << muon_onset_n_B_fm3_
		   << "<n_B<" << sigma_minus_onset_n_B_fm3_
		   << " fm^-3 (strict active npe-mu branch; no threshold smoothing)";
	metadata_ = {
		"track-r-fernandez-reisenegger-2005-free-gas-local",
		"FR2005 noninteracting-Fermi-gas local model; CompactStar local-r1; "
		"R2006 corrected charge-neutral response interpretation",
		"neutrons, protons, electrons, muons (all cold ideal spin-1/2 fermions)",
		"x=(n_B,n_e,n_mu) fm^-3; n_p=n_e+n_mu; n_n=n_B-n_p",
		"T=0 only",
		"total species energy including Zaki neutron/proton/electron/muon rest masses",
		"provider owns analytic free electron and muon contributions exactly once",
		domain.str()};
}

const LocalThermodynamicProviderMetadata &
TrackRFreeGasThermodynamicProvider::Metadata() const noexcept
{
	return metadata_;
}

double TrackRFreeGasThermodynamicProvider::ChargeDensityAtCommonLeptonMuFm3(
	double common_lepton_mu_MeV) const
{
	if (!std::isfinite(common_lepton_mu_MeV) ||
		common_lepton_mu_MeV < electron_.RestMassEnergyMeV())
		throw std::runtime_error(
			"Common lepton chemical potential is outside the cold free-gas domain.");
	const double n_e = electron_.NumberDensityForChemicalPotentialFm3(
		common_lepton_mu_MeV);
	const double n_mu = common_lepton_mu_MeV > muon_.RestMassEnergyMeV()
		? muon_.NumberDensityForChemicalPotentialFm3(common_lepton_mu_MeV)
		: 0.0;
	const double charge_density = n_e + n_mu;
	if (!std::isfinite(charge_density))
		throw std::runtime_error("Lepton charge-density relation is non-finite.");
	return charge_density;
}

double TrackRFreeGasThermodynamicProvider::EquilibriumResidualAtCommonLeptonMu(
	double n_B_fm3, double common_lepton_mu_MeV) const
{
	const double n_p_fm3 = ChargeDensityAtCommonLeptonMuFm3(common_lepton_mu_MeV);
	if (!std::isfinite(n_B_fm3) || n_B_fm3 <= n_p_fm3)
		throw std::runtime_error(
			"Equilibrium common-mu bracket left positive-neutron charge-neutral matter.");
	const double n_n_fm3 = n_B_fm3 - n_p_fm3;
	return neutron_.ChemicalPotentialMeV(n_n_fm3) -
		   proton_.ChemicalPotentialMeV(n_p_fm3) - common_lepton_mu_MeV;
}

ChargeNeutralCompositionState
TrackRFreeGasThermodynamicProvider::SolveActiveEquilibriumUnchecked(double n_B_fm3) const
{
	if (!std::isfinite(n_B_fm3) || n_B_fm3 <= muon_onset_n_B_fm3_)
		throw std::runtime_error(
			"Active npe-mu equilibrium requires n_B strictly above muon onset.");

	const double lower = muon_.RestMassEnergyMeV();
	const double upper = neutron_.ChemicalPotentialMeV(n_B_fm3) -
		proton_.RestMassEnergyMeV();
	if (!(upper > lower) || ChargeDensityAtCommonLeptonMuFm3(upper) >= n_B_fm3)
		throw std::runtime_error(
			"Active npe-mu equilibrium could not construct its physical common-mu bracket.");
	const auto residual = [this, n_B_fm3](double common_mu)
	{
		return EquilibriumResidualAtCommonLeptonMu(n_B_fm3, common_mu);
	};
	const double common_mu = BisectDecreasingRoot(
		lower, upper, residual, 2.0e-12, "Track-R free-gas equilibrium solve");
	const double n_e_fm3 = electron_.NumberDensityForChemicalPotentialFm3(common_mu);
	const double n_mu_fm3 = muon_.NumberDensityForChemicalPotentialFm3(common_mu);
	const auto state = MakeChargeNeutralCompositionState(
		{n_B_fm3, n_e_fm3, n_mu_fm3});
	RequireStrictlyPositive(state.NeutronDensityFm3(), "n_n");
	RequireStrictlyPositive(state.ProtonDensityFm3(), "n_p");
	RequireStrictlyPositive(state.ElectronDensityFm3(), "n_e");
	RequireStrictlyPositive(state.MuonDensityFm3(), "n_mu");
	return state;
}

double TrackRFreeGasThermodynamicProvider::ComputeMuonOnsetBaryonDensityFm3() const
{
	const double common_mu = muon_.RestMassEnergyMeV();
	const double n_e_fm3 = electron_.NumberDensityForChemicalPotentialFm3(common_mu);
	const double mu_p_MeV = proton_.ChemicalPotentialMeV(n_e_fm3);
	const double n_n_fm3 = neutron_.NumberDensityForChemicalPotentialFm3(
		mu_p_MeV + common_mu);
	const double onset = n_n_fm3 + n_e_fm3;
	if (!std::isfinite(onset) || onset <= 0.0)
		throw std::runtime_error("Muon-onset construction left the finite physical domain.");
	return onset;
}

double TrackRFreeGasThermodynamicProvider::ComputeSigmaMinusOnsetBaryonDensityFm3() const
{
	// FR2005 eq. (71): eta_nnSigma p=2 mu_n-mu_Sigma-mu_p.  At
	// first appearance the zero-density Sigma-minus chemical potential is its
	// rest mass; no hyperon thermodynamics is introduced here.
	const double sigma_mass_MeV = Zaki::Physics::SIGMA_MINUS_M_MEV;
	const auto onset_residual = [this, sigma_mass_MeV](double n_B_fm3)
	{
		const auto state = SolveActiveEquilibriumUnchecked(n_B_fm3);
		return sigma_mass_MeV +
			   proton_.ChemicalPotentialMeV(state.ProtonDensityFm3()) -
			   2.0 * neutron_.ChemicalPotentialMeV(state.NeutronDensityFm3());
	};

	double lower = std::nextafter(muon_onset_n_B_fm3_,
								  std::numeric_limits<double>::infinity());
	double upper = 2.0 * lower;
	while (upper < 32.0 && onset_residual(upper) > 0.0)
		upper *= 2.0;
	if (upper >= 32.0 && onset_residual(upper) > 0.0)
		throw std::runtime_error(
			"Sigma-minus onset was not bracketed below 32 fm^-3.");
	return BisectDecreasingRoot(lower, upper, onset_residual, 2.0e-11,
								"Track-R Sigma-minus onset solve");
}

ChargeNeutralCompositionState
TrackRFreeGasThermodynamicProvider::EquilibriumStateAt(double n_B_fm3) const
{
	if (!std::isfinite(n_B_fm3) || n_B_fm3 <= muon_onset_n_B_fm3_ ||
		n_B_fm3 >= sigma_minus_onset_n_B_fm3_)
		throw std::runtime_error(
			"Track-R equilibrium request is outside the strict active npe-mu source domain.");
	return SolveActiveEquilibriumUnchecked(n_B_fm3);
}

std::array<double, 4>
TrackRFreeGasThermodynamicProvider::IntrinsicChemicalPotentialsMeV(
	const ChargeNeutralCoordinates &coordinates) const
{
	const auto state = MakeChargeNeutralCompositionState(coordinates);
	RequireStrictlyPositive(state.NeutronDensityFm3(), "n_n");
	RequireStrictlyPositive(state.ProtonDensityFm3(), "n_p");
	RequireStrictlyPositive(state.ElectronDensityFm3(), "n_e");
	RequireStrictlyPositive(state.MuonDensityFm3(), "n_mu");
	if (state.BaryonDensityFm3() >= sigma_minus_onset_n_B_fm3_)
		throw std::runtime_error(
			"Track-R free-gas npe-mu model ends before beta-equilibrated Sigma-minus onset.");
	return {neutron_.ChemicalPotentialMeV(state.NeutronDensityFm3()),
			proton_.ChemicalPotentialMeV(state.ProtonDensityFm3()),
			electron_.ChemicalPotentialMeV(state.ElectronDensityFm3()),
			muon_.ChemicalPotentialMeV(state.MuonDensityFm3())};
}

LocalThermodynamicEvaluation TrackRFreeGasThermodynamicProvider::Evaluate(
	const ChargeNeutralCoordinates &coordinates) const
{
	const auto state = MakeChargeNeutralCompositionState(coordinates);
	RequireStrictlyPositive(state.NeutronDensityFm3(), "n_n");
	RequireStrictlyPositive(state.ProtonDensityFm3(), "n_p");
	RequireStrictlyPositive(state.ElectronDensityFm3(), "n_e");
	RequireStrictlyPositive(state.MuonDensityFm3(), "n_mu");
	if (state.BaryonDensityFm3() >= sigma_minus_onset_n_B_fm3_)
		throw std::runtime_error(
			"Track-R free-gas npe-mu model ends before beta-equilibrated Sigma-minus onset.");

	const auto neutron = neutron_.Evaluate(state.NeutronDensityFm3());
	const auto proton = proton_.Evaluate(state.ProtonDensityFm3());
	const auto electron = electron_.Evaluate(state.ElectronDensityFm3());
	const auto muon = muon_.Evaluate(state.MuonDensityFm3());
	const double D_n = *neutron.dchemical_potential_dn_MeV_fm3;
	const double D_p = *proton.dchemical_potential_dn_MeV_fm3;
	const double D_e = *electron.dchemical_potential_dn_MeV_fm3;
	const double D_mu = *muon.dchemical_potential_dn_MeV_fm3;

	const double eta_npe = neutron.chemical_potential_MeV -
		proton.chemical_potential_MeV - electron.chemical_potential_MeV;
	const double eta_npmu = neutron.chemical_potential_MeV -
		proton.chemical_potential_MeV - muon.chemical_potential_MeV;
	const NeutralConjugates conjugates{{
		neutron.chemical_potential_MeV, -eta_npe, -eta_npmu}};
	const ChargeNeutralChemicalHessian::Matrix hessian_values{{
		{{D_n, -D_n, -D_n}},
		{{-D_n, D_n + D_p + D_e, D_n + D_p}},
		{{-D_n, D_n + D_p, D_n + D_p + D_mu}},
	}};
	const ChargeNeutralChemicalHessian hessian{hessian_values};
	const double energy_density = neutron.energy_density_MeV_fm3 +
		proton.energy_density_MeV_fm3 + electron.energy_density_MeV_fm3 +
		muon.energy_density_MeV_fm3;
	if (!std::isfinite(energy_density))
		throw std::runtime_error("Track-R free-gas neutral energy is non-finite.");

	return {state, energy_density, conjugates, hessian};
}

} // namespace CompactStar
