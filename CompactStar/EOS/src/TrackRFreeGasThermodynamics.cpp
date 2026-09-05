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
							 " must be finite and strictly positive in the requested active domain.");
}

} // namespace

TrackRFreeGasThermodynamicProvider::TrackRFreeGasThermodynamicProvider()
	: neutron_(ColdRelativisticIdealFermion::Neutron()),
	  proton_(ColdRelativisticIdealFermion::Proton()),
	  electron_(ColdRelativisticIdealFermion::Electron()),
	  muon_(ColdRelativisticIdealFermion::Muon()),
	  neutron_onset_n_B_fm3_(0.0), muon_onset_n_B_fm3_(0.0),
	  sigma_minus_onset_n_B_fm3_(0.0)
{
	neutron_onset_n_B_fm3_ = ComputeNeutronOnsetBaryonDensityFm3();
	muon_onset_n_B_fm3_ = ComputeMuonOnsetBaryonDensityFm3();
	sigma_minus_onset_n_B_fm3_ = ComputeSigmaMinusOnsetBaryonDensityFm3();

	std::ostringstream domain;
	domain << std::setprecision(17)
		   << "Evaluate: n_n,n_p,n_e,n_mu>0 and n_B<"
		   << sigma_minus_onset_n_B_fm3_
		   << " fm^-3; EvaluateNpe: n_n,n_p,n_e>0, n_mu=0, mu_e<m_mu, "
		   << neutron_onset_n_B_fm3_ << "<n_B<" << muon_onset_n_B_fm3_
		   << " fm^-3; EquilibriumAt/EquilibriumStateAt: " << neutron_onset_n_B_fm3_
		   << "<n_B<" << sigma_minus_onset_n_B_fm3_
		   << " fm^-3 (npe / value-only muon onset / npe-mu; p-e branch unavailable; "
		   << "roundoff-unresolved interior solves fail closed; no threshold smoothing)";
	metadata_ = {
		"track-r-fernandez-reisenegger-2005-free-gas-local",
		"FR2005 noninteracting-Fermi-gas local model; CompactStar local-r2-npe; "
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
	// A small residual tolerance must not manufacture a positive muon density
	// when the onset side is unresolved in double precision. Classification is
	// by n_B, without a tolerance; only numerical availability is restricted.
	const double roundoff_bound = 64.0 * std::numeric_limits<double>::epsilon() *
		(neutron_.ChemicalPotentialMeV(n_B_fm3) + proton_.RestMassEnergyMeV() + lower);
	if (residual(lower) <= roundoff_bound)
		throw EquilibriumResolutionError(
			"npe-mu branch classified above onset, but its bracket is roundoff-unresolved.");
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
		const auto state = n_B_fm3 == muon_onset_n_B_fm3_
			? MakeChargeNeutralCompositionState({n_B_fm3,
				electron_.NumberDensityForChemicalPotentialFm3(muon_.RestMassEnergyMeV()), 0.0})
			: SolveActiveEquilibriumUnchecked(n_B_fm3);
		return sigma_mass_MeV +
			   proton_.ChemicalPotentialMeV(state.ProtonDensityFm3()) -
			   2.0 * neutron_.ChemicalPotentialMeV(state.NeutronDensityFm3());
	};

	// Use the exact value-only threshold as the outer bracket endpoint; never
	// ask the inner smooth solver for a one-ULP-above-onset Hessian/state.
	double lower = muon_onset_n_B_fm3_;
	double upper = 2.0 * lower;
	while (upper < 32.0 && onset_residual(upper) > 0.0)
		upper *= 2.0;
	if (upper >= 32.0 && onset_residual(upper) > 0.0)
		throw std::runtime_error(
			"Sigma-minus onset was not bracketed below 32 fm^-3.");
	return BisectDecreasingRoot(lower, upper, onset_residual, 2.0e-11,
								"Track-R Sigma-minus onset solve");
}

double TrackRFreeGasThermodynamicProvider::ComputeNeutronOnsetBaryonDensityFm3() const
{
	// At n_n=0, n_p=n_e=n_B and sqrt(m_p^2+p^2)+sqrt(m_e^2+p^2)=m_n.
	// Factor the mass-square difference to avoid unnecessary cancellation.
	const double mn = neutron_.RestMassEnergyMeV();
	const double mp = proton_.RestMassEnergyMeV();
	const double me = electron_.RestMassEnergyMeV();
	const double mu_e = ((mn - mp) * (mn + mp) + me * me) / (2.0 * mn);
	return electron_.NumberDensityForChemicalPotentialFm3(mu_e);
}

FreeGasEquilibriumDomain
TrackRFreeGasThermodynamicProvider::EquilibriumDomainAt(double n_B_fm3) const
{
	if (!std::isfinite(n_B_fm3) || n_B_fm3 <= 0.0 ||
		n_B_fm3 >= sigma_minus_onset_n_B_fm3_)
		throw std::runtime_error("Track-R equilibrium density is outside the finite source domain.");
	if (n_B_fm3 < neutron_onset_n_B_fm3_)
		return FreeGasEquilibriumDomain::ProtonElectron;
	if (n_B_fm3 == neutron_onset_n_B_fm3_)
		return FreeGasEquilibriumDomain::NeutronOnset;
	if (n_B_fm3 < muon_onset_n_B_fm3_)
		return FreeGasEquilibriumDomain::Npe;
	if (n_B_fm3 == muon_onset_n_B_fm3_)
		return FreeGasEquilibriumDomain::MuonOnset;
	return FreeGasEquilibriumDomain::NpeMuon;
}

ChargeNeutralCompositionState
TrackRFreeGasThermodynamicProvider::SolveNpeEquilibrium(double n_B_fm3) const
{
	if (EquilibriumDomainAt(n_B_fm3) != FreeGasEquilibriumDomain::Npe)
		throw std::runtime_error("Npe equilibrium requires the open neutron-to-muon interval.");
	// At fixed n_B, F(n_e)=mu_n(n_B-n_e)-mu_p(n_e)-mu_e(n_e).
	// F(0)>0 because m_n>m_p+m_e; F(n_B)<0 precisely above neutron onset.
	// F'=-D_n-D_p-D_e<0: [0,n_B] is a physical bracket with one root.
	const auto residual = [this, n_B_fm3](double n_e_fm3)
	{
		return neutron_.ChemicalPotentialMeV(n_B_fm3 - n_e_fm3) -
			proton_.ChemicalPotentialMeV(n_e_fm3) -
			electron_.ChemicalPotentialMeV(n_e_fm3);
	};
	const double roundoff_bound = 64.0 * std::numeric_limits<double>::epsilon() *
		(neutron_.ChemicalPotentialMeV(n_B_fm3) +
		 proton_.ChemicalPotentialMeV(n_B_fm3) + electron_.ChemicalPotentialMeV(n_B_fm3));
	const double muon_onset_n_e = electron_.NumberDensityForChemicalPotentialFm3(
		muon_.RestMassEnergyMeV());
	if (residual(n_B_fm3) >= -roundoff_bound ||
		(n_B_fm3 >= muon_onset_n_e && residual(muon_onset_n_e) >= -roundoff_bound))
		throw EquilibriumResolutionError(
			"npe branch classified in its open interval, but an onset bracket is roundoff-unresolved.");
	double lower = 0.0;
	double upper = n_B_fm3;
	if (!(residual(lower) > 0.0 && residual(upper) < 0.0))
		throw std::runtime_error("Npe equilibrium has no physical decreasing-sign bracket.");
	// Continue to adjacent representable densities, not an absolute density
	// tolerance (which would erase the low-density end of this branch).
	for (int iteration = 0; iteration < 256; ++iteration)
	{
		const double midpoint = lower + 0.5 * (upper - lower);
		if (midpoint == lower || midpoint == upper)
		{
			const double n_e = std::abs(residual(lower)) < std::abs(residual(upper))
				? lower : upper;
			const auto state = MakeChargeNeutralCompositionState({n_B_fm3, n_e, 0.0});
			if (!(state.NeutronDensityFm3() > 0.0 && n_e > 0.0 &&
				  electron_.ChemicalPotentialMeV(n_e) < muon_.RestMassEnergyMeV()) ||
				std::abs(residual(n_e)) > 5.0e-11)
				throw EquilibriumResolutionError("Npe root is not representable in its strict active domain.");
			return state;
		}
		const double value = residual(midpoint);
		if (!std::isfinite(value))
			throw std::runtime_error("Npe equilibrium produced a non-finite residual.");
		if (value > 0.0)
			lower = midpoint;
		else
			upper = midpoint;
	}
	throw std::runtime_error("Npe equilibrium did not converge within 256 bisections.");
}

ChargeNeutralCompositionState
TrackRFreeGasThermodynamicProvider::EquilibriumStateAt(double n_B_fm3) const
{
	switch (EquilibriumDomainAt(n_B_fm3))
	{
	case FreeGasEquilibriumDomain::Npe:
		return SolveNpeEquilibrium(n_B_fm3);
	case FreeGasEquilibriumDomain::MuonOnset:
		return EvaluateMuonThreshold().state;
	case FreeGasEquilibriumDomain::NpeMuon:
		return SolveActiveEquilibriumUnchecked(n_B_fm3);
	case FreeGasEquilibriumDomain::ProtonElectron:
	case FreeGasEquilibriumDomain::NeutronOnset:
		throw std::runtime_error("Track-R p-e branch and neutron-appearance boundary are a separate unimplemented gate.");
	}
	throw std::runtime_error("Unknown Track-R equilibrium domain.");
}

MuonThresholdEvaluation TrackRFreeGasThermodynamicProvider::EvaluateMuonThreshold() const
{
	const double n_e = electron_.NumberDensityForChemicalPotentialFm3(muon_.RestMassEnergyMeV());
	const auto state = MakeChargeNeutralCompositionState({muon_onset_n_B_fm3_, n_e, 0.0});
	const auto n = neutron_.Evaluate(state.NeutronDensityFm3());
	const auto p = proton_.Evaluate(n_e);
	const auto e = electron_.Evaluate(n_e);
	return {state, n.energy_density_MeV_fm3 + p.energy_density_MeV_fm3 + e.energy_density_MeV_fm3,
		{{n.chemical_potential_MeV, -n.chemical_potential_MeV + p.chemical_potential_MeV + e.chemical_potential_MeV}},
		n.chemical_potential_MeV - p.chemical_potential_MeV - muon_.RestMassEnergyMeV()};
}

NpeThermodynamicEvaluation TrackRFreeGasThermodynamicProvider::EvaluateNpe(
	const NpeCoordinates &coordinates) const
{
	const auto state = MakeChargeNeutralCompositionState({coordinates.n_B_fm3, coordinates.n_e_fm3, 0.0});
	RequireStrictlyPositive(state.NeutronDensityFm3(), "n_n");
	RequireStrictlyPositive(state.ProtonDensityFm3(), "n_p=n_e");
	if (EquilibriumDomainAt(coordinates.n_B_fm3) != FreeGasEquilibriumDomain::Npe ||
		electron_.ChemicalPotentialMeV(coordinates.n_e_fm3) >= muon_.RestMassEnergyMeV())
		throw std::runtime_error("Npe response requires the declared open sub-muon domain and mu_e<m_mu.");
	const auto n = neutron_.Evaluate(state.NeutronDensityFm3());
	const auto p = proton_.Evaluate(state.ProtonDensityFm3());
	const auto e = electron_.Evaluate(state.ElectronDensityFm3());
	const double Dn = *n.dchemical_potential_dn_MeV_fm3;
	const double Dp = *p.dchemical_potential_dn_MeV_fm3;
	const double De = *e.dchemical_potential_dn_MeV_fm3;
	// Restrict epsilon to n_mu=0 before differentiation: h=(mu_n,-eta_npe).
	const NpeChemicalHessian H{{{{Dn, -Dn}, {-Dn, Dn + Dp + De}}}};
	return {state, n.energy_density_MeV_fm3 + p.energy_density_MeV_fm3 + e.energy_density_MeV_fm3,
		{{n.chemical_potential_MeV, -n.chemical_potential_MeV + p.chemical_potential_MeV + e.chemical_potential_MeV}},
		H, n.chemical_potential_MeV - p.chemical_potential_MeV - muon_.RestMassEnergyMeV()};
}

ActiveLocalThermodynamicEvaluation TrackRFreeGasThermodynamicProvider::EvaluateActive(
	const ChargeNeutralCoordinates &coordinates) const
{
	// Validate first; no negative/NaN muon density may select an inactive path.
	const auto state = MakeChargeNeutralCompositionState(coordinates);
	if (state.MuonDensityFm3() > 0.0)
		return Evaluate(coordinates);
	if (coordinates.n_B_fm3 == muon_onset_n_B_fm3_ &&
		coordinates.n_e_fm3 == electron_.NumberDensityForChemicalPotentialFm3(muon_.RestMassEnergyMeV()))
		return EvaluateMuonThreshold();
	return EvaluateNpe({coordinates.n_B_fm3, coordinates.n_e_fm3});
}

ActiveLocalThermodynamicEvaluation TrackRFreeGasThermodynamicProvider::EquilibriumAt(double n_B_fm3) const
{
	return EvaluateActive(EquilibriumStateAt(n_B_fm3).Coordinates());
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
