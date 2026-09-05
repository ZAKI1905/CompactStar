// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file TrackRFreeGasThermodynamics.hpp
 * @brief Local cold ideal npe-mu model used by the Track-R free-gas benchmark.
 *
 * This is only the source-model local thermodynamic potential.  It constructs
 * no star and contains no electrostatic projection, stellar integral, or
 * rotochemical-evolution coefficient.
 */

#ifndef CompactStar_TrackRFreeGasThermodynamics_H
#define CompactStar_TrackRFreeGasThermodynamics_H

#include "CompactStar/EOS/LocalThermodynamics.hpp"

#include <array>
#include <stdexcept>

namespace CompactStar
{

/// Physical equilibrium active-set inventory, including the unimplemented p-e gate.
enum class FreeGasEquilibriumDomain
{
	ProtonElectron,
	NeutronOnset,
	Npe,
	MuonOnset,
	NpeMuon
};

/// The requested physical branch exists, but finite precision cannot certify it.
/// No nearby density or different active set is substituted.
class EquilibriumResolutionError : public std::runtime_error
{
  public:
	using std::runtime_error::runtime_error;
};

/**
 * @brief Fernandez--Reisenegger Track-R cold noninteracting npe-mu provider.
 *
 * Evaluate() supports the smooth active-species domain n_n,n_p,n_e,n_mu>0
 * below the beta-equilibrated Sigma-minus onset.  It does not impose beta
 * equilibrium. EvaluateNpe() uses the separate z=(n_B,n_e) active chart.
 * EquilibriumAt() dispatches explicitly between npe, value-only muon onset,
 * and npe-mu. The lower p-e branch is identified but not implemented.
 * AI-authored scientific candidate; see GOVERNANCE.md section 5.
 */
class TrackRFreeGasThermodynamicProvider final : public ILocalThermodynamicProvider
{
  public:
	TrackRFreeGasThermodynamicProvider();

	[[nodiscard]] const LocalThermodynamicProviderMetadata &Metadata() const noexcept override;
	[[nodiscard]] LocalThermodynamicEvaluation
	Evaluate(const ChargeNeutralCoordinates &coordinates) const override;
	[[nodiscard]] ChargeNeutralCompositionState
	EquilibriumStateAt(double n_B_fm3) const override;

	[[nodiscard]] NpeThermodynamicEvaluation EvaluateNpe(const NpeCoordinates &coordinates) const;
	[[nodiscard]] ActiveLocalThermodynamicEvaluation
	EvaluateActive(const ChargeNeutralCoordinates &coordinates) const override;
	[[nodiscard]] ActiveLocalThermodynamicEvaluation EquilibriumAt(double n_B_fm3) const override;
	[[nodiscard]] FreeGasEquilibriumDomain EquilibriumDomainAt(double n_B_fm3) const;
	[[nodiscard]] double NeutronOnsetBaryonDensityFm3() const noexcept
	{
		return neutron_onset_n_B_fm3_;
	}

	[[nodiscard]] double MuonOnsetBaryonDensityFm3() const noexcept
	{
		return muon_onset_n_B_fm3_;
	}
	[[nodiscard]] double SigmaMinusOnsetBaryonDensityFm3() const noexcept
	{
		return sigma_minus_onset_n_B_fm3_;
	}

	/// Optional source-specific species values in order (n,p,e,mu), all in MeV.
	[[nodiscard]] std::array<double, 4>
	IntrinsicChemicalPotentialsMeV(const ChargeNeutralCoordinates &coordinates) const;

  private:
	[[nodiscard]] ChargeNeutralCompositionState SolveNpeEquilibrium(double n_B_fm3) const;
	[[nodiscard]] MuonThresholdEvaluation EvaluateMuonThreshold() const;
	[[nodiscard]] double ComputeNeutronOnsetBaryonDensityFm3() const;
	[[nodiscard]] ChargeNeutralCompositionState
	SolveActiveEquilibriumUnchecked(double n_B_fm3) const;
	[[nodiscard]] double EquilibriumResidualAtCommonLeptonMu(
		double n_B_fm3, double common_lepton_mu_MeV) const;
	[[nodiscard]] double ChargeDensityAtCommonLeptonMuFm3(
		double common_lepton_mu_MeV) const;
	[[nodiscard]] double ComputeMuonOnsetBaryonDensityFm3() const;
	[[nodiscard]] double ComputeSigmaMinusOnsetBaryonDensityFm3() const;

	ColdRelativisticIdealFermion neutron_;
	ColdRelativisticIdealFermion proton_;
	ColdRelativisticIdealFermion electron_;
	ColdRelativisticIdealFermion muon_;
	double neutron_onset_n_B_fm3_;
	double muon_onset_n_B_fm3_;
	double sigma_minus_onset_n_B_fm3_;
	LocalThermodynamicProviderMetadata metadata_;
};

} // namespace CompactStar

#endif /* CompactStar_TrackRFreeGasThermodynamics_H */
