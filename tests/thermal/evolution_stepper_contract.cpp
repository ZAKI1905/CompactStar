// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 * Copyright (c) 2026 Mohammadreza Zakeri
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file evolution_stepper_contract.cpp
 * @brief Guards the evolution stepper contract and thermal-only execution.
 *
 * Self-contained: synthetic star, trivial analytic driver, NO external EOS data.
 *
 * Regressions this catches:
 *   S1 — the default stepper must be a runnable explicit method. A return to
 *        MSBDF-by-default would previously reach GSL with sys.jacobian == nullptr and
 *        SIGSEGV, because gsl_odeiv2.h:69 dereferences (S)->jacobian with no null guard.
 *   S2 — explicitly requesting MSBDF must fail CLEANLY with a diagnostic, never crash and
 *        never be silently substituted by another method.
 *   S3 — a Thermal-only StateVector must execute. EvolutionSystem::operator() used to bind
 *        GetSpin() in a dead debug block, making a registered Spin block mandatory.
 */

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Physics/Driver/IDriver.hpp"
#include "CompactStar/Physics/Evolution/DriverContext.hpp"
#include "CompactStar/Physics/Evolution/EvolutionConfig.hpp"
#include "CompactStar/Physics/Evolution/EvolutionSystem.hpp"
#include "CompactStar/Physics/Evolution/GeometryCache.hpp"
#include "CompactStar/Physics/Evolution/Integrator/GSLIntegrator.hpp"
#include "CompactStar/Physics/Evolution/Run/RunBuilder.hpp"
#include "CompactStar/Physics/Evolution/StarContext.hpp"
#include "CompactStar/Physics/Evolution/StatePacking.hpp"
#include "CompactStar/Physics/State/ThermalState.hpp"

namespace P = CompactStar::Physics;
using CompactStar::Core::StarProfile;

static int g_fail = 0;
static void Report(const std::string &id, bool ok, const std::string &d)
{
	std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << id << " — " << d << "\n";
	if (!ok)
		++g_fail;
}

/// dY/dt = -k Y on the thermal slot; exact solution Y(t) = Y0 exp(-k t).
class DecayDriver final : public P::IDriver
{
  public:
	explicit DecayDriver(double k) : k_(k) {}
	std::string Name() const override { return "DecayDriver"; }
	const std::vector<P::State::StateTag> &DependsOn() const override
	{
		static const std::vector<P::State::StateTag> d{P::State::StateTag::Thermal};
		return d;
	}
	const std::vector<P::State::StateTag> &Updates() const override { return DependsOn(); }
	void AccumulateRHS(double, const P::Evolution::StateVector &Y,
					   P::Evolution::RHSAccumulator &dYdt,
					   const P::Evolution::DriverContext &) const override
	{
		const auto &th = Y.GetThermal();
		if (th.Size() == 0)
			return;
		dYdt.AddTo(P::State::StateTag::Thermal, 0, -k_ * th.LnTinfOverTref());
	}

  private:
	double k_;
};

static void FillProfile(StarProfile &prof, std::size_t N)
{
	auto edit = prof.Edit();
	auto &rad = prof.RadialMutable();
	rad.ClearRows();
	rad.Reserve(8, N);
	const char *labels[] = {"r(km)", "m(km)", "nu_prime(km^-1)", "p(km^-2)",
							"eps(km^-2)", "nB(fm^-3)", "nu", "lambda"};
	using C = StarProfile::Column;
	const C cols[] = {C::Radius, C::Mass, C::MetricNuPrime, C::Pressure,
					  C::EnergyDensity, C::BaryonDensity, C::MetricNu, C::MetricLambda};
	for (int i = 0; i < 8; ++i)
	{
		rad[i].SetLabel(labels[i]);
		prof.SetColumnIndex(cols[i], i);
	}
	const double h = 10.0 / double(N - 1);
	for (std::size_t i = 0; i < N; ++i)
	{
		rad[0].PushBack(double(i) * h);
		for (int c = 1; c <= 4; ++c)
			rad[c].PushBack(0.0);
		rad[5].PushBack(0.3);
		rad[6].PushBack(0.0);
		rad[7].PushBack(0.0);
	}
}

int main()
{
	std::cout << "Evolution stepper contract (synthetic; no external data)\n\n";

	StarProfile prof;
	FillProfile(prof, 65);
	P::Evolution::StarContext starCtx(prof);
	P::Evolution::GeometryCache geo(starCtx);

	const double k = 0.5;
	const double t0 = 0.0, t1 = 2.0;

	auto run = [&](P::Evolution::StepperType st, double &y_out, std::string &err) -> bool {
		P::Evolution::Config cfg;
		cfg.stepper = st;
		cfg.rtol = 1e-9;
		cfg.atol = 1e-12;
		cfg.max_internal_steps = 100000;
		cfg.max_samples = 100000;
		cfg.save_cadence = P::Evolution::SaveCadence::LinearDt;
		cfg.dt_save = 0.1;

		P::Evolution::DriverContext ctx;
		ctx.star = &starCtx;
		ctx.geo = &geo;
		ctx.cfg = &cfg;

		P::State::ThermalState th;
		th.Resize(1);
		th.SetTinf(P::State::ThermalState::Tref_K()); // ln(T/Tref) = 0 -> use 1.0 below
		th.LnTinfOverTref() = 1.0;

		P::Evolution::Run::StateWiring w;
		// THERMAL ONLY — no Spin block registered anywhere.
		w.state_vec.Register(P::State::StateTag::Thermal, th);
		std::vector<P::State::StateTag> tags{P::State::StateTag::Thermal};
		P::Evolution::Run::ConfigureLayout(w, tags);
		P::Evolution::Run::ConfigureRHS(w, tags);

		std::vector<P::Evolution::DriverPtr> drivers;
		drivers.push_back(std::make_shared<DecayDriver>(k));
		P::Evolution::EvolutionSystem sys(ctx, w.state_vec, w.rhs, w.layout,
										  std::move(drivers));
		P::Evolution::GSLIntegrator integ(sys, cfg, w.dim);

		std::vector<double> y(w.dim);
		P::Evolution::PackStateVector(w.state_vec, w.layout, y.data());
		try
		{
			if (!integ.Integrate(t0, t1, y.data()))
			{
				err = "Integrate returned false";
				return false;
			}
		}
		catch (const std::exception &e)
		{
			err = e.what();
			return false;
		}
		y_out = y[0];
		return true;
	};

	// -----------------------------------------------------------------
	std::cout << "S1/S3  default stepper runs a THERMAL-ONLY system\n";
	{
		P::Evolution::Config probe; // default-constructed: the production default
		double y = 0;
		std::string err;
		const bool ok = run(probe.stepper, y, err);
		Report("S1.a default stepper completes without crashing", ok, ok ? "ok" : err);
		if (ok)
		{
			const double exact = std::exp(-k * (t1 - t0));
			const double rel = std::fabs(y - exact) / exact;
			std::cout << "      y(2) = " << y << "  exact = " << exact
					  << "  rel = " << rel << "\n";
			Report("S1.b matches the analytic solution", rel < 1e-6,
				   "rel err " + std::to_string(rel));
		}
		Report("S3 Thermal-only StateVector executes (no Spin registered)", ok,
			   ok ? "no 'Spin is not registered' throw" : err);
	}

	// -----------------------------------------------------------------
	std::cout << "\nS2  explicit MSBDF fails cleanly, never crashes, never substitutes\n";
	{
		double y = 0;
		std::string err;
		const bool ok = run(P::Evolution::StepperType::MSBDF, y, err);
		Report("S2.a MSBDF is rejected rather than run", !ok, err);
		const bool mentions =
			err.find("MSBDF") != std::string::npos && err.find("Jacobian") != std::string::npos;
		Report("S2.b diagnostic names MSBDF and the missing Jacobian", mentions, err);
	}

	// -----------------------------------------------------------------
	std::cout << "\nS4  explicit steppers agree on the analytic solution\n";
	{
		const double exact = std::exp(-k * (t1 - t0));
		for (auto st : {P::Evolution::StepperType::RKF45, P::Evolution::StepperType::RKCK,
						P::Evolution::StepperType::RK8PD})
		{
			double y = 0;
			std::string err;
			const bool ok = run(st, y, err);
			const double rel = ok ? std::fabs(y - exact) / exact : 1.0;
			Report("S4 explicit stepper matches analytic", ok && rel < 1e-6,
				   ok ? ("rel err " + std::to_string(rel)) : err);
		}
	}

	std::cout << "\n"
			  << (g_fail == 0 ? "stepper contract holds"
							  : "CONTRACT FAILURES: " + std::to_string(g_fail))
			  << "\n";
	return g_fail == 0 ? 0 : 1;
}
