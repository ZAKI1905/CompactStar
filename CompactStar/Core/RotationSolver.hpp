// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2023 Mohammadreza Zakeri
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/**
 * @file RotationSolver.hpp
 *
 * @brief Solves Hartle's equation of rotating neutron star.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@eku.edu
 *
 */
#ifndef CompactStar_Core_RotationSolver_H
#define CompactStar_Core_RotationSolver_H

#include <cstdint>
#include <gsl/gsl_spline.h>
#include <vector>

#include <Zaki/Math/Math_Core.hpp>
#include <Zaki/Vector/DataColumn.hpp>

#include "CompactStar/AngularVelocity.hpp"
#include "CompactStar/Core/Prog.hpp"

//==============================================================
namespace CompactStar::Core
{
class NStar;
class MixedStar;

//==============================================================
struct TOVTable
{
  public:
	std::vector<double> r;
	std::vector<double> m;
	std::vector<double> pre;
	std::vector<double> eps;

	size_t Size()
	{
		return r.size();
	}
};
//==============================================================
// Added on Aug 6, 2020
// omega_bar_c: The initial value for bar{omega} at r = 0
// Note that bar{omega}(r) = Omega - omega(r)
//
//  m     :   mass of the star
//  r     :   radius
//  J     :   total angular momentum of the star
//  Omega :   angular velocity of the star
//
struct OmegaSeqPoint
{
	double omega_bar_c, m, r, J, Omega;

	OmegaSeqPoint(const double &in_oc, const double &in_m,
				  const double &in_r, const double &in_J, const double &in_Omega)
		: omega_bar_c(in_oc), m(in_m), r(in_r), J(in_J), Omega(in_Omega)
	{
	}

	std::string Str() const
	{
		std::stringstream ss;
		char tmp[200];
		snprintf(tmp, sizeof(tmp), "%.8e\t %.8e\t %.8e\t %.8e\t %.8e", omega_bar_c, m, r, J, Omega);
		ss << tmp;

		return ss.str();
	}
};
//==============================================================
/// First-order Hartle solution at an explicitly requested physical angular velocity.
///
/// Produced ONLY by `HartleFirstOrderResponse::At(AngularVelocity)`. A star that was never
/// given a physical spin has no object of this type: ADR-0006's binding clarification is that
/// construction does not confer an implicit `Omega`.
///
/// **Storage is geometric, exactly once** (ADR-0006 Q3/P6). `Omega_geom` is the single stored
/// angular velocity; the physical value is *derived* by `OmegaRadPerSecond()` so the two can
/// never drift apart. There is no `Omega_phys` member and there must never be one.
struct PhysicalFirstOrderRotation
{
	/// Angular velocity as seen at infinity, geometric [km^-1]. CANONICAL — the one stored form.
	double Omega_geom = 0.0;

	/// Total angular momentum [km^2], from the exterior matching `omega_bar = Omega - 2J/r^3`
	/// applied to the scaled surface derivative: `J = R^4 omega_bar'(R) / 6`.
	double J = 0.0;

	/// Moment of inertia [km^3]. Scale-free and independent of `Omega`; carried through from the
	/// response unchanged, never recomputed as `J/Omega` (which is 0/0 at zero spin).
	double I = 0.0;

	/// omega_bar(r) = Omega - omega(r) at the requested spin, geometric [km^-1].
	/// Multiply an entry by `c` — i.e. `AngularVelocityGeomToRadPerSecond` — for s^-1.
	Zaki::Vector::DataColumn omega_bar;

	/// d(omega_bar)/dr at the requested spin, geometric [km^-2].
	Zaki::Vector::DataColumn domega_bar;

	/// Non-owning pointer to the radius column [km] this solution is tabulated on.
	/// Lifetime follows the owning `StarProfile`, per the ADR-0003 lifetime rule: this object
	/// must not outlive the star it was materialized from.
	const Zaki::Vector::DataColumn *r_grid = nullptr;

	bool valid = false; ///< True when the underlying first-order response was usable.

	/// The requested angular velocity in physical rad s^-1, derived from `Omega_geom` through
	/// the single geometric -> physical conversion owner (ADR-0006 P2/P6).
	[[nodiscard]] double OmegaRadPerSecond() const noexcept
	{
		return AngularVelocityGeomToRadPerSecond(Omega_geom);
	}
};
//==============================================================
/// The **seed-free** first-order rotational response of a star.
///
/// Every member is a quantity from which the arbitrary internal `omega_bar` normalization
/// cancels analytically, so nothing here depends on the numerical seed (ADR-0006 Q2/Q4):
///
///   - `I = J_raw / Omega_raw` is a ratio of two quantities that scale identically;
///   - `omega_bar_over_Omega(r) = omega_bar_raw(r) / Omega_raw` likewise;
///   - `domega_bar_over_Omega_dr(r) = domega_bar_raw(r) / Omega_raw` likewise.
///
/// This is what ordinary star construction leaves behind. It is a property of the *star*, not
/// of any particular spin — the frame-dragging shape `omega(r)/Omega = 1 - omega_bar/Omega`
/// and the moment of inertia are the same at every angular velocity in the linear theory.
///
/// To obtain a physical solution, supply an explicit `AngularVelocity` to `At()`.
struct HartleFirstOrderResponse
{
	/// Moment of inertia [km^3]. Scale-free; validated as such in Phase 2B-4B.
	double I = 0.0;

	/// omega_bar(r) / Omega — dimensionless.
	Zaki::Vector::DataColumn omega_bar_over_Omega;

	/// [d omega_bar/dr](r) / Omega — [km^-1].
	Zaki::Vector::DataColumn domega_bar_over_Omega_dr;

	/// Non-owning pointer to the radius column [km]. See the lifetime note above.
	const Zaki::Vector::DataColumn *r_grid = nullptr;

	bool valid = false; ///< True after a complete, usable first-order solve.

	/// Materialize the first-order solution at an explicitly requested physical angular
	/// velocity. This is a **scaling of an existing solution, not a new ODE solve**: the
	/// first-order problem is linear and homogeneous, so the response fully determines the
	/// solution at every `Omega` (ADR-0006 P3).
	///
	/// Zero spin is well defined and is not a special case: it yields `Omega = 0`, `J = 0`,
	/// `omega_bar == 0` and the unchanged `I` (ADR-0006 P5).
	///
	/// @throws std::runtime_error if the response is not valid.
	[[nodiscard]] PhysicalFirstOrderRotation At(AngularVelocity omega) const;
};
//==============================================================
/// The O(Omega^2) monopole perturbation of a star at an **explicitly requested** physical
/// angular velocity.
///
/// Produced ONLY by `HartleMonopoleResponse::At(AngularVelocity)`. Computing the normalized
/// response confers no spin: ADR-0006's binding clarification, extended to second order by
/// ADR-0007 P9, is that a physical perturbation exists only after a caller supplies an
/// `AngularVelocity`.
///
/// Every field is the corresponding coefficient times `Omega_geom^2`, so **no ODE is solved
/// here** and materializing many spins costs nothing. Because the quantities are quadratic in
/// `Omega`, `+Omega` and `-Omega` produce bit-identical perturbations, and zero spin produces
/// exact zeros.
struct PhysicalHartleMonopole
{
	/// Angular velocity as seen at infinity, geometric [km^-1]. CANONICAL — the one stored
	/// form; the physical value is derived, never stored alongside (ADR-0006 Q3/P6).
	double Omega_geom = 0.0;

	/// m0(r) — the O(Omega^2) perturbation of the enclosed-mass function [km].
	Zaki::Vector::DataColumn m0;

	/// p0*(r) — Hartle's dimensionless pressure perturbation factor (H67 eqs. 87-88, 99).
	Zaki::Vector::DataColumn p0star;

	/// delta_p0(r) = (eps + p) p0*(r) — the Eulerian monopole pressure perturbation [km^-2].
	Zaki::Vector::DataColumn delta_p0;

	/// xi0(r) = p0*/nu' — the outward displacement of the isobar that sits at r [km].
	Zaki::Vector::DataColumn xi0;

	/// delta_M — the O(Omega^2) change in the star's gravitational mass [km],
	/// `m0(R_*) + 4 pi R_*^2 eps_* xi0(R_*) + J^2/R_*^3` (ADR-0007 P6).
	double delta_M = 0.0;

	/// The surface mass-shell term of `delta_M` alone [km], kept separately because it is
	/// negligible on an EOS-floor star and dominant on a constant-density one.
	double surface_shell_mass = 0.0;

	/// xi0 at the production surface node R_* [km]. **R_* is the EOS-floor surface, not the
	/// exact p = 0 radius** (INV-06; ADR-0007 P7), so this carries that documented systematic.
	double surface_xi0 = 0.0;

	/// Non-owning pointer to the radius column [km] this solution is tabulated on. It must not
	/// outlive the owning `StarProfile`; `source_profile`/`source_version` below let a holder
	/// detect that the source has changed rather than read through a stale grid.
	const Zaki::Vector::DataColumn *r_grid = nullptr;

	/// Provenance of the response this was materialized from (ADR-0003).
	const void *source_profile = nullptr;
	std::uint64_t source_version = 0;

	bool valid = false;

	/// The requested angular velocity in physical rad s^-1, derived through the single
	/// geometric -> physical conversion owner (ADR-0006 P2/P6).
	[[nodiscard]] double OmegaRadPerSecond() const noexcept
	{
		return AngularVelocityGeomToRadPerSecond(Omega_geom);
	}
};
//==============================================================
/// The **seed-free, fixed-central-energy-density** O(Omega^2) monopole (l = 0) structural
/// response of a star, per unit `Omega_geom^2`.
///
/// **Governed by ADR-0007 (ACCEPTED 2026-09-02).** This replaced the retired `HartleResult`
/// candidate of commit `675b4a9`; nothing of that candidate's numerical content survives, and
/// none of its historical output is a reference result (ADR-0007 §9 item 6).
///
/// Every field carries the `_over_Omega2` suffix because that is exactly what it is: a
/// coefficient, not a perturbation. Multiplying by `Omega_geom^2` — which is what `At()` does —
/// is the only way to obtain a physical quantity. The coefficients are independent of `Omega`
/// (the system is linear with sources quadratic in `omega_bar`) and independent of the internal
/// first-order seed, because they are built from the verified seed-free `omega_bar/Omega` and
/// `[d omega_bar/dr]/Omega` and never from the raw stored profile.
///
/// **Fixed epsilon_c.** This is Hartle's particular solution with `m0(0) = p0*(0) = 0` and no
/// homogeneous admixture (H67 p. 1009, p. 1022; HT68 §II f): the rotating star has the same
/// central energy density as the non-rotating one. There is no surface condition, no shooting
/// and no free constant. The regular homogeneous solution — the non-rotating sequence
/// derivative — is deliberately **not** exposed (ADR-0007 P11 as modified at acceptance).
///
/// **EOS source (ADR-0008, ACCEPTED 2026-09-03).** The energy-density contribution to `m0` is
/// the measure `-4 pi r^2 xi0_hat d(eps)` of Hartle's eq. (93), integrated one governed profile
/// segment at a time, with the surface shell the terminal `eps_* -> 0` atom of that same
/// measure. The nodal `d(eps)/dp` column is retained for the regular-centre series and for
/// diagnostics, and is no longer the radial mass source: it cannot represent energy-density
/// variation narrower than the profile spacing, which is why Phase 4D's validation failed
/// (`docs/validation/PHASE4D_MONOPOLE_VALIDATION.md`).
///
/// **Not yet independently validated.** ADR-0007 (as amended by ADR-0008) fixes the contract and
/// this is its conforming implementation; the corrected independent revalidation is still
/// outstanding. No number here is validated physics and no baseline of it may be created
/// (INV-08).
struct HartleMonopoleResponse
{
	/// m0(r)/Omega_geom^2 [km^3].
	Zaki::Vector::DataColumn m0_over_Omega2;

	/// p0*(r)/Omega_geom^2 [km^2]. `p0*` itself is dimensionless.
	Zaki::Vector::DataColumn p0star_over_Omega2;

	/// delta_p0(r)/Omega_geom^2 — dimensionless, `= (eps + p) * p0star_over_Omega2`.
	Zaki::Vector::DataColumn delta_p0_over_Omega2;

	/// xi0(r)/Omega_geom^2 [km^3], `= p0star_over_Omega2 / nu'`.
	Zaki::Vector::DataColumn xi0_over_Omega2;

	/// delta_M/Omega_geom^2 [km^3] — `m0_hat(R_*) + shell_hat + I^2/R_*^3` (ADR-0007 P6).
	double deltaM_over_Omega2 = 0.0;

	/// The surface mass-shell term alone, `4 pi R_*^2 eps_* xi0_hat(R_*)` [km^3].
	double surface_shell_mass_over_Omega2 = 0.0;

	/// xi0_hat at the production surface node [km^3]. See `PhysicalHartleMonopole::surface_xi0`
	/// for the R_* semantics.
	double surface_xi0_over_Omega2 = 0.0;

	/// The production surface radius R_* [km] — the **last profile node**, i.e. the EOS-floor
	/// surface, deliberately not identified with the exact p = 0 radius (INV-06, ADR-0007 P7).
	double R_surface = 0.0;

	/// The moment of inertia [km^3] this response used for the exterior `I^2/R_*^3` term,
	/// carried through from the verified first-order response.
	double I = 0.0;

	/// Non-owning pointer to the radius column [km]. See the lifetime note on
	/// `PhysicalHartleMonopole::r_grid`.
	const Zaki::Vector::DataColumn *r_grid = nullptr;

	/// Provenance (ADR-0003): the identity and `Version()` of the `StarProfile` this response
	/// was computed from. A response whose provenance no longer matches its star is **stale**
	/// and must not be returned as current.
	const void *source_profile = nullptr;
	std::uint64_t source_version = 0;

	bool valid = false; ///< True only after a complete, all-finite solve (published atomically).

	/// True when this response was computed from `profile` at its current version.
	[[nodiscard]] bool MatchesSource(const void *profile, std::uint64_t version) const noexcept
	{
		return valid && source_profile == profile && source_version == version;
	}

	/// Materialize the monopole perturbation at an explicitly requested physical angular
	/// velocity: `Q = Q_hat * Omega_geom^2`. **No ODE is solved** (ADR-0007 P9/§19.3).
	///
	/// Zero spin yields exact zeros; `+Omega` and `-Omega` yield bit-identical results.
	///
	/// @throws std::runtime_error if the response is not valid.
	[[nodiscard]] PhysicalHartleMonopole At(AngularVelocity omega) const;
};
//==============================================================
/// Test-only access seam for `RotationSolver`'s internal numerical state.
///
/// Declared — never defined — in production. Validation harnesses define it to reach the
/// arbitrary `omega_bar` seed and the first-order solve counter, neither of which may become
/// public scientific API (ADR-0006 Q2). Friendship, not a public setter, is what keeps the seed
/// out of the contract while still allowing seed invariance to be *proved* rather than asserted.
struct RotationSolverTestSeam;
//==============================================================
class RotationSolver : public Prog
{
	//--------------------------------------------------------------
	// This struct was added on Aug 6
	struct OmegaPoint
	{
		double r, m, nu, omega_bar, domega_bar, omega;

		OmegaPoint(const double &in_r, const double &in_m,
				   const double &in_nu,
				   const double &in_omega_bar,
				   const double &in_domega_bar)
			: r(in_r), m(in_m), nu(in_nu), omega_bar(in_omega_bar),
			  domega_bar(in_domega_bar), omega(0)
		{
		}

		std::string Str() const
		{
			std::stringstream ss;
			char tmp[150];
			snprintf(tmp, sizeof(tmp), "%-19.8e\t %-19.8e\t %-19.8e\t %-19.8e\t %-19.8e\t %-19.8e", r, m, nu, omega, omega_bar, domega_bar);
			ss << tmp;

			return ss.str();
		}
	};
	//--------------------------------------------------------------
	// This struct was added on Apr 20, 2022
	//  to accommodate the mixed star scenarios
	struct OmegaPointDark
	{
		double r, m, m_d, nu, omega_bar, domega_bar, omega;

		OmegaPointDark(const double &in_r,
					   const double &in_m,
					   const double &in_m_d,
					   const double &in_nu,
					   const double &in_omega_bar,
					   const double &in_domega_bar)
			: r(in_r), m(in_m), m_d(in_m_d), nu(in_nu),
			  omega_bar(in_omega_bar),
			  domega_bar(in_domega_bar), omega(0)
		{
		}

		std::string Str() const
		{
			std::stringstream ss;
			char tmp[150];
			snprintf(tmp, sizeof(tmp), "%-19.8e\t %-19.8e\t %-19.8e\t "
									   "%-19.8e\t %-19.8e\t %-19.8e\t %-19.8e",
					 r, m, m_d, nu, omega, omega_bar, domega_bar);
			ss << tmp;

			return ss.str();
		}
	};

	//--------------------------------------------------------------
  private:
	NStar *nstar_ptr = nullptr;
	MixedStar *mixedstar_ptr = nullptr;

	size_t radial_res = 1000;

	// This workspace stores state variables for interpolation lookups.
	// It caches the previous value of an index lookup.
	// When the subsequent interpolation point falls in the same
	// interval its index value can be returned immediately.
	// gsl_interp_accel *accel = nullptr ;

	// We use cubic spline
	// Cubic spline with natural boundary conditions.
	// The resulting curve is piecewise cubic on each interval,
	// with matching first and second derivatives at the supplied data-points.
	//  The second derivative is chosen to be zero at the first point and last point.

	/// Mass spline
	// gsl_spline *m_spline = nullptr ;
	/// pressure spline
	// gsl_spline *p_spline = nullptr ;
	/// energy denisty spline
	// gsl_spline *e_spline = nullptr ;

	/// Initial bar{omega} (at r = r(i=0)), as used by the most recent solve.
	/// **Arbitrary internal numerical normalization — NOT a physical quantity** (ADR-0006 Q2).
	double init_omega_bar = -1;

	/// The default internal seed for the first-order solve. Its value carries no physical
	/// meaning: the frame-dragging equation is linear and homogeneous, so every seed yields the
	/// same physics once normalized. It is kept at its historical value because changing it
	/// would move `I` at the level of the solver's fixed absolute tolerance for no scientific
	/// gain (`docs/validation/PHASE4_ROTATION_ENTRY.md` §12).
	static constexpr double kDefaultInitOmegaBar = 5e-3;

	/// The seed the next first-order solve will use. Private, with no public setter: it is
	/// reachable only through `RotationSolverTestSeam`, so that seed invariance can be tested
	/// without the seed becoming part of the scientific API (ADR-0006 Q2).
	double seed_omega_bar_ = kDefaultInitOmegaBar;

	/// Number of first-order ODE integrations performed by this solver. Diagnostic only,
	/// readable through `RotationSolverTestSeam`; it exists so that "one raw solve per built
	/// star, and none for a rescaling" can be asserted rather than assumed.
	std::size_t first_order_solve_count_ = 0;

	/// Seed-free first-order response of the attached star, rebuilt by every solve.
	HartleFirstOrderResponse first_order_response_;

	/// Solution to TOV is saved here
	// TOVTable tov_solution ;
	std::vector<OmegaPoint> omega_results;	  // Aug 6, 2020
	std::vector<OmegaSeqPoint> omega_seq_pts; // Aug 6, 2020

	std::vector<OmegaPointDark> omega_results_dark; // Apr 20, 2022

	/// Radius of the star
	// double R_Star = -1 ;

	/// Mass of the star (it will set after importing TOV solution)
	// double M_Star = -1 ;

	// double fast_p_v;
	// double fast_p_d;
	// double fast_e_v;
	// double fast_e_d;
	// double fast_m_tot;
	// -------- Fast Hartle RHS: profile-backed interpolation (thread-safe if solver is not shared) --------
	const Zaki::Vector::DataColumn *fast_r_ = nullptr; // km
	const Zaki::Vector::DataColumn *fast_p_ = nullptr; // geom units (consistent with Hartle equations)
	const Zaki::Vector::DataColumn *fast_e_ = nullptr;
	const Zaki::Vector::DataColumn *fast_m_ = nullptr;

	mutable std::size_t fast_k_ = 0; // cached bracket index for interpolation

	// Mixed-star: total quantities on a single radius grid (we will point these to prebuilt arrays)
	const Zaki::Vector::DataColumn *fast_r_mix_ = nullptr;
	const Zaki::Vector::DataColumn *fast_p_tot_ = nullptr;
	const Zaki::Vector::DataColumn *fast_e_tot_ = nullptr;
	const Zaki::Vector::DataColumn *fast_m_tot_ = nullptr;

	mutable std::size_t fast_k_mix_ = 0;

	void SetFastProfilePtrs_(const Zaki::Vector::DataColumn &r,
							 const Zaki::Vector::DataColumn &p,
							 const Zaki::Vector::DataColumn &e,
							 const Zaki::Vector::DataColumn &m);

	void SetFastMixedPtrs_(const Zaki::Vector::DataColumn &r,
						   const Zaki::Vector::DataColumn &p_tot,
						   const Zaki::Vector::DataColumn &e_tot,
						   const Zaki::Vector::DataColumn &m_tot);

	inline void EvalFastPEM_(double r, double &p, double &e, double &m) const;
	inline void EvalFastMixedPEM_(double r, double &p, double &e, double &m) const;

	/// Raw first-order omega_bar profile at the internal seed, stored per node by
	/// `FindNMomInertia` and divided by `Omega_raw` there to build the seed-free
	/// `first_order_response_`. **First-order machinery** (ADR-0006 Q4) — the O(Omega^2)
	/// solver never reads these; it consumes only the seed-free response.
	Zaki::Vector::DataColumn stored_omega_bar_;
	Zaki::Vector::DataColumn stored_domega_bar_;

	// ---------------------------------------------------------------------------------
	//  O(Omega^2) monopole (l = 0) — ADR-0007. Replaces the retired 675b4a9 candidate.
	// ---------------------------------------------------------------------------------
	/// Background sampled at the ACTUAL radius the ODE driver asks for.
	///
	/// The retired candidate set per-node scalars before each `gsl_odeiv2_driver_apply` and let
	/// the right-hand side use them at whatever internal radius GSL chose — a node value
	/// standing in for an inter-node one. This workspace instead interpolates every input at
	/// the true RHS radius, through **one shared bracket index** so that all eight quantities
	/// refer to the same radial interval, and linearly, which is the interpolation order the
	/// profile itself carries (INV-13).
	struct MonopoleBackground_
	{
		const Zaki::Vector::DataColumn *r = nullptr;	///< [km]
		const Zaki::Vector::DataColumn *p = nullptr;	///< [km^-2]
		const Zaki::Vector::DataColumn *e = nullptr;	///< [km^-2]
		const Zaki::Vector::DataColumn *m = nullptr;	///< [km]
		const Zaki::Vector::DataColumn *nu = nullptr;	///< dimensionless
		const Zaki::Vector::DataColumn *nup = nullptr;	///< [km^-1]
		const Zaki::Vector::DataColumn *dedp = nullptr; ///< dimensionless (ADR-0007 P5); since
														///< ADR-0008 Q8 this is consumed ONLY by
														///< the regular-centre series, never as
														///< the radial mass source
		const Zaki::Vector::DataColumn *s = nullptr;	///< omega_bar/Omega, dimensionless
		const Zaki::Vector::DataColumn *sp = nullptr;	///< [d omega_bar/dr]/Omega, [km^-1]

		mutable std::size_t k = 0; ///< shared cached bracket index

		/// ADR-0008 Q3 — the energy-density MEASURE density of the profile segment the driver
		/// is currently inside: `d(eps)/dr|_segment = (eps_{i+1} - eps_i)/(r_{i+1} - r_i)`
		/// [km^-3]. Piecewise constant by construction, so the segment's total contribution is
		/// exactly `Delta eps_i` whatever the EOS does between the two nodes. Installed by
		/// `ComputeMonopoleResponse()` before each one-segment `gsl_odeiv2_driver_apply`; it is
		/// deliberately NOT interpolated and NOT reconstructed from `dedp * dp/dr`.
		double eps_slope = 0.0;

		/// ADR-0007 P4 regular-centre limit numerator `j_c^2 s_c^2` for
		/// `xi0_hat -> j_c^2 s_c^2 r / [4 pi (eps + 3 p)]`, used by the right-hand side on the
		/// same terms as the derived `xi0` column when `p0*/nu'` is ill-conditioned. Zero until
		/// the centre initialization has run.
		double centre_xi_num = 0.0;

		struct Sample
		{
			double p = 0, e = 0, m = 0, nu = 0, nup = 0, dedp = 0, s = 0, sp = 0;
		};

		[[nodiscard]] bool Complete() const noexcept;
		[[nodiscard]] Sample At(double r_km) const;
	};

	MonopoleBackground_ monopole_bg_;

	/// The published response. Never mutated incrementally: a complete candidate is built in
	/// local storage and assigned here only once every acceptance test has passed.
	HartleMonopoleResponse monopole_response_;

	/// Number of monopole ODE integrations performed by this solver. Read through
	/// `RotationSolverTestSeam` so "one solve per explicit request, none per materialization"
	/// can be proved rather than asserted.
	std::size_t monopole_solve_count_ = 0;

	//--------------------------------------------------------------
  public:
	RotationSolver();
	~RotationSolver();

	/// Copy Constructor
	RotationSolver(const RotationSolver &) = delete;

	/// Assignment operator
	RotationSolver &operator=(const RotationSolver &) = delete;

	// void ImportTOVSolution(const Zaki::String::Directory&) ;
	// void ImportTOVSolution(const TOVTable&) ;

	/// Attaches the NStar pointer to RotationSolver class
	void AttachNStar(NStar *);

	void AttachMixedStar(MixedStar *);

	NStar *GetNStar();

	MixedStar *GetMixedStar();

	/// Returns the mass at a given radius
	double GetMass(const double &);

	/// Returns the pressure at a given radius
	double GetPress(double);

	/// Returns the derivative of pressure at a given radius
	// double GetPressDer(const double&) ;

	/// Returns the energy density at a given radius
	double GetEDens(const double &);

	/// Returns the initial value assumes for omega at r ~ 0
	double GetInitOmegaBar() const;

	/// Returns omega_seq_pts
	std::vector<OmegaSeqPoint> GetOmegaSeq() const;

	/// .......................Aug 6, 2020.........................
	// Returns the metric function nu(r)
	double GetNu(double in_R);

	// nu(r) integrand
	// double NuIntegrand(double r) ;

	// // Returns the value of the j function (Hartle Eq. 6.17, P. 253)
	// double GetHartleJ(double r) ;

	// // Returns the derivative of the j function (Hartle Eq. 6.18, P. 253)
	// double GetHartleJDer(double r) ;

	// /// Coefficient for y[0]
	// double GetHartleOmegaCoeff(const double r) ;

	// /// Coefficient for y[1]
	// double GetHartleDOmegaCoeff(const double r) ;

	/// NStar: Coefficient for y[0]
	double GetHartleOmegaCoeff_N_Fast(const double r);

	/// NStar: Coefficient for y[1]
	double GetHartleDOmegaCoeff_N_Fast(const double r);

	/// MixedStar: Coefficient for y[0]
	double GetHartleOmegaCoeff_Mixed(const double r);

	/// MixedStar: Coefficient for y[1]
	double GetHartleDOmegaCoeff_Mixed(const double r);

	/// MixedStar: Coefficient for y[0]
	double GetHartleOmegaCoeff_Mixed_Fast(const double r);

	/// MixedStar: Coefficient for y[1]
	double GetHartleDOmegaCoeff_Mixed_Fast(const double r);

	/// Input a range of initial omega(0) values
	void Solve(const Zaki::Math::Axis &,
			   const Zaki::String::Directory & = "");

	/// Input an initial omega(0) value
	void Solve(const double &,
			   const Zaki::String::Directory & = "");

	/// MixedStar: Input an initial omega(0) value
	void Solve_Mixed(const double &,
					 const Zaki::String::Directory & = "");

	/// Evaluates the moment of inertia for the neutron star
	void FindNMomInertia();

	/// Evaluates the moment of inertia for the mixed star
	void FindMixedMomInertia();

	/// Compute the governed O(Omega^2) monopole response of the attached star (ADR-0007).
	///
	/// **This performs work**: one ODE integration over the star. It is deliberately not run
	/// by `FindNMomInertia`, and therefore not by ordinary star construction — existing
	/// workflows do not need it and must not pay for it.
	///
	/// Recomputes when the cached response is absent or stale (source profile changed or its
	/// `Version()` moved); reuses it otherwise. Fails closed — leaving no partial response —
	/// if the profile is incomplete, the first-order response is unusable, the profile carries
	/// no authoritative `d(eps)/dp` (ADR-0007 P5 — still required for the regular-centre
	/// series, ADR-0008 Q8), the radial partition is not strictly increasing or carries a
	/// non-finite energy-density increment (ADR-0008 Q3), the integration fails, or any
	/// published value would be non-finite.
	///
	/// @return true if a valid, current response is available afterwards.
	bool ComputeMonopoleResponse();

	/// The current O(Omega^2) monopole response, or nullptr when none has been computed or the
	/// cached one is stale. A cheap read-only accessor: it never integrates. Call
	/// `ComputeMonopoleResponse()` to (re)compute.
	[[nodiscard]] const HartleMonopoleResponse *MonopoleResponse() const;

	/// The seed-free first-order response of the attached star (ADR-0006 Q4).
	///
	/// Valid after `FindNMomInertia()`. Contains no arbitrary normalization and no physical
	/// spin; call `HartleFirstOrderResponse::At()` with an explicit `AngularVelocity` to
	/// materialize a physical first-order solution.
	const HartleFirstOrderResponse &FirstOrderResponse() const;

	/// ..........................................................
	/// Exports the results of solving the rotation equations
	void ExportResults(const Zaki::String::Directory &) const;

	static int ODE(double r, const double y[], double f[], void *params);
	static int ODE_Mixed(double r, const double y[], double f[], void *params);
	static int ODE_Mixed_Out(double r, const double y[], double f[], void *params);
	static int ODE_Mixed_Fast(double r, const double y[], double f[], void *params);
	static int ODE_N_Fast(double r, const double y[], double f[], void *params);

	/// The governed O(Omega^2) monopole right-hand side (ADR-0007 P2 as amended by ADR-0008),
	/// in the normalized variables `y[0] = m0/Omega^2 [km^3]`, `y[1] = p0*/Omega^2 [km^2]`.
	///
	/// The EOS energy-density contribution to `dm0/dr` is the **measure**
	/// `-4 pi r^2 xi0_hat d(eps)` (H67 eq. 93; ADR-0008 Q1/Q3), evaluated inside one governed
	/// profile segment at a time through `MonopoleBackground_::eps_slope`. It is **not** the
	/// pointwise `4 pi r^2 (eps+p)(deps/dp) p0*` rewrite, which cannot represent
	/// energy-density variation narrower than the profile spacing.
	static int ODE_HartleMonopole_(double r, const double y[], double f[], void *params);

	// Resets the containers
	void Reset();

	// Test-only seam (ADR-0006 Q2). Grants no public API; see RotationSolverTestSeam.
	friend struct RotationSolverTestSeam;
};

//==============================================================
} // namespace CompactStar::Core
//==============================================================
#endif /*CompactStar_Core_RotationSolver_H*/
