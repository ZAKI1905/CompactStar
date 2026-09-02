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
/// Result of the **O(Omega^2) second-order candidate** (`SolveHartle2_N`).
///
/// **UNVERIFIED SCIENTIFIC CANDIDATE** under `GOVERNANCE.md` §5 and INV-08: publicly callable,
/// with zero repository callers, and its equations are recorded as defective
/// (`docs/validation/PHASE4_ROTATION_ENTRY.md` §10-§12). Nothing here is validated physics.
///
/// Every field below inherits the **square** of the arbitrary internal first-order seed,
/// because the solve is driven by the raw stored `omega_bar` (measured: entry record §12).
/// ADR-0006 P9 requires any future second-order product to be seed-free; that is Phase-4C work
/// and is deliberately NOT done here.
///
/// The first-order fields this struct used to carry — `Omega`, `J`, `I`, `omega_bar`,
/// `domega_bar` — were removed by ADR-0006: `Omega` was annotated `[s^-1]` while storing km^-1,
/// and all five held seed-normalized values behind a public accessor. First-order results now
/// live in `HartleFirstOrderResponse` and `PhysicalFirstOrderRotation`, which are seed-free by
/// construction.
struct HartleResult
{
	// --- Second-order O(Omega^2) monopole (l=0) --- CANDIDATE, UNVALIDATED ---
	Zaki::Vector::DataColumn m0;  ///< Mass perturbation delta_m(r) [km]
	Zaki::Vector::DataColumn p0;  ///< Pressure perturbation delta_p(r) [km^-2]
	Zaki::Vector::DataColumn xi0; ///< Isobar displacement xi_0(r) [km]

	double p0_c = 0.0;		///< Central pressure perturbation (shooting result) [km^-2]
	double delta_M = 0.0;	///< O(Omega^2) gravitational mass correction [km]

	/// Grid reference (non-owning pointer to star's radius column)
	const Zaki::Vector::DataColumn *r_grid = nullptr;

	bool valid = false; ///< True after a successful second-order solve
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

	// Per-grid-point scalars used ONLY by the second-order O(Omega^2) candidate
	// path (ODE_Hartle2_N_Fast, SolveHartle2_N). Not read by any first-order code,
	// which uses the profile-backed fast_p_/fast_e_/fast_m_ interpolation below.
	double fast_p;
	double fast_e;
	double fast_m;

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

	// Fast interpolation cache for second-order Hartle solver
	double fast_nu;		  ///< Metric nu at current grid point
	double fast_nu_prime; ///< dnu/dr at current grid point
	double fast_dEdP;	  ///< d(eps)/d(p) at current grid point (NStar)
	double fast_dEdP_v;	  ///< d(eps)/d(p) visible sector (MixedStar)
	double fast_dEdP_d;	  ///< d(eps)/d(p) dark sector (MixedStar)

	/// Cached omega_bar profile for second-order integration
	/// (stored as DataColumn for interpolation during m0/p0 solve)
	Zaki::Vector::DataColumn stored_omega_bar_;
	Zaki::Vector::DataColumn stored_domega_bar_;

	// Fast interpolation cache for stored omega_bar profile
	double fast_omega_bar;
	double fast_domega_bar;

	/// Second-order Hartle results
	HartleResult hartle_result_;

	/// Whether to include source terms in the (m0,p0) ODE
	/// (false for homogeneous solution in superposition method)
	bool include_m0p0_source_ = true;

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

	/// Solves the second-order Hartle O(Omega^2) monopole equations
	/// for an NStar. Requires FindNMomInertia() to have been called first
	/// (or uses its own first-order solve). Stores results in hartle_result_.
	void SolveHartle2_N();

	/// Solves the second-order Hartle O(Omega^2) monopole equations
	/// for a MixedStar.
	void SolveHartle2_Mixed();

	/// Returns the second-order Hartle result (const).
	///
	/// **UNVERIFIED SCIENTIFIC CANDIDATE** (`GOVERNANCE.md` §5, INV-08). Its contents are
	/// quadratic in the arbitrary internal seed and are not physical results.
	const HartleResult &GetHartleResult() const;

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

	/// Second-order (m0, p0) ODE for NStar with fast interpolation
	static int ODE_Hartle2_N_Fast(double r, const double y[], double f[], void *params);

	/// Second-order (m0, p0) ODE for MixedStar with fast interpolation
	static int ODE_Hartle2_Mixed_Fast(double r, const double y[], double f[], void *params);

	// Resets the containers
	void Reset();

	// Test-only seam (ADR-0006 Q2). Grants no public API; see RotationSolverTestSeam.
	friend struct RotationSolverTestSeam;
};

//==============================================================
} // namespace CompactStar::Core
//==============================================================
#endif /*CompactStar_Core_RotationSolver_H*/
