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
/// Stores the results of both first-order (omega_bar) and second-order
/// (m0, p0, xi0) Hartle slow-rotation solutions on the radial grid.
struct HartleResult
{
	// --- First-order O(Omega) ---
	double Omega = 0.0;		///< Angular velocity [s^-1]
	double J = 0.0;			///< Angular momentum [km^2]
	double I = 0.0;			///< Moment of inertia [km^3]

	Zaki::Vector::DataColumn omega_bar;	 ///< omega_bar(r) [km^-1]
	Zaki::Vector::DataColumn domega_bar; ///< d(omega_bar)/dr(r) [km^-2]

	// --- Second-order O(Omega^2) monopole (l=0) ---
	Zaki::Vector::DataColumn m0;  ///< Mass perturbation delta_m(r) [km]
	Zaki::Vector::DataColumn p0;  ///< Pressure perturbation delta_p(r) [km^-2]
	Zaki::Vector::DataColumn xi0; ///< Isobar displacement xi_0(r) [km]

	double p0_c = 0.0;		///< Central pressure perturbation (shooting result)
	double delta_M = 0.0;	///< O(Omega^2) gravitational mass correction [km]

	/// Grid reference (non-owning pointer to star's radius column)
	const Zaki::Vector::DataColumn *r_grid = nullptr;

	bool valid = false; ///< True after a successful second-order solve
};
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

	/// Initial bar{omega} (at r = r(i=0))
	double init_omega_bar = -1;

	/// Solution to TOV is saved here
	// TOVTable tov_solution ;
	std::vector<OmegaPoint> omega_results;	  // Aug 6, 2020
	std::vector<OmegaSeqPoint> omega_seq_pts; // Aug 6, 2020

	std::vector<OmegaPointDark> omega_results_dark; // Apr 20, 2022

	/// Radius of the star
	// double R_Star = -1 ;

	/// Mass of the star (it will set after importing TOV solution)
	// double M_Star = -1 ;

	// double fast_p;
	// double fast_e;
	// double fast_m;

	double fast_p_v;
	double fast_p_d;
	double fast_e_v;
	double fast_e_d;
	double fast_m_tot;

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

	/// Returns the second-order Hartle result (const)
	const HartleResult &GetHartleResult() const;

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
};

//==============================================================
} // namespace CompactStar::Core
//==============================================================
#endif /*CompactStar_Core_RotationSolver_H*/
