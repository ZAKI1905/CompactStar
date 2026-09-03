/*
  Last edited in 2026
  RotationSolver class
*/
#include <cmath>
#include <limits>
#include <stdexcept>
#include <gsl/gsl_math.h>
// #include <gsl/gsl_roots.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h> // Aug 6, 2020
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h> // Aug 6, 2020

#include <Zaki/File/CSVIterator.hpp>
#include <Zaki/File/VecSaver.hpp>
#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <Zaki/Util/Instrumentor.hpp>
// #include <Zaki/Vector/DataColumn.hpp>

#include "CompactStar/Core/MixedStar.hpp"
#include "CompactStar/Geometry.hpp"
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/RotationSolver.hpp"

#define R_SOLVER_VERBOSE 0

using namespace CompactStar::Core;
//==============================================================
namespace
{
constexpr double kR_EPS_KM = 1e-6; // 1e-6 km = 0.1 cm; safely nonzero in km units

inline double SafeR0(double r_min_km)
{
	if (std::isfinite(r_min_km) && r_min_km > kR_EPS_KM)
		return r_min_km;
	return kR_EPS_KM;
}
} // anonymous namespace

//==============================================================
//                        RotationSolver class
//==============================================================
// Constructor
RotationSolver::RotationSolver() : Prog("RotationSolver")
{
}

//--------------------------------------------------------------
RotationSolver::~RotationSolver()
{
}

//--------------------------------------------------------------
void RotationSolver::AttachNStar(NStar *ns_ptr)
{
	Z_LOG_INFO("Attaching NStar to RotationSolver class...");

	if (!ns_ptr)
	{
		Z_LOG_ERROR("The NStar pointer is a nullptr!");
		return;
	}

	nstar_ptr = ns_ptr;
}

//--------------------------------------------------------------
void RotationSolver::AttachMixedStar(MixedStar *m_ns_ptr)
{
	Z_LOG_INFO("Attaching MixedStar to RotationSolver class...");

	if (!m_ns_ptr)
	{
		Z_LOG_ERROR("The MixedStar pointer is a nullptr!");
		return;
	}

	mixedstar_ptr = m_ns_ptr;
}

//--------------------------------------------------------------
void RotationSolver::SetFastProfilePtrs_(const Zaki::Vector::DataColumn &r,
										 const Zaki::Vector::DataColumn &p,
										 const Zaki::Vector::DataColumn &e,
										 const Zaki::Vector::DataColumn &m)
{
	fast_r_ = &r;
	fast_p_ = &p;
	fast_e_ = &e;
	fast_m_ = &m;
	fast_k_ = 0;
}

//--------------------------------------------------------------
void RotationSolver::SetFastMixedPtrs_(const Zaki::Vector::DataColumn &r,
									   const Zaki::Vector::DataColumn &p_tot,
									   const Zaki::Vector::DataColumn &e_tot,
									   const Zaki::Vector::DataColumn &m_tot)
{
	fast_r_mix_ = &r;
	fast_p_tot_ = &p_tot;
	fast_e_tot_ = &e_tot;
	fast_m_tot_ = &m_tot;
	fast_k_mix_ = 0;
}

static inline double Lerp_(double x, double x0, double x1, double y0, double y1)
{
	const double t = (x - x0) / (x1 - x0);
	return y0 + t * (y1 - y0);
}

inline void RotationSolver::EvalFastPEM_(double r, double &p, double &e, double &m) const
{
	// Preconditions
	const auto &R = *fast_r_;
	const auto &P = *fast_p_;
	const auto &E = *fast_e_;
	const auto &M = *fast_m_;

	const std::size_t n = R.Size();
	if (n < 2)
	{
		p = P.Empty() ? 0 : P[0];
		e = E.Empty() ? 0 : E[0];
		m = M.Empty() ? 0 : M[0];
		return;
	}

	// Clamp
	if (r <= R.Values().front())
	{
		p = P.Values().front();
		e = E.Values().front();
		m = M.Values().front();
		return;
	}
	if (r >= R.Values().back())
	{
		p = P.Values().back();
		e = E.Values().back();
		m = M.Values().back();
		return;
	}

	// Ensure fast_k_ is a valid left bracket
	std::size_t k = fast_k_;
	if (k >= n - 1)
		k = n - 2;

	// Walk k to bracket r: R[k] <= r <= R[k+1]
	while (k + 1 < n && r > R[k + 1])
		++k;
	while (k > 0 && r < R[k])
		--k;

	fast_k_ = k;

	const double r0 = R[k], r1 = R[k + 1];
	p = Lerp_(r, r0, r1, P[k], P[k + 1]);
	e = Lerp_(r, r0, r1, E[k], E[k + 1]);
	m = Lerp_(r, r0, r1, M[k], M[k + 1]);
}

inline void RotationSolver::EvalFastMixedPEM_(double r, double &p, double &e, double &m) const
{
	const auto &R = *fast_r_mix_;
	const auto &P = *fast_p_tot_;
	const auto &E = *fast_e_tot_;
	const auto &M = *fast_m_tot_;

	const std::size_t n = R.Size();
	if (n < 2)
	{
		p = P.Empty() ? 0 : P[0];
		e = E.Empty() ? 0 : E[0];
		m = M.Empty() ? 0 : M[0];
		return;
	}

	if (r <= R.Values().front())
	{
		p = P.Values().front();
		e = E.Values().front();
		m = M.Values().front();
		return;
	}
	if (r >= R.Values().back())
	{
		p = P.Values().back();
		e = E.Values().back();
		m = M.Values().back();
		return;
	}

	std::size_t k = fast_k_mix_;
	if (k >= n - 1)
		k = n - 2;

	while (k + 1 < n && r > R[k + 1])
		++k;
	while (k > 0 && r < R[k])
		--k;

	fast_k_mix_ = k;

	const double r0 = R[k], r1 = R[k + 1];
	p = Lerp_(r, r0, r1, P[k], P[k + 1]);
	e = Lerp_(r, r0, r1, E[k], E[k + 1]);
	m = Lerp_(r, r0, r1, M[k], M[k + 1]);
}

//--------------------------------------------------------------
// Input radius, output energy density
double RotationSolver::GetEDens(const double &in_R)
{
	return nstar_ptr->GetEnergyDensity(in_R);
}

//--------------------------------------------------------------
// Input radius, output metric function nu(r) (Aug 6, 2020)
// [ Aug 26, 2021 ]: The output is only valid inside the star
// for outside the star this should simply return:
//  0.5*log( 1 - 2*M_Star/R_Star)
//  !!! Radius is in km here !!!
//  Changed on Mar 23, 2022 !
double RotationSolver::GetNu(double in_R)
{
	return nstar_ptr->GetMetricNu(in_R);
}

//--------------------------------------------------------------
// Added on Apr 20, 2022
// Edited on May 3, 2022
int RotationSolver::ODE_Mixed(double r, const double y[], double f[], void *params)
{
	RotationSolver *rot_obj = (RotationSolver *)params;

	// First derivative of bar{omega(r)}
	f[0] = y[1];

	// Second derivative of bar{omega(r)}
	const double r_safe = (r > kR_EPS_KM) ? r : kR_EPS_KM;
	f[1] = -4. * y[1] / r_safe + rot_obj->GetHartleOmegaCoeff_Mixed(r) * y[0] + y[1] * rot_obj->GetHartleDOmegaCoeff_Mixed(r);

	return GSL_SUCCESS;
}

//--------------------------------------------------------------
// Created on June 3, 2022
int RotationSolver::ODE_N_Fast(double r, const double y[], double f[], void *params)
{
	RotationSolver *rot_obj = (RotationSolver *)params;

	// First derivative of bar{omega(r)}
	f[0] = y[1];

	// Second derivative of bar{omega(r)}
	const double r_safe = (r > kR_EPS_KM) ? r : kR_EPS_KM;
	f[1] = -4. * y[1] / r_safe + rot_obj->GetHartleOmegaCoeff_N_Fast(r) * y[0] + y[1] * rot_obj->GetHartleDOmegaCoeff_N_Fast(r);

	return GSL_SUCCESS;
}

//--------------------------------------------------------------
// Created on May 3, 2022
int RotationSolver::ODE_Mixed_Fast(double r, const double y[], double f[], void *params)
{
	RotationSolver *rot_obj = (RotationSolver *)params;

	// First derivative of bar{omega(r)}
	f[0] = y[1];

	// Second derivative of bar{omega(r)}
	const double r_safe = (r > kR_EPS_KM) ? r : kR_EPS_KM;
	f[1] = -4. * y[1] / r_safe + rot_obj->GetHartleOmegaCoeff_Mixed_Fast(r) * y[0] + y[1] * rot_obj->GetHartleDOmegaCoeff_Mixed_Fast(r);

	return GSL_SUCCESS;
}

//--------------------------------------------------------------
// Created on May 3, 2022
int RotationSolver::ODE_Mixed_Out(double r, const double y[], double f[], void *params)
{
	RotationSolver *rot_obj = (RotationSolver *)params;

	// First derivative of bar{omega(r)}
	f[0] = y[1];

	// Second derivative of bar{omega(r)}
	const double r_safe = (r > kR_EPS_KM) ? r : kR_EPS_KM;
	f[1] = -4. * y[1] / r_safe;

	return GSL_SUCCESS;
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[0]
double RotationSolver::GetHartleOmegaCoeff_Mixed(const double r)
{
	return 16. * M_PI * (mixedstar_ptr->GetPress_Visible(r) + mixedstar_ptr->GetPress_Dark(r) + mixedstar_ptr->GetEps_Visible(r) + mixedstar_ptr->GetEps_Dark(r)) / (1. - 2. * mixedstar_ptr->GetMass_Total(r) / r);
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[1]
double RotationSolver::GetHartleDOmegaCoeff_Mixed(const double r)
{
	return 4. * M_PI * (mixedstar_ptr->GetPress_Visible(r) + mixedstar_ptr->GetPress_Dark(r) + mixedstar_ptr->GetEps_Visible(r) + mixedstar_ptr->GetEps_Dark(r)) * r / (1. - 2. * mixedstar_ptr->GetMass_Total(r) / r);
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[0]
double RotationSolver::GetHartleOmegaCoeff_N_Fast(const double r)
{
	double p, e, m;
	EvalFastPEM_(r, p, e, m);
	return 16. * M_PI * (p + e) / (1. - 2. * m / r);
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[1]
double RotationSolver::GetHartleDOmegaCoeff_N_Fast(const double r)
{
	double p, e, m;
	EvalFastPEM_(r, p, e, m);
	return 4. * M_PI * (p + e) * r / (1. - 2. * m / r);
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[0]
double RotationSolver::GetHartleOmegaCoeff_Mixed_Fast(const double r)
{
	// return 16. * M_PI * (fast_p_v + fast_p_d + fast_e_v + fast_e_d) / (1. - 2. * fast_m_tot / r);
	double p, e, m;
	EvalFastMixedPEM_(r, p, e, m);
	return 16. * M_PI * (p + e) / (1. - 2. * m / r);
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[1]
double RotationSolver::GetHartleDOmegaCoeff_Mixed_Fast(const double r)
{
	// return 4. * M_PI * (fast_p_v + fast_p_d + fast_e_v + fast_e_d) * r / (1. - 2. * fast_m_tot / r);
	double p, e, m;
	EvalFastMixedPEM_(r, p, e, m);
	return 4. * M_PI * (p + e) * r / (1. - 2. * m / r);
}

//--------------------------------------------------------------
// bar{omega(r)} not omega
double RotationSolver::GetInitOmegaBar() const
{
	return init_omega_bar;
}

//--------------------------------------------------------------
// Added on Apr 20, 2022
// Input an initial omega(0) value
// Need to fix this method!!!!
void RotationSolver::Solve_Mixed(const double &in_omega_0,
								 const Zaki::String::Directory &in_dir)
{
	PROFILE_FUNCTION();

	// Resolution of the solver (division of the radial distance)
	// unsigned int solver_res = 5000 ;

	// double r_min = mixedstar_ptr->core_region.Min();
	// double r_max = 1.1 * mixedstar_ptr->mantle_region.Max();

	double r_min_raw = mixedstar_ptr->core_region.Min();
	double r0 = SafeR0(r_min_raw);

	double r_max = 1.1 * mixedstar_ptr->mantle_region.Max();

#if R_SOLVER_VERBOSE
	std::cout << "\n\n\t\t ****************************************"
			  << "*******************************" << " \n";
	std::cout << "\t\t *                   "
			  << "Rotation Solver Sequence Results"
			  << "                  * \n";
	std::cout << "\t\t ******************************************"
			  << "*****************************" << "\n\n";
#endif

	init_omega_bar = in_omega_0;
	double r = r0;
	//  double r1 = r_max ;

	double y[2];
	y[0] = init_omega_bar;
	y[1] = 0;

	// ------------------------------------------------
	// Inside the mixed star:
	gsl_odeiv2_system ode_sys = {RotationSolver::ODE_Mixed, nullptr, 2, this};

	gsl_odeiv2_driver *tmp_driver = gsl_odeiv2_driver_alloc_y_new(&ode_sys, gsl_odeiv2_step_rk8pd,
																  1.e-1, 1.e-10, 1.e-10);
	// ------------------------------------------------
	// Outside the mixed star
	gsl_odeiv2_system ode_sys_out = {RotationSolver::ODE_Mixed_Out, nullptr, 2, this};

	gsl_odeiv2_driver *tmp_driver_out = gsl_odeiv2_driver_alloc_y_new(&ode_sys_out, gsl_odeiv2_step_rk8pd,
																	  1.e-1, 1.e-10, 1.e-10);
	// ------------------------------------------------

	//  double min_log_r  = log10(r_min) ;
	//  double max_log_r  = log10(r_max) ;
	// -------------------------------------------------
	// Added on Apr 20, 2022
	double r_surface = mixedstar_ptr->mantle_region.Max();

	// Added on March 23, 2022 to fix the bug
	// in J, Omega, I calculations.
	//  double log_r_star = log10(r_surface) ;

	// Total angular momentum and velocity:
	double ang_mom_J, ang_vel_Omega, mom_inertia;
	// -------------------------------------------------

	//  double step = (max_log_r - min_log_r) / radial_res ;
	double step = (r_max - r0) / radial_res;

	// One extra point for the exact value at r = R_star.
	omega_results.reserve(radial_res + 1);
	//  bool surface_reached = false ;

	// ---------------------------------------------------------------------
	//                         INNER RADIUS LOOP BEGINS
	// ---------------------------------------------------------------------
	for (double r_i = r0; r_i < r_surface; r_i += step)
	{
		//   PROFILE_SCOPE("R-Loop") ;

		//    double ri = pow(10, log_r_i) ;

		// This function evolves the driver system d from t to t1.
		// Initially vector y should contain the values of dependent
		// variables at point t.
		int status = gsl_odeiv2_driver_apply(tmp_driver, &r, r_i, y);

		if (status != GSL_SUCCESS)
		{
#if R_SOLVER_VERBOSE
			printf("error, return value=%d\n", status);
#endif
			break;
		}

		// Units : { [ km ], [ M_Sun ], - , [ km^-1 ], [ km^-2 ] }
		omega_results_dark.emplace_back(r,
										mixedstar_ptr->GetMass_Visible(r) / Zaki::Physics::SUN_M_KM,
										mixedstar_ptr->GetMass_Dark(r) / Zaki::Physics::SUN_M_KM,
										mixedstar_ptr->GetNu(r), y[0], y[1]);

		// March 23, 2022
		//    if ( !surface_reached && abs(log_r_i - log_r_star) < step )
		//    {
		//
		//      status = gsl_odeiv2_driver_apply (tmp_driver, &r, r_surface, y) ;
		//
		//      omega_results_dark.emplace_back(r,
		//                  mixedstar_ptr->GetMass_Visible(r) / Zaki::Physics::SUN_M_KM,
		//                  mixedstar_ptr->GetMass_Dark(r) / Zaki::Physics::SUN_M_KM,
		//                  mixedstar_ptr->GetNu(r), y[0], y[1]) ;
		//
		//      // Total angular momentum and velocity:
		//      ang_mom_J = pow(r_surface,4)*y[1]/6. ;
		//      ang_vel_Omega = y[0] + r_surface*y[1]/3. ;
		//      mom_inertia = ang_mom_J / ang_vel_Omega ;
		//
		//      surface_reached = true ;
		//    }
	}
	// ............ Inner Radius Loop Ends .........................
	gsl_odeiv2_driver_apply(tmp_driver, &r, r_surface, y);

	omega_results_dark.emplace_back(r,
									mixedstar_ptr->GetMass_Visible(r) / Zaki::Physics::SUN_M_KM,
									mixedstar_ptr->GetMass_Dark(r) / Zaki::Physics::SUN_M_KM,
									mixedstar_ptr->GetNu(r), y[0], y[1]);

	// Total angular momentum and velocity:
	ang_mom_J = pow(r_surface, 4) * y[1] / 6.;
	ang_vel_Omega = y[0] + r_surface * y[1] / 3.;
	mom_inertia = ang_mom_J / ang_vel_Omega;

	// ---------------------------------------------------------------------
	//                         OUTER RADIUS LOOP BEGINS
	// ---------------------------------------------------------------------
	for (double r_i = r_surface + step; r_i < r_max; r_i += step)
	{
		if (GSL_SUCCESS != gsl_odeiv2_driver_apply(tmp_driver_out, &r, r_i, y))
			break;

		// Units : { [ km ], [ M_Sun ], - , [ km^-1 ], [ km^-2 ] }
		omega_results_dark.emplace_back(r,
										mixedstar_ptr->GetMass_Visible(r) / Zaki::Physics::SUN_M_KM,
										mixedstar_ptr->GetMass_Dark(r) / Zaki::Physics::SUN_M_KM,
										mixedstar_ptr->GetNu(r), y[0], y[1]);
	}
	// ............ Outer Radius Loop Ends .........................

	// ............................. omega_results_dark Loop Begins .....................................
	// We go back now and update the omega_results_dark with the values for omega(r)
	// We needed ang_vel_Omega for omega(r) since:
	// omega(r) = ang_vel_Omega - bar{omega(r)}
	// The unit for omega(r) is changed from km^{-1} to s^{-1}
	for (size_t i = 0; i < omega_results_dark.size(); i++)
	{
		omega_results_dark[i].omega = (ang_vel_Omega - omega_results_dark[i].omega_bar) * Zaki::Physics::LIGHT_C_KM_S;
	}
	// .............................. omega_results_dark Loop Ends ......................................

#if R_SOLVER_VERBOSE
	// ...............................................................................................
	// Printing the sequence point data in terminal
	printf("\u2554%-40s \u03A9-Sequence %-3zu out of %-3lu %43s\u2557\n",
		   Zaki::String::Multiply("\u2550", 40).c_str(), (size_t)1, (size_t)1,
		   Zaki::String::Multiply("\u2550", 43).c_str());
	printf("\u2551 %14s %14s %17s %16s %14s %12s %12s %9s\n", "R (km)", "M (km)", "\u03B5_c (km^-2)",
		   "\u03C9c_bar (s^-1)", "\u03A9 (s^-1) ", "J (km^2)", "I (km^3)", "\u2551");
	printf("\u2551 %14le %14le %14le %14le %14le %14le %14le %7s\n",
		   r_surface,
		   mixedstar_ptr->mass[-1],
		   mixedstar_ptr->eps[0],
		   init_omega_bar * Zaki::Physics::LIGHT_C_KM_S,
		   ang_vel_Omega * Zaki::Physics::LIGHT_C_KM_S,
		   ang_mom_J,
		   mom_inertia,
		   "\u2551");
	printf("\u255A%110s\u255D\n", Zaki::String::Multiply("\u2550", 110).c_str());
#endif
	// ...............................................................................................

	omega_seq_pts.emplace_back(init_omega_bar * Zaki::Physics::LIGHT_C_KM_S,
							   mixedstar_ptr->sequence.v.m,
							   r_surface, ang_mom_J, ang_vel_Omega * Zaki::Physics::LIGHT_C_KM_S);

	gsl_odeiv2_driver_free(tmp_driver);

	// ------------------------------------------------------------
	//                      Saving to file
	// ------------------------------------------------------------
	if (false && in_dir.Str() != "")
	{
		std::vector<std::string> tmp_fname_v = Zaki::String::Pars(in_dir.Str(), "*");
		std::string tmp_fname;

		char tmp_omeg_str[100];
		snprintf(tmp_omeg_str, sizeof(tmp_omeg_str),
				 "%.3e", init_omega_bar * Zaki::Physics::LIGHT_C_KM_S);

		if (tmp_fname_v.size() < 2)
		{
			Z_LOG_WARNING("File name patter doesn't match '[]*[]'.");
			tmp_fname = tmp_fname_v[0] + std::string(tmp_omeg_str);
		}
		else
		{
			tmp_fname = tmp_fname_v[0] + std::string(tmp_omeg_str) + tmp_fname_v[1];
		}

		Zaki::File::VecSaver vec_saver(wrk_dir_ + "/" + tmp_fname);

		// .............................................
		// Header
		std::stringstream ss;
		char res_header[300];

		snprintf(res_header, sizeof(res_header),
				 "\u2554%-40s Sequence %-3zu out of %-3lu %45s\u2557\n",
				 Zaki::String::Multiply("\u2550", 40).c_str(), (size_t)1, (size_t)1,
				 Zaki::String::Multiply("\u2550", 45).c_str());
		ss << "# " << res_header;

		snprintf(res_header, sizeof(res_header),
				 "\u2551 %14s %12s %17s %16s %14s %12s %12s %11s\n", "R (km)", "M (Sun)", "\u03B5_c (g/cm^3)",
				 "\u03C9c_bar (s^-1)", "\u03A9 (s^-1) ", "J (km^2)", "I (km^3)", "\u2551");
		ss << "# " << res_header;

		snprintf(res_header, sizeof(res_header),
				 "\u2551 %14le %14le %14le %14le %14le %14le %14le %7s\n",
				 mixedstar_ptr->sequence.v.r,
				 mixedstar_ptr->sequence.v.m,
				 mixedstar_ptr->sequence.v.ec,
				 init_omega_bar * Zaki::Physics::LIGHT_C_KM_S,
				 ang_vel_Omega * Zaki::Physics::LIGHT_C_KM_S,
				 ang_mom_J,
				 mom_inertia,
				 "\u2551");
		ss << "# " << res_header;

		snprintf(res_header, sizeof(res_header), "\u255A%110s\u255D\n", Zaki::String::Multiply("\u2550", 110).c_str());
		ss << "# " << res_header;

		snprintf(res_header, sizeof(res_header), "%-19s\t %-19s\t %-19s\t %-19s\t %-19s\t %-19s\t %-19s",
				 "r [km]", "M_V [Sun]", "M_D [Sun]", "nu", "omega [1/s]",
				 "omega_bar [1/km]", "omega'_bar [1/km^2]");

		ss << res_header;

		// std::string tmp_label = res_header ;
		// tmp_label.c_str()
		// .............................................

		vec_saver.SetHeader(ss.str().c_str());
		vec_saver.Export1D(omega_results_dark);
	}
	// ------------------------------------------------------------

	// Testing this for freeing memory
	// omega_results_dark = std::vector<OmegaPointDark>() ;
	// instead of :
	omega_results_dark.clear();
	omega_results_dark.shrink_to_fit();

#if R_SOLVER_VERBOSE
	std::cout << "\n\t\t ************************"
			  << " Rot. Solver Finished *************************" << "\n\n";
#endif
	// ------------------------------------------------------------
}

//--------------------------------------------------------------
// Added on Aug 6, 2020
// Exports the results of solving the rotation equations
void RotationSolver::ExportResults(const Zaki::String::Directory &in_dir) const
{
	Zaki::File::VecSaver vec_saver_2(wrk_dir_ + "/" + in_dir);

	// Header tokens must describe the values actually written (ADR-0006 P7). `omega_seq_pts` is
	// populated ONLY by Solve_Mixed, the legacy two-fluid sequence path, from a caller-supplied
	// omega_bar SEED: `omega_bar_c` is that seed times c and `Omega` is the seed-normalized
	// angular velocity times c. Both are in 1/s, and neither is a requested physical spin — so
	// they are labelled `_seed` / `_seednorm` rather than promoted to physical quantities. `M`,
	// `R` and `J` previously carried no unit at all. The ordinary NStar first-order path never
	// populates this container and does not export through here; its physical results come from
	// NStar::RotationAt().
	char seq_header[300];
	snprintf(seq_header, sizeof(seq_header), "%-22s\t %-14s\t %-14s\t %-18s\t %-20s",
			 "omega_bar_c_seed (1/s)", "M (M_sun)", "R (km)", "J_seednorm (km^2)",
			 "Omega_seednorm (1/s)");

	vec_saver_2.SetHeader(seq_header);
	vec_saver_2.Export1D(omega_seq_pts);
}

//--------------------------------------------------------------
/// Returns omega_seq_pts
std::vector<OmegaSeqPoint> RotationSolver::GetOmegaSeq() const
{
	return omega_seq_pts;
}

//--------------------------------------------------------------
NStar *RotationSolver::GetNStar()
{
	return nstar_ptr;
}

//--------------------------------------------------------------
MixedStar *RotationSolver::GetMixedStar()
{
	return mixedstar_ptr;
}

//--------------------------------------------------------------
// Resets the containers
void RotationSolver::Reset()
{
}

//--------------------------------------------------------------
// Added on June 3, 2022
// Evaluates the moment of inertia for the neutron star.
// Also stores the omega_bar profile for use by the second-order solver.
void RotationSolver::FindNMomInertia()
{
	// double r_min = nstar_ptr->Profile().GetRadius()->operator[](0);
	// double r_surface = nstar_ptr->Profile().GetRadius()->operator[](-1);

	auto *R = nstar_ptr->Profile().GetRadius();
	const std::size_t n = nstar_ptr->Size();

	double r_surface = R->operator[](-1);

	// Find first strictly-positive grid radius (or fall back to epsilon)
	std::size_t i0 = 0;
	while (i0 < n && !(R->operator[](i0) > 0.0))
		++i0;

	double r0 = (i0 < n) ? R->operator[](i0) : kR_EPS_KM;
	if (r0 < kR_EPS_KM)
		r0 = kR_EPS_KM;
	double r = r0;

	auto *P = nstar_ptr->prof_.GetPressure();
	auto *E = nstar_ptr->prof_.GetEnergyDensity();
	auto *M = nstar_ptr->prof_.GetMass();

	// If these are not std::vector<double>, adapt this call to your container;
	// the key point is: store references/pointers to contiguous arrays.
	SetFastProfilePtrs_(*R, *P, *E, *M);

	// Arbitrary internal numerical normalization, NOT a physical spin (ADR-0006 Q2/P4). The
	// equation is linear and homogeneous, so this choice cancels from every quantity the public
	// API exposes. `seed_omega_bar_` defaults to the historical 5e-3 and is reachable only
	// through `RotationSolverTestSeam`, so seed invariance can be proved without the seed
	// becoming scientific API.
	init_omega_bar = seed_omega_bar_;
	// double r = r_min;

	double y[2];
	y[0] = init_omega_bar;
	y[1] = 0;

	double ang_mom_J, ang_vel_Omega, mom_inertia;

	gsl_odeiv2_system fast_ode_sys = {RotationSolver::ODE_N_Fast, nullptr, 2, this};

	gsl_odeiv2_driver *fast_driver = gsl_odeiv2_driver_alloc_y_new(&fast_ode_sys, gsl_odeiv2_step_rk8pd,
																   1.e-1, 1.e-10, 1.e-10);

	// Prepare storage for omega_bar profile
	stored_omega_bar_ = Zaki::Vector::DataColumn("omega_bar", n, 0.0);
	stored_domega_bar_ = Zaki::Vector::DataColumn("domega_bar", n, 0.0);

	// Records whether every grid point was integrated. Used ONLY to decide whether the
	// seed-free response below may be published; the extraction of J, Omega and I is
	// deliberately untouched, so `SeqPoint::I` is unaffected either way.
	bool solve_complete = true;

	// Radius loop inside the core
	for (size_t i = 0; i < n; i++)
	{
		if (GSL_SUCCESS !=
			gsl_odeiv2_driver_apply(fast_driver, &r, R->operator[](i), y))
		{
			solve_complete = false;
			break;
		}

		// Store omega_bar profile at each grid point
		stored_omega_bar_[i] = y[0];
		stored_domega_bar_[i] = y[1];
	}
	// ++++++++++++++++++++++++++++++++++++++++++++++++++

	// Total angular momentum and velocity:
	ang_mom_J = pow(r_surface, 4) * y[1] / 6.;
	ang_vel_Omega = y[0] + r_surface * y[1] / 3.;
	mom_inertia = ang_mom_J / ang_vel_Omega;

	nstar_ptr->MomI = mom_inertia;

	++first_order_solve_count_;

	// ---- Seed-free first-order response (ADR-0006 Q4) -------------------------------------
	// Dividing by Omega_raw removes the arbitrary seed ANALYTICALLY, not approximately: the
	// equation is linear and homogeneous, so omega_bar_raw, domega_bar_raw and Omega_raw all
	// carry the same factor A and every ratio below is a property of the star alone.
	//
	// Nothing published here is a physical spin. Per ADR-0006's binding clarification, an
	// NStar that was never given an explicit angular velocity does NOT acquire an implicit
	// one; a physical Omega, J and omega_bar(r) exist only after a caller supplies an
	// AngularVelocity to HartleFirstOrderResponse::At().
	//
	// A `hartle_result_.Omega/J/I/omega_bar/domega_bar` assignment lived here until ADR-0006.
	// It published Omega_raw behind an accessor annotated `[s^-1]` while storing km^-1, i.e. it
	// exposed the arbitrary seed as physical data. ADR-0006 removed those fields, and ADR-0007
	// retired the rest of that struct with the candidate it belonged to (Phase 4C-I1).
	first_order_response_ = HartleFirstOrderResponse{};
	if (solve_complete && n >= 2 && std::isfinite(ang_vel_Omega) && ang_vel_Omega != 0.0 &&
		std::isfinite(mom_inertia))
	{
		first_order_response_.omega_bar_over_Omega =
			Zaki::Vector::DataColumn("omega_bar_over_Omega", n, 0.0);
		first_order_response_.domega_bar_over_Omega_dr =
			Zaki::Vector::DataColumn("domega_bar_over_Omega_dr", n, 0.0);

		for (size_t i = 0; i < n; i++)
		{
			// Division, not multiplication by a precomputed reciprocal: one correctly-rounded
			// operation per node instead of two.
			first_order_response_.omega_bar_over_Omega[static_cast<int>(i)] =
				stored_omega_bar_[static_cast<int>(i)] / ang_vel_Omega;
			first_order_response_.domega_bar_over_Omega_dr[static_cast<int>(i)] =
				stored_domega_bar_[static_cast<int>(i)] / ang_vel_Omega;
		}

		first_order_response_.I = mom_inertia;
		first_order_response_.r_grid = nstar_ptr->Profile().GetRadius();
		first_order_response_.valid = true;
	}

	// The O(Omega^2) monopole response is NOT computed here. Ordinary star construction runs
	// this function, and ADR-0007's implementation contract keeps the second-order integration
	// off that path: it is work no existing workflow needs. Call
	// `RotationSolver::ComputeMonopoleResponse()` (or `NStar::ComputeHartleMonopoleResponse()`)
	// explicitly. Any previously cached response is left alone; its recorded profile version
	// is what makes it detectably stale.

	gsl_odeiv2_driver_free(fast_driver);
}

// //--------------------------------------------------------------
// // Added on Apr 22, 2022
// void RotationSolver::FindMixedMomInertia()
// {
	// Resolution of the solver (division of the radial distance)
	// unsigned int solver_res = 5000 ;

	// 	// double r_min = mixedstar_ptr->core_region.Min();
	// 	//  double r_mantle = mixedstar_ptr->mantle_region.Min() ;
	// 	//  double r_max = mixedstar_ptr->mantle_region.Max() ;

	// 	double r_min_raw = mixedstar_ptr->core_region.Min();
	// 	double r0 = SafeR0(r_min_raw);
	// 	double r = r0;

	// 	init_omega_bar = 5e-3;
	// 	// double r = r_min;

	// 	double y[2];
	// 	y[0] = init_omega_bar;
	// 	y[1] = 0;

	// 	// -------------------------------------------------
	// 	// Added on Apr 20, 2022
	// 	double r_surface = mixedstar_ptr->mantle_region.Max();

	// 	// Added on March 23, 2022 to fix the bug
	// 	// in J, Omega, I calculations.
	// 	//  double log_r_star = log10(r_surface) ;

	// 	// Total angular momentum and velocity:
	// 	double ang_mom_J, ang_vel_Omega, mom_inertia;
	// 	// -------------------------------------------------

	// 	// double step = (r_surface - r_min) / radial_res ;

	// 	// One extra point for the exact value at r = R_star.
	// 	// bool surface_reached = false ;

	// 	// ---------------------------------------------------------------------
	// 	//                         RADIUS LOOP BEGINS
	// 	// ---------------------------------------------------------------------
	// 	// for (double r_i = r_min;  r_i < r_surface; r_i += step)
	// 	// {
	// 	//   if ( GSL_SUCCESS !=
	// 	//       gsl_odeiv2_driver_apply (tmp_driver, &r, r_i, y) )
	// 	//     break ;
	// 	// }
	// 	// ............  Radius Loop Ends .........................

	// 	// gsl_odeiv2_driver_apply (tmp_driver, &r, r_surface, y) ;
	// 	// ++++++++++++++++++++++++++++++++++++++++++++++++++

	// 	// ++++++++++++++++++++++++++++++++++++++++++++++++++
	// 	// New method:
	// 	gsl_odeiv2_system fast_ode_sys = {RotationSolver::ODE_Mixed_Fast, nullptr, 2, this};

	// 	gsl_odeiv2_driver *fast_driver = gsl_odeiv2_driver_alloc_y_new(&fast_ode_sys, gsl_odeiv2_step_rk8pd,
	// 																   1.e-1, 1.e-10, 1.e-10);

	// 	if (mixedstar_ptr->dark_core)
	// 	{
	// 		auto &Rad_dar = mixedstar_ptr->ds_dar[0];
	// 		const size_t N_dar = Rad_dar.Size();

	// 		size_t i0 = 0;
	// 		while (i0 < N_dar && !(Rad_dar[i0] > r0))
	// 			++i0;
	// 		// Loop inside the core
	// 		for (size_t i = i0; i < N_dar; i++)
	// 		{
	// 			// fast_p_v = mixedstar_ptr->ds_vis[3][i];
	// 			// fast_p_d = mixedstar_ptr->ds_dar[3][i];
	// 			// fast_e_v = mixedstar_ptr->ds_vis[4][i];
	// 			// fast_e_d = mixedstar_ptr->ds_dar[4][i];
	// 			// fast_m_tot = mixedstar_ptr->ds_vis[1][i] + mixedstar_ptr->ds_dar[1][i];

	// 			if (GSL_SUCCESS !=
	// 				gsl_odeiv2_driver_apply(fast_driver, &r, Rad_dar[i], y))
	// 				break;
	// 		}

	// 		auto &Rad_vis = mixedstar_ptr->ds_vis[0];
	// 		const size_t N_vis = Rad_vis.Size();
	// 		// Loop inside the mantle
	// 		for (size_t i = N_dar; i < N_vis; i++)
	// 		{
	// 			// fast_p_v = mixedstar_ptr->ds_vis[3][i];
	// 			// fast_p_d = 0;
	// 			// fast_e_v = mixedstar_ptr->ds_vis[4][i];
	// 			// fast_e_d = 0;
	// 			// fast_m_tot = mixedstar_ptr->ds_vis[1][i];

	// 			if (GSL_SUCCESS !=
	// 				gsl_odeiv2_driver_apply(fast_driver, &r, Rad_vis[i], y))
	// 				break;
	// 		}
	// 	}
	// 	// ++++++++++++++++++++++++++++++++++++++++++++++++++
	// 	else // We have a dark mantle
	// 	{
	// 		auto &Rad_vis = mixedstar_ptr->ds_vis[0];
	// 		const size_t N_vis = Rad_vis.Size();

	// 		size_t i0 = 0;
	// 		while (i0 < N_vis && !(Rad_vis[i0] > r0))
	// 			++i0;
	// 		// Loop inside the core
	// 		for (size_t i = i0; i < N_vis; i++)
	// 		{
	// 			// fast_p_v = mixedstar_ptr->ds_vis[3][i];
	// 			// fast_p_d = mixedstar_ptr->ds_dar[3][i];
	// 			// fast_e_v = mixedstar_ptr->ds_vis[4][i];
	// 			// fast_e_d = mixedstar_ptr->ds_dar[4][i];
	// 			// fast_m_tot = mixedstar_ptr->ds_vis[1][i] + mixedstar_ptr->ds_dar[1][i];

	// 			if (GSL_SUCCESS !=
	// 				gsl_odeiv2_driver_apply(fast_driver, &r, Rad_vis[i], y))
	// 				break;
	// 		}

	// 		auto &Rad_dar = mixedstar_ptr->ds_dar[0];
	// 		const size_t N_dar = Rad_dar.Size();
	// 		// Loop inside the mantle
	// 		for (size_t i = N_vis; i < N_dar; i++)
	// 		{
	// 			// fast_p_v = 0;
	// 			// fast_p_d = mixedstar_ptr->ds_dar[3][i];
	// 			// fast_e_v = 0;
	// 			// fast_e_d = mixedstar_ptr->ds_dar[4][i];
	// 			// fast_m_tot = mixedstar_ptr->ds_dar[1][i];

	// 			if (GSL_SUCCESS !=
	// 				gsl_odeiv2_driver_apply(fast_driver, &r, Rad_dar[i], y))
	// 				break;
	// 		}
	// 	}
	// 	// ++++++++++++++++++++++++++++++++++++++++++++++++++

	// 	// Total angular momentum and velocity:
	// 	ang_mom_J = pow(r_surface, 4) * y[1] / 6.;
	// 	ang_vel_Omega = y[0] + r_surface * y[1] / 3.;
	// 	mom_inertia = ang_mom_J / ang_vel_Omega;

	// 	mixedstar_ptr->MomI = mom_inertia;

	// 	gsl_odeiv2_driver_free(fast_driver);
	// }
	// // ------------------------------------------------------------
	//--------------------------------------------------------------
	// Added on Jan 7, 2026
	void RotationSolver::FindMixedMomInertia()
	{
		PROFILE_FUNCTION();

		// Defensive
		if (!mixedstar_ptr)
		{
			Z_LOG_ERROR("Mixedstar_ptr is null.");
			return;
		}

		// Master-grid totals must be ready (constructed in MixedStar::SurfaceIsReached or import-ctor).
		// If you added a boolean like totals_ready_, use it here; otherwise rely on Size() checks.
		// if (!mixedstar_ptr->totals_ready_) { ... }

		// ---- Master grid + totals (all on the same radius grid) ----
		auto Rdc = mixedstar_ptr->GetRadius_Master();	   // r_master_dc
		auto Pdc = mixedstar_ptr->GetPress_Total_Master(); // pre_tot_dc
		auto Edc = mixedstar_ptr->GetEps_Total_Master();   // eps_tot_dc
		auto Mdc = mixedstar_ptr->GetMass_Total_Master();  // mass_tot_dc (master-grid aligned)

		if (!mixedstar_ptr->HasTotalMasterProfiles())
		{
			Z_LOG_ERROR("Missing master-grid totals.");
			return;
		}
		if (Rdc.Size() < 2 || Pdc.Size() < 2 || Edc.Size() < 2 || Mdc.Size() < 2)
		{
			Z_LOG_ERROR("Master-grid totals are too small.");
			return;
		}

		// Ensure consistent lengths (important for interpolation correctness).
		const std::size_t N = Rdc.Size();
		if (Pdc.Size() != N || Edc.Size() != N || Mdc.Size() != N)
		{
			Z_LOG_ERROR("Inconsistent master-grid total sizes.");
			return;
		}

		// ---- Robust center handling ----
		const double r_min_raw = mixedstar_ptr->core_region.Min();
		const double r0 = SafeR0(r_min_raw);
		double r = r0;

		// Surface is last point of master grid (more consistent than mantle_region.Max()).
		const double r_surface = Rdc.Values().back();

		// ---- Initial conditions (regular at center) ----
		init_omega_bar = 5e-3;
		double y[2];
		y[0] = init_omega_bar; // \bar{\omega}(0)
		y[1] = 0.0;			   // \bar{\omega}'(0)

		// ---- Point the fast mixed interpolation at master-grid totals ----
		SetFastMixedPtrs_(Rdc, Pdc, Edc, Mdc);

		// ---- GSL system/driver ----
		gsl_odeiv2_system sys = {RotationSolver::ODE_Mixed_Fast, nullptr, 2, this};

		gsl_odeiv2_driver *drv = gsl_odeiv2_driver_alloc_y_new(
			&sys, gsl_odeiv2_step_rk8pd,
			1.e-1,	// initial step guess
			1.e-10, // abs tol
			1.e-10	// rel tol
		);

		// ---- Integrate along the master grid (single pass) ----
		// Start from the first grid point strictly greater than r0.
		std::size_t i0 = 0;
		while (i0 < N && !((Rdc)[i0] > r0))
			++i0;

		for (std::size_t i = i0; i < N; ++i)
		{
			const double r_i = (Rdc)[i];
			if (GSL_SUCCESS != gsl_odeiv2_driver_apply(drv, &r, r_i, y))
				break;
		}

		// Ensure we hit the surface exactly (helps J/Omega stability).
		if (r < r_surface)
		{
			(void)gsl_odeiv2_driver_apply(drv, &r, r_surface, y);
		}

		// ---- Compute J, Omega, I ----
		const double ang_mom_J = std::pow(r_surface, 4) * y[1] / 6.0;
		const double ang_vel_Omega = y[0] + r_surface * y[1] / 3.0;
		const double mom_inertia = ang_mom_J / ang_vel_Omega;

		mixedstar_ptr->MomI = mom_inertia;

		gsl_odeiv2_driver_free(drv);
	}
	// ------------------------------------------------------------

	//--------------------------------------------------------------
	const HartleFirstOrderResponse &RotationSolver::FirstOrderResponse() const
	{
		return first_order_response_;
	}

	//--------------------------------------------------------------
	// Materializes the first-order solution at an explicitly requested physical angular
	// velocity (ADR-0006 P3). This is a SCALING of the already-computed seed-free response,
	// never a new ODE solve: the first-order problem is linear and homogeneous, so
	//
	//     omega_bar_phys(r) = [omega_bar/Omega](r) * Omega_geom ,
	//     domega_bar_phys(r) = [domega_bar/(Omega dr)](r) * Omega_geom ,
	//
	// and the exterior matching omega_bar = Omega - 2J/r^3 then gives J from the scaled
	// surface derivative. `I` is carried through unchanged rather than recomputed as J/Omega,
	// which would be 0/0 at zero spin (ADR-0006 P5).
	PhysicalFirstOrderRotation HartleFirstOrderResponse::At(AngularVelocity omega) const
	{
		const int n = static_cast<int>(omega_bar_over_Omega.Size());

		if (!valid || r_grid == nullptr || n < 2 ||
			static_cast<int>(domega_bar_over_Omega_dr.Size()) != n ||
			static_cast<int>(r_grid->Size()) != n)
		{
			throw std::runtime_error(
				"CompactStar::HartleFirstOrderResponse::At: no valid first-order response is "
				"available for this star, so no physical rotation can be materialized. A "
				"response is produced by RotationSolver::FindNMomInertia() on a star with a "
				"complete radial profile. Governed by ADR-0006.");
		}

		// The single physical -> geometric conversion (ADR-0006 P2), owned by AngularVelocity.
		const double Omega_geom = omega.GeomKmInverse();

		PhysicalFirstOrderRotation out;
		out.Omega_geom = Omega_geom;
		out.I = I;
		out.r_grid = r_grid;
		out.omega_bar = Zaki::Vector::DataColumn("omega_bar", static_cast<size_t>(n), 0.0);
		out.domega_bar = Zaki::Vector::DataColumn("domega_bar", static_cast<size_t>(n), 0.0);

		for (int i = 0; i < n; i++)
		{
			out.omega_bar[i] = omega_bar_over_Omega[i] * Omega_geom;
			out.domega_bar[i] = domega_bar_over_Omega_dr[i] * Omega_geom;
		}

		// Exterior matching at the surface, applied to the SCALED derivative. Deriving J from
		// the published array rather than from `I * Omega_geom` keeps the relation J = I Omega
		// a genuine consistency check on the arrays a consumer actually receives.
		const double r_surface = (*r_grid)[-1];
		out.J = pow(r_surface, 4) * out.domega_bar[-1] / 6.;

		out.valid = true;
		return out;
	}

	//==============================================================
	//   O(Omega^2) MONOPOLE (l = 0) — GOVERNED BY ADR-0007 (ACCEPTED 2026-09-02)
	//
	//   This section REPLACES the candidate of commit 675b4a9 (ODE_Hartle2_N_Fast,
	//   SolveHartle2_N, HartleResult), which Phase 4C-G adjudicated invalid against the
	//   primary source and which was deleted in the same commit that added this code
	//   (GOVERNANCE.md 3.1, authorized by ADR-0007 SS9). Nothing here is derived from it: every
	//   term below is transcribed from ADR-0007 P2, whose provenance is the journal scan of
	//   Hartle (1967) ApJ 150, 1005 -- eqs. (97) and (100) -- recorded term by term in
	//   docs/validation/PHASE4C_HARTLE2_DERIVATION.md SS6-SS7.
	//
	//   NOT YET INDEPENDENTLY VALIDATED. Conformance to an accepted contract is not physical
	//   verification; that is Phase 4D (INV-08). No baseline of these numbers may be created
	//   before it.
	//==============================================================

	//--------------------------------------------------------------
	bool RotationSolver::MonopoleBackground_::Complete() const noexcept
	{
		const Zaki::Vector::DataColumn *cols[] = {r, p, e, m, nu, nup, dedp, s, sp};
		for (const auto *c : cols)
			if (c == nullptr)
				return false;

		const std::size_t n = r->Size();
		if (n < 2)
			return false;
		for (const auto *c : cols)
			if (c->Size() != n)
				return false;

		return true;
	}

	//--------------------------------------------------------------
	// Every background input, sampled at the radius the ODE driver is ACTUALLY asking about.
	//
	// The retired candidate assigned per-node scalars before each gsl_odeiv2_driver_apply and
	// let its right-hand side use them at whatever internal radius GSL chose, so a node value
	// stood in for an inter-node one. That is scientifically invalid and is not repeated here.
	//
	// One shared bracket index serves all eight quantities, so they always refer to the same
	// radial interval, and the interpolation is linear -- the order the profile itself carries
	// (INV-13). The bracket walk mirrors EvalFastPEM_, which the validated first-order solver
	// already uses for exactly this purpose.
	RotationSolver::MonopoleBackground_::Sample
	RotationSolver::MonopoleBackground_::At(double r_km) const
	{
		Sample out;

		const auto &R = *r;
		const std::size_t n = R.Size();

		auto gather = [&](std::size_t i) {
			out.p = (*p)[static_cast<int>(i)];
			out.e = (*e)[static_cast<int>(i)];
			out.m = (*m)[static_cast<int>(i)];
			out.nu = (*nu)[static_cast<int>(i)];
			out.nup = (*nup)[static_cast<int>(i)];
			out.dedp = (*dedp)[static_cast<int>(i)];
			out.s = (*s)[static_cast<int>(i)];
			out.sp = (*sp)[static_cast<int>(i)];
		};

		// Outside the tabulated range the driver can only be at an endpoint of its own
		// integration interval; clamping there is the same convention EvalFastPEM_ uses.
		if (r_km <= R.Values().front())
		{
			gather(0);
			return out;
		}
		if (r_km >= R.Values().back())
		{
			gather(n - 1);
			return out;
		}

		std::size_t idx = k;
		if (idx >= n - 1)
			idx = n - 2;
		while (idx + 1 < n && r_km > R[static_cast<int>(idx + 1)])
			++idx;
		while (idx > 0 && r_km < R[static_cast<int>(idx)])
			--idx;
		k = idx;

		const double r0 = R[static_cast<int>(idx)];
		const double r1 = R[static_cast<int>(idx + 1)];
		const double t = (r1 > r0) ? (r_km - r0) / (r1 - r0) : 0.0;

		auto lerp = [&](const Zaki::Vector::DataColumn &c) {
			const double a = c[static_cast<int>(idx)];
			const double b = c[static_cast<int>(idx + 1)];
			return a + t * (b - a);
		};

		out.p = lerp(*p);
		out.e = lerp(*e);
		out.m = lerp(*m);
		out.nu = lerp(*nu);
		out.nup = lerp(*nup);
		out.dedp = lerp(*dedp);
		out.s = lerp(*s);
		out.sp = lerp(*sp);
		return out;
	}

	//--------------------------------------------------------------
	// ADR-0007 P2, verbatim, in the normalized variables
	//
	//     mhat = m0 / Omega_geom^2   [km^3]      phat = p0* / Omega_geom^2   [km^2]
	//
	// driven by the VERIFIED seed-free first-order response s = omega_bar/Omega [1] and
	// s' = [d omega_bar/dr]/Omega [km^-1]. The raw stored_omega_bar_ / stored_domega_bar_
	// arrays are deliberately NOT reachable from here: consuming them would reintroduce the
	// arbitrary internal seed that ADR-0006 P9 forbids in any second-order product.
	//
	//   dmhat/dr = 4 pi r^2 (eps+p) (deps/dp) phat                          [ADR-0007 P2, term 1]
	//            + (1/12) r^4 exp(-2 nu) (1 - 2m/r) s'^2                    [            term 2]
	//            + (8 pi/3) r^4 (eps+p) exp(-2 nu) s^2                      [            term 3]
	//
	//   dphat/dr = - mhat (1 + 8 pi r^2 p)/(r - 2m)^2                       [            term 4]
	//              - 4 pi (eps+p) r^2 phat/(r - 2m)                         [            term 5]
	//              + (1/12) r^3 exp(-2 nu) s'^2                             [            term 6]
	//              + (2/3) r exp(-2 nu) s (s + r s' - r nu' s)              [            term 7]
	//
	// Term 7 is the analytically expanded form of Hartle's (1/3) d/dr[ r^3 j^2 omegabar^2 /
	// (r - 2m) ]; the expansion is derived in PHASE4C_HARTLE2_DERIVATION.md SS7 and is kept in
	// this form so no numerical differentiation of the source appears in the right-hand side.
	//
	// Dimensions, checked term by term: every term of f[0] is km^2 and every term of f[1] is
	// km (PHASE4C_HARTLE2_DERIVATION.md SS6.3). Deliberately NOT algebraically compressed --
	// each term must stay visibly traceable to the contract it came from.
	int RotationSolver::ODE_HartleMonopole_(double r, const double y[], double f[], void *params)
	{
		RotationSolver *rot = static_cast<RotationSolver *>(params);

		// Geometry::MetricDenominator fails closed by throwing (ADR-0004 SS0-Q3). A C++
		// exception must not cross the GSL C callback boundary, so the driver is told instead.
		try
		{
			const auto b = rot->monopole_bg_.At(r);
			const double b_slope = rot->monopole_bg_.eps_slope; // ADR-0008 Q3, this segment only

			const double mhat = y[0];
			const double phat = y[1];

			// ADR-0004: the metric factor has ONE mathematical owner. No clamp, no 1e-15.
			const double D = CompactStar::Geometry::MetricDenominator(r, b.m); // 1 - 2m/r
			const double r_2m = r * D;										  // r - 2m

			const double e_m2nu = std::exp(-2.0 * b.nu);
			const double eps_p = b.e + b.p;

			const double r2 = r * r;
			const double r3 = r2 * r;
			const double r4 = r3 * r;

			// ---- term 1: the EOS energy-density MEASURE (ADR-0008 Q1/Q3) ----------------
			// Hartle's l = 0 field equation (H67 eq. 93) carries the matter source as
			// -8 pi xi dE/dR, i.e. the Lebesgue-Stieltjes measure
			//
			//     dm0_hat|_EOS = -4 pi r^2 xi0_hat d(eps) ,      xi0_hat = p0*_hat / nu' .
			//
			// His eq. (97) rewrites it as 4 pi r^2 (E+P)(dE/dP) p0* through (88) and the TOV
			// relation -- an identity only where eps is differentiable AND the background
			// satisfies dp/dr = -(eps+p) nu' pointwise. On a tabulated background neither holds
			// inside a feature narrower than the node spacing, and a nodal dE/dP column loses
			// whatever energy-density variation falls between its samples: that is the defect
			// Phase 4D measured (PHASE4D_MONOPOLE_VALIDATION.md section 14) and ADR-0008
			// corrects. Within one governed segment the measure density is the segment secant
			// installed in `eps_slope`, so the segment's TOTAL contribution is exactly its
			// Delta eps. `dedp` is NOT read here (ADR-0008 Q8: it remains the centre-series
			// authority only) and there is exactly one active EOS mass source.
			double term1 = 0.0;
			if (b_slope != 0.0)
			{
				// xi0_hat = p0*_hat/nu'. Both vanish at the centre (p0* ~ r^2, nu' ~ r), so the
				// ratio is finite there; where the division is ill-conditioned the SAME
				// regular-centre limit the derived xi0 column uses is applied, never an
				// epsilon regularization and never a fabricated zero (ADR-0008 §11 of the
				// implementation record; PHASE4C_HARTLE2_DERIVATION.md 9.2).
				double xi = phat / b.nup;
				if (!std::isfinite(xi))
				{
					const double denom = 4.0 * M_PI * (b.e + 3.0 * b.p);
					if (!(denom != 0.0) || !std::isfinite(rot->monopole_bg_.centre_xi_num))
						return GSL_EBADFUNC; // fail closed
					xi = rot->monopole_bg_.centre_xi_num * r / denom;
					if (!std::isfinite(xi))
						return GSL_EBADFUNC;
				}
				term1 = -4.0 * M_PI * r2 * xi * b_slope;
			}

			f[0] = term1											 // term 1 (ADR-0008 measure)
				   + (1.0 / 12.0) * r4 * e_m2nu * D * b.sp * b.sp	 // term 2
				   + (8.0 * M_PI / 3.0) * r4 * eps_p * e_m2nu * b.s * b.s; // term 3

			f[1] = -mhat * (1.0 + 8.0 * M_PI * r2 * b.p) / (r_2m * r_2m)  // term 4
				   - 4.0 * M_PI * eps_p * r2 * phat / r_2m				  // term 5
				   + (1.0 / 12.0) * r3 * e_m2nu * b.sp * b.sp			  // term 6
				   + (2.0 / 3.0) * r * e_m2nu * b.s *
						 (b.s + r * b.sp - r * b.nup * b.s);			  // term 7

			if (!std::isfinite(f[0]) || !std::isfinite(f[1]))
				return GSL_EBADFUNC;
		}
		catch (const std::exception &)
		{
			return GSL_EBADFUNC;
		}

		return GSL_SUCCESS;
	}

	//--------------------------------------------------------------
	const HartleMonopoleResponse *RotationSolver::MonopoleResponse() const
	{
		if (nstar_ptr == nullptr)
			return nullptr;

		const auto &prof = nstar_ptr->Profile();
		if (!monopole_response_.MatchesSource(static_cast<const void *>(&prof), prof.Version()))
			return nullptr; // absent or stale -- never returned as current

		return &monopole_response_;
	}

	//--------------------------------------------------------------
	// Computes the governed fixed-eps_c, seed-free O(Omega^2) monopole response.
	//
	// PUBLICATION IS ATOMIC. Everything is built in local storage and `monopole_response_` is
	// assigned exactly once, at the end, only if every acceptance test passed. On any failure
	// the cached response is cleared, so a stale or partial result can never be mistaken for a
	// current one.
	bool RotationSolver::ComputeMonopoleResponse()
	{
		if (nstar_ptr == nullptr)
		{
			Z_LOG_ERROR("ComputeMonopoleResponse: no NStar is attached.");
			monopole_response_ = HartleMonopoleResponse{};
			return false;
		}

		const auto &prof = nstar_ptr->Profile();
		const void *prof_id = static_cast<const void *>(&prof);
		const std::uint64_t prof_ver = prof.Version();

		// Already current for this exact profile state: reuse, and do NOT integrate again.
		if (monopole_response_.MatchesSource(prof_id, prof_ver))
			return true;

		monopole_response_ = HartleMonopoleResponse{};

		// ---- inputs -------------------------------------------------------------------
		const auto *R = prof.GetRadius();
		const auto *P = prof.GetPressure();
		const auto *E = prof.GetEnergyDensity();
		const auto *M = prof.GetMass();
		const auto *NU = prof.GetMetricNu();
		const auto *NUP = prof.GetMetricNuPrime();

		if (R == nullptr || P == nullptr || E == nullptr || M == nullptr || NU == nullptr ||
			NUP == nullptr)
		{
			Z_LOG_ERROR("ComputeMonopoleResponse: the profile is missing a required column.");
			return false;
		}

		// ADR-0007 P5: the ONLY admissible source of d(eps)/dp is the EOS authority, carried
		// on the profile by Phase 4C-I0. No profile finite difference, no fallback value, no
		// EOS object reached from the rotation solver. Absent data fails the computation.
		if (!prof.HasEosDEdP() || prof.GetEosDEdP() == nullptr)
		{
			Z_LOG_ERROR(
				"ComputeMonopoleResponse: this star carries no authoritative d(eps)/dp, so the "
				"governed O(Omega^2) monopole system cannot be integrated. Build the star from "
				"an EOS-backed TOV solve, or supply the derivative explicitly on the TOVPoints "
				"for an analytic star. ADR-0007 P5 forbids substituting a profile finite "
				"difference.");
			return false;
		}
		const auto *DEDP = prof.GetEosDEdP();

		if (!first_order_response_.valid)
		{
			Z_LOG_ERROR("ComputeMonopoleResponse: no valid first-order response; the O(Omega^2) "
						"sources are built from omega_bar/Omega and its derivative.");
			return false;
		}

		const std::size_t n = R->Size();
		if (n < 2 || DEDP->Size() != n ||
			first_order_response_.omega_bar_over_Omega.Size() != n ||
			first_order_response_.domega_bar_over_Omega_dr.Size() != n)
		{
			Z_LOG_ERROR("ComputeMonopoleResponse: input column sizes disagree with the radial "
						"grid.");
			return false;
		}

		// ADR-0008 Q3: the governed measure is carried by the profile partition itself, so the
		// partition must be a partition: strictly increasing radii and finite energy-density
		// increments. A non-monotone or non-finite node makes the segment measure meaningless,
		// and the contract fails closed rather than integrating something else.
		for (std::size_t i = 0; i + 1 < n; ++i)
		{
			const int ii = static_cast<int>(i);
			const double dr = R->operator[](ii + 1) - R->operator[](ii);
			const double de = (*E)[ii + 1] - (*E)[ii];
			if (!(dr > 0.0) || !std::isfinite(dr) || !std::isfinite(de))
			{
				Z_LOG_ERROR("ComputeMonopoleResponse: the radial partition is not strictly "
							"increasing, or carries a non-finite energy-density increment, so "
							"the governed EOS measure (ADR-0008 Q3) is undefined.");
				return false;
			}
		}

		monopole_bg_.r = R;
		monopole_bg_.p = P;
		monopole_bg_.e = E;
		monopole_bg_.m = M;
		monopole_bg_.nu = NU;
		monopole_bg_.nup = NUP;
		monopole_bg_.dedp = DEDP;
		monopole_bg_.s = &first_order_response_.omega_bar_over_Omega;
		monopole_bg_.sp = &first_order_response_.domega_bar_over_Omega_dr;
		monopole_bg_.k = 0;
		monopole_bg_.eps_slope = 0.0;	 // installed per segment below (ADR-0008 Q3)
		monopole_bg_.centre_xi_num = 0.0; // set by the centre initialization

		if (!monopole_bg_.Complete())
		{
			Z_LOG_ERROR("ComputeMonopoleResponse: incomplete background workspace.");
			return false;
		}

		// ---- regular-centre start (ADR-0007 P4) ---------------------------------------
		// The first strictly-positive grid radius, matching the first-order convention
		// (INV-05). Literal {0,0} starts are what the retired candidate did; the governed
		// contract initializes from the regular series instead, which makes "fixed eps_c"
		// exact to rounding rather than to O((r0/R)^2).
		std::size_t i0 = 0;
		while (i0 < n && !(R->operator[](static_cast<int>(i0)) > 0.0))
			++i0;
		double r0 = (i0 < n) ? R->operator[](static_cast<int>(i0)) : kR_EPS_KM;
		if (r0 < kR_EPS_KM)
			r0 = kR_EPS_KM;

		double mhat0 = 0.0, phat0 = 0.0, j2_c = 0.0, s_c = 0.0;
		try
		{
			const auto b0 = monopole_bg_.At(r0);
			// j^2 = exp(-2 nu) (1 - 2m/r), Hartle's j = exp[-(nu+lambda)] in this convention.
			j2_c = std::exp(-2.0 * b0.nu) * CompactStar::Geometry::MetricDenominator(r0, b0.m);
			s_c = b0.s;

			phat0 = (1.0 / 3.0) * j2_c * s_c * s_c * r0 * r0;
			// ADR-0008 Q8: the regular-centre series keeps the POINTWISE EOS derivative of the
			// 4C-I0 authority. It is a local, well-resolved property of the central state, not
			// an integrated measure, and the measure contract does not touch it.
			mhat0 = (4.0 * M_PI / 15.0) * (b0.e + b0.p) * (b0.dedp + 2.0) * j2_c * s_c * s_c *
					std::pow(r0, 5);
			monopole_bg_.centre_xi_num = j2_c * s_c * s_c;
		}
		catch (const std::exception &ex)
		{
			Z_LOG_ERROR(std::string("ComputeMonopoleResponse: centre initialization failed: ") +
						ex.what());
			return false;
		}

		if (!std::isfinite(mhat0) || !std::isfinite(phat0))
		{
			Z_LOG_ERROR("ComputeMonopoleResponse: non-finite regular-centre initialization.");
			return false;
		}

		// ---- integrate -----------------------------------------------------------------
		// PROVISIONAL NUMERICAL IMPLEMENTATION -- the step type and tolerances mirror the
		// first-order driver because that is a reasonable engineering starting point, NOT
		// because the retired candidate used them. Their adequacy, and the convergence of the
		// result under radial refinement, are Phase 4D questions (ADR-0007 SS7 item 9).
		Zaki::Vector::DataColumn mhat_col("m0_over_Omega2", n, 0.0);
		Zaki::Vector::DataColumn phat_col("p0star_over_Omega2", n, 0.0);

		gsl_odeiv2_system sys = {RotationSolver::ODE_HartleMonopole_, nullptr, 2, this};
		gsl_odeiv2_driver *driver =
			gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, 1.e-1, 1.e-10, 1.e-10);

		double r = r0;
		double y[2] = {mhat0, phat0};
		bool ok = true;

		// ADR-0008 Q3: profile-node boundaries are MANDATORY integration boundaries. The driver
		// advances exactly one governed segment per call, and that segment's own energy-density
		// measure density is installed before the call, so no step -- and no adaptive substep,
		// since `gsl_odeiv2_driver_apply` integrates to its target and never past it -- can
		// carry one segment's Delta eps into another. This is what removes the node-placement
		// sensitivity the nodal-derivative source had.
		for (std::size_t i = 0; i < n; ++i)
		{
			const int ii = static_cast<int>(i);
			// Segment ending at node i; node 0 is the start of the integration (no advance on a
			// governed profile), and takes the first segment's measure for the degenerate case
			// where the start radius had to be raised off a non-positive first node.
			const int lo = (i == 0) ? 0 : ii - 1;
			const int hi = (i == 0) ? 1 : ii;
			monopole_bg_.eps_slope =
				((*E)[hi] - (*E)[lo]) / (R->operator[](hi) - R->operator[](lo));

			if (GSL_SUCCESS != gsl_odeiv2_driver_apply(driver, &r, R->operator[](ii), y))
			{
				ok = false;
				break;
			}
			mhat_col[ii] = y[0];
			phat_col[ii] = y[1];
		}
		monopole_bg_.eps_slope = 0.0;

		gsl_odeiv2_driver_free(driver);
		++monopole_solve_count_;

		if (!ok)
		{
			Z_LOG_ERROR("ComputeMonopoleResponse: the O(Omega^2) integration did not complete.");
			return false;
		}

		// ---- derived fields (ADR-0007 P1: p0* is integrated, the rest are derived) ------
		Zaki::Vector::DataColumn dp_col("delta_p0_over_Omega2", n, 0.0);
		Zaki::Vector::DataColumn xi_col("xi0_over_Omega2", n, 0.0);

		for (std::size_t i = 0; i < n; ++i)
		{
			const int ii = static_cast<int>(i);
			const double eps_p = (*E)[ii] + (*P)[ii];

			dp_col[ii] = eps_p * phat_col[ii];

			// xi0 = p0*/nu'. Both vanish at the centre (p0* ~ r^2, nu' ~ r), so the ratio is
			// finite there but the division is not always well conditioned on the first node.
			// Where it is not, the regular-centre series is used -- NOT a zero and not an
			// arbitrary epsilon, which would replace a real value with a fabricated one:
			//
			//     xi0_hat -> j_c^2 s_c^2 r / [4 pi (eps_c + 3 p_c)]
			//
			// (PHASE4C_HARTLE2_DERIVATION.md SS9.2).
			const double nup_i = (*NUP)[ii];
			double xi = phat_col[ii] / nup_i;
			if (!std::isfinite(xi))
			{
				const double denom = 4.0 * M_PI * ((*E)[ii] + 3.0 * (*P)[ii]);
				xi = (denom != 0.0)
						 ? j2_c * s_c * s_c * R->operator[](ii) / denom
						 : std::numeric_limits<double>::quiet_NaN();
			}
			xi_col[ii] = xi;
		}

		// ---- surface (ADR-0007 P6/P7) --------------------------------------------------
		// R_* is the LAST PROFILE NODE -- the production EOS-floor surface (INV-06). It is
		// deliberately not identified with the exact p = 0 radius; the resulting isobar
		// systematic is documented in PHASE4C_HARTLE2_DERIVATION.md SS11 and carried, not hidden.
		const int last = static_cast<int>(n) - 1;
		const double R_star = R->operator[](last);
		const double eps_star = (*E)[last];
		const double xi_star = xi_col[last];
		const double I_first = first_order_response_.I;

		// The surface mass shell is the TERMINAL ATOM of the same energy-density measure
		// (ADR-0008 Q7): the interior measure runs over [r0, R_*) segment by segment above, and
		// the remaining jump eps_* -> 0 at R_* is applied here exactly ONCE,
		//
		//     Delta m0_hat|_{R_*} = 4 pi R_*^2 (eps_* - 0) xi0_hat(R_*) ,
		//
		// the same jump operator a declared internal discontinuity would use. It is never also
		// folded into the last interior segment, whose measure is eps[n-1] - eps[n-2]. It is
		// negligible on an EOS-floor star and DOMINANT on a constant-density one -- so it is
		// computed, never assumed small (ADR-0007 P6).
		const double shell_hat = 4.0 * M_PI * R_star * R_star * eps_star * xi_star;

		// delta_M = m0(R_*) + shell + J^2/R_*^3, with J/Omega = I from the VERIFIED
		// first-order response. No raw J, no seed-dependent quantity, no arbitrary Omega.
		const double deltaM_hat =
			mhat_col[last] + shell_hat + (I_first * I_first) / (R_star * R_star * R_star);

		// ---- acceptance ---------------------------------------------------------------
		for (std::size_t i = 0; i < n; ++i)
		{
			const int ii = static_cast<int>(i);
			if (!std::isfinite(mhat_col[ii]) || !std::isfinite(phat_col[ii]) ||
				!std::isfinite(dp_col[ii]) || !std::isfinite(xi_col[ii]))
			{
				Z_LOG_ERROR("ComputeMonopoleResponse: a non-finite value reached a published "
							"field; the response is discarded.");
				return false;
			}
		}
		if (!std::isfinite(deltaM_hat) || !std::isfinite(shell_hat) || !std::isfinite(xi_star) ||
			!(R_star > 0.0))
		{
			Z_LOG_ERROR("ComputeMonopoleResponse: non-finite surface quantity; discarded.");
			return false;
		}

		// ---- publish, atomically -------------------------------------------------------
		HartleMonopoleResponse out;
		out.m0_over_Omega2 = mhat_col;
		out.p0star_over_Omega2 = phat_col;
		out.delta_p0_over_Omega2 = dp_col;
		out.xi0_over_Omega2 = xi_col;
		out.deltaM_over_Omega2 = deltaM_hat;
		out.surface_shell_mass_over_Omega2 = shell_hat;
		out.surface_xi0_over_Omega2 = xi_star;
		out.R_surface = R_star;
		out.I = I_first;
		out.r_grid = R;
		out.source_profile = prof_id;
		out.source_version = prof_ver;
		out.valid = true;

		monopole_response_ = out;
		return true;
	}

	//--------------------------------------------------------------
	// Materializes the monopole perturbation at an explicitly requested physical angular
	// velocity: Q = Q_hat * Omega_geom^2 (ADR-0007 P9). This is a SCALING, never a new ODE
	// solve -- the coefficients are properties of the star and do not depend on Omega.
	//
	// Zero spin gives exact zeros. Because every quantity is quadratic in Omega, +Omega and
	// -Omega give BIT-IDENTICAL results: (-x)*(-x) and x*x are the same IEEE-754 product.
	PhysicalHartleMonopole HartleMonopoleResponse::At(AngularVelocity omega) const
	{
		const int n = static_cast<int>(m0_over_Omega2.Size());

		if (!valid || r_grid == nullptr || n < 2 ||
			static_cast<int>(p0star_over_Omega2.Size()) != n ||
			static_cast<int>(delta_p0_over_Omega2.Size()) != n ||
			static_cast<int>(xi0_over_Omega2.Size()) != n ||
			static_cast<int>(r_grid->Size()) != n)
		{
			throw std::runtime_error(
				"CompactStar::HartleMonopoleResponse::At: no valid O(Omega^2) monopole response "
				"is available for this star. Call NStar::ComputeHartleMonopoleResponse() first; "
				"it fails closed when the star carries no authoritative d(eps)/dp. Governed by "
				"ADR-0007.");
		}

		// The single physical -> geometric conversion (ADR-0006 P2), owned by AngularVelocity.
		const double Omega_geom = omega.GeomKmInverse();
		const double q = Omega_geom * Omega_geom;

		PhysicalHartleMonopole out;
		out.Omega_geom = Omega_geom;
		out.m0 = Zaki::Vector::DataColumn("m0", static_cast<size_t>(n), 0.0);
		out.p0star = Zaki::Vector::DataColumn("p0star", static_cast<size_t>(n), 0.0);
		out.delta_p0 = Zaki::Vector::DataColumn("delta_p0", static_cast<size_t>(n), 0.0);
		out.xi0 = Zaki::Vector::DataColumn("xi0", static_cast<size_t>(n), 0.0);

		for (int i = 0; i < n; ++i)
		{
			out.m0[i] = m0_over_Omega2[i] * q;
			out.p0star[i] = p0star_over_Omega2[i] * q;
			out.delta_p0[i] = delta_p0_over_Omega2[i] * q;
			out.xi0[i] = xi0_over_Omega2[i] * q;
		}

		out.delta_M = deltaM_over_Omega2 * q;
		out.surface_shell_mass = surface_shell_mass_over_Omega2 * q;
		out.surface_xi0 = surface_xi0_over_Omega2 * q;
		out.r_grid = r_grid;
		out.source_profile = source_profile;
		out.source_version = source_version;
		out.valid = true;
		return out;
	}

	//==============================================================
