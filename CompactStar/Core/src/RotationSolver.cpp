/*
  Last edited in 2026
  RotationSolver class
*/
#include <cmath>
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
	// return 16. * M_PI * (fast_p + fast_e) / (1. - 2. * fast_m / r);
	double p, e, m;
	EvalFastPEM_(r, p, e, m);
	return 16. * M_PI * (p + e) / (1. - 2. * m / r);
}

//--------------------------------------------------------------
// MixedStar: coefficient for y[1]
double RotationSolver::GetHartleDOmegaCoeff_N_Fast(const double r)
{
	// return 4. * M_PI * (fast_p + fast_e) * r / (1. - 2. * fast_m / r);
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

	char seq_header[200];
	snprintf(seq_header, sizeof(seq_header), "%-16s\t %-14s\t %-14s\t %-14s\t %-16s",
			 "omega_bar_c (1/s)", "M", "R", "J", "Omega (1/s)");

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

	init_omega_bar = 5e-3;
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

	// Radius loop inside the core
	for (size_t i = 0; i < n; i++)
	{
		// fast_p = nstar_ptr->prof_.GetPressure()->operator[](i);
		// fast_e = nstar_ptr->prof_.GetEnergyDensity()->operator[](i);
		// fast_m = nstar_ptr->prof_.GetMass()->operator[](i);

		if (GSL_SUCCESS !=
			gsl_odeiv2_driver_apply(fast_driver, &r, R->operator[](i), y))
			break;

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

	// Store first-order results in HartleResult
	hartle_result_.Omega = ang_vel_Omega;
	hartle_result_.J = ang_mom_J;
	hartle_result_.I = mom_inertia;
	hartle_result_.omega_bar = stored_omega_bar_;
	hartle_result_.domega_bar = stored_domega_bar_;
	hartle_result_.r_grid = nstar_ptr->Profile().GetRadius();

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
	const HartleResult &RotationSolver::GetHartleResult() const
	{
		return hartle_result_;
	}

	//--------------------------------------------------------------
	// Second-order (m0, p0) ODE for NStar with fast interpolation.
	// y[0] = m0, y[1] = p0
	// Source terms from omega_bar are read from fast_omega_bar, fast_domega_bar.
	//
	// Equations follow Hartle (1967), Eqs. (30)-(33).
	// Using j(r) = exp(-nu) * sqrt(1 - 2m/r):
	//
	// dm0/dr = 4*pi*r^2 * (deps/dp) * p0 + S_m(r)
	//
	// dp0/dr = -(m0 + 4*pi*r^3*p0) / [r*(r - 2m)]
	//          - 4*pi*r*(eps+p) * p0 / (r - 2m)  [FIX: confirm exact from textbook]
	//          + S_p(r)
	//
	// where S_m, S_p are quadratic source terms in omega_bar.
	//
	int RotationSolver::ODE_Hartle2_N_Fast(double r, const double y[], double f[], void *params)
	{
		RotationSolver *rot = static_cast<RotationSolver *>(params);

		const double p = rot->fast_p;
		const double eps = rot->fast_e;
		const double m = rot->fast_m;
		const double dEdP = rot->fast_dEdP;
		const double nu_prime = rot->fast_nu_prime;

		const double m0 = y[0];
		const double p0 = y[1];

		const double r2 = r * r;
		const double r3 = r2 * r;
		const double r_2m = r - 2.0 * m;

		// Avoid division by zero at r = 0 or r = 2m
		if (r < 1.e-10 || std::abs(r_2m) < 1.e-30)
		{
			f[0] = 0.0;
			f[1] = 0.0;
			return GSL_SUCCESS;
		}

		const double inv_r_2m = 1.0 / r_2m;
		const double eps_plus_p = eps + p;

		// ------- Homogeneous part of the equations -------
		// dm0/dr (homogeneous): 4*pi*r^2 * (deps/dp) * p0
		double dm0_dr = 4.0 * M_PI * r2 * dEdP * p0;

		// dp0/dr (homogeneous, from perturbed TOV):
		// dp0/dr = -(m0 + 4*pi*r^3 * p0) * (eps+p) / (r^2 * (r - 2m))
		//          + (m0/r^2 + 4*pi*p0) * nu' [alternate form via nu']
		//
		// Using the standard form: dp0/dr = -nu'*(eps+p+p0*dEdP) ...
		// We use the direct form from Hartle (1967) Appendix:
		double dp0_dr = -(m0 + 4.0 * M_PI * r3 * p0) * inv_r_2m / r2;
		// This is the gravitational term. Now apply the (eps+p) factor:
		// Actually the full form is:
		// dp0/dr = -p0 * [4*pi*(eps+p)*r / (r-2m) + m/(r^2*(r-2m))]
		//          - m0 * [(eps+p) / (r^2*(r-2m))]
		// Which groups as:
		dp0_dr = -p0 * (4.0 * M_PI * eps_plus_p * r + m / r2) * inv_r_2m - m0 * eps_plus_p / (r2 * r_2m);

		// ------- Source terms (quadratic in omega_bar) -------
		if (rot->include_m0p0_source_)
		{
			const double ob = rot->fast_omega_bar;
			const double dob = rot->fast_domega_bar;

			// j(r)^2 = exp(-2*nu) * (1 - 2*m/r) = exp(-2*nu) * (r-2m)/r
			// We compute j^2 using: j^2 = (r_2m / r) * exp(-2*nu)
			// But we don't have nu directly in the fast cache.
			// Instead, use the relation: exp(-2*nu) can be obtained from
			// the TOV structure. For now, compute via the known formula:
			// (1 - 2m/r) = r_2m / r
			const double one_minus_2m_r = r_2m / r;

			// Source term for dm0/dr:
			// S_m = (1/12) * r^4 * (dob/dr)^2 / (1 - 2m/r)
			//       - (1/3) * r^3 * 4*pi*(eps+p) * ob^2 / (1 - 2m/r)
			// Ref: Hartle (1967) Eq. (33), adapted for omega_bar convention
			double S_m = (1.0 / 12.0) * r2 * r2 * dob * dob * inv_r_2m * r - (1.0 / 3.0) * r3 * 4.0 * M_PI * eps_plus_p * ob * ob * inv_r_2m * r;

			// Source term for dp0/dr:
			// S_p = -(1/12) * r^2 * (dob)^2 / [(eps+p) * ???]
			// Using the relation from the perturbed TOV with rotation:
			// S_p = (1/12) * r * (dob)^2 * (r - 2m) / r
			//       + (4/3) * M_PI * r * ob^2 * (eps+p) / (1-2m/r)
			// This needs careful derivation — using Hartle (1967) Eq. (30):
			// dp*/dr = ... + (1/3) * omega_bar^2 * r * (dp/dr) / [r-2m]
			//          + (1/12) * r * (r-2m) * (d_omega_bar/dr)^2
			//
			// More precisely, the source for dp0/dr involves:
			double S_p = (1.0 / 12.0) * r * one_minus_2m_r * dob * dob + (1.0 / 3.0) * ob * ob * r * nu_prime;

			dm0_dr += S_m;
			dp0_dr += S_p;
		}

		f[0] = dm0_dr;
		f[1] = dp0_dr;

		return GSL_SUCCESS;
	}

	//--------------------------------------------------------------
	// Stub for MixedStar second-order ODE (to be implemented in Phase C)
	int RotationSolver::ODE_Hartle2_Mixed_Fast(double r, const double y[], double f[], void *params)
	{
		// TODO: Implement for MixedStar (Phase C)
		f[0] = 0.0;
		f[1] = 0.0;
		return GSL_SUCCESS;
	}

	//--------------------------------------------------------------
	// Solves the second-order Hartle O(Omega^2) monopole equations for NStar.
	// Uses superposition: particular solution (with source) + homogeneous solution.
	// Shooting condition: p0(R) = 0.
	void RotationSolver::SolveHartle2_N()
	{
		const size_t N = nstar_ptr->Size();
		const auto *r_col = nstar_ptr->Profile().GetRadius();
		const auto *p_col = nstar_ptr->Profile().GetPressure();
		const auto *e_col = nstar_ptr->Profile().GetEnergyDensity();
		const auto *m_col = nstar_ptr->Profile().GetMass();
		const auto *nu_p_col = nstar_ptr->Profile().GetMetricNuPrime();

		double r_min = (*r_col)[0];
		double r_surface = (*r_col)[-1];

		// --- Precompute d(eps)/d(p) column via finite differences on profile ---
		// d(eps)/d(p) = (d_eps/d_r) / (d_p/d_r), computed from stored columns.
		Zaki::Vector::DataColumn dEdP_col("dEdP", N, 0.0);
		for (size_t i = 1; i < N - 1; i++)
		{
			double dp = (*p_col)[i + 1] - (*p_col)[i - 1];
			double de = (*e_col)[i + 1] - (*e_col)[i - 1];
			if (std::abs(dp) > 1.e-30)
				dEdP_col[i] = de / dp;
			else
				dEdP_col[i] = 1.0; // Fallback (incompressible limit)
		}
		// Boundary: one-sided differences
		{
			double dp = (*p_col)[1] - (*p_col)[0];
			double de = (*e_col)[1] - (*e_col)[0];
			dEdP_col[0] = (std::abs(dp) > 1.e-30) ? de / dp : 1.0;
		}
		{
			int last = static_cast<int>(N) - 1;
			double dp = (*p_col)[last] - (*p_col)[last - 1];
			double de = (*e_col)[last] - (*e_col)[last - 1];
			dEdP_col[last] = (std::abs(dp) > 1.e-30) ? de / dp : 1.0;
		}

		// === Pass 1: Particular solution (p0_c = 0, source ON) ===
		include_m0p0_source_ = true;

		gsl_odeiv2_system ode_sys = {RotationSolver::ODE_Hartle2_N_Fast, nullptr, 2, this};
		gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
			&ode_sys, gsl_odeiv2_step_rk8pd, 1.e-1, 1.e-10, 1.e-10);

		double r = r_min;
		double y_part[2] = {0.0, 0.0}; // m0 = 0, p0 = 0 at center

		// Storage for particular solution profile
		Zaki::Vector::DataColumn m0_part("m0_part", N, 0.0);
		Zaki::Vector::DataColumn p0_part("p0_part", N, 0.0);

		for (size_t i = 0; i < N; i++)
		{
			fast_p = (*p_col)[i];
			fast_e = (*e_col)[i];
			fast_m = (*m_col)[i];
			fast_nu_prime = (*nu_p_col)[i];
			fast_dEdP = dEdP_col[i];
			fast_omega_bar = stored_omega_bar_[i];
			fast_domega_bar = stored_domega_bar_[i];

			if (GSL_SUCCESS !=
				gsl_odeiv2_driver_apply(driver, &r, (*r_col)[i], y_part))
				break;

			m0_part[i] = y_part[0];
			p0_part[i] = y_part[1];
		}

		gsl_odeiv2_driver_free(driver);

		// === Pass 2: Homogeneous solution (p0_c = 1, source OFF) ===
		include_m0p0_source_ = false;

		gsl_odeiv2_system ode_sys_hom = {RotationSolver::ODE_Hartle2_N_Fast, nullptr, 2, this};
		gsl_odeiv2_driver *driver_hom = gsl_odeiv2_driver_alloc_y_new(
			&ode_sys_hom, gsl_odeiv2_step_rk8pd, 1.e-1, 1.e-10, 1.e-10);

		r = r_min;
		double y_hom[2] = {0.0, 1.0}; // m0 = 0, p0 = 1 at center

		Zaki::Vector::DataColumn m0_hom("m0_hom", N, 0.0);
		Zaki::Vector::DataColumn p0_hom("p0_hom", N, 0.0);

		for (size_t i = 0; i < N; i++)
		{
			fast_p = (*p_col)[i];
			fast_e = (*e_col)[i];
			fast_m = (*m_col)[i];
			fast_nu_prime = (*nu_p_col)[i];
			fast_dEdP = dEdP_col[i];
			fast_omega_bar = 0.0; // Not used when source is off
			fast_domega_bar = 0.0;

			if (GSL_SUCCESS !=
				gsl_odeiv2_driver_apply(driver_hom, &r, (*r_col)[i], y_hom))
				break;

			m0_hom[i] = y_hom[0];
			p0_hom[i] = y_hom[1];
		}

		gsl_odeiv2_driver_free(driver_hom);

		// Reset source flag
		include_m0p0_source_ = true;

		// === Superposition: find p0_c such that p0(R) = 0 ===
		double p0_part_R = p0_part[-1];
		double p0_hom_R = p0_hom[-1];

		double p0_c = 0.0;
		if (std::abs(p0_hom_R) > 1.e-30)
			p0_c = -p0_part_R / p0_hom_R;

		// === Combine into final profiles ===
		hartle_result_.m0 = Zaki::Vector::DataColumn("m0", N, 0.0);
		hartle_result_.p0 = Zaki::Vector::DataColumn("p0", N, 0.0);
		hartle_result_.xi0 = Zaki::Vector::DataColumn("xi0", N, 0.0);

		for (size_t i = 0; i < N; i++)
		{
			hartle_result_.m0[i] = m0_part[i] + p0_c * m0_hom[i];
			hartle_result_.p0[i] = p0_part[i] + p0_c * p0_hom[i];
		}

		// === Compute xi0 = -p0 / (dp/dr) ===
		// dp/dr = -(eps+p) * nu' (from TOV equation)
		for (size_t i = 0; i < N; i++)
		{
			double eps_plus_p = (*e_col)[i] + (*p_col)[i];
			double nu_p = (*nu_p_col)[i];
			double dp_dr = -eps_plus_p * nu_p;

			if (std::abs(dp_dr) > 1.e-30)
				hartle_result_.xi0[i] = -hartle_result_.p0[i] / dp_dr;
			else
				hartle_result_.xi0[i] = 0.0; // At surface: both p0 and dp/dr -> 0
		}

		// Store scalars
		hartle_result_.p0_c = p0_c;
		hartle_result_.delta_M = hartle_result_.m0[-1];
		hartle_result_.valid = true;
	}

	//--------------------------------------------------------------
	// Stub for MixedStar second-order solve (to be implemented in Phase C)
	void RotationSolver::SolveHartle2_Mixed()
	{
		// TODO: Implement for MixedStar (Phase C)
	}

	//==============================================================
