/*
  Last edited on Dec 29, 2025

  TOVSolver class
*/

// Creating directories
#include <filesystem>
#include <sys/stat.h>

#include <array>
#include <cmath>
#include <limits>

#include <gsl/gsl_const_cgsm.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_roots.h>

#include <Zaki/File/CSVIterator.hpp>
#include <Zaki/File/VecSaver.hpp>
#include <Zaki/Math/GSLFuncWrapper.hpp>
#include <Zaki/Physics/Constants.hpp>
#include <Zaki/Util/Instrumentor.hpp>

#include "CompactStar/Core/Analysis.hpp"
#include "CompactStar/Core/MixedStar.hpp"
#include "CompactStar/Core/NStar.hpp"
#include "CompactStar/Core/TOVSolver.hpp"

#define TOV_SOLVER_VERBOSE 1
#define DarkCore_Analysis 0

using namespace CompactStar::Core;

//==============================================================
//               Safe Spline Evaluation
//==============================================================
static double SafeSplineEval(gsl_spline *sp,
							 gsl_interp_accel *acc,
							 double x,
							 const char *what,
							 bool clamp = true)
{
	if (!sp || !sp->interp)
	{
		Z_LOG_ERROR(std::string("SafeSplineEval: null spline for ") + what);
		return 0.0;
	}

	const double xmin = sp->interp->xmin;
	const double xmax = sp->interp->xmax;

	double x_use = x;

	if (clamp)
	{
		if (x_use < xmin)
			x_use = xmin;
		if (x_use > xmax)
			x_use = xmax;
	}
	else
	{
		if (x_use < xmin || x_use > xmax)
		{
			Z_LOG_ERROR(std::string("SafeSplineEval: x out of range for ") + what +
						"  x=" + std::to_string(x_use) +
						"  [" + std::to_string(xmin) + "," + std::to_string(xmax) + "]");
			return 0.0;
		}
	}

	double y = 0.0;
	const int status = gsl_spline_eval_e(sp, x_use, acc, &y);

	if (status != GSL_SUCCESS)
	{
		Z_LOG_ERROR(std::string("SafeSplineEval: gsl_spline_eval_e failed for ") + what +
					" status=" + std::to_string(status));
		return 0.0;
	}

	if (clamp && x != x_use)
	{
		Z_LOG_WARNING(std::string("SafeSplineEval: clamped x for ") + what +
					  " from " + std::to_string(x) +
					  " to " + std::to_string(x_use));
	}

	return y;
}

//==============================================================
//        Safe Spline Derivative Evaluation — FAIL-CLOSED
//==============================================================
/**
 * A derivative is a *local* physical property of the state being asked about, so the clamping
 * that `SafeSplineEval` applies to a value lookup is not admissible here: clamping would return
 * the derivative at a boundary while pretending to describe the requested state. This helper
 * therefore refuses everything it cannot answer exactly — out of domain, non-finite argument,
 * missing interpolant, GSL error — and returns NaN so that the caller cannot mistake the result
 * for a physical number. Exact endpoints are inside the domain and are evaluated.
 *
 * ADR-0007 P5 (ACCEPTED 2026-09-02); Phase 4C-I0.
 */
static double SafeSplineEvalDeriv(gsl_spline *sp,
								  gsl_interp_accel *acc,
								  double x,
								  const char *what)
{
	const double kNaN = std::numeric_limits<double>::quiet_NaN();

	if (!sp || !sp->interp)
	{
		Z_LOG_ERROR(std::string("SafeSplineEvalDeriv: null spline for ") + what);
		return kNaN;
	}

	if (!std::isfinite(x))
	{
		Z_LOG_ERROR(std::string("SafeSplineEvalDeriv: non-finite x for ") + what);
		return kNaN;
	}

	const double xmin = sp->interp->xmin;
	const double xmax = sp->interp->xmax;

	// Deliberately NOT clamped — see the comment above.
	if (x < xmin || x > xmax)
	{
		Z_LOG_ERROR(std::string("SafeSplineEvalDeriv: x out of range for ") + what +
					"  x=" + std::to_string(x) +
					"  [" + std::to_string(xmin) + "," + std::to_string(xmax) + "]");
		return kNaN;
	}

	double dydx = 0.0;
	const int status = gsl_spline_eval_deriv_e(sp, x, acc, &dydx);

	if (status != GSL_SUCCESS)
	{
		Z_LOG_ERROR(std::string("SafeSplineEvalDeriv: gsl_spline_eval_deriv_e failed for ") +
					what + " status=" + std::to_string(status));
		return kNaN;
	}

	if (!std::isfinite(dydx))
	{
		Z_LOG_ERROR(std::string("SafeSplineEvalDeriv: non-finite derivative for ") + what);
		return kNaN;
	}

	return dydx;
}

//==============================================================
/**
 * Converts the raw `eps(p)` spline derivative from the EOS table's own units
 * (g/cm^3 per dyne/cm^2, i.e. s^2/cm^2) to the dimensionless geometric `d(eps)/dp`.
 *
 * Derived from the *same* two constants `NStar::BuildFromTOV` uses to place `eps` and `p` on
 * the profile:
 *
 *     eps_geom = eps_cgs * (INV_FM4_2_INV_KM2 / INV_FM4_2_G_CM3)
 *     p_geom   = p_cgs   * (INV_FM4_2_INV_KM2 / INV_FM4_2_Dyn_CM2)
 *  => d(eps_geom)/d(p_geom) = (INV_FM4_2_Dyn_CM2 / INV_FM4_2_G_CM3) * d(eps_cgs)/d(p_cgs)
 *
 * The ratio is `c^2` in cgs (cm^2/s^2), because `INV_FM4_2_Dyn_CM2` converts an energy density
 * to dyne/cm^2 = erg/cm^3 while `INV_FM4_2_G_CM3` converts the same energy density to the
 * mass-density equivalent g/cm^3. That identity is asserted against an independent literal
 * `c = 2.99792458e10 cm/s` by `tests/eos/eos_derivative_contract.cpp`; it is not assumed here.
 */
static double EosDerivCgsToDimensionless()
{
	return Zaki::Physics::INV_FM4_2_Dyn_CM2 / Zaki::Physics::INV_FM4_2_G_CM3;
}

//==============================================================
//                        Sequence class
//==============================================================
// Constructor
Sequence::Sequence() : Prog("Sequence")
{
}
//--------------------------------------------------------------
void Sequence::Add(const NStar &in_star)
{
	seq.emplace_back(in_star.GetSequence());
}

//--------------------------------------------------------------
// Exports the star sequence
void Sequence::Export(const Zaki::String::Directory &in_dir)
	const
{
	Zaki::File::VecSaver vec_saver(wrk_dir_ + in_dir);

	char seq_header[400];
	snprintf(seq_header, sizeof(seq_header), "%-14s\t %-14s\t %-14s\t %-14s\t %-14s"
											 "\t %-14s",
			 "ec(g/cm^3)", "M(Sun)", "R(km)", "pc(dyne/cm^2)", "B",
			 "I(km^3)");

	vec_saver.SetHeader(seq_header);
	// std::cout << " ---> in_dir = " << in_dir << "\n";
	// std::cout
	// 	<< "[debug] Exporting sequence to "
	// 	<< (wrk_dir_ + in_dir).Str() << "\n";
	vec_saver.Export1D(seq);
}

//--------------------------------------------------------------
// Combines two sequences
void Sequence::Combine(const Sequence &other)
{
	// std::lock_guard<std::mutex> lock(m_mutex) ;

	seq.reserve(seq.size() + other.seq.size());

	for (auto &&s : other.seq)
	{
		seq.emplace_back(s);
	}
	// std::cout << "[ Thread = " << std::this_thread::get_id()
	//           << " ] " << name << " Size = "
	//           << seq.size() << "\n" ;
}

//--------------------------------------------------------------
// Clears the sequence
void Sequence::Clear()
{
	seq.clear();
}
//--------------------------------------------------------------

//==============================================================
//                        MixedSequence class
//==============================================================
// Constructor
MixedSequence::MixedSequence() : Prog("MixedSequence")
{
}
//--------------------------------------------------------------

void MixedSequence::Add(const MixedStar &in_star)
{
	seq.emplace_back(in_star.sequence);
}

//--------------------------------------------------------------
// Exports the mixed star sequence
void MixedSequence::Export(const Zaki::String::Directory &in_dir)
	const
{
	Zaki::File::VecSaver vec_saver(wrk_dir_ + in_dir);

	char seq_header[400];
	snprintf(seq_header, sizeof(seq_header), "%-6s\t %-6s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s"
											 "\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s",
			 "v_idx", "d_idx", "ec(g/cm^3)", "M(Sun)", "R(km)", "pc(dyne/cm^2)", "B",
			 "I(km^3)", "ec_d(g/cm^3)", "M_d(Sun)", "R_d(km)",
			 "pc_d(dyne/cm^2)", "B_d", "I_d(km^3)");

	vec_saver.SetHeader(seq_header);

	vec_saver.Export1D(seq);
}
//--------------------------------------------------------------
// Combines two sequences
void MixedSequence::Combine(const MixedSequence &other)
{
	// std::lock_guard<std::mutex> lock(m_mutex) ;

	seq.reserve(seq.size() + other.seq.size());

	for (auto &&s : other.seq)
	{
		seq.emplace_back(s);
	}
	// std::cout << "[ Thread = " << std::this_thread::get_id()
	//           << " ] " << name << " Size = "
	//           << seq.size() << "\n" ;
}

//--------------------------------------------------------------
// Clears the sequence
void MixedSequence::Clear()
{
	seq.clear();
}
//--------------------------------------------------------------

//==============================================================

//==============================================================
//                        TOVSolver class
//==============================================================
// Constructor
TOVSolver::TOVSolver() : Prog("TOVSolver")
{
	// prevent GSL from aborting the entire process on domain errors
	gsl_set_error_handler_off();

	mixed_r_accel = gsl_interp_accel_alloc();
	visi_p_accel = gsl_interp_accel_alloc();
	dark_p_accel = gsl_interp_accel_alloc();

	// tov_counter++ ;
}

//--------------------------------------------------------------
TOVSolver::~TOVSolver()
{

	if (visi_eps_p_spline)
		gsl_spline_free(visi_eps_p_spline);

	if (dark_eps_p_spline)
		gsl_spline_free(dark_eps_p_spline);

	if (visi_rho_p_spline)
		gsl_spline_free(visi_rho_p_spline);

	if (dark_rho_p_spline)
		gsl_spline_free(dark_rho_p_spline);

	// if(rho_r_spline)
	//   gsl_spline_free (rho_r_spline);

	// if(rho_r_spline_dark)
	//   gsl_spline_free (rho_r_spline_dark);

	if (nu_der_r_spline)
		gsl_spline_free(nu_der_r_spline);

	if (mixed_r_accel)
		gsl_interp_accel_free(mixed_r_accel);

	// if(dark_accel)
	//   gsl_interp_accel_free (dark_accel);

	if (visi_p_accel)
		gsl_interp_accel_free(visi_p_accel);

	if (dark_p_accel)
		gsl_interp_accel_free(dark_p_accel);

	for (auto sp : visi_rho_i_p_spline)
	{
		if (sp)
			gsl_spline_free(sp);
	}
	for (auto sp : dark_rho_i_p_spline)
	{
		if (sp)
			gsl_spline_free(sp);
	}

	// tov_counter-- ;

	// {
	//   std::lock_guard<std::mutex> lock(m_mutex) ;

	//   mixed_seq_static.Combine(mixed_sequence) ;
	//   // std::cout << "\n\t tov_counter = " << tov_counter << "\n";
	// }

	// if( tov_counter == 0 )
	// {
	//   mixed_seq_static.Export("Full_Mixed_Sequence.txt") ;
	// }
}

//--------------------------------------------------------------
void TOVSolver::SetRadialRes(const size_t &in_radial_res)
{
	radial_res = in_radial_res;
}

//--------------------------------------------------------------
/// Sets maximum value for radius in the solver
/// if a maximum is known, say 15 km for NS,
/// it would increase the radial resolution
/// defined by:
///  delta R = scale(R) * (R_max - R_min) / radial_res
/// Unit must be in km
void TOVSolver::SetMaxRadius(const double &in_r_max)
{
	r_max = in_r_max * 1e5;
	char msg[100];
	snprintf(msg, 100, "Maximum radius changed to %.2f km.", in_r_max);
	Z_LOG_INFO(msg);
}

//--------------------------------------------------------------
/// Sets the work directory for the member objects
void TOVSolver::OnWorkDirChanged(const Zaki::String::Directory &in_dir)
{
	n_star.SetWrkDir(in_dir);
	sequence.SetWrkDir(in_dir);

	mixed_star.SetWrkDir(in_dir);
	mixed_sequence.SetWrkDir(in_dir);

	// {
	//   std::lock_guard<std::mutex> lock(m_mutex) ;

	//   if ( !mixed_seq_static.IsWrkDirSet() )
	//   {
	//     mixed_seq_static.SetWrkDir(in_dir) ;
	//   }
	// }

	// return this;
}

//--------------------------------------------------------------
// void TOVSolver::Hidden_ImportEOS_Vis(const Zaki::String::Directory &f_name)
// {
// 	std::ifstream file((wrk_dir_ + "/" + f_name).Str());

// 	// Error opening the file
// 	if (file.fail())
// 	{
// 		Z_LOG_ERROR("File '" + (wrk_dir_ + "/" + f_name).Str() + "' cannot be opened!");
// 		Z_LOG_ERROR("Importing EOS data failed!");
// 		exit(EXIT_FAILURE);
// 		return;
// 	}

// 	size_t line_num = 0;
// 	// Reading the input file
// 	for (Zaki::File::CSVIterator loop(file, '\t'); loop != Zaki::File::CSVIterator(); ++loop)
// 	{
// 		if ((*loop).size() < 3)
// 		{
// 			std::cout << "\n (*loop)[0]: " << (*loop)[0] << "\n";
// 			std::cout << "(*loop).size(): " << (*loop).size() << "\n";
// 			std::cout << "Line number: " << line_num << "\n";
// 			Z_LOG_ERROR("EOS file is not complete!");
// 			break;
// 		}

// 		// First line
// 		if (line_num == 0)
// 		{
// 			eos_tab.SetLabels(
// 				Zaki::String::Strip((*loop)[0], ' '),
// 				Zaki::String::Strip((*loop)[1], ' '),
// 				Zaki::String::Strip((*loop)[2], ' '));

// 			// We want to do this once only
// 			if ((*loop).size() > 3)
// 			{
// 				for (size_t i = 3; i < (*loop).size(); i++)
// 				{
// 					eos_tab.rho_i.push_back({});
// 					eos_tab.AddExtraLabels(Zaki::String::Strip((*loop)[i], ' '));
// 				}
// 			}
// 		}
// 		else
// 		{
// 			eos_tab.eps.push_back(std::atof((*loop)[0].c_str()));
// 			eos_tab.pre.push_back(std::atof((*loop)[1].c_str()));
// 			eos_tab.rho.push_back(std::atof((*loop)[2].c_str()));

// 			for (size_t i = 3; i < (*loop).size(); i++)
// 				eos_tab.rho_i[i - 3].push_back(
// 					std::atof((*loop)[i].c_str()));
// 		}
// 		line_num++;
// 	}

// 	Z_LOG_INFO("EOS data imported from: " + (wrk_dir_ + f_name).Str() + ".");

// 	visi_eps_p_spline = gsl_spline_alloc(TOV_gsl_interp_type, eos_tab.Size());
// 	visi_rho_p_spline = gsl_spline_alloc(TOV_gsl_interp_type, eos_tab.Size());

// 	for (size_t i = 0; i < eos_tab.pre.size() - 1; i++)
// 	{
// 		if (eos_tab.pre[i] > eos_tab.pre[i + 1])
// 			std::cout << "\n P[ " << i << "] = " << eos_tab.pre[i]
// 					  << " , P[" << i + 1 << "]" << eos_tab.pre[i + 1] << "\n";
// 	}

// 	for (size_t i = 0; i < eos_tab.rho_i.size(); i++)
// 	{
// 		visi_rho_i_p_spline.emplace_back(gsl_spline_alloc(TOV_gsl_interp_type, eos_tab.Size()));
// 		// std::cout << "\n\t eos_tab.pre[0]= " << eos_tab.pre[0] << "\n" << std::flush ;
// 		gsl_spline_init(visi_rho_i_p_spline[i], &eos_tab.pre[0], &eos_tab.rho_i[i][0], eos_tab.Size());
// 	}

// 	Z_LOG_INFO("Initializing the splines for energy density and pressure.");
// 	// This function initializes the interpolation object
// 	// x has to be strictly increasing
// 	gsl_spline_init(visi_eps_p_spline, &eos_tab.pre[0], &eos_tab.eps[0], eos_tab.Size());
// 	gsl_spline_init(visi_rho_p_spline, &eos_tab.pre[0], &eos_tab.rho[0], eos_tab.Size());
// }
//--------------------------------------------------------------
void TOVSolver::Hidden_ImportEOS_Vis(const Zaki::String::Directory &eos_file, const bool absolute_path)
{
	// 0) Resolve the file path EXACTLY as we will open it
	//    (this was ambiguous before)
	Zaki::String::Directory full_path("");
	if (absolute_path)
	{
		full_path = eos_file; // use as-is
	}
	else
	{
		full_path = wrk_dir_ + "/" + eos_file; // relative to work dir
	}

	// std::cout << "\n[TOVSolver::Hidden_ImportEOS_Vis] trying to open: '"
	// 		  << full_path.Str() << "'\n";

	std::ifstream file(full_path.Str());

	// 1) Error opening the file
	if (file.fail())
	{
		Z_LOG_ERROR("File '" + full_path.Str() + "' cannot be opened!");
		Z_LOG_ERROR("Importing EOS data failed!");
		// keep exit for now so we see it clearly
		exit(EXIT_FAILURE);
		return;
	}

	// ---- Free any existing visible splines before re-import ----
	if (visi_eps_p_spline)
	{
		gsl_spline_free(visi_eps_p_spline);
		visi_eps_p_spline = nullptr;
	}
	if (visi_rho_p_spline)
	{
		gsl_spline_free(visi_rho_p_spline);
		visi_rho_p_spline = nullptr;
	}

	for (auto &sp : visi_rho_i_p_spline)
	{
		if (sp)
			gsl_spline_free(sp);
	}
	visi_rho_i_p_spline.clear();

	// clear previous data if any
	eos_tab.eps.clear();
	eos_tab.pre.clear();
	eos_tab.rho.clear();
	eos_tab.rho_i.clear();
	eos_tab.extra_labels.clear();

	size_t line_num = 0;

	// 2) Reading the input file
	for (Zaki::File::CSVIterator loop(file, '\t');
		 loop != Zaki::File::CSVIterator(); ++loop)
	{
		// print first few raw lines to see what’s actually in the file
		if (line_num < 5)
		{
			std::cout << "[EOS:raw] line " << line_num << " has "
					  << (*loop).size() << " fields.\n";
			for (size_t k = 0; k < (*loop).size(); ++k)
			{
				std::cout << "    [" << k << "] = '" << (*loop)[k] << "'\n";
			}
		}

		// if ((*loop).size() < 3)
		// {
		// 	std::cout << "\n[EOS] (*loop)[0]: " << (*loop)[0] << "\n";
		// 	std::cout << "[EOS] (*loop).size(): " << (*loop).size() << "\n";
		// 	std::cout << "[EOS] Line number: " << line_num << "\n";
		// 	Z_LOG_ERROR("EOS file is not complete!");
		// 	break;
		// }
		if ((*loop).size() < 3)
		{
			std::cout << "[EOS] Line " << line_num << " has only "
					  << (*loop).size() << " fields.\n";
			if ((*loop).size() != 0)
				std::cout << "  first token: '" << (*loop)[0] << "'\n";
			Z_LOG_ERROR("EOS file is not complete!");
			continue;
			;
		}

		if (line_num == 0)
		{
			// header
			eos_tab.SetLabels(
				Zaki::String::Strip((*loop)[0], ' '),
				Zaki::String::Strip((*loop)[1], ' '),
				Zaki::String::Strip((*loop)[2], ' '));

			// extra species
			if ((*loop).size() > 3)
			{
				for (size_t i = 3; i < (*loop).size(); i++)
				{
					eos_tab.rho_i.push_back({});
					eos_tab.AddExtraLabels(Zaki::String::Strip((*loop)[i], ' '));
				}
			}

			std::cout << "[EOS] headers: eps='" << (*loop)[0]
					  << "', p='" << (*loop)[1]
					  << "', rho='" << (*loop)[2] << "'\n";
			if (!eos_tab.extra_labels.empty())
			{
				std::cout << "[EOS] extra columns: ";
				for (auto &lbl : eos_tab.extra_labels)
					std::cout << lbl << " ";
				std::cout << "\n";
			}
		}
		else
		{
			const double eps_val = std::atof((*loop)[0].c_str());
			const double pre_val = std::atof((*loop)[1].c_str());
			const double rho_val = std::atof((*loop)[2].c_str());

			eos_tab.eps.push_back(eps_val);
			eos_tab.pre.push_back(pre_val);
			eos_tab.rho.push_back(rho_val);

			// fill extra columns
			for (size_t i = 3; i < (*loop).size(); i++)
			{
				eos_tab.rho_i[i - 3].push_back(std::atof((*loop)[i].c_str()));
			}
		}
		line_num++;
	}

	// 3) After the loop: print sizes
	const std::size_t n = eos_tab.Size();
	std::cout << "[EOS] imported rows (excluding header): " << n << "\n";
	std::cout << "[EOS] eps.size() = " << eos_tab.eps.size() << "\n";
	std::cout << "[EOS] pre.size() = " << eos_tab.pre.size() << "\n";
	std::cout << "[EOS] rho.size() = " << eos_tab.rho.size() << "\n";
	for (size_t i = 0; i < eos_tab.rho_i.size(); ++i)
	{
		std::cout << "[EOS] rho_i[" << i << "].size() = "
				  << eos_tab.rho_i[i].size() << "\n";
	}

	Z_LOG_INFO("EOS data imported from: " + full_path.Str() + ".");

	// 4) If we have < 2 data points, DO NOT build splines
	if (n < 2)
	{
		Z_LOG_ERROR("EOS has too few data points (" + std::to_string(n) +
					") to build GSL splines. Check the path or file format.");
		return; // <- IMPORTANT: don't try to build splines
	}

	// 5) Check monotonicity of pressure – this was already here, just make it louder
	for (size_t i = 0; i + 1 < eos_tab.pre.size(); i++)
	{
		if (eos_tab.pre[i] >= eos_tab.pre[i + 1])
		{
			Z_LOG_ERROR("[EOS][WARN] P[" + std::to_string(i) + "] = " +
						std::to_string(eos_tab.pre[i]) + "  >=  P[" +
						std::to_string(i + 1) + "] = " +
						std::to_string(eos_tab.pre[i + 1]) +
						"  (pressure must be strictly increasing for GSL)\n");
			return;
		}
	}

	// 6) alloc + init (visible)
	visi_eps_p_spline = gsl_spline_alloc(TOV_gsl_interp_type, n);
	visi_rho_p_spline = gsl_spline_alloc(TOV_gsl_interp_type, n);

	Z_LOG_INFO("Initializing main splines for EOS with n = " + std::to_string(n));
	gsl_spline_init(visi_eps_p_spline,
					eos_tab.pre.data(),
					eos_tab.eps.data(),
					n);
	gsl_spline_init(visi_rho_p_spline,
					eos_tab.pre.data(),
					eos_tab.rho.data(),
					n);

	// 7) extra species
	visi_rho_i_p_spline.clear();
	for (size_t i = 0; i < eos_tab.rho_i.size(); i++)
	{
		if (eos_tab.rho_i[i].size() != n)
		{
			Z_LOG_WARNING("Extra column " + std::to_string(i) +
						  " has size " + std::to_string(eos_tab.rho_i[i].size()) +
						  " but expected " + std::to_string(n) +
						  " – skipping spline.");
			visi_rho_i_p_spline.emplace_back(nullptr);
			continue;
		}

		gsl_spline *sp = gsl_spline_alloc(TOV_gsl_interp_type, n);
		gsl_spline_init(sp,
						eos_tab.pre.data(),
						eos_tab.rho_i[i].data(),
						n);
		visi_rho_i_p_spline.emplace_back(sp);
	}

	Z_LOG_INFO("Initializing the splines for energy density and pressure: done.");
}

//--------------------------------------------------------------
void TOVSolver::ImportEOS(const Zaki::String::Directory &f_name,
						  const bool absolute_path)
{
	Hidden_ImportEOS_Vis(f_name, absolute_path);

	n_star.InitFromTOVSolver(this);
}

//--------------------------------------------------------------
void TOVSolver::ImportEOS(const Zaki::String::Directory &vis_eos,
						  const Zaki::String::Directory &dar_eos,
						  const bool absolute_path)
{
	Hidden_ImportEOS_Vis(vis_eos, absolute_path);
	Hidden_ImportEOS_Dar(dar_eos, absolute_path);

	mixed_star.InitVisible(this);
	mixed_star.InitDark(this);
}

//--------------------------------------------------------------
void TOVSolver::Hidden_ImportEOS_Dar(const Zaki::String::Directory &eos_file,
									 const bool absolute_path)
{
	Zaki::String::Directory full_path("");
	if (absolute_path)
	{
		full_path = eos_file; // use as-is
	}
	else
	{
		full_path = wrk_dir_ + "/" + eos_file; // relative to work dir
	}

	std::ifstream file(full_path.Str());

	// Error opening the file
	if (file.fail())
	{
		Z_LOG_ERROR("File '" + full_path.Str() + "' cannot be opened!");
		Z_LOG_ERROR("Importing EOS data failed!");
		exit(EXIT_FAILURE);
		return;
	}

	// ---- Free any existing dark splines before re-import ----
	if (dark_eps_p_spline)
	{
		gsl_spline_free(dark_eps_p_spline);
		dark_eps_p_spline = nullptr;
	}
	if (dark_rho_p_spline)
	{
		gsl_spline_free(dark_rho_p_spline);
		dark_rho_p_spline = nullptr;
	}

	for (auto &sp : dark_rho_i_p_spline)
	{
		if (sp)
			gsl_spline_free(sp);
	}
	dark_rho_i_p_spline.clear();
	// ----------------------------------------
	eos_tab_dark.eps.clear();
	eos_tab_dark.pre.clear();
	eos_tab_dark.rho.clear();
	eos_tab_dark.rho_i.clear();
	eos_tab_dark.extra_labels.clear();
	// ----------------------------------------
	size_t line_num = 0;
	// Reading the input file
	;
	for (Zaki::File::CSVIterator loop(file, '\t'); loop != Zaki::File::CSVIterator(); ++loop)
	{
		if ((*loop).size() < 3)
		{
			std::cout << "\n (*loop)[0]: " << (*loop)[0] << "\n";
			std::cout << "(*loop).size(): " << (*loop).size() << "\n";
			Z_LOG_ERROR("EOS file is not complete!");
			continue;
		}

		// First line
		if (line_num == 0)
		{
			// eos_tab_dark.SetLabels((*loop)[0], (*loop)[1], (*loop)[2]);

			eos_tab_dark.SetLabels(
				Zaki::String::Strip((*loop)[0], ' '),
				Zaki::String::Strip((*loop)[1], ' '),
				Zaki::String::Strip((*loop)[2], ' '));
			// We want to do this once only
			if ((*loop).size() > 3)
			{
				for (size_t i = 3; i < (*loop).size(); i++)
				{
					eos_tab_dark.rho_i.push_back({});
					eos_tab_dark.AddExtraLabels(Zaki::String::Strip((*loop)[i], ' '));
				}
			}
		}
		else
		{
			eos_tab_dark.eps.push_back(std::atof((*loop)[0].c_str()));
			eos_tab_dark.pre.push_back(std::atof((*loop)[1].c_str()));
			eos_tab_dark.rho.push_back(std::atof((*loop)[2].c_str()));

			for (size_t i = 3; i < (*loop).size(); i++)
				eos_tab_dark.rho_i[i - 3].push_back(
					std::atof((*loop)[i].c_str()));
		}
		line_num++;
	}

	Z_LOG_INFO("Dark EOS data imported from: " + full_path.Str() + ".");

	const std::size_t n = eos_tab_dark.Size();
	if (n < 2)
	{
		Z_LOG_ERROR("Dark EOS has too few data points (" + std::to_string(n) + ")");
		return;
	}

	// ----------------------------------------
	// Check monotonicity of pressure
	for (size_t i = 0; i + 1 < eos_tab_dark.pre.size(); ++i)
	{
		if (eos_tab_dark.pre[i] >= eos_tab_dark.pre[i + 1])
		{
			Z_LOG_ERROR("P[" + std::to_string(i) +
						"] = " + std::to_string(eos_tab_dark.pre[i]) +
						" >= P[" + std::to_string(i + 1) + "] = " +
						std::to_string(eos_tab_dark.pre[i + 1]) +
						" (pressure must be strictly increasing for GSL)\n");
			return;
		}
	}
	// ----------------------------------------
	dark_eps_p_spline = gsl_spline_alloc(TOV_gsl_interp_type, n);
	dark_rho_p_spline = gsl_spline_alloc(TOV_gsl_interp_type, n);

	// for (size_t i = 0; i < eos_tab_dark.rho_i.size(); i++)
	// {
	// 	dark_rho_i_p_spline.emplace_back(gsl_spline_alloc(TOV_gsl_interp_type, n));
	// 	gsl_spline_init(dark_rho_i_p_spline[i], eos_tab_dark.pre.data(), eos_tab_dark.rho_i[i].data(), n);
	// }
	dark_rho_i_p_spline.clear();
	for (size_t i = 0; i < eos_tab_dark.rho_i.size(); ++i)
	{
		if (eos_tab_dark.rho_i[i].size() != n)
		{
			Z_LOG_WARNING("[EOS_DARK][WARN] extra column " + std::to_string(i) +
						  " has size " + std::to_string(eos_tab_dark.rho_i[i].size()) +
						  " but expected " + std::to_string(n) + " – skipping spline.");
			dark_rho_i_p_spline.emplace_back(nullptr);
			continue;
		}

		gsl_spline *sp = gsl_spline_alloc(TOV_gsl_interp_type, n);
		gsl_spline_init(sp, eos_tab_dark.pre.data(), eos_tab_dark.rho_i[i].data(), n);
		dark_rho_i_p_spline.emplace_back(sp);
	}

	// This function initializes the interpolation object
	// x has to be strictly increasing
	gsl_spline_init(dark_eps_p_spline, eos_tab_dark.pre.data(), eos_tab_dark.eps.data(), n);
	gsl_spline_init(dark_rho_p_spline, eos_tab_dark.pre.data(), eos_tab_dark.rho.data(), n);
}

//--------------------------------------------------------------
// // Full table (whatever EOSTable::Print() does by default)
// void TOVSolver::PrintEOSTable() const
// {
// 	eos_tab.Print(); // prints header + all
// }

//--------------------------------------------------------------
// Bounded table: header + first `max_rows` data rows
void TOVSolver::PrintEOSTable(const std::size_t max_rows) const
{
	eos_tab.Print(max_rows);
}

//--------------------------------------------------------------
// Compact summary: counts, labels, min/max
void TOVSolver::PrintEOSSummary() const
{
	eos_tab.PrintSummary();
}

//--------------------------------------------------------------
double TOVSolver::GetEOSMinEDens() const
{
	if (eos_tab.eps.empty())
		return 0.0;
	return eos_tab.eps.front();
}

//--------------------------------------------------------------
double TOVSolver::GetEOSMaxEDens() const
{
	if (eos_tab.eps.empty())
		return 0.0;
	return eos_tab.eps.back();
}

//--------------------------------------------------------------
void TOVSolver::SetCentralEDensFloorFactor(double f)
{
	central_eps_floor_factor = f;
}

//--------------------------------------------------------------
//            Added on December 15, 2020
/// Returns the nu_der value given the radius input
// r is in cm!
double TOVSolver::GetNuDerSpline(const double &in_r)
{
	return SafeSplineEval(nu_der_r_spline, mixed_r_accel, in_r, "nu'(r) spline");
}

//--------------------------------------------------------------
//            Added on December 15, 2020
/// Evaluates & exports the nu(r) function table
/// input must be a TOV solution file (with nu' column)
// void TOVSolver::ExportNu(const Zaki::String::Directory& f_name)
// {
//   std::ifstream     file( (wrk_dir_ + "/" + f_name).Str());

//   // Error opening the file
//   if (file.fail())
//   {
//     Z_LOG_ERROR("File '"+(wrk_dir_ + "/" + f_name).Str() +"' cannot be opened!") ;
//     Z_LOG_ERROR("Importing TOV solution failed!") ;
//     exit(EXIT_FAILURE) ;
//     return ;
//   }

//   std::vector<double> tov_radius ;
//   std::vector<double> tov_nu_der ;
//   std::vector<double> tov_mass   ;

//   size_t line_num = 0 ;
//   for(Zaki::File::CSVIterator loop(file, '\t'); loop != Zaki::File::CSVIterator(); ++loop)
//   {
//     // First line
//     if (line_num == 0)
//     {
//       if( Zaki::String::Strip((*loop)[2], ' ') != "nu'")
//       {
//         Z_LOG_ERROR("TOV solution doesn't include nu'(r) !") ;
//         return ;
//       }
//     }
//     // Other lines:
//     else
//     {
//       // Converting the radius into cm --> Why?
//       // Mar-2022: Because nu' is in (1/cm).
//       tov_radius.push_back(std::atof((*loop)[0].c_str()) * 1e5) ;
//       tov_mass.push_back(std::atof((*loop)[1].c_str())) ;
//       tov_nu_der.push_back(std::atof((*loop)[2].c_str())) ;
//     }
//     line_num++ ;
//   }

//   Z_LOG_INFO("TOV soluton imported from: "+ (wrk_dir_+f_name).Str()+".") ;

//   nu_der_r_spline  = gsl_spline_alloc (TOV_gsl_interp_type, tov_radius.size());

//   // This function initializes the interpolation object
//   // x has to be strictly increasing
//   gsl_spline_init (nu_der_r_spline, &tov_radius[0], &tov_nu_der[0], tov_radius.size());

//   // -----------------------------------
//   // Integrate to find nu(r) :
//   // -----------------------------------

//   gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000*tov_radius.size());
//   double err;

//   Zaki::Math::GSLFuncWrapper<TOVSolver, double (TOVSolver::*)( const double& )>
//     Fp(this, &TOVSolver::GetNuDerSpline);

//   gsl_function F = static_cast<gsl_function> (Fp) ;

//   double tmp_nu = 0;
//   std::vector<double> nu_result ;
//   nu_result.reserve(tov_radius.size()) ;

//   for(size_t i = 0 ; i < tov_radius.size() ; ++i)
//   {
//     gsl_integration_qag(&F, tov_radius[0], tov_radius[i], 1e-9, 1e-9, 1000, 1, w, &tmp_nu, &err);
//     nu_result.push_back(tmp_nu) ;
//     tmp_nu = 0 ;
//   }

//   gsl_integration_workspace_free(w);

//   // -------------------------------------------
//   // Saving the results for [ r, M(r), nu(r) ]
//   // -------------------------------------------

//   // Matching the boundary condition for nu(r)
//   double nu_at_R = 0.5*log(1
//                   - 2*Zaki::Physics::SUN_M_KM*tov_mass[tov_mass.size()-1]
//                     / (tov_radius[tov_mass.size()-1]*1e-5)
//                   ) ;

//   double delta_nu_r = nu_result[tov_mass.size()-1] - nu_at_R ;
//   std::vector<TOV_Nu_Point> rmnu_results ;
//   rmnu_results.reserve(tov_radius.size()) ;

//   for (size_t i = 0; i < tov_radius.size(); i++)
//   {
//     rmnu_results.emplace_back(tov_radius[i]*1e-5,
//                     tov_mass[i],
//                     nu_result[i] - delta_nu_r) ;
//   }

//   std::string out_f_name = Zaki::String::Pars((wrk_dir_ + "/" + f_name).Str(), ".tsv")[0] ;

//   Zaki::File::VecSaver vec_saver(out_f_name + "_nu.tsv");

//   char seq_header[200] ;
//   snprintf(seq_header,  sizeof(seq_header), "%-14s\t %-14s\t %-14s",
//           "r(km)", "m",  "nu") ;
//   vec_saver.SetHeader(seq_header) ;
//   vec_saver.Export1D(rmnu_results) ;

// }

//--------------------------------------------------------------
// Input pressure, output energy density
double TOVSolver::GetEDens(const double &in_pres)
{
	return SafeSplineEval(visi_eps_p_spline, visi_p_accel, in_pres, "visible eps(p) spline");
}

//--------------------------------------------------------------
// Input pressure, output energy density (dark sector)
double TOVSolver::GetEDens_Dark(const double &in_pres)
{
	return SafeSplineEval(dark_eps_p_spline, dark_p_accel, in_pres, "dark eps(p) spline");
}

//--------------------------------------------------------------
// EOS thermodynamic derivative — ADR-0007 P5 (Phase 4C-I0).
// The authority is the SAME visi_eps_p_spline that GetEDens reads, through the SAME
// accelerator; there is no second interpolant anywhere in this class.
bool TOVSolver::HasEDensDeriv() const noexcept
{
	return visi_eps_p_spline != nullptr && visi_eps_p_spline->interp != nullptr;
}

//--------------------------------------------------------------
double TOVSolver::EDensDerivPressMin() const noexcept
{
	if (!HasEDensDeriv())
		return std::numeric_limits<double>::quiet_NaN();
	return visi_eps_p_spline->interp->xmin;
}

//--------------------------------------------------------------
double TOVSolver::EDensDerivPressMax() const noexcept
{
	if (!HasEDensDeriv())
		return std::numeric_limits<double>::quiet_NaN();
	return visi_eps_p_spline->interp->xmax;
}

//--------------------------------------------------------------
// Input pressure (dyne/cm^2), output d(eps)/dp — DIMENSIONLESS.
// Fail-closed outside the interpolant's domain: see SafeSplineEvalDeriv.
double TOVSolver::GetEDensDeriv(const double &in_pres)
{
	const double dedp_cgs = SafeSplineEvalDeriv(visi_eps_p_spline, visi_p_accel, in_pres,
											   "visible eps(p) spline derivative");
	if (!std::isfinite(dedp_cgs))
		return std::numeric_limits<double>::quiet_NaN();

	return dedp_cgs * EosDerivCgsToDimensionless();
}

//--------------------------------------------------------------
double TOVSolver::cost_p_of_e(const double in_p)
{
	return GetEDens(in_p) - cost_p_of_e_input;
}

//--------------------------------------------------------------
double TOVSolver::cost_p_of_e_dark(const double in_p)
{
	return GetEDens_Dark(in_p) - cost_p_of_e_input_dark;
}

//--------------------------------------------------------------
// // Inverse function of "GetEDens"
// double TOVSolver::p_of_e(const double &in_e)
// {
// 	double p_min = eos_tab.pre[0];
// 	double p_max = eos_tab.pre[eos_tab.pre.size() - 1];

// 	// std::cout << "p_min = " << p_min
// 	//           << ", p_max = " << p_max ;

// 	// std::cout << ", p_of_e(" <<in_e << ") = " ;

// 	cost_p_of_e_input = in_e;

// 	Zaki::Math::GSLFuncWrapper<TOVSolver, double (TOVSolver::*)(double)>
// 		func(this, &TOVSolver::cost_p_of_e);

// 	gsl_function F = static_cast<gsl_function>(func);

// 	const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
// 	gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

// 	//  gsl_set_error_handler_off() ;
// 	gsl_root_fsolver_set(s, &F, p_min, p_max);

// 	int iter = 0, max_iter = 1000;
// 	int status;
// 	double out_p = 0;
// 	do
// 	{
// 		iter++;
// 		status = gsl_root_fsolver_iterate(s);
// 		out_p = gsl_root_fsolver_root(s);
// 		p_min = gsl_root_fsolver_x_lower(s);
// 		p_max = gsl_root_fsolver_x_upper(s);
// 		// status = gsl_root_test_interval (p_min, p_max,
// 		// 0.0001, 0.0001);
// 		// Edited on Apr 25
// 		status = gsl_root_test_interval(p_min, p_max,
// 										p_of_e_prec, p_of_e_prec);
// 	} while (status == GSL_CONTINUE && iter < max_iter);

// 	gsl_root_fsolver_free(s);

// 	// std::cout << out_p << ", GetEDens("
// 	//           << out_p << ") = "
// 	//           << GetEDens(out_p) << ", iter = " << iter
// 	//           << ", p_min = " << p_min
// 	//           << ", p_max = " << p_max << "\n" ;

// 	return out_p;
// }

//--------------------------------------------------------------
// Inverse function of "GetEDens" with EOS-range guards
double TOVSolver::p_of_e(const double &in_e)
{
	// EOS must be present
	if (eos_tab.eps.empty())
	{
		Z_LOG_ERROR("p_of_e(...) called but EOS table is empty.");
		return 0.0;
	}

	const double eos_e_min = eos_tab.eps.front();
	const double eos_e_max = eos_tab.eps.back();
	const double p_min = eos_tab.pre.front();
	const double p_max = eos_tab.pre.back();

	// ----------------------------------------------------------
	// Outside EOS (below): just return lowest pressure
	// ----------------------------------------------------------
	if (in_e <= eos_e_min)
	{
		Z_LOG_WARNING("p_of_e: requested e = " + std::to_string(in_e) +
					  " is below EOS min = " + std::to_string(eos_e_min) +
					  " -> clamping to p_min.");
		return p_min;
	}

	// ----------------------------------------------------------
	// Outside EOS (above): just return highest pressure
	// ----------------------------------------------------------
	if (in_e >= eos_e_max)
	{
		Z_LOG_WARNING("p_of_e: requested e = " + std::to_string(in_e) +
					  " is above EOS max = " + std::to_string(eos_e_max) +
					  " -> clamping to p_max.");
		return p_max;
	}

	// ----------------------------------------------------------
	// Inside EOS range: do the Brent inversion
	// ----------------------------------------------------------
	cost_p_of_e_input = in_e;

	Zaki::Math::GSLFuncWrapper<TOVSolver, double (TOVSolver::*)(double)>
		func(this, &TOVSolver::cost_p_of_e);

	gsl_function F = static_cast<gsl_function>(func);

	const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

	double a = p_min;
	double b = p_max;

	gsl_root_fsolver_set(s, &F, a, b);

	int iter = 0;
	int status = GSL_CONTINUE;
	double out_p = 0.0;
	const int max_iter = 1000;

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		out_p = gsl_root_fsolver_root(s);
		a = gsl_root_fsolver_x_lower(s);
		b = gsl_root_fsolver_x_upper(s);

		status = gsl_root_test_interval(a, b, p_of_e_prec, p_of_e_prec);
	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);

	return out_p;
}

//--------------------------------------------------------------
// Inverse function of "GetEDens_Dark"
double TOVSolver::p_of_e_dark(const double &in_e)
{
	// EOS must be present
	if (eos_tab_dark.eps.empty())
	{
		Z_LOG_ERROR("p_of_e_dark(...) called but EOS table is empty.");
		return 0.0;
	}

	// double p_min = eos_tab_dark.pre[0];
	// double p_max = eos_tab_dark.pre[eos_tab_dark.pre.size() - 1];

	const double eos_e_min = eos_tab_dark.eps.front();
	const double eos_e_max = eos_tab_dark.eps.back();
	const double p_min = eos_tab_dark.pre.front();
	const double p_max = eos_tab_dark.pre.back();

	// ----------------------------------------------------------
	// Outside EOS (below): just return lowest pressure
	// ----------------------------------------------------------
	if (in_e <= eos_e_min)
	{
		Z_LOG_WARNING("p_of_e_dark: requested e = " + std::to_string(in_e) +
					  " is below EOS min = " + std::to_string(eos_e_min) +
					  " -> clamping to p_min.");
		return p_min;
	}
	// ----------------------------------------------------------
	// Outside EOS (above): just return highest pressure
	// ----------------------------------------------------------
	if (in_e >= eos_e_max)
	{
		Z_LOG_WARNING("p_of_e_dark: requested e = " + std::to_string(in_e) +
					  " is above EOS max = " + std::to_string(eos_e_max) +
					  " -> clamping to p_max.");
		return p_max;
	}
	// ----------------------------------------------------------
	// Inside EOS range: do the Brent inversion
	// ----------------------------------------------------------
	cost_p_of_e_input_dark = in_e;

	Zaki::Math::GSLFuncWrapper<TOVSolver, double (TOVSolver::*)(double)>
		func(this, &TOVSolver::cost_p_of_e_dark);

	gsl_function F = static_cast<gsl_function>(func);

	const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
	gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);

	double a = p_min;
	double b = p_max;

	// gsl_set_error_handler_off(); // Redundant: added to the constructor
	gsl_root_fsolver_set(s, &F, a, b);

	int iter = 0, max_iter = 250;
	int status;
	double out_p = 0;
	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		out_p = gsl_root_fsolver_root(s);
		a = gsl_root_fsolver_x_lower(s);
		b = gsl_root_fsolver_x_upper(s);
		// status = gsl_root_test_interval (p_min, p_max,
		//                                   0.0001, 0.0001);
		// Edited on Apr 25
		status = gsl_root_test_interval(a, b,
										p_of_e_prec, p_of_e_prec);
	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);

	return out_p;
}

//--------------------------------------------------------------
double TOVSolver::PressureCutoff() const
{
	return std::max(1.e-15 * GetInitPress(), eos_tab.pre[0]);
	// return eos_tab.pre[0];
}

//--------------------------------------------------------------
double TOVSolver::PressureCutoff_Dark() const
{
	return std::max(1.e-15 * GetInitPress_Dark(), eos_tab_dark.pre[0]);
	// return eos_tab_dark.pre[0];
}

//--------------------------------------------------------------
// Modified ODE in the presence of a dark core
// ......................
// Dictionary :
// y[0] = visible mass(r)
// y[1] = dark mass(r)
// y[2] = visible pressure(r)
// y[3] = dark pressure(r)
// f[0] = visible m'(r)
// f[1] = dark m'(r)
// f[2] = visible p'(r)
// f[3] = dark p'(r)
// ......................
int TOVSolver::ODE_Dark_Core(double r, const double y[], double f[], void *params)
{
	TOVSolver *tov_obj = (TOVSolver *)params;

	// Visible surface reached
	if (y[2] <= tov_obj->PressureCutoff())
	{
#if TOV_SOLVER_VERBOSE
		printf("\u2554----------------- Visible Core Surface reached ----------------\u2557\n");
		printf("\u2551 %14s %14s %15s %20s\n", "R (km)", "M (M_sun)", "\u03B5_c (g/cm^3)", "\u2551");
		printf("\u2551 %14le %14le %14le %20s\n", r / 1.e+5,
			   y[0] / GSL_CONST_CGSM_SOLAR_MASS, tov_obj->GetEDens(tov_obj->GetInitPress()), "\u2551");
		printf("\u255A---------------------------------------------------------------\u255D\n");

		printf("\n --------------------------------------------------\n");
		std::cout << " Vis. P_0 =" << tov_obj->GetInitPress() << "\n";
		std::cout << " Vis. P Cut-Off =" << tov_obj->PressureCutoff() << "\n";
		std::cout << " Current Vis. P = y[2] = " << y[2] << "\n";
		printf(" --------------------------------------------------\n\n");
#endif
		tov_obj->dark_core = false;

		return GSL_EBADFUNC;
	}

	// Dark core surface reached
	if (y[3] <= tov_obj->PressureCutoff_Dark())
	{
#if TOV_SOLVER_VERBOSE
		printf("\u2554----------------- Dark Core Surface reached ----------------\u2557\n");
		printf("\u2551 %14s %14s %15s %17s\n", "R (km)", "M (M_sun)", "\u03B5_c (g/cm^3)", "\u2551");
		printf("\u2551 %14le %14le %14le %17s\n", r / 1.e+5,
			   y[1] / GSL_CONST_CGSM_SOLAR_MASS, tov_obj->GetEDens_Dark(tov_obj->GetInitPress_Dark()), "\u2551");
		printf("\u255A------------------------------------------------------------\u255D\n");

		printf("\n --------------------------------------------------\n");
		std::cout << " Dark P_0 =" << tov_obj->GetInitPress_Dark() << "\n";
		std::cout << " Dark P Cut-Off =" << tov_obj->PressureCutoff_Dark() << "\n";
		std::cout << " Current Dark P = y[3] = " << y[3] << "\n";
		printf(" --------------------------------------------------\n\n");
#endif
		tov_obj->dark_core = true;

		return GSL_EBADFUNC;
	}

	// Mass Continuity Equations
	f[0] = 4. * M_PI * r * r * tov_obj->GetEDens(y[2]);
	f[1] = 4. * M_PI * r * r * tov_obj->GetEDens_Dark(y[3]);

	// TOV Equation
	// visible pressure derivative
	f[2] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow(r, 2.)) * (tov_obj->GetEDens(y[2]) + (y[2] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))) * (y[0] + y[1] + 4 * M_PI * pow(r, 3.) * (y[2] + y[3]) / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)) / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * (y[0] + y[1]) / (pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));

	// dark pressure derivative
	f[3] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow(r, 2.)) * (tov_obj->GetEDens_Dark(y[3]) + (y[3] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))) * (y[0] + y[1] + 4 * M_PI * pow(r, 3.) * (y[2] + y[3]) / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)) / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * (y[0] + y[1]) / (pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));

	return GSL_SUCCESS;
}

//--------------------------------------------------------------
// Modified ODE in the presence of a dark core
// ......................
// Dictionary :
// y[0] = mantle mass(r)
// y[1] = mantle pressure(r)
// f[0] = mantle m'(r)
// f[1] = mantle p'(r)
// ......................
int TOVSolver::ODE_Dark_Mantle(double r, const double y[], double f[], void *params)
{
	TOVSolver *tov_obj = (TOVSolver *)params;
	double m_c = tov_obj->m_core;

	// Surface reached
	if (tov_obj->dark_core && y[1] <= tov_obj->PressureCutoff())
	{
#if TOV_SOLVER_VERBOSE
		printf("\u2554----------------- Vis. Mantle Surface reached ----------------\u2557\n");
		printf("\u2551 %14s %14s %15s %19s\n", "R (km)", "M (M_sun)", "\u03B5_c (g/cm^3)", "\u2551");
		printf("\u2551 %14le %14le %14le %19s\n", r / 1.e+5,
			   y[0] / GSL_CONST_CGSM_SOLAR_MASS, tov_obj->GetEDens(tov_obj->GetInitPress()), "\u2551");
		printf("\u255A--------------------------------------------------------------\u255D\n");

		printf("\n --------------------------------------------------\n");
		std::cout << " Vis. P_0 =" << tov_obj->GetInitPress() << "\n";
		std::cout << " Vis. P Cut-Off =" << tov_obj->PressureCutoff() << "\n";
		std::cout << " Current Vis. P = y[1] = " << y[1] << "\n";
		printf(" --------------------------------------------------\n\n");
#endif
		return GSL_EBADFUNC;
	}
	if (!(tov_obj->dark_core) && y[1] <= tov_obj->PressureCutoff_Dark())
	{
#if TOV_SOLVER_VERBOSE
		printf("\u2554----------------- Dark Mantle Surface reached ----------------\u2557\n");
		printf("\u2551 %14s %14s %15s %19s\n", "R (km)", "M (M_sun)", "\u03B5_c (g/cm^3)", "\u2551");
		printf("\u2551 %14le %14le %14le %19s\n", r / 1.e+5,
			   y[0] / GSL_CONST_CGSM_SOLAR_MASS, tov_obj->GetEDens_Dark(tov_obj->GetInitPress_Dark()), "\u2551");
		printf("\u255A--------------------------------------------------------------\u255D\n");

		printf("\n --------------------------------------------------\n");
		std::cout << " Dark P_0 =" << tov_obj->GetInitPress_Dark() << "\n";
		std::cout << " Dark P Cut-Off =" << tov_obj->PressureCutoff_Dark() << "\n";
		std::cout << " Current Dark P = y[1] = " << y[1] << "\n";
		printf(" --------------------------------------------------\n\n");
#endif
		return GSL_EBADFUNC;
	}

	double e_den = 0;

	if (tov_obj->dark_core)
		e_den = tov_obj->GetEDens(y[1]);
	else
		e_den = tov_obj->GetEDens_Dark(y[1]);

	// Mass Continuity Equations
	f[0] = 4. * M_PI * r * r * e_den;

	// TOV Equation
	// visible pressure derivative
	f[1] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow(r, 2.)) * (e_den + (y[1] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))) * (y[0] + m_c + 4 * M_PI * pow(r, 3.) * (y[1] + 0) / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)) / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * (y[0] + m_c) / (pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));

	return GSL_SUCCESS;
}

//--------------------------------------------------------------
// ......................
// Dictionary :
// y[0] = mass(r)
// y[1] = pressure(r)
// f[0] = m'(r)
// f[1] = p'(r)
// ......................
int TOVSolver::ODE(double r, const double y[], double f[], void *params)
{
	TOVSolver *tov_obj = (TOVSolver *)params;

	// set a minimum pressure cutoff. if we don't, the ODE solver will wobble all
	// over the surface and crash, or if we make the error tolerance really strict
	// it'll integrate forever
	if (y[1] < tov_obj->PressureCutoff())
	{
#if TOV_SOLVER_VERBOSE
		printf("\u2554----------------- Surface reached ----------------\u2557\n");
		printf("\u2551 %14s %14s %15s %7s\n", "R (km)", "M (M_sun)", "\u03B5_c (g/cm^3)", "\u2551");
		printf("\u2551 %14le %14le %14le %7s\n", r / 1.e+5,
			   y[0] / GSL_CONST_CGSM_SOLAR_MASS, tov_obj->GetEDens(tov_obj->GetInitPress()), "\u2551");
		printf("\u255A--------------------------------------------------\u255D\n");
		std::cout << "\n init_press =" << tov_obj->GetInitPress() << "\n";
		std::cout << "\n Pressure Cut Off =" << tov_obj->PressureCutoff() << "\n";
		std::cout << " y[1] = " << y[1] << "\n";
#endif
		return GSL_EBADFUNC; // this flag tells GSL integrator to quit
	}

	// Mass Continuity Equation
	f[0] = 4. * M_PI * r * r * tov_obj->GetEDens(y[1]);

	// TOV Equation
	//      - [ eps(r) + p/c^2 ]
	// f[1]  = - eos_eps(y[1])
	//       + (y[1] / pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.));

	// //    [ M(r) + 4.pi.r^3 p(r) ]
	// f[1] *= y[0] + 4*M_PI*r*r*r*y[1] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) ;

	// //    / [ r.(r - 2M(r)) ]
	// f[1]  *= 1./r ;
	// f[1]  *= 1./( r/GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT
	//               - 2*y[0]/pow (GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)) ;

	f[1] = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow(r, 2.)) * (tov_obj->GetEDens(y[1]) + (y[1] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))) * (y[0] + 4 * M_PI * pow(r, 3.) * y[1] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)) / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * y[0] / (pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));

	return GSL_SUCCESS;
}

//--------------------------------------------------------------
//               Added on December 15, 2020
/// Returns the derivative of the metric nu(r) function
// r must be in cm!
double TOVSolver::GetNuDer(const double r, const std::vector<double> &y)
{
	// (dp/dr) has units of  [ g / (cm^2 s^2) ]
	double dpdr = -(GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow(r, 2.)) * (GetEDens(y[1]) + (y[1] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.))) * (y[0] + 4 * M_PI * pow(r, 3.) * y[1] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)) / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * y[0] / (pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));

	return -dpdr / (

					   y[1] + pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * GetEDens(y[1])

				   );
}

//--------------------------------------------------------------
double TOVSolver::GetNuDer_Dark(const double r, const std::vector<double> &y)
{

	double out = (GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT / pow(r, 2.)) * (y[0] + 4 * M_PI * pow(r, 3.) * y[1] / pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.)) / (1. - (2. * GSL_CONST_CGSM_GRAVITATIONAL_CONSTANT * y[0] / (pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.) * r)));

	out /= pow(GSL_CONST_CGSM_SPEED_OF_LIGHT, 2.);

	return out;
}

//--------------------------------------------------------------
/// Returns the total baryon number density given pressure
double TOVSolver::GetRho(const double &in_p)
{
	return SafeSplineEval(visi_rho_p_spline, visi_p_accel, in_p, "visible rho(p) spline");
}

//--------------------------------------------------------------
/// Returns the total baryon number density given pressure
double TOVSolver::GetRho_Dark(const double &in_p)
{
	return SafeSplineEval(dark_rho_p_spline, dark_p_accel, in_p, "dark rho(p) spline");
}

//--------------------------------------------------------------
/// Returns the specific number density given pressure
std::vector<double> TOVSolver::GetRho_i(const double &in_p)
{
	std::vector<double> out;
	out.reserve(visi_rho_i_p_spline.size());

	for (auto sp : visi_rho_i_p_spline)
	{
		if (!sp)
		{
			out.push_back(0.0);
			continue;
		} // or NaN if we prefer
		out.push_back(SafeSplineEval(sp, visi_p_accel, in_p, "visible rho_i(p) spline"));
	}

	return out;
}

//--------------------------------------------------------------
/// Returns the specific number density given pressure
std::vector<double> TOVSolver::GetRho_i_Dark(const double &in_p)
{
	std::vector<double> out;
	out.reserve(dark_rho_i_p_spline.size());

	for (auto sp : dark_rho_i_p_spline)
	{
		if (!sp)
		{
			out.push_back(0.0);
			continue;
		} // or NaN if we prefer
		out.push_back(SafeSplineEval(sp, dark_p_accel, in_p, "dark rho_i(p) spline"));
	}

	return out;
}

//--------------------------------------------------------------
double TOVSolver::GetInitPress() const
{
	return init_press;
}

//--------------------------------------------------------------
double TOVSolver::GetInitPress_Dark() const
{
	return init_press_dark;
}

//--------------------------------------------------------------
double TOVSolver::GetInitEDens() const
{
	return init_edens;
}

//--------------------------------------------------------------
// /// Returns the total baryon number density given radius
// double TOVSolver::GetRho_r(const double& in_r)
// {
//   if (results.size() == 0)
//   {
//     Z_LOG_ERROR("Call 'GetRho_r' only after the"
//                 " TOV is fully solved for a star.") ;
//     return 0 ;
//   }

//   std::vector<double> r_list ;
//   std::vector<double> rho_list ;

//   rho_list.reserve(results.size()) ;
//   r_list.reserve(results.size()) ;

//   for (auto &&i : results)
//   {
//     r_list.emplace_back(i.r) ;
//     rho_list.emplace_back(i.rho) ;
//   }

//   rho_r_spline = gsl_spline_alloc (gsl_interp_cspline, results.size());

//   gsl_spline_init (rho_r_spline, &r_list[0], &rho_list[0], results.size()) ;

//   return gsl_spline_eval(rho_r_spline, in_r, accel) ;
// }

//--------------------------------------------------------------
/// Adds the condition for printing the mixed star profile
void TOVSolver::AddMixCondition(bool (*func)(const MixedStar &))
{
	mix_exp_cond_f = func;
}

//--------------------------------------------------------------
/// Adds the condition for printing the star profile
void TOVSolver::AddNCondition(bool (*func)(const NStar &))
{
	n_exp_cond_f = func;
}

//--------------------------------------------------------------
/// Attaches a pointer to analysis
void TOVSolver::AddAnalysis(Analysis *in_analysis)
{
	analysis = in_analysis;
}

//--------------------------------------------------------------
// void TOVSolver::Solve(const Zaki::Math::Axis &in_ax,
// 					  const Zaki::String::Directory &in_dir,
// 					  const Zaki::String::Directory &in_file)
// {
// 	// PROFILE_FUNCTION() ;

// #if TOV_SOLVER_VERBOSE
// 	std::cout << "\n\n\t\t ****************************************"
// 			  << "*******************************" << " \n";
// 	std::cout << "\t\t *                     "
// 			  << "TOV Solver Sequence Results"
// 			  << "                     * \n";
// 	std::cout << "\t\t ******************************************"
// 			  << "*****************************" << "\n\n";
// #endif

// 	// This star object saves the TOV results
// 	// auto tmp_star = std::make_shared<NStar> () ;
// 	// tmp_star->Reserve(1000) ;

// 	for (size_t idx = 0; idx <= in_ax.res; idx++)
// 	{
// 		Z_LOG_INFO("Sequence " + std::to_string(idx + 1) + " out of " + std::to_string(in_ax.res + 1) + ".");

// 		if (idx % 10 == 0)
// 		{
// 			PrintStatus(idx, in_ax.res);
// 		}

// 		// std::cout << "\n\t eps [ " << idx << " ] = "
// 		//           << in_ax[idx] << "\n" ;

// 		// init_edens = in_ax[idx] ;
// 		init_press = p_of_e(in_ax[idx]);

// 		double r = r_min;
// 		double y[2];

// 		y[1] = init_press;
// 		y[0] = (4. / 3.) * M_PI * pow(r, 3.) * GetEDens(y[1]);

// 		// ------------------------------------------
// 		//                RADIUS LOOP
// 		// ------------------------------------------
// 		// ------------------------------------------
// 		// ------------------------------------------
// 		SurfaceIsReached();
// 		// ------------------------------------------

// 		// ----------------------------------------------
// 		if (analysis)
// 			analysis->Analyze(&n_star);
// 		// ----------------------------------------------

// 		// ----------------------------------------------
// 		// Saving the results
// 		// ----------------------------------------------
// 		if (n_exp_cond_f)
// 		{
// 			if (n_exp_cond_f(n_star))
// 			{
// 				ExportNStarProfile(idx, in_dir + "/profiles" + in_file);
// 			}
// 		}
// 		sequence.Add(n_star);
// 		// ----------------------------------------------

// 		// Testing this for freeing memory
// 		// results = std::vector<TOVPoint>();
// 		// results.clear() ;
// 		// results.shrink_to_fit() ;
// 		n_star.Reset();
// 	}
// 	// ----------------------------------------------------------------
// 	//                  TOV Visible sequence loop ends!
// 	// ----------------------------------------------------------------
// #if TOV_SOLVER_VERBOSE
// 	std::cout << "\n\t\t *************************"
// 			  << " TOV Solver Finished *************************" << "\n\n";
// #endif
// 	ExportSequence(in_dir + in_file + "_Sequence.tsv");

// 	// ----------------------------------------------
// 	if (analysis)
// 		analysis->Export(wrk_dir_ + in_dir);
// 	// ----------------------------------------------
// }
//--------------------------------------------------------------
void TOVSolver::Solve(const Zaki::Math::Axis &in_ax,
					  const Zaki::String::Directory &in_dir,
					  const Zaki::String::Directory &in_file)
{
#if TOV_SOLVER_VERBOSE
	std::cout << "\n\n\t\t ****************************************"
			  << "*******************************" << " \n";
	std::cout << "\t\t *                     "
			  << "TOV Solver Sequence Results"
			  << "                     * \n";
	std::cout << "\t\t ******************************************"
			  << "*****************************" << "\n\n";
#endif

	if (eos_tab.eps.empty())
	{
		Z_LOG_ERROR("Solve(...) called but EOS table is empty.");
		return;
	}

	// ------------------------------------------------------------------
	// ADR-0005 (ACCEPTED 2026-09-02): this function is a SEQUENCE / WORKFLOW
	// ORCHESTRATOR, not a numerical authority. The canonical TOV numerical
	// primitive is SingleStarSolveToTOVPoints, and it owns:
	//   the central-density clamp, the p_of_e conversion, the initial
	//   conditions, the GSL RK8PD driver, the radial grid, the step_scale
	//   ladder, the pressure-cutoff termination and TOVPoint construction.
	//
	// What remains here is orchestration only: sweeping the Axis, driving the
	// existing star finalization, the Analysis and export hooks, Sequence
	// accumulation, and the unconditional _Sequence.tsv export.
	//
	// The clamp is deliberately NOT duplicated here any more. The primitive
	// computes floor/ceiling from the same eos_tab with the same
	// central_eps_floor_factor and the same 0.999 margin, so the effective
	// central density is unchanged. The only observable difference is the log
	// line: the warning now originates in SingleStarSolveToTOVPoints and reads
	// "Requested eps(...)" rather than "Solve: requested eps(...)". No
	// numerical behavior changes. See docs/validation/PHASE3E_I1_CANONICAL_TOV.md.
	//
	// Phase 3E-0 measured this delegation as exactly equivalent: all 25 radial
	// columns bit-identical at 14 central densities and radial_res
	// 5000/10000/20000 (docs/validation/TOV_PATH_EQUIVALENCE.md).
	// ------------------------------------------------------------------
	for (size_t idx = 0; idx <= in_ax.res; idx++)
	{
		Z_LOG_INFO("Sequence " + std::to_string(idx + 1) +
				   " out of " + std::to_string(in_ax.res + 1) + ".");

		if (idx % 10 == 0)
			PrintStatus(idx, in_ax.res);

		// The ONE radial integration for this sequence member.
		std::vector<TOVPoint> tov_points;
		SingleStarSolveToTOVPoints(in_ax[idx], tov_points);

		// Feed the points through the EXISTING Path-1 postprocessing,
		// unchanged. Converging this onto BuildFromTOV is ADR-0005 Q3 = P3 and
		// is deliberately NOT part of this increment.
		for (const auto &tp : tov_points)
			n_star.Append(tp);

		SurfaceIsReached();

		if (analysis)
			analysis->Analyze(&n_star);

		if (n_exp_cond_f && n_exp_cond_f(n_star))
			ExportNStarProfile(idx, in_dir + "/profiles" + in_file);

		sequence.Add(n_star);
		n_star.Reset();
	}

#if TOV_SOLVER_VERBOSE
	std::cout << "\n\t\t *************************"
			  << " TOV Solver Finished *************************" << "\n\n";
#endif

	ExportSequence(in_dir + in_file + "_Sequence.tsv");

	if (analysis)
		analysis->Export(wrk_dir_ + in_dir);
}

//--------------------------------------------------------------
// The radius iteration in the neutron star scenario
void TOVSolver::Solve_Mixed(const Zaki::Math::Axis &in_v_ax,
							const Zaki::Math::Axis &in_d_ax,
							const Zaki::String::Directory &in_dir,
							const Zaki::String::Directory &in_file)
{
	PROFILE_FUNCTION();

#if TOV_SOLVER_VERBOSE
	std::cout << "\n\n\t\t ****************************************"
			  << "*******************************" << " \n";
	std::cout << "\t\t *                     "
			  << "TOV Solver Sequence Results"
			  << "                     * \n";
	std::cout << "\t\t ******************************************"
			  << "*****************************" << "\n\n";
#endif

	// ----------------------------------------------------------------
	//                  TOV Dark sequence loop begins
	// ----------------------------------------------------------------
	for (size_t d_idx = 0; d_idx <= in_d_ax.res; d_idx++)
	{
		Z_LOG_INFO("Dark sequence " + std::to_string(d_idx + 1) + " out of " + std::to_string(in_d_ax.res + 1) + ".");

		init_press_dark = p_of_e_dark(in_d_ax[d_idx]);
		// ----------------------------------------------------------------
		//                  TOV Visible sequence loop begins
		// ----------------------------------------------------------------
		for (size_t v_idx = 0; v_idx <= in_v_ax.res; v_idx++)
		{
			if (c_poly.IsExcluded({in_v_ax[v_idx], in_d_ax[d_idx]}))
			{
				ignored_counter++;
				continue;
			}

			Z_LOG_INFO("Mixed sequence (" + std::to_string(d_idx + 1) + ", " + std::to_string(v_idx + 1) + ") out of (" + std::to_string(in_d_ax.res + 1) + ", " + std::to_string(in_v_ax.res + 1) + ").");

			if (v_idx % 10 == 0)
			{
				PrintStatus(v_idx, d_idx, in_v_ax.res, in_d_ax.res);
			}

			init_press = p_of_e(in_v_ax[v_idx]);

			double r = r_min;
			double y[4], y_mantle[2];

			y[2] = init_press;
			y[3] = init_press_dark;
			y[0] = (4. / 3.) * M_PI * pow(r, 3.) * GetEDens(y[2]);
			y[1] = (4. / 3.) * M_PI * pow(r, 3.) * GetEDens_Dark(y[3]);

			// ------------------------------------------
			//                RADIUS LOOP
			// ------------------------------------------
			RadiusLoopMixed(r, y, y_mantle);
			// ------------------------------------------

			// ----------------------------------------------
			// mixed_star.SurfaceIsReached(v_idx, d_idx) ;
			SurfaceIsReached(v_idx, d_idx);
			// ----------------------------------------------

			// ----------------------------------------------
			if (analysis)
				analysis->Analyze(&mixed_star);
			// ----------------------------------------------

#if DarkCore_Analysis
			double m_n = Zaki::Physics::NEUTRON_M_FM;
			double m_chi = 0.8 * Zaki::Physics::NEUTRON_M_FM;

			double critical_rho_val = pow(m_n - m_chi * m_chi / m_n, 3) / (24 * M_PI * M_PI);
			int critical_idx = mixed_star.ds_dar[mixed_star.rho_idx].GetClosestIdx(critical_rho_val);
			std::cout << "Critical Rho = " << critical_rho_val
					  << ", Index = " << critical_idx
					  << ", Size = "
					  << mixed_star.ds_dar.Dim()[0]
					  << ", rho_d [v_idx-1] = " << mixed_star.ds_dar[mixed_star.rho_idx][critical_idx - 1]
					  << ", rho_d [v_idx] = " << mixed_star.ds_dar[mixed_star.rho_idx][critical_idx]
					  << "\n\t r [v_idx-1] = " << mixed_star.ds_dar[mixed_star.r_idx][critical_idx - 1]
					  << ", r [v_idx] = " << mixed_star.ds_dar[mixed_star.r_idx][critical_idx] << "\n";

			// std::map<std::string, std::string> baryons = {
			//   {"10", "Neutron", 6},
			//   {"100", "Lambda", 7}, {"110", "Sigma-", 10}, {"111", "Sigma0"},
			//   {"112", "Sigma+"}, {"120","Xi-"}, {"121", "Xi0"}
			// } ;

			Zaki::Vector::DataColumn b_den = mixed_star.ds_vis[mixed_star.rho_i_v_idx[7]] * mixed_star.ds_vis[mixed_star.rho_idx];
			Zaki::Vector::DataColumn r_set = mixed_star.ds_vis[mixed_star.r_idx];
			Zaki::Vector::DataColumn M_r = mixed_star.ds_vis[mixed_star.m_idx];
			Zaki::Vector::DataColumn nu_r = mixed_star.ds_vis[mixed_star.nu_idx];

			Zaki::Vector::DataColumn integ_fr = 4 * M_PI * b_den;
			integ_fr *= r_set.pow(2) * (1. - 2 * M_r / r_set).pow(-0.5);

			Zaki::Vector::DataColumn integ_b_dot = integ_fr;
			integ_b_dot *= exp(nu_r);

			Zaki::Vector::DataSet integrand({r_set, integ_fr, integ_b_dot});

			std::cout << "\n integ_fr.Size() = " << integ_fr.Size();
			std::cout << "\n integ_b_dot.Size() = " << integ_b_dot.Size();

			integrand.Interpolate(0, {1, 2});
			double fr_result = integrand.Integrate(1, {r_set[0], r_set[-1]});
			double b_dot_result = integrand.Integrate(2, {r_set[0], r_set[-1]});

			double fr_result_choked = integrand.Integrate(1, {r_set[critical_idx], r_set[-1]});
			double b_dot_result_choked = integrand.Integrate(2, {r_set[critical_idx], r_set[-1]});

			std::cout << "\nFraction = " << fr_result * 1e54 / mixed_star.sequence.v.b
					  << ", [B_dot/B] = "
					  << b_dot_result * 1e54 / mixed_star.sequence.v.b
					  << "\n Fraction choked = "
					  << fr_result_choked * 1e54 / mixed_star.sequence.v.b
					  << ", [B_dot/B] choked = "
					  << b_dot_result_choked * 1e54 / mixed_star.sequence.v.b << "\n\n";
#endif
			// ----------------------------------------------
			// Saving the results
			// ----------------------------------------------
			// mixed_star.SetWrkDir( wrk_dir_ ) ;
			if (mix_exp_cond_f)
			{
				// ............ Checking/Making "profiles" directory ............
				// Not needed anymore? Commented on [ Apr 20, 2023 ]
				// (wrk_dir_ + in_dir + "/profiles").Create() ;
				// if (!std::filesystem::is_directory((in_dir + "/profiles").Str()))
				// {
				//   if (mkdir((wrk_dir_ + in_dir + "/profiles").Str().c_str(),
				//               ACCESSPERMS) == -1)
				//   {
				//     Z_LOG_INFO("Directory '"
				//                 + (wrk_dir_ + in_dir + "/profiles").Str()
				//                 + "' wasn't created, because: "
				//                 + strerror(errno)+".") ;
				//   }
				//   else
				//     Z_LOG_INFO("Directory '"
				//                 + (wrk_dir_ + in_dir + "/profiles").Str()
				//                 + "' created.") ;
				// }
				// ..............................................................
				if (mix_exp_cond_f(mixed_star))
					ExportMixedStarProfile(v_idx, d_idx, in_dir + "/profiles" + in_file);
				// mixed_star.Export(in_dir + "/Mixed_" +
				//     std::to_string(d_idx) + "_" +
				//     std::to_string(v_idx) + ".tsv") ;
			}
			mixed_sequence.Add(mixed_star);
			// ----------------------------------------------

#if 0
      // ----------------------
      sequence.emplace_back(
                      GetEDens(init_press), 
                      tmp_m_star->mass[-1]/Zaki::Physics::SUN_M_KM,
                      tmp_m_star->radius[-1], init_press, 
                      tmp_m_star->GetBaryonNum(), 
                      tmp_m_star->GetMomInertia()
                            );
      // ----------------------------------------------

      std::vector<std::string> tmp_fname_v = Zaki::String::Pars(in_file.Str(), "*") ;
      std::string tmp_fname ;
      if(tmp_fname_v.size() <2 )
      {
        Z_LOG_WARNING("File name pattern doesn't match '[]*[]'.") ;
        tmp_fname = tmp_fname_v[0] + std::to_string(v_idx) ;
      }
      else
      {
        tmp_fname = tmp_fname_v[0] + std::to_string(v_idx) + tmp_fname_v[1] ;
      }
      
    // ............ Creating a directory ............
    if (mkdir((wrk_dir_ + in_dir + "/"
                      + std::to_string(d_idx)
                      ).Str().c_str(), ACCESSPERMS) == -1) 
    {
      Z_LOG_INFO("Directory '"+ (wrk_dir_ + in_dir + "/" 
                                      + std::to_string(d_idx)
                                      ).Str()
                                      +"' wasn't created, because: "
                                      +strerror(errno)+".") ;    
    }
    else
      Z_LOG_INFO("Directory '"+(wrk_dir_ + in_dir + "/" 
                                      + std::to_string(d_idx) 
                                      ).Str()
                                      +"' created.") ; 
    // .................................................

      Zaki::File::VecSaver vec_saver(wrk_dir_ + in_dir + "/" 
                                      + std::to_string(d_idx)
                                      + "/" + tmp_fname) ;

      // -------------------- Visible Header Begins --------------------
      char res_header[200] ;
      snprintf(res_header, sizeof(res_header), "%-14s\t %-14s\t %-14s\t %-15s\t %-14s\t %-14s\t %-14s", 
              "r(km)", "m", "nu'(1/cm)", " nu", "p(dyne/cm^2)", "e(g/cm^3)", "rho(1/fm^3)" ) ;

      std::string tmp_label = res_header ;
      for(auto& lab : eos_tab.extra_labels)
      {
        snprintf(res_header, sizeof(res_header), "\t %-14s", Zaki::String::Strip(lab, ' ').c_str() ) ;
        tmp_label += res_header ;
      }
      // -------------------- Visible Header Ends--------------------

      vec_saver.SetHeader(tmp_label.c_str()) ;
      vec_saver.Export1D(results) ;
      
      // -------------------- Dark Header Begins --------------------
      snprintf(res_header, sizeof(res_header), "%-14s\t %-14s\t %-14s\t %-15s\t %-14s\t %-14s\t %-14s", 
              "r(km)", "m", "nu'(1/cm)", " nu", "p(dyne/cm^2)", "e(g/cm^3)", "rho(1/fm^3)" ) ;

      tmp_label = res_header ;
      for(auto& lab : eos_tab_dark.extra_labels)
      {
        snprintf(res_header, sizeof(res_header), "\t %-14s", Zaki::String::Strip(lab, ' ').c_str() ) ;
        tmp_label += res_header ;
      }
      // -------------------- Dark Header Ends--------------------
      vec_saver.SetHeader(tmp_label.c_str()) ;

      tmp_fname_v = Zaki::String::Pars(tmp_fname, ".") ;
      if(tmp_fname_v.size() == 2)
        vec_saver.SetFileName(wrk_dir_ + in_dir +  "/" 
                                      + std::to_string(d_idx)
                                      + "/" + tmp_fname_v[0] 
                                      + "_dark." +  tmp_fname_v[1] ) ;
      else
      {
        vec_saver.SetFileName(wrk_dir_ + in_dir + "/" 
                                      + std::to_string(d_idx)
                                      + "/" + tmp_fname + "_dark" ) ;
      }
      
      vec_saver.Export1D(results_dark) ;

#endif

			mixed_star.Reset();

#if TOV_SOLVER_VERBOSE
			printf("\n========================== %s (%lu) "
				   "==========================\n\n",
				   "End of Seq.", v_idx + 1);
#endif
		}
		// ----------------------------------------------------------------
		//                  TOV Visible sequence loop ends!
		// ----------------------------------------------------------------
#if TOV_SOLVER_VERBOSE
		std::cout << "\n\t\t *************************"
				  << " TOV Solver Finished *************************" << "\n\n";
#endif
	}
	// ----------------------------------------------------------------
	//                  TOV Dark sequence loop ends!
	// ----------------------------------------------------------------
	ExportMixedSequence(in_dir + in_file + "_Sequence.tsv");
	std::cout << "\n\t Number of points ignored: " << ignored_counter << "\n";
	// mixed_sequence.Export(in_dir + "/Mixed_Sequence.tsv") ;

	// ----------------------------------------------
	if (analysis)
		analysis->Export(wrk_dir_ + in_dir);
	// ----------------------------------------------
}

//--------------------------------------------------------------
void TOVSolver::Solve_Mixed(const Contour &eps_cont,
							const Zaki::String::Directory &in_dir,
							const Zaki::String::Directory &in_file)
{
	PROFILE_FUNCTION();

#if TOV_SOLVER_VERBOSE
	std::cout << "\n\n\t\t ****************************************"
			  << "*******************************" << " \n";
	std::cout << "\t\t *                     "
			  << "TOV Solver Sequence Results"
			  << "                     * \n";
	std::cout << "\t\t ******************************************"
			  << "*****************************" << "\n\n";
#endif

	// std::vector<BNV_Rate> bnv_rates ;
	// bnv_rates.reserve(eps_cont.Size());
	// ----------------------------------------------------------------
	//                  TOV sequence loop begins
	// ----------------------------------------------------------------
	for (size_t i = 0; i < eps_cont.Size(); i += 10)
	{
		Z_LOG_INFO("Mixed sequence " + std::to_string(i + 1) + " out of " + std::to_string(eps_cont.Size()) + ".");

		PrintStatus(i, 0, eps_cont.Size(), 0);

		// ------------------ WHILE Loop Starts---------------------------
		std::cout << "\n";
		double tmp_b_tot = 0;
		size_t trial = 0;
		double tmp_e_d = eps_cont.curve[i].y;
		double tmp_e_v = eps_cont.curve[i].x;
		bool b_excess = false;
		bool b_under = false;
		size_t cross = 0;
		double best_rel_err = 1;
		size_t best_idx = 0;
		std::vector<eps_pair> trial_set;
		trial_set.reserve(eps_cont.max_steps);
		bool exp_decrease = false;
		int flip_sign = 1;
		do
		{
			trial++;
			if (trial != 1)
			{
				mixed_star.Reset();
			}

			trial_set.emplace_back(tmp_e_v, tmp_e_d);

			init_press_dark = p_of_e_dark(tmp_e_d);
			init_press = p_of_e(tmp_e_v);

			double r = r_min;
			double y[4], y_mantle[2];

			y[2] = init_press;
			y[3] = init_press_dark;
			y[0] = (4. / 3.) * M_PI * pow(r, 3.) * GetEDens(y[2]);
			y[1] = (4. / 3.) * M_PI * pow(r, 3.) * GetEDens_Dark(y[3]);

			// ------------------------------------------
			//                RADIUS LOOP
			// ------------------------------------------
			RadiusLoopMixed(r, y, y_mantle);
			// ------------------------------------------
			// ------------------------------------------
			SurfaceIsReached(i, 0);
			// ------------------------------------------
			// if(trial != 1)
			// {
			//   if( exp_decrease &&
			//      (mixed_star.sequence.v.b + mixed_star.sequence.d.b) > tmp_b_tot)
			//     {
			//       flip_sign *= -1 ;
			//       std::cout << "Signed flipped!\n" ;
			//     }
			//   if( !exp_decrease &&
			//      (mixed_star.sequence.v.b + mixed_star.sequence.d.b) < tmp_b_tot)
			//     {
			//       flip_sign *= -1 ;
			//       std::cout << "Signed flipped!\n" ;
			//     }
			// }

			tmp_b_tot = mixed_star.sequence.v.b + mixed_star.sequence.d.b;
			double tmp_rel_err = abs(tmp_b_tot - eps_cont.val) / eps_cont.val;
			std::cout.precision(9);
			std::cout << "Trial = " << trial
					  << ": B_tot = " << tmp_b_tot
					  << ", val = " << eps_cont.val
					  << ", rel_err = "
					  << tmp_rel_err
					  << "\n";

			if (tmp_rel_err < best_rel_err)
			{
				best_rel_err = tmp_rel_err;
				best_idx = trial - 1;
			}
			double step_size = 7e-3;

			if (tmp_rel_err < 1e-7)
				step_size = 1e-6;
			else if (tmp_rel_err < 3e-6)
				step_size = 1e-5;
			else if (tmp_rel_err < 1e-5)
				step_size = 0.5e-4;
			else if (tmp_rel_err < 1e-4)
				step_size = 1e-4;
			else if (tmp_rel_err < 5e-4)
				step_size = 3e-3;
			else if (tmp_rel_err < 1e-3)
				step_size = 5e-3;
			// To decrease B, increase e_d and decrease e_v
			if (tmp_b_tot > eps_cont.val)
			{
				exp_decrease = true;
				b_excess = true;
				if (b_under)
				{
					cross++;
					b_under = false;
				}
				tmp_e_d *= 1 + flip_sign * (step_size) / (1. + (1 - tmp_rel_err) * cross);
				tmp_e_v *= 1 - flip_sign * (step_size) / (1. + (1 - tmp_rel_err) * cross);
			}
			// To increase B, decrease e_d and increase e_v
			else if (tmp_b_tot < eps_cont.val)
			{
				exp_decrease = false;
				b_under = true;
				if (b_excess)
				{
					cross++;
					b_excess = false;
				}
				tmp_e_d *= 1 - flip_sign * (step_size) / (1. + (1 - tmp_rel_err) * cross);
				tmp_e_v *= 1 + flip_sign * (step_size) / (1. + (1 - tmp_rel_err) * cross);
			}
		} while (abs(tmp_b_tot - eps_cont.val) / eps_cont.val > eps_cont.precision && trial <= eps_cont.max_steps);
		// ------------------ WHILE Loop Ends ----------------------------
		if (best_idx != trial - 1)
		{
			mixed_star.Reset();
			trial++;
			std::cout << "Best Rel. Err. = " << best_rel_err << "\n";
			tmp_e_d = trial_set[best_idx].y;
			tmp_e_v = trial_set[best_idx].x;

			init_press_dark = p_of_e_dark(tmp_e_d);
			init_press = p_of_e(tmp_e_v);

			double r = r_min;
			double y[4], y_mantle[2];

			y[2] = init_press;
			y[3] = init_press_dark;
			y[0] = (4. / 3.) * M_PI * pow(r, 3.) * GetEDens(y[2]);
			y[1] = (4. / 3.) * M_PI * pow(r, 3.) * GetEDens_Dark(y[3]);

			// ------------------------------------------
			//                RADIUS LOOP
			// ------------------------------------------
			RadiusLoopMixed(r, y, y_mantle);
			// ------------------------------------------
			// ------------------------------------------
			SurfaceIsReached(i, 0);
			// ------------------------------------------
			tmp_b_tot = mixed_star.sequence.v.b + mixed_star.sequence.d.b;
			double tmp_rel_err = abs(tmp_b_tot - eps_cont.val) / eps_cont.val;
			std::cout.precision(9);
			std::cout << "Trial(*) = " << trial
					  << ": B_tot = " << tmp_b_tot
					  << ", val = " << eps_cont.val
					  << ", rel_err = "
					  << tmp_rel_err
					  << "\n";
		}

		trial_set.clear();

		// ----------------------------------------------
		if (analysis)
			analysis->Analyze(&mixed_star);
		// ----------------------------------------------

#if DarkCore_Analysis

		//   double m_n = Zaki::Physics::NEUTRON_M_FM ;
		//   double m_chi = 0.8*Zaki::Physics::NEUTRON_M_FM ;

		//   double critical_rho_val = pow(m_n - m_chi*m_chi/m_n, 3)
		//                             / ( 24 * M_PI*M_PI ) ;
		//   int critical_idx = mixed_star.ds_dar[
		//                           mixed_star.rho_idx
		//                           ].GetClosestIdx(critical_rho_val) ;
		//   // std::cout << "Critical Rho = " << critical_rho_val
		//   //           << ", Index = " << critical_idx
		//   //           << ", Size = "
		//   //           << mixed_star.ds_dar.Dim()[0]
		//   //           << ", rho_d [v_idx-1] = " << mixed_star.ds_dar[
		//   //                         mixed_star.rho_idx
		//   //                         ][critical_idx-1]
		//   //           << ", rho_d [v_idx] = " << mixed_star.ds_dar[
		//   //                         mixed_star.rho_idx
		//   //                         ][critical_idx]
		//   //           << "\n\t r [v_idx-1] = " << mixed_star.ds_dar[
		//   //                         mixed_star.r_idx
		//   //                         ][critical_idx-1]
		//   //           << ", r [v_idx] = " << mixed_star.ds_dar[
		//   //                         mixed_star.r_idx
		//   //                         ][critical_idx] << "\n" ;

		// // std::map<std::string, std::string> baryons = {
		// //   {"10", "Neutron", 6},
		// //   {"100", "Lambda", 7}, {"110", "Sigma-", 10}, {"111", "Sigma0"},
		// //   {"112", "Sigma+"}, {"120","Xi-"}, {"121", "Xi0"}
		// // } ;

		//   Zaki::Vector::DataColumn b_den = mixed_star.ds_vis[
		//                                     mixed_star.rho_i_v_idx[6]
		//                                     ] * mixed_star.ds_vis[
		//                                     mixed_star.rho_idx
		//                                     ];
		//   Zaki::Vector::DataColumn r_set = mixed_star.ds_vis[
		//                                     mixed_star.r_idx
		//                                     ];
		//   Zaki::Vector::DataColumn M_r = mixed_star.ds_vis[
		//                                     mixed_star.m_idx
		//                                     ];
		//   Zaki::Vector::DataColumn nu_r = mixed_star.ds_vis[
		//                                     mixed_star.nu_idx
		//                                     ];

		//   Zaki::Vector::DataColumn integ_fr = 4*M_PI*b_den ;
		//   integ_fr *= r_set.pow(2) * (1. - 2*M_r / r_set ).pow(-0.5) ;

		//   Zaki::Vector::DataColumn integ_b_dot = integ_fr ;
		//   integ_b_dot *= exp( nu_r ) ;

		//   Zaki::Vector::DataSet integrand({r_set, integ_fr, integ_b_dot}) ;

		//   // std::cout << "\n integ_fr.Size() = " << integ_fr.Size() ;
		//   // std::cout << "\n integ_b_dot.Size() = " << integ_b_dot.Size() ;

		//   integrand.Interpolate(0, {1, 2}) ;
		//   double fr_result = integrand.Integrate(1, {r_set[0], r_set[-1]}) ;
		//   double b_dot_result = integrand.Integrate(2, {r_set[0], r_set[-1]}) ;

		//   double fr_result_choked = integrand.Integrate(1, {r_set[critical_idx], r_set[-1]}) ;
		//   double b_dot_result_choked = integrand.Integrate(2, {r_set[critical_idx], r_set[-1]}) ;

		//   // std::cout << "\nFraction = " << fr_result* 1e54 / mixed_star.sequence.v.b
		//   //           << ", [B_dot/B] = "
		//   //           << b_dot_result* 1e54 / mixed_star.sequence.v.b
		//   //           <<"\n Fraction choked = "
		//   //           <<  fr_result_choked* 1e54 / mixed_star.sequence.v.b
		//   //           << ", [B_dot/B] choked = "
		//   //           << b_dot_result_choked* 1e54 / mixed_star.sequence.v.b  << "\n\n" ;
		//   std::cout.precision(9) ;
		//   std::cout << "\n B_dot = "
		//             << b_dot_result* 1e54
		//             << ", B_dot (choked) = "
		//             << b_dot_result_choked* 1e54
		//             <<", B_vis = " << mixed_star.sequence.v.b
		//             << ", B_tot = " << mixed_star.sequence.v.b
		//                 + mixed_star.sequence.d.b
		//             << "\n\n" ;
		//   bnv_rates.emplace_back( b_dot_result* 1e54,
		//                           b_dot_result_choked* 1e54,
		//                           mixed_star.sequence.v.b) ;
#endif
		// ----------------------------------------------
		// Saving the results
		// ----------------------------------------------
		// mixed_star.SetWrkDir( wrk_dir_ ) ;
		if (mix_exp_cond_f)
		{
			if (mix_exp_cond_f(mixed_star))
			{
				// ............ Checking/Making "profiles" directory ............
				// Not needed anymore? Commented on [ Apr 20, 2023 ]
				// (wrk_dir_ + in_dir + "/profiles").Create() ;
				// if (!std::filesystem::is_directory((in_dir + "/profiles").Str()))
				// {
				//   if (mkdir((wrk_dir_ + in_dir + "/profiles").Str().c_str(),
				//               ACCESSPERMS) == -1)
				//   {
				//     Z_LOG_INFO("Directory '"
				//                 + (wrk_dir_ + in_dir + "/profiles").Str()
				//                 + "' wasn't created, because: "
				//                 + strerror(errno)+".") ;
				//   }
				//   else
				//     Z_LOG_INFO("Directory '"
				//                 + (wrk_dir_ + in_dir + "/profiles").Str()
				//                 + "' created.") ;
				// }
				// ..............................................................
				ExportMixedStarProfile(i, 0, in_dir + "/profiles" + in_file);
			}
			// mixed_star.Export(in_dir + "/Mixed_" +
			//     std::to_string(d_idx) + "_" +
			//     std::to_string(v_idx) + ".tsv") ;
		}
		mixed_sequence.Add(mixed_star);
		// ----------------------------------------------

		mixed_star.Reset();

#if TOV_SOLVER_VERBOSE
		printf("\n========================== %s (%lu) "
			   "==========================\n\n",
			   "End of Seq.", i + 1);
#endif
	}
	// ----------------------------------------------------------------
	//                  TOV Visible sequence loop ends!
	// ----------------------------------------------------------------
#if TOV_SOLVER_VERBOSE
	std::cout << "\n\t\t *************************"
			  << " TOV Solver Finished *************************" << "\n\n";
#endif

	ExportMixedSequence(in_dir + in_file + "_Sequence.tsv");

	// ----------------------------------------------
	if (analysis)
		analysis->Export(wrk_dir_ + in_dir);
	// ----------------------------------------------
	// Zaki::File::VecSaver vec_saver(wrk_dir_ + in_dir + "/BNV_rates.tsv") ;
	// char bnv_header[100] ;
	// snprintf(bnv_header, sizeof(bnv_header), "%-14s\t %-14s\t %-14s",
	//             "B_dot", "B_dot_choke", "B" ) ;
	// vec_saver.SetHeader(bnv_header) ;
	// vec_saver.Export1D(bnv_rates) ;
}
//--------------------------------------------------------------
// Single-star TOV solve → vector<TOVPoint>
//--------------------------------------------------------------
int TOVSolver::SingleStarSolveToTOVPoints(double ec_central,
										  std::vector<TOVPoint> &out_tov)
{
	PROFILE_FUNCTION();

	out_tov.clear();

	if (eos_tab.eps.empty())
	{
		Z_LOG_ERROR("EOS table is empty.");
		return 0;
	}

	// ----------------------------------------------------------
	// 0) Clamp central energy density to EOS range (same as Solve)
	// ----------------------------------------------------------
	const double eos_e_min = eos_tab.eps.front();
	const double eos_e_max = eos_tab.eps.back();

	const double floor_e = central_eps_floor_factor * eos_e_min;
	const double ceil_e = 0.999 * eos_e_max;

	const double ec_req = ec_central;
	double ec = ec_req;

	if (ec < floor_e)
	{
		Z_LOG_WARNING("Requested eps(" +
					  std::to_string(ec_req) + ") < floor(" +
					  std::to_string(floor_e) + ") -> clamping.");
		ec = floor_e;
	}
	else if (ec > ceil_e)
	{
		Z_LOG_WARNING("Requested eps(" +
					  std::to_string(ec_req) + ") > ceil(" +
					  std::to_string(ceil_e) + ") -> clamping.");
		ec = ceil_e;
	}

	// ----------------------------------------------------------
	// 1) Convert central ε to central pressure, set initial y[]
	// ----------------------------------------------------------
	init_press = p_of_e(ec);

	double r = r_min; // cm
	double y[2];

	// y[1] = p(r), y[0] = m(r) in cgs (g)
	y[1] = init_press;
	y[0] = (4.0 / 3.0) * M_PI * std::pow(r, 3.0) * GetEDens(y[1]);

	// ----------------------------------------------------------
	// 2) GSL ODE setup -- the ONE ordinary-star radial driver (ADR-0005)
	// ----------------------------------------------------------
	gsl_odeiv2_system ode_sys = {TOVSolver::ODE, nullptr, 2, this};

	gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new(
		&ode_sys,
		gsl_odeiv2_step_rk8pd,
		1.e-1,	// initial step
		1.e-10, // abs tol
		1.e-10	// rel tol
	);

	double min_log_r = r_min;
	double max_log_r = r_max;

	double step = (max_log_r - min_log_r) / radial_res;
	double step_scale = 1.0;

	const double p_cut = PressureCutoff();

	// ----------------------------------------------------------
	// 3) Radius loop — the canonical ordinary-star radial integration.
	//    Historically this was a copy of TOVSolver::RadiusLoop; that duplicate was retired
	//    in Phase 3E-I4 and this is now the sole implementation.
	// ----------------------------------------------------------
	for (double log_r_i = min_log_r;
		 log_r_i <= max_log_r;
		 log_r_i += step * step_scale)
	{
		double ri = log_r_i;

		double tmp_delta_p = y[1]; // kept for potential future step-control tweaks

		int status = gsl_odeiv2_driver_apply(driver, &r, ri, y);

		if (status != GSL_SUCCESS)
		{
#if TOV_SOLVER_VERBOSE
			printf("\t-------------------%s-------------------\n", "GSL");
			printf("error, return value=%d\n.", status);
			printf("Pressure = %2.2e.\n", y[1]);
#endif
			break;
		}

		// ------------------------------------------------------
		// Step scaling
		// r is in cm; thresholds 100, 1000, ... are cm as well.
		// ------------------------------------------------------
		if (ri < 100.0) // < 1 m
		{
			step_scale = 0.005;
		}
		else if (ri < 1000.0) // 1 m - 10 m
		{
			step_scale = 0.025;
		}
		else if (ri < 10000.0) // 10 m - 100 m
		{
			step_scale = 0.05;
		}
		else if (ri < 100000.0) // 100 m - 1 km
		{
			step_scale = 0.25;
		}
		else // > 1 km
		{
			step_scale = 1.0;
		}

		// ------------------------------------------------------
		// Build and store TOVPoint at this radius
		//
		// Conventions:
		//  - r in km         → r / 1e5 (cm → km)
		//  - m in solar mass → y[0] / GSL_CONST_CGSM_SOLAR_MASS
		//  - ν' from GetNuDer
		//  - ν = 0 (rebuilt later by NStar::BuildFromTOV)
		//  - p as is (cgs)
		//  - e = GetEDens(p)
		//  - ρ = GetRho(p)
		//  - ρ_i = GetRho_i(p) (vector)
		// ------------------------------------------------------
		const double r_km = r / 1.e5;
		const double m_msun = y[0] / GSL_CONST_CGSM_SOLAR_MASS;

		const double nu_prime = GetNuDer(r, {y[0], y[1]});
		const double e_here = GetEDens(y[1]);
		const double rho_here = GetRho(y[1]);
		const std::vector<double> rho_i_here = GetRho_i(y[1]);

		// ADR-0007 P5: the authoritative barotropic derivative, taken from the same eps(p)
		// interpolant that produced `e_here` one line above. Dimensionless. NaN if the
		// pressure sits outside the interpolant's domain, which is exactly the state in
		// which no derivative may be claimed; a profile carrying any NaN simply publishes
		// no EOS-derivative data and every consumer of it fails closed.
		const double dedp_here = GetEDensDeriv(y[1]);

		out_tov.emplace_back(
			r_km,
			m_msun,
			nu_prime,
			0.0,	    // ν will be reconstructed later from ν'(r)
			y[1],	    // p
			e_here,	    // ε
			rho_here,   // total baryon density
			rho_i_here, // species densities
			dedp_here   // d(eps)/dp, dimensionless (ADR-0007 P5)
		);

		// ------------------------------------------------------
		// Termination condition: pressure below cutoff
		// (same notion as used throughout TOVSolver)
		// ------------------------------------------------------
		if (y[1] <= p_cut)
			break;
	}

	gsl_odeiv2_driver_free(driver);

	if (out_tov.empty())
		return 0;

	return static_cast<int>(out_tov.size());
}

//--------------------------------------------------------------
// SolveToProfile: Find e_c such that M(e_c) ≈ target_M_solar
//   using root-finding in central energy density.
//--------------------------------------------------------------
int TOVSolver::SolveToProfile(double target_M_solar,
							  std::vector<TOVPoint> &out_tov,
							  std::vector<std::string> *out_species_labels)
{
	PROFILE_FUNCTION();

	// (void)model_name; // currently unused; EOS assumed already imported

	out_tov.clear();

	if (eos_tab.eps.empty())
	{
		Z_LOG_ERROR("EOS table is empty. Call ImportEOS(...) first.");
		return 0;
	}

	// ----------------------------------------------------------
	// 0) Define allowed central ε range (same logic as Solve)
	// ----------------------------------------------------------
	const double eos_e_min = eos_tab.eps.front();
	const double eos_e_max = eos_tab.eps.back();

	const double floor_e = central_eps_floor_factor * eos_e_min;
	const double ceil_e = 0.999 * eos_e_max;

	if (floor_e <= 0.0 || ceil_e <= floor_e)
	{
		Z_LOG_ERROR("SolveToProfile: invalid EOS energy-density range.");
		return 0;
	}

	// ----------------------------------------------------------
	// 1) Coarse sampling in log ε_c to bracket the target mass
	// ----------------------------------------------------------
	const int N_coarse = 24; // number of coarse samples (configurable)

	std::vector<double> ec_grid;
	std::vector<double> M_grid;
	ec_grid.reserve(N_coarse + 1);
	M_grid.reserve(N_coarse + 1);

	int best_idx = -1;
	double best_mass_diff = std::numeric_limits<double>::infinity();

	const double log_e_lo = std::log10(floor_e);
	const double log_e_hi = std::log10(ceil_e);

	for (int i = 0; i <= N_coarse; ++i)
	{
		const double t = static_cast<double>(i) / static_cast<double>(N_coarse);
		const double log_e = log_e_lo + t * (log_e_hi - log_e_lo);
		const double ec = std::pow(10.0, log_e);

		std::vector<TOVPoint> tmp;
		const int npts = SingleStarSolveToTOVPoints(ec, tmp);

		if (npts <= 0 || tmp.empty())
		{
			Z_LOG_ERROR("SolveToProfile: SingleStarSolveToTOVPoints failed at ec = " +
						std::to_string(ec));
			return 0;
		}

		const double M_here = tmp.back().m; // Msun (by construction in SingleStarSolveToTOVPoints)

		ec_grid.push_back(ec);
		M_grid.push_back(M_here);

		const double diff = std::fabs(M_here - target_M_solar);
		if (diff < best_mass_diff)
		{
			best_mass_diff = diff;
			best_idx = static_cast<int>(ec_grid.size()) - 1;
		}
	}

	if (ec_grid.empty())
	{
		Z_LOG_ERROR("SolveToProfile: coarse sampling produced no valid points.");
		return 0;
	}

	// ----------------------------------------------------------
	// 2) Find a *stable-branch* bracket where M crosses the target
	//    Stable branch ≡ dM/de_c > 0 → M_{i+1} > M_i
	// ----------------------------------------------------------
	int idx_lo = -1;
	int idx_hi = -1;

	for (int i = 0; i < static_cast<int>(ec_grid.size()) - 1; ++i)
	{
		const double Mi = M_grid[static_cast<std::size_t>(i)];
		const double Mi1 = M_grid[static_cast<std::size_t>(i + 1)];

		// require positive slope (stable segment)
		if (Mi1 <= Mi)
			continue;

		// Check if target lies between Mi and Mi1
		if (Mi <= target_M_solar && target_M_solar <= Mi1)
		{
			idx_lo = i;
			idx_hi = i + 1;
			break;
		}
	}

	// If we didn't find a stable-branch bracket, we fall back to "closest"
	if (idx_lo < 0 || idx_hi < 0)
	{
		Z_LOG_WARNING("SolveToProfile: could not bracket target mass on a stable branch. "
					  "Falling back to closest coarse sample.");

		if (best_idx < 0)
			return 0;

		const double ec_best = ec_grid[static_cast<std::size_t>(best_idx)];

		std::vector<TOVPoint> tmp;
		const int npts = SingleStarSolveToTOVPoints(ec_best, tmp);

		if (npts <= 0 || tmp.empty())
		{
			Z_LOG_ERROR("SolveToProfile: fallback SingleStarSolveToTOVPoints failed.");
			return 0;
		}

		out_tov = std::move(tmp);

		if (out_species_labels)
			*out_species_labels = eos_tab.extra_labels;

		return static_cast<int>(out_tov.size());
	}

	// ----------------------------------------------------------
	// 3) Bisection in e_c on the monotonic (stable) segment
	// ----------------------------------------------------------
	double ec_lo = ec_grid[static_cast<std::size_t>(idx_lo)];
	double ec_hi = ec_grid[static_cast<std::size_t>(idx_hi)];

	double M_lo = M_grid[static_cast<std::size_t>(idx_lo)];
	double M_hi = M_grid[static_cast<std::size_t>(idx_hi)];

	// enforce M_lo < M_hi on the bracket
	if (M_lo > M_hi)
	{
		std::swap(M_lo, M_hi);
		std::swap(ec_lo, ec_hi);
	}

	const double mass_tol = 1e-4; // Msun absolute tolerance
	const int max_iter = 40;

	std::vector<TOVPoint> best_profile;
	double best_M = std::numeric_limits<double>::quiet_NaN();
	best_mass_diff = std::numeric_limits<double>::infinity();

	for (int iter = 0; iter < max_iter; ++iter)
	{
		const double ec_mid = 0.5 * (ec_lo + ec_hi);

		std::vector<TOVPoint> tmp;
		const int npts = SingleStarSolveToTOVPoints(ec_mid, tmp);

		if (npts <= 0 || tmp.empty())
		{
			Z_LOG_ERROR("SolveToProfile: SingleStarSolveToTOVPoints failed at ec_mid = " +
						std::to_string(ec_mid));
			break;
		}

		const double M_mid = tmp.back().m;
		const double diff = std::fabs(M_mid - target_M_solar);

		if (diff < best_mass_diff)
		{
			best_mass_diff = diff;
			best_profile = std::move(tmp);
			best_M = M_mid;
		}

		if (diff < mass_tol)
		{
			// good enough
			break;
		}

		// since we are on a monotonic branch (M increasing with ec),
		// decide which side to keep
		if (M_mid < target_M_solar)
		{
			ec_lo = ec_mid;
			M_lo = M_mid;
		}
		else
		{
			ec_hi = ec_mid;
			M_hi = M_mid;
		}
	}

	if (best_profile.empty())
	{
		Z_LOG_ERROR("SolveToProfile: bisection failed to produce a valid profile.");
		return 0;
	}

	out_tov = std::move(best_profile);

	if (out_species_labels)
		*out_species_labels = eos_tab.extra_labels;

	return static_cast<int>(out_tov.size());
}

//--------------------------------------------------------------
/// Sets the printing precision for the NStar profiles
void TOVSolver::SetProfilePrecision(const int &in_prec)
{
	profile_precision = in_prec;
}

//--------------------------------------------------------------
// Exports the neutron star profile
void TOVSolver::ExportNStarProfile(const size_t &idx,
								   const Zaki::String::Directory &in_dir)
{
	n_star.SetProfilePrecision(profile_precision);
	n_star.Export(in_dir + "_" + std::to_string(idx) + ".tsv");
}

//--------------------------------------------------------------
void TOVSolver::ExportMixedStarProfile(const size_t &v_idx, const size_t &d_idx,
									   const Zaki::String::Directory &in_dir)
{
	mixed_star.Export(in_dir + "_" +
					  std::to_string(d_idx) + "_" +
					  std::to_string(v_idx) + ".tsv");
}

//--------------------------------------------------------------
void TOVSolver::ExportMixedSequence(const Zaki::String::Directory &in_dir)
{
	mixed_sequence.Export(in_dir);
}

//--------------------------------------------------------------
void TOVSolver::SurfaceIsReached(const size_t &v_idx,
								 const size_t &d_idx)
{
	mixed_star.SurfaceIsReached(v_idx, d_idx);
}

//--------------------------------------------------------------
void TOVSolver::SurfaceIsReached()
{
	n_star.FinalizeSurface();
}

//--------------------------------------------------------------
// The radius iteration in the mixed star scenario
void TOVSolver::RadiusLoopMixed(double &in_r, double *in_y,
								double *in_y_mantle)
{
	PROFILE_FUNCTION();

	//----------------------------------------
	//          GSL ODE SYSTEM SETUP
	//----------------------------------------
	gsl_odeiv2_system ode_sys_core =
		{TOVSolver::ODE_Dark_Core, nullptr, 4, this};

	gsl_odeiv2_driver *tmp_driver_core = gsl_odeiv2_driver_alloc_y_new(&ode_sys_core, gsl_odeiv2_step_rk8pd,
																	   1.e-1, 1.e-10, 1.e-10);

	gsl_odeiv2_system ode_sys_mantle =
		{TOVSolver::ODE_Dark_Mantle, nullptr, 2, this};

	gsl_odeiv2_driver *tmp_driver_mantle = gsl_odeiv2_driver_alloc_y_new(&ode_sys_mantle, gsl_odeiv2_step_rk8pd,
																		 1.e-1, 1.e-10, 1.e-10);
	//----------------------------------------

	// double min_log_r = log10(r_min) ;
	// double max_log_r = log10(r_max) ;

	// double min_log_r = log(r_min) ;
	// double max_log_r = log(r_max) ;

	double min_log_r = r_min;
	double max_log_r = r_max;

	double step = (max_log_r - min_log_r) / radial_res;
	double step_scale = 1; // Adaptive steps (Aug 6, 2020)
	// mixed_star.Reserve(radial_res) ;
	bool CORE_REGION = true;

	for (double log_r_i = min_log_r; log_r_i <= max_log_r; log_r_i += step * step_scale)
	{
		// double ri = pow(10, log_r_i) ;
		// double ri = exp(log_r_i) ;
		double ri = log_r_i;

		if (ri < 100) // < 1 m
		{
			step_scale = 0.005;
		}
		else if (ri < 1000) // 1 m - 10 m
		{
			step_scale = 0.025;
		}
		else if (ri < 10000) // 10 m - 100 m
		{
			step_scale = 0.05;
		}
		else if (ri < 100000) // 100 m - 1 km
		{
			step_scale = 0.25;
		}
		else // > 1 km
		{
			step_scale = 1;
		}

		// double tmp_delta_p = 0 ;     // Adaptive steps (Aug 6, 2020)

		int status = -1;

		if (CORE_REGION)
		{
			// tmp_delta_p = in_y[2] ;
			status = gsl_odeiv2_driver_apply(tmp_driver_core, &in_r, ri, in_y);
			// tmp_delta_p = (tmp_delta_p - in_y[2]) / in_y[2] ;
		}
		else
		{
			// tmp_delta_p = in_y_mantle[1] ;

			status = gsl_odeiv2_driver_apply(tmp_driver_mantle, &in_r, ri, in_y_mantle);

			// tmp_delta_p = (tmp_delta_p - in_y_mantle[1]) / in_y_mantle[1] ;
		}

		// 		if (status != GSL_SUCCESS)
		// 		{
		// #if TOV_SOLVER_VERBOSE
		// 			printf("\t-------------------%s-------------------\n", "GSL");
		// 			printf("\t GSL error, return value=%d\n", status);
		// #endif
		// 			if (CORE_REGION)
		// 			{
		// #if TOV_SOLVER_VERBOSE
		// 				printf("\t Visible Pressure = %2.2e.\n", in_y[2]);
		// 				printf("\t Dark Pressure    = %2.2e.\n", in_y[3]);
		// #endif
		// 				if (dark_core)
		// 				{
		// 					m_core = in_y[1];

		// 					// Initial condition for ODE_Mantle
		// 					in_y_mantle[0] = in_y[0];
		// 					in_y_mantle[1] = in_y[2];
		// 				}
		// 				else
		// 				{
		// 					m_core = in_y[0];

		// 					// Initial condition for ODE_Mantle
		// 					in_y_mantle[0] = in_y[1];
		// 					in_y_mantle[1] = in_y[3];
		// #if TOV_SOLVER_VERBOSE
		// 					std::cout << "\n\t m_core = " << in_y[0] / GSL_CONST_CGSM_SOLAR_MASS << "\n";
		// 					std::cout << "\t in_y_mantle[0] = " << in_y_mantle[0] << "\n";
		// 					std::cout << "\t in_y_mantle[1] = " << in_y_mantle[1] << "\n";
		// #endif
		// 				}

		// 				CORE_REGION = false;
		// 				gsl_odeiv2_driver_free(tmp_driver_core);
		// 			}
		// 			else // Mantle's surface reached!
		// 			{
		// #if TOV_SOLVER_VERBOSE
		// 				printf("\t  Surface Pressure = %2.2e.\n", in_y_mantle[1]);
		// 				printf("\t-------------------%s-------------------\n", "GSL");
		// #endif
		// 				break;
		// 			}
		// #if TOV_SOLVER_VERBOSE
		// 			printf("\t-------------------%s-------------------\n", "GSL");
		// #endif

		// 			// Experimental : NOT SURE !!!!!!!!!!!!!!!!!!!
		// 			continue; // Jump over the boundary to avoid duplicate values
		// 		}
		if (status != GSL_SUCCESS)
		{
			// Expected "event" from our RHS: pressure crossed cutoff
			if (status == GSL_EBADFUNC)
			{
				if (CORE_REGION)
				{
					// ODE_Dark_Core sets tov_obj->dark_core to tell us which component ended first
					if (dark_core)
					{
						// Dark pressure hit cutoff first => dark is the core, visible continues as mantle
						m_core = in_y[1];

						in_y_mantle[0] = in_y[0]; // visible mass at boundary radius
						in_y_mantle[1] = in_y[2]; // visible pressure at boundary radius
					}
					else
					{
						// Visible pressure hit cutoff first => visible is the core, dark continues as mantle
						m_core = in_y[0];

						in_y_mantle[0] = in_y[1]; // dark mass at boundary radius
						in_y_mantle[1] = in_y[3]; // dark pressure at boundary radius
					}

					CORE_REGION = false;
					gsl_odeiv2_driver_free(tmp_driver_core);

					// Important: We do NOT treat this like a generic failure; just move on in mantle mode
					continue;
				}
				else
				{
					// Mantle surface reached
					break;
				}
			}

			// Anything else is a real integration failure
#if TOV_SOLVER_VERBOSE
			printf("\t-------------------%s-------------------\n", "GSL");
			printf("\t GSL error, return value=%d\n", status);
			printf("\t-------------------%s-------------------\n", "GSL");
#endif
			break; // or return;
		}

		if (CORE_REGION)
		{
			mixed_star.Append_Core(
				{in_r / 1.e+5, in_y[0] / GSL_CONST_CGSM_SOLAR_MASS,
				 GetNuDer_Dark(in_r, {in_y[0] + in_y[1], in_y[2] + in_y[3]}), 0,
				 in_y[2], GetEDens(in_y[2]),
				 GetRho(in_y[2]), GetRho_i(in_y[2])},
				{in_r / 1.e+5, in_y[1] / GSL_CONST_CGSM_SOLAR_MASS,
				 GetNuDer_Dark(in_r, {in_y[0] + in_y[1], in_y[2] + in_y[3]}), 0,
				 in_y[3], GetEDens_Dark(in_y[3]),
				 GetRho_Dark(in_y[3]), GetRho_i_Dark(in_y[3])});
		}
		else if (dark_core) // dark core with a visible mantle
		{
			mixed_star.Append_Visible_Mantle(
				{in_r / 1.e+5, in_y_mantle[0] / GSL_CONST_CGSM_SOLAR_MASS,
				 GetNuDer_Dark(in_r, {in_y_mantle[0] + m_core, in_y_mantle[1]}), 0,
				 in_y_mantle[1], GetEDens(in_y_mantle[1]),
				 GetRho(in_y_mantle[1]), GetRho_i(in_y_mantle[1])});
		}
		else // visible core, with a dark mantle
		{
			mixed_star.Append_Dark_Mantle(
				{in_r / 1.e+5, in_y_mantle[0] / GSL_CONST_CGSM_SOLAR_MASS,
				 GetNuDer_Dark(in_r, {in_y_mantle[0] + m_core, in_y_mantle[1]}), 0,
				 in_y_mantle[1], GetEDens_Dark(in_y_mantle[1]),
				 GetRho_Dark(in_y_mantle[1]), GetRho_i_Dark(in_y_mantle[1])});
		}
	}

	gsl_odeiv2_driver_free(tmp_driver_mantle);
}

//--------------------------------------------------------------
// Exports the star sequence
void TOVSolver::ExportSequence(const Zaki::String::Directory &in_dir) const
{
	// Zaki::File::VecSaver vec_saver_2(wrk_dir_ + "/" + in_dir);

	// char seq_header[200] ;
	// snprintf(seq_header, sizeof(seq_header), "%-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s",
	//         "ec(g/cm^3)", "M",  "R(km)", "pc(dyne/cm^2)", "B", "I(km^3)" ) ;
	// vec_saver_2.SetHeader(seq_header) ;
	// vec_saver_2.Export1D(sequence) ;

	sequence.Export(in_dir);
}
//--------------------------------------------------------------
// Exports the mixed star sequence
// void TOVSolver::ExportMixedSequence(const Zaki::String::Directory& in_dir) const
// {
//   Zaki::File::VecSaver vec_saver_2(wrk_dir_ + "/" + in_dir);

//   char seq_header[400] ;
//   snprintf(seq_header, sizeof(seq_header), "%-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s"
//                   "\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s\t %-14s",
//           "ec(g/cm^3)", "M",  "R(km)", "pc(dyne/cm^2)", "B", "I(km^3)",
//           "ec_d(g/cm^3)", "M_d",  "R_d(km)", "pc_d(dyne/cm^2)", "B_d", "I_d(km^3)" ) ;

//   vec_saver_2.SetHeader(seq_header) ;

//   vec_saver_2.Export1D(mixed_sequence) ;
// }
//--------------------------------------------------------------
void TOVSolver::SetExclusionRegion(const Zaki::Math::Cond_Polygon &in_c_poly)
{
	c_poly = in_c_poly;
}
//--------------------------------------------------------------
void TOVSolver::PrintStatus(const size_t &in_v_idx,
							const size_t &in_d_idx, const size_t &in_v_res, const size_t &in_d_res)
{
	char tmp_term[150];
	snprintf(tmp_term, sizeof(tmp_term), "Mixed sequence (%3.lu, %3.lu) out "
										 "of (%3.lu, %3.lu).\r",
			 in_v_idx + 1, in_d_idx + 1,
			 in_v_res + 1, in_d_res + 1);
	std::cout << tmp_term << std::flush;
}

//--------------------------------------------------------------
void TOVSolver::PrintStatus(const size_t &in_idx,
							const size_t &in_res)
{
	char tmp_term[100];
	snprintf(tmp_term, sizeof(tmp_term), "Sequence %3.lu out "
										 "of %3.lu.\r",
			 in_idx + 1, in_res + 1);
	std::cout << tmp_term << std::flush;
}
//--------------------------------------------------------------
// Empties the sequence
void TOVSolver::ClearSequence()
{
	sequence.Clear();
}
//--------------------------------------------------------------
class TOVSolver_Tests
{
  protected:
	size_t test_size;

  public:
	size_t Size() const
	{
		return test_size;
	}
	TOVSolver_Tests(const size_t &in_size) : test_size(in_size) {}
	virtual void Modify(TOVSolver *tov_ptr, const size_t &idx) = 0;
};

//--------------------------------------------------------------
class radial_res_test : public TOVSolver_Tests
{
  private:
	const std::vector<size_t> radial_res_set = {10000, 15000, 20000,
												25000, 30000, 40000,
												50000, 55000, 60000,
												65000, 70000, 75000,
												80000, 85000, 90000,
												100000};

  public:
	radial_res_test() : TOVSolver_Tests(16) {}

	virtual void Modify(TOVSolver *tov_ptr, const size_t &idx) override
	{
		tov_ptr->SetRadialRes(radial_res_set[idx]);
	}
};
//--------------------------------------------------------------

//--------------------------------------------------------------
// It generates a sequence of NS by varying radial resolution
// as a function of central energy density
void TOVSolver::GenTestSequence(const double &in_e_c,
								const Zaki::String::Directory &in_dir,
								const Zaki::String::Directory &in_file)
{

//   std::cout << "\n Dir = " << in_dir << "\t file = " << in_file << "\n" ;
// return ;
#if TOV_SOLVER_VERBOSE
	std::cout << "\n\n\t\t ****************************************"
			  << "*******************************" << " \n";
	std::cout << "\t\t *                     "
			  << "TOV Solver Sequence Results"
			  << "                     * \n";
	std::cout << "\t\t ******************************************"
			  << "*****************************" << "\n\n";
#endif
	std::cout << "\n\t TOV_gsl_interp_type = '" << TOV_gsl_interp_type->name << "'\n";
	// size_t test_seq_nums = 10 ;
	// std::vector<size_t> radial_res_set = {10000, 15000, 20000, 25000, 30000} ;
	// std::vector<double> p_prec_set = {1e-4, 1e-5, 1e-6, 1e-7, 1e-8} ;

	radial_res_test test;
	p_of_e_prec = 1e-9;

	for (size_t idx = 0; idx < test.Size(); idx++)
	{
		Z_LOG_WARNING("Sequence " + std::to_string(idx + 1) + " out of " + std::to_string(test.Size()) + ".");

		// SetRadialRes(radial_res_set[idx]) ;
		// p_of_e_prec = p_prec_set[idx] ;

		test.Modify(this, idx);

		if (idx % 1 == 0)
		{
			PrintStatus(idx, test.Size());
		}

		// ------------------------------------------
		//        CANONICAL RADIAL INTEGRATION
		// ------------------------------------------
		// ADR-0005 (ACCEPTED 2026-09-02): SingleStarSolveToTOVPoints is the ONE ordinary-star
		// radial numerical authority. This diagnostic is a resolution-sweep ORCHESTRATOR over
		// it, exactly as Solve(Axis) is a sequence orchestrator over it. Phase 3E-I4 replaced
		// the former inline setup + RadiusLoop call; RadiusLoop itself is gone.
		//
		// Bit-identical on the ordinary valid domain: measured 16/16 bitwise on
		// (ec, M, R, pc, B, I) at every one of the 16 resolutions before the migration
		// (docs/validation/PHASE3E_I4_RADIUSLOOP_RETIREMENT.md).
		//
		// DELIBERATE BEHAVIOR CHANGE, out-of-domain only: the former code called
		// p_of_e(in_e_c) directly, bypassing the central-density floor/ceiling clamp that
		// Solve() and SolveToProfile() both apply. Measured, that path did not degrade
		// gracefully -- a request below the EOS energy-density minimum ABORTED the process
		// (SIGABRT). The canonical primitive clamps into the valid band instead, so such a
		// request now yields the clamped star rather than a crash. Compatibility is guaranteed
		// on the ordinary valid domain; the out-of-domain change is an improvement, recorded
		// rather than claimed as preservation.
		//
		// p_of_e_prec = 1e-9 (set above, and as before never restored) still applies: the
		// primitive calls p_of_e.
		std::vector<TOVPoint> tov_points;
		SingleStarSolveToTOVPoints(in_e_c, tov_points);

		for (const auto &tp : tov_points)
			n_star.Append(tp);
		// ------------------------------------------
		SurfaceIsReached();
		// ------------------------------------------

		// ----------------------------------------------
		if (analysis)
			analysis->Analyze(&n_star);
		// ----------------------------------------------

		// ----------------------------------------------
		// Saving the results
		// ----------------------------------------------
		if (n_exp_cond_f)
		{
			if (n_exp_cond_f(n_star))
			{
				ExportNStarProfile(idx, in_dir + "/profiles" + in_file);
			}
		}
		sequence.Add(n_star);
		// ----------------------------------------------
		n_star.Reset();
	}
	// ----------------------------------------------------------------
	//                  TOV Visible sequence loop ends!
	// ----------------------------------------------------------------
#if TOV_SOLVER_VERBOSE
	std::cout << "\n\t\t *************************"
			  << " TOV Solver Finished *************************" << "\n\n";
#endif
	ExportSequence(in_dir + in_file + "_TestSequence.tsv");

	// ----------------------------------------------
	if (analysis)
		analysis->Export(wrk_dir_ + in_dir);
	// ----------------------------------------------
}

//==============================================================
