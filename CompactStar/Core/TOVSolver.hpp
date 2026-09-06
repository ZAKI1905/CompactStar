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
 * @file TOVSolver.hpp
 *
 * @brief Solves TOV equations.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@eku.edu
 *
 */
// Last edit on Dec 15 2020
#ifndef CompactStar_Core_TOVSolver_H
#define CompactStar_Core_TOVSolver_H

#include <gsl/gsl_spline.h>
#include <limits>
#include <vector>

#include <Zaki/Math/Math_Core.hpp>

#include "CompactStar/Core/MixedStar.hpp"
#include "CompactStar/Core/Prog.hpp"
#include "CompactStar/Core/StarProfile.hpp"
//==============================================================
namespace CompactStar::Core
{

class Analysis;
//==============================================================
//                      eps_pair Class
//==============================================================
/**
 * @class eps_pair
 * @brief Class representing a pair of energy density and pressure
 * values for the dark and visible components of a mixed star.
 *
 */
class eps_pair : public Zaki::Math::Coord2D
{
  public:
	/// Constructor from two double values
	/// @param in_e_v Energy density of the visible component (x-coordinate)
	/// @param in_e_d Energy density of the dark component (y-coordinate)
	/// @details The constructor initializes the x and y coordinates of the
	/// eps_pair object with the provided energy density values.
	eps_pair(const double &in_e_v, const double &in_e_d)
		: Coord2D(in_e_v, in_e_d) {}

	/// Constructor from a Coord2D object
	/// @param in_c Coord2D object containing the energy density values
	/// @details The constructor initializes the x and y coordinates of the
	/// eps_pair object with the x and y coordinates of the provided Coord2D object.
	eps_pair(const Coord2D &in_c)
		: Coord2D(in_c.x, in_c.y) {}

	/// @brief Getter for the visble energy density (x-coordinate)
	/// @return The visible energy density value
	double e_v() const
	{
		return x;
	}

	/// @brief Getter for the dark energy density (y-coordinate)
	/// @return The dark energy density value
	double e_d() const
	{
		return y;
	}
};
//==============================================================

//==============================================================
//                      Contour Class
//==============================================================
/**
 * @class Contour
 * @brief Class representing a contour in the TOV solution space.
 *
 */

class Contour
{
  public:
	Zaki::Math::Curve2D curve;
	double val = 0;
	double precision = 1e-8;
	size_t max_steps = 35;

	/// @brief Default constructor
	/// @details Initializes an empty contour object.
	/// @note The curve is empty and the value is set to 0.
	Contour() {}

	/// @brief Constructor with label
	/// @param in_label Label for the contour
	/// @details Initializes a contour object with the specified label.
	/// @note The curve is empty and the value is set to 0.
	Contour(const std::string &in_label) : curve(in_label) {}

	/// @brief Getter for the size of the curve
	/// @details Returns the number of points in the curve.
	/// @return The size of the curve
	size_t Size() const
	{
		return curve.Size();
	}

	/// @brief Imports the contour from a file
	/// @param in_file
	void Import(const Zaki::String::Directory &in_file)
	{
		curve.Import(in_file);
		val = std::atof(curve.GetLabel().c_str());
	}
};

//==============================================================
//                      EOSTable Struct
//==============================================================
/**
 * @struct EOSTable
 * @brief Class representing a table of equation of state (EOS) data.
 *
 * @details The EOSTable class stores the energy density, pressure, and baryon number density
 * values for a given EOS. It also provides methods to set labels, add extra labels, and print the table.
 */
struct EOSTable
{
  private:
	std::string eps_label;
	std::string pre_label;
	std::string rho_label;

  public:
	std::vector<double> eps;
	std::vector<double> pre;
	std::vector<double> rho;

	std::vector<std::vector<double>> rho_i;
	std::vector<std::string> extra_labels;

	/// @brief Getter for the size of the table
	/// @details Returns the number of rows in the table.
	/// @return The size of the table
	size_t Size()
	{
		return eps.size();
	}

	/// @brief Setter for the labels
	/// @details Sets the labels for the energy density, pressure, and baryon number density.
	/// @param in_eps_label Label for the energy density
	/// @param in_pre_label Label for the pressure
	/// @param in_rho_label Label for the baryon number density
	/// @note The labels are stripped of leading and trailing spaces.
	/// @note The labels are stored as strings.
	void SetLabels(const std::string &in_eps_label,
				   const std::string &in_pre_label,
				   const std::string &in_rho_label)
	{
		eps_label = Zaki::String::Strip(in_eps_label, ' ');
		pre_label = Zaki::String::Strip(in_pre_label, ' ');
		rho_label = Zaki::String::Strip(in_rho_label, ' ');
	}

	/// @brief Setter for the extra labels
	/// @param in_label Label for the extra data
	void AddExtraLabels(const std::string &in_label)
	{
		extra_labels.emplace_back(in_label);
	}

	/// @brief Printer for the table
	/// @details Prints the table to the standard output.
	// void Print() const
	// {
	// 	std::cout << " *-------------------------------------------* " << "\n";
	// 	std::cout << " | " << eps_label << "   | " << pre_label
	// 			  << ")   | " << rho_label << "  |\n";
	// 	std::cout << " *-------------------------------------------* " << "\n";
	// 	std::string tmp_str;
	// 	for (size_t i = 0; i < eps.size(); i++)
	// 	{
	// 		std::cout << " | " << eps[i] << "\t      | "
	// 				  << pre[i] << "\t      | " << rho[i] << "\n";
	// 	}
	// 	std::cout << " *--------------------------------------------* " << "\n";
	// }

	/// @brief Prints the EOS table to the standard output.
	/// @details
	///  - By default, prints only the header and the first few lines (5).
	///  - If `n_lines` is specified and positive, prints up to that many rows.
	///  - If `n_lines` exceeds the table size, prints the entire table.
	///  - Negative values for `n_lines` (e.g. -1) force printing *all* rows.
	/// @param n_lines Optional number of data rows to print (default = 5).
	void Print(const int n_lines = 5) const
	{
		// Guard: no data
		if (eps.empty() || pre.empty() || rho.empty())
		{
			std::cout << "\n[EOS::Print] Warning: table is empty — nothing to print.\n";
			return;
		}

		std::cout << "\n *-------------------------------------------* \n";
		std::cout << " | " << eps_label << "   | " << pre_label
				  << "   | " << rho_label << "  |\n";
		std::cout << " *-------------------------------------------* \n";

		// Determine number of rows to print
		size_t total = eps.size();
		size_t limit = total;

		if (n_lines >= 0)
		{
			limit = std::min<size_t>(n_lines, total);
		} // if n_lines < 0, print all rows

		// Loop over selected rows
		for (size_t i = 0; i < limit; i++)
		{
			std::cout << " | "
					  << std::scientific << std::setprecision(6)
					  << eps[i] << "   | "
					  << pre[i] << "   | "
					  << rho[i] << "\n";
		}

		// If we didn’t print all lines, indicate how many were omitted
		if (limit < total)
		{
			std::cout << " | ... (" << (total - limit)
					  << " more rows omitted) ...\n";
		}

		std::cout << " *-------------------------------------------* \n";
	}

	/// @brief Prints a compact summary of the EOS table.
	/// @details
	///  Shows:
	///   - Number of data rows.
	///   - Column labels.
	///   - Minimum and maximum values of ε, P, and ρ.
	///   - Extra species columns, if present.
	///  Intended for quick sanity checks before running the TOV solver.
	void PrintSummary() const
	{
		std::cout << "\n================ EOS Table Summary ================\n";

		if (eps.empty() || pre.empty() || rho.empty())
		{
			std::cout << "[EOS::PrintSummary] Table is empty.\n";
			std::cout << "===================================================\n";
			return;
		}

		const size_t n = eps.size();

		auto [eps_min, eps_max] = std::minmax_element(eps.begin(), eps.end());
		auto [pre_min, pre_max] = std::minmax_element(pre.begin(), pre.end());
		auto [rho_min, rho_max] = std::minmax_element(rho.begin(), rho.end());

		std::cout << std::scientific << std::setprecision(6);
		std::cout << "  # of rows: " << n << "\n";
		std::cout << "  Columns: " << eps_label << " | " << pre_label << " | " << rho_label << "\n";
		std::cout << "  ε (energy density):  [" << *eps_min << ", " << *eps_max << "]\n";
		std::cout << "  P (pressure):         [" << *pre_min << ", " << *pre_max << "]\n";
		std::cout << "  ρ (baryon density):   [" << *rho_min << ", " << *rho_max << "]\n";

		if (!extra_labels.empty())
		{
			std::cout << "  Additional species: ";
			for (size_t i = 0; i < extra_labels.size(); ++i)
			{
				std::cout << extra_labels[i];
				if (i + 1 < extra_labels.size())
					std::cout << ", ";
			}
			std::cout << "\n";
		}

		std::cout << "===================================================\n";
	}
};
//==============================================================
//            nu_der added on December 15, 2020
// Struct representing TOV solution points
struct TOVPoint
{
	// The canonical integrator holds raw mass in grams internally. Publication:
	// r [km], m = m_grams/GSL_CONST_CGSM_SOLAR_MASS (literal public ratio),
	// nu_der [cm^-1], p [dyn/cm^2], e [g/cm^3], rho [fm^-3].
	// Ordinary NStar maps m,p,e,nu_der to one GSL spacetime via RelativityUnits
	// (ADR-0012); m here is neither geometric km nor an IAU nominal GM ratio.
	double r, m, nu_der, nu, p, e, rho;
	std::vector<double> rho_i;

	/**
	 * @brief The barotropic EOS derivative d(eps)/dp at this node — **DIMENSIONLESS**.
	 *
	 * ADR-0007 P5 (ACCEPTED 2026-09-02) makes the EOS the sole authority for this quantity:
	 * it is the derivative of the *same* `eps(p)` interpolant that constructed the star
	 * (`TOVSolver::GetEDensDeriv`), never a finite difference of the radial profile.
	 * `d(eps)/dp = 1/c_s^2` in `c = 1` units, so causal matter has `d(eps)/dp >= 1` and
	 * incompressible matter has `d(eps)/dp = 0` as a formal limit.
	 *
	 * Unlike `p` and `e`, which are carried across this boundary in cgs, this member is
	 * dimensionless and therefore unit-system independent.
	 *
	 * **NaN means "not supplied".** That is the default, and it is deliberately distinct from
	 * an explicitly supplied `0.0` (the correct value for incompressible matter). A profile
	 * built from points that do not all carry a finite value simply has no EOS-derivative
	 * data, and any consumer that requires one fails closed — see
	 * `StarProfile::HasEosDEdP()`.
	 */
	double dedp;

	TOVPoint(const double &in_r, const double &in_m,
			 const double &in_nu_der, const double &in_nu,
			 const double &in_p, const double &in_e,
			 const double &in_rho, const std::vector<double> &in_rho_i,
			 const double &in_dedp = std::numeric_limits<double>::quiet_NaN())
		: r(in_r), m(in_m), nu_der(in_nu_der), nu(in_nu),
		  p(in_p), e(in_e), rho(in_rho),
		  rho_i(in_rho_i), dedp(in_dedp)
	{
	}

	std::string Str() const
	{
		std::stringstream ss;
		char tmp[200];
		snprintf(tmp, sizeof(tmp), "%.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e\t %.8e",
				 r, m, nu_der, nu, p, e, rho);
		ss << tmp;
		for (double tmp_rho : rho_i)
		{
			snprintf(tmp, sizeof(tmp), "\t %.8e", tmp_rho);
			ss << tmp;
		}

		return ss.str();
	}
};
//==============================================================
//            nu_der added on December 15, 2020
// Struct representing TOV solution points
struct TOV_Nu_Point
{
	double r, m, nu;

	TOV_Nu_Point(const double &in_r, const double &in_m, const double &in_nu)
		: r(in_r), m(in_m), nu(in_nu)
	{
	}

	std::string Str() const
	{
		std::stringstream ss;
		char tmp[150];
		snprintf(tmp, sizeof(tmp), "%.8e\t %.8e\t %.8e", r, m, nu);
		ss << tmp;

		return ss.str();
	}
};

//==============================================================
//                        Sequence class
//==============================================================
/**
 * @class Sequence
 * @brief Class representing a sequence of TOV solution points.
 *
 * @details The Sequence class stores a vector of SeqPoint objects,
 * which represent individual TOV solution points. It provides methods
 * to add points, export the sequence to a file, combine sequences,
 * and clear the sequence.
 */

class Sequence : public Prog
{
  private:
	/// @brief Vector of SeqPoint objects representing the TOV solution points
	/// @details The vector stores the sequence of TOV solution points.
	/// Each SeqPoint object contains the energy density, mass, radius,
	/// central pressure, baryon number, and moment of inertia.
	std::vector<SeqPoint> seq;

  public:
	/// @brief Default constructor
	/// @details Initializes an empty sequence object.
	Sequence();

	/// @brief Destructor
	~Sequence()
	{
		// std::cout << "[ Thread = " << std::this_thread::get_id()
		//           << "] Destructor called for: "
		//           << name << "\n" ;
	}

	/// @brief Adds a neutron star point to the sequence
	/// @param in_star the neutron star point to be added
	void Add(const NStar &in_star);

	/// @brief Exports the star sequence to a file
	/// @param in_dir the directory to export the sequence to
	void Export(const Zaki::String::Directory &in_dir = "") const;

	/// @brief Combines two sequences
	/// @param other the other sequence to combine with
	void Combine(const Sequence &other);

	/// @brief Clears the sequence
	void Clear();
};

//==============================================================
//                        MixedSequence class
//==============================================================
class MixedSequence : public Prog
{
  private:
	std::vector<MixedSeqPoint> seq;
	// std::mutex m_mutex ;

  public:
	MixedSequence();
	// MixedSequence(const std::string& in_name )
	// : Prog(in_name) {}

	~MixedSequence()
	{
		// std::cout << "[ Thread = " << std::this_thread::get_id()
		//           << "] Destructor called for: "
		//           << name << "\n" ;
	}

	void Add(const MixedStar &in_star);

	// Exports the mixed star sequence
	void Export(const Zaki::String::Directory &in_dir = "") const;

	// Combines two sequences
	void Combine(const MixedSequence &other);

	// Clears the sequence
	void Clear();
};

//==============================================================
/// Completion of the ordinary-star solve (ADR-0009).
enum class TOVSolveStatus
{
	SURFACE_REACHED,
	GSL_FAILURE,
	EOS_DOMAIN_FAILURE,
	R_MAX_EXHAUSTED,
	INVALID_INITIAL_STATE,
	PARTITION_INVALID,
	MASS_TARGET_NOT_REACHED
};

class TOVSolver : public Prog
{
	friend class MixedStar;
	friend class NStar;

  private:
	TOVSolveStatus last_solve_status = TOVSolveStatus::INVALID_INITIAL_STATE;
	static int PressureCoordinateODE(double p, const double y[], double f[], void *params);
	bool LocateSurface(double p_above, double r_above, double m_above,
					   double r_below, double &r_surface, double &m_surface);
	//--------------------------------------------------------------
  protected:
	/// The number of TOV object instances
	// static inline std::atomic<size_t> tov_counter = 0;

	// std::mutex m_mutex ;

	/// The radial resolution for the solver
	size_t radial_res = 10000;

	EOSTable eos_tab;
	EOSTable eos_tab_dark;

	MixedStar mixed_star;
	NStar n_star;

	// This workspace stores state variables for interpolation lookups.
	// It caches the previous value of an index lookup.
	// When the subsequent interpolation point falls in the same
	// interval its index value can be returned immediately.
	// There has to be separate accelerators for variables
	// with different domains.
	//
	// EOS splines ( Domain = pressure )
	gsl_interp_accel *visi_p_accel = nullptr;
	gsl_interp_accel *dark_p_accel = nullptr;
	gsl_interp_accel *mixed_r_accel = nullptr;

	const gsl_interp_type *TOV_gsl_interp_type = gsl_interp_steffen;
	// const gsl_interp_type* TOV_gsl_interp_type = gsl_interp_linear ;

	// gsl_interp_accel *accel = nullptr ;
	// gsl_interp_accel *dark_accel = nullptr ;

	// We use cubic spline
	// Cubic spline with natural boundary conditions.
	// The resulting curve is piecewise cubic on each interval,
	// with matching first and second derivatives at the supplied data-points.
	//  The second derivative is chosen to be zero at the first point and last point.
	gsl_spline *visi_eps_p_spline = nullptr;
	gsl_spline *dark_eps_p_spline = nullptr;

	/// Total baryon number density (function of p)
	gsl_spline *visi_rho_p_spline = nullptr;
	gsl_spline *dark_rho_p_spline = nullptr;

	/// Total baryon number density (function of r)
	// gsl_spline *visi_rho_r_spline = nullptr ;
	// gsl_spline *dark_rho_r_spline = nullptr ;

	/// This is for the extra rho's (individual species)
	std::vector<gsl_spline *> visi_rho_i_p_spline;
	std::vector<gsl_spline *> dark_rho_i_p_spline;

	// Added on December 15, 2020
	/// nu(r) derivative spline
	gsl_spline *nu_der_r_spline = nullptr;

	/// Initial pressure (at r = r(i=0))
	double init_press = -1;
	double init_press_dark = -1;

	/// Initial energy density (at r = r(i=0))
	double init_edens = -1;
	double init_edens_dark = -1;

	/// Solution to TOV are saved here
	// std::vector<TOVPoint> results ;
	// std::vector<TOVPoint> results_dark ;

	// std::vector<SeqPoint> sequence ;
	Sequence sequence;
	// std::vector<MixedSeqPoint> mixed_sequence ;
	MixedSequence mixed_sequence;

	// static inline MixedSequence mixed_seq_static ;

	/// Returns the pressure corresponding to in_e
	/// It's the inverse function of "GetEDens"
	double p_of_e(const double &in_e);
	double p_of_e_dark(const double &in_e);

	/// Cost function for finding the pressure
	/// for "p_of_e" method.
	double cost_p_of_e(const double in_e);
	double cost_p_of_e_dark(const double in_e);
	double cost_p_of_e_input = 0;
	double cost_p_of_e_input_dark = 0;

	/// The value of pressure cut-off is the pressure
	/// at the surface of the star
	/// Theoretically it's zero, but we choose
	/// 1e-5 times smaller than the central pressure
	/// or the lowest possible pressure given by the Eq.
	/// of state.
	double PressureCutoff() const;
	double PressureCutoff_Dark() const;

	/// Minimum radius in the solver (cm)
	double r_min = 1.;

	/// Maximum radius in the solver (cm)
	// double r_max = 20e+5 ; // 20 km
	// double r_max = 50e+5 ; // 50 km
	double r_max = 70e+5; // 70 km

	// dark core or visible core :
	bool dark_core = true;

	double m_core = -1;
	enum class TOVSolverStatus
	{
		Vis_Surf_Reached = 100,
		Dark_Surf_Reached = 101,
		Mantle_Surf_Reached = 102
	};

	/**
	 * @brief Propagate the new working directory to internal members.
	 *
	 * Called automatically by SetWrkDir() after this object's own
	 * working directory is updated and created on disk.
	 *
	 * Derived classes override this method to push the directory
	 * into member objects that also produce output files.
	 *
	 * @param dir Newly assigned working directory.
	 */
	void OnWorkDirChanged(const Zaki::String::Directory &dir) override;

	/// Pointer to the condition for exporting
	///  the mixed star profile
	bool (*mix_exp_cond_f)(const MixedStar &) = nullptr;

	/// Pointer to the condition for exporting
	/// the star profile
	bool (*n_exp_cond_f)(const NStar &) = nullptr;

	/// Pointer to the analysis object
	Analysis *analysis = nullptr;

	/// Exclusion region in (ec_v, ec_d) space
	Zaki::Math::Cond_Polygon c_poly;

	/// Number of points that are ignored
	///  because they are in the exclusion region.
	size_t ignored_counter = 0;

	void Hidden_ImportEOS_Vis(const Zaki::String::Directory &eos_file,
							  const bool absolute_path = false);
	void Hidden_ImportEOS_Dar(const Zaki::String::Directory &eos_file,
							  const bool absolute_path = false);

	// The precision in printing the profiles
	int profile_precision = 9;

	/// The precision in evaluation pressure as a function of density
	double p_of_e_prec = 1e-4;

	/**
	 * @brief Multiplier applied to the EOS minimum energy density
	 *        to build a lower bound for central energy density.
	 *
	 * @details If the user provides an axis with values below
	 *          @c central_eps_floor_factor * eos_tab.eps.front(),
	 *          the solver will clamp to that floor.
	 *
	 * Default: 10.0.
	 */
	double central_eps_floor_factor = 10.0;
	//--------------------------------------------------------------
  public:
	/**
	 * @brief Construct a new TOVSolver object
	 *
	 */
	TOVSolver();

	/**
	 * @brief Destroy the TOVSolver object
	 *
	 */
	~TOVSolver();

	/**
	 * @brief Deleted copy constructor to prevent copying of TOVSolver objects.
	 */
	TOVSolver(const TOVSolver &) = delete;

	/**
	 * @brief Deleted assignment operator to prevent copying of TOVSolver objects.
	 */
	TOVSolver &operator=(const TOVSolver &) = delete;

	/**
	 * @brief Import the equation of state (EOS) from a file.
	 * @param eos_file Path to the EOS file.
	 * @param absolute_path If true, the provided path is treated as an absolute path.
	 */
	void ImportEOS(const Zaki::String::Directory &eos_file, const bool absolute_path = false);

	/**
	 * @brief Import the visible and dark equation of state (EOS) from files.
	 * @param vis_eos Path to the visible EOS file.
	 * @param dar_eos Path to the dark EOS file.
	 * @param absolute_path If true, the provided paths are treated as absolute paths.
	 */
	void ImportEOS(const Zaki::String::Directory &vis_eos,
				   const Zaki::String::Directory &dar_eos,
				   const bool absolute_path = false);

	// ----------------------------------------------------------
	// EOS inspection / debugging
	// ----------------------------------------------------------
	/**
	 * @brief Print (at most) @p max_rows rows of the EOS table.
	 *
	 * @param max_rows Maximum number of data rows (after the header)
	 *                 to print. If @p max_rows == 0 the function
	 *                 only prints the header and a short footer.
	 *
	 * @details This is useful for large CompOSE tables: we can
	 * quickly verify that the file was read correctly, that labels
	 * were parsed, and that the first few rows are monotonic in
	 * pressure, without flooding the terminal.
	 */
	void PrintEOSTable(const std::size_t max_rows) const;

	/**
	 * @brief Print a compact EOS summary.
	 *
	 * @details Shows number of rows, column labels, and min/max for
	 * energy density, pressure, and baryon density. This is what we
	 * typically want right after ImportEOS(...).
	 */
	void PrintEOSSummary() const;

	/**
	 * @brief Return the minimum energy density available in the EOS.
	 *
	 * @return Minimum energy density (first row) or 0.0 if the EOS
	 *         table is empty.
	 */
	double GetEOSMinEDens() const;

	/**
	 * @brief Return the maximum energy density available in the EOS.
	 *
	 * @return Maximum energy density (last row) or 0.0 if the EOS
	 *         table is empty.
	 */
	double GetEOSMaxEDens() const;

	/**
	 * @brief Set the floor multiplier for central energy density.
	 *
	 * @param f Factor by which the minimum EOS energy density is
	 *          multiplied to form the *lowest* allowed central
	 *          energy density in Solve(...).
	 *
	 * @details If the EOS minimum energy density is
	 * \f$\varepsilon_{\min}\f$, then the solver will clamp user
	 * inputs to at least
	 * \f$ \varepsilon_{\text{floor}} = f \times \varepsilon_{\min} \f$.
	 *
	 * This is to avoid trying to build stars whose *central* energy
	 * density is already too close to the surface value of the EOS,
	 * which typically makes the TOV integration stop immediately.
	 *
	 * @note A typical value is 10.0.
	 */
	void SetCentralEDensFloorFactor(double f);
	// ----------------------------------------------------------

	/**
	 * @brief Get the energy density corresponding to a given pressure.
	 * @param in_pressure value for which to find the energy density.
	 * @return Energy density corresponding to the given pressure.
	 */
	double GetEDens(const double &in_pressure);

	// ----------------------------------------------------------
	//   EOS thermodynamic derivative — the single scientific authority
	//   (ADR-0007 P5, ACCEPTED 2026-09-02; Phase 4C-I0)
	// ----------------------------------------------------------
	/**
	 * @brief Whether the visible `eps(p)` interpolant exists and can supply a derivative.
	 *
	 * False before `ImportEOS`, or when the table was too short / non-monotonic for GSL.
	 */
	bool HasEDensDeriv() const noexcept;

	/// Lower end of the `eps(p)` interpolant's pressure domain, in the EOS table's own
	/// units (dyne/cm^2). NaN when no interpolant exists.
	double EDensDerivPressMin() const noexcept;

	/// Upper end of the `eps(p)` interpolant's pressure domain, in the EOS table's own
	/// units (dyne/cm^2). NaN when no interpolant exists.
	double EDensDerivPressMax() const noexcept;

	/**
	 * @brief The barotropic EOS derivative `d(eps)/dp` at a given pressure — **DIMENSIONLESS**.
	 *
	 * This is the quantity Hartle's O(Omega^2) `l = 0` system needs (H67 eq. 97, `dE/dP`), and
	 * ADR-0007 P5 fixes its authority: it is the analytic derivative of the **same**
	 * `visi_eps_p_spline` from which `GetEDens` reads `eps`, evaluated through the same
	 * accelerator. There is no second interpolant, no profile finite difference, and no
	 * fallback value.
	 *
	 * **Units.** The table stores `eps` in g/cm^3 against `p` in dyne/cm^2, so the raw spline
	 * derivative carries s^2/cm^2. The profile conversion in `NStar::BuildFromTOV` scales
	 * `eps` by `INV_FM4_2_INV_KM2 / INV_FM4_2_G_CM3` and `p` by
	 * `INV_FM4_2_INV_KM2 / INV_FM4_2_Dyn_CM2`, so the geometric (and hence dimensionless)
	 * derivative is the raw one times `INV_FM4_2_Dyn_CM2 / INV_FM4_2_G_CM3` — which is `c^2` in
	 * cgs. That factor is applied here, and this function is its only owner.
	 *
	 * **Domain semantics — stricter than a value lookup, deliberately.** A derivative is a
	 * *local* property of the state, so a clamped out-of-range query would silently describe a
	 * different state than the one asked about. `GetEDens` clamps and warns; this function
	 * **fails closed**: outside `[EDensDerivPressMin(), EDensDerivPressMax()]`, on a non-finite
	 * argument, with no interpolant, or on any GSL error it logs and returns NaN. Exact
	 * endpoints are evaluated, not rejected. Callers must test `std::isfinite`.
	 *
	 * @param in_pressure pressure in the EOS table's own units (dyne/cm^2).
	 * @return `d(eps)/dp`, dimensionless, or NaN if it cannot be evaluated for this argument.
	 */
	double GetEDensDeriv(const double &in_pressure);

	/**
	 * @brief Get the energy density of the dark component corresponding to a given pressure.
	 * @param in_pressure value for which to find the dark energy density.
	 *
	 * @return Dark energy density corresponding to the given pressure.
	 */
	double GetEDens_Dark(const double &in_pressure);

	/**
	 * @brief Get the total baryon number density corresponding to a given pressure.
	 * @param in_pressure value for which to find the baryon number density.
	 * @return Total baryon number density corresponding to the given pressure.
	 */
	double GetRho(const double &in_pressure);

	/**
	 * @brief Get the total baryon number density of the dark component corresponding to a given pressure.
	 * @param in_pressure value for which to find the dark baryon number density.
	 *
	 * @return Dark total baryon number density corresponding to the given pressure.
	 */
	double GetRho_Dark(const double &in_pressure);

	//........... NOV 3, 2021 Begins ...........
	/// Returns the total baryon number density given radius
	double GetRho_r(const double &);
	//........... NOV 3, 2021 Ends  ............

	/// Returns the specific number density given pressure
	std::vector<double> GetRho_i(const double &in_p);
	std::vector<double> GetRho_i_Dark(const double &in_p);

	double GetInitPress() const;
	double GetInitEDens() const;
	double GetInitPress_Dark() const;

	//               Added on December 15, 2020
	/// Returns the derivative of the metric nu(r) function
	double GetNuDer(const double r, const std::vector<double> &y);
	double GetNuDer_Dark(const double r, const std::vector<double> &y);

	//            Added on December 15, 2020
	// Returns the nu_der value given the radius input
	double GetNuDerSpline(const double &in_r);

	/// Adds the condition for printing the mixed star profile
	void AddMixCondition(bool (*func)(const MixedStar &));

	/// Adds the condition for printing the star profile
	void AddNCondition(bool (*func)(const NStar &));

	/// Attaches a pointer to analysis
	void AddAnalysis(Analysis *);

	/**
	 * @brief Solve TOV equations over a range of central energy densities.
	 * @param in_ax Axis defining the range of central energy densities.
	 * @param dir Directory to export the results to.
	 * @param file_name File name for the results.
	 */
	void Solve(const Zaki::Math::Axis &,
			   const Zaki::String::Directory &,
			   const Zaki::String::Directory &file_name);

	void Solve_Mixed(const Zaki::Math::Axis &vis_ax,
					 const Zaki::Math::Axis &dark_ax,
					 const Zaki::String::Directory &dir,
					 const Zaki::String::Directory &file_name);

	void Solve_Mixed(const Contour &eps_cont,
					 const Zaki::String::Directory &dir,
					 const Zaki::String::Directory &file_name);

	// The radius iteration in the mixed star scenario
	void RadiusLoopMixed(double &r, double *y_core,
						 double *y_mantle);

	/**
	 * @brief Export the generated sequence of neutron stars.
	 * @param dir Directory to export the sequence to.
	 *
	 */
	void ExportSequence(const Zaki::String::Directory &) const;

	/// Exports the mixed sequence
	virtual void ExportMixedSequence(const Zaki::String::Directory &);

	/// Exports the mixed star profile
	virtual void ExportMixedStarProfile(const size_t &v_idx,
										const size_t &d_idx, const Zaki::String::Directory &);

	/// Exports the neutron star profile
	virtual void ExportNStarProfile(const size_t &idx,
									const Zaki::String::Directory &);

	virtual void SurfaceIsReached(const size_t &v_idx,
								  const size_t &d_idx);

	/**
	 * @brief Handle the event when the surface is reached.
	 *
	 */
	void SurfaceIsReached();

	// For mixed stars
	virtual void PrintStatus(const size_t &v_idx,
							 const size_t &d_idx,
							 const size_t &v_res,
							 const size_t &d_res);
	// For neutron stars
	void PrintStatus(const size_t &in_idx, const size_t &in_res);

	//            Added on December 15, 2020
	/// Evaluates & exports the nu(r) function table
	// void ExportNu(const Zaki::String::Directory& in_dir) ;

	static int ODE(double r, const double y[], double f[], void *params);
	static int ODE_Dark_Core(double r, const double y[], double f[], void *params);
	static int ODE_Dark_Mantle(double r, const double y[], double f[], void *params);

	/**
	 * @brief Set the Exclusion Region object
	 * @param polygon The polygon defining the exclusion region in (ec_v, ec_d) space
	 * @details This method sets the exclusion region for the TOV solver.
	 *
	 */
	void SetExclusionRegion(const Zaki::Math::Cond_Polygon &);

	/**
	 * @brief Set the Radial Res object
	 * @param in_res The radial resolution for the solver
	 * @details This method sets the radial resolution for the TOV solver.
	 *
	 */
	void SetRadialRes(const size_t &);

	/**
	 * @brief Set the Profile Precision object
	 *
	 * @param prec The printing precision for the star profiles
	 * @details This method sets the printing precision for the star profiles.
	 */
	void SetProfilePrecision(const int &prec);

	/// Sets maximum value for radius in the solver
	/// if a maximum is known, say 15 km for NS,
	/// it would increase the radial resolution
	/// defined by:
	///  delta R = scale(R) * (R_max - R_min) / radial_res
	/**
	 * @brief Set the Max Radius object
	 * @param The maximum radius for the solver
	 * @details This method sets the maximum radius for the TOV solver.
	 */
	void SetMaxRadius(const double &);

	/**
	 * @brief Clear the generated sequence of neutron stars.
	 *
	 */
	void ClearSequence();

	/** @brief It generates a sequence of NS by varying radial resolution
	 * @param e_c central energy density
	 * @param dir directory for saving the results
	 * @param file result file name
	 */
	void GenTestSequence(const double &e_c,
						 const Zaki::String::Directory &dir,
						 const Zaki::String::Directory &file);

	/**
	 * @brief Solve TOV for a single star with given central energy density
	 *        and return the radial TOV points.
	 *
	 * Contract:
	 *  - @p ec_central is in the same units as eos_tab.eps (g/cm^3).
	 *  - Only SURFACE_REACHED returns >0 and fills @p out_tov with (r, m, nu', nu=0, p, e, rho, rho_i).
	 *  - On failure, returns 0 and leaves @p out_tov empty.
	 *
	 * Output r is km, m is solar masses, pressure and energy density are cgs,
	 * and nu_prime is cm^-1. Baryon-number / I integration is left to NStar::BuildFromTOV.
	 */
	TOVSolveStatus LastSolveStatus() const noexcept { return last_solve_status; }

	int SingleStarSolveToTOVPoints(double ec_central,
								   std::vector<TOVPoint> &out_tov);

	/**
	 * @brief Solve the TOV equations for a *single neutron star* specified by
	 *        a target gravitational mass, returning the full radial structure
	 *        as a vector of @ref TOVPoint.
	 *
	 * This routine performs an internal search over central energy density
	 * \( \varepsilon_c \) to find the value that yields the requested
	 * gravitational mass (in solar masses). The EOS must have been imported
	 * beforehand via @ref ImportEOS.
	 *
	 * ### Algorithm summary
	 * 1. Construct a coarse, logarithmically spaced grid in \( \varepsilon_c \)
	 *    between the allowed EOS range (with the usual floor/ceiling safety margins).
	 * 2. For each sampled \( \varepsilon_c \):
	 *      - Integrate the TOV equations using
	 *        @ref SingleStarSolveToTOVPoints,
	 *      - Record mass only for SURFACE_REACHED.
	 * 3. Identify a *stable branch* interval where
	 *      \f$ M(\varepsilon_{c,i+1}) > M(\varepsilon_{c,i}) \f$
	 *    and where \( M(\varepsilon_c) \) brackets the target mass.
	 * 4. If such a bracket exists, perform a bisection search in
	 *    \( \varepsilon_c \) on that monotonic segment until a desired
	 *    mass tolerance is achieved.
	 * 5. If a stable-branch bracketing interval cannot be found (e.g. EOS
	 *    with no monotonic segment covering the target), fail closed. No nearest-sample fallback.
	 *
	 * ### Output
	 * - On **success**, @p out_tov is filled with the radial grid
	 *   \f$ (r, m, \nu', \nu, p, \epsilon, n_B, \rho_i) \f$ produced by TOV
	 *   integration, and the return value is the number of radial points.
	 *
	 * - If @p out_species_labels is non-null, it is filled with the species
	 *   labels inferred from the EOS (typically @ref EOSTable::extra_labels).
	 *
	 * - On **failure**, the function returns `0` and leaves @p out_tov empty.
	 *
	 * @param[in]  target_M_solar
	 *     Target gravitational mass in units of solar masses.
	 *
	 * @param[out] out_tov
	 *     Vector that receives the final TOV radial structure.
	 *     It is cleared on entry.
	 *
	 * @param[out] out_species_labels
	 *     Optional pointer to a vector that receives the species labels
	 *     corresponding to the \f$ \rho_i(r) \f$ columns. May be nullptr.
	 *
	 * @return
	 *     Number of TOV radial points on success, or `0` on failure.
	 *
	 * @pre An EOS must have been loaded via @ref ImportEOS.
	 *
	 * @note This method performs no unit conversions; raw TOV values are
	 *       returned exactly as produced by @ref SingleStarSolveToTOVPoints.
	 *
	 * @see SingleStarSolveToTOVPoints
	 * @see ImportEOS
	 * @see TOVPoint
	 */
	int SolveToProfile(double target_M_solar,
					   std::vector<TOVPoint> &out_tov,
					   std::vector<std::string> *out_species_labels = nullptr);
};

//==============================================================
} // namespace CompactStar::Core
//==============================================================
// std::ostream& operator << ( std::ostream &, const Coord2D&);
// std::ostream& operator << ( std::ostream &, const Segment&);
//==============================================================
#endif /*CompactStar_Core_TOVSolver_H*/
