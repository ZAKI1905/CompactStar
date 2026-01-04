#pragma once
/**
 * @file CompOSE_Thermo.hpp
 * @brief Reader + interpolator for CompOSE finite-temperature thermodynamic tables (eos.thermo),
 *        with entropy-based heat capacity support.
 *
 * CompOSE thermodynamic tables parameterize the EOS on a 3D grid:
 *   - temperature        T   [MeV]
 *   - baryon density     nB  [fm^-3]
 *   - charge fraction    Yq  [dimensionless], net strong-sector electric charge per baryon
 *
 * The required thermo columns include:
 *   Q2 = s/nB  (entropy per baryon), dimensionless.
 *
 * From s = nB * Q2 and c_V = T (∂s/∂T)_{nB,Yq}, we obtain:
 *
 *   c_V(T,nB,Yq) = T * nB * (∂Q2/∂T)_{nB,Yq}.
 *
 * This class:
 *   - reads eos.t, eos.nb, eos.yq, eos.thermo,
 *   - stores Q2 on each temperature plane as a dense (nB,Yq) grid,
 *   - evaluates Q2 on planes via bilinear interpolation,
 *   - computes dQ2/dT using temperature-plane finite differences,
 *   - provides CvDensity and CvPerBaryon convenience methods.
 *
 * Equilibrium note:
 *   This object does NOT impose beta equilibrium. Callers provide Yq (e.g., from a solved star
 *   composition profile). Beta equilibrium corresponds to a trajectory Yq^β(T,nB), not a constant Yq.
 */

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace CompactStar
{
namespace EOS
{

/**
 * @brief Parser/interpolator for CompOSE "eos.thermo" tables with emphasis on entropy-based c_V.
 */
class CompOSE_Thermo
{
  public:
	/**
	 * @brief Numeric policy controlling temperature differentiation and bounds handling.
	 *
	 * This is deliberately a simple POD-like type with an explicit default constructor.
	 * That avoids compiler corner-cases with nested default member initializers + default arguments.
	 */
	struct Options
	{
		bool use_central_difference;	///< Use central stencils for dQ2/dT when possible.
		double Tmin_for_derivative_MeV; ///< If >0, snap very small T to the first-interval derivative.
		bool clamp_to_domain;			///< Clamp queries to grid bounds instead of throwing.

		// --- Low-T cooling policy (new) -----------------------------------------
		bool enable_lowT_fit;		 ///< If true, use multi-point fit to get dQ2/dT|_{T->0}.
		int lowT_fit_points;		 ///< Number of first temperature grid points to use (including T=0).
		double lowT_switch_MeV;		 ///< Switch temperature for using fit-based slope.
		double lowT_blend_width_MeV; ///< Blend width around switch to smooth derivative transition.

		/// @brief Default options (robust and conservative).
		Options()
			: use_central_difference(true),
			  Tmin_for_derivative_MeV(0.0),
			  clamp_to_domain(true),
			  enable_lowT_fit(true),
			  lowT_fit_points(3),		// uses T = 0,2,4 MeV by default
			  lowT_switch_MeV(2.0),		// use fit slope up to 2 MeV (configurable)
			  lowT_blend_width_MeV(1.0) // blend from (switch - w) to (switch + w)
		{
		}
	};

  public:
	/// @brief Construct empty; call LoadFromDirectory().
	CompOSE_Thermo() = default;

	/// @brief Construct and load immediately.
	explicit CompOSE_Thermo(const std::string &directory, const Options &opt = Options());

	/// @brief Load/reload from directory containing eos.t, eos.nb, eos.yq, eos.thermo.
	void LoadFromDirectory(const std::string &directory, const Options &opt = Options());

	/// @brief Whether data have been loaded successfully.
	bool IsLoaded() const noexcept { return m_loaded; }

	/// @name Grid accessors
	/// @{
	const std::vector<double> &TGrid_MeV() const { return m_T; }
	const std::vector<double> &NbGrid_fm3() const { return m_nb; }
	const std::vector<double> &YqGrid() const { return m_Yq; }

	std::size_t NT() const noexcept { return m_T.size(); }
	std::size_t NNb() const noexcept { return m_nb.size(); }
	std::size_t NYq() const noexcept { return m_Yq.size(); }
	/// @}

	/// @name Thermodynamic queries
	/// @{

	/**
	 * @brief Evaluate entropy per baryon Q2 = s/nB at (T,nB,Yq).
	 * @return Q2 (dimensionless).
	 */
	double Q2(double T_MeV, double nb_fm3, double Yq) const;

	/**
	 * @brief Entropy per baryon Q2 = s/nB for cooling use.
	 *
	 * For T <= lowT_switch_MeV, uses a low-T linear model
	 *   Q2(T) = a0 * T
	 * where a0 is obtained from a constrained multi-point fit near T=0.
	 *
	 * Above the switch temperature, this blends smoothly to the table-interpolated Q2().
	 */
	double Q2_ForCooling(double T_MeV, double nb_fm3, double Yq) const;

	/**
	 * @brief Evaluate (∂Q2/∂T) at fixed (nB,Yq).
	 * @return dQ2/dT in units of 1/MeV.
	 */
	double dQ2dT(double T_MeV, double nb_fm3, double Yq) const;

	/**
	 * @brief Volumetric heat capacity: c_V = T*nB*dQ2/dT.
	 * @return c_V in native CompOSE (natural) units (MeV * fm^-3 * (1/MeV) = fm^-3).
	 */
	double CvDensity_Natural(double T_MeV, double nb_fm3, double Yq) const;

	/**
	 * @brief Volumetric heat capacity: c_V = T*nB*dQ2/dT.
	 * @return c_V in cV in cgs: erg cm^-3 K^-1.
	 */
	double CvDensity_cgs(double T_MeV, double nb_fm3, double Yq) const;

	/**
	 * @brief Heat capacity per baryon: (c_V/nB) = T*dQ2/dT.
	 */
	double CvPerBaryon(double T_MeV, double nb_fm3, double Yq) const;

	/// @}

	/// @name Low-level utilities
	/// @{

	/**
	 * @brief Evaluate Q2 on a specific temperature plane iT using bilinear interpolation in (nB,Yq).
	 */
	double Q2_OnPlane(std::size_t iT, double nb_fm3, double Yq) const;

	/// @brief Options used by this instance.
	const Options &GetOptions() const noexcept { return m_opt; }

	/// @}

	/**
	 * @brief Low-T improved derivative dQ2/dT for cooling use.
	 *
	 * Uses a multi-point constrained fit near T=0 (Q2(0)=0) to infer the T->0 slope,
	 * and blends to the standard table derivative above lowT_switch_MeV.
	 */
	double dQ2dT_ForCooling(double T_MeV, double nb_fm3, double Yq) const;

	/**
	 * @brief Low-T improved c_V density for cooling use.
	 *
	 * c_V = T * nB * dQ2/dT_ForCooling
	 * c_V in native CompOSE (natural) units (MeV * fm^-3 * (1/MeV) = fm^-3).
	 */
	double CvDensity_Natural_ForCooling(double T_MeV, double nb_fm3, double Yq) const;

	/**
	 * @brief Low-T improved c_V density for cooling use.
	 *
	 * c_V = T * nB * dQ2/dT_ForCooling
	 * c_V in cgs: erg cm^-3 K^-1.
	 */
	double CvDensity_cgs_ForCooling(double T_MeV, double nb_fm3, double Yq) const;

  private:
	void ReadAxes_(const std::string &directory);
	void ReadThermoQ2_(const std::string &directory);

	static std::size_t BracketIndex_(const std::vector<double> &grid, double x);
	void ClampToDomain_(double &T_MeV, double &nb_fm3, double &Yq) const;

	double LowT_Slope_dQ2dT0_(double nb_fm3, double Yq) const;
	static double SmoothStep01_(double x); // x in [0,1]
	double BlendLowT_(double T, double lowT, double highT, double w) const;

  private:
	// Axes
	std::vector<double> m_T;  ///< Temperature grid [MeV]
	std::vector<double> m_nb; ///< Baryon density grid [fm^-3]
	std::vector<double> m_Yq; ///< Charge fraction grid [dimensionless]

	// Q2 values per temperature plane:
	// m_Q2_plane_values[iT] is a flat array sized NNb*NYq with indexing [nb_idx * NYq + yq_idx].
	std::vector<std::vector<double>> m_Q2_plane_values;

	// Metadata (from first line of eos.thermo)
	double m_mn_MeV = 0.0; ///< Neutron mass [MeV]
	double m_mp_MeV = 0.0; ///< Proton mass [MeV]
	int m_Il = 0;		   ///< Lepton flag (1 if leptons included)

	Options m_opt;
	bool m_loaded = false;
};

} // namespace EOS
} // namespace CompactStar