/**
 * @file StarBuilder.hpp
 * @brief Utilities to construct a neutron-star structural profile from
 *        sequence TSVs (the files we already export from the TOV pipeline).
 *
 * This is the code that used to live in `Pulsar::FindProfile(...)`, but made
 * standalone so that:
 *
 *  1. Pulsar doesn't have to know how profiles are built.
 *  2. MixedStar / NStar / BNV / Thermal code can all reuse the same logic.
 *  3. We can later add other builders (from TOVSolver directly, from JSON,
 *     from MHD runs, etc.) without touching the core containers.
 *
 * The builder does the following:
 *
 *  - load `<wrk_dir><rel_dir><model_name>_Sequence.tsv`
 *  - find the closest mass to the requested one
 *  - if we're at the high-mass end, just load that profile
 *  - otherwise, linearly interpolate *between* the two neighboring profiles
 *    (our logic: pick the shorter profile's radial grid, interpolate
 *     the longer one to it, then blend column-by-column)
 *  - compute `eta_I` using the derivatives in the sequence file
 *  - locate the blanket radius (by energy density)
 *  - build a DUrca “mask” column using the Fermi-momentum condition
 *
 *
 * @ingroup Core
 */

#ifndef CompactStar_Core_StarBuilder_H
#define CompactStar_Core_StarBuilder_H

#include <string>
#include <utility>

#include "CompactStar/Core/SeqPoint.hpp"
#include <Zaki/Physics/Constants.hpp>
#include <Zaki/String/Directory.hpp>
#include <Zaki/Vector/DataColumn.hpp>
#include <Zaki/Vector/DataSet.hpp>

namespace CompactStar
{
namespace Core
{
/**
 * @namespace CompactStar::Core::StarBuilder
 * @brief Namespace for profile-construction helpers.
 */
namespace StarBuilder
{

/**
 * @brief Options controlling how the builder treats blanket / DUrca, etc.
 */
struct Options
{
	/**
	 * @brief Energy density used to define the heat-blanket radius.
	 *
	 * In the old code this was
	 * \f$ \varepsilon \simeq 10^{10}\ \text{g cm}^{-3} \approx 7.4237\times10^{-9}\ \text{km}^{-2} \f$.
	 */
	double blanket_energy_density_km2 = 7.4237e-9;

	/**
	 * @brief Whether to try to build a DUrca mask from composition.
	 *
	 * This requires columns:
	 *  - total baryon density      → column index 5 in our old profile
	 *  - neutron fraction          → label "10"
	 *  - proton fraction           → label "11"
	 *  - electron fraction         → label "0"
	 *
	 * If any of these are missing, the builder will simply leave the mask empty.
	 */
	bool compute_durca_mask = true;
};

/**
 * @brief Full output of the builder.
 *
 * This is everything our old `Pulsar::FindProfile(...)` produced, but
 * gathered into a single struct.
 */
struct Output
{
	/**
	 * @brief The *radial* profile (radius, mass, energy density, composition...).
	 *
	 * This is what we will most often move into `StarProfile` or hand to
	 * `StarProfileView`.
	 */
	Zaki::Vector::DataSet profile;

	/**
	 * @brief The sequence point we landed on (or interpolated to).
	 */
	SeqPoint seq_point;

	/**
	 * @brief Index of the sequence entry we took as “closest”.
	 *
	 * If interpolation was needed, this is still the upper index.
	 */
	int seq_index = -1;

	/**
	 * @brief η_I as defined in [arXiv:2201.02637] logic:
	 *  \f$ \eta_I = \frac{b}{dB/d\varepsilon_c} \frac{dI/d\varepsilon_c}{I} \f$.
	 */
	double eta_I = 0.0;

	/**
	 * @brief Radius (km) where the heat-blanket density threshold is reached.
	 */
	double r_blanket_km = 0.0;

	/**
	 * @brief Index in the radial grid corresponding to @ref r_blanket_km.
	 */
	int r_blanket_idx = -1;

	/**
	 * @brief Radius (km) where the DUrca threshold is reached.
	 *
	 * Only set if we could actually compute the DUrca mask.
	 */
	double r_durca_km = 0.0;

	/**
	 * @brief A 0/1 column that is 1 for r < r_DUrca and 0 outside.
	 *
	 * Same object we used to call `r_durca_cond` in the old code.
	 */
	Zaki::Vector::DataColumn durca_mask;
};

/**
 * @brief Build a neutron-star(-like) profile from a precomputed sequence.
 *
 * This is a direct extraction of our previous logic:
 *
 *  1. load the sequence TSV
 *  2. find the closest mass
 *  3. if we’re at the upper-mass end: just load that profile
 *  4. else: load the neighbor and linearly interpolate
 *  5. compute η_I from sequence derivatives
 *  6. locate blanket radius
 *  7. (optional) build DUrca mask
 *
 * @param wrk_dir    base working directory (what used to be `wrk_dir` in Pulsar)
 * @param rel_dir    relative directory where sequence and profiles live
 *                   (what used to be `in_dir`)
 * @param model_name model/EOS label, e.g. `"DD2"`, `"SLy4"`, ...
 * @param target_mass_Msun target gravitational mass in solar masses
 * @param out        filled on success with profile + seq_point + extras
 * @param opt        optional behavior flags / thresholds
 *
 * @return index in the sequence that was chosen as “closest” (same as old code)
 *
 * @throws std::runtime_error on obvious I/O or interpolation problems
 */
int BuildFromSequence(const Zaki::String::Directory &wrk_dir,
					  const Zaki::String::Directory &rel_dir,
					  const std::string &model_name,
					  double target_mass_Msun,
					  Output &out,
					  const Options &opt = Options());

} // namespace StarBuilder
} // namespace Core
} // namespace CompactStar

#endif // CompactStar_Core_StarBuilder_H