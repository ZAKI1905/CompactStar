/**
 * @file Pulsar.cpp
 * @brief Source for the lightweight, modular Pulsar class.
 *
 * This implementation corresponds to the **new** Pulsar interface in
 * `CompactStar/Core/Pulsar.hpp`, i.e. the one that:
 *
 *  - owns a `StarProfile` (`prof_`)
 *  - exposes a `StarProfileView` (`view_`)
 *  - stores spin/thermal/BNV state in tiny PODs
 *  - delegates *all actual physics* to the headers in `CompactStar/Physics/`
 *
 * The goal is to keep this file SHORT. Anything that starts looking like
 * "derive emissivity as a function of radius" goes into
 * `Physics/Thermal.cpp` (or the header-only version of it).
 */

#include "CompactStar/Core/Pulsar.hpp"
#include "CompactStar/Core/StarBuilder.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Physics/Driver/Spin/MagneticDipole.hpp"
#include "CompactStar/Physics/Spin.hpp"

namespace PhysDriv = CompactStar::Physics::Driver;
namespace CompactStar::Core
{

//======================================================================
//                         Pulsar Class
//======================================================================
// -------------------------------------------------------
//                         CTORS / DTOR
// -------------------------------------------------------
/**
 * @brief Default constructor.
 *
 * We give a generic name to the base `Prog` so logging still works,
 * and we immediately point the view at our internal profile.
 * (The view will be invalid until we actually load a profile.)
 */
Pulsar::Pulsar()
	: Prog("pulsar")
{
	// view_ is a non-owning handle — just point it to our internal storage.
	view_.p = &prof_;
}

// ------------------------------------------------------------
Pulsar::Pulsar(const std::string &name,
			   const Zaki::Math::Quantity &mass,
			   const Physics::State::SpinState &spin)
	: Prog(name),
	  mp(mass),
	  spin_(spin)
{
}

// ------------------------------------------------------------
/// backward-compatible ctor
Pulsar::Pulsar(const std::string &name,
			   const Zaki::Math::Quantity &mass,
			   const Zaki::Math::Quantity &spin_p,
			   const Zaki::Math::Quantity &spin_pdot)
	: Pulsar(
		  name,
		  mass,
		  [spin_p, spin_pdot]()
		  {
			  CompactStar::Physics::State::SpinState s;
			  s.P = spin_p;
			  s.Pdot = spin_pdot;
			  s.mu = {0.0, 0.0};
			  s.d = {0.0, 0.0};

			  s.Resize(1);
			  s.Omega() = 2.0 * M_PI / spin_p.val;

			  return s;
		  }())
{
}

// -------------------------------------------------------
/**
 * @brief Named constructor.
 *
 * @param name Name for logging / workspace.
 */
Pulsar::Pulsar(const std::string &name)
	: Prog(name)
{
	view_.p = &prof_;
}

//==============================================================
//                     Pulsar::FindProfile
//==============================================================

/**
 * @brief Construct the pulsar’s internal structural profile from precomputed
 *        sequence and profile files on disk.
 *
 * This replaces the old, inline logic that used to perform file import,
 * interpolation, and DUrca/blanket calculations directly inside Pulsar.
 * It now delegates those operations to @ref CompactStar::StarBuilder::BuildFromSequence,
 * which handles:
 *  - loading the sequence file `<model_name>_Sequence.tsv`
 *  - identifying the profile index closest to the pulsar mass (`mp`)
 *  - interpolating between neighboring profiles if needed
 *  - computing the structural point @ref seq_point and η_I parameter
 *  - locating the heat-blanket radius (`r_blanket`)
 *  - optionally building a Direct-Urca mask (`r_durca_cond`)
 *
 * The resulting data set is stored in the Pulsar’s internal members for later use
 * by thermodynamic, BNV, or spin modules.
 *
 * @param model_name  Name of the EOS/sequence model (e.g. `"DD2"`, `"SLy4"`).
 * @param in_dir      Relative directory path containing the sequence and profile files.
 *                    Typically this is the same as the one used for the TOV sequence exports.
 *
 * @return Index of the sequence point used (same as returned by
 *         @ref StarBuilder::BuildFromSequence).
 *
 * @throws std::runtime_error if required files are missing or interpolation fails.
 *
 * @note This function automatically exports the resulting profile as
 *       `<wrk_dir>/<in_dir>/<name>.tsv` for consistency with older workflows.
 *       Disable that export if we prefer in-memory use only.
 *
 * @see CompactStar::StarBuilder::BuildFromSequence
 * @see CompactStar::SeqPoint
 * @see CompactStar::StarProfile
 */
int Pulsar::FindProfile(const std::string &model_name,
						const Zaki::String::Directory &in_dir)
{
	// ------------------------------------------------------------
	// 1. Build the star profile using the generic StarBuilder
	// ------------------------------------------------------------
	StarBuilder::Output out;
	StarBuilder::Options opt; // default: blanket_energy_density_km2 = 7.4237e-9, compute_durca_mask = true

	const int seq_idx = StarBuilder::BuildFromSequence(
		wrk_dir_,	// base working directory (inherited from Prog)
		in_dir,		// relative subdirectory holding sequence & profiles
		model_name, // EOS / model label
		mp.val,		// target gravitational mass [M_sun]
		out,		// filled output structure
		opt);

	// ------------------------------------------------------------
	// 2. Copy results into Pulsar members
	// ------------------------------------------------------------
	// profile = std::move(out.profile);
	// seq_point = out.seq_point;
	// eta_I = out.eta_I;
	// r_blanket = out.r_blanket_km;
	// r_blanket_idx = out.r_blanket_idx;
	// r_durca_thresh = out.r_durca_km;
	// r_durca_cond = out.durca_mask;

	// ------------------------------------------------------------
	// 3. Optional: export the resulting profile
	// ------------------------------------------------------------
	out.profile.SetWrkDir(wrk_dir_);
	out.profile.SetPrecision(12);
	out.profile.Export(in_dir + name_ + ".tsv");

	return seq_idx;
}

// -------------------------------------------------------
/**
 * @brief Import an already computed profile.
 *
 * Semantics are the same as `FindProfile(...)`, except we do **not**
 * attempt to run a TOV solve — we just read what’s on disk.
 *
 * @param model_name model/EOS tag
 * @param in_dir     directory relative to work dir
 */
void Pulsar::ImportProfile(const std::string &model_name,
						   const Zaki::String::Directory &in_dir)
{
	const Zaki::String::Directory full_path = (wrk_dir_ + in_dir) + model_name + ".tsv";

	auto edit = prof_.Edit(); // batches all touches below into one bump
	prof_.RadialMutable().Import(full_path);
	view_.p = &prof_;

	// If we later set scalars here, it remains one bump:
	// prof_.SetSurfaceScalars(M, R, z);

	// optional: same comment as in FindProfile(...)
	// const std::string seq_path = (wrk_dir_ + in_dir) + model_name + "_Sequence.tsv";
	// (void)seq_path;
}
//------------------------------------------------------
// 					   Spin interface
//------------------------------------------------------
// Returns the estimated equatorial dipole field [G].
double Pulsar::GetDipoleFieldEstimate() const
{
	return view_.valid() ? Physics::Spin::DipoleFieldEstimate(spin_, view_) : 0.0;
}
//------------------------------------------------------
// Convenience wrapper for the characteristic age \(\tau_c = P/(2\dot{P})\).
// Unit is seconds.
// Returns 0.0 if \(\dot{P}=0\).
double Pulsar::GetCharacteristicAge() const
{
	return Physics::Spin::CharacteristicAge(spin_);
}

//======================================================================
//                     END NAMESPACE
//======================================================================
} // namespace CompactStar