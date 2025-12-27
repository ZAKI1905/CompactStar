/**
 * @file StarBuilder.cpp
 * @brief Implementation of profile-building helpers.
 */

#include "CompactStar/Core/StarBuilder.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

#include <Zaki/Physics/Constants.hpp>
// #include <Zaki/String/String.hpp> // for Multiply(...) if we want logging
// #include <Zaki/Vector/DataSet.hpp>

namespace CompactStar::Core
{
namespace StarBuilder
{

/**
 * @brief Small helper: do we have a column with this label?
 */
static bool HasLabel(const Zaki::Vector::DataSet &ds, const std::string &label)
{
	for (const auto &col : ds)
	{
		if (col.Label() == label)
			return true;
	}
	return false;
}

/**
 * @brief Build DUrca mask using our old Fermi-momentum condition.
 *
 * Condition:
 *   kF_n - kF_p - kF_e = 0  → threshold
 *   mask = 1 for r < r_threshold
 *
 * This requires:
 *  - total baryon density  = ds[5]
 *  - neutron fraction      = ds["10"]
 *  - proton fraction       = ds["11"]
 *  - electron fraction     = ds["0"]
 *
 * If any are missing, we just leave the mask empty.
 */
static void BuildDurcaMask(const Zaki::Vector::DataSet &prof,
						   Output &out)
{
	// Basic presence check
	if (prof.Dim().size() < 6)
		return; // no baryon-density column

	if (!HasLabel(prof, "10") || !HasLabel(prof, "11") || !HasLabel(prof, "0"))
		return;

	// Grab what we need
	const Zaki::Vector::DataColumn &r = prof[0];
	const Zaki::Vector::DataColumn &nB = prof[5]; // total baryon density

	// The label-based columns
	Zaki::Vector::DataColumn nFrac = prof["10"];
	Zaki::Vector::DataColumn pFrac = prof["11"];
	Zaki::Vector::DataColumn eFrac = prof["0"];

	// Fermi momenta
	Zaki::Vector::DataColumn kF_n =
		(3 * M_PI * M_PI * nB * nFrac).pow(1.0 / 3.0);
	Zaki::Vector::DataColumn kF_p =
		(3 * M_PI * M_PI * nB * pFrac).pow(1.0 / 3.0);
	Zaki::Vector::DataColumn kF_e =
		(3 * M_PI * M_PI * nB * eFrac).pow(1.0 / 3.0);

	Zaki::Vector::DataColumn kF_diff = kF_n - kF_p - kF_e;

	// threshold index
	int durca_idx = kF_diff.GetClosestIdx(0.0);

	out.r_durca_km = r[durca_idx];

	// build 0/1 column
	out.durca_mask.Resize(r.Size());
	out.durca_mask.Fill(0.0);
	for (std::size_t i = 0; i < static_cast<std::size_t>(durca_idx); ++i)
		out.durca_mask[i] = 1.0;
}

/**
 * @brief Core implementation.
 */
int BuildFromSequence(const Zaki::String::Directory &wrk_dir,
					  const Zaki::String::Directory &rel_dir,
					  const std::string &model_name,
					  double target_mass_Msun,
					  Output &out,
					  const Options &opt)
{
	// ------------------------------------------------------------
	// 1. Load the sequence
	// ------------------------------------------------------------
	const Zaki::String::Directory seq_path =
		(wrk_dir + rel_dir) + model_name + "_Sequence.tsv";

	Zaki::Vector::DataSet seq_ds;
	seq_ds.Import(seq_path);

	if (seq_ds.Dim().empty())
		throw std::runtime_error("StarBuilder::BuildFromSequence: sequence file is empty: " + seq_path.Str());

	// According to our old code:
	//   seq_ds[0] = central energy density
	//   seq_ds[1] = mass [Msun]
	//   seq_ds[4] = total baryon number B(ec)
	//   seq_ds[5] = moment of inertia I(ec)
	//
	// We'll keep that assumption.
	const int mass_col = 1;

	// ------------------------------------------------------------
	// 2. Find closest mass index
	// ------------------------------------------------------------
	const int tmp_idx = seq_ds[mass_col].GetClosestIdx(target_mass_Msun);

	// full row → SeqPoint
	SeqPoint seq_i(seq_ds.GetDataRows()[tmp_idx].vals);

	// ------------------------------------------------------------
	// 3. Load the profile(s)
	//    Profiles live in: <wrk_dir><rel_dir>/profiles/<model>_<idx>.tsv
	// ------------------------------------------------------------
	const Zaki::String::Directory profiles_dir = wrk_dir + rel_dir + "/profiles";

	// current
	Zaki::Vector::DataSet prof_i(profiles_dir,
								 model_name + "_" + std::to_string(tmp_idx) + ".tsv");

	if (prof_i.Dim().empty())
		throw std::runtime_error("StarBuilder::BuildFromSequence: profile file not found or empty for index " + std::to_string(tmp_idx));

	// ------------------------------------------------------------
	// 4. Decide whether we need interpolation
	// ------------------------------------------------------------
	const double closest_mass = seq_ds[mass_col][tmp_idx];
	const double max_mass_in_seq = seq_ds[mass_col].Max();

	if (closest_mass <= target_mass_Msun &&
		tmp_idx == seq_ds[mass_col].MaxIdx())
	{
		// --------------------------------------------------------
		// 4a. We are at/highest mass → just take it
		// --------------------------------------------------------
		out.profile = prof_i;
		out.seq_point = seq_i;
		out.seq_index = tmp_idx;
	}
	else
	{
		// --------------------------------------------------------
		// 4b. We need to interpolate between (tmp_idx - 1) and tmp_idx
		// --------------------------------------------------------
		if (tmp_idx == 0)
			throw std::runtime_error("StarBuilder::BuildFromSequence: cannot interpolate below index 0");

		SeqPoint seq_i_1(seq_ds.GetDataRows()[tmp_idx - 1].vals);

		const double m_i_1 = seq_i_1.m;
		const double m_i = seq_i.m;

		// same formula we had:
		//    x = (m_i - m_pulsar) / (m_i - m_i_1)
		const double x = (m_i - target_mass_Msun) / (m_i - m_i_1);

		// show / log if we want
		std::cout << "----------------------------------------------------\n";
		std::cout << " StarBuilder: closest index = " << tmp_idx << "\n";
		std::cout << "   M_{" << tmp_idx - 1 << "} = " << m_i_1 << "\n";
		std::cout << "   M_{" << tmp_idx << "}   = " << m_i << "\n";
		std::cout << "   x = " << x << "\n";
		std::cout << "----------------------------------------------------\n";

		// The “ideal” sequence point
		out.seq_point = seq_i_1 * x + seq_i * (1.0 - x);
		out.seq_index = tmp_idx;

		// load the neighbor profile
		Zaki::Vector::DataSet prof_i_1(profiles_dir,
									   model_name + "_" + std::to_string(tmp_idx - 1) + ".tsv");

		if (prof_i_1.Dim().empty())
			throw std::runtime_error("StarBuilder::BuildFromSequence: neighbor profile not found for index " + std::to_string(tmp_idx - 1));

		// --------------------------------------------------------
		// Interpolate and blend profiles:
		// pick the shorter radial grid → copy it
		// interpolate the longer one onto it → blend
		// --------------------------------------------------------
		Zaki::Vector::DataSet short_prof;
		Zaki::Vector::DataSet long_prof;
		double x_long = 0.0;
		double x_short = 0.0;

		if (prof_i[0][-1] < prof_i_1[0][-1])
		{
			// i is shorter
			short_prof = prof_i;
			long_prof = prof_i_1;
			x_long = x;
			x_short = 1.0 - x;
		}
		else
		{
			// (i-1) is shorter
			short_prof = prof_i_1;
			long_prof = prof_i;
			x_long = 1.0 - x;
			x_short = x;
		}

		// start output with the short profile's radius column
		out.profile = {short_prof[0]};

		// now blend each column
		for (std::size_t c = 1; c < short_prof.Dim().size(); ++c)
		{
			long_prof.Interpolate(0, static_cast<int>(c));

			Zaki::Vector::DataColumn dc;
			dc.SetLabel(short_prof[c].Label());
			dc.Reserve(short_prof[c].Size());

			for (std::size_t i = 0; i < short_prof[0].Size(); ++i)
			{
				const double r_km = short_prof[0][i];
				const double v_short = short_prof[c][i];
				const double v_long = long_prof.Evaluate(static_cast<int>(c), r_km);

				dc.PushBack(x_long * v_long + x_short * v_short);
			}

			out.profile.AddColumn(std::move(dc));
		}
	}

	// ------------------------------------------------------------
	// 5. Compute eta_I exactly like old code
	//    eta_I = (b / dB/dεc) * (dI/dεc / I)
	// ------------------------------------------------------------
	{
		// we assume:
		// seq_ds[0] = ε_c
		// seq_ds[4] = B(ε_c)
		// seq_ds[5] = I(ε_c)
		seq_ds.Interpolate(0, 4);
		const double dB_over_deps = seq_ds.Derivative(4, out.seq_point.ec);

		seq_ds.Interpolate(0, 5);
		const double dI_over_deps = seq_ds.Derivative(5, out.seq_point.ec);

		double eta_I = out.seq_point.b / dB_over_deps;
		eta_I *= dI_over_deps / out.seq_point.I;

		out.eta_I = eta_I;
	}

	// ------------------------------------------------------------
	// 6. Blanket radius
	//    defined by energy density ≃ 1e10 g/cm^3 → 7.4237e-9 km^-2
	//    This is column 4 in our old profile.
	// ------------------------------------------------------------
	{
		const Zaki::Vector::DataColumn &eps = out.profile[4];
		const Zaki::Vector::DataColumn &r = out.profile[0];

		const int idx = eps.GetClosestIdx(opt.blanket_energy_density_km2);

		out.r_blanket_idx = idx;
		out.r_blanket_km = r[idx];
	}

	// ------------------------------------------------------------
	// 7. DUrca mask (optional)
	// ------------------------------------------------------------
	if (opt.compute_durca_mask)
	{
		BuildDurcaMask(out.profile, out);
	}

	// ------------------------------------------------------------
	// 8. (optional) export — we do NOT export here.
	//    Caller (Pulsar, NStar, …) can do:
	//        out.profile.SetWrkDir(wrk_dir);
	//        out.profile.Export(rel_dir + pulsar_name + ".tsv");
	// ------------------------------------------------------------

	return out.seq_index;
}

} // namespace StarBuilder
} // namespace CompactStar