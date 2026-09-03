// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025
 *   Mohammadreza Zakeri
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
 * @file StarProfile.hpp
 * @brief Unified structural / metric / composition profile for compact stars.
 *
 * This is the “fully populated” phase-2 version. It mirrors the column layout
 * we already used in `NStar`:
 *
 *   0: r         (radius)
 *   1: m         (enclosed mass)
 *   2: nu'       (dν/dr)
 *   3: p         (pressure)
 *   4: eps       (energy density)
 *   5: rho       (total baryon number density)
 *   6: nu        (metric exponent)
 *
 * plus optional:
 *   7: lambda    (the other metric exponent / g_rr)
 *
 * and then **per-species** density columns after that.
 *
 * The idea is: if the TOV solver wrote it, this struct should be able to carry it.
 *
 * @ingroup Core
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@eku.edu
 */

#ifndef CompactStar_Core_StarProfile_H
#define CompactStar_Core_StarProfile_H

#include "CompactStar/Core/SeqPoint.hpp"
#include <Zaki/Vector/DataColumn.hpp>
#include <Zaki/Vector/DataSet.hpp>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include <stdexcept>
#include <string>
#include <vector>
namespace CompactStar::Core
{

//==============================================================
//                        StarProfile
//==============================================================
/**
 * @struct StarProfile
 * @brief Structural / metric / composition data for a compact-star configuration.
 *
 * One instance = one star.
 *
 * It contains:
 *  - full radial data as a DataSet;
 *  - sequence metadata (M, R, central values) via SeqPoint;
 *  - explicit indices for all “core” columns (r, m, ν′, p, ε, ρ, ν, λ);
 *  - a typed `enum class Column` so new code can use enum and old code can use ints;
 *  - per-species density columns, with names and indices.
 */
struct StarProfile
{
  public:
	// ------------------------------------------------------------
	// 0) Versioning / mutation tracking
	// ------------------------------------------------------------
	/**
	 * @brief Monotonic version counter for this profile.
	 *
	 * Increments whenever the profile's structure, column mappings, species registry,
	 * or any other data that cached derived quantities might depend on is modified.
	 *
	 * Consumers (e.g., StarContext caches) can store the last-seen version and
	 * invalidate/rebuild derived caches when Version() changes.
	 */
	std::uint64_t Version() const noexcept { return m_version; }

	/**
	 * @brief Mark the profile as modified (versioning hook).
	 *
	 * This function is the single choke-point for cache invalidation.
	 *
	 * Behavior:
	 *  - If no edit scope is active, it bumps Version() immediately.
	 *  - If an edit scope is active, it only marks the scope "dirty";
	 *    the version is bumped exactly once when the outermost scope ends.
	 *
	 * Use cases:
	 *  - Call from all mutating accessors/mutators (recommended).
	 *  - If external code mutates the profile in a way that bypasses StarProfile
	 *    (should be rare), wrap the mutation inside an EditScope (preferred),
	 *    or call Touch() once after the mutation (fallback).
	 */
	void Touch() noexcept { TouchImpl_(); }

	// ------------------------------------------------------------
	// 0.1) Batched editing (RAII transaction)
	// ------------------------------------------------------------
	/**
	 * @class EditScope
	 * @brief RAII scope that batches multiple StarProfile mutations into one version bump.
	 *
	 * Typical usage:
	 * @code
	 * {
	 *     StarProfile::EditScope edit(prof); // begin edit transaction
	 *     prof.RadialMutable().Import(path); // may call Touch() internally
	 *     prof.SetSurfaceScalars(M, R, z);
	 * } // end transaction -> Version() increments once if anything changed
	 * @endcode
	 *
	 * Notes:
	 *  - Scopes may be nested. Only the outermost scope commits.
	 *  - If no mutation occurs, the version is not bumped.
	 */
	class EditScope
	{
	  public:
		/// Begin an edit scope for @p prof. Supports nesting.
		explicit EditScope(StarProfile &prof) noexcept : m_prof(prof)
		{
			++m_prof.m_editDepth;
		}

		/// End the scope. If this is the outermost scope and anything changed,
		/// bump the version exactly once.
		~EditScope() noexcept
		{
			// Defensive: should never underflow.
			if (m_prof.m_editDepth == 0)
				return;

			--m_prof.m_editDepth;

			// Commit only when leaving the outermost scope.
			if (m_prof.m_editDepth == 0 && m_prof.m_editDirty)
			{
				++m_prof.m_version;
				m_prof.m_editDirty = false;
			}
		}

		EditScope(const EditScope &) = delete;
		EditScope(EditScope &&) noexcept = default;

		EditScope &operator=(const EditScope &) = delete;
		EditScope &operator=(EditScope &&) = delete;

	  private:
		StarProfile &m_prof;
	};

	/**
	 * @brief Convenience creator for an edit scope.
	 *
	 * Usage:
	 * @code
	 * auto edit = prof.Edit();
	 * // ... mutate prof ...
	 * @endcode
	 */
	EditScope Edit() noexcept { return EditScope(*this); }

  private:
	/// Monotonic modification counter used for cache invalidation.
	std::uint64_t m_version = 0;

	/// Active edit-scope nesting depth (0 => no transaction active).
	std::uint32_t m_editDepth = 0;

	/// True if something changed during the current (outermost) edit scope.
	bool m_editDirty = false;

	/** @brief Precision used when exporting profile values. */
	int profile_precision = 9;

	/**
	 * @brief Internal implementation for Touch().
	 *
	 * This keeps the public Touch() semantics stable while allowing
	 * version batching inside EditScope.
	 */
	void TouchImpl_() noexcept
	{
		if (m_editDepth > 0)
		{
			// Defer version bump until end of outermost scope.
			m_editDirty = true;
			return;
		}
		++m_version;
	}

  public:
	// ------------------------------------------------------------
	// 1) Strongly-typed column identifiers
	// ------------------------------------------------------------
	/**
	 * @enum Column
	 * @brief Scoped identifiers for common profile columns.
	 *
	 * We match the legacy layout:
	 *  0: Radius
	 *  1: Mass
	 *  2: MetricNuPrime
	 *  3: Pressure
	 *  4: EnergyDensity
	 *  5: BaryonDensity
	 *  6: MetricNu
	 *  7: MetricLambda
	 */
	enum class Column : int
	{
		Radius = 0,	   ///< r [km]
		Mass,		   ///< m(r) [km]
		MetricNuPrime, ///< dν/dr or ν'(r)
		Pressure,	   ///< p(r)
		EnergyDensity, ///< ε(r)
		BaryonDensity, ///< n_B(r)
		MetricNu,	   ///< ν(r) such that g_tt = -e^{2ν}
		MetricLambda   ///< λ(r) such that g_rr = e^{2λ}
	};

  private:
	// ------------------------------------------------------------
	// 2) Raw radial data
	// ------------------------------------------------------------
	/**
	 * @brief Radial profile as produced by the TOV solver or file importer.
	 *
	 * Column order is expected to **at least** follow the first 7 entries above.
	 */
	Zaki::Vector::DataSet radial;

	// ------------------------------------------------------------
	// 3) Sequence metadata
	// ------------------------------------------------------------
	SeqPoint seq_point{}; ///< Sequence point this profile belongs to.

	// ------------------------------------------------------------
	// 4) Global star properties (surface values)
	// ------------------------------------------------------------
	double M = 0.0;		 ///< Gravitational mass at surface (km).
	double R = 0.0;		 ///< Circumferential radius at surface (km).
	double z_surf = 0.0; ///< Surface redshift factor e^{ν(R)} (dimensionless).

	// ------------------------------------------------------------
	// 5) Backward-compatible integer indices
	// ------------------------------------------------------------
	/**
	 * @name Backward-compatible integer indices
	 * @{
	 */
	int idx_r = static_cast<int>(Column::Radius);
	int idx_m = static_cast<int>(Column::Mass);
	int idx_nuprime = static_cast<int>(Column::MetricNuPrime);
	int idx_p = static_cast<int>(Column::Pressure);
	int idx_eps = static_cast<int>(Column::EnergyDensity);
	int idx_nb = static_cast<int>(Column::BaryonDensity);
	int idx_nu = static_cast<int>(Column::MetricNu);
	int idx_lambda = static_cast<int>(Column::MetricLambda); // may be invalid for older files
	/** @} */

	// ------------------------------------------------------------
	// 6) Per-species densities
	// ------------------------------------------------------------
	/**
	 * @brief Labels for per-species density columns (same order as species_idx).
	 *
	 * Example:
	 *   species_labels = {"n", "p", "Lambda", "Sigma-"}
	 *   species_idx    = { 8,   9,   10,       11      }
	 *
	 * so that radial[ species_idx[i] ] is the column for species_labels[i].
	 */
	std::vector<std::string> species_labels;

	/**
	 * @brief Column indices in `radial` for each species.
	 *
	 * Must be the same length as `species_labels`.
	 */
	std::vector<int> species_idx;

	// ------------------------------------------------------------
	// 7) EOS-derived thermodynamic derivative (ADR-0007 P5, Phase 4C-I0)
	// ------------------------------------------------------------
	/**
	 * @brief `d(eps)/dp` at each radial node — **dimensionless**, EOS-owned.
	 *
	 * Deliberately **not** a column of `radial`. `radial` holds what the TOV integration
	 * produced along the star; this is a property of the *equation of state*, evaluated by its
	 * owner (`TOVSolver::GetEDensDeriv`) on the same `eps(p)` interpolant that built the star
	 * and merely carried here. Keeping it out of `radial` also means the profile's column
	 * count, its species-column indices, its export layout and every file derived from them
	 * are untouched by this addition — the schema is unchanged, and there is still exactly one
	 * authoritative copy of the derivative per profile.
	 *
	 * Empty and `has_eos_dedp_ == false` unless a complete, all-finite set was installed.
	 */
	Zaki::Vector::DataColumn eos_dedp_;

	/// Whether `eos_dedp_` holds a complete, validated, all-finite set for this profile.
	bool has_eos_dedp_ = false;

  public:
	/// Read-only access to the radial dataset (does not change version).
	const Zaki::Vector::DataSet &Radial() const noexcept { return radial; }

	/**
	 * @brief Mutable access to the radial dataset(bumps version).
	 * @return Reference to radial profile data.
	 * @note It will call Touch()
	 */
	Zaki::Vector::DataSet &RadialMutable()
	{
		Touch();
		return radial;
	}

	/**
	 * @brief Get the gravitational mass at the stellar surface.
	 *
	 * @return double Mass at surface (km).
	 */
	double MassSurface() const noexcept { return M; }

	/**
	 * @brief Get the circumferential radius at the stellar surface.
	 *
	 * @return double Radius at surface (km).
	 */
	double RadiusSurface() const noexcept { return R; }

	/**
	 * @brief Get the exponential of the metric function ν at the stellar surface.
	 *
	 * @return double e^{ν(R)} at surface (dimensionless).
	 */
	double ExpNuSurface() const noexcept { return z_surf; }

	/**
	 * @brief Get the sequence point metadata.
	 *
	 * @return const SeqPoint& Reference to the SeqPoint.
	 */
	const SeqPoint &Seq() const noexcept { return seq_point; }

	/**
	 * @brief Mutable access to the SeqPoint (bumps version).
	 * @return Reference to the SeqPoint.
	 * @note It will call Touch()
	 */
	SeqPoint &SeqMutable()
	{
		Touch();
		return seq_point;
	}

	/**
	 * @brief Set the surface properties (bumps version if changed).
	 *
	 * @param M_km Mass at surface (km).
	 * @param R_km Radius at surface (km).
	 * @param zsurf e^{ν(R)} at surface (dimensionless).
	 */
	void SetSurfaceScalars(double M_km, double R_km, double zsurf)
	{
		if (M != M_km || R != R_km || z_surf != zsurf)
		{
			M = M_km;
			R = R_km;
			z_surf = zsurf;
			Touch();
		}
	}

	// ------------------------------------------------------------
	// 7) Basic queries
	// ------------------------------------------------------------
	/**
	 * @brief Test whether the profile has any radial samples.
	 * @return true if there is no column or first column is empty.
	 */
	bool empty() const
	{
		return radial.Empty();
		// auto dims = radial.Dim();
		// if (dims.empty())
		// 	return true;
		// return radial[0].Size() == 0;
	}

	/**
	 * @brief Number of radial grid points.
	 */
	std::size_t size() const
	{
		auto dims = radial.Dim();
		if (dims.empty())
			return 0;
		return radial[0].Size();
	}

	/**
	 * @brief Total number of columns available in the radial profile.
	 */
	std::size_t column_count() const
	{
		return radial.Dim().size();
	}

	/**
	 * @brief Check if a raw column index is valid.
	 */
	bool IsValidColumnIndex(int idx) const noexcept
	{
		if (idx < 0)
			return false;
		return static_cast<std::size_t>(idx) < column_count();
	}

	/**
	 * @brief Check whether the given Column identifier maps to a valid column.
	 */
	bool HasColumn(Column col) const noexcept
	{
		return IsValidColumnIndex(GetColumnIndex(col));
	}

	/**
	 * @brief Whether λ(r) data is available.
	 */
	bool HasMetricLambda() const noexcept
	{
		return HasColumn(Column::MetricLambda);
	}

	// ------------------------------------------------------------
	// 8) Column index management
	// ------------------------------------------------------------
	/**
	 * @brief Get the current integer index that corresponds to a given Column.
	 *
	 * These can be overriden with SetColumnIndex(...) if the file layout
	 * is different.
	 */
	int GetColumnIndex(Column col) const
	{
		switch (col)
		{
		case Column::Radius:
			return idx_r;
		case Column::Mass:
			return idx_m;
		case Column::MetricNuPrime:
			return idx_nuprime;
		case Column::Pressure:
			return idx_p;
		case Column::EnergyDensity:
			return idx_eps;
		case Column::BaryonDensity:
			return idx_nb;
		case Column::MetricNu:
			return idx_nu;
		case Column::MetricLambda:
			return idx_lambda;
		default:
			return -1;
		}
	}

	/**
	 * @brief Override the integer index for a given Column.
	 * @param col Column enum value.
	 * @param idx New integer index.
	 * @note This will call Touch() if the index actually changes.
	 */
	void SetColumnIndex(Column col, int idx)
	{
		bool changed = false;
		switch (col)
		{
		case Column::Radius:
			changed = (idx_r != idx);
			idx_r = idx;
			break;
		case Column::Mass:
			changed = (idx_m != idx);
			idx_m = idx;
			break;
		case Column::MetricNuPrime:
			changed = (idx_nuprime != idx);
			idx_nuprime = idx;
			break;
		case Column::Pressure:
			changed = (idx_p != idx);
			idx_p = idx;
			break;
		case Column::EnergyDensity:
			changed = (idx_eps != idx);
			idx_eps = idx;
			break;
		case Column::BaryonDensity:
			changed = (idx_nb != idx);
			idx_nb = idx;
			break;
		case Column::MetricNu:
			changed = (idx_nu != idx);
			idx_nu = idx;
			break;
		case Column::MetricLambda:
			changed = (idx_lambda != idx);
			idx_lambda = idx;
			break;
		default:
			break;
		}

		if (changed)
			Touch();
	}
	// void SetColumnIndex(Column col, int idx)
	// {
	// 	switch (col)
	// 	{
	// 	case Column::Radius:
	// 		idx_r = idx;
	// 		break;
	// 	case Column::Mass:
	// 		idx_m = idx;
	// 		break;
	// 	case Column::MetricNuPrime:
	// 		idx_nuprime = idx;
	// 		break;
	// 	case Column::Pressure:
	// 		idx_p = idx;
	// 		break;
	// 	case Column::EnergyDensity:
	// 		idx_eps = idx;
	// 		break;
	// 	case Column::BaryonDensity:
	// 		idx_nb = idx;
	// 		break;
	// 	case Column::MetricNu:
	// 		idx_nu = idx;
	// 		break;
	// 	case Column::MetricLambda:
	// 		idx_lambda = idx;
	// 		break;
	// 	}
	// }

	// ------------------------------------------------------------
	// 9) Enum-based data access
	// ------------------------------------------------------------
	Zaki::Vector::DataColumn Get(Column col) const
	{
		int idx = GetColumnIndex(col);
		if (!IsValidColumnIndex(idx))
			throw std::out_of_range("Column index out of range: " +
									std::to_string(static_cast<int>(col)));
		return radial[static_cast<std::size_t>(idx)]; // return by value
	}

	// keep a NON-const version only if we really need mutation;
	// otherwise just delete it. Right now, keep it like this:
	/**
	 * @brief Mutable access to a column by Column enum (bumps version).
	 *
	 * @param col Column enum value.
	 * @return Reference to the column.
	 * @note It will call Touch()
	 */
	Zaki::Vector::DataColumn &Get(Column col)
	{
		Touch();
		int idx = GetColumnIndex(col);
		if (!IsValidColumnIndex(idx))
			throw std::out_of_range("Column index out of range: " +
									std::to_string(static_cast<int>(col)));
		return radial[static_cast<std::size_t>(idx)]; // this is OK if non-const operator[] returns ref
	}

	/**
	 * @brief Get a pointer to the column for a given logical column.
	 *
	 * @param col Column enum value.
	 * This uses const_cast to get the non-const operator[] so we do NOT
	 * return a reference to a temporary.
	 *
	 * @return pointer to column or nullptr if invalid.
	 */
	const Zaki::Vector::DataColumn *GetPtr(Column col) const
	{
		int idx = GetColumnIndex(col);
		if (!IsValidColumnIndex(idx))
			return nullptr;

		// DataSet::operator[] on const returns by value, so we must go through
		// the non-const operator[] to get a stable reference.
		auto &nonconst_ds = const_cast<Zaki::Vector::DataSet &>(radial);
		return &nonconst_ds[static_cast<std::size_t>(idx)];
	}
	// ------------------------------------------------------------
	// 10) Convenience accessors (core)
	// ------------------------------------------------------------
	const Zaki::Vector::DataColumn *GetRadius() const { return GetPtr(Column::Radius); }
	const Zaki::Vector::DataColumn *GetMass() const { return GetPtr(Column::Mass); }
	const Zaki::Vector::DataColumn *GetMetricNuPrime() const { return GetPtr(Column::MetricNuPrime); }
	const Zaki::Vector::DataColumn *GetPressure() const { return GetPtr(Column::Pressure); }
	const Zaki::Vector::DataColumn *GetEnergyDensity() const { return GetPtr(Column::EnergyDensity); }
	const Zaki::Vector::DataColumn *GetBaryonDensity() const { return GetPtr(Column::BaryonDensity); }
	const Zaki::Vector::DataColumn *GetMetricNu() const { return GetPtr(Column::MetricNu); }

	/**
	 * @brief Get λ(r) column, if present.
	 *
	 */
	const Zaki::Vector::DataColumn *GetMetricLambda() const
	{
		if (!HasMetricLambda())
			return nullptr;
		return GetPtr(Column::MetricLambda);
	}

	// -----------------------------------------------------------------
	// LOW-LEVEL POINTER ACCESSORS (for code that needs an actual column)
	// -----------------------------------------------------------------
	/**
	 * @brief Get a pointer to the underlying column, or nullptr if invalid.
	 *
	 * This uses const_cast to get the non-const operator[] so we do NOT
	 * return a reference to a temporary.
	 */
	const Zaki::Vector::DataColumn *GetColumnPtr(int idx) const
	{
		if (!IsValidColumnIndex(idx))
			return nullptr;

		// we know radial is the owner; we just need the non-const operator[]
		auto &ds = const_cast<Zaki::Vector::DataSet &>(radial);
		return &ds[static_cast<std::size_t>(idx)];
	}

	/**
	 * @brief Mutable version of GetColumnPtr.
	 *
	 * Returns a direct pointer to the underlying column if valid,
	 * or nullptr if the index is out of range.
	 */
	Zaki::Vector::DataColumn *GetColumnPtr(int idx)
	{
		if (!IsValidColumnIndex(idx))
			return nullptr;

		Touch();

		return &radial[static_cast<std::size_t>(idx)];
	}

	/**
	 * @brief Get a pointer to the column for a given logical column.
	 */
	const Zaki::Vector::DataColumn *GetColumnPtr(Column col) const
	{
		return GetColumnPtr(GetColumnIndex(col));
	}

	/**
	 * @brief Get a pointer to the column for a given logical column.
	 *
	 * @param col Column enum value.
	 *
	 * Mutable version.
	 * @return pointer to column or nullptr if invalid.
	 */
	Zaki::Vector::DataColumn *GetColumnPtr(Column col)
	{
		return GetColumnPtr(GetColumnIndex(col));
	}

	/**
	 * @brief Get a pointer to a species column by label.
	 */
	const Zaki::Vector::DataColumn *GetSpeciesPtr(const std::string &label) const
	{
		int li = SpeciesLocalIndex(label);
		if (li < 0 || static_cast<std::size_t>(li) >= species_idx.size())
			return nullptr;

		int col_idx = species_idx[static_cast<std::size_t>(li)];
		return GetColumnPtr(col_idx);
	}

	/**
	 * @brief Mutable version of GetSpeciesPtr.
	 * @param label species name, e.g. "n", "p", "Lambda"
	 * @return pointer to column or nullptr if not found.
	 */
	Zaki::Vector::DataColumn *GetSpeciesPtr(const std::string &label)
	{
		int li = SpeciesLocalIndex(label);
		if (li < 0 || static_cast<std::size_t>(li) >= species_idx.size())
			return nullptr;
		// Touch(); // Will be called inside GetColumnPtr(int)

		int col_idx = species_idx[static_cast<std::size_t>(li)];
		return GetColumnPtr(col_idx); // Calls Touch() inside GetColumnPtr(int)
	}
	// ------------------------------------------------------------
	// 11) Per-species helpers
	// ------------------------------------------------------------
	/**
	 * @brief Number of per-species density columns attached to this profile.
	 */
	// std::size_t SpeciesCount() const noexcept
	// {
	// 	return species_labels.size();
	// }
	std::size_t SpeciesCount() const noexcept
	{
		return std::min(species_labels.size(), species_idx.size());
	}

	/**
	 * @brief Check if a species label exists.
	 *
	 * @param label e.g. "n", "p", "Lambda"
	 * @return true if present.
	 */
	bool HasSpecies(const std::string &label) const
	{
		auto it = std::find(species_labels.begin(), species_labels.end(), label);
		return it != species_labels.end();
	}

	/**
	 * @brief Get index in `species_labels` for a given label, or -1 if not found.
	 */
	int SpeciesLocalIndex(const std::string &label) const
	{
		for (std::size_t i = 0; i < species_labels.size(); ++i)
			if (species_labels[i] == label)
				return static_cast<int>(i);
		return -1;
	}

	// j-th species column index (in the radial DataSet), or -1 if out of range.
	int SpeciesColumnIndex(std::size_t j) const noexcept
	{
		if (j >= species_idx.size())
			return -1;
		return species_idx[j];
	}

	/**
	 * @brief Get radial column for a given species label (const).
	 *
	 * @param label species name, e.g. "n", "p", "Lambda"
	 *
	 * Returns by value to avoid returning a reference to a temporary,
	 * since DataSet::operator[] on a const DataSet returns by value.
	 * @throws std::out_of_range if label not found or index invalid.
	 * @return DataColumn for the species.
	 */
	Zaki::Vector::DataColumn GetSpecies(const std::string &label) const
	{
		int li = SpeciesLocalIndex(label);
		if (li < 0 || static_cast<std::size_t>(li) >= species_idx.size())
			throw std::out_of_range("Unknown species label: " + label);

		int col_idx = species_idx[static_cast<std::size_t>(li)];
		if (!IsValidColumnIndex(col_idx))
			throw std::out_of_range("Column index out of range for species: " + label);

		// return by value — SAFE
		return radial[static_cast<std::size_t>(col_idx)];
	}

	/**
	 * @brief Mutable access to species column.
	 */
	Zaki::Vector::DataColumn &GetSpecies(const std::string &label)
	{
		Touch();
		int li = SpeciesLocalIndex(label);
		if (li < 0 || static_cast<std::size_t>(li) >= species_idx.size())
			throw std::out_of_range("Unknown species label: " + label);
		int col_idx = species_idx[static_cast<std::size_t>(li)];
		if (!IsValidColumnIndex(col_idx))
			throw std::out_of_range("Column index out of range for species: " + label);
		return radial[static_cast<std::size_t>(col_idx)];
	}

	/**
	 * @brief Add/register a species column.
	 *
	 * Call this from our TOV importer after pushing the column into `radial`.
	 *
	 * @param label   species name
	 * @param col_idx column index in `radial`
	 * @note This will call Touch().
	 */
	void AddSpecies(const std::string &label, int col_idx)
	{
		species_labels.push_back(label);
		species_idx.push_back(col_idx);
		Touch();
	}

	/**
	 * @brief Set/update the column index for an existing species label.
	 *
	 * If the label is not found, it will be added as a new species.
	 *
	 * @param label   species name
	 * @param col_idx column index in `radial`
	 * @note This will call Touch() if the index changes or if a new species is added.
	 */
	void SetSpeciesColumn(const std::string &label, int col_idx)
	{
		for (std::size_t i = 0; i < species_labels.size(); ++i)
		{
			if (species_labels[i] == label)
			{
				if (species_idx[i] != col_idx)
				{
					species_idx[i] = col_idx;
					Touch();
				}
				// species_idx[i] = col_idx;
				return;
			}
		}
		// fallback: if not found, register new
		AddSpecies(label, col_idx); // Calls Touch()
	}

	// ------------------------------------------------------------
	// 12) Profile precision
	// ------------------------------------------------------------
	/**
	 * @brief Set the Profile Precision object
	 *
	 * @param profile_precision Number of digits after decimal point
	 * @note This will call Touch() if the precision changes.
	 */
	void SetProfilePrecision(const int &profile_precision)
	{
		if (this->profile_precision != profile_precision)
		{
			this->profile_precision = profile_precision;
			Touch();
		}
	}

	// ------------------------------------------------------------
	// 13) Export (in-place; no copying)
	// ------------------------------------------------------------
	void Export(const Zaki::String::Directory &out_dir, int precision = -1);

	// ------------------------------------------------------------
	// 14) Reset
	// ------------------------------------------------------------
	// 12) EOS thermodynamic derivative (ADR-0007 P5, ACCEPTED 2026-09-02)
	// ------------------------------------------------------------
	/**
	 * @brief Whether this profile carries an authoritative `d(eps)/dp` for every node.
	 *
	 * False is the normal state for a profile built without EOS-derivative data — e.g. a
	 * point-constructed analytic star whose caller did not supply one. Consumers that require
	 * the derivative (the governed O(Omega^2) monopole solver) must test this and fail closed;
	 * consumers that do not (all first-order and TOV physics) are unaffected.
	 */
	bool HasEosDEdP() const noexcept { return has_eos_dedp_; }

	/**
	 * @brief The `d(eps)/dp` column — **dimensionless** — or nullptr when absent.
	 *
	 * One value per radial node, in node order. `d(eps)/dp = 1/c_s^2` in `c = 1` units:
	 * causal matter gives `>= 1`, incompressible matter `0` as a formal limit.
	 */
	const Zaki::Vector::DataColumn *GetEosDEdP() const noexcept
	{
		return has_eos_dedp_ ? &eos_dedp_ : nullptr;
	}

	/**
	 * @brief Install a complete EOS-derivative set (bulk path).
	 *
	 * Validated before it is accepted: the size must match the radius column and every value
	 * must be finite. Anything else clears the data and returns false, so that a partial or
	 * corrupt set can never be published as authoritative.
	 *
	 * @return true if the data was accepted.
	 * @note Calls Touch() (ADR-0003: any cache keyed on Version() sees this).
	 */
	bool SetEosDEdP(const Zaki::Vector::DataColumn &in_dedp)
	{
		const int rcol = GetColumnIndex(Column::Radius);
		const std::size_t n_nodes =
			IsValidColumnIndex(rcol)
				? const_cast<Zaki::Vector::DataSet &>(radial)[static_cast<std::size_t>(rcol)].Size()
				: 0;

		if (n_nodes == 0 || in_dedp.Size() != n_nodes)
		{
			ClearEosDEdP();
			return false;
		}
		for (std::size_t i = 0; i < n_nodes; ++i)
		{
			if (!std::isfinite(in_dedp[static_cast<int>(i)]))
			{
				ClearEosDEdP();
				return false;
			}
		}

		eos_dedp_ = in_dedp;
		eos_dedp_.SetLabel("dEdP");
		has_eos_dedp_ = true;
		Touch();
		return true;
	}

	/**
	 * @brief Append one node's `d(eps)/dp` (incremental path).
	 *
	 * Accumulates only; the data does not become authoritative until FinalizeEosDEdP()
	 * validates the completed set. Mirrors how the incremental construction path builds the
	 * radial columns one appended point at a time.
	 * @note Calls Touch().
	 */
	void AppendEosDEdP(double in_dedp)
	{
		has_eos_dedp_ = false; // not authoritative until finalized
		eos_dedp_.PushBack(in_dedp);
		Touch();
	}

	/**
	 * @brief Validate and publish an incrementally accumulated EOS-derivative set.
	 *
	 * Applies exactly the same acceptance test as SetEosDEdP(): complete and all-finite, or
	 * discarded. Returns false when the profile simply has no derivative data, which is a
	 * normal outcome, not an error.
	 */
	bool FinalizeEosDEdP()
	{
		if (eos_dedp_.Size() == 0)
		{
			ClearEosDEdP();
			return false;
		}
		const Zaki::Vector::DataColumn staged = eos_dedp_;
		return SetEosDEdP(staged);
	}

	/**
	 * @brief Drop any EOS-derivative data. @note Calls Touch().
	 */
	void ClearEosDEdP() noexcept
	{
		eos_dedp_ = Zaki::Vector::DataColumn();
		has_eos_dedp_ = false;
		Touch();
	}

	// ------------------------------------------------------------
	/**
	 * @brief Reset profile to an empty state. Invalidates all views.
	 * @note This will call Touch().
	 */
	void Reset()
	{
		radial.ClearRows();
		seq_point.clear();

		M = 0.0;
		R = 0.0;
		z_surf = 0.0;
		// species_labels.clear();
		// species_idx.clear();

		// EOS-derivative data belongs to the profile that was just discarded (ADR-0007 P5).
		eos_dedp_ = Zaki::Vector::DataColumn();
		has_eos_dedp_ = false;

		Touch();
	}

	// Clear + reserve in one place (single semantic mutation)
	void ResetSpecies(std::size_t n_reserve = 0)
	{
		species_labels.clear();
		species_idx.clear();
		if (n_reserve > 0)
		{
			species_labels.reserve(n_reserve);
			species_idx.reserve(n_reserve);
		}
		Touch();
	}
};

//==============================================================
//                        StarProfileView
//==============================================================
/**
 * @struct StarProfileView
 * @brief Non-owning view into a `StarProfile`.
 *
 * Use this in algorithm headers (thermal, BNV, rotochemical) to avoid copying
 * large datasets. We can check `valid()` before accessing.
 */
struct StarProfileView
{
	const StarProfile *p = nullptr;

	bool valid() const { return (p != nullptr) && !p->empty(); }

	int GetColumnIndex(StarProfile::Column col) const
	{
		if (!valid() || !p->HasColumn(col))
			return -1;
		return p->GetColumnIndex(col);
	}

	// const Zaki::Vector::DataColumn *Get(StarProfile::Column col) const
	// {
	// 	if (!valid())
	// 		return nullptr;

	// 	const int idx = p->GetColumnIndex(col);
	// 	if (!p->IsValidColumnIndex(idx))
	// 		return nullptr;

	// 	return p->GetColumnPtr(idx);
	// }

	/**
	 * @brief Get a pointer to the column for a given logical column.
	 *
	 * @param col Column enum value.
	 * @return pointer to column or nullptr if invalid.
	 */
	const Zaki::Vector::DataColumn *GetPtr(StarProfile::Column col) const
	{
		if (!valid())
			return nullptr;

		// get the int index for this column in the underlying profile
		const int idx = p->GetColumnIndex(col);
		if (!p->IsValidColumnIndex(idx))
			return nullptr;

		// use the profile's pointer-returning helper
		return p->GetColumnPtr(idx);
	}

	// Convenience forwards (so thermal/BNV code can stay short)
	// const Zaki::Vector::DataColumn &Radius() const { return Get(StarProfile::Column::Radius); }
	// const Zaki::Vector::DataColumn &Mass() const { return Get(StarProfile::Column::Mass); }
	// const Zaki::Vector::DataColumn &Pressure() const { return Get(StarProfile::Column::Pressure); }
	// const Zaki::Vector::DataColumn &EnergyDensity() const { return Get(StarProfile::Column::EnergyDensity); }
	// const Zaki::Vector::DataColumn &Baryon() const { return Get(StarProfile::Column::BaryonDensity); }
	// const Zaki::Vector::DataColumn &Nu() const { return Get(StarProfile::Column::MetricNu); }

	const Zaki::Vector::DataColumn *GetRadius() const { return GetPtr(StarProfile::Column::Radius); }
	const Zaki::Vector::DataColumn *GetMass() const { return GetPtr(StarProfile::Column::Mass); }
	const Zaki::Vector::DataColumn *GetMetricNuPrime() const { return GetPtr(StarProfile::Column::MetricNuPrime); }
	const Zaki::Vector::DataColumn *GetPressure() const { return GetPtr(StarProfile::Column::Pressure); }
	const Zaki::Vector::DataColumn *GetEnergyDensity() const { return GetPtr(StarProfile::Column::EnergyDensity); }
	const Zaki::Vector::DataColumn *GetBaryonDensity() const { return GetPtr(StarProfile::Column::BaryonDensity); }
	const Zaki::Vector::DataColumn *GetMetricNu() const { return GetPtr(StarProfile::Column::MetricNu); }
};

} // namespace CompactStar::Core

#endif /* CompactStar_Core_StarProfile_H */
