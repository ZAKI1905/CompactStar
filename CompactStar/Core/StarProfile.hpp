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

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "CompactStar/Core/SeqPoint.hpp"
#include <Zaki/Vector/DataColumn.hpp>
#include <Zaki/Vector/DataSet.hpp>
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
  private:
	/** @brief Precision used when exporting profile values. */
	int profile_precision = 9;

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
	 */
	void SetColumnIndex(Column col, int idx)
	{
		switch (col)
		{
		case Column::Radius:
			idx_r = idx;
			break;
		case Column::Mass:
			idx_m = idx;
			break;
		case Column::MetricNuPrime:
			idx_nuprime = idx;
			break;
		case Column::Pressure:
			idx_p = idx;
			break;
		case Column::EnergyDensity:
			idx_eps = idx;
			break;
		case Column::BaryonDensity:
			idx_nb = idx;
			break;
		case Column::MetricNu:
			idx_nu = idx;
			break;
		case Column::MetricLambda:
			idx_lambda = idx;
			break;
		}
	}

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
	Zaki::Vector::DataColumn &Get(Column col)
	{
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

		int col_idx = species_idx[static_cast<std::size_t>(li)];
		return GetColumnPtr(col_idx);
	}
	// ------------------------------------------------------------
	// 11) Per-species helpers
	// ------------------------------------------------------------
	/**
	 * @brief Number of per-species density columns attached to this profile.
	 */
	std::size_t SpeciesCount() const noexcept
	{
		return species_labels.size();
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
	 */
	void AddSpecies(const std::string &label, int col_idx)
	{
		species_labels.push_back(label);
		species_idx.push_back(col_idx);
	}

	/**
	 * @brief Set/update the column index for an existing species label.
	 *
	 * If the label is not found, it will be added as a new species.
	 *
	 * @param label   species name
	 * @param col_idx column index in `radial`
	 */
	void SetSpeciesColumn(const std::string &label, int col_idx)
	{
		for (std::size_t i = 0; i < species_labels.size(); ++i)
		{
			if (species_labels[i] == label)
			{
				species_idx[i] = col_idx;
				return;
			}
		}
		// fallback: if not found, register new
		AddSpecies(label, col_idx);
	}

	// ------------------------------------------------------------
	// 12) Profile precision
	// ------------------------------------------------------------
	void SetProfilePrecision(const int &profile_precision)
	{
		this->profile_precision = profile_precision;
	}

	// ------------------------------------------------------------
	// 13) Export (in-place; no copying)
	// ------------------------------------------------------------
	void Export(const Zaki::String::Directory &out_dir, int precision = -1);

	// ------------------------------------------------------------
	// 14) Reset
	// ------------------------------------------------------------
	/// @brief Reset profile to an empty state. Invalidates all views.
	void Reset()
	{
		radial.ClearRows();
		// radial.data_set.clear();
		seq_point.clear();

		M = 0.0;
		R = 0.0;
		z_surf = 0.0;
		// species_labels.clear();
		// species_idx.clear();
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
