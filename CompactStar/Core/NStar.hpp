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
 * @file NStar.hpp
 * @brief Neutron-star class (phase 2): legacy DataSet API + unified, SAFE StarProfile.
 *
 * This version:
 *  - uses our latest StarProfile (with safety checks, HasColumn, HasSpecies, etc.),
 *  - keeps the legacy ds-based API but comments out the ones we want to migrate away from,
 *  - prefers StarProfile (profile path) when it is non-empty,
 *  - falls back to legacy (ds + SeqPoint) otherwise.
 *
 * @ingroup Core
 * @author Mohammadreza Zakeri
 * @contact M.Zakeri@eku.edu
 */
// Last edit Oct 31, 2025

#ifndef CompactStar_Core_NStar_H
#define CompactStar_Core_NStar_H

#include <Zaki/Physics/Constants.hpp>
#include <Zaki/String/Directory.hpp>
#include <Zaki/Vector/DataSet.hpp>

#include <string_view>
#include <vector>

#include "CompactStar/Core/Prog.hpp"
#include "CompactStar/Core/RotationSolver.hpp"
#include "CompactStar/Core/SeqPoint.hpp"
#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/EOS/CompOSE_EOS.hpp"

namespace CompactStar::Core
{

// Forward declarations
struct TOVPoint;
class TOVSolver;
class RotationSolver;
class Sequence;

//==============================================================
//                        NStar Class
//==============================================================
/**
 * @class NStar
 * @brief Neutron-star container for TOV solutions, rotation, and export.
 *
 * Transitional “phase 2” class:
 *  - **new** path: TOV → StarProfile (preferred)
 *  - **old** path: TOV → ds (kept, but marked for replacement)
 */
class NStar : public Prog
{
	friend class RotationSolver;
	friend class TOVSolver;
	friend class Sequence;

	// Test-only seam (ADR-0006 Q2): lets a validation harness read the owned solver's internal
	// numerical state — the arbitrary seed and the first-order solve counter — without any of
	// it becoming public scientific API. Declared, never defined, in production.
	friend struct RotationSolverTestSeam;

  private:
	// ------------------------------------------------------------
	// 1) Legacy storage (pre-StarProfile)
	// ------------------------------------------------------------
	/** @brief Legacy dataset holding radius, mass, pressure, etc. */
	// Zaki::Vector::DataSet ds;

	/** @brief Integrand dataset for baryon number calculation (legacy). */
	Zaki::Vector::DataSet B_integrand;

	/** @brief Whether legacy surface info has been finalized. */
	bool surface_ready = false;

	/** @brief Rotation solver instance. */
	RotationSolver rot_solver;

	/** @brief Cached sequence point after surface is reached/imported. */
	// SeqPoint sequence;

	/** @brief Cached moment of inertia. */
	double MomI = 0.0;

	/** @brief Precision used when exporting profile values. */
	// int profile_precision = 9;

	/**
	 * @brief Legacy column layout in `ds`:
	 *
	 *  [0] r
	 *  [1] m
	 *  [2] nu'
	 *  [3] p
	 *  [4] eps
	 *  [5] rho
	 *  [6] nu
	 */
	// enum class Col : int
	// {
	// 	R = 0,
	// 	M = 1,
	// 	NuPrime = 2,
	// 	P = 3,
	// 	Eps = 4,
	// 	Rho = 5,
	// 	Nu = 6
	// };
	// static constexpr std::size_t kFixedCols = 7;
	// static constexpr int idx(Col c) noexcept { return static_cast<int>(c); }

	// const Zaki::Vector::DataColumn &col(Col c) const noexcept { return ds[idx(c)]; }
	// Zaki::Vector::DataColumn &col(Col c) noexcept { return ds[idx(c)]; }

	/** @brief Indices for individual species densities (legacy ds style). */
	// std::vector<int> rho_i_idx;

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

	/**
	 * @brief Legacy helper that builds `ds` from TOV points (pre-profile path).
	 */
	void BuildFromTOV(const std::vector<TOVPoint> &in_tov,
					  const std::vector<std::string> *species_labels = nullptr);

	// ------------------------------------------------------------
	// 2) New unified structural / metric / composition profile
	// ------------------------------------------------------------
	/** @brief Unified, safe StarProfile (with p, nu', lambda, species). */
	StarProfile prof_;

	/** @brief EOS pointer used by profile-based solves (not owned). */
	// CompOSE_EOS *eos_ = nullptr;

	/**
	 * @brief Build / register all interpolants on the current StarProfile.
	 *
	 * We must call DataSet::Interpolate
	 * on the profile’s dataset, just like we did on `ds`.
	 *
	 * Safe to call multiple times.
	 */
	void InitInterpolantsFromProfile_();

  public:
	// ------------------------------------------------------------
	// 3) Ctors / Dtor
	// ------------------------------------------------------------
	NStar();
	explicit NStar(const std::vector<TOVPoint> &in_tov_results);
	NStar(const std::vector<TOVPoint> &in_tov,
		  const std::vector<std::string> &species_labels);
	~NStar();

	NStar(const NStar &) = delete;
	NStar &operator=(const NStar &) = delete;

	// ------------------------------------------------------------
	// 4) EOS management
	// ------------------------------------------------------------
	// void AttachEOS(CompOSE_EOS *eos) { eos_ = eos; }
	// CompOSE_EOS *GetEOS() const { return eos_; }

	// ------------------------------------------------------------
	// 5) Legacy TOV init path (keep)
	// ------------------------------------------------------------
	void InitFromTOVSolver(const TOVSolver *in_tov_solver);
	void Append(const TOVPoint &in_tov_pts);
	/**
	 * @brief Debug helper: print column labels and sizes.
	 */
	void PrintProfileColumnSizes() const;

	void FinalizeSurface();
	[[nodiscard]] bool IsSurfaceFinalized() const noexcept { return surface_ready; }

	/// @brief Reset B_integrand, sequence, flags, and profile to an empty state. Invalidates all profile views.
	void Reset();

	// ------------------------------------------------------------
	// 6) profile-based solving / importing
	// ------------------------------------------------------------
	//! Solve TOV for a single star and populate this NStar's StarProfile.
	/*!
	 * @param eos_file Full path (directory + filename) of the EOS table,
	 *                 e.g. `eos_root + eos_name + "/" + eos_name + ".eos"`.
	 * @param target_M_solar Target gravitational mass in units of M_sun.
	 * @param rel_out_dir Optional subdirectory (relative to this NStar's
	 *                    working directory) where the TOV solver can dump
	 *                    any sequence / debug output if it wants.
	 *
	 * @return Number of radial points in the final profile on success,
	 *         or 0 on failure.
	 *
	 * This is a high–level one-shot constructor:
	 *  - builds an internal TOVSolver,
	 *  - imports the EOS from @p eos_file,
	 *  - calls TOVSolver::SolveToProfile(...),
	 *  - then calls BuildFromTOV(...) to populate @c prof_ and all
	 *    derived quantities (B, I, z_surf, etc.).
	 */
	int SolveTOV_Profile(const Zaki::String::Directory &eos_file,
						 double target_M_solar,
						 const Zaki::String::Directory &rel_out_dir = "");

	/**
	 * @brief Import a precomputed StarProfile from disk and store internally.
	 */
	void ImportProfile_Profile(const std::string &model_name,
							   const Zaki::String::Directory &dir);

	/**
	 * @brief Access the owned structural profile (safe).
	 */
	const StarProfile &Profile() const { return prof_; }

	/**
	 * @brief Get a non-owning view into the structural profile.
	 */
	StarProfileView View() const { return {&prof_}; }

	// ------------------------------------------------------------
	// 7) High-level properties (prefer profile, fallback to legacy)
	// ------------------------------------------------------------

	/**
	 * @brief Retrieve the computed sequence point (legacy).
	 */
	// [[nodiscard]] SeqPoint GetSequence() const { return sequence; }

	/// brief Return const reference to the sequence point (from the profile).
	/// details If the profile is empty, returns a static empty SeqPoint.
	const SeqPoint &GetSequence() const noexcept;

	/// brief Mutable access to the sequence point (from the profile).
	/// details. We use this to tweak B, I, etc. after finalize.
	SeqPoint &GetSequence() noexcept;

	/**
	 * @brief Get radius at star surface.
	 */
	[[nodiscard]] double RadiusSurface() const noexcept
	{
		if (!prof_.empty() && prof_.RadiusSurface() > 0.0)
			return prof_.RadiusSurface();
		return GetSequence().r;
	}

	/**
	 * @brief Get mass at star surface in M_sun units.
	 */
	[[nodiscard]] double MassSurface() const noexcept
	{
		if (!prof_.empty() && prof_.MassSurface() > 0.0)
			return prof_.MassSurface() / Zaki::Physics::SUN_M_KM; // in M_Sun units
		return GetSequence().m;									  // in M_Sun units
	}

	/**
	 * @brief Get number of radial grid points.
	 *
	 * Uses the safe helpers we added: first we check `prof_.empty()`,
	 * then we fall back to the legacy `ds` if necessary.
	 */
	[[nodiscard]] std::size_t Size() const noexcept
	{
		if (!prof_.empty())
			return prof_.size();

		// legacy branch
		// auto ds_dims = ds.Dim();
		// if (!ds_dims.empty())
		// 	return ds[0].Size();

		return 0;
	}

	// ------------------------------------------------------------
	// 8) Per-species handling
	// ------------------------------------------------------------
	/**
	 * @brief Check if density data exists for a labeled species (profile-aware).
	 */
	[[nodiscard]] bool HasRho_i(std::string_view label) const noexcept
	{
		// profile-aware check (preferred)
		if (!prof_.empty() && prof_.HasSpecies(std::string(label)))
			return true;

		// legacy: we don’t have label→index here, so we can’t be exact
		return false;
	}

	/**
	 * @brief Get baryon density DataColumn for a labeled species (profile-aware).
	 *
	 * @return pointer to column or nullptr if not found.
	 */
	// 1) CONST VERSION
	// profile-aware, read-only
	const Zaki::Vector::DataColumn *
	GetRho_i(const std::string_view &label) const;
	// {
	// 	if (prof_.empty())
	// 		return nullptr;

	// 	// assume StarProfile has: const DataColumn *GetSpeciesPtr(const std::string &) const
	// 	return prof_.GetSpeciesPtr(std::string(label));
	// }

	/**
	 * @brief Mutable version of species access.
	 */
	// 2) MUTABLE VERSION
	// profile-aware, returns non-const column
	Zaki::Vector::DataColumn *GetRho_i(const std::string &label)
	{
		if (prof_.empty())
			return nullptr;

		// assume StarProfile has: DataColumn *GetSpeciesPtr(const std::string &)
		return prof_.GetSpeciesPtr(label);
	}

	// ------------------------------------------------------------
	// 9) Rotation / moment of inertia (keep)
	// ------------------------------------------------------------
	[[nodiscard]] double Find_MomInertia();

	/**
	 * @brief The star's **seed-free** first-order (Hartle O(Omega)) rotational response.
	 *
	 * Governed by **ADR-0006** (ACCEPTED 2026-09-02), Q4. Available for every star built
	 * through the ordinary construction path, because `Find_MomInertia()` runs there already.
	 *
	 * Contains `I` and the dimensionless shapes `omega_bar(r)/Omega` and
	 * `[d omega_bar/dr](r)/Omega`. The arbitrary internal normalization of the frame-dragging
	 * solve cancels analytically from each of them, so **nothing here depends on the numerical
	 * seed** and nothing here is a physical spin.
	 *
	 * **This star has no implicit angular velocity.** Per ADR-0006's binding clarification, a
	 * physical `Omega`, `J` and `omega_bar(r)` come into existence only when a caller supplies
	 * an explicit `AngularVelocity` — see `RotationAt()`.
	 *
	 * @return A reference owned by this `NStar`; it must not outlive the star, and its
	 *         `r_grid` points into this star's profile (the ADR-0003 lifetime rule).
	 */
	[[nodiscard]] const HartleFirstOrderResponse &RotationResponse() const;

	/**
	 * @brief Materialize the first-order rotation of this star at an explicitly requested
	 *        physical angular velocity.
	 *
	 * Governed by **ADR-0006** Q1/P3. `omega` is a physical angular velocity in rad s^-1,
	 * carried by an explicit type so the unit cannot be lost at the call site:
	 *
	 * @code
	 *   const auto rot = star.RotationAt(AngularVelocity::FromHz(716.0));
	 *   const double J = rot.J;                       // [km^2]
	 *   const double W = rot.OmegaRadPerSecond();     // [s^-1], as requested
	 * @endcode
	 *
	 * This is a **scaling of the response, not a new ODE solve** — the first-order problem is
	 * linear, so calling it at many angular velocities costs no additional integration.
	 * Zero spin is well defined: it yields `Omega = 0`, `J = 0`, `omega_bar == 0` and the
	 * unchanged `I`.
	 *
	 * @throws std::runtime_error if this star has no valid first-order response.
	 */
	[[nodiscard]] PhysicalFirstOrderRotation RotationAt(AngularVelocity omega) const;

	// ------------------------------------------------------------
	void EvaluateNu();

	// ------------------------------------------------------------
	// 10) Export (keep)
	// ------------------------------------------------------------
	void SetProfilePrecision(const int &profile_precision);
	void Export(const Zaki::String::Directory &in_dir);

	// =========================================================
	// ===== INTERPOLATED ACCESSORS (profile-first) ============
	// =========================================================

	/// Metric function ν(r)
	double GetMetricNu(const double &in_r) const;

	/// Mass m(r)
	double GetMass(const double &in_r) const;

	/// Total baryon density n_B(r)
	double GetBaryonDensity(const double &in_r) const;

	/// Energy density ε(r)
	double GetEnergyDensity(const double &in_r) const;

	/// Pressure p(r)
	double GetPressure(const double &in_r) const;

	/// Baryon number integrand 4π r^2 n_B / sqrt(1-2M/r)
	double BaryonNumIntegrand(double in_r) const;
	// ------------------------------------------------------------
	// 11) ===== LEGACY ds-BASED GETTERS (commented, to be replaced) =====
	// ------------------------------------------------------------
	// These are the ones we can grep in VSCode and replace with
	// Profile().Get(...), Profile().GetPressure(), etc.

	// [[nodiscard]] double GetMass(const double& in_r) const;
	// Zaki::Vector::DataColumn*       GetMass() noexcept;
	// const Zaki::Vector::DataColumn* GetMass() const noexcept;

	// [[nodiscard]] double GetRho(const double& in_r) const noexcept;
	// Zaki::Vector::DataColumn*       GetRho() noexcept;
	// const Zaki::Vector::DataColumn* GetRho() const noexcept;

	// [[nodiscard]] double GetEps(const double& in_r) const noexcept;
	// Zaki::Vector::DataColumn*       GetEps() noexcept;
	// const Zaki::Vector::DataColumn* GetEps() const noexcept;

	// [[nodiscard]] double GetPress(const double& in_r) const noexcept;
	// Zaki::Vector::DataColumn*       GetPress() noexcept;
	// const Zaki::Vector::DataColumn* GetPress() const noexcept;

	// [[nodiscard]] double GetNu(const double& in_r) const;
	// Zaki::Vector::DataColumn*       GetNu() noexcept;
	// const Zaki::Vector::DataColumn* GetNu() const noexcept;

	// Zaki::Vector::DataColumn*       GetRadius() noexcept;
	// const Zaki::Vector::DataColumn* GetRadius() const noexcept;

	// double BaryonNumIntegrand(double in_r) const;
};

//==============================================================
//              Candidate removals / replacements
//==============================================================
/*
 * Replace these calls in our codebase:
 *
 * 1) GetMass(...) / GetMass()
 *      → Profile().GetMass()
 *      → Profile().M   (surface)
 *
 * 2) GetRho(...) / GetRho()
 *      → Profile().GetBaryonDensity()
 *
 * 3) GetEps(...) / GetEps()
 *      → Profile().GetEnergyDensity()
 *
 * 4) GetPress(...) / GetPress()
 *      → Profile().GetPressure()
 *
 * 5) GetNu(...) / GetNu()
 *      → Profile().GetMetricNu()
 *
 * 6) GetRadius()
 *      → Profile().GetRadius()
 *
 * 7) BaryonNumIntegrand(...)
 *      → move to a thermal/BNV helper that takes StarProfileView and uses only
 *        safe accessors (HasColumn(...) before Get(...)).
 *
 * NOTE: species
 *  - old code can stay on NStar::GetRho_i(...)
 *  - new code should do: nstar.Profile().GetSpecies("n") (now safe)
 */

} // namespace CompactStar::Core

#endif /* CompactStar_Core_NStar_H */