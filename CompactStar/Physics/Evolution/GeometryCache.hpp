// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2025 Mohammadreza Zakeri
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
 * @file GeometryCache.hpp
 * @brief Precomputed geometric factors for fast radial integration (DataColumn-native).
 *
 * Metric convention assumed:
 * \f[
 *   ds^2 = -e^{2\nu(r)} dt^2 + e^{2\Lambda(r)} dr^2 + r^2 d\Omega^2
 * \f]
 *
 * where:
 *  - \f$\nu(r)\f$ is the time metric exponent
 *  - \f$\Lambda(r)\f$ is the radial metric exponent
 *
 * This cache stores primitive, reusable columns:
 *  - \f$r\f$ (km)
 *  - \f$A(r)=4\pi r^2\f$
 *  - \f$e^{\nu}\f$, \f$e^{-\nu}\f$, \f$e^{2\nu}\f$
 *  - \f$e^{\Lambda}\f$, \f$e^{-\Lambda}\f$
 *  - \f$e^{\nu-\Lambda}\f$, \f$e^{-(\nu+\Lambda)}\f$
 *  - \f$w_V(r)=4\pi r^2 e^{\Lambda}\f$     (proper-volume shell weight)
 *  - \f$w_V e^{\nu}\f$, \f$w_V e^{2\nu}\f$ (common redshifted variants)
 *
 * @ingroup PhysicsEvolution
 */

#ifndef CompactStar_Physics_Evolution_GeometryCache_H
#define CompactStar_Physics_Evolution_GeometryCache_H

#include <Zaki/Vector/DataColumn.hpp>
#include <cstddef>

// namespace Zaki
// {
// namespace Vector
// {
// class DataColumn;
// } // namespace Vector
// } // namespace Zaki

namespace CompactStar
{
namespace Physics
{
namespace Evolution
{

class StarContext;

/**
 * @class GeometryCache
 * @brief Geometry-only cached columns for repeated integrals.
 *
 * This class is intended to be:
 *  - Cheap to copy/move (columns store vectors internally),
 *  - Purely geometric (no microphysics),
 *  - DataColumn-native (supports algebra and Mathematica-like indexing).
 */
class GeometryCache
{
  public:
	//--------------------------------------------------------------
	/**
	 * @brief Construct and build cache columns from a StarContext.
	 *
	 * Requirements:
	 *  - ctx.Radius() != nullptr and Size() > 0
	 *  - ctx.Nu()     != nullptr and same Size()
	 *
	 * Lambda handling:
	 *  - If ctx.Lambda() exists and matches Size(), use it.
	 *  - Otherwise derive Lambda from (m,r):
	 *    \f[
	 *      \Lambda = -\tfrac12\ln(1-2m/r)
	 *    \f]
	 *    with a small clamp if (1-2m/r) <= 0.
	 */
	explicit GeometryCache(const StarContext &ctx);

	//--------------------------------------------------------------
	/** @brief Number of radial samples (N). */
	[[nodiscard]] std::size_t Size() const;

	//--------------------------------------------------------------
	/** @name Primitive grids / exponentials */
	///@{
	const Zaki::Vector::DataColumn &R() const;				///< r(km)
	const Zaki::Vector::DataColumn &Mass() const;			///< m(km)
	const Zaki::Vector::DataColumn &Area() const;			///< 4*pi*r^2
	const Zaki::Vector::DataColumn &ExpNu() const;			///< exp(nu)
	const Zaki::Vector::DataColumn &ExpMinusNu() const;		///< exp(-nu)
	const Zaki::Vector::DataColumn &Exp2Nu() const;			///< exp(2*nu)
	const Zaki::Vector::DataColumn &ExpLambda() const;		///< exp(Lambda)
	const Zaki::Vector::DataColumn &ExpMinusLambda() const; ///< exp(-Lambda)
	///@}

	//--------------------------------------------------------------
	/** @name Common mixed metric products */
	///@{
	const Zaki::Vector::DataColumn &ExpNuMinusLambda() const;	   ///< exp(nu - Lambda)
	const Zaki::Vector::DataColumn &ExpMinusNuMinusLambda() const; ///< exp(-(nu + Lambda))
	///@}

	//--------------------------------------------------------------
	/** @name Canonical weights */
	///@{
	const Zaki::Vector::DataColumn &WV() const;		  ///< 4*pi*r^2*exp(Lambda)
	const Zaki::Vector::DataColumn &WVExpNu() const;  ///< WV*exp(nu)
	const Zaki::Vector::DataColumn &WVExp2Nu() const; ///< WV*exp(2nu)
													  ///@}

  private:
	//--------------------------------------------------------------
	/** @brief Build all cached columns (called in ctor). */
	void Build_(const StarContext &ctx);

	//--------------------------------------------------------------
	/**
	 * @brief Derive Lambda from mass and radius, row-by-row, with safety clamp.
	 *
	 * Computes:
	 * \f[
	 *   \Lambda_i = -\tfrac12\ln(\max(1-2m_i/r_i,\epsilon))
	 * \f]
	 *
	 * @param r Radius column (km)
	 * @param m Mass column (km)
	 * @param eps Clamp value for denom <= 0
	 */
	static Zaki::Vector::DataColumn
	DeriveLambdaFromMR_(const Zaki::Vector::DataColumn &r,
						const Zaki::Vector::DataColumn &m,
						double eps = 1e-15);

	//--------------------------------------------------------------
	// Cached columns
	Zaki::Vector::DataColumn m_r;
	Zaki::Vector::DataColumn m_mass; // in km
	Zaki::Vector::DataColumn m_area;

	Zaki::Vector::DataColumn m_expNu;
	Zaki::Vector::DataColumn m_expMinusNu;
	Zaki::Vector::DataColumn m_exp2Nu;

	Zaki::Vector::DataColumn m_expLam;
	Zaki::Vector::DataColumn m_expMinusLam;

	Zaki::Vector::DataColumn m_expNuMinusLam;
	Zaki::Vector::DataColumn m_expMinusNuMinusLam;

	Zaki::Vector::DataColumn m_wV;
	Zaki::Vector::DataColumn m_wVExpNu;
	Zaki::Vector::DataColumn m_wVExp2Nu;
};

} // namespace Evolution
} // namespace Physics
} // namespace CompactStar

#endif /* CompactStar_Physics_Evolution_GeometryCache_H */