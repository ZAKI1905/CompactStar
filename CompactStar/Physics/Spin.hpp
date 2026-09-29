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
 * @file Spin.hpp
 * @brief Spin-related derived quantities for pulsars.
 * @ingroup Physics
 *
 * @author Mohammadreza Zakeri
 * Contact: M.Zakeri@eku.edu
 */
#ifndef CompactStar_Physics_Spin_H
#define CompactStar_Physics_Spin_H

#include "CompactStar/Core/StarProfile.hpp"
#include "CompactStar/Physics/State/SpinState.hpp"

/**
 * @namespace CompactStar::Physics::Spin
 * @brief Spin-related physics utilities (characteristic age, dipole fields).
 *
 * This namespace holds lightweight, stateless functions for spin diagnostics
 * using SpinState + structural profiles. Evolution drivers live separately
 * under Physics::Driver::Spin.
 */
namespace CompactStar::Physics::Spin
{
/**
 * @brief Characteristic age \(\tau_c = P/(2\dot{P})\).
 * @ingroup Physics
 * @param s Spin state (uses P and Pdot).
 * @throws std::logic_error This declared diagnostic is not yet implemented/qualified.
 */
double CharacteristicAge(const State::SpinState &s);

/**
 * @brief Magnetic field estimate for a dipole-in-vacuum model.
 * @ingroup Physics
 *
 * Commonly \(B \sim 3.2\times10^{19} \sqrt{P\,\dot{P}}\ \mathrm{G}\) with model factors.
 * This function is a placeholder for our preferred normalization and moments of inertia.
 *
 * @param s    Spin state.
 * @param view Structural profile (for I(M,R) if we include it).
 * @throws std::logic_error The normalization is not yet specified/qualified.
 */
double DipoleFieldEstimate(const State::SpinState &s, Core::StarProfileView view);

} // namespace CompactStar::Physics::Spin

#endif /* CompactStar_Physics_Spin_H */