// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// Copyright (C) 2014-2020 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#pragma once

// Reaktoro includes
#include <Reaktoro/Core/ReactionThermoProps.hpp>

namespace Reaktoro {

/// Return a function that calculates thermodynamic properties of a reaction using PHREEQC's analytical expression.
///
/// In this model, the equilibrium constant of the reaction
/// is computed at a given temperature with:
///
/// @eqc{\log_{10}K=A_{1}+A_{2}T+A_{3}T^{-1}+A_{4}\log_{10}T+A_{5}T^{-2}+A_{6}T^{2},}
///
/// where @eq{A_i} are given coefficients. From this model, we can calculate
/// the standard Gibbs energy of reaction using:
///
/// @eqc{\Delta G^{\circ}=-RT\left(A_{1}+A_{2}T+A_{3}T^{-1}+A_{4}\log_{10}T+A_{5}T^{-2}+A_{6}T^{2}\right)\ln_{10}}
///
/// and the standard enthalpy of reaction using:
///
/// @eqc{\Delta H^{\circ}=R\left(A_{2}T^{2}-A_{3}+\frac{A_{4}}{\ln10}T-2A_{5}T^{-1}+2A_{6}T^{3}\right)\ln_{10},}
///
/// considering that:
///
/// @eqc{\Delta G^{\circ}=-RT\ln K}
///
/// and
///
/// @eqc{\Delta H^{\circ}\equiv T^{2}\frac{\partial}{\partial T}\left(-\frac{\Delta G^{\circ}}{T}\right)=RT^{2}\frac{\partial\ln K}{\partial T}.}
///
/// @param A1 The coefficient @eq{A_1} in the reaction thermodynamic model.
/// @param A2 The coefficient @eq{A_2} in the reaction thermodynamic model.
/// @param A3 The coefficient @eq{A_3} in the reaction thermodynamic model.
/// @param A4 The coefficient @eq{A_4} in the reaction thermodynamic model.
/// @param A5 The coefficient @eq{A_5} in the reaction thermodynamic model.
/// @param A6 The coefficient @eq{A_6} in the reaction thermodynamic model.
/// @return ReactionThermoModel
auto ReactionThermoModelAnalyticalPHREEQC(real A1, real A2, real A3, real A4, real A5, real A6) -> ReactionThermoModel;

} // namespace Reaktoro