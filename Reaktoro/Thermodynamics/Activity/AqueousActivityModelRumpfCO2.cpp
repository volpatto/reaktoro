// Reaktoro is a C++ library for computational reaction modelling.
//
// Copyright (C) 2014 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "AqueousActivityModelRumpfCO2.hpp"

// Reaktoro includes
#include <Reaktoro/Common/ChemicalScalar.hpp>
#include <Reaktoro/Common/ThermoScalar.hpp>
#include <Reaktoro/Thermodynamics/Mixtures/AqueousMixture.hpp>

namespace Reaktoro {

auto aqueousActivityModelRumpfCO2(const AqueousMixture& mixture) -> AqueousActivityModel
{
    // The number of speciesn and charged species
    const unsigned nspecies = mixture.numSpecies();
    const unsigned nions = mixture.numChargedSpecies();

    // The local indices of some charged species among all charged species
    const Index iNa  = mixture.indexChargedSpecies("Na+");
    const Index iK   = mixture.indexChargedSpecies("K+");
    const Index iCa  = mixture.indexChargedSpecies("Ca++");
    const Index iMg  = mixture.indexChargedSpecies("Mg++");
    const Index iCl  = mixture.indexChargedSpecies("Cl-");

    AqueousActivityModel f = [=](const AqueousMixtureState& state)
    {
        // Extract temperature from the parameters
        const ThermoScalar T = ThermoScalar::Temperature(state.T);

        // The stoichiometric molalities of the ions in the aqueous mixture and their molar derivatives
        const ChemicalVector& ms = state.ms;

        // Extract the stoichiometric molalities of the specific ions and their molar derivatives
        ChemicalScalar zero(nspecies);
        ChemicalScalar mNa = (iNa < nions) ? ms[iNa] : zero;
        ChemicalScalar mK  = (iK  < nions) ? ms[iK]  : zero;
        ChemicalScalar mCa = (iCa < nions) ? ms[iCa] : zero;
        ChemicalScalar mMg = (iMg < nions) ? ms[iMg] : zero;
        ChemicalScalar mCl = (iCl < nions) ? ms[iCl] : zero;

        // The Pitzer's parameters of the Rumpf et al. (1994) model
        const ThermoScalar B = 0.254 - 76.82/T - 10656.0/(T*T) + 6312.0e+3/(T*T*T);
        const double Gamma = -0.0028;

        // The ln activity coefficient of CO2(aq)
        ChemicalScalar ln_gCO2 = 2*B*(mNa + mK + 2*mCa + 2*mMg) + 3*Gamma*(mNa + mK + mCa + mMg)*mCl;

        return ln_gCO2;
    };

    return f;
}

} // namespace Reaktoro
