// Reaktor is a C++ library for computational reaction modelling.
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

#include "PyChemicalSystem.hpp"

// Boost includes
#include <boost/python.hpp>
#include <boost/smart_ptr.hpp>
namespace py = boost::python;

// Reaktor includes
#include <Reaktor/Core/ChemicalSystem.hpp>
#include <Reaktor/Interfaces/Gems.hpp>
#include <Reaktor/Interfaces/Phreeqx.hpp>

namespace Reaktor {
namespace {

auto createChemicalSystemGems(const Gems& gems) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(gems);
}

auto createChemicalSystemPhreeqx(const Phreeqx& phreeqx) -> boost::shared_ptr<ChemicalSystem>
{
    return boost::make_shared<ChemicalSystem>(phreeqx);
}

} // namespace

auto export_ChemicalSystem() -> void
{
    py::class_<ChemicalSystemModels>("ChemicalSystemModels")
        .def_readwrite("phases", &ChemicalSystemModels::phases)
        .def_readwrite("standard_gibbs_energy_fn", &ChemicalSystemModels::standard_gibbs_energy_fn)
        .def_readwrite("standard_enthalpy_fn", &ChemicalSystemModels::standard_enthalpy_fn)
        .def_readwrite("standard_helmholtz_energy_fn", &ChemicalSystemModels::standard_helmholtz_energy_fn)
        .def_readwrite("standard_entropy_fn", &ChemicalSystemModels::standard_entropy_fn)
        .def_readwrite("standard_volume_fn", &ChemicalSystemModels::standard_volume_fn)
        .def_readwrite("standard_internal_energy_fn", &ChemicalSystemModels::standard_internal_energy_fn)
        .def_readwrite("standard_heat_capacity_fn", &ChemicalSystemModels::standard_heat_capacity_fn)
        .def_readwrite("concentration_fn", &ChemicalSystemModels::concentration_fn)
        .def_readwrite("activity_coefficient_fn", &ChemicalSystemModels::activity_coefficient_fn)
        .def_readwrite("activity_fn", &ChemicalSystemModels::activity_fn)
        .def_readwrite("chemical_potential_fn", &ChemicalSystemModels::chemical_potential_fn)
        .def_readwrite("phase_molar_volume_fn", &ChemicalSystemModels::phase_molar_volume_fn)
        ;

    py::class_<Connectivity>("Connectivity")
        .def_readwrite("element_to_species", &Connectivity::element_to_species)
        .def_readwrite("species_to_elements", &Connectivity::species_to_elements)
        .def_readwrite("species_to_phase", &Connectivity::species_to_phase)
        .def_readwrite("phase_to_species", &Connectivity::phase_to_species)
        .def_readwrite("element_to_phases", &Connectivity::element_to_phases)
        .def_readwrite("phase_to_elements", &Connectivity::phase_to_elements)
        ;

    auto element1 = static_cast<const Element&(ChemicalSystem::*)(Index) const>(&ChemicalSystem::element);
    auto element2 = static_cast<const Element&(ChemicalSystem::*)(std::string) const>(&ChemicalSystem::element);

    auto species1 = static_cast<const std::vector<Species>&(ChemicalSystem::*)() const>(&ChemicalSystem::species);
    auto species2 = static_cast<const Species&(ChemicalSystem::*)(Index) const>(&ChemicalSystem::species);
    auto species3 = static_cast<const Species&(ChemicalSystem::*)(std::string) const>(&ChemicalSystem::species);

    auto phase1 = static_cast<const Phase&(ChemicalSystem::*)(Index) const>(&ChemicalSystem::phase);
    auto phase2 = static_cast<const Phase&(ChemicalSystem::*)(std::string) const>(&ChemicalSystem::phase);

    auto indicesElementsInSpecies1 = static_cast<Indices(ChemicalSystem::*)(Index) const>(&ChemicalSystem::indicesElementsInSpecies);
    auto indicesElementsInSpecies2 = static_cast<Indices(ChemicalSystem::*)(const Indices&) const>(&ChemicalSystem::indicesElementsInSpecies);


    py::class_<ChemicalSystem>("ChemicalSystem")
        .def(py::init<>())
        .def(py::init<const ChemicalSystemModels&>())
        .def("__init__", py::make_constructor(createChemicalSystemGems))
        .def("__init__", py::make_constructor(createChemicalSystemPhreeqx))
        .def("numElements", &ChemicalSystem::numElements)
        .def("numSpecies", &ChemicalSystem::numSpecies)
        .def("numSpeciesInPhase", &ChemicalSystem::numSpeciesInPhase)
        .def("numPhases", &ChemicalSystem::numPhases)
        .def("elements", &ChemicalSystem::elements, py::return_value_policy<py::copy_const_reference>())
        .def("species", species1, py::return_value_policy<py::copy_const_reference>())
        .def("phases", &ChemicalSystem::phases, py::return_value_policy<py::copy_const_reference>())
        .def("formulaMatrix", &ChemicalSystem::formulaMatrix, py::return_value_policy<py::copy_const_reference>())
        .def("connectivity", &ChemicalSystem::connectivity, py::return_value_policy<py::copy_const_reference>())
        .def("element", element1, py::return_value_policy<py::copy_const_reference>())
        .def("element", element2, py::return_value_policy<py::copy_const_reference>())
        .def("species", species2, py::return_value_policy<py::copy_const_reference>())
        .def("species", species3, py::return_value_policy<py::copy_const_reference>())
        .def("phase", phase1, py::return_value_policy<py::copy_const_reference>())
        .def("phase", phase2, py::return_value_policy<py::copy_const_reference>())
        .def("indexElement", &ChemicalSystem::indexElement)
        .def("indexElementWithError", &ChemicalSystem::indexElementWithError)
        .def("indexSpecies", &ChemicalSystem::indexSpecies)
        .def("indexSpeciesWithError", &ChemicalSystem::indexSpeciesWithError)
        .def("indexPhase", &ChemicalSystem::indexPhase)
        .def("indexPhaseWithError", &ChemicalSystem::indexPhaseWithError)
        .def("indexPhaseWithSpecies", &ChemicalSystem::indexPhaseWithSpecies)
        .def("indicesElements", &ChemicalSystem::indicesElements)
        .def("indicesSpecies", &ChemicalSystem::indicesSpecies)
        .def("indicesPhases", &ChemicalSystem::indicesPhases)
        .def("indicesPhasesWithSpecies", &ChemicalSystem::indicesPhasesWithSpecies)
        .def("indicesElementsInSpecies", indicesElementsInSpecies1)
        .def("indicesElementsInSpecies", indicesElementsInSpecies2)
        .def("indexFirstSpeciesInPhase", &ChemicalSystem::indexFirstSpeciesInPhase)
        .def("standardGibbsEnergies", &ChemicalSystem::standardGibbsEnergies)
        .def("standardEnthalpies", &ChemicalSystem::standardEnthalpies)
        .def("standardHelmholtzEnergies", &ChemicalSystem::standardHelmholtzEnergies)
        .def("standardEntropies", &ChemicalSystem::standardEntropies)
        .def("standardVolumes", &ChemicalSystem::standardVolumes)
        .def("standardInternalEnergies", &ChemicalSystem::standardInternalEnergies)
        .def("standardHeatCapacities", &ChemicalSystem::standardHeatCapacities)
        .def("concentrations", &ChemicalSystem::concentrations)
        .def("activityCoefficients", &ChemicalSystem::activityCoefficients)
        .def("activity_fn", &ChemicalSystem::activities)
        .def("chemicalPotentials", &ChemicalSystem::chemicalPotentials)
        .def("phaseMolarVolumes", &ChemicalSystem::phaseMolarVolumes)
        .def("phaseDensities", &ChemicalSystem::phaseDensities)
        .def("phaseTotalAmounts", &ChemicalSystem::phaseTotalAmounts)
        .def("phaseVolumes", &ChemicalSystem::phaseVolumes)
        .def("elementAmounts", &ChemicalSystem::elementAmounts)
        .def("elementAmountsInPhase", &ChemicalSystem::elementAmountsInPhase)
        .def("elementAmountsInSpecies", &ChemicalSystem::elementAmountsInSpecies)
        .def("elementAmount", &ChemicalSystem::elementAmount)
        .def("elementAmountInPhase", &ChemicalSystem::elementAmountInPhase)
        .def("elementAmountsInSpecies", &ChemicalSystem::elementAmountsInSpecies)
        .def(py::self_ns::str(py::self_ns::self));
        ;
}

} // namespace Reaktor
