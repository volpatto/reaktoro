import numpy as np
import pandas as pd
import pathlib

import reaktoro

import matplotlib.pyplot as plt


def get_test_data_dir():
    from pathlib import Path
    import os
    return Path(os.path.abspath(__file__)).parents[0] / "data"


def _add_hydrocarbons_to_database(db):
    hydrocarbons_db = reaktoro.Database(str(get_test_data_dir() / 'hydrocarbons.xml'))

    element_names_in_db = set(e.name() for e in db.elements())
    for element in hydrocarbons_db.elements():
        if element.name() not in element_names_in_db:
            db.addElement(element)

    for species in hydrocarbons_db.gaseousSpecies():
        db.addGaseousSpecies(species)

    for species in hydrocarbons_db.liquidSpecies():
        db.addLiquidSpecies(species)


temperature = 280  # degF
pressure = 500  # psi

composition = np.array([0.5, 0.42, 0.08])

# db = reaktoro.Database('supcrt07.xml')
# _add_hydrocarbons_to_database(db)

db = reaktoro.Database(str(get_test_data_dir() / 'hydrocarbons.xml'))

editor = reaktoro.ChemicalEditor(db)

gaseous_species = ["C1(g)", "C4(g)", "C10(g)"]
oil_species = ["C1(liq)", "C4(liq)", "C10(liq)"]

eos_params = reaktoro.CubicEOSParams(
    model=reaktoro.CubicEOSModel.PengRobinson,
    phase_identification_method=reaktoro.PhaseIdentificationMethod.GibbsEnergyAndEquationOfStateMethod,
)

editor.addGaseousPhase(gaseous_species).setChemicalModelCubicEOS(eos_params)
editor.addLiquidPhase(oil_species).setChemicalModelCubicEOS(eos_params)

system = reaktoro.ChemicalSystem(editor)

problem = reaktoro.EquilibriumProblem(system)

problem.setTemperature(temperature, 'degF')
problem.setElementAmount('C1', composition[0])
problem.setElementAmount('C4', composition[1])
problem.setElementAmount('C10', composition[2])

options = reaktoro.EquilibriumOptions()
options.hessian = reaktoro.GibbsHessian.Exact
options.optimum.max_iterations = 1000
options.optimum.tolerance = 1e-18
options.optimum.output.active = False

solver = reaktoro.EquilibriumSolver(system)
solver.setOptions(options)

state = reaktoro.ChemicalState(system)
pressure_values = np.linspace(50, 2000)
phase_fractions_liquid = list()
phase_fractions_gas = list()
pressure_values_converged = list()
for P in pressure_values:
    problem.setPressure(P, 'psi')
    solver.solve(state, problem)

    molar_base = state.phaseAmount('Gaseous') + state.phaseAmount('Liquid')

    gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
    phase_fractions_gas.append(gas_phase_molar_fraction)

    liquid_phase_molar_fraction = state.phaseAmount('Liquid') / molar_base
    phase_fractions_liquid.append(liquid_phase_molar_fraction)

    pressure_values_converged.append(P)

    # Use previous solution as Initial Guess
    previous_composition_x0 = state.speciesAmount("C1(liq)")
    previous_composition_x1 = state.speciesAmount("C4(liq)")
    previous_composition_x2 = state.speciesAmount("C10(liq)")
    state.setSpeciesAmount(oil_species[0], previous_composition_x0)
    state.setSpeciesAmount(oil_species[1], previous_composition_x1)
    state.setSpeciesAmount(oil_species[2], previous_composition_x2)

    previous_composition_y0 = state.speciesAmount("C1(g)")
    previous_composition_y1 = state.speciesAmount("C4(g)")
    previous_composition_y2 = state.speciesAmount("C10(g)")
    state.setSpeciesAmount(gaseous_species[0], previous_composition_y0)
    state.setSpeciesAmount(gaseous_species[1], previous_composition_y1)
    state.setSpeciesAmount(gaseous_species[2], previous_composition_y2)

phase_fractions_gas = np.array(phase_fractions_gas)
phase_fractions_liquid = np.array(phase_fractions_liquid)
pressure_values_converged = np.array(pressure_values_converged)

data_dir = pathlib.Path(__file__).parent.absolute() / "data"
df_pvtlib_result = pd.read_csv(data_dir / "pvtlib_table_whitson18_fixed_T.csv")

pvtlib_phase_fractions_gas = df_pvtlib_result[df_pvtlib_result.phase_id == 1].F.values
pvtlib_pressure_values_gas = df_pvtlib_result[df_pvtlib_result.phase_id == 1].P.values

pvtlib_phase_fractions_liquid = df_pvtlib_result[df_pvtlib_result.phase_id == 0].F.values
pvtlib_pressure_values_liquid = df_pvtlib_result[df_pvtlib_result.phase_id == 0].P.values

plt.figure(figsize=(8, 6))
plt.plot(pressure_values_converged, phase_fractions_liquid, "-x", label="Liquid (Reaktoro)")
plt.plot(pvtlib_pressure_values_liquid, pvtlib_phase_fractions_liquid, "-x", label="Liquid (pvtlib)")

plt.plot(pressure_values_converged, phase_fractions_gas, "-o", label="Gas (Reaktoro)")
plt.plot(pvtlib_pressure_values_gas, pvtlib_phase_fractions_gas, "-o", label="Gas (pvtlib)")

plt.xlabel("Pressure [psi]")
plt.ylabel("Phase molar fraction [mol / mol]")
plt.title(f"Fixed T = {temperature} degF")
plt.legend(shadow=True)

plt.grid(True)

plt.savefig("reaktoro_pvtlib_using_prev.png", dpi=300)
plt.show()
