import numpy as np
import pandas as pd

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

molar_base = 1
composition = molar_base * np.array([0.5, 0.42, 0.08])

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
options.optimum.tolerance = 1e-10
options.optimum.output.active = False

state = reaktoro.equilibrate(problem, options)
solver = reaktoro.EquilibriumSolver(system)
solver.setOptions(options)

pressure_values = np.linspace(50, 2000)
phase_fractions_liquid = list()
phase_fractions_gas = list()
pressure_values_converged = list()
for P in pressure_values:
    problem.setPressure(P, 'psi')
    try:
        solver.solve(state, problem)

        molar_base = state.phaseAmount('Gaseous') + state.phaseAmount('Liquid')

        gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
        phase_fractions_gas.append(gas_phase_molar_fraction)

        liquid_phase_molar_fraction = state.phaseAmount('Liquid') / molar_base
        phase_fractions_liquid.append(liquid_phase_molar_fraction)

        pressure_values_converged.append(P)
    except RuntimeError:
        pass

phase_fractions_gas = np.array(phase_fractions_gas)
phase_fractions_liquid = np.array(phase_fractions_liquid)
pressure_values_converged = np.array(pressure_values_converged)

df_pvtlib_result = pd.read_csv("./data/pvtlib_table_whitson18_fixed_T.csv")

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

plt.savefig("reaktoro_pvtlib_comparison.png", dpi=300)
plt.show()
