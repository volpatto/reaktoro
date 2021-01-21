import numpy as np
import pandas as pd

import reaktoro

import matplotlib.pyplot as plt


def get_test_data_dir():
    from pathlib import Path
    import os
    return Path(os.path.abspath(__file__)).parent.absolute() / "data"


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


# Reaktoro basic setup
temperature = 280  # degF
pressure = 500  # psi

molar_base = 1
composition = molar_base * np.array([0.5, 0.42, 0.08])

# db = reaktoro.Database('supcrt07.xml')
# _add_hydrocarbons_to_database(db)

data_dir = get_test_data_dir()
db = reaktoro.Database(str(data_dir / 'hydrocarbons.xml'))

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

state = reaktoro.ChemicalState(system)

problem = reaktoro.EquilibriumProblem(system)

problem.setTemperature(temperature, 'degF')
problem.setElementAmount('C1', composition[0])
problem.setElementAmount('C4', composition[1])
problem.setElementAmount('C10', composition[2])

options = reaktoro.EquilibriumOptions()
options.hessian = reaktoro.GibbsHessian.Exact
options.optimum.max_iterations = 1000
options.optimum.tolerance = 1e-12
options.optimum.output.active = False

solver = reaktoro.EquilibriumSolver(system)
solver.setOptions(options)

# Gathering pvtlib results and some preprocessing
df_pvtlib_result = pd.read_csv(str(data_dir / "pvtlib_table_whitson18_fixed_T.csv"))
columns_to_use = [
    "P", 
    "T", 
    "F", 
    "x0", 
    "x1", 
    "x2", 
    "fugacity_0", 
    "fugacity_1", 
    "fugacity_2",
    "phase_id"
]
df_pvtlib_result = df_pvtlib_result[columns_to_use]

df_pvtlib_result_gas = df_pvtlib_result[df_pvtlib_result.phase_id == 1]
new_composition_names_gas = {
    "x0": "y0", 
    "x1": "y1", 
    "x2": "y2",
    "fugacity_0": "fugacity_0_g", 
    "fugacity_1": "fugacity_1_g", 
    "fugacity_2": "fugacity_2_g",
    "F": "Fv"
}
df_pvtlib_result_gas.rename(columns=new_composition_names_gas, inplace=True)
df_pvtlib_result_gas.reset_index(drop=True, inplace=True)

df_pvtlib_result_liq = df_pvtlib_result[df_pvtlib_result.phase_id == 0]
new_composition_names_liq = {
    "fugacity_0": "fugacity_0_l", 
    "fugacity_1": "fugacity_1_l", 
    "fugacity_2": "fugacity_2_l",
    "F": "Fl"
}
df_pvtlib_result_liq.rename(columns=new_composition_names_liq, inplace=True)
df_pvtlib_result_liq.reset_index(drop=True, inplace=True)

# Merging pvtlib liq and gas results
dfs_to_merge = [df_pvtlib_result_liq, df_pvtlib_result_gas]
df_pvtlib_merged = df_pvtlib_result_gas.merge(df_pvtlib_result_liq, on=["P", "T"], how="outer")
df_pvtlib_merged.fillna(value=0.0, inplace=True)

# Solving using pvtlib results as Initial Guess
phase_fractions_liquid = list()
phase_fractions_gas = list()
pressure_values_converged = df_pvtlib_merged.P.values
for index, df_row in df_pvtlib_merged.iterrows():
    P = df_row.P
    composition_0_g = df_row.y0
    composition_1_g = df_row.y1
    composition_2_g = df_row.y2
    state.setSpeciesAmount(gaseous_species[0], composition_0_g)
    state.setSpeciesAmount(gaseous_species[1], composition_1_g)
    state.setSpeciesAmount(gaseous_species[2], composition_2_g)
    
    composition_0_l = df_row.x0
    composition_1_l = df_row.x1
    composition_2_l = df_row.x2
    state.setSpeciesAmount(oil_species[0], composition_0_l)
    state.setSpeciesAmount(oil_species[1], composition_1_l)
    state.setSpeciesAmount(oil_species[2], composition_2_l)
    
    problem.setPressure(P, 'psi')

    solver.solve(state, problem)

    molar_base = state.phaseAmount('Gaseous') + state.phaseAmount('Liquid')

    gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
    phase_fractions_gas.append(gas_phase_molar_fraction)

    liquid_phase_molar_fraction = state.phaseAmount('Liquid') / molar_base
    phase_fractions_liquid.append(liquid_phase_molar_fraction)

phase_fractions_gas = np.array(phase_fractions_gas)
phase_fractions_liquid = np.array(phase_fractions_liquid)
pressure_values_converged = np.array(pressure_values_converged)

# Plotting the results
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

plt.savefig("reaktoro_pvtlib_with initial_guess.png", dpi=300)
plt.show()