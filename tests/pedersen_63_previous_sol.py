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


temperature = -42.0  # degC

db = reaktoro.Database('supcrt98.xml')
# _add_hydrocarbons_to_database(db)

# db = reaktoro.Database(str(get_test_data_dir() / 'hydrocarbons.xml'))

editor = reaktoro.ChemicalEditor(db)

gaseous_species = ["CH4(g)", "CO2(g)"]
oil_species = ["CH4(liq)", "CO2(liq)"]
composition = np.array([0.40, 0.60])

def calculate_bips(T):
    k = np.zeros((2, 2))
    k[0, 1] = 0.12
    k[1, 0] = 0.12
    bips = reaktoro.BinaryInteractionParams(k)
    return bips

eos_params = reaktoro.CubicEOSParams(
    model=reaktoro.CubicEOSModel.SoaveRedlichKwong,
    phase_identification_method=reaktoro.PhaseIdentificationMethod.GibbsEnergyAndEquationOfStateMethod,
    binary_interaction_values=calculate_bips
)

editor.addGaseousPhase(gaseous_species).setChemicalModelCubicEOS(eos_params)
editor.addLiquidPhase(oil_species).setChemicalModelCubicEOS(eos_params)

system = reaktoro.ChemicalSystem(editor)
problem = reaktoro.EquilibriumProblem(system)

problem.setTemperature(temperature, 'degC')
problem.add('CH4', composition[0], "mol")
problem.add('CO2', composition[1], "mol")

options = reaktoro.EquilibriumOptions()
options.hessian = reaktoro.GibbsHessian.Exact  # required change for this case

solver = reaktoro.EquilibriumSolver(system)
solver.setOptions(options)

# In this case, we should provide a better initial guess instead of the default.
state = reaktoro.ChemicalState(system)
state.setSpeciesAmounts(0.001)
state.setSpeciesAmount("CH4(g)", composition[0])
state.setSpeciesAmount("CO2(g)", composition[1])

# Retrieve pvtlib solution
data_dir = pathlib.Path(__file__).parent.absolute() / "data"
df_pvtlib_result = pd.read_csv(
    data_dir / "pvtlib-pedersen63-phase-compositions.csv",
)

df_pressures = df_pvtlib_result["Pressure [Pa]"].copy()
df_pressures.drop_duplicates(inplace=True)

df_pvtlib_result_gas = df_pvtlib_result[df_pvtlib_result.Phase == 1]
pvtlib_phase_fractions_gas = df_pvtlib_result_gas["Molar Fraction"].values
pvtlib_pressure_values_gas = df_pvtlib_result_gas["Pressure [Pa]"].values
pvtlib_pressure_values_gas = pvtlib_pressure_values_gas / 100000  # to bar

df_pvtlib_result_liq = df_pvtlib_result[df_pvtlib_result.Phase == 2]
pvtlib_phase_fractions_liquid = df_pvtlib_result_liq["Molar Fraction"].values
pvtlib_pressure_values_liquid = df_pvtlib_result_liq["Pressure [Pa]"].values
pvtlib_pressure_values_liquid = pvtlib_pressure_values_liquid / 100000  # to bar

# Phase molar fractions
pressure_values = df_pressures.values / 100000  # to bar
num_of_components = len(gaseous_species)
num_of_phases = 2
phase_fractions_liquid = list()
phase_fractions_gas = list()
composition_liq = list()
composition_gas = list()
pressure_values_converged = list()
for P in pressure_values:
    problem.setPressure(P, 'bar')
    has_converged = solver.solve(state, problem)
    print(f"P = {P:.3f} bar; T = {temperature} degC; Converged? {has_converged.optimum.succeeded}")

    molar_base = state.phaseAmount('Gaseous') + state.phaseAmount('Liquid')

    gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
    phase_fractions_gas.append(gas_phase_molar_fraction)

    liquid_phase_molar_fraction = state.phaseAmount('Liquid') / molar_base
    phase_fractions_liquid.append(liquid_phase_molar_fraction)

    mixture_properties = state.properties()
    x_gas = mixture_properties.moleFractions().val[0:num_of_components]
    composition_gas.append(x_gas)
    x_liq = mixture_properties.moleFractions().val[num_of_components:num_of_phases * num_of_components]
    composition_liq.append(x_liq)

    pressure_values_converged.append(P)

phase_fractions_gas = np.array(phase_fractions_gas)
phase_fractions_liquid = np.array(phase_fractions_liquid)
composition_gas = np.array(composition_gas)
composition_liq = np.array(composition_liq)
pressure_values_converged = np.array(pressure_values_converged)

plt.figure(figsize=(8, 6))
plt.plot(pressure_values_converged, phase_fractions_liquid, "-x", label="Liquid (Reaktoro)")
plt.plot(pvtlib_pressure_values_liquid, pvtlib_phase_fractions_liquid, "-x", label="Liquid (pvtlib)")

plt.plot(pressure_values_converged, phase_fractions_gas, "-o", label="Gas (Reaktoro)")
plt.plot(pvtlib_pressure_values_gas, pvtlib_phase_fractions_gas, "-o", label="Gas (pvtlib)")

plt.xlabel("Pressure [bar]")
plt.ylabel("Phase molar fraction [mol / mol]")
plt.title(f"Fixed T = {temperature} degC")
plt.legend(shadow=True)

plt.grid(True)

plt.savefig("reaktoro_pvtlib_using_prev_pedersen.png", dpi=300)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(pressure_values_converged, composition_gas[:, 0], "-x", label="C1(g) - Reaktoro")
plt.plot(pvtlib_pressure_values_gas, df_pvtlib_result_gas["C1"], "-x", label="C1(g) - pvtlib")
plt.plot(pressure_values_converged, composition_gas[:, 1], "-x", label="CO2(g) - Reaktoro")
plt.plot(pvtlib_pressure_values_gas, df_pvtlib_result_gas["CO2"], "-x", label="CO2(g) - pvtlib")

plt.ylim([0, 1])
plt.grid(True)

plt.xlabel("Pressure [bar]")
plt.ylabel("Molar fraction [mol / mol]")
plt.title(f"Fixed T = {temperature} degC")
plt.legend(shadow=True, ncol=2)

plt.savefig("reaktoro_pvtlib_compositions_pedersen_gas.png", dpi=300)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(pressure_values_converged, composition_liq[:, 0], "-x", label="C1(liq) - Reaktoro")
plt.plot(pvtlib_pressure_values_liquid, df_pvtlib_result_liq["C1"], "-x", label="C1(liq) - pvtlib")
plt.plot(pressure_values_converged, composition_liq[:, 1], "-x", label="CO2(liq) - Reaktoro")
plt.plot(pvtlib_pressure_values_liquid, df_pvtlib_result_liq["CO2"], "-x", label="CO2(liq) - pvtlib")

plt.ylim([0, 1])
plt.grid(True)

plt.xlabel("Pressure [bar]")
plt.ylabel("Molar fraction [mol / mol]")
plt.title(f"Fixed T = {temperature} degC")
plt.legend(shadow=True, ncol=2)

plt.savefig("reaktoro_pvtlib_compositions_pedersen_liq.png", dpi=300)
plt.show()
