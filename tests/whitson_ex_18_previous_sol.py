import numpy as np
import pandas as pd
import pathlib

import reaktoro

import matplotlib.pyplot as plt


def get_test_data_dir():
    from pathlib import Path
    import os
    return Path(os.path.abspath(__file__)).parents[0] / "data"


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
options.hessian = reaktoro.GibbsHessian.Exact  # required change for this case

solver = reaktoro.EquilibriumSolver(system)
solver.setOptions(options)

# In this case, we should provide a better initial guess instead of the default.
state = reaktoro.ChemicalState(system)
state.setSpeciesAmounts(0.001)  # start will all having 0.001 moles
state.setSpeciesAmount("C1(g)", composition[0])  # overwrite amount of C1(g) (same below)
state.setSpeciesAmount("C4(g)", composition[1])  
state.setSpeciesAmount("C10(g)", composition[2])

# Retrieve pvtlib solution
data_dir = pathlib.Path(__file__).parent.absolute() / "data"
df_pvtlib_result = pd.read_csv(data_dir / "pvtlib_table_whitson18_fixed_T.csv")

df_pressures = df_pvtlib_result["P"].copy()
df_pressures.drop_duplicates(inplace=True)

df_pvtlib_result_gas = df_pvtlib_result[df_pvtlib_result.phase_id == 1]
pvtlib_phase_fractions_gas = df_pvtlib_result_gas.F.values
pvtlib_pressure_values_gas = df_pvtlib_result_gas.P.values

df_pvtlib_result_liq = df_pvtlib_result[df_pvtlib_result.phase_id == 0]
pvtlib_phase_fractions_liquid = df_pvtlib_result_liq.F.values
pvtlib_pressure_values_liquid = df_pvtlib_result_liq.P.values

pressure_values = df_pressures.values
num_of_components = len(gaseous_species)
num_of_phases = 2
phase_fractions_liquid = list()
phase_fractions_gas = list()
composition_liq = list()
composition_gas = list()
pressure_values_converged = list()
for P in pressure_values:
    problem.setPressure(P, 'psi')
    has_converged = solver.solve(state, problem)
    print(f"P = {P:.3f} psi; T = {temperature} degF; Converged? {has_converged.optimum.succeeded}")

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

plt.xlabel("Pressure [psi]")
plt.ylabel("Phase molar fraction [mol / mol]")
plt.title(f"Fixed T = {temperature} degF")
plt.legend(shadow=True)

plt.grid(True)

plt.savefig("reaktoro_pvtlib_using_prev.png", dpi=300)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(pressure_values_converged, composition_gas[:, 0], "-x", label="C1(g) - Reaktoro")
plt.plot(pvtlib_pressure_values_gas, df_pvtlib_result_gas.x0, "-x", label="C1(g) - pvtlib")
plt.plot(pressure_values_converged, composition_gas[:, 1], "-x", label="C4(g) - Reaktoro")
plt.plot(pvtlib_pressure_values_gas, df_pvtlib_result_gas.x1, "-x", label="C4(g) - pvtlib")
plt.plot(pressure_values_converged, composition_gas[:, 2], "-x", label="C10(g) - Reaktoro")
plt.plot(pvtlib_pressure_values_gas, df_pvtlib_result_gas.x2, "-x", label="C10(g) - pvtlib")

plt.ylim([0, 1])
plt.grid(True)

plt.xlabel("Pressure [bar]")
plt.ylabel("Molar fraction [mol / mol]")
plt.title(f"Fixed T = {temperature} degC")
plt.legend(shadow=True, ncol=3)

plt.savefig("reaktoro_pvtlib_compositions_whitson_gas.png", dpi=300)
plt.show()

plt.figure(figsize=(8, 6))
plt.plot(pressure_values_converged, composition_liq[:, 0], "-x", label="C1(liq) - Reaktoro")
plt.plot(pvtlib_pressure_values_liquid, df_pvtlib_result_liq.x0, "-x", label="C1(liq) - pvtlib")
plt.plot(pressure_values_converged, composition_liq[:, 1], "-x", label="C4(liq) - Reaktoro")
plt.plot(pvtlib_pressure_values_liquid, df_pvtlib_result_liq.x1, "-x", label="C4(liq) - pvtlib")
plt.plot(pressure_values_converged, composition_liq[:, 2], "-x", label="C10(liq) - Reaktoro")
plt.plot(pvtlib_pressure_values_liquid, df_pvtlib_result_liq.x2, "-x", label="C10(liq) - pvtlib")

plt.ylim([0, 1])
plt.grid(True)

plt.xlabel("Pressure [bar]")
plt.ylabel("Molar fraction [mol / mol]")
plt.title(f"Fixed T = {temperature} degC")
plt.legend(shadow=True, ncol=3)

plt.savefig("reaktoro_pvtlib_compositions_whitson_liq.png", dpi=300)
plt.show()
