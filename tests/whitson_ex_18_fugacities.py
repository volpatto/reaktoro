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

state = reaktoro.ChemicalState(system)
state.setTemperature(temperature, 'degF')

# Gathering pvtlib results
data_dir = pathlib.Path(__file__).parent.absolute() / "data"
df_pvtlib_result = pd.read_csv(data_dir / "pvtlib_table_whitson18_fixed_T.csv")
df_pvtlib_result_gas = df_pvtlib_result[df_pvtlib_result.phase_id == 1]
df_pvtlib_result_liq = df_pvtlib_result[df_pvtlib_result.phase_id == 0]

# Properties' calculation loop for gas phase
activities_gas = list()
pvtlib_fugacities_gas = list()
for index, df_row in df_pvtlib_result_gas.iterrows():
    P = df_row.P
    composition_0 = df_row.x0
    composition_1 = df_row.x1
    composition_2 = df_row.x2
    state.setPressure(P, 'psi')
    state.setSpeciesAmount(gaseous_species[0], composition_0)
    state.setSpeciesAmount(gaseous_species[1], composition_1)
    state.setSpeciesAmount(gaseous_species[2], composition_2)

    fugacity_0 = df_row.fugacity_0
    fugacity_1 = df_row.fugacity_1
    fugacity_2 = df_row.fugacity_2
    fugacities = np.array([fugacity_0, fugacity_1, fugacity_2])

    gas_properties = state.properties()
    activities = np.exp(gas_properties.lnActivities().val)
    activities_gas.append(activities)
    pvtlib_fugacities_gas.append(fugacities)

activities_gas = np.array(activities_gas)
pvtlib_fugacities_gas = np.array(pvtlib_fugacities_gas)
pressure_values = df_pvtlib_result_gas.P.values

plt.figure(figsize=(8, 6))

plt.plot(pressure_values, activities_gas[:, 0], "-x", label="C1(g) - Reaktoro")
plt.plot(pressure_values, pvtlib_fugacities_gas[:, 0], "-o", label="C1(g) - pvtlib")
plt.plot(pressure_values, activities_gas[:, 1], "-x", label="C4(g) - Reaktoro")
plt.plot(pressure_values, pvtlib_fugacities_gas[:, 1], "-o", label="C4(g) - pvtlib")
plt.plot(pressure_values, activities_gas[:, 2], "-x", label="C10(g) - Reaktoro")
plt.plot(pressure_values, pvtlib_fugacities_gas[:, 2], "-o", label="C10(g) - pvtlib")

plt.xlabel("Pressure [psi]")
plt.ylabel("Fugacities [psi]")
plt.title(f"Fixed T = {temperature} degF")
plt.legend(shadow=True)

plt.grid(True)

plt.savefig("reaktoro_pvtlib_fugacities_gas.png", dpi=300)
plt.show()

# Properties' calculation loop for liq phase
activities_liq = list()
pvtlib_fugacities_liq = list()
for index, df_row in df_pvtlib_result_liq.iterrows():
    P = df_row.P
    composition_0 = df_row.x0
    composition_1 = df_row.x1
    composition_2 = df_row.x2
    state.setPressure(P, 'psi')
    state.setSpeciesAmount(oil_species[0], composition_0)
    state.setSpeciesAmount(oil_species[1], composition_1)
    state.setSpeciesAmount(oil_species[2], composition_2)

    fugacity_0 = df_row.fugacity_0
    fugacity_1 = df_row.fugacity_1
    fugacity_2 = df_row.fugacity_2
    fugacities = np.array([fugacity_0, fugacity_1, fugacity_2])

    gas_properties = state.properties()
    activities = np.exp(gas_properties.lnActivities().val)
    activities_liq.append(activities)
    pvtlib_fugacities_liq.append(fugacities)

activities_liq = np.array(activities_liq)
pvtlib_fugacities_liq = np.array(pvtlib_fugacities_liq)
pressure_values = df_pvtlib_result_liq.P.values

plt.figure(figsize=(8, 6))

plt.plot(pressure_values, activities_liq[:, 0], "-x", label="C1(liq) - Reaktoro")
plt.plot(pressure_values, pvtlib_fugacities_liq[:, 0], "-o", label="C1(liq) - pvtlib")
plt.plot(pressure_values, activities_liq[:, 1], "-x", label="C4(liq) - Reaktoro")
plt.plot(pressure_values, pvtlib_fugacities_liq[:, 1], "-o", label="C4(liq) - pvtlib")
plt.plot(pressure_values, activities_liq[:, 2], "-x", label="C10(liq) - Reaktoro")
plt.plot(pressure_values, pvtlib_fugacities_liq[:, 2], "-o", label="C10(liq) - pvtlib")

plt.xlabel("Pressure [psi]")
plt.ylabel("Fugacities [psi]")
plt.title(f"Fixed T = {temperature} degF")
plt.legend(shadow=True)

plt.grid(True)

plt.savefig("reaktoro_pvtlib_fugacities_liq.png", dpi=300)
plt.show()
