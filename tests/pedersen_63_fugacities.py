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

db = reaktoro.Database('supcrt07.xml')
# _add_hydrocarbons_to_database(db)

# db = reaktoro.Database(str(get_test_data_dir() / 'hydrocarbons.xml'))

editor = reaktoro.ChemicalEditor(db)

gaseous_species = ["CH4(g)", "CO2(g)"]
oil_species = ["CH4(liq)", "CO2(liq)"]

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

state = reaktoro.ChemicalState(system)
state.setTemperature(temperature, 'degC')

## Gathering pvtlib results ##
# Fugacities
data_dir = pathlib.Path(__file__).parent.absolute() / "data"
df_pvtlib_fugacities = pd.read_csv(data_dir / "pvtlib-pedersen63-phase-fugacities.csv")

# Composition
df_pvtlib_compositions = pd.read_csv(data_dir / "pvtlib-pedersen63-phase-compositions.csv")

# Normalizing fugacities
P_1_atm_as_Pa = 101325.0
df_pvtlib_fugacities["C1"] = df_pvtlib_fugacities["C1"] / P_1_atm_as_Pa
df_pvtlib_fugacities["CO2"] = df_pvtlib_fugacities["CO2"] / P_1_atm_as_Pa

# Retrieve results by phase
df_pvtlib_fugacities_gas = df_pvtlib_fugacities[df_pvtlib_fugacities.Phase == 1]
df_pvtlib_fugacities_gas.reset_index(drop=True, inplace=True)
df_pvtlib_fugacities_liq = df_pvtlib_fugacities[df_pvtlib_fugacities.Phase == 2]
df_pvtlib_fugacities_liq.reset_index(drop=True, inplace=True)
df_pvtlib_composition_gas = df_pvtlib_compositions[df_pvtlib_compositions.Phase == 1]
df_pvtlib_composition_gas.reset_index(drop=True, inplace=True)
df_pvtlib_composition_liq = df_pvtlib_compositions[df_pvtlib_compositions.Phase == 2]
df_pvtlib_composition_liq.reset_index(drop=True, inplace=True)

# Properties' calculation loop for gas phase
activities_gas = list()
pvtlib_fugacities_gas = list()
for index, df_row in df_pvtlib_composition_gas.iterrows():
    P = df_row["Pressure [Pa]"]
    composition_0 = df_row["C1"]
    composition_1 = df_row["CO2"]
    state.setPressure(P, 'Pa')
    state.setSpeciesAmount(gaseous_species[0], composition_0)
    state.setSpeciesAmount(gaseous_species[1], composition_1)

    df_row_fugacities = df_pvtlib_fugacities_gas.iloc[index, :]
    fugacity_0 = df_row_fugacities["C1"]
    fugacity_1 = df_row_fugacities["CO2"]
    fugacities = np.array([fugacity_0, fugacity_1])

    gas_properties = state.properties()
    activities = np.exp(gas_properties.lnActivities().val)
    activities_gas.append(activities)
    pvtlib_fugacities_gas.append(fugacities)

activities_gas = np.array(activities_gas)
pvtlib_fugacities_gas = np.array(pvtlib_fugacities_gas)
pressure_values = df_pvtlib_composition_gas["Pressure [Pa]"].values / 100000  # to bar

plt.figure(figsize=(8, 6))

plt.plot(pressure_values, activities_gas[:, 0], "-x", label="C1(g) - Reaktoro")
plt.plot(pressure_values, pvtlib_fugacities_gas[:, 0], "-o", label="C1(g) - pvtlib")
plt.plot(pressure_values, activities_gas[:, 1], "-x", label="CO2(g) - Reaktoro")
plt.plot(pressure_values, pvtlib_fugacities_gas[:, 1], "-o", label="CO2(g) - pvtlib")

plt.xlabel("Pressure [bar]")
plt.ylabel("Fugacities [bar]")
plt.title(f"Fixed T = {temperature} degC")
plt.legend(shadow=True)

plt.grid(True)

plt.savefig("reaktoro_pvtlib_fugacities_gas_pedersen.png", dpi=300)
plt.show()

# Properties' calculation loop for gas phase
activities_liq = list()
pvtlib_fugacities_liq = list()
for index, df_row in df_pvtlib_composition_liq.iterrows():
    P = df_row["Pressure [Pa]"]
    composition_0 = df_row["C1"]
    composition_1 = df_row["CO2"]
    state.setPressure(P, 'Pa')
    state.setSpeciesAmount(oil_species[0], composition_0)
    state.setSpeciesAmount(oil_species[1], composition_1)

    df_row_fugacities = df_pvtlib_fugacities_liq.iloc[index, :]
    fugacity_0 = df_row_fugacities["C1"]
    fugacity_1 = df_row_fugacities["CO2"]
    fugacities = np.array([fugacity_0, fugacity_1])

    liq_properties = state.properties()
    activities = np.exp(liq_properties.lnActivities().val)
    activities_liq.append(activities)
    pvtlib_fugacities_liq.append(fugacities)

activities_liq = np.array(activities_liq)
pvtlib_fugacities_liq = np.array(pvtlib_fugacities_liq)
pressure_values = df_pvtlib_composition_liq["Pressure [Pa]"].values / 100000  # to bar

plt.figure(figsize=(8, 6))

plt.plot(pressure_values, activities_liq[:, 2], "-x", label="C1(liq) - Reaktoro")
plt.plot(pressure_values, pvtlib_fugacities_liq[:, 0], "-o", label="C1(liq) - pvtlib")
plt.plot(pressure_values, activities_liq[:, 3], "-x", label="CO2(liq) - Reaktoro")
plt.plot(pressure_values, pvtlib_fugacities_liq[:, 1], "-o", label="CO2(liq) - pvtlib")

plt.xlabel("Pressure [bar]")
plt.ylabel("Fugacities [bar]")
plt.title(f"Fixed T = {temperature} degC")
plt.legend(shadow=True)

plt.grid(True)

plt.savefig("reaktoro_pvtlib_fugacities_liq_pedersen.png", dpi=300)
plt.show()
