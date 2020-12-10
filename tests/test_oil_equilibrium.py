import numpy as np
import pytest
import reaktoro


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


def test_ternary_c1_c4_c10_mixture():
    """
    This is a ternary example from Whitson monograph. Retrieved from Problem 18 in
    Appendix B.

    .. reference:
        Whitson, C. H., & Brulé, M. R. (2000). Phase behavior (Vol. 20).
        Richardson, TX: Henry L. Doherty Memorial Fund of AIME, Society of Petroleum Engineers.
    """
    temperature = 280  # degF
    pressure = 500  # psi

    molar_base = 1
    composition = molar_base * np.array([0.5, 0.42, 0.08])

    db = reaktoro.Database('supcrt07.xml')
    _add_hydrocarbons_to_database(db)

    editor = reaktoro.ChemicalEditor(db)

    gaseous_species = ["CH4(g)", "C4H10(g)", "C10H22(g)"]
    oil_species = ["CH4(liq)", "C4H10(liq)", "C10H22(liq)"]

    eos_params = reaktoro.CubicEOSParams(
        model=reaktoro.CubicEOSModel.PengRobinson,
        phase_identification_method=reaktoro.PhaseIdentificationMethod.GibbsEnergyAndEquationOfStateMethod,
    )

    editor.addGaseousPhase(gaseous_species).setChemicalModelCubicEOS(eos_params)
    editor.addLiquidPhase(oil_species).setChemicalModelCubicEOS(eos_params)

    system = reaktoro.ChemicalSystem(editor)

    problem = reaktoro.EquilibriumProblem(system)

    problem.setTemperature(temperature, 'degF')
    problem.setPressure(pressure, 'psi')
    problem.add('CH4(g)', composition[0], 'mol')
    problem.add('C4H10(g)', composition[1], 'mol')
    problem.add('C10H22(g)', composition[2], 'mol')

    state = reaktoro.ChemicalState(system)  # for debug purposes only
    solver = reaktoro.EquilibriumSolver(system)

    options = reaktoro.EquilibriumOptions()
    options.hessian = reaktoro.GibbsHessian.Exact
    options.nonlinear.max_iterations = 1000
    options.optimum.max_iterations = 1000
    options.optimum.ipnewton.step = reaktoro.StepMode.Conservative
    options.optimum.tolerance = 1e-18
    solver.setOptions(options)

    state = reaktoro.ChemicalState(system)

    result = solver.solve(state, problem)
    converged = result.optimum.succeeded
    assert converged

    gas_phase_molar_fraction = state.phaseAmount('Gaseous') / molar_base
    liquid_molar_fraction = state.phaseAmount('Liquid') / molar_base
    phase_fractions = np.array([gas_phase_molar_fraction, liquid_molar_fraction])
    phase_fractions_expected = np.array([0.853401, 1 - 0.853401])  # Fv, Fl
    assert phase_fractions == pytest.approx(phase_fractions_expected, rel=1e-3)

    c1_gas_fraction = state.speciesAmount('CH4(g)') / molar_base
    c4_gas_fraction = state.speciesAmount('C4H10(g)') / molar_base
    c10_gas_fraction = state.speciesAmount('C10H22(g)') / molar_base
    equilibrium_gas_composition = np.array([
        c1_gas_fraction,
        c4_gas_fraction,
        c10_gas_fraction
    ])
    expected_equilibrium_gas_composition = np.array([
        0.57114,
        0.41253,
        0.01633
    ])

    c1_liq_fraction = state.speciesAmount('CH4(liq)') / molar_base
    c4_liq_fraction = state.speciesAmount('C4H10(liq)') / molar_base
    c10_liq_fraction = state.speciesAmount('C10H22(liq)') / molar_base
    equilibrium_liq_composition = np.array([
        c1_liq_fraction,
        c4_liq_fraction,
        c10_liq_fraction
    ])
    expected_equilibrium_liq_composition = np.array([
        0.08588,
        0.46349,
        0.45064
    ])

    assert equilibrium_gas_composition == pytest.approx(
        expected_equilibrium_gas_composition,
        rel=1e-2
    )
    assert equilibrium_liq_composition == pytest.approx(
        expected_equilibrium_liq_composition,
        rel=1e-2
    )

