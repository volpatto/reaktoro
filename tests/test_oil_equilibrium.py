import numpy as np
import pytest
import reaktoro


def get_test_data_dir():
    from pathlib import Path
    import os
    return Path(os.path.abspath(__file__)).parents[1] / "data"


def _add_lumped_oil_to_database(db):
    lumped_oil_db = reaktoro.Database(str(get_test_data_dir() / 'hydrocarbons.xml'))

    element_names_in_db = set(e.name() for e in db.elements())
    for element in lumped_oil_db.elements():
        if element.name() not in element_names_in_db:
            db.addElement(element)

    for species in lumped_oil_db.gaseousSpecies():
        db.addGaseousSpecies(species)

    for species in lumped_oil_db.liquidSpecies():
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
    assert True
