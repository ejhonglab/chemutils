#!/usr/bin/env python3

import numpy as np
import pytest
from rdkit import Chem

import chemutils as cu


def _inchi(smiles: str) -> str:
    return Chem.MolToInchi(Chem.MolFromSmiles(smiles))


# Abraham L solute descriptors (log10 K hexadecane/air at 25 C, from measurements).
# All have vapor pressures from thermo's fits to measurements (so PubChem is not used).
abraham_l = {
    # hexane
    'CCCCCC': 2.668,
    # toluene
    'Cc1ccccc1': 3.325,
    # benzene
    'c1ccccc1': 2.786,
    # ethanol
    'CCO': 1.485,
    # methanol
    'CO': 0.970,
    # acetone
    'CC(C)=O': 1.696,
    # ethyl acetate
    'CCOC(C)=O': 2.314,
    # 1-butanol
    'CCCCO': 2.601,
    # benzaldehyde
    'O=Cc1ccccc1': 4.008,
    # diethyl ether
    'CCOCC': 2.015,
}

@pytest.mark.parametrize('smiles,expected_log_k', abraham_l.items())
def test_k_henry_pfo_unifac(smiles, expected_log_k):
    ret = cu.inchi2k_henry_pfo(_inchi(smiles))

    assert ret['k_henry_pfo_method'] == cu.pfo_unifac_method
    assert ret['k_henry_pfo_measured_partitioning'] == False
    assert ret['vapor_pressure_measured'] == True
    assert ret['vapor_pressure_source'].startswith('thermo')
    assert ret['unifac_groups_source'] == 'DDBST'
    assert ret['solvent_surrogate'] == 'hexadecane'

    assert ret['log10_k_oil_air'] == pytest.approx(expected_log_k, abs=0.3)

    assert ret['k_henry_pfo_units'] == 'mol / m^3 / Pa'
    k_oil_air = cu.ureg.Quantity(ret['k_henry_pfo'], ret['k_henry_pfo_units']) * (
        cu.ureg.Quantity(cu.gas_constant, 'J / mol / K') *
        cu.ureg.Quantity(298.15, 'K')
    )
    assert np.log10(k_oil_air.to('dimensionless').magnitude) == pytest.approx(
        ret['log10_k_oil_air']
    )


def test_k_henry_pfo_squalane():
    inchi = _inchi('CCOC(C)=O')
    hexadecane = cu.inchi2k_henry_pfo(inchi)
    squalane = cu.inchi2k_henry_pfo(inchi, solvent='squalane')
    assert squalane['solvent_surrogate'] == 'squalane'
    assert squalane['log10_k_oil_air'] == pytest.approx(hexadecane['log10_k_oil_air'],
        abs=0.3
    )


def test_k_henry_pfo_temperature():
    inchi = _inchi('CCOC(C)=O')
    # (less soluble at higher temperatures)
    assert (cu.inchi2k_henry_pfo(inchi, temp_c=35)['k_henry_pfo'] <
        cu.inchi2k_henry_pfo(inchi)['k_henry_pfo']
    )


def test_k_henry_pfo_measured():
    inchi = _inchi('CCOC(C)=O')
    rt = cu.gas_constant * 298.15

    ret = cu.inchi2k_henry_pfo(inchi, abraham_l=2.314)
    assert ret['k_henry_pfo_method'] == cu.pfo_abraham_l_method
    assert ret['k_henry_pfo_measured_partitioning'] == True
    assert ret['gamma_inf'] is None and ret['vapor_pressure_kpa'] is None
    assert ret['k_henry_pfo'] == pytest.approx(10**2.314 / rt)

    ret = cu.inchi2k_henry_pfo(inchi, log_k_oil_air=2.5, abraham_l=2.314)
    assert ret['k_henry_pfo_method'] == cu.pfo_log_k_method
    assert ret['k_henry_pfo_measured_partitioning'] == True
    assert ret['k_henry_pfo'] == pytest.approx(10**2.5 / rt)

    ret = cu.inchi2k_henry_pfo(inchi, temp_c=35, abraham_l=2.314)
    assert 'not corrected for temperature' in ret['notes']


def test_k_henry_pfo_vapor_pressure_passed_in():
    inchi = _inchi('CCOC(C)=O')
    ret = cu.inchi2k_henry_pfo(inchi)
    psat_kpa = ret['vapor_pressure_kpa']

    ret2 = cu.inchi2k_henry_pfo(inchi, vapor_pressure_kpa=psat_kpa / 2)
    assert ret2['vapor_pressure_source'] == 'passed in'
    assert ret2['k_henry_pfo_method'] == cu.pfo_unifac_method
    assert ret2['k_henry_pfo'] == pytest.approx(2 * ret['k_henry_pfo'])


def test_k_henry_pfo_null():
    assert cu.inchi2k_henry_pfo(None).isnull().all()
