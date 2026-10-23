#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pubchempy as pcp
import pytest

from chemutils import is_one2one, compound2cas, allowed_kwarg
import chemutils as cu


def _is_one2one(xs, ys, **kwargs):
    cols = ['x', 'y']
    df = pd.DataFrame(dict(zip(cols, [xs, ys])))
    r1 = is_one2one(df, *cols, **kwargs)
    r2 = is_one2one(df, *cols[::-1], **kwargs)
    assert r1 == r2
    return r1


def test_is_one2one():
    assert _is_one2one([1,2], [2,4])

    # if either col has >1 unique value for any 1 value of the other col, is_one2one
    # should return False
    assert not _is_one2one([1, 1], [2, 3])

    # duplicate combos of [x,y] should not prevent is_one2one from returning True
    assert _is_one2one([1, 1, 2], [2, 2, 4])


def test_is_one2one_with_null():
    xs1 = [1,None,None]
    xs2 = [1,3,None]
    ys1 = [2,4,8]
    ys2 = [2,0,0]
    for dropna in (True, False):
        assert _is_one2one(xs1, ys2, dropna=dropna)
        assert _is_one2one(xs2, ys1, dropna=dropna)

        assert dropna == _is_one2one(xs1, ys1, dropna=dropna)
        assert dropna == _is_one2one(xs2, ys2, dropna=dropna)


@pytest.mark.network
def test_compound2cas():
    # alpha-pinene
    # https://pubchem.ncbi.nlm.nih.gov/compound/6654
    cid = 6654
    compound = pcp.Compound.from_cid(cid)
    # first one i can see under "Synonyms" section on page above is "80-56-8", and this
    # one also seems to list the greatest number of sources under "Other Identifiers" >
    # "CAS" section, though there are also:
    # - 25766-18-1 (CAS Common Chemistry; ChemIDplus)
    # - 2437-95-8 (ECHA)
    # - 7785-70-8 (EPA)
    cas = compound2cas(compound)
    assert cas == '80-56-8'


def test_allowed_kwarg():
    # compound2cas DOES have a verbose=True kwarg (specified as keyword only, via
    # preceding `*, ` in def
    assert allowed_kwarg(compound2cas, 'verbose')

    # compound2name takes no arguments other than a single required positional arg
    assert not allowed_kwarg(cu.compound2name, 'verbose')

    def f1(x, verbose=True):
        return x

    def f2(x, **kwargs):
        return x

    # despite this being as positional argument, can still pass as f3(verbose=True)
    def f3(verbose):
        pass

    # here verbose is positional only, so can not specify it as kwarg
    def f4(verbose, /):
        pass

    assert allowed_kwarg(f1, 'verbose')
    assert allowed_kwarg(f2, 'verbose')
    assert allowed_kwarg(f3, 'verbose')
    assert not allowed_kwarg(f4, 'verbose')



@pytest.mark.parametrize('value,expected', [
    (60000.0, 60000.0),
    ('60000.', 60000.0),
    ('6.0×10+8', 6.0e8),
    ('1.2×10−3', 1.2e-3),
    ('>6.0×10+8', np.nan),
    ('<4.0×10+11', np.nan),
])
def test_parse_nist_number(value, expected):
    assert cu._parse_nist_number(value) == pytest.approx(expected, nan_ok=True)


@pytest.mark.parametrize('value,expected', [
    ('1.9', 1.9),
    (0.2, 0.2),
    ('1.2×10−3', 1.2e-3),
    ('9.9×1023', 9.9e23),
])
def test_parse_henrys_law_org_number(value, expected):
    assert cu._parse_henrys_law_org_number(value) == pytest.approx(expected)


def test_inchi2est_k_henry(monkeypatch):
    monkeypatch.setattr(cu, 'inchi2water_solubility', lambda inchi, **kwargs:
        cu._property_series('water_solubility', 0.5, 'x', 'g / L')
    )
    monkeypatch.setattr(cu, 'inchi2vapor_pressure', lambda inchi, **kwargs:
        cu._property_series('vapor_pressure', 2.0, 'y', 'kPa', measured=True)
    )
    monkeypatch.setattr(cu, 'inchi2molar_mass', lambda inchi, **kwargs:
        cu._property_series('molar_mass', 100.0, 'z', 'g/mol')
    )
    ret = cu.inchi2est_k_henry('1S/CH4/h1H4')
    # 0.5 g/L / (100 g/mol) = 5 mol/m^3, / 2000 Pa
    assert ret['k_henry'] == pytest.approx(5 / 2000)
    assert ret['k_henry_units'] == 'mol / m^3 / Pa'
    assert 'VP/AS' in ret['k_henry_sources']

    monkeypatch.setattr(cu, 'inchi2vapor_pressure', lambda inchi, **kwargs:
        cu._property_series('vapor_pressure', measured=None)
    )
    assert cu.inchi2est_k_henry('1S/CH4/h1H4').isnull().all()


@pytest.mark.parametrize('measured', [True, False, None])
def test_pfo_vapor_pressure_measured(monkeypatch, measured):
    # (PubChem value, as there is no thermo Chemical)
    values = {'vapor_pressure_measured': measured}
    if measured is None:
        # (cached before this key was added)
        values = {}
    monkeypatch.setattr(cu, 'inchi2vapor_pressure', lambda inchi, **kwargs:
        cu._property_series('vapor_pressure', 2.0, 'url', 'kPa', **{
            k.replace('vapor_pressure_', ''): v for k, v in values.items()
        })
    )
    psat, source, ret_measured = cu._pfo_vapor_pressure('1S/CH4/h1H4', None, 298.15, None)
    assert psat == pytest.approx(2000)
    assert source.startswith('PubChem')
    assert ret_measured == measured


@pytest.mark.network
def test_compound2title():
    compound = pcp.Compound.from_cid(6549)
    assert cu.compound2title(compound) == 'Linalool'
