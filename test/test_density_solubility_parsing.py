#!/usr/bin/env python3
"""
Tests parsing of the strings in PubChem's 'Density' and 'Solubility' sections (as
`inchi2density` and `inchi2water_solubility` use them), with examples from PubChem pages.
"""

import pubchemprops.pubchemprops as pp
import pytest

import chemutils as cu


def _in_units(parsed, units):
    return cu.ureg.Quantity(parsed['quantity'], parsed['units']).to(units).magnitude


@pytest.mark.parametrize('string,expected_g_ml', [
    # Relative densities (to water), as most of PubChem's are
    ('Relative density (water = 1): 0.87', 0.87),
    ('0.876 at 68 °F (USCG, 1999) - Less dense than water; will float', 0.876),
    ('1.174 (USCG, 1999) - Denser than water; will sink', 1.174),
    ('0.876 at 15 °C/4 °C', 0.876),
    ('0.8053 @ 20 °C/4 °C', 0.8053),
    ('0.868-0.878', 0.873),
    ('Sp gr: 0.884-0.891', 0.8875),
    ('(57 °F): 0.81', 0.81),
    ('Density: 0.870-0.890 at 15 °C', 0.88),
    ('0.805 at 20 °C/4 °C; 0.7997 at 25 °C/4 °C; 0.8255 at 0 °C/4 °C', 0.7997),
    ('Density: 1.266 at 16.60 °C (solid); 1.053 at 16.67 °C (liquid); contracts '
        'slightly on freezing', 1.053
    ),
    # With units
    ('0.9003 g/cu cm at 20 °C', 0.9003),
    ('1.0328 g.cu cm at 20 °C', 1.0328),
    ('1.02 g/cm³', 1.02),
    ('1.02 g/cm3', 1.02),
    # Not of the compound as a liquid (or solid)
    ('DENSITY OF SATURATED AIR (AIR= 1) 1.02; CONVERSION FACTORS: 1 MG/L= 278 PPM; '
        '1 PPM= 3.60 MG/CU M', None
    ),
    ('Relative density of the vapour/air-mixture at 20 °C (air = 1): 1', None),
    ('Density (gas): 1.00045 g/l', None),
    ('1.639 g/L', None),
    ('1.27(relative gas density)', None),
    ('Bulk density: 6.71 lb/gal at 20 °C', None),
    ('1.05 at 59 °F for 10.17% weight/weight solution (EPA, 1998) - Denser than water; '
        'will sink', None
    ),
    ('Coefficient of cubical expansion: 9.56X10-7 at 0-94 °C; critical density: 0.2734',
        None
    ),
    ('< 0.810', None),
    ('DENSITY OF NATURAL ESTER IS ABOUT 1.180', None),
])
def test_parse_density(string, expected_g_ml):
    parsed = cu.parse_pubchemprops_string(string, expected_units='density')
    if expected_g_ml is None:
        assert parsed['quantity'] is None
    else:
        assert _in_units(parsed, 'g / mL') == pytest.approx(expected_g_ml)


@pytest.mark.parametrize('string,expected_g_l', [
    ('In water, 8.0X10+4 mg/L at 25 °C', 80),
    ('Solubility in water, g/100ml at 20 °C: 8.7 (poor)', 87),
    ('2.5 g in 100 mL water', 25),
    ('1.3 G SOL IN 100 ML WATER', 13),
    ('1 gram dissolves in 1000 mL water, 80 mL boiling water, 0.8 mL alcohol', 1),
    # (would need the density)
    ('1 ml in 40 ml water', None),
    # In other solvents
    ('Solubility in methanol: 54.6 g/100 g solution at -10 °C', None),
    ('1 ml in 5 ml 50% alcohol (in ethanol)', None),
])
def test_parse_water_solubility(string, expected_g_l):
    parsed = cu.parse_pubchemprops_string(string, expected_units='concentration')
    if expected_g_l is None:
        assert parsed['quantity'] is None
    else:
        assert _in_units(parsed, 'g / L') == pytest.approx(expected_g_l)


def _pubchem_values(monkeypatch, heading, values):
    """Has PubChem lookups give `values` (str for strings, else as given) for `heading`.
    """
    results = [
        {'Value': {'StringWithMarkup': [{'String': v}]} if isinstance(v, str) else v}
        for v in values
    ]
    monkeypatch.setattr(pp, 'get_second_layer_props',
        lambda name, props, cid=False: {heading: results}
    )
    monkeypatch.setattr(cu, 'convert', lambda *args, **kwargs: 1)


def test_inchi2density(monkeypatch):
    # Ethyl acetate's (https://pubchem.ncbi.nlm.nih.gov/compound/8857), the third of
    # which used to raise ValueError, and a table (as hydrogen chloride's page has)
    _pubchem_values(monkeypatch, 'Density', [
        '0.902 at 68 °F (USCG, 1999) - Less dense than water; will float',
        '0.9003 g/cu cm at 20 °C',
        'DENSITY OF SATURATED AIR (AIR= 1) 1.02; CONVERSION FACTORS: 1 MG/L= 278 PPM; '
            '1 PPM= 3.60 MG/CU M',
        'Relative density (water = 1): 0.9', '0.894-0.898', '0.9', '0.9003 @ 20°C',
        '0.90', {'Binary': ['eA=='], 'MimeType': 'text/csv'},
    ])
    ret = cu.inchi2density('1S/CH4/h1H4')
    # (the mean of those not from ranges)
    assert ret['density'] == pytest.approx((0.902 + 0.9003 + 0.9 + 0.9 + 0.9003 + 0.9) / 6)
    assert ret['density_units'] == 'g / mL'


def test_inchi2water_solubility(monkeypatch):
    _pubchem_values(monkeypatch, 'Solubility', [
        'In water, 8.0X10+4 mg/L at 25 °C',
        {'Number': [100], 'Unit': 'g/l at 20 °C (68 °F)'},
        {'Number': [8.7], 'Unit': 'g/100 mL'},
        # (not parseable, so not used)
        '1 g /o-methoxyphenol/ is soluble in mL water, in 1 mL glycerol',
        'Miscible with ethanol, ethyl ether; very soluble in acetone, benzene',
    ])
    ret = cu.inchi2water_solubility('1S/CH4/h1H4')
    assert ret['water_solubility'] == pytest.approx((80 + 100 + 87) / 3)
    assert ret['water_solubility_units'] == 'g / L'
