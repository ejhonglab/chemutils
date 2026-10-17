#!/usr/bin/env python3
"""
Checks lookups (against the real websites) return correct values, and that the values
are then served from the cache.

Run only these with `pytest -m network`, or skip them with `pytest -m "not network"`.
"""

from chemspipy.errors import ChemSpiPyHTTPError
import pytest

import chemutils as cu


pytestmark = pytest.mark.network

# Ethyl acetate: https://pubchem.ncbi.nlm.nih.gov/compound/8857
cid = 8857
inchi = '1S/C4H8O2/c1-3-6-4(2)5/h3H2,1-2H3'
inchikey = 'XEKOWRVHYACXOJ-UHFFFAOYSA-N'
cas = '141-78-6'


def check_second_call_cached(network_guard, fn, *args, **kwargs):
    """Returns `fn(*args, **kwargs)` after resetting cache state and blocking network.
    """
    cu._reset_cache_state()
    network_guard.blocking = True
    value = fn(*args, **kwargs)
    assert network_guard.attempts == []
    return value


@pytest.mark.parametrize('from_type,chem_id,second_chem_id,to_type,expected', [
    ('name', 'ethyl acetate', 'Ethyl Acetate (CAS)', 'cid', cid),
    ('name', 'Ethyl Acetate', 'ethyl acetate', 'inchi', inchi),
    ('name', 'ethyl acetate', ' ETHYL ACETATE', 'inchikey', inchikey),
    ('name', 'ethyl acetate', 'ethyl acetate', 'cas', cas),
    ('cas', cas, f' "{cas}"', 'cid', cid),
    ('cas', cas, cas, 'inchi', inchi),
    ('inchi', f'InChI={inchi}', inchi, 'cid', cid),
    ('inchi', inchi, f'InChI=InChI={inchi}', 'inchikey', inchikey),
    ('inchikey', inchikey, inchikey.lower(), 'cid', cid),
    ('inchikey', inchikey, f'InChIKey={inchikey}', 'inchi', inchi),
    ('cid', cid, float(cid), 'name', 'ethyl acetate'),
    ('cid', float(cid), str(cid), 'inchi', inchi),
    ('cid', cid, cid, 'cas', cas),
    ('cid', cid, cid, 'smiles', 'CCOC(=O)C'),
    # (E)-2-hexenal, which has stereochemistry
    ('name', '(E)-2-hexenal', 'E2-hexenal', 'cid', 5281168),
    ('name', '(E)-2-hexenal', 'E2-hexenal', 'smiles', 'CCCC=CC=O'),
    ('name', '(E)-2-hexenal', 'E2-hexenal', 'isomeric_smiles', 'CCC/C=C/C=O'),
    ('cid', 5281168, 5281168.0, 'inchi',
        '1S/C6H10O/c1-2-3-4-5-6-7/h4-6H,2-3H2,1H3/b5-4+'
    ),
])
def test_convert(network_guard, from_type, chem_id, second_chem_id, to_type, expected):
    value = cu.convert(chem_id, from_type=from_type, to_type=to_type)
    assert value == expected
    if to_type == 'cid':
        assert type(value) is int

    assert check_second_call_cached(network_guard, cu.convert, second_chem_id,
        from_type=from_type, to_type=to_type
    ) == expected


def in_range(prop, units, low, high):
    def check(value):
        assert value[f'{prop}_units'] == units
        assert low < value[prop] < high, f'{prop}={value[prop]} not in ({low}, {high})'
    return check


@pytest.mark.parametrize('fn_name,chem_id,second_chem_id,check', [
    ('cid2molar_mass', cid, float(cid), in_range('molar_mass', 'g/mol', 88.0, 88.2)),
    ('inchi2molar_mass', inchi, f'InChI={inchi}',
        in_range('molar_mass', 'g/mol', 88.0, 88.2)
    ),
    # Literature values below (for ~25 C) from PubChem / NIST / henrys-law.org pages
    # for ethyl acetate.
    # ~12.4 kPa (93.2 mmHg)
    ('inchi2vapor_pressure', f'InChI={inchi}', inchi,
        in_range('vapor_pressure', 'kPa', 8, 16)
    ),
    # ~80 g/L
    ('inchi2water_solubility', inchi, f'InChI={inchi}',
        in_range('water_solubility', 'g / L', 50, 120)
    ),
    # ~0.90 g/mL
    pytest.param('inchi2density', inchi, f'InChI={inchi}',
        in_range('density', 'g / mL', 0.85, 0.95),
        # pubchemprops_lookup gives 0.9003 g/mL with skip_unit_parse_err=True, but
        # inchi2density does not pass that, so it raises.
        marks=pytest.mark.xfail(raises=ValueError, strict=True, reason='parsing of '
            'PubChem string "CONVERSION FACTORS: 1 MG/L= 278 PPM" raises ValueError'
        )
    ),
    # ~6-9 mol/(kg*bar)
    ('inchi2nist_k_henry', inchi, f'InChI={inchi}',
        in_range('k_henry', 'mol / kg / bar', 3, 15)
    ),
    # ~6e-2 mol/(m^3*Pa) (same as ~6 mol/(kg*bar))
    ('inchi2maxplanck_k_henry', inchi, f'InChI={inchi}',
        in_range('k_henry', 'mol / m^3 / Pa', 3e-2, 0.15)
    ),
])
def test_cached_lookup(network_guard, fn_name, chem_id, second_chem_id, check):
    fn = getattr(cu, fn_name)
    value = fn(chem_id)
    check(value)

    assert check_second_call_cached(network_guard, fn, second_chem_id).equals(value)


@pytest.mark.skipif(cu.chemspipy_api_key is None, reason='no ChemSpider API key set')
def test_inchi2chemspider_id(network_guard):
    # http://www.chemspider.com/Chemical-Structure.8525.html
    try:
        csid = cu.inchi2chemspider_id(f'InChI={inchi}')
    except ChemSpiPyHTTPError as err:
        pytest.skip(f'ChemSpider API request failed (check API key): {err!r}')

    assert csid == 8525

    assert check_second_call_cached(network_guard, cu.inchi2chemspider_id, inchi
        ) == csid
