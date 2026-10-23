#!/usr/bin/env python3

import numpy as np
import pandas as pd
import pytest

import chemutils as cu
from test_cache import ethyl_acetate


# Including inputs that (some versions of) normalize_name handled incorrectly.
name2normalized = {
    'ethyl acetate': 'ethyl acetate',
    '  Ethyl  Acetate (CAS) ': 'ethyl acetate',
    'E2-hexenal': '(E)-2-hexenal',
    'e2-hexenal': '(E)-2-hexenal',
    '(E)-2-hexenal': '(E)-2-hexenal',
    '(e)-2-hexenal': '(E)-2-hexenal',
    'Z3-hexenol': '(Z)-3-hexenol',
    'E3-hexenol': '(E)-hex-3-en-1-ol',
    '2,4-hexadienal, (E,E)-': '(E,E)-2,4-hexadienal',
    '(E,E)-2,4-hexadienal': '(E,E)-2,4-hexadienal',
    '(e,e)-2,4-hexadienal': '(E,E)-2,4-hexadienal',
    'ethyl (E)-3-hexenoate': 'ethyl (E)-3-hexenoate',
    'a -pinene': 'alpha-pinene',
    'g-hexalactone': 'gamma-hexalactone',
    '1-octen- 3-ol': '1-octen-3-ol',
    'linalool oxide': 'trans-linalool oxide',
    '4-ethyl guaiacol': '4-ethylguaiacol',
    'docosane (branched)': 'docosane',
    # Old normalize_name joined these words, changing which compound PubChem returned
    'ethyl 3-hydroxyhexanoate': 'ethyl 3-hydroxyhexanoate',
    'ethyl (methylthio)acetate': 'ethyl (methylthio)acetate',
    'hexyl 2-methylbutanoate': 'hexyl 2-methylbutanoate',
    'methyl 1-propenyl trisulfide': 'methyl 1-propenyl trisulfide',
    # Old normalize_name was not idempotent for these
    'x (abcd) (efgh)': 'x',
    'ethyl ester hexanoic acid': 'hexanoic acid, ethyl ester',
    # Prefix correction needed after ester/acid reordering
    'ethyl ester g-hydroxybutyric acid': 'gamma-hydroxybutyric acid, ethyl ester',
    '(1S)-(+)-3-carene': '(1s)-(+)-3-carene',
    # Not stripping the only parenthetical (which used to leave nothing)
    '(trimethylsilyl)': '(trimethylsilyl)',
    ' (Trimethylsilyl) (CAS)': '(trimethylsilyl)',
}

@pytest.mark.parametrize('name,expected', name2normalized.items())
def test_normalize_name(name, expected):
    assert cu.normalize_name(name) == expected


def hardcoded_names():
    names = set(cu.hardcoded_type2null_keys['name'])
    for f2t_dict in cu.hardcoded['name'].values():
        names.update(f2t_dict)
    return sorted(names)


@pytest.mark.parametrize('name',
    sorted(set(name2normalized) | set(name2normalized.values())) + hardcoded_names()
)
def test_normalize_name_idempotent(name):
    normed = cu.normalize_name(name)
    assert cu.normalize_name(normed) == normed


inchi = ethyl_acetate['inchi']
inchikey = ethyl_acetate['inchikey']

# chem_id_type -> list of (raw input, expected normalized)
type2cases = {
    'cas': [
        ('141-78-6', '141-78-6'),
        (' "141-78-6" ', '141-78-6'),
        ('0-00-0', None),
    ],
    'inchi': [
        (inchi, inchi),
        (f'InChI={inchi}', inchi),
        (f' InChI=InChI={inchi} ', inchi),
        (f'inchi={inchi}', inchi),
    ],
    'inchikey': [
        (inchikey, inchikey),
        (f' {inchikey.lower()} ', inchikey),
        (f'InChIKey={inchikey}', inchikey),
    ],
    'cid': [
        (8857, 8857),
        (8857.0, 8857),
        (np.int64(8857), 8857),
        (np.float64(8857.0), 8857),
        (' 8857 ', 8857),
    ],
    'smiles': [
        (' CCOC(=O)C ', 'CCOC(=O)C'),
    ],
}

@pytest.mark.parametrize('chem_id_type,raw,expected',
    [(t, r, e) for t, cases in type2cases.items() for r, e in cases]
)
def test_normalize_ids(chem_id_type, raw, expected):
    normed = cu.normalize_id(chem_id_type, raw)
    assert normed == expected
    if chem_id_type == 'cid':
        assert type(normed) is int

    # Idempotent
    assert cu.normalize_id(chem_id_type, normed) == normed


@pytest.mark.parametrize('chem_id_type', ['name', 'cas', 'inchi', 'inchikey', 'cid',
    'smiles'
])
@pytest.mark.parametrize('null', [None, np.nan])
def test_normalize_null(chem_id_type, null):
    assert cu.isnull(cu.normalize_id(chem_id_type, null))


def test_normalize_cid_noninteger():
    with pytest.raises(ValueError):
        cu.normalize_cid(8857.5)


def test_normalize_id_equivalent_col_names():
    assert cu.normalize_id('cas_number', ' "141-78-6" ') == '141-78-6'
    assert cu.normalize_id('odor', 'E2-hexenal') == '(E)-2-hexenal'
    # Types without a normalize_<type> fn are returned unchanged
    assert cu.normalize_id('not_a_type', ' x ') == ' x '


def test_is_series():
    assert cu.is_series(pd.Series({'a': 1}))
    assert not cu.is_series(pd.DataFrame({'a': [1]}))
    assert not cu.is_series(None)

    assert cu.isnull(pd.Series({'a': None, 'b': None}))
    assert not cu.isnull(pd.Series({'a': None, 'b': 1.0}))
    assert cu.isnull(None) and cu.isnull(np.nan)
    assert not cu.isnull('x')
    # (arrays are not null, rather than raising ValueError)
    assert not cu.isnull(np.array([np.nan, np.nan]))


@pytest.mark.parametrize('name', ['', '  ', '(CAS)'])
def test_normalize_name_empty(name):
    assert cu.normalize_name(name) is None
