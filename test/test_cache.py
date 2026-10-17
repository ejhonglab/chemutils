#!/usr/bin/env python3
"""
Tests of chemutils caching, against a fake (offline) PubChem.

Each test uses a new, empty cache (see `chemutils_test_cache` in ../conftest.py).
`cu._reset_cache_state()` is used to simulate a new process, so later lookups can only
come from the saved cache.
"""

import sqlite3

import numpy as np
import pandas as pd
import pubchempy as pcp
import pytest

import chemutils as cu


ethyl_acetate = dict(
    cid=8857,
    name='ethyl acetate',
    iupac_name='ethyl acetate',
    inchi='1S/C4H8O2/c1-3-6-4(2)5/h3H2,1-2H3',
    inchikey='XEKOWRVHYACXOJ-UHFFFAOYSA-N',
    cas='141-78-6',
    smiles='CCOC(=O)C',
    isomeric_smiles='CCOC(=O)C',
    molecular_weight='88.11',
)
hexenal = dict(
    cid=5281168,
    name='(E)-2-hexenal',
    iupac_name='(E)-hex-2-enal',
    inchi='1S/C6H10O/c1-2-3-4-5-6-7/h4-6H,2-3H2,1H3/b5-4+',
    inchikey='MBDOYVRWFFCFHM-SNAWJCMRSA-N',
    cas='6728-26-3',
    smiles='CCCC=CC=O',
    isomeric_smiles='CCC/C=C/C=O',
    molecular_weight='98.14',
)

# Names of PubChem SMILES properties for ('smiles', 'isomeric_smiles'), since 2025.
new_smiles_names = ('Connectivity', 'Absolute')

def make_compound(c, smiles_names=new_smiles_names) -> pcp.Compound:
    def prop(label, name, value):
        urn = {'label': label}
        if name is not None:
            urn['name'] = name
        return {'urn': urn, 'value': {'sval': value}}

    return pcp.Compound({
        'id': {'id': {'cid': c['cid']}},
        'atoms': {'aid': [], 'element': []},
        'props': [
            prop('IUPAC Name', 'Preferred', c['iupac_name']),
            prop('InChI', 'Standard', 'InChI=' + c['inchi']),
            prop('InChIKey', 'Standard', c['inchikey']),
            prop('SMILES', smiles_names[0], c['smiles']),
            prop('SMILES', smiles_names[1], c['isomeric_smiles']),
            prop('Molecular Weight', None, c['molecular_weight']),
        ],
    })


class FakePubChem:
    """Replaces the pubchempy fns chemutils uses, counting calls.

    Name searches are case-insensitive (as PubChem's seem to be), but otherwise exact.
    """
    def __init__(self, compounds):
        self.cid2compound = {c['cid']: c for c in compounds}
        self.names = dict()
        self.inchis = dict()
        self.inchikeys = dict()
        for c in compounds:
            self.names[c['name'].lower()] = c['cid']
            # (cas2cid also searches names)
            self.names[c['cas']] = c['cid']
            self.inchis['InChI=' + c['inchi']] = c['cid']
            self.inchikeys[c['inchikey']] = c['cid']

        self.calls = []
        # If True, any call raises (to check things come from the cache)
        self.fail = False

    def _call(self, *args):
        if self.fail:
            raise AssertionError(f'unexpected PubChem call: {args}')
        self.calls.append(args)

    def get_compounds(self, identifier, namespace='cid', **kwargs):
        self._call('get_compounds', identifier, namespace)
        namespace2ids = {
            'name': (self.names, identifier.lower()),
            'inchi': (self.inchis, identifier),
            'inchikey': (self.inchikeys, identifier),
        }
        ids, identifier = namespace2ids[namespace]
        if identifier not in ids:
            return []
        return [make_compound(self.cid2compound[ids[identifier]])]

    def from_cid(self, cid, **kwargs):
        self._call('from_cid', cid)
        # Should always be normalized to int before this.
        assert type(cid) is int
        return make_compound(self.cid2compound[cid])

    def synonyms(self, compound):
        self._call('synonyms', compound.cid)
        c = self.cid2compound[compound.cid]
        return [c['name'], c['cas']]

    def n_calls(self, fn_name=None):
        return len([c for c in self.calls if fn_name is None or c[0] == fn_name])


@pytest.fixture
def pubchem(monkeypatch):
    fake = FakePubChem([ethyl_acetate, hexenal])
    monkeypatch.setattr(pcp, 'get_compounds', fake.get_compounds)
    monkeypatch.setattr(pcp.Compound, 'from_cid', fake.from_cid)
    monkeypatch.setattr(pcp.Compound, 'synonyms', property(fake.synonyms))
    return fake


def db_rows(from_type, to_type):
    """Returns {key: decoded value} saved in cache database (bypassing chemutils).
    """
    if not cu.cache_path().exists():
        return dict()

    conn = sqlite3.connect(cu.cache_path())
    rows = conn.execute('SELECT key, value FROM entries WHERE from_type = ? AND '
        'to_type = ?', (from_type, to_type)
    ).fetchall()
    conn.close()
    return {k: cu._decode_value(v) for k, v in rows}


def new_process():
    cu._reset_cache_state()


# (from_type, to_type, list of raw (non-normalized) inputs, normalized input, expected,
# value to seed cache with)
convert_cases = [
    ('name', 'inchi', ['  Ethyl  Acetate (CAS) ', 'ETHYL ACETATE'], 'ethyl acetate',
        ethyl_acetate['inchi'], 'seeded-inchi'
    ),
    ('name', 'inchi', ['E2-hexenal', '(e)-2-hexenal'], '(E)-2-hexenal',
        hexenal['inchi'], 'seeded-inchi'
    ),
    ('cas', 'inchi', [' "141-78-6" '], '141-78-6', ethyl_acetate['inchi'],
        'seeded-inchi'
    ),
    ('inchi', 'cid', ['InChI=' + ethyl_acetate['inchi'],
        ' InChI=InChI=' + ethyl_acetate['inchi']], ethyl_acetate['inchi'], 8857, 999999
    ),
    ('inchikey', 'cid', [' ' + ethyl_acetate['inchikey'].lower(),
        'InChIKey=' + ethyl_acetate['inchikey']], ethyl_acetate['inchikey'], 8857, 999999
    ),
    ('cid', 'inchi', [8857.0, '8857', np.int64(8857)], 8857, ethyl_acetate['inchi'],
        'seeded-inchi'
    ),
    ('cid', 'cas', [8857.0], 8857, ethyl_acetate['cas'], 'seeded-cas'),
    ('name', 'cid', ['Ethyl Acetate'], 'ethyl acetate', 8857, 999999),
]
convert_params = [
    pytest.param(ft, tt, raw, normed, expected, seed_value,
        id=f'{ft}->{tt}:{raw!r}'
    )
    for ft, tt, raws, normed, expected, seed_value in convert_cases for raw in raws
]

@pytest.mark.parametrize('from_type,to_type,raw,normed,expected,seed_value',
    convert_params
)
@pytest.mark.parametrize('seed', [None, 'raw', 'normed'])
@pytest.mark.parametrize('first', ['raw', 'normed'])
@pytest.mark.parametrize('second', ['raw', 'normed'])
def test_convert_cache(pubchem, from_type, to_type, raw, normed, expected, seed_value,
    seed, first, second):
    """Checks raw and normalized IDs share one cache entry, whichever is used to look
    up / fill it, and that the cache is used whenever the entry exists.
    """
    inputs = {'raw': raw, 'normed': normed}

    if seed is not None:
        # The cache view normalizes keys, so this should be equivalent to setting it via
        # the normalized ID.
        cu.cache[from_type][to_type][inputs[seed]] = seed_value
        new_process()
        assert db_rows(from_type, to_type) == {normed: seed_value}
        expected = seed_value

    assert cu.convert(inputs[first], from_type=from_type, to_type=to_type) == expected
    if seed is None:
        assert pubchem.n_calls() > 0
    else:
        assert pubchem.n_calls() == 0

    # Value should be saved under the normalized ID only.
    assert db_rows(from_type, to_type) == {normed: expected}
    if seed is None and from_type != 'cid':
        assert db_rows(from_type, 'cid') == {normed: ethyl_acetate['cid']
            if expected in (ethyl_acetate['inchi'], 8857) else hexenal['cid']
        }

    new_process()
    pubchem.fail = True
    assert cu.convert(inputs[second], from_type=from_type, to_type=to_type) == expected

    # Should also be used within the same process (from in-memory cache).
    assert cu.convert(inputs[second], from_type=from_type, to_type=to_type) == expected


def test_convert_uses_cached_intermediates(pubchem):
    raw = 'Ethyl Acetate'
    cu.cache['name']['cid'][raw] = 8857
    new_process()

    # CID from cache, but Compound needs lookup.
    assert cu.convert(raw, from_type='name', to_type='inchikey') == \
        ethyl_acetate['inchikey']
    assert pubchem.n_calls('get_compounds') == 0
    assert pubchem.n_calls('from_cid') == 1

    # Compound now also cached.
    new_process()
    pubchem.fail = True
    assert cu.convert('ethyl acetate', from_type='name', to_type='inchi') == \
        ethyl_acetate['inchi']
    assert cu.convert(8857.0, from_type='cid', to_type='name') == 'ethyl acetate'

    compound = db_rows('cid', 'compound')[8857]
    assert isinstance(compound, pcp.Compound)
    assert compound == make_compound(ethyl_acetate)


@pytest.mark.parametrize('smiles_names', [
    new_smiles_names,
    # Before 2025
    ('Canonical', 'Isomeric'),
])
def test_compound2smiles(smiles_names):
    compound = make_compound(hexenal, smiles_names=smiles_names)
    assert cu.compound2smiles(compound) == 'CCCC=CC=O'
    assert cu.compound2isomeric_smiles(compound) == 'CCC/C=C/C=O'


def test_compound2smiles_record_variants():
    compound = make_compound(hexenal, smiles_names=('Canonical', 'Isomeric'))
    # Some records from around the renaming had 'Absolute', but not 'Connectivity'.
    compound.record['props'].append({'urn': {'label': 'SMILES', 'name': 'Absolute'},
        'value': {'sval': 'from-absolute'}
    })
    assert cu.compound2smiles(compound) == 'CCCC=CC=O'
    assert cu.compound2isomeric_smiles(compound) == 'from-absolute'

    compound.record['props'] = [p for p in compound.record['props']
        if p['urn']['label'] != 'SMILES'
    ]
    assert cu.compound2smiles(compound) is None
    assert cu.compound2isomeric_smiles(compound) is None


@pytest.mark.parametrize('to_type,expected', [
    ('smiles', hexenal['smiles']),
    ('isomeric_smiles', hexenal['isomeric_smiles']),
])
def test_convert_smiles(pubchem, to_type, expected):
    assert cu.convert('E2-hexenal', from_type='name', to_type=to_type) == expected
    assert db_rows('name', to_type) == {'(E)-2-hexenal': expected}

    new_process()
    pubchem.fail = True
    assert cu.convert('(e)-2-hexenal', from_type='name', to_type=to_type) == expected
    assert cu.convert(hexenal['cid'], from_type='cid', to_type=to_type) == expected


def test_convert_failed_lookup_cached(pubchem):
    raw = ' Not A Real Chemical '
    assert cu.convert(raw, from_type='name', to_type='inchi', allow_nan=True) is None
    n_calls = pubchem.n_calls()
    assert n_calls > 0

    # Null values are saved for all ID types, under the normalized key.
    normed = 'not a real chemical'
    for tt in cu.chem_id_types:
        assert db_rows('name', tt) == {normed: None}

    new_process()
    pubchem.fail = True
    for chem_id in (raw, normed):
        for tt in ('inchi', 'cas'):
            assert cu.convert(chem_id, from_type='name', to_type=tt, allow_nan=True
                ) is None

    # ...unless ignore_cache='if_null'
    pubchem.fail = False
    assert cu.convert(normed, from_type='name', to_type='inchi', allow_nan=True,
        ignore_cache='if_null') is None
    assert pubchem.n_calls() > n_calls


def test_convert_try_non_normalized(pubchem):
    # Only findable by the name before normalization.
    raw = 'Weird Name (CAS)'
    pubchem.names[raw.lower()] = 8857
    normed = cu.normalize_name(raw)
    assert normed == 'weird name'

    assert cu.convert(raw, from_type='name', to_type='cid', try_non_normalized=False,
        allow_nan=True) is None
    cu.clear_cache('name', 'cid')

    assert cu.convert(raw, from_type='name', to_type='cid') == 8857
    names_searched = [c[1] for c in pubchem.calls if c[0] == 'get_compounds']
    assert names_searched[-2:] == [normed, raw]

    # Saved under normalized name.
    assert db_rows('name', 'cid') == {normed: 8857}
    new_process()
    pubchem.fail = True
    assert cu.convert(normed, from_type='name', to_type='cid') == 8857


def test_convert_series(pubchem):
    names = pd.Series(['Ethyl Acetate', 'ethyl acetate', ' ETHYL ACETATE (CAS)',
        'E2-hexenal', None], name='name'
    )
    inchis = cu.convert(names, to_type='inchi', dropna=False, allow_nan=True)
    assert list(inchis[:4]) == [ethyl_acetate['inchi']] * 3 + [hexenal['inchi']]
    assert pd.isnull(inchis[4])
    # One search per normalized name.
    assert pubchem.n_calls('get_compounds') == 2

    new_process()
    pubchem.fail = True
    assert cu.convert(names, to_type='inchi', dropna=False, allow_nan=True
        ).equals(inchis)


def test_convert_dataframe(pubchem):
    df = pd.DataFrame({'name': ['Ethyl Acetate', 'E2-hexenal'],
        'cas': ['141-78-6', None]
    })
    converted = cu.convert(df, to_type='inchi')
    assert list(converted.inchi) == [ethyl_acetate['inchi'], hexenal['inchi']]
    # Input ID columns are normalized in output.
    assert list(converted.name) == ['ethyl acetate', '(E)-2-hexenal']

    new_process()
    pubchem.fail = True
    assert cu.convert(df, to_type='inchi').equals(converted)


def test_hardcoded_values(pubchem):
    pubchem.fail = True

    # Hardcoded with raw key 'a -pinene' and 'E2-hexenal', but should apply to
    # anything normalizing to the same name.
    for name in ('a -pinene', 'alpha-pinene', 'Alpha-Pinene'):
        assert cu.convert(name, from_type='name', to_type='cas') == '102640-64-2'

    for name in ('E2-hexenal', '(E)-2-hexenal', '(e)-2-hexenal'):
        assert cu.convert(name, from_type='name', to_type='cas') == '1335-39-3'

    assert cu.convert('odor', from_type='name', to_type='inchi', allow_nan=True) is None

    # Hardcoded values take precedence over saved values, and still apply with
    # ignore_cache=True.
    cu.cache['name']['cas']['alpha-pinene'] = '80-56-8'
    assert cu.convert('a -pinene', from_type='name', to_type='cas') == '102640-64-2'
    assert cu.convert('a -pinene', from_type='name', to_type='cas', ignore_cache=True
        ) == '102640-64-2'

    # ...and are never written to the database.
    assert db_rows('name', 'cas') == {'alpha-pinene': '80-56-8'}
    assert db_rows('name', 'inchi') == dict()


def test_ignore_cache(pubchem):
    cu.cache['name']['inchi']['ethyl acetate'] = 'seeded-inchi'
    new_process()

    assert cu.convert('Ethyl Acetate', from_type='name', to_type='inchi',
        ignore_cache=True) == ethyl_acetate['inchi']
    assert pubchem.n_calls() > 0
    # Looked up value overwrites saved one.
    assert db_rows('name', 'inchi') == {'ethyl acetate': ethyl_acetate['inchi']}

    # Values looked up earlier in this process are still used with ignore_cache=True.
    pubchem.fail = True
    assert cu.convert('ethyl acetate', from_type='name', to_type='inchi',
        ignore_cache=True) == ethyl_acetate['inchi']


def test_clear_cache(pubchem):
    for name in ('ethyl acetate', 'E2-hexenal'):
        cu.convert(name, from_type='name', to_type='inchi')
    assert set(db_rows('name', 'inchi')) == {'ethyl acetate', '(E)-2-hexenal'}

    # Raw value should delete normalized key.
    cu.clear_cache('name', 'inchi', 'E2-Hexenal')
    assert set(db_rows('name', 'inchi')) == {'ethyl acetate'}

    cu.clear_cache('name', 'inchi')
    assert db_rows('name', 'inchi') == dict()
    assert len(db_rows('name', 'cid')) == 2

    # No args: saved values are ignored for the rest of this process, but not deleted,
    # and values looked up after are still saved and used.
    cu.cache['name']['cas']['ethyl acetate'] = 'seeded-cas'
    cu.clear_cache()
    n_calls = pubchem.n_calls()
    assert cu.convert('ethyl acetate', from_type='name', to_type='cas') == \
        ethyl_acetate['cas']
    assert pubchem.n_calls() > n_calls

    pubchem.fail = True
    assert cu.convert('ethyl acetate', from_type='name', to_type='cas') == \
        ethyl_acetate['cas']
    assert db_rows('name', 'cas') == {'ethyl acetate': ethyl_acetate['cas']}


def test_delete_cache(pubchem, chemutils_test_cache):
    cu.convert('ethyl acetate', from_type='name', to_type='inchi')
    assert chemutils_test_cache.exists()
    cu.delete_cache()
    assert not chemutils_test_cache.exists()

    pubchem.fail = True
    with pytest.raises(AssertionError, match='unexpected PubChem call'):
        cu.convert('ethyl acetate', from_type='name', to_type='inchi')


def test_cache_view():
    sub_cache = cu.cache['inchi']['cid']
    sub_cache['InChI=' + ethyl_acetate['inchi']] = 8857
    assert len(sub_cache) == 1
    assert list(sub_cache) == [ethyl_acetate['inchi']]
    assert sub_cache[ethyl_acetate['inchi']] == 8857
    assert 'InChI=' + ethyl_acetate['inchi'] in sub_cache
    assert 'not an inchi' not in sub_cache

    del sub_cache[' InChI=' + ethyl_acetate['inchi']]
    assert len(sub_cache) == 0
    with pytest.raises(KeyError):
        sub_cache[ethyl_acetate['inchi']]


def test_cache_values_roundtrip():
    values = [
        None,
        1,
        2.5,
        float('nan'),
        'x',
        pd.Series({'density': 0.9, 'density_sources': 'x', 'density_units': 'g / mL'}),
        pd.Series({'density': None, 'density_sources': None, 'density_units': None}),
        pd.Series({'a': 1.0, 'b': 2.0}),
        make_compound(ethyl_acetate),
        # (falls back to pickle)
        {'a': [1, 2]},
    ]
    for i, value in enumerate(values):
        cu.cache['inchi']['test'][str(i)] = value

    new_process()
    for i, value in enumerate(values):
        decoded = cu.cache['inchi']['test'][str(i)]
        assert type(decoded) is type(value)
        if cu.is_series(value):
            assert decoded.equals(value)
            assert decoded.dtype == value.dtype
        else:
            assert cu._values_equal(decoded, value)


def test_cache_persists_across_processes(pubchem, chemutils_test_cache):
    import subprocess
    import sys

    cu.convert('E2-hexenal', from_type='name', to_type='inchi')
    cu._reset_cache_state()

    code = ('import socket\n'
        'def blocked(*args, **kwargs): raise RuntimeError("network")\n'
        'socket.getaddrinfo = blocked\n'
        'import chemutils as cu\n'
        'print(cu.convert("(e)-2-hexenal", from_type="name", to_type="inchi"))\n'
    )
    out = subprocess.run([sys.executable, '-c', code], check=True,
        capture_output=True, text=True,
        env={**__import__('os').environ, 'CHEMUTILS_CACHE': str(chemutils_test_cache)}
    ).stdout.strip()
    assert out == hexenal['inchi']


# Using a fn defined in here (rather than one of the real lookup fns), to test cached
# behavior without needing to fake any websites.
fake_prop_calls = []

@cu.cached
def inchi2fake_prop(inchi, verbose=False):
    fake_prop_calls.append(inchi)
    if cu.isnull(inchi) or inchi == 'null':
        return pd.Series({'fake_prop': None, 'fake_prop_units': None})
    return pd.Series({'fake_prop': float(len(inchi)), 'fake_prop_units': 'kPa'})


@pytest.fixture
def fake_prop():
    fake_prop_calls.clear()
    return fake_prop_calls


def fake_prop_value(inchi):
    return pd.Series({'fake_prop': float(len(inchi)), 'fake_prop_units': 'kPa'})


@pytest.mark.parametrize('seed', [None, 'raw', 'normed'])
@pytest.mark.parametrize('first', ['raw', 'normed'])
@pytest.mark.parametrize('second', ['raw', 'normed'])
def test_cached(fake_prop, seed, first, second):
    normed = ethyl_acetate['inchi']
    inputs = {'raw': ' InChI=' + normed, 'normed': normed}

    expected = fake_prop_value(normed)
    if seed is not None:
        expected = pd.Series({'fake_prop': -1.0, 'fake_prop_units': 'kPa'})
        cu.cache['inchi']['fake_prop'][inputs[seed]] = expected
        new_process()

    assert inchi2fake_prop(inputs[first]).equals(expected)
    # Wrapped fn should only ever get normalized input.
    assert fake_prop == ([] if seed is not None else [normed])
    assert list(db_rows('inchi', 'fake_prop')) == [normed]

    new_process()
    assert inchi2fake_prop(inputs[second]).equals(expected)
    assert inchi2fake_prop(inputs[second], verbose=True).equals(expected)
    assert len(fake_prop) == (0 if seed is not None else 1)


def test_cached_null_input(fake_prop):
    for null in (None, np.nan):
        inchi2fake_prop(null)
    # Passed to fn without caching
    assert len(fake_prop) == 2
    assert db_rows('inchi', 'fake_prop') == dict()


def test_cached_ignore_cache(fake_prop):
    inchi = ethyl_acetate['inchi']
    cu.cache['inchi']['fake_prop'][inchi] = 'seeded'
    new_process()

    assert inchi2fake_prop(inchi, ignore_cache=True).equals(fake_prop_value(inchi))
    assert len(fake_prop) == 1
    # Value computed earlier in this process is still used with ignore_cache=True
    inchi2fake_prop(inchi, ignore_cache=True)
    assert len(fake_prop) == 1

    assert db_rows('inchi', 'fake_prop')[inchi].equals(fake_prop_value(inchi))


def test_cached_ignore_cache_null(fake_prop):
    inchi = ethyl_acetate['inchi']
    null_value = pd.Series({'fake_prop': None, 'fake_prop_units': None})
    cu.cache['inchi']['fake_prop'][inchi] = null_value
    new_process()

    assert inchi2fake_prop(inchi).equals(null_value)
    assert len(fake_prop) == 0

    # Recomputed, since cached value is null
    assert inchi2fake_prop(inchi, ignore_cache_null=True).equals(fake_prop_value(inchi))
    assert len(fake_prop) == 1

    # Non-null cached values are still used
    inchi2fake_prop(inchi, ignore_cache_null=True)
    assert len(fake_prop) == 1


def test_cached_leave_nonnull_cache(fake_prop):
    # Fake fn returns null for this input
    key = 'null'
    cu.cache['inchi']['fake_prop'][key] = fake_prop_value('x')
    new_process()

    ret = inchi2fake_prop(key, ignore_cache=True, leave_nonnull_cache=True)
    assert cu.isnull(ret)
    assert db_rows('inchi', 'fake_prop')[key].equals(fake_prop_value('x'))

    inchi2fake_prop(key, ignore_cache=True)
    assert cu.isnull(db_rows('inchi', 'fake_prop')[key])


def test_cached_checks_property_units():
    @cu.cached
    def inchi2vapor_pressure(inchi):
        return pd.Series({'vapor_pressure': 1.0, 'vapor_pressure_sources': 'x',
            'vapor_pressure_units': 'g / mL'
        })

    with pytest.raises(ValueError, match='units'):
        inchi2vapor_pressure(ethyl_acetate['inchi'])

    assert db_rows('inchi', 'vapor_pressure') == dict()


def property_value(prop, value, units):
    return pd.Series({prop: value, f'{prop}_sources': 'seeded', f'{prop}_units': units})


inchi = ethyl_acetate['inchi']
# (lookup fn name, from_type, raw input, normalized input, value to seed cache with)
cached_fn_cases = [
    ('cid2molar_mass', 'cid', 8857.0, 8857, property_value('molar_mass', 1, 'g/mol')),
    ('inchi2molar_mass', 'inchi', f'InChI={inchi}', inchi,
        property_value('molar_mass', 1, 'g/mol')
    ),
    ('inchi2water_solubility', 'inchi', f'InChI={inchi}', inchi,
        property_value('water_solubility', 1, 'g / L')
    ),
    ('inchi2vapor_pressure', 'inchi', f'InChI={inchi}', inchi,
        property_value('vapor_pressure', 1, 'kPa')
    ),
    ('inchi2density', 'inchi', f'InChI={inchi}', inchi,
        property_value('density', 1, 'g / mL')
    ),
    ('inchi2nist_k_henry', 'inchi', f'InChI={inchi}', inchi,
        property_value('k_henry', 1, 'mol / kg / bar')
    ),
    ('inchi2maxplanck_k_henry', 'inchi', f'InChI={inchi}', inchi,
        property_value('k_henry', 1, 'mol / m^3 / Pa')
    ),
    ('inchi2chemspider_id', 'inchi', f'InChI={inchi}', inchi, 12345),
]

@pytest.mark.parametrize('fn_name,from_type,raw,normed,seed_value', cached_fn_cases)
@pytest.mark.parametrize('seed', ['raw', 'normed'])
def test_lookup_fns_use_cache(pubchem, fn_name, from_type, raw, normed, seed_value,
    seed):
    """Checks each @cached lookup fn returns saved values, for raw / normalized input.

    Network access is blocked (see conftest.py) and fake PubChem fails on any call.
    See test_network_lookups.py for checks that the lookups compute correct values.
    """
    fn = getattr(cu, fn_name)
    to_type = fn_name.split('2', 1)[1]

    cu.cache[from_type][to_type][{'raw': raw, 'normed': normed}[seed]] = seed_value
    new_process()
    pubchem.fail = True

    for chem_id in (raw, normed):
        assert cu._values_equal(fn(chem_id), seed_value)


@pytest.mark.parametrize('fn_name,chem_id', [
    ('cid2molar_mass', 8857.0),
    ('inchi2molar_mass', f'InChI={inchi}'),
])
def test_molar_mass_lookups(pubchem, fn_name, chem_id):
    fn = getattr(cu, fn_name)
    value = fn(chem_id)
    assert value['molar_mass'] == 88.11
    assert value['molar_mass_units'] == 'g/mol'

    new_process()
    pubchem.fail = True
    assert fn(cu.normalize_id(fn_name.split('2')[0], chem_id)).equals(value)


def test_add_properties_to_cache():
    df = pd.DataFrame({
        'inchi': [f'InChI={inchi}'],
        'density': [900.0],
        'density_units': ['g / L'],
        'density_sources': ['manual'],
    })
    cu.add_properties_to_cache(df)

    conn = sqlite3.connect(cu.cache_path())
    rows = conn.execute("SELECT key, manual FROM entries WHERE to_type = 'density'"
        ).fetchall()
    conn.close()
    assert rows == [(inchi, 1)]

    new_process()
    value = cu.inchi2density(inchi)
    assert value['density'] == pytest.approx(0.9)
    assert value['density_units'] == 'g / mL'
