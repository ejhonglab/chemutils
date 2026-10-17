#!/usr/bin/env python3

import importlib.util
from pathlib import Path
import pickle

import pandas as pd
import pytest

import chemutils as cu
from test_cache import db_entries, ethyl_acetate, make_compound, new_process, pubchem


script = Path(__file__).resolve().parents[1] / 'scripts' / 'migrate_cache.py'
spec = importlib.util.spec_from_file_location('migrate_cache', script)
migrate_cache = importlib.util.module_from_spec(spec)
spec.loader.exec_module(migrate_cache)

inchi = ethyl_acetate['inchi']
inchikey = ethyl_acetate['inchikey']

vapor_pressure = pd.Series({'vapor_pressure': 12.4, 'vapor_pressure_sources': 'manual',
    'vapor_pressure_units': 'kPa'
})


def make_old_cache():
    # Same structure as the old pickle cache
    c = {ft: {tt: dict() for tt in cu.chem_id_types + cu.properties}
        for ft in cu.chem_id_types
    }
    c['name']['cid'].update({
        # Old normalize_name joined these words. Lookup with joined name failed.
        'hexyl 2-methylbutanoate': 24838,
        'hexyl2-methylbutanoate': None,
        # ...and here joined name gave a different (wrong) compound.
        'ethyl 3-hydroxyhexanoate': 61293,
        'ethyl3-hydroxyhexanoate': 736357,
        # Now normalize to same name
        '(E)-2-hexenal': 5281168,
        '(e)-2-hexenal': 5281168,
        'docosane (branched)': None,
        'docosane': 12405,
        # Hardcoded null key
        'odor': None,
    })
    c['name']['inchi'].update({
        'hexyl 2-methylbutanoate': 'hexyl-2-methylbutanoate-inchi',
        'hexyl2-methylbutanoate': None,
        'ethyl 3-hydroxyhexanoate': 'ethyl-3-hydroxyhexanoate-inchi',
        'ethyl3-hydroxyhexanoate': 'wrong-inchi',
        'docosane (branched)': None,
        'docosane': 'docosane-inchi',
    })
    # Cached when CID lookup of 'docosane (branched)' failed. Now that it normalizes to
    # 'docosane' (which has a CID), it should be dropped.
    c['name']['inchikey']['docosane (branched)'] = None
    # Same as hardcoded value
    c['name']['cas']['a -pinene'] = '102640-64-2'

    c['inchi']['cid']['InChI=' + inchi] = 8857.0
    c['inchi']['vapor_pressure'][inchi] = vapor_pressure
    c['inchi']['vapor_pressure'][None] = pd.Series({'vapor_pressure': None,
        'vapor_pressure_sources': None, 'vapor_pressure_units': None
    })
    c['inchikey']['cid'][inchikey + ' '] = 8857
    c['cid']['compound'] = {8857.0: make_compound(ethyl_acetate)}
    c['cid']['name'][8857] = 'ethyl acetate'

    is_manual = {('inchi', 'vapor_pressure', inchi): True}
    return c, is_manual


@pytest.fixture
def old_pickle(tmp_path):
    path = tmp_path / 'old_cache.p'
    old_cache, is_manual = make_old_cache()
    with open(path, 'wb') as f:
        pickle.dump({'cache': old_cache, 'is_manual': is_manual}, f)
    return path


def test_migrate(old_pickle, pubchem):
    pickle_bytes = old_pickle.read_bytes()

    _, report = migrate_cache.migrate(old_pickle, cu.cache_path(), verbose=False)

    # Old cache should not be modified.
    assert old_pickle.read_bytes() == pickle_bytes

    entries = db_entries()
    for (ft, tt, key), (value, manual) in entries.items():
        assert cu.normalize_id(ft, key) == key

    assert report['wrong_cid_dropped'] == [
        ('ethyl 3-hydroxyhexanoate', 61293, 'ethyl3-hydroxyhexanoate', 736357)
    ]
    assert report['conflicts'] == []

    name_keys = {k for ft, tt, k in entries if ft == 'name'}
    # Incorrectly normalized keys, with one-sided null or wrong CID, are removed.
    assert 'hexyl2-methylbutanoate' not in name_keys
    assert 'ethyl3-hydroxyhexanoate' not in name_keys
    # Keys normalizing to the same thing are merged.
    assert '(e)-2-hexenal' not in name_keys
    assert 'docosane (branched)' not in name_keys
    # Stale null from failed CID lookup, now that name has a CID.
    assert ('name', 'inchikey', 'docosane') not in entries
    # Hardcoded values (which are applied separately) are not stored.
    assert 'odor' not in name_keys
    assert 'a -pinene' not in name_keys and 'alpha-pinene' not in name_keys

    assert entries[('inchi', 'vapor_pressure', inchi)][1] == 1
    assert not any(k is None for _, _, k in entries)

    pubchem.fail = True
    new_process()

    name2cid = {
        'hexyl 2-methylbutanoate': 24838,
        'ethyl 3-hydroxyhexanoate': 61293,
        '(E)-2-hexenal': 5281168,
        '(e)-2-hexenal': 5281168,
        'E2-Hexenal': 5281168,
        'docosane': 12405,
        'docosane (branched)': 12405,
        'odor': None,
    }
    for name, cid in name2cid.items():
        assert cu.convert(name, from_type='name', to_type='cid', allow_nan=True) == cid

    assert cu.convert('docosane (branched)', from_type='name', to_type='inchi'
        ) == 'docosane-inchi'
    assert cu.convert('a -pinene', from_type='name', to_type='cas') == '102640-64-2'

    for chem_id in ('InChI=' + inchi, inchi):
        cid = cu.convert(chem_id, from_type='inchi', to_type='cid')
        assert cid == 8857 and type(cid) is int
        assert cu.inchi2vapor_pressure(chem_id).equals(vapor_pressure)

    assert cu.convert(inchikey.lower(), from_type='inchikey', to_type='cid') == 8857

    # Compound from cache
    assert cu.convert(8857, from_type='cid', to_type='inchikey') == inchikey
    assert entries[('cid', 'compound', 8857)][0] == make_compound(ethyl_acetate)


def test_migrate_dry_run(old_pickle):
    resolved, report = migrate_cache.migrate(old_pickle, cu.cache_path(),
        dry_run=True, verbose=False
    )
    assert len(resolved) > 0
    assert not cu.cache_path().exists()


def test_migrate_existing(old_pickle):
    cu.cache['name']['cid']['x'] = 1
    new_process()

    with pytest.raises(RuntimeError, match='--force'):
        migrate_cache.migrate(old_pickle, cu.cache_path(), verbose=False)
    assert ('name', 'cid', 'x') in db_entries()

    migrate_cache.migrate(old_pickle, cu.cache_path(), force=True, verbose=False)
    assert ('name', 'cid', 'x') not in db_entries()
