import base64
from collections import Counter
from collections.abc import MutableMapping
import contextvars
import functools
import inspect
from io import StringIO
import json
import logging
import os
from os.path import split, join, exists
from pathlib import Path
import pickle
from pprint import pformat
import re
import sqlite3
import time
from typing import Optional, Callable
import traceback
import urllib.error
import urllib.parse

from bs4 import BeautifulSoup
import numpy as np
import pubchempy as pcp
from pubchempy import Compound
import pandas as pd
import pint
# Should be included with pint. Some pint fn can throw this.
from tokenize import TokenError
from rdkit import Chem
from rdkit.Chem.inchi import InchiToInchiKey
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from chemspipy import ChemSpider
from chemspipy.search import Results


session = requests.Session()
# https://stackoverflow.com/questions/15431044
# (429 is what PubChem responds with when rate limiting)
retry = Retry(total=5, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504])
session.mount('http://', HTTPAdapter(max_retries=retry))
session.mount('https://', HTTPAdapter(max_retries=retry))

logger = logging.getLogger(__name__)
warn = logger.warning

# Warnings are appended to this file (opened when first needed).
log_path = Path(os.environ.get('CHEMUTILS_LOG', join('~', '.chemutils.log'))).expanduser()
shared_log_format = '%(levelname)s[%(filename)s:%(lineno)d]: %(message)s'

_file_handler = logging.FileHandler(log_path, mode='a', delay=True)
_file_handler.setFormatter(logging.Formatter(f'%(asctime)s {shared_log_format}\n'))
logger.addHandler(_file_handler)

# https://stackoverflow.com/questions/52161735/logging-warn-add-stacktrace
class WarnWithStackHandler(logging.StreamHandler):
    """Prints warnings along with the stack of the code that issued them.
    """
    def emit(self, record):
        if record.levelno == logging.WARNING:
            print()
            # (without frames from the logging module, or this one)
            stack = traceback.extract_stack()[:-1]
            while stack and Path(stack[-1].filename) == Path(logging.__file__):
                stack.pop()

            for line in traceback.format_list(stack):
                print(line, end='')

        super().emit(record)

# To stop printing warnings (they still go to log_path): logger.removeHandler(sh)
sh = WarnWithStackHandler()
sh.setFormatter(logging.Formatter(shared_log_format))
logger.addHandler(sh)

pkg_data_dir = split(__file__)[0]

# Can't use default xlrd now that it has dropped support for XLSX files, and although
# pandas >= 1.2 might default to openpyxl if it's installed, it seems that some earlier
# pandas versions (e.g. 1.1.5) may not. Note that writing still uses xlsxwriter, which
# itself should be using openpyxl.
pandas_excel_engine = 'openpyxl'

# Lookups only return magnitudes and unit strings (not pint objects), so other code
# should convert units with this registry (`chemutils.ureg`), as pint quantities from
# different registries can not be combined.
ureg = pint.UnitRegistry()
# (so pint objects unpickled elsewhere use this registry)
pint.set_application_registry(ureg)

# Persistent cache of lookup results (an SQLite database). If the CHEMUTILS_CACHE
# environment variable is set, it is used as the path instead (e.g. in tests).
default_cache_file = os.path.expanduser(join('~', '.chemutils_cache.sqlite'))
# Cache format used before the SQLite one. Migrate with scripts/migrate_cache.py
old_pickle_cache_file = os.path.expanduser(join('~', '.chemutils_cache.p'))

chemspipy_api_key = os.environ.get('CHEMSPIDER_API_KEY')

# (env var takes precedence over ~/.chemspipy_api_key)
if chemspipy_api_key is None:
    chemspipy_api_key_file = os.path.expanduser(join('~', '.chemspipy_api_key'))
    if exists(chemspipy_api_key_file):
        with open(chemspipy_api_key_file, 'r') as f:
            chemspipy_api_key = f.read().strip()

# In descending order of the certainty each ID provides as to whether we found
# the right PubChem entry.
chem_id_types = [
    'cid',
    'inchikey',
    'inchi',
    # Without stereochemistry (see compound2smiles). 'isomeric_smiles' (with
    # stereochemistry) is also supported as a `convert` `to_type`.
    'smiles',
    'cas',
    'name',
]
equivalent_col_names = {
    'odor': 'name',
    'cas_number': 'cas',
}
# No matter which `chem_id_type` these are to be converted to, they will always
# return a null value. (maybe handle some other way?)
# NOTE: keys in `properties` are currently except from this. maybe they should
# not be?
hardcoded_type2null_keys = {
    'name': {
        'spontaneous firing rate',
        'odor',
    }
}
# If you add a property here, you probably also want to add preferred units for
# it in `property2preferred_units` below.
properties  = [
    'density',
    'k_henry',
    'vapor_pressure',
    'water_solubility',
    'molar_mass',
]

concentration_dimensionalities = [
    '[length]^-3 [mass]',
    # This seems to be how pint formats molarity (M).
    '[length]^-3 [substance]',
]
preferred_concentration_unit = 'g / L'
preferred_density_unit = 'g / mL'
dimensionality2preferred_units = {
    # (densities have the same dimensionality, so pubchemprops_lookup converts those to
    # preferred_density_unit instead)
    '[length]^-3 [mass]': preferred_concentration_unit,
    '[length]^-3 [substance]': 'molar',
    '[length]^-1 [mass] [time]^-2': 'kPa',
}
# Checking that pint parses all preferred units to the expected
# dimensionality.
for dim, preferred_unit_str in dimensionality2preferred_units.items():
    assert ureg(preferred_unit_str).check(dim)

property2preferred_units = {
    'concentration': preferred_concentration_unit,
    'density': preferred_density_unit,
    # Just what NIST uses. May want to change.
    'k_henry': 'mol / kg / bar',
    'vapor_pressure': 'kPa',
    'water_solubility': 'g / L',
    'molar_mass': 'g / mol',
}


def _property_series(prop: str, value=None, sources=None, units=None, **kwargs
    ) -> pd.Series:
    """Returns property lookup output (all null by default).

    Has keys `prop`, '<prop>_sources', and '<prop>_units', and '<prop>_<key>' for each
    of `kwargs`.
    """
    return pd.Series({prop: value, f'{prop}_sources': sources, f'{prop}_units': units,
        **{f'{prop}_{k}': v for k, v in kwargs.items()}
    })


def is_concentration_unit(units):
    """
    Returns whether a `pint` `Unit`/`Quantity` has dimensionality we accept as a
    concentration.
    """
    # This is to detect the pint Unit class, as opposed to Quantity.
    if not hasattr(units, 'u'):
        # Since units don't have the same check methods as quantities.
        quantity = ureg.Quantity(1, units)
    else:
        quantity = units

    return any([quantity.check(d) for d in concentration_dimensionalities])


def in_expected_units(expected_units, test_units):
    """Returns whether `expected_units` and `test_units` are equivalent.

    `expected_units` must be one of `None` or a `str` that is either
        'concentration', 'density', or something `pint` can parse.

    `test_units` is a `pint` `Quantity` or `Unit` object.

    Equivalent in the sense that they share the same number of factors of each
    base unit type.
    """
    # Lone numbers, for instance, parse to dimensionless units.
    if expected_units is None:
        return not test_units.dimensionless

    elif expected_units == 'concentration':
        return is_concentration_unit(test_units)

    elif expected_units == 'density':
        return (test_units.dimensionality ==
            ureg(preferred_density_unit).dimensionality
        )

    else:
        return test_units.dimensionality == ureg(expected_units).dimensionality


def change_series_units(series, new_unit_str, new_units=None):
    """Converts units of Series with at least the keys ['value','units']
    """
    if new_units is None:
        # implement only if useful to refactor other stuff to use this
        # (calculate from new_unit_str)
        raise NotImplementedError

    unit_str = series['units']
    if not new_units.check(unit_str):
        raise ValueError(f'units {unit_str} inconsistent with {new_unit_str}')

    m = (series['value'] * ureg(unit_str)).to(new_units).magnitude

    new_series = series.copy()
    new_series['value'] = m
    new_series['units'] = new_unit_str
    return new_series


def change_frame_units(df, prop=None, prop_col='property', value_col='value',
    units_col='units'):
    """Returns single-property DataFrame with property in preferred units.

    Assuming that input is not null in any of the columns this fn uses, all
    unit strings are parseasble, and all property values are numeric.
    """
    if prop is None:
        assert prop_col in df.columns, 'define prop or prop_col'
        unique_props = df[prop_col].unique()
        assert len(unique_props) == 1
        prop = unique_props[0]

    assert prop in property2preferred_units, \
        'not supporting manual units being passed in for now'

    preferred_unit_str = property2preferred_units[prop]
    preferred_units = ureg(preferred_unit_str)

    series_fn = lambda s: change_series_units(s, preferred_unit_str,
        preferred_units
    )
    return df.apply(series_fn, axis=1)


hardcoded_properties_xlsx = join(pkg_data_dir, 'hardcoded_properties.xlsx')
def load_manual_properties(to_update):
    """Updates passed in `dict` with chemical properties in spreadsheet."""
    if not exists(hardcoded_properties_xlsx):
        return

    df = pd.read_excel(hardcoded_properties_xlsx, engine=pandas_excel_engine)
    df.dropna(how='all', inplace=True)

    expected_cols = [
        'property',
        'inchi',
        'value',
        'units',
        'sources'
    ]
    assert all([e in df.columns for e in expected_cols])

    # So you don't have to keep re-entering these for consecutive rows that
    # would have the same value.
    df['property'] = df['property'].ffill()
    df['inchi'] = df['inchi'].ffill()

    cols_have_null = df[expected_cols].isnull().any()
    if cols_have_null.any():
        nc_names = list(np.array(expected_cols)[cols_have_null])
        warn(f'{hardcoded_properties_xlsx} has null values in columns: {nc_names}\n\n'
            f'Rows with null values in ANY of these columns will not be used!'
        )
        df.dropna(subset=expected_cols, inplace=True)

    assert not df.duplicated(subset=['property','inchi','sources']).any()

    name_cols = [c for c in df.columns if c.startswith('name')]
    if len(name_cols) > 0:
        assert len(name_cols) == 1
        name_col = name_cols[0]
        # Since we are ffilling 'inchi', to avoid null values causing false
        # positive duplicates.
        df[name_col] = df[name_col].ffill()
        assert not df.duplicated(subset=['property', name_col, 'sources']).any()

    # NOTE: not using groupby(...).apply(change_frame_units), since pandas>=3 excludes
    # the grouping column ('property') from the frames passed to the fn.
    df = pd.concat([
        change_frame_units(gdf, prop=prop) for prop, gdf in df.groupby('property')
    ])
    assert (df.groupby('property').units.nunique() == 1).all()

    # So far, this is the only "to_type" I've used for any property lookups.
    from_type = 'inchi'
    if from_type not in to_update:
        to_update[from_type] = dict()
    to_update_f = to_update[from_type]

    # This averages across sources, if there are multiple.
    df = df.groupby(['property','inchi'], sort=False).agg({
        'value': 'mean', 'units': 'first', 'sources': '\n'.join
    })
    for row in df.itertuples():
        to_type, inchi = row.Index
        assert to_type in properties

        if to_type not in to_update_f:
            to_update_f[to_type] = dict()

        to_update_f[to_type][inchi] = _property_series(to_type, row.value, row.sources,
            row.units
        )


# Values used instead of lookup results (these take precedence over the cache, even
# with ignore_cache=True), unless `ignore_hardcoded` is True. The name -> CAS values
# differ from what looking up the name gives (as of 2026-10-23). Properties from
# hardcoded_properties.xlsx are added below.
hardcoded = {
    'name': {
        'cas': {
            'g-hexalactone': '57129-70-1',
            'g-decalactone': '1336-42-1',
            'a -pinene': '102640-64-2',
            # might not be right
            '(1S)-(+)-3-carene': '1352058-87-7',
            'a -terpineol': '10482-56-1',
            # trans? racemic? (what to do w/ cas number in that case?)
            'linalool oxide': '10448-29-0',
            'E2-hexenal': '1335-39-3',
            'E2-hexenol': '114411-82-4',
            'E3-hexenol': '544-12-7',
            'Z3-hexenol': '82658-58-0',
            'E2-hexenyl acetate': '10094-40-3'
        },
        'cid': {
            # This should be the name after normalize name
            # (it's a hardcoded correction).
            '(E,E)-2,4-hexadienal': 11829369
        }
    }
}
load_manual_properties(hardcoded)

# If True, `hardcoded` values are not used (e.g. to check whether they still differ
# from lookup results).
ignore_hardcoded = False


def is_series(obj) -> bool:
    """Returns whether `obj` is a pandas Series.
    """
    # NOTE: previously checked type(obj).__module__.startswith('pandas.'), but in
    # pandas>=3 Series.__module__ is just 'pandas', so that was False for all Series.
    return isinstance(obj, pd.Series)


def isnull(x):
    """
    Returns True if either:
    1) `pd.isnull(val)` would return True without raising an error.
    2) `val` is a Series, and all values are null.
    ...otherwise returns False.
    """
    # We specifically only want this behavior for pandas Series,
    # not DataFrame or numpy ndarray objects.
    if is_series(x):
        return x.isnull().all()
    else:
        # Forcing to bool first so ValueError is triggered here, rather than
        # actually causing an error in the calling code. DataFrames or numpy
        # arrays will raise this ValueError, with the message:
        # "The truth value of an array with more than one element is ambiguous.
        # Use a.any() or a.all()"
        is_null = pd.isnull(x)
        try:
            return bool(is_null)

        except ValueError:
            return False


def _values_equal(a, b) -> bool:
    """Returns whether two cache values are equal (handling Series and nulls)."""
    if is_series(a) or is_series(b):
        return is_series(a) and is_series(b) and a.equals(b)

    if isnull(a) or isnull(b):
        return isnull(a) and isnull(b)

    return bool(a == b)


@functools.lru_cache(maxsize=None)
def _get_hardcoded():
    """Returns `hardcoded` and `hardcoded_type2null_keys` values, with normalized keys.

    Returned dict is of the form {(from_type, to_type): {normalized_key: value}}.

    Built on first use, since the normalization fns are defined below. These values
    take precedence over (and are never written to) the persistent cache.
    """
    normed = dict()
    for ft, tts in hardcoded.items():
        for tt, f2t_dict in tts.items():
            normed_f2t = normed.setdefault((ft, tt), dict())
            for k, v in f2t_dict.items():
                nk = normalize_id(ft, k)
                if nk in normed_f2t and not _values_equal(normed_f2t[nk], v):
                    raise ValueError(f'hardcoded {ft}->{tt} values for keys that '
                        f'normalize to {nk!r} conflict'
                    )
                normed_f2t[nk] = v

    for ft, null_keys in hardcoded_type2null_keys.items():
        for nk in null_keys:
            nk = normalize_id(ft, nk)
            for tt in chem_id_types:
                normed.setdefault((ft, tt), dict()).setdefault(nk, None)

    return normed


_type_tag = '__chemutils_type__'

def _to_json_obj(value):
    if value is None or isinstance(value, (bool, str)):
        return value

    if isinstance(value, (int, np.integer)):
        return int(value)

    if isinstance(value, (float, np.floating)):
        return float(value)

    if is_series(value):
        return {_type_tag: 'Series',
            'dtype': str(value.dtype),
            'name': _to_json_obj(value.name),
            'index': [_to_json_obj(x) for x in value.index],
            'values': [_to_json_obj(x) for x in value],
        }

    if isinstance(value, Compound):
        return {_type_tag: 'Compound', 'record': value.record}

    return {_type_tag: 'pickle',
        'data': base64.b64encode(pickle.dumps(value)).decode('ascii')
    }


def _from_json_obj(obj):
    if not isinstance(obj, dict):
        return obj

    value_type = obj[_type_tag]
    if value_type == 'Series':
        return pd.Series([_from_json_obj(x) for x in obj['values']],
            index=[_from_json_obj(x) for x in obj['index']],
            dtype=obj['dtype'], name=_from_json_obj(obj['name'])
        )

    elif value_type == 'Compound':
        return Compound(obj['record'])

    elif value_type == 'pickle':
        return pickle.loads(base64.b64decode(obj['data']))

    raise ValueError(f'unrecognized cache value type: {value_type}')


def _encode_value(value) -> str:
    """Returns JSON str for storing `value` in the cache database."""
    return json.dumps(_to_json_obj(value))


def _decode_value(text: str):
    return _from_json_obj(json.loads(text))


def cache_path() -> Path:
    """Returns path to the persistent cache (CHEMUTILS_CACHE env var, if set).

    Only read when the cache is first used in a process (or after
    `_reset_cache_state()`).
    """
    return Path(os.environ.get('CHEMUTILS_CACHE', default_cache_file)).expanduser()


_cache_schema_version = 1

def _connect(path) -> sqlite3.Connection:
    """Returns connection to cache database at `path`, creating it if needed.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    # isolation_level=None for autocommit, so each write is saved (and visible to
    # other processes) immediately.
    conn = sqlite3.connect(path, timeout=30, isolation_level=None)
    # WAL + synchronous=NORMAL make each commit cheap (no fsync per write).
    conn.execute('PRAGMA journal_mode=WAL')
    conn.execute('PRAGMA synchronous=NORMAL')

    version = conn.execute('PRAGMA user_version').fetchone()[0]
    if version not in (0, _cache_schema_version):
        conn.close()
        raise RuntimeError(f'cache at {path} has unsupported schema version {version}')

    # key has no declared type, so CIDs are stored as integers and other IDs as text.
    # value is JSON (see _encode_value).
    conn.execute('CREATE TABLE IF NOT EXISTS entries ('
        'from_type TEXT NOT NULL, '
        'to_type TEXT NOT NULL, '
        'key NOT NULL, '
        'value TEXT NOT NULL, '
        'manual INTEGER NOT NULL DEFAULT 0, '
        'updated REAL NOT NULL, '
        'PRIMARY KEY (from_type, to_type, key)'
        ') WITHOUT ROWID'
    )
    # Only setting if needed, as it requires a write lock (which would make opening
    # wait for any other process writing).
    if version == 0:
        conn.execute(f'PRAGMA user_version = {_cache_schema_version}')

    return conn


_insert_entry_sql = 'INSERT OR REPLACE INTO entries VALUES (?, ?, ?, ?, ?, ?)'

def _entry_row(from_type: str, to_type: str, key, value, *, manual: bool = False,
    updated: Optional[float] = None) -> tuple:
    """Returns row of `entries` table (see `_connect`), for `_insert_entry_sql`.
    """
    if updated is None:
        updated = time.time()
    return (from_type, to_type, key, _encode_value(value), int(manual), updated)


def _iter_entries(conn: sqlite3.Connection, from_type: Optional[str] = None,
    to_type: Optional[str] = None):
    """Yields (from_type, to_type, key, value, manual) for saved entries.

    Optionally only for entries with the specified `from_type` and/or `to_type`.
    """
    conditions = {'from_type': from_type, 'to_type': to_type}
    conditions = {k: v for k, v in conditions.items() if v is not None}
    where = ''
    if len(conditions) > 0:
        where = ' WHERE ' + ' AND '.join(f'{k} = ?' for k in conditions)

    for ft, tt, key, value, manual in conn.execute('SELECT from_type, to_type, key, '
        f'value, manual FROM entries{where}', list(conditions.values())):

        yield ft, tt, key, _decode_value(value), manual


# Sentinel returned by _cache_get when there is no entry.
_missing = object()

_conn = None

# (from_type, to_type, key) -> (value [or _missing], whether value was written in this
# process)
_memo = dict()
# Set by clear_cache(). If True, only values written in this process are used.
_session_only = False

def _reset_cache_state():
    """Closes cache database and clears in-memory state (e.g. to change cache path).
    """
    global _conn, _session_only

    if _conn is not None:
        _conn.close()

    _conn = None
    _memo.clear()
    _session_only = False


def _forget_parent_conn():
    global _conn
    # Not closing connection inherited from parent process, as that can break the
    # parent's locks.
    _conn = None

# (not available on Windows, which also can't fork)
if hasattr(os, 'register_at_fork'):
    os.register_at_fork(after_in_child=_forget_parent_conn)


def _get_conn() -> sqlite3.Connection:
    global _conn
    if _conn is not None:
        return _conn

    path = cache_path()
    _conn = _connect(path)

    if (path == Path(default_cache_file) and exists(old_pickle_cache_file) and
        _conn.execute('SELECT 1 FROM entries LIMIT 1').fetchone() is None):

        warn(f'cache at {path} is empty, but old cache exists at '
            f'{old_pickle_cache_file}. migrate it with: '
            f'python {join(pkg_data_dir, "scripts", "migrate_cache.py")}'
        )

    return _conn


def _cache_get(from_type: str, to_type: str, key, *, session_only: bool = False,
    use_hardcoded: bool = True):
    """Returns cached value for `key` (normalized here), or `_missing`.

    Hardcoded values take precedence over everything else.

    Args:
        session_only: if True, only values written during this process are used
            (to implement `ignore_cache=True`).
    """
    key = normalize_id(from_type, key)

    if use_hardcoded and not ignore_hardcoded:
        hardcoded_f2t = _get_hardcoded().get((from_type, to_type))
        if hardcoded_f2t is not None and key in hardcoded_f2t:
            return hardcoded_f2t[key]

    memo_key = (from_type, to_type, key)
    session_only = session_only or _session_only
    memoized = _memo.get(memo_key)
    if memoized is not None and (memoized[1] or not session_only):
        return memoized[0]

    if session_only:
        return _missing

    row = _get_conn().execute('SELECT value FROM entries WHERE from_type = ? AND '
        'to_type = ? AND key = ?', memo_key
    ).fetchone()
    value = _missing if row is None else _decode_value(row[0])
    _memo[memo_key] = (value, False)
    return value


def _cache_set(from_type: str, to_type: str, key, value, *, manual: bool = False
    ) -> None:
    """Saves `value` for `key` (normalized here) to the cache.

    Raises ValueError (without saving) if `value` is a property with unexpected units.
    """
    key = normalize_id(from_type, key)
    _check_property_value(to_type, value)
    _get_conn().execute(_insert_entry_sql,
        _entry_row(from_type, to_type, key, value, manual=manual)
    )
    _memo[(from_type, to_type, key)] = (value, True)


def _cache_delete(from_type: str, to_type: str, keys=None) -> None:
    """Deletes `keys` (normalized here) from cache. Deletes all keys if None.
    """
    conn = _get_conn()
    if keys is None:
        conn.execute('DELETE FROM entries WHERE from_type = ? AND to_type = ?',
            (from_type, to_type)
        )
        memo_keys = [k for k in _memo if k[:2] == (from_type, to_type)]
    else:
        memo_keys = [(from_type, to_type, normalize_id(from_type, k)) for k in keys]
        conn.executemany('DELETE FROM entries WHERE from_type = ? AND to_type = ? '
            'AND key = ?', memo_keys
        )

    for memo_key in memo_keys:
        _memo.pop(memo_key, None)


class _SubCache(MutableMapping):
    """Persistent cache values for one (from_type, to_type) pair.

    Keys are normalized on access, so raw and normalized IDs refer to the same entry.
    Does not include hardcoded values.
    """
    def __init__(self, from_type: str, to_type: str):
        self.from_type = from_type
        self.to_type = to_type

    def __getitem__(self, key):
        value = _cache_get(self.from_type, self.to_type, key, use_hardcoded=False)
        if value is _missing:
            raise KeyError(key)
        return value

    def __setitem__(self, key, value):
        _cache_set(self.from_type, self.to_type, key, value)

    def __delitem__(self, key):
        self[key]
        _cache_delete(self.from_type, self.to_type, [key])

    def __iter__(self):
        rows = _get_conn().execute('SELECT key FROM entries WHERE from_type = ? AND '
            'to_type = ?', (self.from_type, self.to_type)
        ).fetchall()
        return iter([r[0] for r in rows])

    def __len__(self):
        return _get_conn().execute('SELECT COUNT(*) FROM entries WHERE from_type = ? '
            'AND to_type = ?', (self.from_type, self.to_type)
        ).fetchone()[0]

    # (overriding these MutableMapping defaults, which query each key separately)
    def clear(self):
        _cache_delete(self.from_type, self.to_type)

    def items(self):
        return [(k, v) for _, _, k, v, _ in
            _iter_entries(_get_conn(), self.from_type, self.to_type)
        ]


class _Cache:
    """View of persistent cache, as `cache[from_type][to_type][chem_id]`.
    """
    def __getitem__(self, from_type: str):
        return _FromTypeCache(from_type)

class _FromTypeCache:
    def __init__(self, from_type: str):
        self.from_type = from_type

    def __getitem__(self, to_type: str) -> _SubCache:
        return _SubCache(self.from_type, to_type)

cache = _Cache()


def clear_cache(*args):
    """
    clear_cache(): ignores all previously saved values for the rest of this process,
        without changing the saved cache (values looked up after this are still saved,
        and used).

    clear_cache(from_type, to_type): deletes all saved values for this conversion.

    clear_cache(from_type, to_type, vals): deletes saved values for input IDs `vals`.
    """
    global _session_only

    if len(args) == 0:
        _memo.clear()
        _session_only = True
        return

    if args[0] not in chem_id_types:
        raise ValueError(f'from_type must be one of {chem_id_types}')

    if len(args) == 2:
        _cache_delete(*args)

    elif len(args) == 3:
        ft, tt, vals = args
        if not pd.api.types.is_list_like(vals):
            vals = [vals]
        _cache_delete(ft, tt, vals)

    else:
        raise ValueError('can only call with 0, 2, or 3 arguments')


def _clear_null_cache_vals():
    """Deletes all null values (including all-null Series) from cache.
    """
    to_delete = dict()
    for from_type, to_type, key, value, _ in _iter_entries(_get_conn()):
        if isnull(value):
            to_delete.setdefault((from_type, to_type), []).append(key)

    for (from_type, to_type), keys in to_delete.items():
        _cache_delete(from_type, to_type, keys)


def delete_cache():
    """Deletes the cache database at `cache_path()`.
    """
    path = cache_path()
    _reset_cache_state()
    for p in (path, Path(f'{path}-wal'), Path(f'{path}-shm')):
        p.unlink(missing_ok=True)


@functools.lru_cache(maxsize=None)
def _preferred_units(prop: str):
    return ureg(property2preferred_units[prop])


def _check_property_value(prop: str, value) -> None:
    """Raises ValueError if `value` (a lookup output for `prop`) has unexpected units.
    """
    if prop not in property2preferred_units:
        return

    unit_key = f'{prop}_units'
    if not is_series(value) or prop not in value or unit_key not in value:
        if isnull(value):
            return

        raise ValueError(f'{prop} value should be a Series with {prop} and {unit_key}:'
            f'\n{value}'
        )

    if pd.isnull(value[prop]):
        return

    preferred_units = _preferred_units(prop)
    units = value[unit_key]
    if pd.isnull(units) or not preferred_units.check(units):
        raise ValueError(f'{prop} units {units} inconsistent with {preferred_units}')


def allowed_kwarg(fn: Callable, kwarg_name: str) -> bool:
    """Returns whether `kwarg_name` can be passed as a keyword argument to `fn`.
    """
    args, varargs, varkw, defaults, kwonlyargs, kwonlydefaults, annotations = \
        inspect.getfullargspec(fn)

    # https://stackoverflow.com/questions/71106783
    signature = inspect.signature(fn)
    positional_only_args = [
        x.name for x in signature.parameters.values() if x.kind == x.POSITIONAL_ONLY
    ]

    # if fn takes **kwargs, it can take any keyword argument (unless maybe also
    # specified as positional only? don't think i care about that case)
    if varkw is not None:
        return True

    # I guess something doesn't need a default to be passable as a kwarg...
    return (
        (kwarg_name in args and kwarg_name not in positional_only_args)
        or kwarg_name in kwonlyargs
    )


# Whether lookups (including nested ones, e.g. of a CID within a property lookup) are
# currently ignoring the cache, as set by calling a @cached fn with ignore_cache=True.
_ignoring_cache = contextvars.ContextVar('chemutils_ignoring_cache', default=False)

def cached(fn):
    """Caches outputs of lookup fn `<from_type>2<to_type>` in the persistent cache.

    The input is normalized (with `normalize_id(from_type, ...)`) before the cache
    lookup and before being passed to `fn`, so raw and normalized inputs share the same
    cache entries. Null inputs are passed to `fn` without any caching.

    Adds these kwargs to the wrapped fn:
        ignore_cache: if True, values saved before this process are ignored (values
            computed earlier in this process, and hardcoded values, are still used).
            Also applies to any lookups (cached fns or `convert`) `fn` makes.

        ignore_cache_null: if True, cached null values are re-computed.

        leave_nonnull_cache: if True, null values will not overwrite cached non-null
            values.
    """
    # Would need to modify if I wanted kwargs for the decorator that explicitly
    # specify the from / to types.
    # https://stackoverflow.com/questions/627501
    from_type, to_type = fn.__name__.split('2', 1)

    # So this cached decorator doesnt change function __name__
    @functools.wraps(fn)
    def cached_fn(key, ignore_cache=False, ignore_cache_null=False,
        leave_nonnull_cache=False, **kwargs):

        verbose = bool(kwargs.get('verbose', False))
        if 'verbose' in kwargs and not allowed_kwarg(fn, 'verbose'):
            kwargs.pop('verbose')

        key = normalize_id(from_type, key)
        if isnull(key):
            return fn(key, **kwargs)

        ignore_cache = ignore_cache or _ignoring_cache.get()
        value = _cache_get(from_type, to_type, key, session_only=ignore_cache)
        if value is not _missing and not (ignore_cache_null and isnull(value)):
            if verbose:
                print(f"{fn.__name__}('{key}') returning {value} from cache")
            return value

        token = _ignoring_cache.set(ignore_cache)
        try:
            value = fn(key, **kwargs)
        finally:
            _ignoring_cache.reset(token)

        if leave_nonnull_cache and isnull(value):
            cached_value = _cache_get(from_type, to_type, key)
            if cached_value is not _missing and not isnull(cached_value):
                return value

        _cache_set(from_type, to_type, key, value)
        return value

    return cached_fn


# From Wikipedia page on InChI
inchi_prefix2layer = {
    'c': 'atom connections',
    'h': 'hydrogen atoms',
    'p': 'proton sublayer',
    'q': 'charge sublayer',
    'b': 'double bonds and cumulenes',
    't': 'tetrahedral stereochemistry of atoms and allenes',
    # (whether the actual structure is the mirror image of the one /t describes)
    'm': 'inverted tetrahedral stereochemistry',
    's': 'type of stereochemistry information',
    # b,t,m,s[,h] can all also have isotopic information for
    # "isotopic stereochemistry"
    'i': 'isotopic layer',
    # 'may end with "o" sublayer; never included in standard InChI'
    # what is "o" sublayer?
    'f': 'fixed-H layer',
    # also never in standard InChI
    'r': 'reconnected metal atoms layer'
}


def basic_inchi(inchi, no_stereochem=False):
    """Returns `inchi` with only connectivity (and double bond stereochemistry) layers.

    Keeps the formula, connections (/c), hydrogens (/h), and double bond
    stereochemistry (/b) layers, so stereoisomers that only differ in tetrahedral
    stereochemistry (e.g. enantiomers) have the same output. If `no_stereochem=True`,
    the /b layer is also removed, so E/Z isomers also have the same output.
    """
    parts = inchi.split('/')
    if no_stereochem:
        keep = {'c','h'}
    else:
        keep = {'c','h','b'}
    return '/'.join(parts[:2] + [p for p in parts[2:] if p[0] in keep])


def inchi_layers(inchi):
    # First part should just be 1 or 1S (for "standard InChI")
    # Second part is the chemical formula w/o a prefix.
    prefixes = [x[0] for x in inchi.split('/')[2:]]
    layers = [inchi_prefix2layer[p] for p in prefixes]
    return layers


def inchi_layer_set(inchis: pd.Series) -> set:
    """Returns set of all layers (descriptions in `inchi_prefix2layer`) in `inchis`.
    """
    return set().union(*inchis.apply(lambda x: set(inchi_layers(x))))


def is_one2one(df: pd.DataFrame, x: str, y: str, *, dropna: bool = True) -> bool:
    """Returns True if x and y are 1:1 with each other in df.

    Handling of null still broken.
    """
    unique_combos = df[[x, y]].drop_duplicates()

    # needed this *in addition to* passing dropna to groupby and dropna
    if dropna:
        unique_combos = unique_combos.dropna()

    # https://stackoverflow.com/questions/50643386
    ny_for_each_x = unique_combos.groupby(x, dropna=dropna)[y].nunique(dropna=dropna)
    if (ny_for_each_x != 1).any():
        return False

    nx_for_each_y = unique_combos.groupby(y, dropna=dropna)[x].nunique(dropna=dropna)
    if (nx_for_each_y != 1).any():
        return False

    return True


def print_full_df(df, index=False):
    with pd.option_context('display.max_colwidth', None):
        print(df.to_string(index=index))


def inchi_diff_in_details(df, no_stereochem=False):
    """Prints InChIs with several combinations of other IDs, and InChIs that only differ
    in stereochemistry.

    Args:
        no_stereochem: passed to `basic_inchi` (if True, InChIs that only differ in
            double bond stereochemistry are also printed).
    """
    assert 'inchi' in df.columns
    other_keys = [x for x in df.columns if x in chem_id_types]
    other_keys += [x for x in df.columns if x in equivalent_col_names]
    cols_to_show = other_keys + ['inchi']
    ncdf = df.drop_duplicates(subset=cols_to_show)

    print_header = f'InChI with multiple combinations of {other_keys}:'
    printed = False
    for gn, gdf in ncdf.groupby('inchi'):
        if len(gdf) <= 1:
            continue

        if not printed:
            print(print_header)
            printed = True

        print_full_df(gdf[cols_to_show])
        print('')

    if printed:
        print('')

    inchi_counts = df.inchi.value_counts(sort=False)

    df = df.drop_duplicates(subset='inchi').copy()

    df['inchi_counts'] = df.inchi.apply(lambda i: inchi_counts.at[i])

    df['basic_inchi'] = df.inchi.apply(basic_inchi, no_stereochem=no_stereochem)

    cols_to_show += ['inchi_counts', 'basic_inchi']

    print_header = \
        'Multiple standard InChI that map to InChI w/o chirality info:'
    printed = False
    for gn, gdf in df.groupby('basic_inchi'):
        if len(gdf) <= 1:
            continue

        if not printed:
            print(print_header)
            printed = True

        print_full_df(gdf[cols_to_show])
        for c in convert(gdf.inchi, to_type='cid'):
            print(pubchem_url(c))
        print('')

    if printed:
        print('')


def count_inchi_layers(inchis):
    counts = Counter(inchis.apply(inchi_layers).agg('sum'))
    return counts


def convertable(chem_id_type: str) -> bool:
    """Returns whether the string input is recognized as a convertable type.

    (valid for the `to_type` and `from_type` arguments to the `convert` function)
    """
    if chem_id_type in chem_id_types or chem_id_type in equivalent_col_names:
        return True
    return False


def _cache_failed_lookup(from_type: str, chem_id, skip=tuple()) -> None:
    """Caches null values for `chem_id`, for each of `chem_id_types` not in `skip`.
    """
    for tt in chem_id_types:
        # To not overwrite existing (e.g. hardcoded) values for other types.
        if tt not in skip and _cache_get(from_type, tt, chem_id) is _missing:
            _cache_set(from_type, tt, chem_id, None)


def convert(chem_id, from_type=None, to_type='inchi', *, verbose: bool = False,
    allow_nan: bool = False, allow_conflicts: bool = True, ignore_cache: bool = False,
    exclude_cols=tuple(), dropna: bool = True,
    report_missing: bool = True, report_conflicts: bool = True,
    report_many2one: bool = True, try_non_normalized: bool = True,
    check_one2one: bool = False, keep_originals: bool = False, orig_prefix: str = 'orig_',
    used_id_col: Optional[str] = None):
    """Converts a chemical ID (or a Series / Index / DataFrame of them) to `to_type`.

    Input IDs are always normalized (with `normalize_<from_type>`, if defined) before
    cache lookup, and values are cached under the normalized ID, so raw and normalized
    inputs share the same cache entries.

    For a Series / Index, `from_type` defaults to its name. For a DataFrame, its index
    and columns are converted if their names are ID types. Otherwise, all of its columns
    named as ID types are converted, and a `to_type` column is added with the output of
    the highest priority (in `chem_id_types` order) column with an output.

    Args:
        dropna: whether to drop input subset where lookups fail (if setting this to
            False, probably want `allow_nan=True`)

        allow_nan: whether to err if any lookups fail (does not depend on
            `dropna=False`)

        ignore_cache: if True, values saved before this process are ignored (values
            looked up earlier in this process, and hardcoded values, are still used).
            If 'if_null', cached null values are looked up again.

        report_many2one: for Series / Index inputs, whether to print any outputs that
            several (normalized) inputs convert to.

        try_non_normalized: if lookup of the normalized ID fails, also try looking up
            the input ID as passed (result is still cached under the normalized ID).

        check_one2one: for Series / Index inputs, if True, raises ValueError if several
            (normalized) inputs convert to the same output.

        keep_originals: for DataFrame inputs, whether to add a column with the original
            values of each converted index / column (named with `orig_prefix`).

        used_id_col: for DataFrame inputs converted from several ID columns, adds a
            column of this name, with which column the output came from.
    """
    def conversion_fail_errmsg():
        return f'conversion from {from_type} to {to_type} failed for'

    def conversion_fail_err():
        msg = conversion_fail_errmsg() + f' {chem_id}'

        if not allow_nan:
            raise ValueError(msg)

        elif report_missing:
            print(msg)


    valid_ignore_cache_vals = (True, False, 'if_null')
    if ignore_cache not in valid_ignore_cache_vals:
        raise ValueError(
            f'ignore_cache must be one of: {valid_ignore_cache_vals}'
        )

    could_keep_originals = False
    if hasattr(chem_id, 'shape') and len(chem_id.shape) > 0:
        to_type = equivalent_col_names.get(to_type, to_type)

        # This covers both Index and Series (converting their values).
        if len(chem_id.shape) == 1:
            if from_type is None:
                if chem_id.name is None:
                    raise ValueError('either pass from_type or name axes')

                from_type = equivalent_col_names.get(chem_id.name, chem_id.name)

            # (allow_nan=True and report_missing=False, so failures can all be reported
            # together below)
            fn = lambda x: convert(x, from_type=from_type, to_type=to_type,
                verbose=verbose, ignore_cache=ignore_cache, allow_nan=True,
                report_missing=False, try_non_normalized=try_non_normalized
            )
            converted = chem_id.map(fn)
            converted.name = to_type

            missing = converted.isnull()
            if missing.any():
                err_str = conversion_fail_errmsg() + ':\n'
                # maybe don't dropna if not allow_nan (especially if that's
                # all that would be in the err message after the header...)
                missing = chem_id[missing].dropna().unique()
                for m in missing:
                    err_str += str(m) + '\n'

                if not allow_nan:
                    raise ValueError(err_str)
                # Second condition is to check it wasn't all null
                elif report_missing and len(missing) > 0:
                    print(err_str)
            del missing

            if check_one2one or report_many2one:
                pairs = pd.DataFrame({
                    'input': [normalize_id(from_type, x) for x in chem_id],
                    'output': list(converted),
                }).dropna().drop_duplicates()
                many2one = pairs[pairs.output.duplicated(keep=False)]
                if len(many2one) > 0:
                    err_str = f'multiple {from_type} converted to the same {to_type}:\n'
                    for output, gdf in many2one.groupby('output', sort=False):
                        err_str += f'{output}: {list(gdf.input)}\n'

                    if check_one2one:
                        raise ValueError(err_str)
                    print(err_str)

            if dropna:
                converted = converted.dropna()

            return converted

        # This covers DataFrames
        elif len(chem_id.shape) == 2:
            could_keep_originals = True
            # So that if there are already an cols that have this prefix that we
            # did not add, those aren't printed in some circumstances where we
            # only want to print the ones modified in this function.
            orig_cols_added = []

            if from_type is not None:
                raise NotImplementedError('only conversion of named indices'
                    ' or columns supported for DataFrames'
                )

            df = chem_id.copy()

            def add_original(name: str, values) -> None:
                orig_col = orig_prefix + name
                values = pd.Series(np.asarray(values), index=df.index)
                if orig_col in df.columns:
                    if not df[orig_col].astype(object).equals(values.astype(object)):
                        raise ValueError(f'{orig_col} already in columns, with different '
                            f'values from {name}'
                        )
                else:
                    df[orig_col] = values
                orig_cols_added.append(orig_col)

            # allow_conflicts and exclude_cols are only relevant for DataFrame
            # case, and therefore do not need to be passed.
            # from_type should always be inferred in DataFrame case as well.
            kwargs = {
                'to_type': to_type,
                'allow_nan': allow_nan,
                'ignore_cache': ignore_cache,
                'dropna': dropna,
                'report_missing': report_missing,
                'report_many2one': report_many2one,
                'try_non_normalized': try_non_normalized,
                'check_one2one': check_one2one,
                'verbose': verbose
            }
            converted_an_index = False
            if convertable(df.index.name):
                if keep_originals:
                    add_original(df.index.name, df.index)

                df.index = convert(df.index, **kwargs)
                converted_an_index = True

            if convertable(df.columns.name):
                if keep_originals:
                    add_original(df.columns.name, df.columns)

                df.columns = convert(df.columns, **kwargs)
                converted_an_index = True

            if not converted_an_index:
                kwargs['allow_nan'] = True
                kwargs['dropna'] = False
                kwargs['report_missing'] = False

                attempts = dict()
                for c in df.columns:
                    if c not in exclude_cols and convertable(c):
                        ft = equivalent_col_names.get(c, c)

                        if keep_originals:
                            add_original(c, df[c])

                        if verbose:
                            print(f'applying normalize_{ft} to input column')

                        df[c] = df[c].map(functools.partial(normalize_id, ft))
                        attempts[c] = convert(df[c], **kwargs)

                attempts = pd.DataFrame(attempts)

                missing = attempts.isnull().all(axis=1)
                if missing.any():
                    cols = orig_cols_added + list(attempts.columns)
                    err_str = (f'conversion from {[c for c in cols]} to {to_type} '
                        'failed for:\n'
                    )
                    missing = df[missing].drop_duplicates(subset=cols).dropna(
                        subset=cols, how='all')[cols]

                    with pd.option_context('display.max_colwidth', None):
                        err_str += missing.to_string()

                    if not allow_nan:
                        raise ValueError(err_str)
                    elif report_missing and len(missing) > 0:
                        print(err_str)
                del missing

                input_priorities = {k: p for p, k in enumerate(chem_id_types)}

                def priority(id_type):
                    return input_priorities[equivalent_col_names.get(id_type, id_type)]

                cols_in_order = sorted(attempts.columns, key=priority)
                # (unique() keeps the order values first appear in, so the first value
                # is from the highest priority column with an output)
                values = attempts[cols_in_order].apply(lambda x: x.dropna().unique(),
                    axis=1
                )
                conflicts = values.apply(lambda x: len(x) > 1)

                df[to_type] = values.apply(
                    lambda x: None if len(x) == 0 else x[0]
                )
                if used_id_col is not None:
                    df[used_id_col] = attempts[cols_in_order].apply(lambda x:
                        next((c for c in cols_in_order if pd.notnull(x[c])), None),
                        axis=1
                    )

                if dropna:
                    df.dropna(subset=[to_type], inplace=True)

                if conflicts.any():
                    if not allow_conflicts or report_conflicts:
                        print('conflicting lookup results:')
                        printed = set()
                        for i, conflict in attempts[conflicts].iterrows():
                            conflict = conflict.dropna()
                            inputs = df.loc[i, conflict.keys()]
                            key = (tuple(inputs), tuple(conflict))
                            if key in printed:
                                continue
                            printed.add(key)

                            cdf = pd.DataFrame({'input': inputs,
                                to_type: conflict
                            })
                            with pd.option_context('display.max_colwidth', None):
                                print(cdf.to_string(justify='left'))

                    if not allow_conflicts:
                        raise ValueError('conflicting lookup results')

            return df

        else:
            raise ValueError('unexpected number of dimensions')

    elif from_type is None:
        raise ValueError('specify from_type if not passing pandas object w/ '
            'named axes'
        )

    if keep_originals and not could_keep_originals:
        raise ValueError('input must be a DataFrame to keep original IDs')

    from_type = equivalent_col_names.get(from_type, from_type)

    if verbose:
        print(f'trying to convert {chem_id} from {from_type} to {to_type}')

    # To short-circuit normalization in the more-common case where null exists
    # before normalization.
    if pd.isnull(chem_id):
        return chem_id

    # Cache keys are always normalized IDs (and normalization is idempotent), so raw
    # and normalized inputs share the same cache entries.
    old_chem_id = chem_id
    chem_id = normalize_id(from_type, chem_id)
    if verbose and old_chem_id != chem_id:
        print(f'normalize_{from_type}({old_chem_id}) -> {chem_id}')

    # Since sometimes (just normalize_cas, for now) normalize fns can return
    # null.
    if pd.isnull(chem_id):
        return chem_id

    # ignore_cache=True still uses values looked up earlier in this process.
    session_only = ignore_cache == True or _ignoring_cache.get()

    val = _cache_get(from_type, to_type, chem_id, session_only=session_only)
    if val is not _missing and not (ignore_cache == 'if_null' and isnull(val)):
        if verbose:
            print(f'returning {val} from cache')
        return val

    if from_type == 'cid':
        cid = chem_id
    else:
        cid = _cache_get(from_type, 'cid', chem_id, session_only=session_only)
        if cid is not _missing:
            if verbose:
                print(f'{chem_id} of type {from_type} had CID {cid} in cache')

            if isnull(cid) and ignore_cache != 'if_null':
                if verbose:
                    print('CID in cache was None!')
                conversion_fail_err()
                return None

    if cid is _missing or isnull(cid):
        f2cid_fn_name = from_type + '2cid'

        if f2cid_fn_name not in globals():
            raise NotImplementedError(f'define function {f2cid_fn_name} to support '
                f'conversion from type {from_type}'
            )

        if verbose:
            print('calling CID lookup function', f2cid_fn_name)

        cid = globals()[f2cid_fn_name](chem_id)

        if cid is None:
            if verbose:
                print(f'looking up CID for {chem_id} of type {from_type} failed!')

            if try_non_normalized and old_chem_id != chem_id:
                if verbose:
                    print('falling back to non-normalized chem_id')
                cid = globals()[f2cid_fn_name](old_chem_id)

                if cid is None and verbose:
                    print(f'looking up CID for {old_chem_id} of type {from_type} '
                        'failed!'
                    )

        if cid is None:
            _cache_failed_lookup(from_type, chem_id)
            conversion_fail_err()
            return None

        cid = normalize_cid(cid)
        if verbose:
            print(f'CID={cid}')

        _cache_set(from_type, 'cid', chem_id, cid)

    # (the <type>2cid fns cache the Compound, if their search returned it)
    compound = _cache_get('cid', 'compound', cid)
    if compound is _missing:
        t0 = time.time()
        compound = pcp.Compound.from_cid(cid)
        if verbose:
            print(f'pcp.Compound.from_cid({cid}) took {(time.time() - t0):.3f}s')

        _cache_set('cid', 'compound', cid, compound)

    if compound is None:
        if verbose:
            print(f'creating Compound from CID {cid} failed!')

        # Skipping 'cid' just because we technically did get a non-null cid before.
        _cache_failed_lookup(from_type, chem_id, skip=('cid',))
        conversion_fail_err()
        return None

    compound2t_fn_name = 'compound2' + to_type
    if compound2t_fn_name not in globals():
        raise NotImplementedError(f'define function {compound2t_fn_name} to support '
            f'conversion to type {to_type}'
        )

    if verbose:
        print(f'calling function {compound2t_fn_name} to get {to_type} from Compound')

    compound2t_fn = globals()[compound2t_fn_name]

    if allowed_kwarg(compound2t_fn, 'verbose'):
        to_type_val = compound2t_fn(compound, verbose=verbose)
    else:
        to_type_val = compound2t_fn(compound)

    to_type_val = _normalize_value(to_type, to_type_val)
    _cache_set(from_type, to_type, chem_id, to_type_val)

    if verbose and to_type_val is None:
        print(f'conversion of {chem_id} from Compound to {to_type} failed!')
        conversion_fail_err()

    return to_type_val


def normalize_id(chem_id_type: str, chem_id):
    """Returns `chem_id` normalized with `normalize_<chem_id_type>`, if defined.

    All normalization fns are idempotent and return null for null input.
    """
    norm_fn = _normalize_fn(chem_id_type)
    if norm_fn is None:
        return chem_id

    return norm_fn(chem_id)


@functools.lru_cache(maxsize=None)
def _normalize_fn(chem_id_type: str) -> Optional[Callable]:
    chem_id_type = equivalent_col_names.get(chem_id_type, chem_id_type)
    return globals().get(f'normalize_{chem_id_type}')


# Types whose values (not just keys) are normalized before caching.
_value_normalized_types = ('cid', 'inchi', 'inchikey', 'cas')

def _normalize_value(chem_id_type: str, value):
    """Returns lookup output `value`, normalized if its type is in
    `_value_normalized_types`.
    """
    if chem_id_type not in _value_normalized_types:
        return value

    return normalize_id(chem_id_type, value)


# E/Z stereo descriptors (e.g. '(e)-', '(e,z)-'), which normalize_name keeps uppercase,
# despite lowercasing everything else.
_ez_descriptor_re = re.compile(r'\(([ez](?:,[ez])*)\)-')

def _null_preserving(fn):
    """Decorates a fn of one argument to return null inputs as-is.
    """
    @functools.wraps(fn)
    def wrapped(x):
        if pd.isnull(x):
            return x
        return fn(x)

    return wrapped


@_null_preserving
def normalize_name(name):
    """Returns lowercase name, with some corrections to how names are written.

    For example, trailing parentheticals are removed, words are joined across hyphens,
    and prefixes like 'a-' are expanded to 'alpha-'. Returns None for empty names.
    """
    return _normalize_name(name)


@functools.lru_cache(maxsize=2**16)
def _normalize_name(name: str) -> Optional[str]:
    name = name.replace('(CAS)', '').strip().lower()

    # Removing trailing parentheticals of at least 4 characters (e.g. ' (branched)',
    # ' (99%)'). All of them (not just the last), so this is idempotent. Not if that
    # would leave nothing (e.g. for a name that is all one parenthetical).
    stripped = re.sub(r'(\s*\([^\(\)]{4,}\))+$', '', name)
    if len(stripped.strip()) > 0:
        name = stripped

    parts = name.split()
    if len(parts) == 0:
        return None

    normed_name = parts[0]

    for a, b in zip(parts, parts[1:]):
        # Only joining across a hyphen (e.g. 'a -pinene' -> 'a-pinene'), and never
        # after a comma. Joining in other cases changed which compound some names refer
        # to (e.g. 'ethyl 3-hydroxyhexanoate' -> 'ethyl3-hydroxyhexanoate', which
        # PubChem resolves to a different CID).
        if a[-1] == ',' or not (a[-1] == '-' or b[0] == '-'):
            b = ' ' + b
        normed_name += b

    corrections = {
        'linalool oxide': 'trans-linalool oxide',
        '4-ethyl guaiacol': '4-ethylguaiacol',
        # (copied from normalize_name in natural_odors)
        'dihydrazide ethanediimidic acid':
            'ethanediimidic acid, dihydrazide',
        # Since search of 'E3-hexenol' fails and '(E)-3-hexenol' returns
        # compound for the (Z)-3 version.
        'e3-hexenol': '(E)-hex-3-en-1-ol',
        # (this one doesn't actually yield the stereochem i want, which is
        # CID 11829369, not 8901 or 637564, so I'm going to hardcode the CID
        # for now anyway)
        '2,4-hexadienal, (e,e)-': '(E,E)-2,4-hexadienal'
    }
    if normed_name in corrections:
        return corrections[normed_name]

    prefix_corrections = {
        'a-': 'alpha-',
        'b-': 'beta-',
        'g-': 'gamma-',
        'e2-': '(E)-2-',
        'e3-': '(E)-3-',
        'z2-': '(Z)-2-',
        'z3-': '(Z)-3-',
    }
    def correct_prefix(n):
        for p, c in prefix_corrections.items():
            if n.startswith(p):
                n = c + n[len(p):]
        return n

    normed_name = correct_prefix(normed_name)

    def ester_acid(n):
        if ' ester ' in n and ' acid' in n:
            return True
        return False

    def rev_ester_acid(n):
        parts = n.split()
        return ' '.join([parts[2], parts[3] + ',', parts[0], parts[1]])

    indicator2correction = (
        (ester_acid, rev_ester_acid),
    )
    for i, c in indicator2correction:
        if i(normed_name):
            print(f'indicator {i.__name__} was True for {normed_name}')
            normed_name = c(normed_name)
            print(f'applying correction {c.__name__} to get {normed_name}.')
            # Since correction can move a part needing prefix correction to the start.
            normed_name = correct_prefix(normed_name)

    # So '(E)-'/'(Z)-' (from corrections above, or in input) are not lowercased.
    normed_name = _ez_descriptor_re.sub(lambda m: m.group(0).upper(), normed_name)

    return normed_name


@_null_preserving
def normalize_cas(cas):
    """Returns CAS number without quotes or whitespace (as some sources have).

    Returns None for '0-00-0', which some sources use for a missing CAS number.
    """
    normed_cas = ''.join(cas.replace('"','').split())

    if normed_cas == '0-00-0':
        return None

    return normed_cas


def _normalize_str(chem_id, prefix: str = ''):
    """Strips whitespace and (any repeats of) case-insensitive `prefix`.

    Returns null input as-is, and None for empty strings.
    """
    if pd.isnull(chem_id):
        return chem_id

    chem_id = chem_id.strip()
    while prefix and chem_id.lower().startswith(prefix.lower()):
        chem_id = chem_id[len(prefix):].strip()

    return chem_id if len(chem_id) > 0 else None


def normalize_inchi(inchi):
    """Returns InChI without whitespace or 'InChI=' prefix (as used in cache).

    >>> normalize_inchi('InChI=1S/CH4/h1H4')
    '1S/CH4/h1H4'
    """
    return _normalize_str(inchi, 'InChI=')


def _inchi_with_prefix(inchi: str) -> str:
    """Returns InChI with exactly one 'InChI=' prefix (as websites expect).
    """
    return 'InChI=' + normalize_inchi(inchi)


def normalize_inchikey(inchikey):
    """Returns uppercase InChIKey without whitespace or 'InChIKey=' prefix.
    """
    inchikey = _normalize_str(inchikey, 'InChIKey=')
    return inchikey.upper() if isinstance(inchikey, str) else inchikey


def normalize_cid(cid) -> Optional[int]:
    """Returns PubChem compound ID as an int (e.g. from float in Series w/ NaN).

    >>> normalize_cid(6654.0)
    6654
    """
    if pd.isnull(cid):
        return cid

    if isinstance(cid, str):
        cid = float(cid.strip())

    int_cid = int(cid)
    if int_cid != cid:
        raise ValueError(f'CID must be an integer. got: {cid!r}')

    return int_cid


def normalize_smiles(smiles):
    """Returns SMILES without surrounding whitespace."""
    return _normalize_str(smiles)


pubchem_root_url = 'https://pubchem.ncbi.nlm.nih.gov'
@_null_preserving
def pubchem_url(cid: int) -> str:
    # .0f to work w/ floats (e.g. from a Series that is a mix of ints and NaN)
    return f'{pubchem_root_url}/compound/{cid:.0f}'


def pubchem_search_url(text: str) -> str:
    return f'{pubchem_root_url}/#query={urllib.parse.quote(text)}'


def inchi2pubchem_url(inchi, **kwargs):
    """Returns PubChem URL for `inchi`, or None if it can not be found.

    `kwargs` are passed to `convert`.
    """
    if pd.isnull(inchi):
        return inchi

    kwargs.setdefault('allow_nan', True)
    cid = convert(inchi, from_type='inchi', to_type='cid', **kwargs)
    return pubchem_url(cid)


def _search_cid(identifier, namespace: str, verbose: bool = False) -> Optional[int]:
    """Returns CID of PubChem's first (best) match for `identifier`, or None if none.

    Also caches the Compound for that CID (so `convert` does not need to request it).
    Identifiers PubChem can not parse (e.g. malformed InChIs) are taken to have no match.
    """
    try:
        results = pcp.get_compounds(identifier, namespace)
    except pcp.BadRequestError as err:
        if verbose:
            print(f'PubChem could not parse {namespace}={identifier!r}: {err}')
        return None
    except urllib.error.URLError as err:
        warn(f'{err}\nReturning None')
        return None

    # (unknown InChIs give a result without a CID)
    results = [r for r in results if r.cid is not None]
    if len(results) == 0:
        return None

    if len(results) > 1 and verbose:
        print(f'got {len(results)} results from {namespace}={identifier!r} (using '
            'first):'
        )
        for r in results:
            print(f'{pubchem_url(r.cid)} {r.inchi}')
        print()

    compound = results[0]
    cid = normalize_cid(compound.cid)
    if _cache_get('cid', 'compound', cid) is _missing:
        _cache_set('cid', 'compound', cid, compound)

    return cid


def name2cid(name, verbose=False):
    return _search_cid(name, 'name', verbose=verbose)


def cas2cid(cas, verbose=False):
    # (PubChem has CAS numbers as synonyms, so they are searched as names)
    return _search_cid(cas, 'name', verbose=verbose)


def inchikey2cid(inchikey, verbose=False):
    return _search_cid(inchikey, 'inchikey', verbose=verbose)


def inchi2cid(inchi, verbose=False):
    return _search_cid(_inchi_with_prefix(inchi), 'inchi', verbose=verbose)


def compound2cid(compound: Compound) -> int:
    """Returns PubChem compound ID.

    PubChem has separate CIDs for each stereoisomer (and another for the compound with
    stereochemistry unspecified), so the CID does reflect stereochemistry.
    """
    return compound.cid


def compound2name(compound: Compound) -> str:
    """Returns PubChem's preferred IUPAC name.

    Includes stereochemistry, when defined (e.g. '(E)-hex-2-enal' for CID 5281168). See
    `compound2title` for the (usually common) name PubChem shows at the top of its page.
    """
    return compound.iupac_name


def compound2title(compound: Compound) -> Optional[str]:
    """Returns name at the top of PubChem's page for compound (e.g. 'Linalool').

    Usually a common name, unlike `compound2name`'s IUPAC name (e.g.
    '3,7-dimethylocta-1,6-dien-3-ol'). Available as `convert(..., to_type='title')`.
    Requires another PubChem request.
    """
    properties = pcp.get_properties(['Title'], compound.cid, 'cid')
    if len(properties) == 0:
        return None
    return properties[0].get('Title')


def compound2inchi(compound: Compound) -> str:
    """Returns standard InChI, without the 'InChI=' prefix.

    Includes stereochemistry layers (/b, /t, /m, /s), when defined (e.g.
    '1S/C6H10O/c1-2-3-4-5-6-7/h4-6H,2-3H2,1H3/b5-4+' for CID 5281168). See
    `basic_inchi(..., no_stereochem=True)` to remove them.
    """
    inchi = compound.inchi
    assert inchi.startswith('InChI=')
    inchi = inchi[6:]
    assert len(inchi) > 0
    return inchi


def compound2inchikey(compound: Compound) -> str:
    """Returns standard InChIKey.

    Includes stereochemistry, when defined, as it is a hash of the full standard InChI
    (the second block is 'UHFFFAOYSA' when there are no stereochemistry or other
    non-connectivity layers).
    """
    return compound.inchikey


def compound2cas(compound: Compound, *, verbose: bool = False) -> Optional[str]:
    """Returns first CAS number in PubChem synonyms (None if there are none).

    CAS numbers are assigned per substance, so may or may not be specific to one
    stereoisomer (some PubChem compounds also list CAS numbers for other forms).
    """
    cas_numbers = [
        x for x in compound.synonyms if re.match(r'(\d{2,7}-\d\d-\d)', x)
    ]
    cid = compound.cid
    if len(cas_numbers) == 0:
        if verbose:
            print(f'no CAS numbers found in PubChem for {cid=}')
            print(pubchem_url(cid))

        cas = None
    else:
        if verbose and len(cas_numbers) > 1:
            print(f'multiple PubChem CAS for {cid=}: {cas_numbers}')
            print('returning first of these CAS')

        # had previously sorted and returned first, but at least for my alpha-pinene
        # (cid=6654) test, the first CAS in synonyms (80-56-8) is what we want, not the
        # first in sorted order (2437-95-8).
        cas = cas_numbers[0]

    return cas


def _compound_smiles(compound: Compound, names) -> Optional[str]:
    """Returns value of first SMILES property of `compound` with name in `names`.

    Reading record directly (rather than via pubchempy attributes), since PubChem renamed
    its SMILES properties in 2025, and records saved before only have the old names.
    """
    smiles_props = {p['urn'].get('name'): p['value']['sval']
        for p in compound.record.get('props', []) if p['urn'].get('label') == 'SMILES'
    }
    for name in names:
        if name in smiles_props:
            return smiles_props[name]

    return None


def compound2smiles(compound: Compound) -> Optional[str]:
    """Returns PubChem SMILES *without* stereochemistry (or isotope) information.

    Only describes connectivity, so all stereoisomers of a compound have the same value
    (e.g. 'CCCC=CC=O' for CID 5281168, (E)-2-hexenal). See `compound2isomeric_smiles`
    for SMILES with stereochemistry.

    PubChem calls this the "Connectivity" SMILES ("Canonical" before 2025).
    """
    return _compound_smiles(compound, ('Connectivity', 'Canonical'))


def compound2isomeric_smiles(compound: Compound) -> Optional[str]:
    """Returns PubChem SMILES *with* stereochemistry (and isotope) information.

    E.g. 'CCC/C=C/C=O' for CID 5281168, (E)-2-hexenal. Same as `compound2smiles` output
    for compounds without defined stereochemistry. Available as
    `convert(..., to_type='isomeric_smiles')`.

    PubChem calls this the "Absolute" SMILES ("Isomeric" before 2025).
    """
    return _compound_smiles(compound, ('Absolute', 'Isomeric'))


def fmt_id_type(chem_id_type):
    """Returns how to write ID type (e.g. 'InChI' for 'inchi') in plot labels, etc.
    """
    if chem_id_type in equivalent_col_names:
        chem_id_type = equivalent_col_names[chem_id_type]

    if chem_id_type not in chem_id_types:
        raise KeyError(f'unrecognized chem_id type. options are: {chem_id_types}')

    fmt_dict = {
        'cas': 'CAS',
        'inchi': 'InChI',
        'inchikey': 'InChIKey',
        'smiles': 'SMILES',
        'cid': 'PubChem compound ID'
    }
    if chem_id_type in fmt_dict:
        return fmt_dict[chem_id_type]
    else:
        return chem_id_type


def name2compound(name, verbose=False) -> Optional[Compound]:
    """Returns PubChem Compound for best match for `name` (None if none).
    """
    cid = name2cid(name, verbose=verbose)
    if cid is None:
        return None

    compound = _cache_get('cid', 'compound', cid)
    if compound is _missing:
        compound = pcp.Compound.from_cid(cid)
    return compound


def normalize_chem_ids(df, chem_id, keep_originals=False, **kwargs):
    """Normalize on chem_id then back convert to names.
    """
    df = convert(df, to_type=chem_id, keep_originals=keep_originals, **kwargs)
    df['name'] = convert(df[chem_id], to_type='name', **kwargs)
    return df


@_null_preserving
def chemspider_url(csid):
    assert np.isclose(csid, int(csid))

    # ':.0f' just so if input comes from a a Series where NaN was added (thus converting
    # any prior int types to float), the URLs still work.
    return f'http://www.chemspider.com/Chemical-Structure.{csid:.0f}.html'


# NOTE: chemspider's policy is that scraping is not allowed, so may want to take further
# steps to rate limit / appear as a normal browser (mechanize w/ similar setup to
# VCF scraper?).
def chemspider_experimental_density(csid, verbose=False):
    url = chemspider_url(csid)
    if verbose:
        print(url)
    response = session.get(url)
    # crude rate limit
    time.sleep(0.5)

    html = response.text
    soup = BeautifulSoup(html, features='html.parser')
    tag = soup.find(text='Experimental Density:',
        attrs={'class': 'user_data_property_name'}
    )
    if tag is None:
        return None

    table = tag.parent.parent.next_sibling.next_sibling
    density = None
    for row in table.find_all('td'):
        parts = row.text.split()
        try:
            curr_density = float(parts[0].strip())
        except ValueError:
            if verbose:
                print(url)
                print('Failed to parse density from row text:', row.text)
            #raise
            return None

        units = ureg(parts[1])
        assert units.check(preferred_density_unit)

        curr_density = (curr_density * units).to(preferred_density_unit
            ).magnitude

        # To filter out what seem to be plainly wrong densities for e.g.
        # http://www.chemspider.com/Chemical-Structure.4444608.html
        # (20 g/mL)
        # Lowered threshold b/c 6.4 g/mL from Biosynth at this source:
        # http://www.chemspider.com/Chemical-Structure.66391.html
        if curr_density < 2 and curr_density > 0.1:
            density = curr_density

    # So we don't have a dict w/ density null and other data not-null.
    if density is None:
        return None

    return pd.Series({
        'density': density,
        'density_sources': url,
        'density_units': preferred_density_unit,
    })


nist_url_prefix = 'http://webbook.nist.gov'
def inchi2nist_webbook_url(inchi):
    query = urllib.parse.urlencode({'InChI': normalize_inchi(inchi), 'Units': 'SI'})
    return f'{nist_url_prefix}/cgi/cbook.cgi?{query}'


def inchi2nist_henrys_law_url(inchi):
    return inchi2nist_webbook_url(inchi) + '&Mask=10'


def inchi2nist_vapor_pressure_url(inchi):
    # NOTE: often there are tables here that aren't the antoinne equation parameters we
    # want (and are there simpler records of vapor pressure, that aren't intended to
    # work across temperatures, that exist in the NIST webbook? maybe in a place i'm not
    # looking?)
    return inchi2nist_webbook_url(inchi) + '&Mask=4'


@cached
def cid2molar_mass(cid):
    if isnull(cid):
        return _property_series('molar_mass')

    compound = _cache_get('cid', 'compound', cid)
    if compound is _missing:
        compound = pcp.Compound.from_cid(cid)

    return _property_series('molar_mass', float(compound.molecular_weight),
        pubchem_url(cid), 'g/mol'
    )


@cached
def inchi2molar_mass(inchi):
    cid = convert(inchi, from_type='inchi', to_type='cid', allow_nan=True)
    return cid2molar_mass(cid)


def try_parse_units(string):
    """Returns either pint units or None. Raises ValueError when pint does.

    pint may also raise other errors I do not specifically anticipate in here.
    """
    curr_units = None
    try:
        curr_units = ureg.parse_units(string)

    except pint.errors.UndefinedUnitError:
        pass

    # This will happen if there is an unmatched lparen.
    except TokenError:
        pass

    # This happens if string == ':'
    except AttributeError:
        pass

    # '/ext/' yields this error, w/ text 'missing unary operator "/"'
    except pint.errors.DefinitionSyntaxError:
        pass

    # 'experimentally-derived', yields this w/ text:
    # "unsupported operand type(s) for -: 'ParserHelper' and 'ParserHelper'"
    except TypeError:
        pass

    # '()' yields this
    except AssertionError:
        pass

    return curr_units


def degree_str2unit(temp_unit_str):
    if temp_unit_str == 'C':
        temp_unit = ureg.degC

    elif temp_unit_str == 'F':
        temp_unit = ureg.degF

    elif temp_unit_str == 'K':
        temp_unit = ureg.degK

    # Although since the regex specifies there must be one of the CFK
    # characters, this should not be reached unless the regex changes.
    else:
        raise ValueError('unexpected temperature unit str')

    return temp_unit


def _match_groups_range_mean(groups):
    """Returns mean of range when given match groups.

    Assumes first two groups have str min and max of range.
    """
    range_min = float(groups[0])
    range_max = float(groups[1])
    assert range_min < range_max, 'range not in expected order'
    quantity = range_min + (range_max - range_min) / 2
    return quantity


float_re = r'(\d+(?:,\d{3})*(?:\.\d*)?)'
float_range_re = float_re + r'\s?(?:-|to)\s?' + float_re
# (see `_stated_temps_c` for other formats temperatures may be in)
_at_re = r'(?:@|at)\s*'
degree_re = r'\s*°\s*([CFK])'
temp_re = _at_re + float_re + degree_re
temp_range_re = _at_re + float_range_re + degree_re

# Densities in the 'Density' section of PubChem pages that are not of the compound as a
# liquid (or solid), e.g. 'Relative density of the vapour/air-mixture at 20 °C (air = 1):
# 1', 'Density (gas): 1.00045 g/l', 'Bulk density = 6.9 wt/gal', '1.05 at 59 °F for
# 10.17% weight/weight solution', 'critical density: 0.2734', and conversion factors
# ('1 MG/L= 278 PPM').
_not_liquid_density_re = re.compile(
    r'\b(?:air|vapou?r|gas|gases|critical|bulk|solutions?|conversion|ppm)\b',
    re.IGNORECASE
)
# What may come before a relative density (its value then the first thing after it):
# 'Density:', 'Relative density (water = 1):', 'Sp gr:', and a temperature, e.g. '(57
# °F): 0.81'.
_density_label_re = re.compile(
    r'^\s*(?:(?:relative\s+)?density|specific\s+gravity|sp\.?\s*gr\.?)?\s*'
    r'(?:\([^)]*\))?\s*:?\s*',
    re.IGNORECASE
)
# The temperature of the water a relative density is to, e.g. '/4 °C' in '0.876 at
# 15 °C/4 °C'.
_relative_to_temp_re = re.compile(r'/\s*-?\d+(?:\.\d+)?\s*°\s*[CF]')
_cu_cm_re = re.compile(r'[/.]\s*cu\.?\s*cm\b', re.IGNORECASE)
# Solvents other than water, which solubilities in are not used as water solubilities.
_other_solvent_re = re.compile(
    r'\b(?:alcohols?|ethanol|methanol|ethers?|acetone|benzene|chloroform|oils?|'
    r'glycerol|glycerin|glycols?|hexane|toluene|DMSO|(?:organic\s+)?solvents)\b',
    re.IGNORECASE
)
# An amount in a volume of water, e.g. '2.5 g in 100 mL water', '1 ml in 40 ml water',
# '1 gram dissolves in 1000 mL water, 80 mL boiling water, 0.8 mL alcohol'.
_amount_in_volume_re = re.compile(float_re + r'\s*(mg|g|grams?|kg|ml|l)\b[^;]*?\bin\s+'
    r'(?:about\s+)?' + float_re + r'\s*(ml|l)\b(?:\s+of)?\s+water\b', re.IGNORECASE
)
# e.g. 'Soluble in about 720 parts water' (1 g in 720 mL water, as pharmacopeias use)
_parts_water_re = re.compile(r'\bin\s+(?:about\s+)?' + float_re +
    r'\s+parts\s+(?:of\s+)?water\b', re.IGNORECASE
)
# e.g. 'In water, 6,700 ppm at 25 °C'
_ppm_in_water_re = re.compile(r'\bin\s+water\b[^;]*?' + float_re + r'\s*ppm\b',
    re.IGNORECASE
)
# e.g. '0.43% (by wt) in water'
_weight_percent_in_water_re = re.compile(float_re +
    r'\s*%\s*\(?\s*(?:by\s+)?(?:wt|weight)\.?\s*\)?[^;]*?\bin\s+water\b', re.IGNORECASE
)

def parse_temperature_c(string):
    """Returns None or scalar temparture in celsius and part of string not used.
    """
    match = re.search(temp_re, string)
    if match is None:
        match = re.search(temp_range_re, string)
        if match is None:
            return None, string
        else:
            groups = match.groups()
            temp_unit = degree_str2unit(groups[-1])
            rmean = _match_groups_range_mean(groups)
            temperature = ureg.Quantity(rmean, temp_unit).to(ureg.degC)

            string = re.sub(temp_range_re, '', string)
    else:
        groups = match.groups()
        temp_unit = degree_str2unit(groups[-1])
        temperature = ureg.Quantity(float(groups[0]), temp_unit).to(ureg.degC)

        # Remove matched part of the string before search for units and
        # quantities below (with the goal of making the search below easier).
        string = re.sub(temp_re, '', string)

    return temperature.magnitude, string


def _parse_float(string: str) -> Optional[float]:
    """Returns number in `string` (e.g. '1,000', '1.5X10+4', '1.5X10-4'), or None.
    """
    try:
        return float(string.replace(',',''))
    except ValueError:
        # Assuming there will be exactly one digit in the exponent.
        match = re.fullmatch(float_re + r'X10((?:\+|-)\d)', string)
        if match is None:
            return None
        groups = match.groups()
        return float(groups[0].replace(',','') + 'e' + groups[-1])


def parse_pubchemprops_string(string, expected_units=None, target_temp_c=None,
    url=None):
    """
    Returns a dict with quantities, units (uses `pint`), conditions they were
    measured at (e.g. temperature), and sources, from a string returned via
    `pubchemprops`.

    Any values not parseable will be have `None` values in the returned `dict`.

    If the `quantity` or `units` would be `None`, all values in returned `dict`
    are set to `None`.

    If there is data parseable at multiple temperatures, the data from the
    temperature closest to `target_temp_c` will be used.

    Assumes any values that would otherwise return exactly 0.0 return `None`,
    to not cause calculations to fail for lack of the real small value.

    Raises ValueError if there are several numbers (other than a range or temperature)
    or units that could be the value.
    """
    orig_str = string

    if target_temp_c is None:
        target_temp_c = 25.0

    keys_to_parse = ['quantity', 'units', 'temperature_c', 'source',
        'from_range'
    ]
    null_dict = {k: None for k in keys_to_parse}

    # In my experience so far, only things that would either not be parseable
    # or that have data at two temperatures (and should be parseable when
    # split) have a semicolon (from solubility test data).
    if ';' in string:
        smallest_diff = None
        data = None
        # Parts with a value but no stated temperature, e.g. the first of
        # '1.174 (USCG, 1999) - Denser than water; will sink'
        undated = []
        error = None
        for part in string.split(';'):
            try:
                curr_data = parse_pubchemprops_string(part,
                    expected_units=expected_units, target_temp_c=target_temp_c,
                    url=url
                )
            # (other parts may still have a value, e.g. the density in
            # 'DENSITY OF "SATURATED" AIR: ...; 1 MG/L= 408 PPM & 1 PPM= 2.45 MG/CU M')
            except ValueError as err:
                if error is None:
                    error = err
                continue

            if curr_data['quantity'] is None:
                continue

            curr_temp = curr_data['temperature_c']
            if curr_temp is not None:
                curr_diff = abs(curr_temp - target_temp_c)
                if smallest_diff is None or curr_diff < smallest_diff:
                    data = curr_data
                    smallest_diff = curr_diff
            else:
                undated.append(curr_data)

        # Only if it is the one value (several without temperatures could be of
        # different things).
        if data is None and len(undated) == 1:
            data = undated[0]

        # This coverts strings that have no parseable portions, but still
        # have a semicolon, such as:
        # "Miscible in ethanol, ethyl ether; slightly soluble in carbon
        #  tetrachloride"
        if data is None:
            if error is not None:
                raise error
            return null_dict

        return data

    parsed_dict = {k: None for k in keys_to_parse}

    relative_density = False
    if expected_units == 'density':
        if _not_liquid_density_re.search(string):
            return null_dict

        # The temperature of the water a relative density is to, e.g. '/4 °C' in
        # '0.876 at 15 °C/4 °C' (as at 4 C water is 1.000 g/mL, and 0.997 g/mL at 25 C,
        # these are taken as g/mL).
        string = _relative_to_temp_re.sub('', string)

        # A value with no units, as PubChem's are when relative to water (specific
        # gravity), e.g. 'Relative density (water = 1): 0.87', '0.868-0.878', or
        # '0.876 at 68 °F (USCG, 1999) - Less dense than water'.
        relative_density = bool(re.match(r'\d', _density_label_re.sub('', string, count=1)))

        # e.g. '1.02 g/cm3' (pint parses 'cm³', but not 'cm3')
        string = re.sub(r'\b(c?m)3\b', r'\1**3', string)

    elif expected_units == 'concentration':
        # Solubilities in other solvents, e.g. 'Solubility in alcohol: 1 mL in 1 mL 95%
        # ethanol'
        if _other_solvent_re.search(string) and not re.search(r'\bwater\b', string,
            re.IGNORECASE):

            return null_dict

    temp_c, string = parse_temperature_c(string)
    if temp_c is not None:
        parsed_dict['temperature_c'] = temp_c

    # Assuming there are always four digits for the year.
    source_re = r'\((([^(]+),\s(\d{4}))\)'
    match = re.search(source_re, string)
    if match is not None:
        parsed_dict['source'] = match.groups()[0]

    # So that we don't have to worry about parsing the years as our quantity.
    string = re.sub(source_re, '', string)

    if expected_units == 'concentration':
        # e.g. '2.5 g in 100 mL water' (as g/mL), or '1 ml in 40 ml water' (which needs
        # the density, so is not used)
        match = _amount_in_volume_re.search(string)
        if match is not None:
            amount, amount_unit, volume, volume_unit = match.groups()
            # (as pint takes e.g. 'ML' to be megaliters)
            amount_unit, volume_unit = amount_unit.lower(), volume_unit.lower()
            if not ureg.Quantity(1, amount_unit).check('[mass]'):
                return null_dict

            parsed_dict['quantity'] = (
                float(amount.replace(',', '')) / float(volume.replace(',', ''))
            )
            parsed_dict['units'] = ureg.parse_units(f'{amount_unit} / {volume_unit}')
            parsed_dict['from_range'] = False
            return parsed_dict

        # ppm (by mass) as mg/L, and weight % as 10 g/L per %, as dilute solutions are
        # about as dense as water
        for pattern, scale, units in (
            (_parts_water_re, None, 'g / mL'),
            (_ppm_in_water_re, 1, 'mg / L'),
            (_weight_percent_in_water_re, 10, 'g / L'),
            ):
            match = pattern.search(string)
            if match is None:
                continue

            value = _parse_float(match.group(1))
            if not value:
                return null_dict

            parsed_dict['quantity'] = 1 / value if scale is None else scale * value
            parsed_dict['units'] = ureg.parse_units(units)
            parsed_dict['from_range'] = False
            return parsed_dict

    quantity = None

    parsed_dict['from_range'] = False

    # First part means match either whitespace or the start of the string.
    range_re = r'(?:^|\s)' + float_range_re
    match = re.search(range_re, string)
    if match is not None:
        parsed_dict['from_range'] = True

        groups = match.groups()
        quantity = _match_groups_range_mean(groups)

        # Assuming we won't need to check for inequalities in this range case.
        parts_before_quantity = []

        # So the components of the range are not also matched in loop over
        # parts.
        string = re.sub(range_re, '', string)

    # This is to fix a common sorce of vapor pressure measurements
    # (https://haz-map.com), which seems to have all of its data in PubChem in this
    # format (e.g. '1.5 [mmHg]'.
    string = string.replace('[mmHg]', 'mmHg').replace('[mm Hg]', 'mmHg')

    # e.g. '1.0446 g/cu cm', and '1.0328 g.cu cm' (a typo on one PubChem page)
    string = _cu_cm_re.sub('/mL', string)

    # Assuming we can remove all parenthetical stuff that remains
    # (source was already parsed and removed above) without issue. Also removes empty
    # parentheses left after removing temperature, e.g. from '11 mmHg (@ 77 °F)'.
    string = re.sub(r'\(.*\)', '', string)

    parts = string.strip().split()

    remaining_parts = []
    for p in parts:
        curr_quantity = _parse_float(p)
        if curr_quantity is None:
            remaining_parts.append(p)
            continue

        if quantity is None:
            parts_before_quantity = list(remaining_parts)
            quantity = curr_quantity
        else:
            err_msg = (f'several numbers that could be the value in "{orig_str}" '
                f'(after removing temperature, source, and parentheticals: "{string}")'
            )
            if url is not None:
                err_msg = f'{url}\n{err_msg}'

            raise ValueError(err_msg)

    # (strings without a number are common, e.g. 'Miscible with ethanol')
    if quantity is None:
        return null_dict

    # If it seems the quantity is only defined as part of an inequality,
    # don't return it (because we have no way to indicate that).
    for p in parts_before_quantity:
        p = p.lower()
        if '>' in p or '<' in p or 'less' in p or 'greater' in p:
            return null_dict

    # This will make this unit of pressure something pint parses correctly.
    remaining_str = ' '.join(remaining_parts).replace('mm Hg', 'mmHg')
    remaining_parts = remaining_str.split()

    units = None
    last_part = None
    for p in remaining_parts:
        p = p.strip()

        # This is a hack to handle cases like this:
        # 'In water, 0.6 wt% (6000 mg/L) at 20 °C'
        # Assuming no need to check for inequality qualifiers in this case.
        last_quantity = None
        if (p[-1] == ')' and '(' not in p and
            last_part[0] == '(' and ')' not in last_part):

            last_quantity = _parse_float(last_part[1:])
            # If previous part had a quantiy, then use second part of
            # this parenthetical expression as a unit.
            # Need to remove the end paren or it won't parse.
            if last_quantity is not None:
                p = p[:-1]
        last_part = p

        # (e.g. 'water,' in 'In water, 8.0X10+4 mg/L at 25 °C' is not a unit)
        if 'water' in p:
            continue

        try:
            curr_units = try_parse_units(p)
            if curr_units is None:
                continue

        # This is one case that prompts this: 'g/100ml'
        # try_parse_units lets this error propagate up from ureg.parse_units
        except ValueError as e:
            # Only case encountered so far.
            assert str(e) == 'Unit expression cannot have a scaling factor.'
            ps = p.split('/')
            if len(ps) != 2:
                continue

            match = re.match(float_re, ps[1])
            if not match:
                continue

            scale_from_denom = float(match.groups()[0])
            ps[1] = re.sub(float_re, '', ps[1])
            p = '/'.join(ps)
            curr_units = try_parse_units(p)
            if curr_units is None:
                continue
            quantity = quantity / scale_from_denom

        if not in_expected_units(expected_units, curr_units):
            continue

        if last_quantity is not None:
            quantity = last_quantity

        if units is None:
            units = curr_units
        else:
            err_msg = (f'several units that could be those of the value in '
                f'"{orig_str}" (parts not parsed as numbers: "{remaining_str}")'
            )
            if url is not None:
                err_msg = f'{url}\n{err_msg}'

            raise ValueError(err_msg)

    if units is None and relative_density:
        units = ureg.parse_units('g / mL')

    # One case that triggers this: "2 % (NIOSH, 2016)"
    if units is None:
        return null_dict

    # Densities of a gas, e.g. '1.639 g/L' (hydrogen chloride's), as no liquid or solid
    # is this light.
    if (expected_units == 'density' and
        ureg.Quantity(quantity, units).to('g / mL').magnitude < 0.1):

        return null_dict

    # Since the measurement / reporting probably just doesn't have enough
    # sigfigs to report the real value, and 0.0 will cause any calculations
    # with multiplications to probably behave not as desired.
    if quantity == 0.0:
        return null_dict

    parsed_dict['units'] = units
    parsed_dict['quantity'] = quantity
    return parsed_dict


_num_re = r'(-?\d+(?:\.\d+)?)'
# Temperatures (or ranges of them) in the formats PubChem vapor pressure strings use,
# e.g. 'at 25 °C', '(77°F)', '68° F', '-1.5 °C', '294 deg K', '300 K', '20-25 °C'
_stated_temp_re = re.compile(_num_re + rf'(?:\s*(?:-|to)\s*{_num_re})?'
    r'\s*(?:(?:°|deg)\s*([CFK])|(K))\b'
)
# Number after 'at' / '@' (once stated temperatures are removed), which may be a
# temperature without units (e.g. '0.0135 mm Hg at 25 /extrapolated/').
_unitless_temp_re = re.compile(r'(?:\bat|@)\s*-?\d', re.IGNORECASE)

def _stated_temps_c(string: str) -> Optional[list]:
    """Returns temperatures [C] stated in `string`, including both ends of ranges.

    Returns None if there may also be a temperature without units.
    """
    temps_c = []
    for match in _stated_temp_re.finditer(string):
        low, high, unit_str, kelvin = match.groups()
        unit = ureg.degK if kelvin else degree_str2unit(unit_str)
        for temp in (low, high):
            if temp is not None:
                temps_c.append(ureg.Quantity(float(temp), unit).to(ureg.degC).magnitude)

    if _unitless_temp_re.search(_stated_temp_re.sub('', string)):
        return None

    return temps_c


def _pubchemprops_temp_rank(part: str, string: str, target_temp_c: float,
    max_temp_diff_c: float) -> Optional[tuple]:
    """Returns rank of value in `part` of PubChem `string`, by temperature and source.

    Lower ranks are preferred. Measured values are preferred to estimates (see
    `_is_estimate`). Then, values stated at temperatures within 1 C of `target_temp_c`
    are preferred to those within `max_temp_diff_c`, which are preferred to those with
    no stated temperature. NTP values are ranked below all other measured values (for 3
    of 28 compounds with NTP and other vapor pressures at 25 C, NTP's were >10x off).

    Returns None if the value should not be used: any temperature stated is more than
    `max_temp_diff_c` from `target_temp_c` (or may be stated without units), or no
    temperature is stated and the value is from USCG (most of which were well above
    other values at 25 C).
    """
    temps_c = _stated_temps_c(part)
    if temps_c is None:
        return None

    if len(temps_c) == 0:
        if 'USCG' in string:
            return None
        temp_rank = 2
    else:
        max_diff = max(abs(t - target_temp_c) for t in temps_c)
        if max_diff > max_temp_diff_c:
            return None
        temp_rank = 0 if max_diff <= 1 else 1

    # (NTP strings can have several values, with source only after the last)
    return (_is_estimate(string), 'NTP' in string, temp_rank)


# PubChem values that are estimates rather than measurements, e.g. '0.39 mm Hg at 25 °C
# /est/'. Values calculated from fits to measurements (e.g. '0.453 mm Hg @ 25 °C,
# calculated from experimentally derived coefficients') count as measurements.
_estimate_re = re.compile(r'\best\b|estimat|calculat|predict', re.IGNORECASE)
_from_measurements_re = re.compile(r'\bexper', re.IGNORECASE)

def _is_estimate(string: str) -> bool:
    return bool(_estimate_re.search(string)) and not _from_measurements_re.search(string)


# e.g. '5.28X10+9 mm Hg at 25 °C /extrapolated/', which is the only vapor pressure data
# on the page for Lysine (https://pubchem.ncbi.nlm.nih.gov/compound/5962), but is absurd.
# (and 'exptrapolated', a typo on some pages)
_extrapolated_re = re.compile(r'exp?trapolat', re.IGNORECASE)

_pug_view_min_interval_s = 0.2
_last_pug_view_request = None

def pubchem_heading_information(cid: int, heading: str) -> list:
    """Returns 'Information' entries in `heading` section of PubChem page for `cid`.

    `heading` is a section (e.g. 'Vapor Pressure') under 'Computed Descriptors', 'Other
    Identifiers', or 'Experimental Properties' (under 'Names and Identifiers' or
    'Chemical and Physical Properties'). Each entry is a dict with a 'Value' key (with
    either a 'StringWithMarkup' or 'Number' and 'Unit'). Returns an empty list if the
    page has no such section.

    Makes at most 5 requests per second, as PubChem asks.
    """
    global _last_pug_view_request

    if _last_pug_view_request is not None:
        wait_s = _pug_view_min_interval_s - (time.monotonic() - _last_pug_view_request)
        if wait_s > 0:
            time.sleep(wait_s)
    _last_pug_view_request = time.monotonic()

    response = session.get(f'{pubchem_root_url}/rest/pug_view/data/compound/{cid}/JSON',
        params={'heading': heading}, timeout=60
    )
    # (if the page has no section with this heading)
    if response.status_code == 404:
        return []
    response.raise_for_status()

    information = []
    for section in response.json()['Record'].get('Section', []):
        if section['TOCHeading'] not in ('Names and Identifiers',
            'Chemical and Physical Properties'):
            continue

        for subsection in section.get('Section', []):
            if subsection['TOCHeading'] not in ('Computed Descriptors',
                'Other Identifiers', 'Experimental Properties'):
                continue

            for prop_section in subsection.get('Section', []):
                if prop_section['TOCHeading'] == heading:
                    information = prop_section['Information']

    return information


def pubchemprops_lookup(inchi, prop_name, expected_units=None,
    target_temp_c=None, max_temp_diff_c=None, skip_unit_parse_err=False,
    verbose=False) -> Optional[pd.Series]:
    """Returns mean of values in `prop_name` section of PubChem page for `inchi`.

    Values (parsed with `parse_pubchemprops_string`) are converted to the same units
    before averaging. If any values are not from ranges (e.g. '0.868-0.878'), those from
    ranges are not used. Extrapolated values are never used.

    Args:
        target_temp_c: for strings with values at several temperatures, the value at the
            temperature closest to this is used (default 25 C).

        max_temp_diff_c: if not None, only values stated to be within this many C of
            `target_temp_c` are used, and values with no stated temperature are only
            used if there are none of those. Estimates are only used if there are no
            measured values (see `_pubchemprops_temp_rank` for details). All values in
            strings with several are considered. Values for solutions are also not
            used.

        skip_unit_parse_err: if True, values that can not be parsed (e.g. with several
            numbers or units) are skipped with a warning, rather than raising
            ValueError.

    Returns:
        None if there are no values to use. Otherwise a Series with the mean
        ('quantity'), its 'units', 'sources' (the PubChem URL), and 'measured' (whether
        the values used were measured, rather than estimates, or None if
        `max_temp_diff_c` is None, as estimates are not checked for then).
    """
    cid = convert(inchi, from_type='inchi', to_type='cid')
    if cid is None:
        return None

    url = pubchem_url(cid)

    if verbose:
        print(f'looking up {prop_name} of {inchi}')
        print(url)

    try:
        results = pubchem_heading_information(cid, prop_name)
    # (warning so which lookups failed is recorded in log_path)
    except (requests.RequestException, ValueError) as err:
        warn(f'looking up {prop_name} for {url} failed: {err}')
        raise

    if len(results) == 0:
        return None

    parsed_dicts = []
    n_parse_errors = 0
    for result in results:
        value_dict = result['Value']

        if 'Number' in value_dict:
            # I've only seen it have 'Unit' and 'Number' keys.
            assert len(value_dict) <= 2

            keys = ['quantity', 'units', 'temperature_c', 'source', 'from_range']
            parsed_dict = {k: None for k in keys}

            unit_str = value_dict['Unit']

            if max_temp_diff_c is not None:
                parsed_dict['temp_rank'] = _pubchemprops_temp_rank(unit_str, unit_str,
                    target_temp_c, max_temp_diff_c
                )
                if parsed_dict['temp_rank'] is None:
                    continue

            temp_c, unit_str = parse_temperature_c(unit_str)
            # `temp_c` can be `None`, but that will just overwrite w/ existing value.
            parsed_dict['temperature_c'] = temp_c

            # (as for strings above, e.g. '(68 °F)' left of 'g/l at 20 °C (68 °F)')
            unit_str = re.sub(r'\(.*\)', '', unit_str).strip()

            # This is to special case one instance where the unit_str was 'mg/mL, clear,
            # colorless'. Probably no cases where I want everything surrounding a comma,
            # though this may not always be right...
            if ',' in unit_str:
                unit_str = unit_str.split(',')[0].strip()
                assert len(unit_str) > 0, ('hack to fix unit str with extra'
                    ' info following comma failed'
                )

            unit_str = _cu_cm_re.sub('/mL', unit_str)
            # e.g. 'g/100 mL' as g/(100 mL) (pint would take it as g/100 * mL)
            unit_str = re.sub(r'/\s*(\d+(?:\.\d+)?)\s*([a-zA-Z]+)', r'/(\1 \2)', unit_str)
            try:
                # (a Quantity, so that any scaling factor is kept, e.g. 'g/100 mL')
                units = ureg.Quantity(unit_str)
            except Exception:
                units = None

            if units is None or not in_expected_units(expected_units, units):
                err_msg = f'{url}\nunexpected units: {value_dict}'
                if skip_unit_parse_err:
                    warn(err_msg)
                    n_parse_errors += 1
                    continue
                raise ValueError(err_msg)

            parsed_dict['units'] = units.u

            nums = value_dict['Number']
            assert len(nums) == 1
            parsed_dict['quantity'] = nums[0] * units.magnitude

            parsed_dict['from_range'] = False

            parsed_dicts.append(parsed_dict)

        elif 'StringWithMarkup' in value_dict:
            assert len(value_dict) == 1
            pstr = value_dict['StringWithMarkup'][0]['String']

            if _extrapolated_re.search(pstr):
                warn(f'excluding StringWithMarkup="{pstr}" from consideration because '
                    'it is extrapolated (as some absurd values were)'
                )
                continue

            if max_temp_diff_c is None:
                parts = [pstr]
            else:
                if 'solution' in pstr:
                    continue
                # (rather than letting parse_pubchemprops_string pick the value
                # closest to target_temp_c, so all within range can be used)
                parts = pstr.split(';')

            for part in parts:
                temp_rank = None
                if max_temp_diff_c is not None:
                    temp_rank = _pubchemprops_temp_rank(part, pstr, target_temp_c,
                        max_temp_diff_c
                    )
                    if temp_rank is None:
                        if verbose:
                            print(f'excluding "{part}" (temperature may not be within '
                                f'{max_temp_diff_c} C of {target_temp_c} C)'
                            )
                        continue

                if verbose:
                    print(f'trying to parse "{part}"')

                try:
                    parsed_dict = parse_pubchemprops_string(part,
                        expected_units=expected_units, target_temp_c=target_temp_c,
                        url=url
                    )

                # e.g. from strings with several numbers or units
                except ValueError as err:
                    if skip_unit_parse_err:
                        warn(str(err))
                        n_parse_errors += 1
                        continue
                    else:
                        raise

                if verbose:
                    print(f'parsed: {pformat(parsed_dict)}')

                # (strings without a value in the expected units, e.g. 'Miscible with
                # ethanol')
                if parsed_dict['quantity'] is None:
                    continue

                parsed_dict['temp_rank'] = temp_rank
                parsed_dicts.append(parsed_dict)

            assert len(value_dict['StringWithMarkup']) == 1
            # If there is an extra Key under the first element, it seems to have always
            # just been a reference to water...
        # e.g. {'Binary': [...], 'MimeType': 'text/csv'}, a table (as for hydrogen
        # chloride's density)
        else:
            warn(f'{url}\nskipping value without a number or string: {sorted(value_dict)}')
            continue

    if len(parsed_dicts) == 0:
        if n_parse_errors > 0:
            warn(f'{url}\nall {n_parse_errors} {prop_name} values that could have been '
                'used failed to parse'
            )
        return None

    measured = None
    if max_temp_diff_c is not None:
        best_rank = min(d['temp_rank'] for d in parsed_dicts)
        parsed_dicts = [d for d in parsed_dicts if d['temp_rank'] == best_rank]
        measured = not best_rank[0]

    if any([d['from_range'] == False for d in parsed_dicts]):
        parsed_dicts = [d for d in parsed_dicts if d['from_range'] == False]

    for parsed_dict in parsed_dicts:
        m = parsed_dict['quantity']
        u = parsed_dict['units']
        uq = ureg.Quantity(m, u)

        found_dims = False
        if expected_units == 'density':
            preferred_unit = preferred_density_unit
        else:
            preferred_unit = None
            for dims, curr_preferred_unit in dimensionality2preferred_units.items():

                if uq.check(dims):
                    found_dims = True
                    preferred_unit = curr_preferred_unit
                    # Assuming check will only be True once (as it should).
                    break

            assert found_dims, 'no matching dimensionality with preferred units'

        parsed_dict['quantity'] = uq.to(preferred_unit).magnitude
        # Changing type to a string here.
        # Consider keeping as a pint obj (though that may be tricky,
        # which is why I opted to do it this way).
        parsed_dict['units'] = preferred_unit

        # Molar concentrations (assuming they are either those or mass / volume), to
        # mass / volume.
        if expected_units == 'concentration' and (
            '[substance]' in uq.dimensionality.keys()):

            molar_mass = ureg.Quantity(cid2molar_mass(cid)['molar_mass'], 'g/mol')
            parsed_dict['quantity'] = (uq * molar_mass
                ).to(preferred_concentration_unit).magnitude
            parsed_dict['units'] = preferred_concentration_unit

    quantities = [d['quantity'] for d in parsed_dicts]

    units_set = {d['units'] for d in parsed_dicts if pd.notnull(d['units'])}
    if len(units_set) == 0:
        units = None
    else:
        # May need to handle cases beyond this?
        assert len(units_set) == 1
        units = units_set.pop()

    return pd.Series({
        'quantity': np.mean(quantities),
        'sources': url,
        'units': units,
        'measured': measured,
    })


@cached
def inchi2water_solubility(inchi, verbose=False):
    """Returns mean PubChem solubility in water [g/L].

    Solubilities in other solvents, and strings that cannot be parsed, are not used
    (with a warning for the latter). Not checked to be at any particular temperature.
    See `inchi2est_k_henry` to estimate Henry's law constants from this and the vapor
    pressure.
    """
    ret = pubchemprops_lookup(inchi, 'Solubility',
        expected_units='concentration', skip_unit_parse_err=True, verbose=verbose
    )
    if ret is None:
        return _property_series('water_solubility')

    return _property_series('water_solubility', ret['quantity'], ret['sources'],
        ret['units']
    )


vapor_pressure_temp_c = 25.0

@cached
def inchi2vapor_pressure(inchi, skip_unit_parse_err=True, verbose=False):
    """Returns mean PubChem vapor pressure [kPa] at (or near) 25 C.

    Only values stated to be within 5 C of 25 C are used, preferring those within 1 C.
    If there are none, values with no stated temperature are used (except USCG
    values). NTP values are only used if there are no others, and estimates only if
    there are no measured values. See `_pubchemprops_temp_rank` for details.

    Returns:
        Series with 'vapor_pressure', 'vapor_pressure_sources', 'vapor_pressure_units',
        and 'vapor_pressure_measured' (False if only estimates were used, None if there
        is no value). Values cached before 'vapor_pressure_measured' was added may lack
        it.
    """
    ret = pubchemprops_lookup(inchi, 'Vapor Pressure', expected_units='kPa',
        target_temp_c=vapor_pressure_temp_c, max_temp_diff_c=5.0,
        skip_unit_parse_err=skip_unit_parse_err, verbose=verbose
    )
    if ret is None:
        return _property_series('vapor_pressure', measured=None)

    return _property_series('vapor_pressure', ret['quantity'], ret['sources'],
        ret['units'], measured=ret['measured']
    )


# SIMPOL.1 coefficients (B_k,1, ..., B_k,4), from Table 5 of Pankow & Asher (2008)
simpol_doi_url = 'https://doi.org/10.5194/acp-8-2773-2008'
_simpol_coefs = {
    'zeroeth': (-4.26938e+02, 2.89223e-01, 4.42057e-03, 2.92846e-01),
    'carbon': (-4.11248e+02, 8.96919e-01, -2.48607e-03, 1.40312e-01),
    'carbon_acid_side_amide': (-1.46442e+02, 1.54528e+00, 1.71021e-03, -2.78291e-01),
    'aromatic_ring': (3.50262e+01, -9.20839e-01, 2.24399e-03, -9.36300e-02),
    'nonaromatic_ring': (-8.72770e+01, 1.78059e+00, -3.07187e-03, -1.04341e-01),
    'c_c': (5.73335e+00, 1.69764e-02, -6.28957e-04, 7.55434e-03),
    'c_c_c_o_in_ring': (-2.61268e+02, -7.63282e-01, -1.68213e-03, 2.89038e-01),
    'hydroxyl': (-7.25373e+02, 8.26326e-01, 2.50957e-03, -2.32304e-01),
    'aldehyde': (-7.29501e+02, 9.86017e-01, -2.92664e-03, 1.78077e-01),
    'ketone': (-1.37456e+01, 5.23486e-01, 5.50298e-04, -2.76950e-01),
    'carboxylic_acid': (-7.98796e+02, -1.09436e+00, 5.24132e-03, -2.28040e-01),
    'ester': (-3.93345e+02, -9.51778e-01, -2.19071e-03, 3.05843e-01),
    'ether': (-1.44334e+02, -1.85617e+00, -2.37491e-05, 2.88290e-01),
    'ether_alicyclic': (4.05265e+01, -2.43780e+00, 3.60133e-03, 9.86422e-02),
    'ether_aromatic': (-7.07406e+01, -1.06674e+00, 3.73104e-03, -1.44003e-01),
    'nitrate': (-7.83648e+02, -1.03439e+00, -1.07148e-03, 3.15535e-01),
    'nitro': (-5.63872e+02, -7.18416e-01, 2.63016e-03, -4.99470e-02),
    'aromatic_hydroxyl': (-4.53961e+02, -3.26105e-01, -1.39780e-04, -3.93916e-02),
    'amine_primary': (3.71375e+01, -2.66753e+00, 1.01483e-03, 2.14233e-01),
    'amine_secondary': (-5.03710e+02, 1.04092e+00, -4.12746e-03, 1.82790e-01),
    'amine_tertiary': (-3.59763e+01, -4.08458e-01, 1.67264e-03, -9.98919e-02),
    'amine_aromatic': (-6.09432e+02, 1.50436e+00, -9.09024e-04, -1.35495e-01),
    'amide_primary': (-1.02367e+02, -7.16253e-01, -2.90670e-04, -5.88556e-01),
    'amide_secondary': (-1.93802e+03, 6.48262e-01, 1.73245e-03, 3.47940e-02),
    'amide_tertiary': (-5.26919e+00, 3.06435e-01, 3.25397e-03, -6.81506e-01),
    'carbonylperoxynitrate': (-2.84042e+02, -6.25424e-01, -8.22474e-04, -8.80549e-02),
    'peroxide': (1.50093e+02, 2.39875e-02, -3.37969e-03, 1.52789e-02),
    'hydroperoxide': (-2.03387e+01, -5.48718e+00, 8.39075e-03, 1.07884e-01),
    'carbonylperoxyacid': (-8.38064e+02, -1.09600e+00, -4.24385e-04, 2.81812e-01),
    'nitrophenol': (-5.27934e+01, -4.63689e-01, -5.11647e-03, 3.84965e-01),
    'nitroester': (-1.61520e+03, 9.01669e-01, 1.44536e-03, 2.66889e-01),
}
# Temperature range of the data SIMPOL.1 was fit to.
_simpol_temp_range_c = (0.0, 120.0)

# Carbons not double bonded to O/N/S (i.e. not acyl or imine carbons).
_simpol_c = '[#6;!$([#6]=[O,N,S])]'
# N single bonded only to carbons, none of which are acyl/imine carbons.
_simpol_amine_n = '[NX3;+0;!$(N~[!#6]);!$(N-[#6]=[O,N,S])'
# (group, SMARTS, index of match atom that identifies each group). Matches with no
# index are unique atom sets.
_simpol_smarts = [
    ('hydroxyl', '[OX2H1]-[C;!$(C=[O,N,S])]', 0),
    ('aromatic_hydroxyl', '[OX2H1]-c', 0),
    ('aldehyde', '[CX3;H1,H2;!$(C-[!#6])]=[OX1]', 0),
    ('ketone', '[#6]-[CX3](=[OX1])-[#6]', 1),
    ('carboxylic_acid', '[CX3](=[OX1])-[OX2H1]', 0),
    ('ester', '[CX3](=[OX1])-[OX2]-[#6]', 0),
    # (split into the 3 SIMPOL.1 ether groups in _simpol_groups)
    ('ether', f'{_simpol_c}-,:[#8X2]-,:{_simpol_c}', 1),
    ('nitrate', '[#6]-[OX2]-[NX3](~[OX1])~[OX1]', 2),
    ('nitro', '[#6]-[NX3](~[OX1])~[OX1]', 1),
    ('amine_primary', f'{_simpol_amine_n};H2;!$(N-a)]', 0),
    ('amine_secondary', f'{_simpol_amine_n};H1;!$(N-a)]', 0),
    ('amine_tertiary', f'{_simpol_amine_n};H0;!$(N-a)]', 0),
    ('amine_aromatic', f'{_simpol_amine_n};$(N-a)]', 0),
    ('amide_primary', '[NX3;+0;H2]-[CX3]=[OX1]', 0),
    ('amide_secondary', '[NX3;+0;H1](-[#6])-[CX3]=[OX1]', 0),
    ('amide_tertiary', '[NX3;+0;H0](-[#6])(-[#6])-[CX3]=[OX1]', 0),
    ('carbonylperoxynitrate', '[CX3](=[OX1])-[OX2]-[OX2]-[NX3](~[OX1])~[OX1]', 0),
    ('peroxide', f'{_simpol_c}-[OX2]-[OX2]-{_simpol_c}', 1),
    ('hydroperoxide', f'{_simpol_c}-[OX2]-[OX2H1]', 1),
    ('carbonylperoxyacid', '[CX3](=[OX1])-[OX2]-[OX2H1]', 0),
    ('c_c', 'C=C', None),
    ('c_c_c_o_in_ring', 'C=;@C-;@[C;R]=[OX1]', None),
]
_simpol_patterns = None

def _simpol_groups(mol: Chem.Mol) -> Optional[dict]:
    """Returns counts of each SIMPOL.1 group in `mol`, or None if it is not supported.

    Only single (neutral, non-radical) molecules made of C, H, O and N are supported,
    and every O and N must be in one of the SIMPOL.1 functional groups (so e.g.
    nitriles, imines, and aromatic nitrogens are not supported).
    """
    global _simpol_patterns
    if _simpol_patterns is None:
        _simpol_patterns = [(name, Chem.MolFromSmarts(smarts), i)
            for name, smarts, i in _simpol_smarts
        ]

    if len(Chem.GetMolFrags(mol)) > 1:
        return None

    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ('C', 'O', 'N') or atom.GetNumRadicalElectrons() > 0:
            return None

    counts = {k: 0 for k in _simpol_coefs}
    counts['zeroeth'] = 1
    counts['carbon'] = sum(a.GetSymbol() == 'C' for a in mol.GetAtoms())

    aromatic_rings = []
    for ring in Chem.GetSSSR(mol):
        if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings.append(set(ring))
        else:
            counts['nonaromatic_ring'] += 1
    counts['aromatic_ring'] = len(aromatic_rings)

    # group -> list of matches (one per group)
    matches = dict()
    covered = set()
    for name, pattern, key_index in _simpol_patterns:
        group_matches = mol.GetSubstructMatches(pattern)
        if key_index is not None:
            by_key = dict()
            for match in group_matches:
                by_key.setdefault(match[key_index], match)
            group_matches = list(by_key.values())

        matches[name] = group_matches
        for match in group_matches:
            covered.update(match)

    # Checking all N/O (and charged atoms, which are only allowed in nitro / nitrate
    # groups) are in some group.
    for atom in mol.GetAtoms():
        if (atom.GetSymbol() != 'C' or atom.GetFormalCharge() != 0) and (
            atom.GetIdx() not in covered):
            return None

    for name, group_matches in matches.items():
        if name != 'ether':
            counts[name] = len(group_matches)

    for _, o_idx, _ in matches['ether']:
        oxygen = mol.GetAtomWithIdx(o_idx)
        if any(n.GetIsAromatic() for n in oxygen.GetNeighbors()):
            counts['ether_aromatic'] += 1
        elif oxygen.IsInRing():
            counts['ether_alicyclic'] += 1
        else:
            counts['ether'] += 1

    nitro_n = {match[1] for match in matches['nitro']}

    def has_nitro(atom_idx):
        return any(n.GetIdx() in nitro_n
            for n in mol.GetAtomWithIdx(atom_idx).GetNeighbors()
        )

    def connected_atoms(start_idx, exclude_idx):
        """Returns indices of atoms connected to `start_idx`, not via `exclude_idx`.
        """
        seen = {start_idx, exclude_idx}
        to_visit = [start_idx]
        while to_visit:
            for n in mol.GetAtomWithIdx(to_visit.pop()).GetNeighbors():
                if n.GetIdx() not in seen:
                    seen.add(n.GetIdx())
                    to_visit.append(n.GetIdx())

        seen.remove(exclude_idx)
        return seen

    # Nitrophenol group replaces aromatic hydroxyl if a nitro is on the same ring.
    for _, c_idx in matches['aromatic_hydroxyl']:
        if any(c_idx in r and any(has_nitro(i) for i in r) for r in aromatic_rings):
            counts['aromatic_hydroxyl'] -= 1
            counts['nitrophenol'] += 1

    # Nitroester group replaces ester if a nitro is on the acid side of the ester.
    for c_idx, _, o_idx, _ in matches['ester']:
        if any(has_nitro(i) for i in connected_atoms(c_idx, o_idx)):
            counts['ester'] -= 1
            counts['nitroester'] += 1

    # Carbons on the acid side of each amide, including the carbonyl carbon (e.g. 4 for
    # N-propyl-butyramide).
    for name in ('amide_primary', 'amide_secondary', 'amide_tertiary'):
        for match in matches[name]:
            n_idx, c_idx = match[0], match[-2]
            counts['carbon_acid_side_amide'] += sum(
                mol.GetAtomWithIdx(i).GetSymbol() == 'C'
                for i in connected_atoms(c_idx, n_idx)
            )

    return counts


# Compounds SIMPOL.1 groups can be assigned to, but for which estimates were far from
# PubChem values (see scripts/validate_est_vapor_pressure.py).
_simpol_excluded_smarts = {
    # (~30-300x too high)
    'lactone / cyclic anhydride': '[#6;R](=O)@[#8;R]',
}
_simpol_excluded_patterns = None

def _simpol_exclusion_reason(mol: Chem.Mol, counts: dict) -> Optional[str]:
    """Returns which excluded class `mol` is in, or None if not excluded.
    """
    global _simpol_excluded_patterns
    if _simpol_excluded_patterns is None:
        _simpol_excluded_patterns = {name: Chem.MolFromSmarts(smarts)
            for name, smarts in _simpol_excluded_smarts.items()
        }

    for name, pattern in _simpol_excluded_patterns.items():
        if mol.HasSubstructMatch(pattern):
            return name

    # (would be zwitterions, which are far less volatile than estimated)
    if counts['carboxylic_acid'] > 0 and any(counts[f'amine_{x}'] > 0
        for x in ('primary', 'secondary', 'tertiary')):
        return 'amino acid'

    return None


def _simpol_log10_p_atm(counts: dict, temp_k: float) -> float:
    """Returns SIMPOL.1 log10 of (subcooled) liquid vapor pressure [atm].
    """
    log10_p = 0.0
    for name, n in counts.items():
        b1, b2, b3, b4 = _simpol_coefs[name]
        log10_p += n * (b1 / temp_k + b2 + b3 * temp_k + b4 * np.log(temp_k))

    return log10_p


def inchi2est_vapor_pressure(inchi, temp_c: float = 25.0, verbose: bool = False
    ) -> pd.Series:
    """Returns vapor pressure estimated from structure, for when no measurements found.

    Intended for compounds where `inchi2vapor_pressure` returns null values. Uses the
    SIMPOL.1 group contribution method (Pankow & Asher, 2008), which only depends on the
    structure, so this requires no network access (and is not cached).

    Estimates are of the liquid vapor pressure, which for compounds that are solid at
    `temp_c` is that of the subcooled liquid (higher than that of the solid).

    Compared to PubChem values stated at 25 C (scripts/validate_est_vapor_pressure.py),
    the median error was ~2x (65% of 384 compounds within 3x, and 88% within 10x).
    Errors were lower for simple aliphatic esters, aldehydes, ketones, alcohols, and
    hydrocarbons (often within 1.5x), but these were typically underestimated:
    - furans: ~10x
    - tertiary alcohols (e.g. linalool): ~6x
    - phenols with an O-containing group ortho to the OH (e.g. guaiacol, eugenol,
      methyl salicylate), where intramolecular hydrogen bonds increase volatility:
      ~10-50x
    - ethers: ~3x

    Estimates are less accurate away from ~25 C (e.g. n-hexane at its boiling point is
    overestimated ~1.7x).

    Args:
        temp_c: temperature [C] to estimate vapor pressure at. SIMPOL.1 was fit to
            data from 0-120 C.

        verbose: if True, prints SIMPOL.1 group counts.

    Returns:
        Series with same keys as `inchi2vapor_pressure` output ('vapor_pressure_measured'
        is False). Values are null (and 'vapor_pressure_measured' None) if the
        compound is not supported, which includes compounds with atoms other than C, H,
        O and N, compounds with N/O not in one of the SIMPOL.1 groups (e.g. nitriles,
        imines, aromatic nitrogens), lactones (estimates ~30-300x too high), and amino
        acids.
    """
    null_ret = _property_series('vapor_pressure', measured=None)

    if isnull(inchi):
        return null_ret

    lo, hi = _simpol_temp_range_c
    if not lo <= temp_c <= hi:
        warn(f'{temp_c=} outside temperature range SIMPOL.1 was fit to ({lo}-{hi} C)')

    mol = Chem.MolFromInchi(_inchi_with_prefix(inchi))
    if mol is None:
        warn(f'could not parse InChI {inchi}')
        return null_ret

    counts = _simpol_groups(mol)
    if counts is None:
        warn(f'can not estimate vapor pressure of {inchi} ({Chem.MolToSmiles(mol)}): '
            'has atoms / functional groups not supported by SIMPOL.1'
        )
        return null_ret

    exclusion_reason = _simpol_exclusion_reason(mol, counts)
    if exclusion_reason is not None:
        warn(f'not estimating vapor pressure of {inchi} ({Chem.MolToSmiles(mol)}): '
            f'SIMPOL.1 estimates for class ({exclusion_reason}) are far from measured '
            'values'
        )
        return null_ret

    if verbose:
        print(f'SIMPOL.1 groups for {inchi}:')
        print(pformat({k: v for k, v in counts.items() if v != 0}))

    log10_p_atm = _simpol_log10_p_atm(counts, temp_c + 273.15)
    units = property2preferred_units['vapor_pressure']
    return _property_series('vapor_pressure',
        ureg.Quantity(10**log10_p_atm, 'atm').to(units).magnitude,
        f'SIMPOL.1 estimate ({simpol_doi_url})', units, measured=False
    )


_chemspider_session = None
@cached
def inchi2chemspider_id(inchi):
    """Returns lowest ChemSpider ID of compounds with `inchi` (None if none).

    Requires `chemspipy_api_key` to be set. Raises ChemSpiPyRateError if the API key's
    limit is exceeded (each search seems to use several API calls).
    """
    global _chemspider_session

    # If chemspipy_api_key is None, trying to construct ChemSpider would give us:
    # chemspipy.errors.ChemSpiPyAuthError: Unauthorised; Invalid API Key;
    if chemspipy_api_key is None:
        raise RuntimeError('chemspipy_api_key must be set to use inchi2chemspider_id')

    if _chemspider_session is None:
        _chemspider_session = ChemSpider(chemspipy_api_key)

    cs = _chemspider_session

    # Docs seem to claim inchi could be passed just as if it were name,
    # but my testing did not seem to bear that out.
    inchi_with_prefix = _inchi_with_prefix(inchi)

    # (sleeping doesn't seem to really have an effect on how many calls we can get)
    # 1st key: (sleep=0.2) -319
    # 2nd key: (sleep=0.2) 319-653
    # 3rd key: (no sleep) 653-989
    #time.sleep(0.2)

    results = Results(cs, cs.filter_inchi, (inchi_with_prefix,), raise_errors=True)
    results.wait()
    assert results.status == 'Complete'

    if len(results) == 0:
        return None

    csid = min([r.record_id for r in results])
    return csid


def inchi2chemspider_url(inchi, **kwargs):
    # NOTE: this is @cached, so no need to cache inchi2chemspider_url
    csid = inchi2chemspider_id(inchi, **kwargs)
    if csid is None:
        return None

    return chemspider_url(csid)


_warned_no_chemspider_key = False

@cached
def inchi2density(inchi, verbose=False):
    """Returns mean PubChem density [g/mL] of the compound (else ChemSpider's, if
    `chemspipy_api_key` is set).

    PubChem's relative densities (to water, as most of its values are) are taken as g/mL
    (water is 0.997-1.000 g/mL from 4-25 C). Densities of the vapor, of air saturated
    with it, of solutions and bulk densities are not used, nor are strings that cannot be
    parsed (with a warning).
    """
    global _warned_no_chemspider_key

    null_ret = _property_series('density')
    if isnull(inchi):
        return null_ret

    pubchem_dict = pubchemprops_lookup(inchi, 'Density',
        expected_units='density', skip_unit_parse_err=True, verbose=verbose
    )
    if pubchem_dict is not None:
        return _property_series('density', pubchem_dict['quantity'],
            pubchem_dict['sources'], pubchem_dict['units']
        )

    if chemspipy_api_key is None:
        if not _warned_no_chemspider_key:
            warn('could not lookup densities on chemspider b/c no chemspipy_api_key set '
                '(only warning once)'
            )
            _warned_no_chemspider_key = True
        return null_ret

    # (lowest compound ID, if there are several)
    csid = inchi2chemspider_id(inchi)

    if csid is None:
        if verbose:
            print('COULD NOT FIND INCHI IN CHEMSPIDER')

        return null_ret

    # Since I couldn't find a way to get arbitrary properties like
    # "Experimental Density" through the API.
    cs_dict = chemspider_experimental_density(csid, verbose=verbose)
    if verbose:
        if cs_dict is not None:
            print('CHEMSPIDER DENSITY:', cs_dict['density'])
        else:
            print('NO CHEMSPIDER DENSITY')

    if cs_dict is None:
        return null_ret
    else:
        return cs_dict


def _parse_nist_number(value) -> float:
    """Returns number in NIST WebBook table (e.g. '6.0×10+8'), or NaN if an inequality.
    """
    if not isinstance(value, str):
        return float(value)

    value = value.strip().replace('−', '-')
    if value.startswith(('<', '>')):
        return np.nan

    if '×10' in value:
        mantissa, power = value.split('×10')
        return float(mantissa) * 10**int(power)

    return float(value)


nist_k_henry_col = 'k°H (mol/(kg*bar))'

@cached
def inchi2nist_k_henry(inchi, verbose=False):
    """Returns mean Henry's law solubility constant in water from the NIST WebBook.

    Values (in mol / kg / bar, at 298.15 K) are from the "Henry's Law data" section of
    the compound's page (https://webbook.nist.gov). Values given as inequalities are
    not used.
    """
    null_ret = _property_series('k_henry')
    if inchi is None:
        return null_ret

    url = inchi2nist_henrys_law_url(inchi)
    if verbose:
        print(url)

    response = session.get(url)
    response.raise_for_status()

    try:
        # (pandas>=3 requires literal HTML be wrapped in StringIO)
        # (flavor='lxml', as with the default, pages without the table make pandas try
        # another parser, which needs html5lib)
        dfs = pd.read_html(StringIO(response.text), flavor='lxml',
            attrs={'aria-label': "Henry's Law constant (water solution)"}
        )
    except ValueError as err:
        if 'No tables found' not in str(err):
            raise
        return null_ret

    assert len(dfs) == 1
    df = dfs[0]

    kh = df[nist_k_henry_col].map(_parse_nist_number).mean()
    if pd.isnull(kh):
        return null_ret

    return _property_series('k_henry', kh, url, property2preferred_units['k_henry'])


# (named for where henrys-law.org used to be hosted. names kept, as they are also cache
# keys)
mp_henry_prefix = 'https://henrys-law.org/henry/'
def inchikey2maxplanck_k_henry_url(inchikey):
    return f'{mp_henry_prefix}search_identifier.html?search={inchikey}'


# Types of values on henrys-law.org (Sander, 2015: https://doi.org/10.5194/acp-15-4399-2015)
henrys_law_org_types = {
    'L': 'literature review',
    'M': 'measured',
    'V': 'VP/AS = vapor pressure/aqueous solubility',
    'R': 'recalculation',
    'T': 'thermodynamical calculation',
    'X': 'original paper not available',
    'C': 'citation',
    'Q': 'QSPR',
    'E': 'estimate',
    '?': 'unknown',
    'W': 'wrong',
}

def _parse_henrys_law_org_number(value) -> float:
    """Returns number in henrys-law.org table (e.g. '1.2×10−3' for 1.2e-3).
    """
    try:
        return float(value)

    except ValueError:
        # float(...) doesn't work with strings containing the more exotic left character
        value = value.replace('−', '-')

        exp_prefix = '×10'
        assert exp_prefix in value
        # Strings where the exponent is x10^1 should appear here as 'x101'
        # The original HTML has <sup>...</sup> tags around exponent.
        assert not value.endswith(exp_prefix)

        mantissa, power = value.split(exp_prefix)
        return float(f'{float(mantissa)}e{int(power)}')


@cached
def inchi2maxplanck_k_henry(inchi, verbose=False):
    """Returns Henry's law solubility constant in water from https://henrys-law.org

    The site is Rolf Sander's compilation (it was hosted at satellite.mpic.de/henry).
    Values (in mol / m^3 / Pa) are at 298.15 K, and data from temperatures above 373 K
    were not used. The first (most reliable) value for the compound is used, other than
    those the site marks as wrong, and those given as inequalities (except Altschuh et
    al. (1999) values, given as lower bounds, which are used as values).
    """
    null_ret = _property_series('k_henry')
    if inchi is None:
        return null_ret

    # Returns None on error
    inchikey = InchiToInchiKey(_inchi_with_prefix(inchi))
    assert inchikey is not None

    search_url = inchikey2maxplanck_k_henry_url(inchikey)

    response = session.get(search_url)
    html = response.text

    # This parser is case insensitive (useful here).
    # Requires: pip install lxml
    soup = BeautifulSoup(html, 'xml')
    tags = soup.find_all('h1')
    # Website must have changed if this happens.
    assert len(tags) == 1
    tag = tags[0]
    n_results = int(tag.text.split()[0])
    assert n_results >= 0
    if n_results == 0:
        return null_ret

    if n_results > 1:
        raise ValueError('>1 results for this inchikey')

    # not sure why it requires two next_sibling calls rather than 1...
    # from looking at source it seems like maybe it should be 1
    table = tag.next_sibling.next_sibling
    # assuming all links in table would go to same place
    # (and that there are links) (both assumptions seem true for now)
    url_suffix = '/'.join(table.find('a')['href'].split('/')[-2:])
    url = mp_henry_prefix + url_suffix

    if verbose:
        print(url)

    response = session.get(url)
    html = response.text

    # (pandas>=3 requires literal HTML be wrapped in StringIO)
    dfs = pd.read_html(StringIO(html), flavor='lxml')

    # Looking for one DataFrame with an index like this:
    # [MultiIndex([(                  'Hscp',       '[mol/(m3Pa)]'),
    #             ('d ln Hs  cp /  d (1/T)',                '[K]'),
    #             (            'References', 'Unnamed: 2_level_1'),
    #             (                  'Type', 'Unnamed: 3_level_1'),
    #             (                 'Notes', 'Unnamed: 4_level_1')],
    #            )]
    target_cols = ['Hscp', 'd ln Hs  cp /  d (1/T)', 'References', 'Type', 'Notes']

    dfs = [x for x in dfs
        if x.shape[1] == len(target_cols) and
        (x.columns.get_level_values(0) == target_cols).all()
    ]
    assert len(dfs) == 1
    df = dfs[0]
    assert len(df) > 0

    # Remainder in this level should be empty in HTML table (-> 'Unnamed: <x> ' here)
    assert (df.columns.get_level_values(1)[:2] == ['[mol/(m3Pa)]', '[K]']).all()
    # pint seems to parse this succcessfully
    units = 'mol / m^3 / Pa'

    df = df.droplevel(1, axis='columns')

    type_set = set(df.Type)
    assert all(x in henrys_law_org_types.keys() for x in type_set), \
        f'unknown types in {type_set}'

    # Rows are in decreasing order of reliability.
    for _, row in df[df.Type != 'W'].iterrows():
        Hscp = row.Hscp
        short_ref = row.References.replace('\xa0', ' ')

        # For Altschuh et al. (1999) data (https://doi.org/10.1016/S0045-6535(99)00082-X),
        # other values are often lower despite the '>', and of other types and close,
        # e.g. for https://henrys-law.org/henry/casrn/60-12-8,
        # https://henrys-law.org/henry/casrn/3360-41-6, and
        # https://henrys-law.org/henry/casrn/122-97-4
        if type(Hscp) is str and ('<' in Hscp or '>' in Hscp):
            gt_ref = 'Altschuh et al. (1999)'
            if short_ref == gt_ref and '<' not in Hscp:
                warn(f"dropping '>' from {gt_ref} value='{Hscp}'\n{url}\n")
                Hscp = Hscp.strip('> ')
            else:
                if verbose:
                    print(f"not using data with inequality: '{Hscp}'")
                continue

        source = f'{short_ref} ({row.Type})'
        notes = row.at['Notes']
        if not pd.isna(notes):
            source += f'\nnotes: {notes}'
        source += f'\n{url}'

        return _property_series('k_henry', _parse_henrys_law_org_number(Hscp), source,
            units
        )

    warn(f'not using henrys-law.org data, as all is wrong or inequalities\n{url}\n')
    return null_ret


def inchi2k_henry(inchi, verbose=False, **kwargs):
    """Returns Henry's law solubility constant in water, from henrys-law.org or NIST.

    The henrys-law.org value (from `inchi2maxplanck_k_henry`, in mol / m^3 / Pa) is used
    if there is one, otherwise the NIST WebBook's (from `inchi2nist_k_henry`, in mol / kg
    / bar). Which was used can be seen from 'k_henry_sources' (and 'k_henry_units').
    `kwargs` are passed to both lookups. See `inchi2est_k_henry` for an estimate when
    neither has a value.
    """
    keys = ['k_henry', 'k_henry_sources', 'k_henry_units']

    for fn in (inchi2maxplanck_k_henry, inchi2nist_k_henry):

        ret = fn(inchi, verbose=verbose, **kwargs)
        assert set(ret.keys()) == set(keys)
        assert is_series(ret)

        # Taking first non-null value (so order of fns about is priority).
        if pd.notna(ret['k_henry']):
            break

    # Last call should give us null Series in this case
    return ret


def inchi2est_k_henry(inchi, verbose=False) -> pd.Series:
    """Returns Henry's law solubility constant in water, from solubility / vapor pressure.

    Estimated with the vapor pressure / aqueous solubility (VP/AS) method: C_sat / p_sat
    [mol / m^3 / Pa], with C_sat the solubility in water (`inchi2water_solubility`) and
    p_sat the vapor pressure of the pure compound (`inchi2vapor_pressure`). This assumes
    the compound's activity coefficient in water is the same at saturation as at
    infinite dilution, which is only reasonable for compounds that are not very soluble,
    and is meaningless for those miscible with water. The PubChem solubilities are not
    checked to be at any particular temperature. For toluene and hexane, estimates were
    within 1% of `inchi2k_henry` values, and 1.6x and 2.5x those for benzaldehyde and
    linalool (~1.6 g/L in water).

    Returns:
        Series with same keys as `inchi2k_henry`. Values are null if either the
        solubility or the vapor pressure are not available.
    """
    null_ret = _property_series('k_henry')
    if isnull(inchi):
        return null_ret

    solubility = inchi2water_solubility(inchi, verbose=verbose)
    vapor_pressure = inchi2vapor_pressure(inchi, verbose=verbose)
    if isnull(solubility['water_solubility']) or isnull(vapor_pressure['vapor_pressure']):
        return null_ret

    molar_mass = inchi2molar_mass(inchi)
    c_sat = ureg.Quantity(solubility['water_solubility'],
        solubility['water_solubility_units']
    ) / ureg.Quantity(molar_mass['molar_mass'], molar_mass['molar_mass_units'])
    p_sat = ureg.Quantity(vapor_pressure['vapor_pressure'],
        vapor_pressure['vapor_pressure_units']
    )
    units = 'mol / m^3 / Pa'
    return _property_series('k_henry', (c_sat / p_sat).to(units).magnitude,
        'VP/AS estimate from water solubility '
        f'({solubility["water_solubility_sources"]}) and vapor pressure '
        f'({vapor_pressure["vapor_pressure_sources"]})', units
    )


gas_constant = 8.314462618  # J / (mol K)
# n-alkanes used as surrogates for paraffin oil (pfo) in inchi2k_henry_pfo, by CAS
pfo_surrogate_cas = {'hexadecane': '544-76-3', 'squalane': '111-01-3'}
# thermo vapor pressure methods that are estimates (rather than fit to measurements)
_thermo_estimated_psat_methods = {'BOILING_CRITICAL', 'LEE_KESLER_PSAT', 'AMBROSE_WALTON',
    'SANJARI', 'EDALAT', 'EOS'
}
_thermo_estimated_tm_methods = {'JOBACK'}

pfo_unifac_method = 'UNIFAC estimate (no measured solubility in oil or alkanes used)'
pfo_log_k_method = 'measured log10 K oil/air (passed in)'
pfo_abraham_l_method = 'Abraham L (passed in; log10 K hexadecane/air at 25 C)'

# (temp_k, solvent) -> (molar density [mol/m^3], UNIFAC groups)
_pfo_surrogates = dict()

def _pfo_surrogate(solvent: str, temp_k: float) -> tuple:
    """Returns molar density [mol / m^3] and UNIFAC groups of oil surrogate `solvent`.
    """
    from chemicals.identifiers import search_chemical
    from thermo import Chemical
    from thermo import unifac

    key = (solvent, temp_k)
    if key not in _pfo_surrogates:
        cas = pfo_surrogate_cas[solvent]
        unifac.load_group_assignments_DDBST()
        groups = dict(unifac.DDBST_MODIFIED_UNIFAC_assignments[
            search_chemical(cas).InChI_key
        ])
        molar_volume = Chemical(cas, T=temp_k).VolumeLiquid(temp_k, 101325)
        _pfo_surrogates[key] = (1 / molar_volume, groups)

    return _pfo_surrogates[key]


def _thermo_search(inchikey: str):
    """Returns `chemicals` metadata for compound with `inchikey`, or None if not found.
    """
    from chemicals.identifiers import search_chemical
    try:
        return search_chemical(f'InChIKey={inchikey}')
    except ValueError:
        return None


def _unifac_groups(inchikey: str, smiles: str) -> tuple:
    """Returns modified UNIFAC (Dortmund) groups and where they are from.

    Groups are DDBST's published assignments (bundled with thermo) where available,
    otherwise from ugropy. Groups are None if neither has any.
    """
    from thermo import unifac

    unifac.load_group_assignments_DDBST()
    table = unifac.DDBST_MODIFIED_UNIFAC_assignments
    if inchikey in table:
        return dict(table[inchikey]), 'DDBST'

    # Same skeleton (e.g. other stereoisomers). Groups don't depend on stereochemistry.
    prefix = inchikey.split('-')[0] + '-'
    matches = sorted(k for k in table if k.startswith(prefix))
    if len(matches) > 0:
        match = next((k for k in matches if k.endswith('-UHFFFAOYSA-N')), matches[0])
        return dict(table[match]), f'DDBST (for {match}, same skeleton)'

    from ugropy import dortmund

    try:
        groups = dortmund.get_groups(smiles, 'smiles').subgroups_num
    # (ugropy can fail many ways on unusual structures)
    except Exception as err:
        return None, f'ugropy failed ({err})'

    if not groups:
        return None, 'ugropy found no complete Dortmund fragmentation'

    return {int(k): int(v) for k, v in groups.items()}, 'ugropy (automatic)'


def _pfo_vapor_pressure(inchi, chemical, temp_k: float, solid: Optional[bool]
    ) -> tuple:
    """Returns liquid vapor pressure [Pa] (or None), source, and whether measured.

    Whether measured can be None if not known.

    For compounds solid at `temp_k`, this is the subcooled liquid vapor pressure, so
    PubChem values (which would be of the solid) are not used for them.
    """
    psat = None
    estimated_psat = None
    if chemical is not None:
        try:
            psat = chemical.VaporPressure(temp_k)
        # (thermo raises various errors outside a method's valid range)
        except Exception:
            psat = None

        if psat is not None and psat > 0:
            method = chemical.VaporPressure.method
            source = f'thermo {method}'
            if method not in _thermo_estimated_psat_methods:
                return psat, source, True
            estimated_psat = (psat, source)

    temp_c = temp_k - 273.15
    if abs(temp_c - vapor_pressure_temp_c) <= 1 and not solid:
        ret = inchi2vapor_pressure(inchi)
        if not isnull(ret['vapor_pressure']):
            psat = ureg.Quantity(ret['vapor_pressure'], ret['vapor_pressure_units']
                ).to('Pa').magnitude
            # (None if cached before this was added)
            measured = ret.get('vapor_pressure_measured')
            return psat, f'PubChem ({ret["vapor_pressure_sources"]})', measured

    if estimated_psat is not None:
        return (*estimated_psat, False)

    ret = inchi2est_vapor_pressure(inchi, temp_c=temp_c)
    if not isnull(ret['vapor_pressure']):
        psat = ureg.Quantity(ret['vapor_pressure'], ret['vapor_pressure_units']
            ).to('Pa').magnitude
        return psat, ret['vapor_pressure_sources'], False

    return None, None, None


def inchi2k_henry_pfo(inchi, temp_c: float = 25.0, solvent: str = 'hexadecane',
    log_k_oil_air: Optional[float] = None, abraham_l: Optional[float] = None,
    vapor_pressure_kpa: Optional[float] = None, verbose: bool = False) -> pd.Series:
    """Returns estimated Henry's law solubility constant in paraffin oil (pfo).

    Estimated as in oil_headspace.py: the oil/air partition coefficient
    K = C_oil / C_air = c_solv * R * T / (gamma_inf * p_sat), where gamma_inf is the
    infinite-dilution activity coefficient in an n-alkane surrogate for the oil, from
    modified UNIFAC (Dortmund), c_solv is the surrogate's molar density, and p_sat is the
    solute's (subcooled) liquid vapor pressure. The returned constant is
    C_oil / p = K / (R T), in mol / m^3 / Pa. This is the same definition (and units) as
    `inchi2k_henry` values from henrys-law.org (for water, the solvent there). Its NIST
    values are in mol / kg / bar (~100x larger in water than mol / m^3 / Pa).

    No measured solubilities or partition coefficients in oil or alkanes (e.g. in
    squalane or hexadecane) are used, unless passed in (see `log_k_oil_air` and
    `abraham_l`). gamma_inf is predicted by UNIFAC from the solute's functional groups
    alone ('squalane' and 'hexadecane' are UNIFAC models of those solvents, not data for
    them). UNIFAC's group interaction parameters were fit to measurements on many
    compounds (possibly including the solute), but no measurements for the particular
    solute are used. The 'k_henry_pfo_method' and 'k_henry_pfo_measured_partitioning'
    outputs say which was used for each compound.

    For simple functional groups, UNIFAC K is typically within ~2x of measured values
    (for 10 simple compounds with measured Abraham L, hexadecane log10 K estimates were
    within 0.26 of it, median 0.06). Carboxylic acids (which dimerize) are less certain, and
    solutes that self-associate (e.g. alcohols) are less soluble than predicted at
    higher concentrations. Vapor pressures estimated with SIMPOL.1 can be biased for some
    classes (see `inchi2est_vapor_pressure`), which biases K by the same factor.

    UNIFAC groups are DDBST's published assignments (bundled with thermo) where
    available, otherwise from ugropy.

    Args:
        temp_c: temperature [C].

        solvent: n-alkane surrogate for the oil ('hexadecane' or 'squalane'). The two
            typically give log10 K within ~0.2 of each other.

        log_k_oil_air: measured log10 K (C_oil / C_air) for paraffin oil at `temp_c`. If
            passed, it is used instead of the UNIFAC estimate.

        abraham_l: Abraham L solute descriptor (log10 K hexadecane/air at 25 C, e.g. from
            the UFZ-LSER database), used instead of the UNIFAC estimate if passed (and
            `log_k_oil_air` is not). Not corrected for temperature.

        vapor_pressure_kpa: liquid (or subcooled liquid) vapor pressure at `temp_c`, for
            the UNIFAC estimate. Otherwise, in order of preference: from thermo's fits to
            measurements, from PubChem (`inchi2vapor_pressure`, only at 25 C, and not
            for compounds thermo has as solid), estimated by thermo, or estimated by
            SIMPOL.1 (`inchi2est_vapor_pressure`). PubChem values for compounds that are
            solid (but that thermo has no melting point for) would be of the solid, and
            would make K too high.

        verbose: if True, prints the returned values.

    Returns:
        Series with:
        - 'k_henry_pfo', 'k_henry_pfo_sources', 'k_henry_pfo_units': as for
          `inchi2k_henry` (with sources describing everything the value depends on).
        - 'k_henry_pfo_method': one of `pfo_unifac_method`, `pfo_log_k_method`, or
          `pfo_abraham_l_method`.
        - 'k_henry_pfo_measured_partitioning': whether measured partitioning between
          oil (or an alkane) and air was used (only if passed in).
        - 'log10_k_oil_air': log10 of the oil/air partition coefficient K.
        - 'solvent_surrogate', 'gamma_inf', 'unifac_groups', 'unifac_groups_source',
          'unifac_parameters': for UNIFAC estimates.
        - 'vapor_pressure_kpa', 'vapor_pressure_source', 'vapor_pressure_measured': the
          vapor pressure used for UNIFAC estimates, and whether it is from measurements
          (None if passed in, or not known, as for PubChem values cached before
          `inchi2vapor_pressure` checked).
        - 'notes': caveats for this compound.

        Values are null (with the reason in 'notes') if K could not be estimated.
    """
    # (imported here, as importing thermo takes ~0.5s)
    from thermo import Chemical
    from thermo import unifac

    temp_k = temp_c + 273.15
    ret = {k: None for k in ('k_henry_pfo', 'k_henry_pfo_sources', 'k_henry_pfo_units',
        'k_henry_pfo_method', 'k_henry_pfo_measured_partitioning', 'log10_k_oil_air',
        'solvent_surrogate', 'gamma_inf', 'unifac_groups', 'unifac_groups_source',
        'unifac_parameters', 'vapor_pressure_kpa', 'vapor_pressure_source',
        'vapor_pressure_measured', 'notes'
    )}
    notes = []

    def output():
        ret['notes'] = '; '.join(notes) if notes else None
        series = pd.Series(ret)
        if verbose:
            print(series.to_string())
        return series

    if isnull(inchi):
        return output()

    if log_k_oil_air is not None:
        log_k = log_k_oil_air
        ret['k_henry_pfo_method'] = pfo_log_k_method
        sources = f'{pfo_log_k_method}: {log_k_oil_air:.3g}'

    elif abraham_l is not None:
        log_k = abraham_l
        ret['k_henry_pfo_method'] = pfo_abraham_l_method
        sources = f'{pfo_abraham_l_method}: {abraham_l:.3g}'
        if abs(temp_c - 25) > 0.5:
            notes.append('Abraham L is for 25 C (not corrected for temperature)')

    else:
        ret['k_henry_pfo_method'] = pfo_unifac_method
        ret['solvent_surrogate'] = solvent

        mol = Chem.MolFromInchi(_inchi_with_prefix(inchi))
        if mol is None:
            notes.append('could not parse InChI')
            return output()
        inchikey = Chem.MolToInchiKey(mol)

        meta = _thermo_search(inchikey)
        if meta is None:
            # (e.g. stereoisomers thermo doesn't have, with the same properties)
            meta = _thermo_search(inchikey.split('-')[0] + '-UHFFFAOYSA-N')
            if meta is not None:
                notes.append(f'thermo properties from {meta.common_name} (same '
                    'skeleton, without stereochemistry)'
                )

        chemical = None
        solid = None
        if meta is not None:
            try:
                chemical = Chemical(meta.CASs, T=temp_k)
            # (thermo can fail to load properties many ways)
            except Exception as err:
                notes.append(f'thermo could not load properties ({err})')

        if chemical is not None and chemical.Tm is not None and (
            chemical.Tm_source not in _thermo_estimated_tm_methods):
            solid = temp_k < chemical.Tm
            if solid:
                notes.append(f'thermo melting point is {chemical.Tm - 273.15:.0f} C '
                    f'({chemical.Tm_source}), so solid at {temp_c:g} C: K only applies '
                    'below its solubility in the oil'
                )

        if vapor_pressure_kpa is not None:
            psat = vapor_pressure_kpa * 1000
            ret['vapor_pressure_source'] = 'passed in'
        else:
            psat, ret['vapor_pressure_source'], ret['vapor_pressure_measured'] = (
                _pfo_vapor_pressure(inchi, chemical, temp_k, solid)
            )

        if psat is None:
            notes.append('no vapor pressure')
            return output()
        ret['vapor_pressure_kpa'] = psat / 1000

        groups, ret['unifac_groups_source'] = _unifac_groups(inchikey,
            Chem.MolToSmiles(mol)
        )
        if groups is None:
            notes.append(f'no UNIFAC groups: {ret["unifac_groups_source"]}')
            return output()
        ret['unifac_groups'] = ' '.join(f'{unifac.DOUFSG[k].group}:{n}'
            for k, n in sorted(groups.items())
        )

        solvent_conc, solvent_groups = _pfo_surrogate(solvent, temp_k)

        parameter_sets = {'Dortmund 2016': unifac.DOUFIP2016,
            'modified UNIFAC 2.0': unifac.DOUF2IP
        }
        main_groups = sorted({unifac.DOUFSG[k].main_group_id
            for g in (groups, solvent_groups) for k in g
        })
        for name, params in parameter_sets.items():
            if all(b in params.get(a, {}) for a in main_groups for b in main_groups
                if a != b):
                ret['unifac_parameters'] = name
                break
        else:
            notes.append('no UNIFAC interaction parameters for some groups')
            return output()

        model = unifac.UNIFAC.from_subgroups(T=temp_k, xs=[0.5, 0.5],
            chemgroups=[groups, solvent_groups], version=1,
            interaction_data=params, subgroups=unifac.DOUFSG
        )
        ret['gamma_inf'] = float(model.gammas_infinite_dilution()[0])
        log_k = np.log10(solvent_conc * gas_constant * temp_k / (ret['gamma_inf'] * psat))

        main_group_names = {unifac.DOUFSG[k].main_group for k in groups}
        if main_group_names & {'COOH', 'HCOOH'}:
            notes.append('carboxylic acid: UNIFAC ignores acid dimerization, so K is '
                'uncertain'
            )
        elif main_group_names & {'OH', 'CH3OH', 'ACOH'}:
            notes.append('OH group: self-association lowers solubility at higher '
                'concentrations'
            )

        if ret['vapor_pressure_source'] == 'passed in':
            measured = 'passed in'
        else:
            measured = {True: 'measured', False: 'not measured', None:
                'not known whether measured'
            }[ret['vapor_pressure_measured']]
        sources = (f'{pfo_unifac_method}: gamma_inf in {solvent} (modified UNIFAC, '
            f'{ret["unifac_parameters"]} parameters, groups from '
            f'{ret["unifac_groups_source"]}), vapor pressure from '
            f'{ret["vapor_pressure_source"]} ({measured})'
        )

    ret['k_henry_pfo_measured_partitioning'] = (
        ret['k_henry_pfo_method'] != pfo_unifac_method
    )
    ret['log10_k_oil_air'] = log_k
    ret['k_henry_pfo'] = 10**log_k / (gas_constant * temp_k)
    ret['k_henry_pfo_units'] = 'mol / m^3 / Pa'
    ret['k_henry_pfo_sources'] = sources
    return output()


def add_properties(df, props=None, verbose=False):
    """Adds chemicals properties to DataFrame `df` with `inchi` in columns.
    """
    # Checking it's a DataFrame
    assert len(df.shape) == 2
    assert 'inchi' in df.columns

    if props is None:
        props = properties

    if verbose:
        print('Trying to add properties:', props)

    for prop in props:
        new_cols = [prop, f'{prop}_sources', f'{prop}_units']
        inchi2prop_fn_name = f'inchi2{prop}'
        if verbose:
            print(f'Using function {inchi2prop_fn_name} for lookup')

        # This will fail if the function is not defined.
        inchi2prop_fn = globals()[inchi2prop_fn_name]

        # Since `cached` wrapper currently will make `allowed_kwarg` return
        # True for all inputs.
        if hasattr(inchi2prop_fn, '__wrapped__'):
            kwarg_check_fn = inchi2prop_fn.__wrapped__
        else:
            kwarg_check_fn = inchi2prop_fn

        # This will probably fail if either kwargs are not accepted,
        # or columns returned do not match.
        df[new_cols] = df.inchi.apply(inchi2prop_fn, verbose=verbose)

    return df


# https://stackoverflow.com/questions/19153462
_LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
def index2excel_col_name(i):
    """ Convert given column number (0-indexed) to an Excel-style cell name. """
    i = i + 1
    result = []
    while i:
        i, rem = divmod(i - 1, 26)
        result[:0] = _LETTERS[rem]
    return ''.join(result)


def load_or_save_excel_properties(df, excel_fname, name='', props=None,
    overwrite_excel=False, verbose=False, all_urls=True):
    """Saves `df` with properties to `excel_fname` for missing value entry.

    All properties are looked up before the spreadsheet is written.

    If `excel_fname` already exists, it is loaded instead of being written to,
    and any new values are added to the cache, for all future lookups.

    `df` is assumed to already have ['name','inchi'] among its columns.

    `df` argument is ignored if spreadsheet will be loaded.
    """
    if props is None:
        props = properties
    else:
        for p in props:
            if p not in properties:
                raise ValueError(f'props must all be in {properties}')

    if name != '':
        name += ' '

    if overwrite_excel or not exists(excel_fname):
        assert 'name' in df.columns and 'inchi' in df.columns

        if 'pubchem_url' not in df.columns:
            df['pubchem_url'] = df.inchi.map(inchi2pubchem_url)

        if all_urls:
            # (as ChemSpider can only be searched with an API key)
            if 'chemspider_url' not in df.columns and chemspipy_api_key is not None:
                df['chemspider_url'] = df.inchi.map(lambda x:
                    inchi2chemspider_url(x, ignore_cache_null=True,
                    leave_nonnull_cache=True)
                )

            if 'nist_webbook_url' not in df.columns:
                df['nist_webbook_url'] = df.inchi.map(inchi2nist_webbook_url)

        add_properties(df, props, verbose=verbose)

        property_cols = [c for p in props for c in
            [p, f'{p}_units', f'{p}_sources']
        ]
        url_cols = [c for c in df.columns if c.endswith('_url')]
        col_order = ['name'] + property_cols + url_cols + ['inchi']
        other_cols = [c for c in df.columns if c not in col_order]
        col_order += other_cols

        df = df[col_order].copy()

        print(f'Writing {name}to {excel_fname} for manual entry of missing '
            'properties.'
        )
        writer = pd.ExcelWriter(excel_fname, engine='xlsxwriter')
        df.to_excel(writer, index=False)

        # A dict of sheet name -> Worksheet objects.
        assert len(writer.sheets) == 1
        # Default sheet name is Sheet1
        worksheet = writer.sheets['Sheet1']

        # Property values from lookups are gray, so they can be distinguished from those
        # entered manually (in empty cells, which have the column format).
        number_format = {'num_format': '0.000'}
        manual_format = writer.book.add_format(number_format)
        lookup_format = writer.book.add_format({**number_format, 'font_color': 'gray'})
        for i, col in enumerate(df.columns):
            if col in props:
                for j, value in enumerate(df[col]):
                    if pd.notnull(value):
                        worksheet.write_number(j + 1, i, value, lookup_format)

        fixed_col_widths = {
            'chemspider_url': len('chemspider_url'),
            'nist_webbook_url': len('nist_webbook_url'),
        }
        # (using length of str repr of values as proxy for width, which may be off for
        # floats)
        max_col_widths = df.map(lambda x: len(str(x))).max(axis=0)
        for i, w in enumerate(max_col_widths):
            n = df.columns[i]
            if n in fixed_col_widths:
                w = fixed_col_widths[n]
            else:
                # So column names are also considered. Could maybe do once
                # before loop...
                w = max(w, len(n))

            c = index2excel_col_name(i)
            r = f'{c}:{c}'
            # I thought integers were supposed to work for first argument,
            # but they did not.
            worksheet.set_column(r, w, manual_format if n in props else None)

        writer.close()

    else:
        print(f'Reading {name}property data from', excel_fname)
        df = pd.read_excel(excel_fname, engine=pandas_excel_engine)
        add_properties_to_cache(df, verbose=verbose)

    return df


def add_properties_to_cache(df, overwrite_close=False, dropna=True,
    verbose=False):
    """Adds chemical properties in `df` to cache under 'inchi' (marked as manual).

    Values must have units and sources (not just whitespace), and are converted to
    units in `property2preferred_units`. Values close to those already cached are not
    added.
    """
    df = df.assign(inchi=df.inchi.map(normalize_inchi))

    # So we don't have to check whether values behind each instance
    # of same inchi are consistent.
    assert not df.inchi.duplicated().any(), 'inchi can not have duplicates'

    df = df.set_index('inchi')

    add_count = 0
    for prop in properties:
        if prop in df.columns:
            preferred_unit_str = property2preferred_units[prop]
            preferred_units = ureg(preferred_unit_str)

            def add_to_cache(i, p, row):
                nonlocal add_count

                sources = getattr(row, f'{prop}_sources')
                if pd.isnull(sources) or len(str(sources).strip()) == 0:
                    raise ValueError(f'{row.Index} {prop} missing source')

                new_entry = _property_series(prop, p, sources, preferred_unit_str)
                _cache_set('inchi', prop, i, new_entry, manual=True)
                add_count += 1

            if verbose:
                print(f'adding {prop}')

            cols = [prop, f'{prop}_sources', f'{prop}_units']
            assert all([c in df.columns for c in cols[1:]])

            for row in df[cols].dropna(subset=[prop]).itertuples():
                # Since usual dict key indexing doesn't work w/ itertuples
                # objects.
                unit_str = getattr(row, f'{prop}_units')
                if pd.isnull(unit_str):
                    raise ValueError(f'{row.Index} {prop} missing units')

                if not preferred_units.check(unit_str):
                    raise ValueError(f'{row.Index} {prop} units of {unit_str} '
                        f'inconsistent with {preferred_unit_str}'
                    )

                m = getattr(row, prop)
                try:
                    float(m)
                except ValueError:
                    raise ValueError(f'{row.Index} {prop} not a number')

                q = (m * ureg(unit_str)).to(preferred_unit_str).magnitude

                inchi = row.Index
                cached_value = _cache_get('inchi', prop, inchi)
                if cached_value is _missing:
                    if verbose:
                        print(f'adding {inchi} not already in cache')
                    add_to_cache(inchi, q, row)
                    continue

                if isnull(cached_value) or pd.isnull(cached_value[prop]):
                    if verbose:
                        print(f'adding {inchi} null in cache')
                    add_to_cache(inchi, q, row)
                    continue

                if not np.isclose(cached_value[prop], q):
                    if verbose:
                        print(f'adding {inchi} b/c not close to cache value')
                    add_to_cache(inchi, q, row)
                else:
                    if verbose:
                        print(f'not adding {inchi} b/c close to cache value')

    if add_count == 0:
        print('Nothing new found to add to cache.')


_odor_inventory_gsheet = None
def odor_inventory_gsheet(use_cache=False, verbose=False):
    """Returns a DataFrame with data odor inventory data from Google sheet.
    """
    global _odor_inventory_gsheet
    if _odor_inventory_gsheet is not None:
        return _odor_inventory_gsheet

    gsheet_cache_file = '.odor_inventory_gsheet_cache.p'
    if use_cache and exists(gsheet_cache_file):
        print(f'loading odor inventory sheet data from cache at {gsheet_cache_file}')
        df = pd.read_pickle(gsheet_cache_file)
    else:
        # Falling back to os.getcwd() prefix as I don't have the data installed
        # correctly with conda right now. Might want to allow setting via env
        # var or something like that too.
        odor_inventory_link_fname = 'odor_inventory_sheet_link.txt'
        dirs_to_try = (os.getcwd(), pkg_data_dir)
        paths_to_try = [join(d, odor_inventory_link_fname) for d in dirs_to_try]
        link_filename = None
        for p in paths_to_try:
            if exists(p):
                link_filename = p
                break

        if link_filename is None:
            raise IOError(f'{odor_inventory_link_fname} not found in any of '
                f'{dirs_to_try}'
            )

        with open(link_filename, 'r') as f:
            gsheet_url = \
                f.readline().split('/edit')[0] + '/export?format=csv&gid='

        gid = '0'
        df = pd.read_csv(gsheet_url + gid)

        df.dropna(how='all', inplace=True)
        df.to_pickle(gsheet_cache_file)

    df.drop(columns=['Quantity', 'Recieved', 'Purity', 'Aliquots', 'Notes'],
        inplace=True
    )

    df.rename(columns={
        'Chemical': 'name',
        'CAS #': 'cas',
        'Ionization potential (eV)': 'ionization_v',
        'Abbreviation': 'abbrev',
        'InChI Key': 'inchikey'
    }, inplace=True)
    df.rename(columns=lambda s: s.lower().replace(' ','_'), inplace=True)

    df.abbrev = df.abbrev.apply(lambda s: s.strip() if pd.notnull(s) else s)

    df['original_name'] = df.name.copy()

    df.name = df.name.apply(normalize_name)

    df = convert(df, to_type='inchi', allow_nan=False, verbose=verbose,
        report_conflicts=False
    )

    _odor_inventory_gsheet = df
    return df


def inchi2abbrev_dict(df, allow_orphan_abbrevs=False):
    abbrev_notnull = df.abbrev.notnull()
    abbrevs_without_inchi = abbrev_notnull & df.inchi.isnull()
    if abbrevs_without_inchi.any():
        print('Abbreviations without inchi:')
        print(df[abbrevs_without_inchi])
        print('')
        if not allow_orphan_abbrevs:
            raise ValueError('all abbreviations must have inchi')

    return dict(zip(df.inchi[abbrev_notnull], df.abbrev[abbrev_notnull]))


def odor_is_mix(odor_name):
    """Returns True if odor is the mixture in my complex mix expts, else False.
    For ordering odors in plots.
    """
    lo = odor_name.lower()
    return 'approx' in lo or 'mix' in lo


hardcoded_odor2abbrevs = {
    'paraffin': 'pfo',
    'd3 kiwi': 'kiwi',
    'fly food b': 'fly food'
}
hardcoded_not_to_abbrev = {
    'pfo',
    'water'
}
_inchi2abbrev = None
def odor2abbrev(odor_name, *args, allow_orphan_abbrevs=False,
    skip=set(), **kwargs):
    """
    Uses abbreviation column in odor inventory gsheet to lookup abbreviation for
    arbitrary chemical names (doesn't have to match name in gsheet exactly).
    """
    global _inchi2abbrev
    global _odor_inventory_gsheet

    if odor_name in skip or odor_name in hardcoded_not_to_abbrev:
        return odor_name

    if odor_name in hardcoded_odor2abbrevs:
        return hardcoded_odor2abbrevs[odor_name]

    if odor_is_mix(odor_name):
        return 'mix'

    if len(args) == 1:
        inchi2abbrev = args[0]

    elif len(args) == 0:
        if _inchi2abbrev is None:
            if _odor_inventory_gsheet is None:
                odor_inventory_gsheet(**kwargs)

            _inchi2abbrev = inchi2abbrev_dict(_odor_inventory_gsheet,
                allow_orphan_abbrevs=allow_orphan_abbrevs)

        inchi2abbrev = _inchi2abbrev

    verbose = False
    if 'verbose' in kwargs:
        verbose = kwargs['verbose']

    inchi = convert(odor_name, from_type='name', verbose=verbose)
    if pd.isnull(inchi):
        print(f'could not find inchi for odor {odor_name}!')
        return inchi

    if pd.isnull(inchi):
        return inchi

    if inchi not in inchi2abbrev:
        return None

    return inchi2abbrev[inchi]


def odor2abbrev_dict(odors, single_letter_abbrevs=True):
    """Takes a unique list of odors and returns abbreviations.

    If single_letter_abbrevs is True, abbreviations are assigned from the
    alphabet. Otherwise, abbreviations are looked up from abbreviation
    column in odor inventory Google sheet.
    """
    if single_letter_abbrevs:
        # Just to have consistent case w/ capital alphabetical labels.
        mix_abbrev = 'MIX'
    else:
        mix_abbrev = 'mix'

    o2a = dict()
    found_mix = False
    i = 0
    for o in odors:
        if odor_is_mix(o):
            if found_mix:
                raise ValueError('only expected one mix per expt')
            else:
                found_mix = True
                o2a[o] = mix_abbrev
        else:
            if single_letter_abbrevs:
                o2a[o] = chr(ord('A') + i)
                i += 1
            else:
                abbrev = odor2abbrev(o)
                if abbrev is None:
                    warn(f'no abbreviation in odor inventory gsheet for odor "{o}"!')
                    # consistent case be damned here.
                    abbrev = chr(ord('A') + i)
                    i += 1

                o2a[o] = abbrev
    return o2a

