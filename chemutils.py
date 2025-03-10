
# TODO delete and use python3 version if there is one (or delete below comment if not).
# with current escape sequence in `num_regex` def, this also required for py3
#
# So that raw string with unicode degree symbol for parsing temperatures
# could also work in Python 2.
from __future__ import unicode_literals

import atexit
from collections import Counter
import functools
import inspect
import logging
from logging import warning as warn
import os
from os.path import split, join, exists
from pathlib import Path
import pickle
from pprint import pformat
import re
import time
from typing import Optional, Callable
import traceback
import urllib.error

import numpy as np
import pubchempy as pcp
from pubchempy import Compound
import pandas as pd
# TODO maybe also install and use the pandas support module for this.
# may need to do conversion of whole pandas objects in here than, rather
# than just defining conversion fns (that use pint) on single elements...
import pint
# Should be included with pint. Some pint fn can throw this.
from tokenize import TokenError
# TODO add rdkit to requirements if not already there
# (does it still need special build for inchi support? pypi version have it?)
from rdkit.Chem.inchi import InchiToInchiKey
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
# TODO add to setup.py
from chemspipy import ChemSpider
from chemspipy.search import Results
from chemspipy.errors import ChemSpiPyRateError


session = requests.Session()
# https://stackoverflow.com/questions/15431044
retry = Retry(total=5, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
session.mount('http://', HTTPAdapter(max_retries=retry))
session.mount('https://', HTTPAdapter(max_retries=retry))

# TODO prob delete this if i'm not gonna bother silencing other annoying debug warnings
# too (and just change level to warning, which is the default)
# To get rid of annoying 'DEBUG:selector_events.py:59 Using selector: EpollSelector'
logging.getLogger('asyncio').setLevel(logging.WARNING)

# (assuming __file__ available)
# mainly makes sense being used in an editable install, but w/e
# TODO maybe move to ./log/, so as to stop screwing with tab complete when trying to
# edit this .py file...
_this_file = Path(__file__).resolve()
log_path = _this_file.with_suffix('.log')

# TODO TODO TODO finish configuring to also *APPEND* to a central file by default
# https://docs.python.org/3/howto/logging-cookbook.html#logging-to-multiple-destinations

# NOTE: this needs to come before addition of handler below, or else this config will
# not apply (still?)!
#
# TODO what is '-12s's? name?
#shared_log_format = '%(name)-12s: %(levelname)s [%(filename)s:%(lineno)d] %(message)s'
shared_log_format = '%(levelname)s[%(filename)s:%(lineno)d]: %(message)s'
logging.basicConfig(
    # TODO happy with this?
    format=f'%(asctime)s {shared_log_format}\n', filename=log_path, filemode='a'
)

# https://stackoverflow.com/questions/52161735/logging-warn-add-stacktrace
class WarnWithStackHandler(logging.StreamHandler):
    def emit(self, record):
        if record.levelno == logging.WARNING:
            print()
            stack = traceback.extract_stack()
            # TODO maybe show a bit less? any significance to `:-7`?
            # skip logging internal stacks
            stack = stack[:-7]
            for line in traceback.format_list(stack):
                print(line, end='')

        super().emit(record)

sh = WarnWithStackHandler()

formatter = logging.Formatter(shared_log_format)
sh.setFormatter(formatter)

# TODO still works w/ __name__, right (NO!)?
# TODO TODO why does cookbook add to root handler ('') instead?
# https://stackoverflow.com/questions/50714316
#logging.getLogger(__name__).addHandler(sh)
logging.getLogger('').addHandler(sh)

# TODO TODO delete?
# doesn't seem to be doing what i wanted.... (is this other direction?)
logging.captureWarnings(True)
# TODO color console log output (according to level)
# https://stackoverflow.com/questions/384076

pkg_data_dir = split(__file__)[0]

# Can't use default xlrd now that it has dropped support for XLSX files, and although
# pandas >= 1.2 might default to openpyxl if it's installed, it seems that some earlier
# pandas versions (e.g. 1.1.5) may not. Note that writing still uses xlsxwriter, which
# itself should be using openpyxl.
pandas_excel_engine = 'openpyxl'

# See notes in pint docs about concerns when pickling pint quantities or
# trying to use multiple unit registries...
# (may want to just convert everything to some standard units, and return
# only the scalar magnitudes from chemutils)

# TODO if I call set_application_registry here and in some other script using
# units originally created / loading from pickles here, does that have the
# effect of enforcing there is only one true shared registry?
# (for now, will just get the ureg defined in here via `chemutils.ureg`)
ureg = pint.UnitRegistry()
# TODO if this is required for correct pickling / unpickling behavior,
# does pickling / unpickling fail if this is not set? PR to make it fail?
# not exactly sure what this does...
pint.set_application_registry(ureg)

cache_file = os.path.expanduser(join('~', '.chemutils_cache.p'))

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
    # TODO maybe handle this some other way / provide some way of indicating
    # whether we want units + sources for each property, if i decide i don't
    # want to also store those things for this (since it is always from pubchem
    # w/ same units)
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
    # TODO i guess this has the same dimensionality as 'density'??
    # maybe don't just lookup via dimensionality then, if i want
    # different answers for density and concentration (g/L) inputs?
    # TODO or some way to have two different substances in pint unit repr?
    # (they'd be same in density case vs diff in conc case)
    '[length]^-3 [mass]': preferred_concentration_unit,
    '[length]^-3 [substance]': 'molar',
    '[length]^-1 [mass] [time]^-2': 'kPa',
}
# Checking that pint parses all preferred units to the expected
# dimensionality.
for dim, preferred_unit_str in dimensionality2preferred_units.items():
    assert ureg(preferred_unit_str).check(dim)

# TODO maybe this should completely replace the above?
property2preferred_units = {
    'concentration': preferred_concentration_unit,
    'density': preferred_density_unit,
    # Just what NIST uses. May want to change.
    'k_henry': 'mol / kg / bar',
    'vapor_pressure': 'kPa',
    'water_solubility': 'g / L',
    # TODO want this here?
    'molar_mass': 'g / mol',
}


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

    # TODO maybe also use new_units here to avoid having to parse a second
    # time??
    m = (series['value'] * ureg(unit_str)).to(new_unit_str).magnitude

    new_series = series.copy()
    new_series['value'] = m
    new_series['units'] = new_unit_str
    return new_series


# TODO maybe the pint pandas module can do some / all of this?
# TODO refactor add_properties / add_properties_to_cache to use this or
# something like this?
def change_frame_units(df, prop=None, prop_col='property', value_col='value',
    units_col='units'):
    """Returns single-property DataFrame with property in preferred units.

    Assuming that input is not null in any of the columns this fn uses, all
    unit strings are parseasble, and all property values are numeric.
    """
    # TODO delete
    '''
    print(type(df))
    print(df)
    return df
    '''
    #
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


# TODO maybe move the whole hardcode thing to some external data file like this
# TODO facilities for deployment specific additions to the hardcoded stuff
# (out of package dir)
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

    # TODO add a 'weight' column that can be filled in for values with multiple
    # sources
    assert not df.duplicated(subset=['property','inchi','sources']).any()

    # TODO maybe check inchi against name (but otherwise don't use name)
    # (to catch some data entry errors, where inchi and name are misaligned)
    name_cols = [c for c in df.columns if c.startswith('name')]
    if len(name_cols) > 0:
        assert len(name_cols) == 1
        name_col = name_cols[0]
        # Since we are ffilling 'inchi', to avoid null values causing false
        # positive duplicates.
        df[name_col] = df[name_col].ffill()
        assert not df.duplicated(subset=['property', name_col, 'sources']).any()

    # NOTE: I only added group_keys=False because when running w/ pandas==1.5.3, I got a
    # FutureWarning saying that would preserve old behavior. not sure it matters for me
    # here...
    df = df.groupby('property', group_keys=False).apply(change_frame_units)
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

        # TODO again, factor out this prop series creation stuff?
        to_update_f[to_type][inchi] = pd.Series({
            # TODO maybe refactor everything so it doesn't duplicate the
            # property name in each key??? not sure why i ever thought that'd be
            # a good idea... was it just for the few cases where i use a wide
            # spreadsheet, using the prefixes to keep the column names unique?
            to_type: row.value,
            f'{to_type}_units': row.units,
            f'{to_type}_sources': row.sources
        })


# TODO TODO TODO flag to ignore hardcoded / manual overrides, to the extent that
# they are causing problems. maybe phase them out altogether?

# TODO TODO try deleting / ignoring this hardcoded stuff and see if it still
# works. my changes to normalize_name probably fixed a few of these cases.
# (and if it can be ignored, delete it)
# TODO TODO TODO modify so these all still apply when ignore_cache=True
# (have that just apply to values that originally came from a lookup)
# (maybe have another flag / another level of ignore_cache that also applies to
# these values)
hardcoded = {
    'name': {
        'cas': {
            # TODO deal w/ linalool (or maybe at cid level?)
            # (it was some spelling that's not the one below)
            'g-hexalactone': '57129-70-1',
            'g-octalactone': '104-50-7',
            'g-decalactone': '1336-42-1',
            'a -pinene': '102640-64-2',
            'b -pinene': '127-91-3',
            # might not be right
            '(1S)-(+)-3-carene': '1352058-87-7',
            'a -terpineol': '10482-56-1',
            # trans? racemic? (what to do w/ cas number in that case?)
            'linalool oxide': '10448-29-0',
            'E2-hexenal': '1335-39-3',
            '4-ethyl guaiacol': '2785-89-9',
            'E2-hexenol': '114411-82-4',
            'Z2-hexenol': '928-94-9',
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


# TODO TODO conversion to names with a priority on various sources:
# [0) hardcoded stuff]
# 1) names in inventory
# 2) hallem
# 3) other (good default available in pubchempy? probably not iupac...)


# TODO TODO use this in all the other places that do something similar
# TODO TODO unit test this!
def is_series(obj):
    # TODO are my fears that the full path the the series object might change /
    # have changed at some point justified, or should i just check the fully
    # qualified path (or comparison of type to that of initialized empty one?)
    # (if not, maybe just use isinstance(obj, pd.Series)? or issubclass?)

    # Adapted from: https://stackoverflow.com/questions/2020014
    _type = type(obj)
    # The intermediate part (of __module__) that is not matched is:
    # 'core.series' for me, w/ pandas 0.25.1
    if _type.__name__ == 'Series' and _type.__module__.startswith('pandas.'):
        return True

    # TODO subclasses? (not that it really matters for any of my use cases
    # now... it just might if i wanted to solve the is_series matter once and
    # for all...)

    return False


# TODO TODO use this in all the other places that do something similar
# TODO TODO unit test this!
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

        # TODO test cases w/ other things that raise value their own valueerrors
        # (or other errors) on casting to bool? check the str message is the
        # numpy array one? if there are other ValueErrors, and this is the only
        # error casting to bool should raise, maybe this is the behavior i want
        # though?
        except ValueError as e:
            return False


def add_hardcoded_overrides(_cache, overwrite_existing=True):
    """
    Args:
    overwrite_existing (bool): If True, non-null keys will not be overwritten by
        values in `hardcoded` (though values in `hardcoded_type2null_keys` will
        still overwrite with `None` unconditionally).
    """
    for ft, null_keys in hardcoded_type2null_keys.items():
        for nk in null_keys:
            for tt in chem_id_types:
                _cache[ft][tt][nk] = None

    # TODO delete
    first = True
    #
    for ft, tts in hardcoded.items():
        for tt in tts:
            if overwrite_existing:
                _cache[ft][tt].update(tts[tt])

            # TODO TODO add a test that ensures this is equivalent to the update
            # call in the case where there are no overlap between existing
            # values and those to be hardcoded
            else:
                f2t_dict = _cache[ft][tt]
                hardcoded_f2t_dict = tts[tt]
                for hk, hv in hardcoded_f2t_dict.items():
                    if hk not in f2t_dict or isnull(f2t_dict[hk]):
                        # TODO delete
                        '''
                        if first:
                            first = False
                            if hk not in f2t_dict:
                                print('Not already in cache')
                            else:
                                print('Existing was null')
                            print("_cache['" + "']['".join([ft, tt, hk]) + "']")
                            print('Hardcoding:', hk, hv)
                            import ipdb; ipdb.set_trace()
                        '''
                        #
                        f2t_dict[hk] = hv

    if not overwrite_existing and not first:
        import ipdb; ipdb.set_trace()


def init_cache():
    # TODO if load_state only adds properties (or whatever i rename that to)
    # behind into, should probably be consistent here
    # (same w/ other things not added behind all from types)
    # TODO maybe modify this fn so it can init only sub-caches not in current
    # cache (and then always call init_cache some way like that on loading)?
    _cache = {ft: {tt: dict() for tt in chem_id_types + properties}
        for ft in chem_id_types
    }
    add_hardcoded_overrides(_cache)
    return _cache


def _clear_one_cache(ft, tt):
    assert ft in cache and tt in cache[ft]
    cache[ft][tt] = dict()
    # TODO need to do something w/ is_manual here too?
    save_cache(merge_with_saved=False)


# TODO maybe relax str restriction on values so this works in compound ID case
# too?
# allow null as val? probably not? null allowed as keys, as-is?
def _clear_vals_in_cache(ft, tt, vals):
    assert ft in cache and tt in cache[ft]

    if type(vals) is str:
        vals = [vals]
    else:
        try:
            iter(vals)
        except TypeError:
            raise ValueError('vals must be str or iterable of strs')

    c = cache[ft][tt]
    for v in vals:
        assert type(v) is str
        try:
            del c[v]
        # Does nothing and moves on to next value if any given val not present.
        except KeyError:
            pass

    # TODO need to do something w/ is_manual here too?
    save_cache(merge_with_saved=False)


# TODO TODO if any of these clearing fns delete stuff that has is_manual set,
# also delete that part of is_manual
def clear_cache(*args):
    global cache
    if len(args) == 0:
        cache = init_cache()

    elif len(args) == 2:
        # args: (from_type, to_type)
        _clear_one_cache(*args)

    elif len(args) == 3:
        _clear_vals_in_cache(*args)

    else:
        raise ValueError('can only call with 0, 2, or 3 arguments')


# TODO options to clear null in specific cache?
# and do i really want this write arg? (why didn't i have it above, in old
# 0-arg version of clear_cache?)
def _clear_null_cache_vals(write=True):
    global cache
    for from_type, to_type_to_cache_dict in cache.items():
        for to_type, f2t_dict in to_type_to_cache_dict.items():
            # Need list() since .items() is an iterator that expects contents
            # not to change.
            for f, t in list(f2t_dict.items()):
                # TODO TODO use this commented bit to replace the below after
                # verifying isnull works
                #if isnull(t):
                #    del f2t_dict[f]

                if not hasattr(t, 'shape') and pd.isnull(t):
                    del f2t_dict[f]

                # TODO maybe error handling in case it had shape but not .all()
                # (or more specific way of detecting series)
                elif hasattr(t, 'shape') and t.isnull().all():
                    del f2t_dict[f]

    if write:
        # TODO and why merge_with_saved=False here? document that param...
        save_cache(merge_with_saved=False)


def delete_cache():
    """Deletes and clears the cache at cache_file
    """
    clear_cache()
    if exists(cache_file):
        os.remove(cache_file)


def load_state():
    """Loads `cache` and `is_manual` from `cache_file`.
    """
    if exists(cache_file):
        try:
            with open(cache_file, 'rb') as f:
                data = pickle.load(f)
                _cache = data['cache']
                _is_manual = data['is_manual']

            # TODO delete (was planning to use to load inverted maps here)
            '''
            for from_type, to_type_dict in _cache.items():
                for to_type, f2t_cache in to_type_dict.items():
            '''
            #

            add_hardcoded_overrides(_cache, overwrite_existing=False)

            # TODO shouldn't this logic (also, w/o refactoring) be in
            # init_cache?
            for k in properties:
                if k not in _cache['inchi']:
                    _cache['inchi'][k] = dict()

            return _cache, _is_manual

        # TODO TODO TODO what raises this ValueError????
        # (pickle load? put it only around the potentially offending section!!!)
        except ValueError as e:
            raise
            # TODO TODO probably just print warning, maybe print a command
            # that can be manually run to delete cache / prompt about whether
            # cache should be deleted.
            # TODO maybe all changes to the cache should be round tripped
            # through this load, to check this doesn't happen?
            '''
            print('Cache was in unreadable format. Deleting.')
            delete_cache()
            # TODO cache guaranteed to be defined here...? not sure it is. test.
            # Will have already been init'd in delete_cache call.
            return cache, dict()
            '''
    else:
        return init_cache(), dict()


def save_cache(merge_with_saved=True):
    """
    Load old cache first, so that clear_cache doesn't have to worry about
    screwing up on-disk cache when atexit call to save_cache triggers.
    """
    # TODO why did i ever want merge_with_saved again??
    if merge_with_saved:
        # TODO need "global cache", as long as cache def is below fn def, or
        # what?  if that's all, move "cache = load_state()" just above this, and
        # declare as global in load_state, cause circularity (?)
        # TODO also want to hae merge_with_saved apply to is_manual somewhat
        # similarly (that's the second returned value currently ignored here)?
        old_cache, _ = load_state()

        # So any current values overwrite old_cache values.
        assert old_cache.keys() == cache.keys()
        for from_type, to_type_dicts in cache.items():
            old_cache_to_type_dicts = old_cache[from_type]
            # Needed to comment to modify cache to add 'density' key behind
            # 'inchi'.
            #assert old_cache_to_type_dicts.keys() == to_type_dicts.keys()
            for to_type, f2t_dict in to_type_dicts.items():
                if to_type not in old_cache_to_type_dicts:
                    old_cache_to_type_dicts[to_type] = dict()

                # TODO how to allow deleting keys (from terminal dict)?

                old_cache_f2t_dict = old_cache_to_type_dicts[to_type]
                old_cache_f2t_dict.update(f2t_dict)

        to_save = old_cache
    else:
        to_save = cache

    # TODO TODO TODO change so it write to a temp file and then moves over old thing
    # when done (to be more atomic, so being interrupted isn't as bad)
    with open(cache_file, 'wb') as f:
        # TODO TODO TODO is is_manual ever curently defined?
        data = {'cache': to_save, 'is_manual': is_manual}
        pickle.dump(data, f)


# TODO fn to manually enter persistent overrides?

# TODO test
cache, is_manual = load_state()

# TODO maybe factor into load_state / some other check fn called there
# TODO maybe add to cache fn to check units of these if they are properties?
# TODO TODO why did i decide to do this on load? why not just do it before each addition
# to cache / on write?
# TODO TODO TODO change to checking in a wrapper of each lookup fn / similar
for prop, preferred_unit_str in property2preferred_units.items():
    prop_units = ureg(preferred_unit_str)
    prop_unit_key = f'{prop}_units'
    #print('checking cache units for:', prop)
    for ft, ft_cache in cache.items():
        if prop in ft_cache.keys():
            f2p_cache = ft_cache[prop]
            for f, p in f2p_cache.items():
                assert prop in p and prop_unit_key in p
                # TODO delete
                # a crude check that the object is a pd.Series
                #assert hasattr(p, 'shape'), 'not a series'
                #
                assert is_series(p), 'not a series'

                if isnull(p):
                    # TODO delete after it is clear it is equiv to my own
                    # `isnull`
                    # This was the old way I calculated it.
                    assert pd.isnull(p[prop])
                    #
                    continue

                # TODO should it be an error for units to not be specified?
                # (probably, that's how i handle it in some lookup fns)
                # (-> so delete any entries in cache that have null units!!!)
                c_units = p[prop_unit_key]
                assert prop_units.check(c_units), f'{prop_units} != {c_units}'

atexit.register(save_cache)


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


_nonpersistent_cache = dict()
# TODO use this for basically every lookup fn?
# or still let convert handle caching for some of the functions?
# way to have it refer to the unwrapped versions in those cases?
# maybe just don't use decorator syntax at that point, and define new
# cached / un-cached versions of same fns, using same wrapper as decorator?
# TODO ideally refactor so that lookup (e.g. fetching html) and parsing (of the lookup
# output) can be separated (e.g. so we can try diff parsing strategies on cached html)
def cached(fn):
    # TODO TODO document what this adds (does it add ignore_cache option?)

    # TODO TODO rename all f->from_type, t->to_type / similar

    # Would need to modify if I wanted kwargs for the decorator that explicitly
    # specify the from / to types.
    # https://stackoverflow.com/questions/627501
    f, t = fn.__name__.split('2')

    #assert f in cache and t in cache[f], f'{f}->{t} not in cache'
    # TODO maybe just handle correctly in init_cache / load_state?
    # (would probably need to either hardcode new things to add alongside
    # any new conversion fn defs, or enumerate types from possible conversion
    # fn defs)
    assert f in cache
    if t not in cache[f]:
        warn(f'trying to cache {fn.__name__}, but {t} not in {f} cache. adding empty '
            f'dict for {t}.'
        )
        cache[f][t] = dict()

    # TODO maybe modify so manual hardcoding into `properties`
    # not required for fns referencing new "to" types?
    # (above does this, but maybe prefer to enumerate fns rather than waiting
    # for an add)

    # TODO maybe specifically only add kwargs that wrapped fn takes somehow,
    # rather than taking any kwargs, and assuming user is calling wrapped
    # fn w/ only kwargs it can take?
    # (or how to lookup kwargs of wrapped fn post hoc?)
    # So this cached decorator doesnt change function __name__
    # TODO TODO TODO modify ignore_cache so it still caches calls w/ same input *within
    # a run* (so we don't need to make a set of input IDs before applying lookup fn to a
    # pandas series, for instance)
    @functools.wraps(fn)
    def cached_fn(key, ignore_cache=False, ignore_cache_null=False,
        leave_nonnull_cache=False, **kwargs):
        """
        leave_nonnull_cache: if True, null values will not overwrite
            cached non-null values.
        """
        np_key = (f, t, key)
        # ignore_cache shouldn't apply to this one, as it is within each run
        if np_key in _nonpersistent_cache:
            return _nonpersistent_cache[np_key]

        # TODO TODO fix surrounding expectations so this doesn't break stuff?
        # don't want to have to manually create a series w/ all null keys in all those
        # separate places ideally...
        # (see inchi2vapor_pressure usage in vcf_scraper/make_csv.py for example where
        # this broke something. 'vapor_pressure' column was no longer defined)
        # TODO is this decorator also broken if fn happens to return None?
        '''
        if pd.isnull(key):
            _nonpersistent_cache[np_key] = key
            return key
        '''

        ft_cache = cache[f][t]

        verbose = False
        if 'verbose' in kwargs:
            if kwargs['verbose']:
                verbose = True

            if not allowed_kwarg(fn, 'verbose'):
                # TODO delete. trying to test my changes to inspect module usage.
                # TODO test kwonly args too
                try:
                    value = fn(key, **kwargs)
                    assert False, f'expected kwarg=verbose to cause error...'
                # TODO what is right kind of error again?
                except:
                    pass
                #
                kwargs.pop('verbose')

        if not ignore_cache and key in ft_cache:
            value = ft_cache[key]
            # TODO maybe also support same null checking that works
            # w/ all null series here? (currently just using for single elements
            # where chemspider id is null though...) (comment is probably no
            # longer relevant now that i'm using my own `isnull`, but check more
            # closely)
            if not ignore_cache_null or isnull(value):
                if verbose:
                    print(f"{fn.__name__}('{key}') returning {value} "
                        "from cache"
                    )
                return value

        value = fn(key, **kwargs)

        _nonpersistent_cache[np_key] = value

        # TODO clean this conditional up
        # TODO + test it
        if not (leave_nonnull_cache and isnull(value) and
            key in ft_cache and not isnull(ft_cache[key])):

            ft_cache[key] = value

        return value

    return cached_fn


# From Wikipedia page on InChI
_t_ster = 'tetrahedral stereochemistry of atoms and allenes'
inchi_prefix2layer = {
    'c': 'atom connections',
    'h': 'hydrogen atoms',
    'p': 'proton sublayer',
    'q': 'charge sublayer',
    'b': 'double bonds and cumulenes',
    # TODO do t and m differ in any particular way?
    't': _t_ster,
    'm': _t_ster,
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
    parts = inchi.split('/')
    if no_stereochem:
        keep = {'c','h'}
    else:
        # NOTE: 'b' layer is only one of the stereochemistry layers.
        # It is specifically "double bonds and cumulenes", whereas "t"/"m"/"s" layers
        # can have other types of stereochemistry information.
        # TODO might want to handle those other stereochem layers consistently8 w/ 'b'
        # layer... (any complications though?)
        keep = {'c','h','b'}
    return '/'.join(parts[:2] + [p for p in parts[2:] if p[0] in keep])


def inchi_layers(inchi):
    # First part should just be 1 or 1S (for "standard InChI")
    # Second part is the chemical formula w/o a prefix.
    prefixes = [x[0] for x in inchi.split('/')[2:]]
    layers = [inchi_prefix2layer[p] for p in prefixes]
    return layers


def inchi_layer_set(inchis):
    from functools import reduce
    # TODO why again can list 'sum' be taken directly by agg but
    # set.union can't? way to reference list sum as set.union is ref'd?
    union = lambda x: reduce(set.union, x)
    return inchis.apply(lambda x: set(inchi_layers(x))).agg(union)


# TODO move to / dupe to hong2p? (along w/ unit test i now have for this)
def is_one2one(df: pd.DataFrame, x: str, y: str, *, dropna: bool = True) -> bool:
    """Returns True if x and y are 1:1 with each other in df.

    Handling of null still broken.
    """
    # TODO maybe modify this to print violations of 1:1-ness, when there are any
    # (otherwise write another fn for this)?

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
    with pd.option_context('display.max_colwidth', -1):
        print(df.to_string(index=False))


# TODO TODO modify this to actually use basic_inchi w/ no_stereochem=True
# (at least provide option to make that comparision here)
# TODO better name?
def inchi_diff_in_details(df):
    assert 'inchi' in df.columns
    other_keys = [x for x in df.columns if x in chem_id_types]
    other_keys += [x for x in df.columns if x in equivalent_col_names]
    cols_to_show = other_keys + ['inchi']
    # TODO also count name/cas/(name,cas) here, as w/ inchi_counts?
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

    # 208
    #print(len(df.drop_duplicates(subset=['name','cas_number','inchi'])))
    #
    df = df.drop_duplicates(subset='inchi').copy()
    # 201
    #print(len(df))

    df['inchi_counts'] = df.inchi.apply(lambda i: inchi_counts.at[i])

    df['basic_inchi'] = df.inchi.apply(basic_inchi)

    cols_to_show += ['inchi_counts', 'basic_inchi']

    # TODO maybe just print part after common prefix for each of these?
    # (part after basic inchi / basic inchi w/o h)

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


# TODO why does "tetrahedral stereochemistry of atoms and allenes" layer
# show up more than "double bonds..."?
# TODO maybe draw structures with / without this info (if possible) to see what
# information it contains?
def count_inchi_layers(inchis):
    counts = Counter(inchis.apply(inchi_layers).agg('sum'))
    return counts


# TODO rename to is_chem_id[_var] or something?
def convertable(chem_id_type: str) -> bool:
    """Returns whether the string input is recognized as a convertable type.

    (valid for the `to_type` and `from_type` arguments to the `convert` function)
    """
    if chem_id_type in chem_id_types or chem_id_type in equivalent_col_names:
        return True
    return False


# TODO TODO tqdm in all but single element case (way to make it work in case
# where something is iterated over w/ call to convert at each iteration?)

cid2compound_cache = dict()
# TODO in ignore_cache case, still make a interpreter run / call specific
# cache, or like de-dupe and re-dupe, so as to still test new behavior, but also
# not waste time. (not as much of a priority if clear_cache approach works)
# TODO add flag to also return conversion table (so that input not suitable for
# inserting new values along original can still have access)
# (like w/ hallem data, where chem id is the sole row index)
# TODO TODO add flag to keep original identifiers (when possible) (and err
# if this flag is passed and it's not possible to insert them in the input
# object) ('orig_*')
# TODO add kwarg flag to enable checks that output id is 1:1 w/ input ids
# (unless theres already some reason to expect this to be true?)
# TODO TODO try to refactor use of cache across here+cached decorator to share the
# same logic (+ to share logic for having a non-persistent (within run) cache, that is
# used no matter what
def convert(chem_id, from_type=None, to_type='inchi', *, verbose: bool = False,
    allow_nan: bool = False, allow_conflicts: bool = True, ignore_cache: bool = False,
    exclude_cols=tuple(), already_normalized: bool = False, dropna: bool = True,
    report_missing: bool = True, report_conflicts: bool = True,
    try_non_normalized: bool = True, check_one2one: bool = False,
    keep_originals: bool = False, orig_prefix: str = 'orig_'):
    """
    Args:
        dropna: whether to drop input subset where lookups fail (if setting this to
            False, probably want `allow_nan=True`)

        allow_nan: whether to err if any lookups fail (does not depend on
            `dropna=False`)
    """
    def conversion_fail_errmsg():
        return f'conversion from {from_type} to {to_type} failed for'

    def conversion_fail_err():
        msg = conversion_fail_errmsg() + f' {chem_id}'

        if not allow_nan:
            raise ValueError(msg)

        elif report_missing:
            print(msg)


    if check_one2one:
        raise NotImplementedError

    valid_ignore_cache_vals = (True, False, 'if_null')
    if ignore_cache not in valid_ignore_cache_vals:
        raise ValueError(
            f'ignore_cache must be one of: {valid_ignore_cache_vals}'
        )

    could_keep_originals = False
    # TODO test w/ index/series/df (w/ index / columns of matching name)
    if hasattr(chem_id, 'shape') and len(chem_id.shape) > 0:
        if to_type in equivalent_col_names:
            to_name = to_type
            to_type = equivalent_col_names[to_type]
        else:
            to_name = to_type

        # TODO TODO TODO whenever this thing maps >1 of old id to 1 of new type,
        # report! (series/index case may be sufficient)
        # flag for it indep of verbose. maybe flag to err?

        # This covers both Index and Series.
        if len(chem_id.shape) == 1:
            if from_type is None:
                if chem_id.name is None:
                    raise ValueError('either pass from_type or name axes')

                # TODO can series index also have a name? if i'm trying to
                # support the index renaming case, would also have to check that
                # here (might want to enforce only either the index or series
                # itself has a name indicating it can be converted)
                from_type = chem_id.name
                if from_type in equivalent_col_names:
                    from_type = equivalent_col_names[from_type]

            # TODO why are allow_nan and report_missing not just threaded thru?
            # + is this all the current kwargs?
            #
            # Forcing allow_nan to True so we can report each of the failing
            # lookups.
            fn = lambda x: convert(x, from_type=from_type, to_type=to_type,
                verbose=verbose, ignore_cache=ignore_cache, allow_nan=True,
                already_normalized=already_normalized, report_missing=False
            )

            # TODO test each of these cases
            # TODO also test w/ passing / failing allow_nan case for each
            # This means it was a Series, not an Index.
            # TODO TODO which case do i actually want to support:
            # series w/ values to convert in index, or w/ them as the value of
            # the series???? either?
            if hasattr(chem_id, 'index'):
                '''
                raise NotImplementedError
                converted = chem_id.rename(fn).rename(to_name)
                to_check = chem_id.index
                '''
                # TODO maybe get rid of this whole if/else, if i don't also plan
                # to support the converting-index-of-series case
                converted = chem_id.map(fn)
                converted.name = to_type
                to_check = converted
            else:
                converted = chem_id.map(fn)
                converted.name = to_type
                to_check = converted

            missing = to_check.isnull()
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

            # TODO also support checking names of single rows/columns and
            # converting those (adding a column)?
            if from_type is not None:
                raise NotImplementedError('only conversion of named indices'
                    ' or columns supported for DataFrames'
                )

            df = chem_id.copy()

            # TODO this all the current kwargs?
            #
            # allow_conflicts and exclude_cols are only relevant for DataFrame
            # case, and therefore do not need to be passed.
            # from_type should always be inferred in DataFrame case as well.
            kwargs = {
                'to_type': to_type,
                'allow_nan': allow_nan,
                'ignore_cache': ignore_cache,
                'dropna': dropna,
                'already_normalized': already_normalized,
                'report_missing': report_missing,
                'verbose': verbose
            }
            converted_an_index = False
            if convertable(df.index.name):
                if keep_originals:
                    # TODO probably factor our this orig col handling logic to
                    # a nested fn
                    orig_col = orig_prefix + df.index.name

                    # TODO TODO replace w/ asserting it's either not there OR equal to
                    # the column already there
                    assert orig_col not in df.columns

                    orig_cols_added.append(orig_col)
                    df[orig_col] = df.index.copy()

                df.index = convert(df.index, **kwargs)
                converted_an_index = True

            if convertable(df.columns.name):
                if keep_originals:
                    orig_col = orig_prefix + df.columns.name

                    # TODO TODO replace w/ asserting it's either not there OR equal to
                    # the column already there
                    assert orig_col not in df.columns

                    orig_cols_added.append(orig_col)
                    df[orig_col] = df.columns.copy()

                df.columns = convert(df.columns, **kwargs)
                converted_an_index = True

            if not converted_an_index:
                kwargs['allow_nan'] = True
                kwargs['dropna'] = False
                kwargs['already_normalized'] = True
                kwargs['report_missing'] = False

                attempts = dict()
                for c in df.columns:
                    if c not in exclude_cols and convertable(c):
                        # TODO TODO factor this into a fn / have convert
                        # recursively handle this w/ code in same place somehow
                        if c in equivalent_col_names:
                            ft = equivalent_col_names[c]
                        else:
                            ft = c

                        if keep_originals:
                            orig_col = orig_prefix + c

                            # TODO TODO replace w/ asserting it's either not there OR
                            # equal to the column already there
                            assert orig_col not in df.columns

                            orig_cols_added.append(orig_col)
                            df[orig_col] = df[c].copy()

                        # TODO TODO test other normalize case + other fn name
                        # lookup cases in c in equivalent_col_names case, since
                        # that caused the failure here
                        norm_fn_name = 'normalize_' + ft
                        if norm_fn_name in globals():
                            if verbose:
                                print(f'applying {norm_fn_name} to input column')

                            # TODO still print same kind of verbose stuff as
                            # other normalization section does (refactor?)
                            norm_fn = globals()[norm_fn_name]
                            df[c] = df[c].apply(norm_fn)
                        #
                        attempts[c] = convert(df[c], **kwargs)

                attempts = pd.DataFrame(attempts)

                missing = attempts.isnull().all(axis=1)
                if missing.any():
                    cols = orig_cols_added + list(attempts.columns)
                    # TODO why not using conversion_fail_err[msg] (presumably b/c that
                    # seems set up for single chem_ids, not Series of them? update?)?
                    # delete one or the other? refactor?
                    err_str = (f'conversion from {[c for c in cols]} to {to_type} '
                        'failed for:\n'
                    )
                    missing = df[missing].drop_duplicates(subset=cols).dropna(
                        subset=cols, how='all')[cols]

                    with pd.option_context('display.max_colwidth', -1):
                        err_str += missing.to_string()

                    if not allow_nan:
                        raise ValueError(err_str)
                    elif report_missing and len(missing) > 0:
                        print(err_str)
                del missing

                input_priorities = {k: p for p, k in enumerate(chem_id_types)}

                def priority(id_type):
                    if id_type in equivalent_col_names:
                        id_type = equivalent_col_names[id_type]

                    return input_priorities[id_type]

                cols_in_order = sorted(attempts.columns, key=priority)
                values = attempts[cols_in_order].apply(lambda x: x.dropna().unique(),
                    axis=1
                )
                # TODO use that str len fn?
                conflicts = values.apply(lambda x: len(x) > 1)

                # TODO maybe in the case that to_type=inchi, say which layers
                # differ in conflict?

                # TODO TODO TODO maybe provide some indicator column that says
                # which ID was ultimately used (either highest priority of those
                # that agreed, or the whole set that agreed), for
                # troubleshooting other normalization problems

                # TODO TODO may want to find set of CIDs that are common to
                # other lookups (though maybe require they are in the highest
                # priority), and then try to resolve from there

                # TODO TODO test that this is actually respecting priority
                df[to_type] = values.apply(
                    lambda x: None if len(x) == 0 else x[0]
                )
                if dropna:
                    df.dropna(subset=[to_type], inplace=True)

                if conflicts.any():
                    if not allow_conflicts or report_conflicts:
                        # TODO make sure only unique conflicts are printed
                        print('conflicting lookup results:')
                        for i, conflict in attempts[conflicts].iterrows():
                            conflict = conflict.dropna()
                            inputs = df.loc[i, conflict.keys()]
                            cdf = pd.DataFrame({'input': inputs,
                                to_type: conflict
                            })
                            with pd.option_context('display.max_colwidth', -1):
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

    # TODO duplicated above (in at least one case?)? clean up / refactor, if so
    if from_type in equivalent_col_names:
        from_type = equivalent_col_names[from_type]

    if verbose:
        print(f'trying to convert {chem_id} from {from_type} to {to_type}')

    # To short-circuit normalization in the more-common case where null exists
    # before normalization.
    if pd.isnull(chem_id):
        return chem_id

    # TODO unit test each of these cases, including None handling

    # TODO caching that can also short circuit this step? use same cache?
    # (seems to be taking up most of time in odor2abbrev case)
    old_chem_id = chem_id
    if not already_normalized:
        norm_fn_name = 'normalize_' + from_type
        if norm_fn_name in globals():
            chem_id = globals()[norm_fn_name](chem_id)
            if verbose and old_chem_id != chem_id:
                print(f'{norm_fn_name}({old_chem_id}) -> {chem_id}')

    # Since sometimes (just normalize_cas, for now) normalize fns can return
    # null.
    if pd.isnull(chem_id):
        return chem_id

    # TODO TODO for relationships that are *guaranteed* to be one-to-one,
    # could also populate the reverse direction in the cache when the other
    # direction is filled (CID <-> inchi? definitely inchikey <-> inchi, right?)
    # ...and do they have to be 1:1?

    #if not ignore_cache:
    #if ignore_cache == False or ignore_cache == 'if_null':
    if ignore_cache != True:
        # TODO should this fail into elif if cached value is None?
        # (be consistent w/ all branches on try_non_normalized)
        if chem_id in cache[from_type][to_type]:
            val = cache[from_type][to_type][chem_id]
            # TODO need other handling of ignore_cache == False case?
            if not (ignore_cache == 'if_null' and pd.isnull(val)):
                if verbose:
                    print(f'returning {val} from cache')
                return val

        elif try_non_normalized and old_chem_id in cache[from_type][to_type]:
            val = cache[from_type][to_type][old_chem_id]
            # TODO need other handling of ignore_cache == False case?
            if not (ignore_cache == 'if_null' and pd.isnull(val)):
                if verbose:
                    # TODO here and in other similar places, replace
                    # 'chem_id' w/ value of `from_type`
                    print('falling back to non-normalized chem_id')
                    print(f'returning {val} from cache')
                return val

    # Added 2020-06-25 as part of trying to add ignore_cache == 'if_null'
    # support.
    if from_type.lower() != 'cid':
        cid = None
    else:
        # TODO have next conditionals skipped in this branch anyway?
        cid = chem_id

    #if not ignore_cache and chem_id in cache[from_type]['cid']:
    if ignore_cache != True and chem_id in cache[from_type]['cid']:
        cid = cache[from_type]['cid'][chem_id]
        if verbose:
            print(f'{chem_id} of type {from_type} had CID {cid} in cache')

        if cid is None:
            if ignore_cache != 'if_null':
                if verbose:
                    print('CID in cache was None!')
                conversion_fail_err()
                return None

    #elif (not ignore_cache and try_non_normalized and
    elif (ignore_cache != True and try_non_normalized and
        old_chem_id in cache[from_type]['cid']):

        cid = cache[from_type]['cid'][old_chem_id]
        if verbose:
            print('falling back to non-normalized chem_id')
            print(f'{old_chem_id} of type {from_type} had CID {cid} in cache')

        if cid is None:
            # TODO could probably handle this once below in `if cid is None`
            # case, rather than here and in `if` above...
            if ignore_cache != 'if_null':
                if verbose:
                    print('CID in cache was None!')
                conversion_fail_err()
                return None

    # Changed 2020-06-25 as part of trying to add ignore_cache == 'if_null'
    # support.
    #else:
    if cid is None:
        f2cid_fn_name = from_type + '2cid'

        if f2cid_fn_name not in globals():
            # TODO delete. why is this failing in cid->name case
            import ipdb; ipdb.set_trace()
            #
            raise NotImplementedError(f'define function {f2cid_fn_name} to support '
                f'conversion from type {from_type}'
            )

        if verbose:
            print('calling CID lookup function', f2cid_fn_name)

        cid = globals()[f2cid_fn_name](chem_id)

        if cid is None:
            if verbose:
                print(f'looking up CID for {chem_id} of type {from_type} failed!')

            if try_non_normalized:
                if verbose:
                    print('falling back to non-normalized chem_id')
                cid = globals()[f2cid_fn_name](old_chem_id)

        if cid is None:
            if try_non_normalized and verbose:
                print(f'looking up CID for {old_chem_id} of type {from_type} failed!')

            for tt in chem_id_types:
                to_type_cache = cache[from_type][tt]
                # To not overwrite hardcoded values for other types.
                if chem_id not in to_type_cache:
                    to_type_cache[chem_id] = None

            conversion_fail_err()
            return None

        if verbose:
            print(f'CID={cid}')

        cache[from_type]['cid'][chem_id] = cid

    # TODO TODO should i also support some way of directly going from from_type
    # to to_type, w/o requiring a set of fns about a Compound intermediate?
    # somethings i want that can't be pushed through that first?

    # TODO way to go direct to Compound?
    # (__init__ takes something called a "record dict from the PubChem PUG REST
    # service", which we might already have in the results...) idk...
    # TODO TODO TODO possible to cache this? Compounds pickleable?

    if 'compound' not in cache['cid']:
        cache['cid']['compound'] = dict()

    elif cid in cache['cid']['compound']:
        # TODO delete
        if verbose:
            print('using cid->compound cache')
        #
        compound = cache['cid']['compound'][cid]

    else:
        # TODO delete
        if verbose:
            print('NOT using cid->compound cache')
        #

        # TODO delete
        t0 = time.time()
        #

        compound = pcp.Compound.from_cid(cid)

        # TODO delete
        print(f'pcp.Compound.from_cid({cid}) took {(time.time() - t0):.3f}s')
        #import ipdb; ipdb.set_trace()
        #

        # TODO delete???
        # TODO this being saved by atexit?
        cache['cid']['compound'][cid] = compound
        #


    if compound is None:
        # TODO should i just assert false here or something?
        if verbose:
            print(f'creating Compound from CID {cid} failed!')

        for tt in chem_id_types:
            # TODO if we didn't need this though, we could factor this into a fn
            # and use both here and above in cid == None case...
            # Just because we technically did get a non-null cid before.
            if tt == 'cid':
                continue
            to_type_dict = cache[from_type][tt]
            # To not overwrite hardcoded values for other types.
            if chem_id not in to_type_dict:
                to_type_dict[chem_id] = None

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

    cache[from_type][to_type][chem_id] = to_type_val

    if verbose and to_type_val is None:
        print(f'conversion of {chem_id} from Compound to {to_type} failed!')
        conversion_fail_err()

    return to_type_val


def normalize_name(name):
    # TODO maybe remove first null check in convert if this is the expected
    # behavior of the norm fns?
    if pd.isnull(name):
        return name

    # TODO don't lowercase letters appearing by themselves?
    name = name.replace('(CAS)', '').strip().lower()

    # TODO add test cases for things this re.sub and following for loop
    # were trying to address

    # TODO test this doesn't break any of my natural_odors usage
    # TODO maybe similar length condition on stuff before open paren?
    # .{3,} means "at least 3 of any character"
    # Square bracket part means "any character but parentheses"
    # TODO TODO TODO was this causing the parts[0] IndexError on Virginie's
    # compiled data?
    # Added $ (end of string) 2020-06-17
    name = re.sub(r'\([^\(\)]{4,}\)$', '', name)

    parts = name.split()

    # TODO TODO delete try / except, after ensuring it won't happen
    try:
        normed_name = parts[0]
    except IndexError:
        warn('hacky fix to IndexError in normalize_name still there!')
        # TODO TODO TODO why were these names empty? cas? am i throwing away
        # important stuff?
        '''
        print(name)
        print(type(name))
        raise
        '''
        return None

    for a, b in zip(parts, parts[1:]):
        # We don't want to join adjacent words, or after comma.
        if a[-1] == ',' or (a[-1].isalpha() and b[0].isalpha()):
            b = ' ' + b
        normed_name += b

    corrections = {
        'linalool oxide': 'trans-linalool oxide',
        '4-ethyl guaiacol': '4-ethylguaiacol',
        # TODO check that this one is still triggered by something in
        # natural_odors (copied from normalize_name there)
        'dihydrazide ethanediimidic acid':
            'ethanediimidic acid, dihydrazide',
        # Since search of 'E3-hexenol' fails and '(E)-3-hexenol' returns
        # compound for the (Z)-3 version.
        'e3-hexenol': '(E)-hex-3-en-1-ol',
        # TODO implement more general handling of this type of sterochemistry
        # nomenclature (don't want lookup to return version w/ sterochem not
        # specified)
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
        # TODO is search case sensitive? (i.e. is initial lower() OK?)
        # (seems not to matter here at least)
        'e2-': '(E)-2-',
        'e3-': '(E)-3-',
        'z2-': '(Z)-2-',
        'z3-': '(Z)-3-',
    }
    for p, c in prefix_corrections.items():
        if normed_name.startswith(p):
            normed_name = c + normed_name[len(p):]

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

    return normed_name


def normalize_cas(cas):
    # TODO maybe remove first null check in convert if this is the expected
    # behavior of the norm fns?
    # TODO make a null_preserving decorator and use for each fn that has similar
    # logic?
    if pd.isnull(cas):
        return cas

    # TODO include comment saying which kind of inputs this is meant to process
    # (i think it was some libraries that had unnecessary quotes in their CAS values?)
    normed_cas = ''.join(cas.replace('"','').split())

    # TODO library seems to return 0-00-0 sometimes... but this is incorrect,
    # right? replace w/ NaN?
    if normed_cas == '0-00-0':
        return None

    return normed_cas


pubchem_root_url = 'https://pubchem.ncbi.nlm.nih.gov'
def pubchem_url(cid: int) -> Optional[str]:
    # TODO return None always (would prob make Optional[str] type hint more accurate,
    # but could break some stuff...)?
    if pd.isnull(cid):
        return cid

    # .0f to work w/ floats (e.g. from a Series that is a mix of ints and NaN)
    return f'{pubchem_root_url}/compound/{cid:.0f}'


# TODO any params i could add to make it only search compounds?
# maybe '&selected_id_type=cid'? unclear that helped in browser tho
# TODO or &tab=compound maybe?
def pubchem_search_url(text: str) -> str:
    # TODO need to do any encoding of text? 'β' seems to show up in url in browser
    # as-is, so maybe not (or at least not for these chars)?
    return f'{pubchem_root_url}/#query={text}'


# TODO worth caching this too (w/ decorator)?
def inchi2pubchem_url(inchi, **kwargs):
    if pd.isnull(inchi):
        return inchi
    # TODO TODO TODO fix cause of this failing in input data
    # (in natural_odors/literature_data.py inputs)
    try:
        cid = convert(inchi, from_type='inchi', to_type='cid', **kwargs)
    except pcp.BadRequestError as e:
        print(f'Error converting {inchi} to pubchem_url')
        return None

    # TODO TODO TODO make this conditional on some input kwarg / fix cases
    except ValueError as e:
        warn(str(e))
        return None
    #

    return pubchem_url(cid)


# TODO maybe also take fns <type>2results or something? which always gets
# CID as here? probably not worth it...
def name2cid(name, verbose=False):
    try:
        # TODO TODO should this have been get_compounds, as in cas2name above?
        # what's the difference?
        #results = pcp.get_synonyms(name, 'name')
        results = pcp.get_compounds(name, 'name')
    except urllib.error.URLError as e:
        traceback.print_exc()
        warn(f'{e}\nReturning None')
        return None

    if len(results) == 0:
        # TODO delete
        # TODO are there substances w/o cids? search for those too?
        '''
        substances = pcp.get_substances(name, 'name')
        substance_cids = pcp.get_cids(name, 'name', 'substance')
        #assert len(substance_cids) == 0
        import ipdb; ipdb.set_trace()
        '''
        #
        return None

    # TODO delete / pick some other strategy if this isn't always true
    # (yes, this does fail sometimes)
    ######assert len(results) == 1
    if len(results) > 1:
        # TODO TODO TODO focus on cases where there are multiple
        # fewest-InChI-layer-results, and figure out whether / what
        # other strategies are necessary to resolve results
        print(f'got multiple results from name={name}:')
        n_inchi_parts = [len(r.inchi.split('/')) for r in results]
        fewest_inchi_parts = sorted(n_inchi_parts)[0]
        counts = Counter(n_inchi_parts)
        print(f'fewest InChI parts: {fewest_inchi_parts}')
        print(f'# InChIs w/ that many parts: {counts[fewest_inchi_parts]}')
        print(counts)
        for r in results:
            print(r.inchi)
            print(pubchem_url(r.cid))
        print('')
    #

    #return results[0]['CID']
    return results[0].cid


# TODO perhaps just call name2cid here, unless we have to take different
# strategies to deal w/ len(results) > 1 case, if and when that happens
def cas2cid(cas, verbose=False):
    try:
        # TODO and maybe if i'm using get_compounds, make this a cas2compound
        # fn or something, and skip the compound.from_cid step?
        results = pcp.get_compounds(cas, 'name')
    except urllib.error.URLError as e:
        warn(f'{e}\nReturning None.')
        return None

    if len(results) == 0:
        # TODO delete
        '''
        r2 = pcp.get_synonyms(cas, 'name')
        assert len(r2) == 0
        # TODO are there substances w/o cids? search for those too?
        substance_cids = pcp.get_cids(cas, 'name' 'substance')
        assert len(substance_cids) == 0
        '''
        #
        return None

    # TODO delete / pick some other strategy if this isn't always true
    # damn... even this one fails sometimes
    #assert len(results) == 1
    if len(results) > 1:
        # TODO TODO TODO focus on cases where there are multiple
        # fewest-InChI-layer-results, and figure out whether / what
        # other strategies are necessary to resolve results
        print(f'got multiple results from CAS={cas}')
        n_inchi_parts = [len(r.inchi.split('/')) for r in results]
        fewest_inchi_parts = sorted(n_inchi_parts)[0]
        counts = Counter(n_inchi_parts)
        print(f'fewest InChI parts: {fewest_inchi_parts}')
        print(f'# InChIs w/ that many parts: {counts[fewest_inchi_parts]}')
        print(counts)
        for r in results:
            print(r.inchi)
            print(pubchem_url(r.cid))
        print('')
    #
    return results[0].cid


def inchikey2cid(inchikey, verbose=False):
    try:
        results = pcp.get_compounds(inchikey, 'inchikey')
    except urllib.error.URLError as e:
        warn(f'{e}\nReturning None.')
        return None
    assert len(results) == 1
    return results[0].cid


def inchi2cid(inchi, verbose=False):
    inchi = 'InChI=' + inchi
    try:
        results = pcp.get_compounds(inchi, 'inchi')
    except urllib.error.URLError as e:
        warn(f'{e}\nReturning None.')
        return None
    assert len(results) == 1
    return results[0].cid


def compound2cid(compound: Compound) -> int:
    return compound.cid


# TODO TODO TODO scrape name at the top of pubchem page if there is no other way
# of getting this. most of the time this would probably be preferable to IUPAC
# name (as long as it doesn't depend on how the page is accessed...).
# e.g. 'linalool' vs IUPAC '3,7-dimethylocta-1,6-dien-3-ol'
def compound2name(compound: Compound) -> str:
    # TODO use name displayed on pubchem page (accessible thru pubchempy?)?
    return compound.iupac_name


def compound2inchi(compound: Compound) -> str:
    inchi = compound.inchi
    # TODO sometimes, is it just the prefix?
    assert inchi.startswith('InChI=')
    inchi = inchi[6:]
    # TODO actually check format meets some minimum of the inchi standard
    assert len(inchi) > 0
    return inchi


def compound2inchikey(compound: Compound) -> str:
    return compound.inchikey


def compound2cas(compound: Compound, *, verbose: bool = False) -> Optional[str]:
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
        #
        # TODO sorting improve ability to repro cas in yuansheng excel file?
        # can probably resolve those things be re-looking up cas though, so shouldn't
        # need to change this
        cas = cas_numbers[0]

    return cas


def compound2smiles(compound):
    return compound.canonical_smiles


# TODO doc (/ better name). what is this for (only used in natural_odors/odors.py)?
def fmt_id_type(chem_id_type):
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


# TODO still used, or a remnant from old style cache stuff?
def name2compound(name, verbose=False):
    cid = name2cid(name, verbose=verbose)
    compound = pcp.Compound.from_cid(cid)
    return compound


# TODO support hardcoded odor->canonical names
# TODO pick chem_id automtically in order from df columns avail +
# err if necessary
# TODO maybe move back to natural_odors/odors.py now that i'm not actually using
# this elsewhere?
def normalize_chem_ids(df, chem_id, keep_originals=False, **kwargs):
    """Normalize on chem_id then back convert to names.
    """
    df = convert(df, to_type=chem_id, keep_originals=keep_originals, **kwargs)
    # TODO check that introduction of allow_nan to this conversion as well
    # (was previously only passed to line above) doesn't break any old
    # uses of this function
    df['name'] = convert(df[chem_id], to_type='name', **kwargs)
    return df


# TODO add decorator (used elsewhere), where null input gets reflected back
def chemspider_url(csid):
    if not pd.isna(csid):
        assert np.isclose(csid, int(csid))

    # ':.0f' just so if input comes from a a Series where NaN was added (thus converting
    # any prior int types to float), the URLs still work.
    return f'http://www.chemspider.com/Chemical-Structure.{csid:.0f}.html'


# NOTE: chemspider's policy is that scraping is not allowed, so may want to take further
# steps to rate limit / appear as a normal browser (mechanize w/ similar setup to
# VCF scraper?).
def chemspider_experimental_density(csid, verbose=False):
    import requests
    from bs4 import BeautifulSoup

    url = chemspider_url(csid)
    if verbose:
        print(url)
    # TODO implement some requests rate limiting
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
        # TODO maybe valid second part (more if includes whitespace)
        # is parseble as valid density unit? use `pint` for that?

        parts = row.text.split()
        try:
            curr_density = float(parts[0].strip())
        except ValueError:
            if verbose:
                print(url)
                print('Failed to parse density from row text:', row.text)
            #raise
            return None

        # TODO ever need to parse this from multiple parts?
        # TODO need to handle null case?
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
        # TODO may want to use None instead if units can not be parsed
        # (or just return everything as null in that case, as i do in at
        # least one other place...)
        'density_units': preferred_density_unit,
    })


nist_url_prefix = 'http://webbook.nist.gov'
def inchi2nist_webbook_url(inchi):
    # TODO replace w/ general url escaping library call? stdlib?
    inchi_str = 'InChI=' + inchi.replace('/','%2F').replace(',','%2C').replace(
        '+', '%2B'
    )
    # TODO need to explicitly ask for the data we want though? links should show up in
    # "Other data available" section if we have it anyway...
    # in particular, might want 'Phase Change' (for Antoinne equation parameters for
    # vapor pressure) and "Henry's Law" checkboxes.
    return f'{nist_url_prefix}/cgi/cbook.cgi?{inchi_str}&Units=SI'


def inchi2nist_henrys_law_url(inchi):
    return inchi2nist_webbook_url(inchi) + '&Mask=10'


def inchi2nist_vapor_pressure_url(inchi):
    # NOTE: often there are tables here that aren't the antoinne equation parameters we
    # want (and are there simpler records of vapor pressure, that aren't intended to
    # work across temperatures, that exist in the NIST webbook? maybe in a place i'm not
    # looking?)
    return inchi2nist_webbook_url(inchi) + '&Mask=4'


# TODO do i even use the pint= kwarg in any of the wrapped fns rn?
# TODO maybe also add kwargs to specify maximum allowable temperature
# deviations? or do none of the property fns actually return that?
# TODO maybe use pint-pandas instead?
# (only place i currently see that uses this is old chemutils/scripts/watersol.py)
'''
def opt_pint_return(property_fn):
    """
    Decorator for property lookup functions that adds a boolean `pint` kwarg to
    return a pint Unit object instead of the default pandas Series.

    If the function to be wrapped returns `None`, so does the wrapped version.
    If it returns a pandas Series with all null values, `None` is returned.

    The function must take an InChI input (without the 'InChI=' prefix), and it
    must return a pandas Series with '<property>' and '<property>_units' as
    keys (with numeric and a `str` with parseable units, respectively).

    The pandas Series may have slightly more information, such as the source of
    the lookup information.
    """
    _, prop_name = property_fn.__name__.split('2')

    @functools.wraps(property_fn)
    def wrapped_prop_fn(inchi, pint=False, **kwargs):
        series = property_fn(inchi, **kwargs)
        if not pint:
            return series

        if isnull(series):
            return None

        assert prop_name in series
        # TODO check to ensure property is numeric?
        unit_key = prop_name + '_units'
        assert unit_key in series

        unit_str = series[unit_key]
        assert type(unit_str) is str

        quantity = ureg.Quantity(series[prop_name], unit_str)
        return quantity

    return wrapped_prop_fn
'''


@cached
def cid2molar_mass(cid):
    # TODO null handling?
    compound = pcp.Compound.from_cid(cid)

    # TODO maybe delete sources + units from this fn. just to be consistent w/
    # other properties for now, so nothing breaks
    url = pubchem_url(cid)
    units = 'g/mol'

    return pd.Series({
        'molar_mass': float(compound.molecular_weight),
        'molar_mass_sources': url,
        'molar_mass_units': 'g/mol'
    })


# TODO TODO should try to pass ignore_cache of this wrapped fn (wrapper adds the
# kwarg) through to the cid2molar_mass lookup, right?
#@opt_pint_return
@cached
def inchi2molar_mass(inchi):
    cid = convert(inchi, from_type='inchi', to_type='cid')
    # TODO delete after handling null case
    assert type(cid) is int
    #
    return cid2molar_mass(cid)


def try_parse_units(string):
    """Returns either pint units or None. Raises ValueError when pint does.

    pint may also raise other errors I do not specifically anticipate in here.
    """
    curr_units = None
    try:
        curr_units = ureg.parse_units(string)

    # TODO are these the only errors parse_units could throw?
    # (from playing around with it, it seems so, but check docs / source)
    except pint.errors.UndefinedUnitError:
        pass

    # This will happen if there is an unmatched lparen.
    except TokenError:
        pass

    # This happens if p == ':' (p?)
    except AttributeError as e:
        assert str(e) == "'NoneType' object has no attribute 'evaluate'"

    # '/ext/' yields this error, w/ text 'missing unary operator "/"'
    except pint.errors.DefinitionSyntaxError:
        pass

    # 'experimentally-derived', yields this w/ text:
    # "unsupported operand type(s) for -: 'ParserHelper' and 'ParserHelper'"
    except TypeError:
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
# TODO also need to match case where degree sign is not there?
# case where F/C not there (maybe C assumed?)?
# TODO maybe check there were not multiple matches to this regex
_at_re = r'(?:@|at)\s*'
degree_re = r'\s*°\s*([CFK])'
temp_re = _at_re + float_re + degree_re
temp_range_re = _at_re + float_range_re + degree_re

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


# TODO maybe have "Insoluble in water" / "Insol in water..." return some small
# value that may be typical (for solubility case)?
# TODO see if cases that list things like "2%" have more information in other
# things returned by pubchemprops (seems they don't, at least not in what
# that library returns)
# TODO maybe if we would return 0, we should either return None or a positive
# value smaller than can be represented with the number of digits in the string
# representation it was parsed from (i.e. '0.00 M' -> 0.001)?
# TODO maybe define an acceptable temperature deviation?
# (if we only have data at like 100C, do we want to use that?)
# TODO TODO TODO should i always split on first '(' / '/' character, if there are any),
# and throw away everything after?
# TODO TODO TODO keep track of inputs strings that produce failed parsing (persistently
# too, mabye), so i can check what the problem cases are easily
# TODO see if https://chemistry-tools.readthedocs.io/ has any useful parsing stuff
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
        for part in string.split(';'):
            curr_data = parse_pubchemprops_string(part,
                expected_units=expected_units, target_temp_c=target_temp_c,
                url=url
            )
            curr_temp = curr_data['temperature_c']
            if curr_temp is not None:
                curr_diff = abs(curr_temp - target_temp_c)
                if smallest_diff is None or curr_diff < smallest_diff:
                    data = curr_data
                    smallest_diff = curr_diff

        # This coverts strings that have no parseable portions, but still
        # have a semicolon, such as:
        # "Miscible in ethanol, ethyl ether; slightly soluble in carbon
        #  tetrachloride"
        if data is None:
            return null_dict

        assert data['quantity'] is not None
        return data

    parsed_dict = {k: None for k in keys_to_parse}

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
    string = string.replace('[mmHg]', 'mmHg')

    # TODO maybe use regex for this. implement more general correctly handling
    # somewhere else? didn't i already have one hack for something like this?
    string = string.replace('/cu cm', '/mL')

    # Assuming we can remove all parenthetical stuff that remains
    # (source was already parsed and removed above) without issue.
    string = re.sub(r'\(.+\)', '', string)

    parts = string.strip().split()

    remaining_parts = []
    for p in parts:
        # TODO may need to strip adjacent paren / colon from part w/ quantity?
        try:
            curr_quantity = float(p.replace(',',''))
        except ValueError:
            # To match format like: 1.5X10+4
            # Assuming there will be exactly one digit in the exponent.
            sn_re = float_re + r'X10((?:\+|-)\d)'
            match = re.match(sn_re, p)
            if match:
                groups = match.groups()
                curr_quantity = float(groups[0].replace(',','') + 'e' + groups[-1])
            else:
                remaining_parts.append(p)
                continue

        if quantity is None:
            parts_before_quantity = list(remaining_parts)
            quantity = curr_quantity
        else:
            # TODO how to handle cases like:
            # https://pubchem.ncbi.nlm.nih.gov/compound/22311
            # "0.25 to 0.67 kPa at 20 °C (1.9 to 5 mm Hg)"
            # where there is redundant parentheical information?
            # can i always throw out parenthetical stuff, or would that
            # screw up some cases that previously worked?

            # TODO implement fn to check changes against all cached lookup values! (?)

            # TODO improve error message (esp clarify wrt similar error below, to avoid
            # need for '(1)' / '(2)' thing)
            # TODO what even is "remaining"? improve error message
            err_msg = (f'(1) multiple unit definitions:\n"{orig_str}"\n'
                f'remaining: {string}'
            )
            if url is not None:
                err_msg = f'{url}\n{err_msg}'

            raise ValueError(err_msg)

    # TODO TODO why not just raise here (and in other cases where null_dict is
    # returned)? (might simplify collecting issues with the parser, esp as i try to do
    # so for a greater diversity of inputs/properties) (w/ error message including which
    # other parts *had* been parsed, which might include the quantity that had been
    # mistaken for something else?)
    if quantity is None:
        return null_dict
        # TODO implement
        #raise ValueError('')

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

            try:
                # TODO may want to wrap all float parseing
                # (which includes scientific notation handling)
                # into a fn, so it can also handle scientific notation here!
                last_quantity = float(last_part[1:])

                # If previous part had a quantiy, then use second part of
                # this parenthetical expression as a unit.
                # Need to remove the end paren or it won't parse.
                p = p[:-1]

            except ValueError:
                pass
        last_part = p

        # TODO way to get pint to not parse things like this as concentration
        # units? some way to test if it was something like this? other things i
        # will have to specifically exclude?
        # TODO TODO what was this excluding?
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
            # TODO maybe warn / print if verbose here
            continue

        if last_quantity is not None:
            quantity = last_quantity

        # TODO TODO check dimensions of units match one of expected things
        # (or are at least not dimensionless, if no expected units provided)
        # TODO maybe store cache in chemutils so units can be checked consistent
        # against other lookups of same property, in case expected units not
        # specified

        if units is None:
            units = curr_units
        else:
            # TODO what even is "remaining"? improve error message
            err_msg = (f'(2) multiple unit definitions:\n"{orig_str}"\n'
                f'remaining: {remaining_str}'
            )
            if url is not None:
                err_msg = f'{url}\n{err_msg}'

            raise ValueError(err_msg)

    # One case that triggers this: "2 % (NIOSH, 2016)"
    if units is None:
        return null_dict

    # Since the measurement / reporting probably just doesn't have enough
    # sigfigs to report the real value, and 0.0 will cause any calculations
    # with multiplications to probably behave not as desired.
    if quantity == 0.0:
        return null_dict

    parsed_dict['units'] = units
    parsed_dict['quantity'] = quantity

    # TODO if expected_units is not None, convert to quantity to those units
    # (and also change / do not return unit str then)
    # (if i want to support sequence of expected units, would need to decide
    # what to do in that case. probably don't convert?)

    return parsed_dict


# TODO delete
# str -> dict
_str_parsings = dict()
# inchi -> list[str]
_parse_failures1 = dict()
_parse_failures2 = dict()

_no_data = set()
_total_failures = set()
# (successes)
_had_num_lookup = set()
_had_str_lookup = set()
#
def pubchemprops_lookup(inchi, prop_name, expected_units=None,
    target_temp_c=None, skip_unit_parse_err=False, verbose=False):
    # TODO TODO doc whether this averages over sources / how it picks one / etc
    """
    Requires my fork of pubchemprops which adds the `cid` kwarg to
    `get_second_layer_props`.
    """
    # TODO TODO still err if *no* strings can be parsed successfully, even if
    # skip_unit_parse_err=True (or some other value for that?)

    from pubchemprops.pubchemprops import get_second_layer_props

    # TODO also pass ignore_cache added by `cached` through so it's usable here?
    # TODO TODO also handle case where this is null
    cid = convert(inchi, from_type='inchi', to_type='cid')
    if cid is None:
        return None

    url = pubchem_url(cid)

    # TODO delete
    if verbose:
        print(f'LOOKING UP {prop_name} OF {inchi}')
        print(url)
    #

    # TODO TODO TODO catch and persistently log / store in pickle which inputs produce
    # which error[ codes]. also catch which produce JSONDecodeError / assert none do
    #
    # TODO where will i have to check for / return None in case of failed
    # lookup?
    ret = get_second_layer_props(cid, [prop_name], cid=True)

    if len(ret) == 0:
        # TODO delete
        _no_data.add(inchi)
        #

        # TODO or do i want to return null_dict? as long as this works w/ pandas
        # fns...
        return None

    # might not be true. if not, would need to also return None in this case
    assert len(ret[prop_name]) > 0

    parsed_dicts = []
    # TODO also aggregate and (optionally?) return sources of the data
    # TODO TODO maybe take closest to room temp / average across all
    for result in ret[prop_name]:
        value_dict = result['Value']

        # TODO if verbose, print which of these cases at least?

        if 'Number' in value_dict:
            # TODO use result['Reference'] for source in this case

            # I've only seen it have 'Unit' and 'Number' keys.
            assert len(value_dict) <= 2

            # TODO maybe just share from parsing fn above?
            keys = ['quantity', 'units', 'temperature_c', 'source', 'from_range']
            parsed_dict = {k: None for k in keys}

            # TODO TODO are there cases where i'd also want to share all of the unit
            # parsing logic above here?
            unit_str = value_dict['Unit']

            temp_c, unit_str = parse_temperature_c(unit_str)
            # `temp_c` can be `None`, but that will just overwrite w/ existing value.
            parsed_dict['temperature_c'] = temp_c

            unit_str = unit_str.replace('()','').strip()

            # This is to special case one instance where the unit_str was 'mg/mL, clear,
            # colorless'. Probably no cases where I want everything surrounding a comma,
            # though this may not always be right...
            if ',' in unit_str:
                unit_str = unit_str.split(',')[0].strip()
                assert len(unit_str) > 0, ('hack to fix unit str with extra'
                    ' info following comma failed'
                )

            units = ureg(unit_str)
            assert in_expected_units(expected_units, units)
            parsed_dict['units'] = units.u

            nums = value_dict['Number']
            assert len(nums) == 1
            parsed_dict['quantity'] = nums[0]

            parsed_dict['from_range'] = False

            parsed_dicts.append(parsed_dict)
            # TODO maybe also inspect these to sanity check this path isn't mangling the
            # quantities / units?

            # TODO delete
            _had_num_lookup.add(inchi)
            #

        elif 'StringWithMarkup' in value_dict:
            assert len(value_dict) == 1
            pstr = value_dict['StringWithMarkup'][0]['String']

            # TODO but add flag to include (some/all) of these (rather than doing our
            # own estimates)?
            # TODO TODO TODO also deal w/ 'calculated' same way (and check for other
            # things like 'est', etc)
            # '/Estimated/' exists in some places too. check 'est' in lower()?

            # This was to filter out '5.28X10+9 mm Hg at 25 °C /extrapolated/' and
            # anything like it, which is the only vapor pressure data on the page for
            # Lysine, but is absurd. Not sure if the actual reference actually has this
            # value (couldn't find online access).
            # https://pubchem.ncbi.nlm.nih.gov/compound/5962
            # TODO other strings to blacklist?
            if 'extrapolated' in pstr:
                warn(f'excluding StringWithMarkup="{pstr}" from consideration because '
                    'of sub-string "extrapolated", which was present in some other '
                    'problem data'
                )
                continue

            if verbose:
                print(f'trying to parse "{pstr}"')

            try:
                parsed_dict = parse_pubchemprops_string(pstr,
                    expected_units=expected_units, target_temp_c=target_temp_c,
                    url=url
                )

            # ValueError: multiple unit definitions. remaining:...
            # TODO TODO be more clear about what is happening here
            except ValueError as err:
                # TODO delete
                # TODO include the specific error message alongside?
                if inchi in _parse_failures1:
                    if pstr not in _parse_failures1:
                        _parse_failures1[inchi].append(pstr)
                else:
                    _parse_failures1[inchi] = [pstr]
                #
                if skip_unit_parse_err:
                    warn(str(err))
                    continue
                else:
                    raise

            if verbose:
                print(f'parsed: {pformat(parsed_dict)}')

            # TODO TODO why would parse_pubchemprops_string not just err in this case?
            # is this actually happening? (modify it so it does, and in all other cases
            # it is effectively returning nothing (`null_dict`))
            if parsed_dict['quantity'] is None:
                # TODO delete
                if inchi in _parse_failures2:
                    if pstr not in _parse_failures2:
                        _parse_failures2[inchi].append(pstr)
                else:
                    _parse_failures2[inchi] = [pstr]
                #
                continue
            else:
                parsed_dicts.append(parsed_dict)

                # TODO delete
                _had_str_lookup.add(inchi)
                # copying so that when we mutate dict to change units below, we still
                # have original units here (which is what we want to check the parsing)
                _str_parsings[pstr] = dict(parsed_dict)
                #

            assert len(value_dict['StringWithMarkup']) == 1
            # If there is an extra Key under the first element, it seems to have always
            # just been a reference to water...
        else:
            raise ValueError(f'{url}\nunexpected value dict')

    if len(parsed_dicts) == 0:
        # TODO delete
        if len(ret[prop_name]) > 0:
            _total_failures.add(inchi)
        #
        return None

    # TODO maybe consider this earlier to not continue parsing
    # anything that seems like a range if we already have non-range data
    if any([d['from_range'] == False for d in parsed_dicts]):
        parsed_dicts = [d for d in parsed_dicts if d['from_range'] == False]

    # TODO TODO maybe also convert all concentration units to same
    # dimensionality for water_solubility case, given known properties of water
    # (and I suppose also the chemicals... so maybe best left to be done
    # outside)
    for parsed_dict in parsed_dicts:
        m = parsed_dict['quantity']
        u = parsed_dict['units']
        uq = ureg.Quantity(m, u)

        found_dims = False
        if expected_units == 'density':
            preferred_unit = preferred_density_unit
        else:
            preferred_unit = None
            # TODO should i allow property2preferred_units to override this?
            # (would need to translate from prop_name here to my naming conventions for
            # same properties throughout this file)
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

    # TODO maybe check that dimensionality of input doesn't otherwise indicate it's a
    # concentration, if i'm only going to handle concentrations as i want when
    # expected_units is explicitly 'concentration'?
    # TODO move inside loop above?
    if expected_units == 'concentration':
        molar_mass = None
        for parsed_dict in parsed_dicts:
            m = parsed_dict['quantity']
            u = parsed_dict['units']
            uq = ureg.Quantity(m, u)

            # It must be a molar concentration then (assuming it's just between
            # that option and mass/vol)
            if '[substance]' in uq.dimensionality.keys():
                # TODO TODO if i cache this part, probably still thread the
                # ignore_cache flag through... (or maybe it global-ish somehow?)
                # (and do the same anywhere else there are nested cached fns)
                if molar_mass is None:
                    molar_mass = ureg.Quantity(
                        cid2molar_mass(cid)['molar_mass'], 'g/mol'
                    )

                parsed_dict['quantity'] = (uq * molar_mass
                    ).to(preferred_concentration_unit).magnitude
                parsed_dict['units'] = preferred_concentration_unit

    # TODO TODO outlier / within-order-of-magnitude detection here
    # TODO TODO + averaging / selection based on temperatures

    # TODO TODO maybe also prefer stuff that has a temperature to stuff that
    # doesn't... (assuming it's in a reasonable range about temp we want)

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
    })


# TODO hardcode these values:
# "Soluble in about 720 parts water, in many organic solvents" (l117 in test
# data)
# "0.43% (by wt) in water" (l71)
# "In water, 6,700 ppm at 25 °C" (l24)
# "Soluble in 120 parts water at 25 °C" (l25)
# maybe this and all the other NIOSH stuff
# "5 % (NIOSH, 2016)" (l39)
#@opt_pint_return
@cached
def inchi2water_solubility(inchi, verbose=False):
    ret_keys = [
        'water_solubility',
        'water_solubility_sources',
        'water_solubility_units',
    ]

    # TODO TODO maybe convert to str repr before caching so that i don't
    # have to worry about that pint app registry thing?
    # (see notes in docs about comparing stuff across ureg instances)
    # (or just test stuff as i'm doing it now... ?)

    # TODO also pass ignore_cache added by `cached` through so it's usable here?
    ret = pubchemprops_lookup(inchi, 'Solubility',
        expected_units='concentration', verbose=verbose
    )
    if ret is None:
        return pd.Series({k: None for k in ret_keys})

    # TODO TODO is this something that can be retrieved via chemspipy??
    # fallback to that if so.
    # TODO TODO NIST?

    return pd.Series({
        'water_solubility': ret['quantity'],
        'water_solubility_sources': ret['sources'],
        'water_solubility_units': ret['units'],
    })


#@opt_pint_return
@cached
# TODO include skip_unit_parse_err in wrapper to be applied to other fns using
# pubchemprops_lookup
# TODO TODO TODO try replacing w/ pubchempy / bare REST call to pubchem
# (wanna get rid of my pubchemprops usage / associated parsing code v badly)
# TODO TODO TODO try looking for other libraries / APIs for this
# TODO available in SDF files from pubchem? python library to parse those?
# (even if so, i doubt it's in an easier to parse format...)
def inchi2vapor_pressure(inchi, skip_unit_parse_err=True, verbose=False):
    ret_keys = [
        'vapor_pressure',
        'vapor_pressure_sources',
        'vapor_pressure_units',
    ]
    # TODO TODO TODO modify this so it fails fully if there are any input strings, but
    # ALL fail parsing (-> fix parsing for those cases).
    # TODO TODO TODO +count/report number of such cases
    ret = pubchemprops_lookup(inchi, 'Vapor Pressure', expected_units='kPa',
        skip_unit_parse_err=skip_unit_parse_err, verbose=verbose
    )
    # TODO delete (after experimenting w/ some more direct pubchem/rest/etc calls)
    _debug = False
    if _debug and inchi is not None:
        print()
        print(f'{inchi=}')
        print(f'url={inchi2pubchem_url(inchi)}')
    #

    # TODO TODO TODO try chemicalbook.com?
    # e.g. https://www.chemicalbook.com/ChemicalProductProperty_EN_CB5347184.htm
    # where are they getting their data though? at least it seems they indicate when
    # it's estimated?

    if ret is None:
        # TODO delete
        # TODO TODO TODO try to focus on ones specifically where parsing has an issue
        if _debug and inchi is not None:
            print('current lookup failed!')
            import ipdb; ipdb.set_trace()
        #
        return pd.Series({k: None for k in ret_keys})

    # TODO delete
    if _debug and inchi is not None:
        print(ret.to_string())
        import ipdb; ipdb.set_trace()
    #

    # TODO TODO is this something that can be retrieved via chemspipy??
    # fallback to that if so.
    # TODO TODO TODO NIST? (might need to use antoine equation parameters, and check
    # validity range first -> compute vapor pressure at stp)
    # TODO TODO other good databases available?

    return pd.Series({
        'vapor_pressure': ret['quantity'],
        'vapor_pressure_sources': ret['sources'],
        'vapor_pressure_units': ret['units'],
    })


_chemspider_session = None
# TODO TODO so is list of results just empty after i pass api call limit?
# some way to detect and notify about that? (b/c all but first few of
# a series of url lookups seemed to be null...)
@cached
def inchi2chemspider_id(inchi):
    global _chemspider_session

    # If chemspipy_api_key is None, trying to construct ChemSpider would give us:
    # chemspipy.errors.ChemSpiPyAuthError: Unauthorised; Invalid API Key;
    if chemspipy_api_key is None:
        raise RuntimeError('chemspipy_api_key must be set to use inchi2chemspider_id')

    if _chemspider_session is None:
        _chemspider_session = ChemSpider(chemspipy_api_key)

    cs = _chemspider_session

    # TODO factor assertion that it doesn't already start w/ inchi into @cached(/wrapper
    # that expands on what @cached does)
    assert not inchi.startswith('InChI=')

    # Docs seem to claim inchi could be passed just as if it were name,
    # but my testing did not seem to bear that out.
    inchi_with_prefix = 'InChI=' + inchi

    # (sleeping doesn't seem to really have an effect on how many calls we can get)
    # 1st key: (sleep=0.2) -319
    # 2nd key: (sleep=0.2) 319-653
    # 3rd key: (no sleep) 653-989
    #time.sleep(0.2)

    results = Results(cs, cs.filter_inchi, (inchi_with_prefix,), raise_errors=True)

    try:
        # TODO do the query status calls this makes also count toward api limit? that
        # would suck. max # calls should be time it takes query / 0.2 seconds, at least
        # (+ intial call)
        # TODO should this be in try block too?
        results.wait()

    # TODO this only error possible? can Results constructor above also throw or nah?
    except ChemSpiPyRateError as err:
        raise

    # TODO options here? should it always be 'Complete' here?
    assert results.status == 'Complete'

    if len(results) == 0:
        return None

    # TODO TODO TODO ran into "Too Many Requests" at 319/2618 InChI. does this mean it's
    # making ~3 API calls per search? can i modify ChemSpiPy to make it just one call
    # per (did old synchronous code do that?)?

    # TODO TODO warn if there are multiple? store all somewhere else (to use others w/o
    # having to make more previous API calls, if others would be useful)
    csid = min([r.record_id for r in results])
    return csid


# TODO TODO maybe return url that achieves search for inchi if csid is null
# (particularly if that can happen just b/c api call limit exceeded or
# something)? possible?
def inchi2chemspider_url(inchi, **kwargs):
    # NOTE: this is @cached, so no need to cache inchi2chemspider_url
    csid = inchi2chemspider_id(inchi, **kwargs)
    if csid is None:
        return None

    return chemspider_url(csid)


#@opt_pint_return
@cached
def inchi2density(inchi, verbose=False):
    ret_keys = ['density', 'density_sources', 'density_units']
    null_ret = pd.Series({k: None for k in ret_keys})

    # TODO is this going to be in consistent units? convert?
    # pass expected_units?
    pubchem_dict = pubchemprops_lookup(inchi, 'Density',
        expected_units='density', verbose=verbose
    )
    if pubchem_dict is not None:
        return pd.Series({
            'density': pubchem_dict['quantity'],
            'density_sources': pubchem_dict['sources'],
            'density_units': pubchem_dict['units'],
        })

    # TODO maybe cache these separately if we can't tell where None came from?
    if chemspipy_api_key is None:
        warn('could not lookup density on chemspider b/c no chemspipy_api_key set')
        return null_ret

    # TODO in the future, only warn if any are missing (?)

    # TODO this inchi guaranteed to be non-null?
    # TODO pass ignore_cache through to this if applicable
    # TODO TODO again maybe ignore_cache should just be set globally to make
    # this possible?
    # TODO TODO in future, look for availability of properties in any associated
    # compounds, and average anything found.
    # (for now, just using lowest compound ID as a hack)
    csid = inchi2chemspider_id(inchi)

    if csid is None:
        if verbose:
            print('COULD NOT FIND INCHI IN CHEMSPIDER')

        return null_ret


    # TODO could enable behind some flag for extra checks
    #result = [r for r in results if r.csid == csid][0]
    #assert result.inchi == inchi_with_prefix

    # Since I couldn't find a way to get arbitrary properties like
    # "Experimental Density" through the API.
    cs_dict = chemspider_experimental_density(csid, verbose=verbose)
    if verbose:
        if cs_dict is not None:
            print('CHEMSPIDER DENSITY:', cs_dict['density'])
        else:
            print('NO CHEMSPIDER DENSITY')

    # TODO why does an inchi search sometimes return multiple?
    # possible to just restrict to exact match after the fact? what else to do?
    # (there can apparently be multiple compounds w/ same inchi)
    # e.g. 11709 and 21231965
    # TODO maybe assert that if there are multiple results, they all have same
    # inchi? (could add lookup time, but would see whether exact matching on
    # inchi would ever be useful)

    if cs_dict is None:
        return null_ret
    else:
        return cs_dict


# TODO get temp dependence coeff if available
# TODO assert temp it's defined at is same-ish (or is it always same?)?

#@opt_pint_return
@cached
def inchi2nist_k_henry(inchi, verbose=False):
    ret_keys = ['k_henry', 'k_henry_sources', 'k_henry_units']

    # TODO TODO refactor
    null_dict = pd.Series({k: None for k in ret_keys})

    if inchi is None:
        return null_dict

    # TODO try just getting from inchi2nist_webbook_url (table should be below
    # anyway, if it exists. then i'd need to follow it though...). i dont think
    # inchi2nist_henrys_law_url works anymore
    url = inchi2nist_henrys_law_url(inchi)
    if verbose:
        print(url)

    # TODO TODO TODO add retry logic (see one of my changes to pubchemprops)
    response = session.get(url)

    response.raise_for_status()

    # TODO TODO TODO probably def my own read_html fn (wrapping pd.read_html), that uses
    # retry logic and response.text, rather than passing URL directly

    # TODO test in case where there shouldn't be any tables
    #dfs = pd.read_html(url,
    try:
        dfs = pd.read_html(response.text,
            attrs={'aria-label': "Henry's Law constant (water solution)"}
        )

    # ValueError: No tables found
    except ValueError:
        # TODO assert error message matches expected

        # TODO TODO TODO refactor so that i can just return None here
        # (rather than needing to make a Series w/ the right keys)
        return null_dict

    '''
    if len(dfs) == 0:
        # TODO TODO TODO refactor so that i can just return None here
        # (rather than needing to make a Series w/ the right keys)
        #return None
        return null_dict
    '''

    assert len(dfs) == 1
    df = dfs[0]

    # TODO TODO average only w/in highest quality 'M' measurements (in 'Method' column)?
    # prioritize those that include temperature dependence, if available?  other
    # considerations?

    # TODO check these units (copied 2023-05-30) are still reflected in old preferred
    # units i had (+ that curr units are assigned correctly before any conversion):
    # "k°H (mol/(kg*bar))"

    # TODO TODO throw out stuff w/ 'missing [citation]' in df.Comment?

    # TODO delete try/except
    try:
        kh = df['k°H (mol/(kg*bar))'].mean()
    # TypeError: Could not convert >6.0×10+8<4.0×10+1160000. to numeric
    except TypeError:
        # TODO TODO deal with: https://webbook.nist.gov/cgi/cbook.cgi?InChI=1S%2FC3H8O3%2Fc4-1-3(6)2-5%2Fh3-6H%2C1-2H2&Units=SI&Mask=10
        # 'M' should probably also not be used here b/c '...unreliable'
        # also missing citation.
        # TODO TODO TODO fix!
        #import ipdb; ipdb.set_trace()
        print('NIST HENRY CONSTANT GETTING IGNORED! FIX!')
        return null_dict
    #

    return pd.Series({
        'k_henry': kh,
        # TODO try to include actual source?
        'k_henry_sources': url,
        # TODO check pint would parse this string correctly
        'k_henry_units': 'mol / kg / bar',
    })
    # TODO TODO TODO delete below

    response = session.get(url)
    # TODO maybe factor all rate limiting into another decorator
    # (would want to factor out loop portion below to a fn then probably...)
    # crude rate limit
    time.sleep(0.5)
    html = response.text
    soup = BeautifulSoup(html, features='html.parser')

    table = None
    # TODO check this condition is only True if table would be None
    if '<h2>Matches found</h2>' in html:
        inchi_str = 'InChI=' + inchi
        if verbose:
            print(len(soup.find('ol').find_all('a')), 'matches')

        # Should only be one of these <ol> tags in the html.
        for match in soup.find('ol').find_all('a'):
            href = match.attrs['href']
            url = f'{nist_url_prefix}{href}'
            if verbose:
                print(url)

            response = session.get(url)
            # crude rate limit
            time.sleep(0.5)
            html = response.text
            soup = BeautifulSoup(html, features='html.parser')
            inchi_matches = soup.find_all('span', string=inchi_str)
            if len(inchi_matches) == 0:
                continue

            if verbose:
                print(len(inchi_matches))

            assert len(inchi_matches) == 1
            inchi_match = inchi_matches[0]
            # TODO test substrings dont match
            if verbose:
                print(inchi_str)
                print(inchi_match.text)

            table = soup.find('table',
                attrs={'aria-label': "Henry's Law constant (water solution)"}
            )
            if table is not None:
                print('NEW URL:', url)
                break

    # TODO TODO identify + take user input / alt ID to select URL in this case
    # "There are no matching entries in the database for this IUPAC
    #  International Chemical Identifier"
    # TODO it seems even if InChI is among options on page returned above,
    # it still sometimes yields that page for disambiguation...
    # way to force the inchi i want?
    # seems like stuff ending in '+' might not have a permanent url for some
    # reason... related?
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C6728263&Units=SI and cf
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C1335393&Units=SI

    else:
        if verbose:
            print('MATCH WAS UNAMBIGUOUS')

        # TODO units always mol/kg*bar (probably)?
        table = soup.find('table',
            attrs={'aria-label': "Henry's Law constant (water solution)"}
        )

    if table is None:
        if verbose:
            print('NO HENRY LAW DATA\n')

        # TODO TODO check this doesn't break any possible places
        # i assumed things w/ series output in some cases always have series
        # output...
        return None

    # TODO ever have multiple w/ best method and presence of ref?
    # need to resolve across them?
    kh = None
    last_had_ref = False
    last_measured = False
    for row in table.find_all('tr', class_='exp'):
        cols = row.find_all('td')
        # First column will have Kh (second is temp. dep. part, if there)
        kh_str = cols[0].text
        try:
            curr_kh = float(kh_str)
        except ValueError as e:
            # This branch triggered by at least 'Methyl Isobutyl Ketone'
            parts = kh_str.split()
            if parts[1] == '-':
                kh_min = float(parts[0])
                kh_max = float(parts[2])
                curr_kh = kh_min + (kh_max - kh_min) / 2
            else:
                raise

        measured = cols[2].text == 'M'
        had_ref = cols[3].text != 'N/A'

        # Just taking first measured (not estimated / unclear origin)
        # value with a reference (if we have one).
        if measured and had_ref:
            kh = curr_kh
            break

        elif measured:
            kh = curr_kh
            last_measured = True
            last_had_ref = False

        # Being measured takes priority over having a reference.
        elif not last_measured and had_ref:
            kh = curr_kh
            last_had_ref = True

        elif not last_measured and not last_had_ref:
            kh = curr_kh

    # TODO (optionally?) also return temp (+P?) Kh defined at, any temp
    # dependence coefficient, method, and reference (use dict or pd.Series)
    # + service / database if using options besides NIST
    # (also implement in inchi2density)

    if verbose:
        print('Kh:', kh)
        print()

    # TODO if possible, maybe factor return of all these prop lookup fns
    # to fill in property name and prefix sources and units cols w/
    # that property name?
    return pd.Series({
        'k_henry': kh,
        'k_henry_sources': url,
        # TODO check pint would parse this string correctly
        'k_henry_units': 'mol / kg / bar',
    })


# TODO update to remove 'mp'/max planck references (now that it's just henrys-law.org)
mp_henry_prefix = 'https://henrys-law.org/henry/'
def inchikey2maxplanck_k_henry_url(inchikey):
    # TODO TODO TODO follow link in search results?
    return f'{mp_henry_prefix}search_identifier.html?search={inchikey}'


# TODO when doing my own henry's law constant unit conversions, check them against the
# calculator on this website
# TODO conversion to inchikey best done in here or outside?
# TODO uncomment decorators after getting this to work
# TODO TODO also consider using one of the two downloads available on the website (.sql
# or .f90, both kind of weird formats, to generate my own CSV / similar representation)
# (could be useful for fitting a model...)
# TODO add exclude_types kwarg for e.g. excluding estimates / QSPR / VP/AS / etc
#@opt_pint_return
@cached
def inchi2maxplanck_k_henry(inchi, verbose=False):
    """
    Looks up Henry's law constant at: https://www.henrys-law.org
    http://satellite.mpic.de/henry (redirected from http://www.henrys-law.org).

    2023-05-30:
    it seems henrys-law.org is no longer redirecting to satellite.mpic.de/henry

    The site was made by Rolf Sander.

    It says on its homepage that data from temperatures above 373K (~100C)
    were not used, but that might still mean some data is at temperatures
    somewhat different from what we want.
    """
    # TODO the PDF compilations he publish have anything the site doesn't?
    # TODO anyway to get a download of the data? request?
    import requests
    # TODO move up top + add to setup.py
    from bs4 import BeautifulSoup

    ret_keys = ['k_henry', 'k_henry_sources', 'k_henry_units']

    # TODO TODO refactor
    null_dict = pd.Series({k: None for k in ret_keys})

    # TODO deal with this in a wrapper ideally (maybe just in cached? rename at that
    # point?)
    if inchi is None:
        return null_dict

    assert not inchi.startswith('InChI=')
    # Returns None on error
    inchikey = InchiToInchiKey(f'InChI={inchi}')
    assert inchikey is not None

    # TODO cf
    # https://henrys-law.org/henry/search_identifier.html?search=PEDCQBHIVMGVHV-UHFFFAOYSA-N
    # -> https://henrys-law.org/henry/casrn/56-81-5
    # w/ http://webbook.nist.gov/cgi/cbook.cgi?InChI=1S%2FC3H8O3%2Fc4-1-3(6)2-5%2Fh3-6H%2C1-2H2&Units=SI&Mask=10

    search_url = inchikey2maxplanck_k_henry_url(inchikey)

    response = session.get(search_url)
    html = response.text

    # This parser is case insensitive (useful here).
    # Requires: pip install lxml
    soup = BeautifulSoup(html, 'xml')
    tags = soup.findAll('h1')
    # Website must have changed if this happens.
    assert len(tags) == 1
    tag = tags[0]
    n_results = int(tag.text.split()[0])
    assert n_results >= 0
    if n_results == 0:
        return null_dict

    if n_results > 1:
        #import ipdb; ipdb.set_trace()
        # TODO prob convert to error
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

    # TODO want to sleep between this request and the above?

    # TODO TODO use his (what seems like) temperature dependence data (but how?) (same
    # w/ nist stuff still...)
    # TODO how do i know what temp Hcsp is defined at? vary? standardized somehow?
    # (the "d ln Hcp / d (1/T)" column)

    # TODO maybe download the whole page to some webpage cache, and then
    # reprocess stuff from there (so i can change my criteria / average across
    # different subsets of table rows) w/o having to make more requests?
    # (maybe some idiomatic way to make these kinds of caches?)
    # (and would probably want to do w/ all things i scrape manually w/
    # beautifulsoup)

    response = session.get(url)
    html = response.text

    dfs = pd.read_html(html)

    # Looking for one DataFrame with an index like this:
    # [MultiIndex([(                  'Hscp',       '[mol/(m3Pa)]'),
    #             ('d ln Hs  cp /  d (1/T)',                '[K]'),
    #             (            'References', 'Unnamed: 2_level_1'),
    #             (                  'Type', 'Unnamed: 3_level_1'),
    #             (                 'Notes', 'Unnamed: 4_level_1')],
    #            )]
    # TODO normalize this second column / drop (i think it has some weird whitespace in
    # it)
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

    # NOTE: see R. Sander, 2015 for more details.
    # https://doi.org/10.5194/acp-15-4399-2015
    type2description = {
        'L': 'literature review',
        'M': 'measured',

        # TODO TODO just use this method myself to estimate missing values?
        'V': 'VP/AS = vapor pressure/aqueous solubility',

        # TODO how is this done? in any of our VCF chemicals?
        'R': 'recalculation',

        # TODO TODO how is this done?
        'T': 'thermodynamical calculation',

        'X': 'original paper not available',

        # TODO TODO what does this mean???
        'C': 'citation',

        'Q': 'QSPR',
        'E': 'estimate',
        '?': 'unknown',
        # TODO TODO TODO exclude at least this!
        'W': 'wrong',
    }

    type_set = set(df.Type)
    assert all(x in type2description.keys() for x in type_set), \
        f'unknown types in {type_set}'

    # TODO go through notes page and find any references that seem like they shouldn't
    # be used? https://henrys-law.org/henry/notes.html

    # TODO TODO generally warn if entries differ by orders of magnitude?
    # (will probably often be the case)

    # TODO TODO select first row? avg? should be listed in decreasing order of
    # reliability
    # TODO maybe avg within first value for 'Type' col?
    row = df.iloc[0]
    Hscp = row.Hscp
    # TODO there might be other places where there are '\xa0' characters (non-breaking
    # space) that we might like to remove too
    short_ref = row.References.replace('\xa0', ' ')

    # TODO TODO TODO TODO where is these values coming from? can't be right...
    # TODO delete
    # TODO check for other big/crazy values (inspect biggest / smallest)
    sus_strs = [
        # https://henrys-law.org/henry/casrn/87-89-8
        # Estimate from Saxena and Hildemann (1996)
        # "using the group contribution method"
        # https://doi.org/10.1007/BF00053823
        '9.9×1023',

        '1.2×1014',

        # TODO TODO which to use here? either?
        # https://henrys-law.org/henry/casrn/56-40-6
        # just V -> 1.2e11, E (estimate) -> 8.9e5
        '1.2×1011',

        # https://henrys-law.org/henry/casrn/51-43-4
        '1.4×1013',
    ]
    if Hscp in sus_strs:
        print()
        print(url)
        print(f'sussy {Hscp=}\n')
    #

    # TODO see if any non-ambiguous data available in other rows for any of these
    # TODO TODO does rolph say why he indicates '>' for data from: Altschuh et al.
    # (1999) (seems like that's where all this stuff is coming from)
    # https://doi.org/10.1016/S0045-6535(99)00082-X
    # ('Altschuh' is not in notes page. any note # references in other places Altschuh
    # data is mentioned?)
    #
    # seems like it might make most sense to just discard '>'...
    #
    # For many of these, other values are lower despite '>' sign and higher
    # "reliability" on value with inequality...
    #
    # Triggered by:
    # https://henrys-law.org/henry/casrn/60-12-8
    # (other data, of diff types, available and close)
    #
    # https://henrys-law.org/henry/casrn/3360-41-6
    # (other data, of diff types, available and close)
    #
    # https://henrys-law.org/henry/casrn/122-97-4
    # (other data, of diff types, available and close)
    #
    # https://henrys-law.org/henry/casrn/100-51-6
    # measurements AFTER the next two (type V from Mackay), are similar
    # TODO blacklist mackay 1995/2006c sources / mackay in general?
    # how else to handle?
    if type(Hscp) is str and ('<' in Hscp or '>' in Hscp):

        gt_ref = 'Altschuh et al. (1999)'
        if short_ref == gt_ref:
            assert '<' not in Hscp
            warn(f"dropping '>' from {gt_ref} value='{Hscp}'\n{url}\n")
            Hscp = Hscp.strip('> ')
        else:
            warn(f"not using data with inequality: '{Hscp}'\n{url}\n")
            return null_dict

    # TODO TODO should i throw out all type='V'
    # ("vapor pressure/aqueous solubility") estimates?

    try:
        # This will handle both exising floats as well as strings immediately parseable
        Hscp = float(Hscp)

    except ValueError:
        # float(...) doesn't work with strings containing the more exotic left character
        Hscp = Hscp.replace('−', '-')

        # TODO TODO sanity check that stuff w/ negative exponent should actually be
        # negative. really small ones especially.
        exp_prefix = '×10'
        assert exp_prefix in Hscp
        # Strings where the exponent is x10^1 should appear here as 'x101'
        # The original HTML has <sup>...</sup> tags around exponent.
        assert not Hscp.endswith(exp_prefix)

        mantissa, power = Hscp.split(exp_prefix)
        mantissa = float(mantissa)
        power = int(power)

        Hscp = float(f'{mantissa}e{power}')

        Hscp2 = mantissa * 10**power
        # NOTE: seems precision issues can make this not exactly equal, and i figure
        # using the float(...) method is slightly more what we want in those cases
        assert np.isclose(Hscp, Hscp2), f'{Hscp} != {Hscp2}'

    # TODO get full reference from list below table?
    source = f'{short_ref} ({row.Type})'
    notes = row.at['Notes']
    if not pd.isna(notes):
        # TODO TODO TODO regenerate with this!!!
        source += f'\nnotes: {notes}'
    source += f'\n{url}'

    return pd.Series({
        'k_henry': Hscp,
        'k_henry_sources': source,
        # TODO check pint would parse this string correctly
        'k_henry_units': units,
    })


# TODO TODO TODO come up with a general solution for handling multiple sources
# separately in their own fns, but having one central fn that can select between them.
# property2preferred_units / etc should all be agnostic to the source of the data.
# pretty sure i want the caching to happen on each source-specific lookup, so that i can
# more flexibly change the priority of the sources, etc.
#
# TODO TODO TODO also parse this from pubchem now. it is sometimes available.
# (prob need to edit some things inside my pubchemprops fork tho)
# TODO TODO k_henry data ever available in chemspipy (or chemspider more
# generally...)? doesn't seem so, BUT it does have estimates in 'Predicted - EPISuite'
# tab (VP here too)
def inchi2k_henry(inchi, verbose=False, **kwargs):
    # TODO TODO TODO add column to show which of the two lookup fns is used
    # (and maybe build some machinery to handle this in general, when there are multiple
    # source-specific lookup fns for a given property)
    keys = ['k_henry', 'k_henry_sources', 'k_henry_units']

    # TODO TODO TODO also make a source-specific fn for pubchem
    # (modifying pubchemprops fork to handle this), and use here
    # (ever have data when one of these two doesn't?)

    # TODO TODO option to average all outputs?
    # TODO if verbose, warn if they differ much?
    for fn in (inchi2maxplanck_k_henry, inchi2nist_k_henry):

        ret = fn(inchi, verbose=verbose, **kwargs)
        assert set(ret.keys()) == set(keys)
        assert isinstance(ret, pd.Series)

        # Taking first non-null value (so order of fns about is priority).
        if pd.notna(ret['k_henry']):
            break

    # Last call should give us null Series in this case
    return ret


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
            # TODO uncomment after lockout expires / after fixing issue
            '''
            if 'chemspider_url' not in df.columns:
                #df['chemspider_url'] = df.inchi.map(inchi2chemspider_url)
                df['chemspider_url'] = df.inchi.map(lambda x:
                    inchi2chemspider_url(x, ignore_cache_null=True,
                    leave_nonnull_cache=True)
                )
            '''

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

        def color_values_from_lookup(val):
            # Empty cells would otherwise still have the text this color
            # when it *is* manually entered.
            color = 'black' if pd.isnull(val) else 'gray'
            return f'color: {color}'

        print(f'Writing {name}to {excel_fname} for manual entry of missing '
            'properties.'
        )
        writer = pd.ExcelWriter(excel_fname, engine='xlsxwriter')

        styled = df.style.applymap(color_values_from_lookup, subset=props)
        #styled = df
        styled.to_excel(writer, index=False)

        # A dict of sheet name -> Worksheet objects.
        assert len(writer.sheets) == 1
        # Default sheet name is Sheet1
        worksheet = writer.sheets['Sheet1']

        workbook = writer.book
        fmt = workbook.add_format()

        # TODO some way to say scientific only for values small / large enough?
        # i want that?
        # 'Scientific' does not produce meaningful number outputs, though
        # docs seemed to suggest that at least something like that would work
        #fmt.set_num_format('Scientific')

        # TODO TODO some way to only show a certain number of sigfigs???
        # (with or without scientific format)
        # TODO also don't set this format for columns just that have integer
        # values (like in_dweck_active)
        # TODO TODO why is this not working anymore?
        # (i think it just will not apply to any stuff pandas has already styled
        # (bug?), so i'll need to implement my formatting logic inside calls to
        # write individual cells... (possible?) (or missing value conditional
        # text color and still use set_column?))
        fmt.set_num_format('0.000')

        # TODO potentially see worksheet.write_url as for a way to eliminate
        # default red formatting of urls
        # TODO and the FAQ says it can't autofit column widths, but it does
        # seem to indicate you can determine width of text after writing:
        # how to do this with pandas???
        # (for now, using write of str repr of input as proxy, which may fail
        # badly for floats)

        # TODO also turn off going beyond bounds of column for these?
        # (like when next cols in same row are null)
        # (and then move url_cols back before property_cols)
        fixed_col_widths = {
            'chemspider_url': len('chemspider_url'),
            'nist_webbook_url': len('nist_webbook_url'),
        }
        max_col_widths = df.applymap(lambda x: len(str(x))).max(axis=0)
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
            worksheet.set_column(r, w, fmt)

        writer.save()

        # TODO TODO automated checks that all things (non-nan) are either equal
        # or np.allclose (i checked manually and it's true)
        #rt = pd.read_excel(excel_fname, engine=pandas_excel_engine)

    else:
        print(f'Reading {name}property data from', excel_fname)
        df = pd.read_excel(excel_fname, engine=pandas_excel_engine)
        add_properties_to_cache(df, verbose=verbose)

    return df


# TODO may want to store record of what all was manually entered,
# so as to be able to treat that stuff differently
# (maybe ignore_cache shouldn't always apply to that stuff?)
def add_properties_to_cache(df, overwrite_close=False, dropna=True,
    verbose=False):
    """Adds chemical properties in `df` to cache under 'inchi'.
    """
    # So we don't have to check whether values behind each instance
    # of same inchi are consistent.
    assert not df.inchi.duplicated().any(), 'inchi can not have duplicates'

    df = df.set_index('inchi')

    add_count = 0
    for prop in properties:
        if prop in df.columns:
            preferred_unit_str = property2preferred_units[prop]
            preferred_units = ureg(preferred_unit_str)
            inchi2prop_cache = cache['inchi'][prop]
            # TODO this work? test
            def add_to_cache(i, p, row):
                nonlocal add_count

                # TODO maybe move unit conversion (currently in loop below)
                # up here, so it's clear we are converting, and not just setting
                # the unit_str to something different...

                # TODO maybe a fn for making these? (including null ones)
                new_entry = pd.Series({
                    prop: p,
                    f'{prop}_sources': getattr(row, f'{prop}_sources'),
                    f'{prop}_units': preferred_unit_str
                })

                # TODO should also make sure that whitespace only strings aren't
                # accepted (maybe check it's a url?)
                # If I don't want this error condition, should support adding
                # sources after previous entry of a value.
                if pd.isnull(getattr(row, f'{prop}_sources')):
                    raise ValueError(f'{row.Index} {prop} missing source')

                inchi2prop_cache[i] = new_entry
                is_manual[('inchi', prop, i)] = True
                add_count += 1
            #

            if verbose:
                print(f'adding {prop}')

            cols = [prop, f'{prop}_sources', f'{prop}_units']
            assert all([c in df.columns for c in cols[1:]])

            # TODO probably factorize out conversion of units on a df
            # (also want to do something similar in fn to load hardcoded manual
            # properties, though that spreadsheet is in a different format)

            for row in df[cols].dropna(subset=[prop]).itertuples():
                # TODO TODO maybe change so that i can more easily clear things
                # already entered (by deleting that row in spreadsheet and
                # re-adding) (delete existing stuff in cache if valueerror,
                # before raising (and if `is_manual` is True for it)?)

                # TODO for all these errors, try to use name rather than inchi
                # (at least if available)

                # Since usual dict key indexing doesn't work w/ itertuples
                # objects.
                unit_str = getattr(row, f'{prop}_units')
                if pd.isnull(unit_str):
                    raise ValueError(f'{row.Index} {prop} missing units')

                if not preferred_units.check(unit_str):
                    raise ValueError(f'{row.Index} {prop} units of {unit_str} '
                        f'inconsistent with {preferred_unit_str}'
                    )

                # TODO TODO may need extra lookup of molar mass if want to
                # convert between concentration units of different
                # dimensionality...

                # TODO probably vectorize this numeric check
                m = getattr(row, prop)
                try:
                    float(m)
                except ValueError:
                    raise ValueError(f'{row.Index} {prop} not a number')

                q = (m * ureg(unit_str)).to(preferred_unit_str).magnitude

                inchi = row.Index
                if inchi not in inchi2prop_cache:
                    if verbose:
                        print(f'adding {inchi} not already in cache')
                    #inchi2prop_cache[inchi] = new_entry
                    #add_count += 1
                    add_to_cache(inchi, q, row)
                    continue

                # TODO replace w/ isnull(inchi2prop_cache[inchi]), or no?
                if pd.isnull(inchi2prop_cache[inchi][prop]):
                    if verbose:
                        print(f'adding {inchi} null in cache')
                    add_to_cache(inchi, q, row)

                if not np.isclose(inchi2prop_cache[inchi][prop], q):
                    if verbose:
                        print(f'adding {inchi} b/c not close to cache value')
                    add_to_cache(inchi, q, row)
                    #inchi2prop_cache[inchi] = new_entry
                    #add_count += 1
                else:
                    if verbose:
                        print(f'not adding {inchi} b/c close to cache value')

    if add_count == 0:
        print('Nothing new found to add to cache.')

    # So that calls clearing cache before atexit trigger / problems with
    # atexit trigger do not prevent these values from being saved.
    # TODO want merge_with_saved=False here? (may not)
    save_cache()


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
        # TODO use pandas interface (if not factoring out whole gsheet reading
        # thing)
        with open(gsheet_cache_file, 'rb') as f:
            df = pickle.load(f)
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

        # TODO use pandas interface (if not factoring out whole gsheet reading
        # thing)
        with open(gsheet_cache_file, 'wb') as f:
            pickle.dump(df, f)

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

    # TODO probably just strip all string entries in df?
    df.abbrev = df.abbrev.apply(lambda s: s.strip() if pd.notnull(s) else s)

    df['original_name'] = df.name.copy()

    df.name = df.name.apply(normalize_name)

    df = convert(df, to_type='inchi', allow_nan=False, verbose=verbose,
        report_conflicts=False
    )

    # TODO could copy name to original_name and normalize ids (+ name)
    # as in natural_odors/odors.py
    # but may want to preserve

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

    # TODO check this doesn't break other place i use odor2abbrev
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
    # TODO maybe support other type inputs for convenience:
    # - non-unique list to be unique'd while preserving order
    # - others? ask lab mates?

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
                    # TODO make sure this prints for each odor,
                    # rather than somehow getting surpressed
                    warn(f'no abbreviation in odor inventory gsheet for odor "{o}"!')
                    # consistent case be damned here.
                    abbrev = chr(ord('A') + i)
                    i += 1

                o2a[o] = abbrev
    return o2a

