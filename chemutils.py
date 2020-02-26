
# So that raw string with unicode degree symbol for parsing temperatures
# could also work in Python 2.
from __future__ import unicode_literals

import re
import os
from os.path import split, join, exists
import pickle
import atexit
import urllib.error
import warnings
from collections import Counter
import traceback
import functools
import inspect

import numpy as np
import pubchempy as pcp
import pandas as pd
# TODO maybe also install and use the pandas support module for this.
# may need to do conversion of whole pandas objects in here than, rather
# than just defining conversion fns (that use pint) on single elements...
import pint
# Should be included with pint. Some pint fn can throw this.
from tokenize import TokenError


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

cache_file = os.path.expanduser('~/.chemutils_cache.p')

# In descending order of the certainty each ID provides as to whether we found
# the right PubChem entry.
chem_id_types = [
    'cid',
    'inchikey',
    'inchi',
    'smiles',
    'cas',
    'name'
]
equivalent_col_names = {
    'odor': 'name',
    'cas_number': 'cas'
}
# No matter the type these are to be converted to, they will always return a
# null value. (maybe handle some other way?)
manual_type2null_keys = {
    'name': {
        'spontaneous firing rate',
        'odor'
    }
}
extra_lookup_keys = [
    'density',
    'k_henry',
    'vapor_pressure',
    'water_solubility'
]

# TODO TODO TODO flag to ignore manual overrides, to the extent that they are
# causing problems. maybe phase them out altogether?

# TODO TODO try deleting / ignoring this hardcoded stuff and see if it still
# works. my changes to normalize_name probably fixed a few of these cases.
# (and if it can be ignored, delete it)
# TODO TODO TODO modify so these all still apply when ignore_cache=True
# (have that just apply to values that originally came from a lookup)
# (maybe have another flag / another level of ignore_cache that also applies to
# these values)
manual = {
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


# TODO TODO conversion to names with a priority on various sources:
# [0) hardcoded stuff]
# 1) names in inventory
# 2) hallem
# 3) other (good default available in pubchempy? probably not iupac...)


def add_manual_overrides(_cache):
    for ft, null_keys in manual_type2null_keys.items():
        for nk in null_keys:
            for tt in chem_id_types:
                _cache[ft][tt][nk] = None

    for ft, tts in manual.items():
        for tt in tts:
            _cache[ft][tt].update(tts[tt])


def init_cache():
    _cache = {ft: {tt: dict() for tt in chem_id_types + extra_lookup_keys}
        for ft in chem_id_types
    }
    add_manual_overrides(_cache)
    return _cache


def clear_cache():
    global cache
    cache = init_cache()


def _clear_null_cache_vals(write=True):
    global cache
    for from_type, to_type_to_cache_dict in cache.items():
        for to_type, f2t_dict in to_type_to_cache_dict.items():
            # Need list() since .items() is an iterator that expects contents
            # not to change.
            for f, t in list(f2t_dict.items()):
                if pd.isnull(t):
                    del f2t_dict[f]

    if write:
        save_cache(merge_with_saved=False)


def delete_cache():
    """Deletes and clears the cache at cache_file
    """
    clear_cache()
    if exists(cache_file):
        os.remove(cache_file)


def load_cache():
    if exists(cache_file):
        try:
            with open(cache_file, 'rb') as f:
                _cache = pickle.load(f)

            # TODO delete (was planning to use to load inverted maps here)
            '''
            for from_type, to_type_dict in _cache.items():
                for to_type, f2t_cache in to_type_dict.items():
            '''
            #

            # TODO shouldn't this logic (also, w/o refactoring) be in
            # init_cache?
            for k in extra_lookup_keys:
                if k not in _cache['inchi']:
                    _cache['inchi'][k] = dict()
            return _cache
        except ValueError:
            print('Cache was in unreadable format. Deleting.')
            delete_cache()
            # Will have already been init'd in delete_cache call.
            return cache
    else:
        return init_cache()


def save_cache(merge_with_saved=True):
    """
    Load old cache first, so that clear_cache doesn't have to worry about
    screwing up on-disk cache when atexit call to save_cache triggers.
    """
    if merge_with_saved:
        # TODO need "global cache", as long as cache def is below fn def, or
        # what?  if that's all, move "cache = load_cache()" just above this, and
        # declare as global in load_cache, cause circularity
        old_cache = load_cache()

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

    with open(cache_file, 'wb') as f:
        pickle.dump(to_save, f)


# TODO provide a "private" fn to delete single entries from cache
# TODO and another to manually enter persistent overrides

# TODO test
cache = load_cache()
atexit.register(save_cache)


# TODO maybe cache the result of this fn? when most values are cached,
# this might constitute a big relative portion of lookup time...
def allowed_kwarg(fn, kwarg_name):
    """Returns whether `kwarg_name` can be passed as a keyword argument to `fn`.
    """
    args, varargs, varkw, defaults = inspect.getargspec(fn)
    if varkw:
        return True

    # I guess something doesn't need a default to be passable as a kwarg...
    return kwarg_name in args


# TODO use this for basically every lookup fn?
# or still let convert handle caching for some of the functions?
# way to have it refer to the unwrapped versions in those cases?
# maybe just don't use decorator syntax at that point, and define new
# cached / un-cached versions of same fns, using same wrapper as decorator?
def cached(fn):
    # Would need to modify if I wanted kwargs for the decorator that explicitly
    # specify the from / to types.
    # https://stackoverflow.com/questions/627501
    f, t = fn.__name__.split('2')
    assert f in cache and t in cache[f]

    # TODO maybe modify so manual hardcoding into extra_lookup_keys
    # not required for fns referencing new "to" types?

    # So this cached decorator doesnt change function __name__
    @functools.wraps(fn)
    def cached_fn(f_val, ignore_cache=False, **kwargs):
        ft_cache = cache[f][t]

        verbose = False
        if 'verbose' in kwargs:
            if kwargs['verbose']:
                verbose = True
            
            if not allowed_kwarg(fn, 'verbose'):
                kwargs.pop('verbose')

        if not ignore_cache and f_val in ft_cache:
            t_val = ft_cache[f_val]
            if verbose:
                print(f"{fn.__name__}('{f_val}') returning {t_val} from cache")
            return t_val

        t_val = fn(f_val, **kwargs)
        ft_cache[f_val] = t_val
        return t_val

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


# TODO unit test this
def is_one2one(df, col1, col2):
    # TODO TODO maybe modify this to print violations of 1:1-ness, when there
    # are any (otherwise write another fn for this)

    # from answer by zipa @ stackoverflow.com/questions/50643386
    # might be a faster way to do this...
    first = df.drop_duplicates([col1, col2]).groupby(col1)[col2].count().max()
    second = df.drop_duplicates([col1, col2]).groupby(col2)[col1].count().max()
    return first + second == 2


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

    print_header = 'InChI with multiple combinations of {}:'.format(other_keys)
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


def convertable(chem_id_type):
    """Returns whether the string input is recognized as a convertable type.

    (valid for the `to_type` and `from_type` arguments to the `convert`
    function)
    """
    if chem_id_type in chem_id_types or chem_id_type in equivalent_col_names:
        return True
    return False


# TODO TODO tqdm in all but single element case (way to make it work in case
# where something is iterated over w/ call to convert at each iteration?)

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
def convert(chem_id, from_type=None, to_type='inchi', verbose=False,
    allow_nan=False, allow_conflicts=True, ignore_cache=False, exclude_cols=[],
    already_normalized=False, drop_na=True, report_missing=True,
    report_conflicts=True, try_non_normalized=True, check_one2one=False,
    keep_originals=False, orig_prefix='orig_'):
    """
    allow_nan (bool): whether to err if any lookups fail
    """

    def conversion_fail_errmsg():
        return 'Conversion from {} to {} failed for'.format(from_type, to_type)

    def conversion_fail_err():
        msg = conversion_fail_errmsg() + ' {}'.format(chem_id)
        if not allow_nan:
            raise ValueError(msg)
        elif report_missing:
            print(msg)

    if check_one2one:
        raise NotImplementedError

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

            if drop_na:
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

            # allow_conflicts and exclude_cols are only relevant for DataFrame
            # case, and therefore do not need to be passed.
            # from_type should always be inferred in DataFrame case as well.
            kwargs = {
                'to_type': to_type,
                'allow_nan': allow_nan,
                'ignore_cache': ignore_cache,
                'drop_na': drop_na,
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
                    assert orig_col not in df.columns
                    orig_cols_added.append(orig_col)
                    df[orig_col] = df.index.copy()

                df.index = convert(df.index, **kwargs)
                converted_an_index = True

            if convertable(df.columns.name):
                if keep_originals:
                    orig_col = orig_prefix + df.columns.name
                    assert orig_col not in df.columns
                    orig_cols_added.append(orig_col)
                    df[orig_col] = df.columns.copy()

                df.columns = convert(df.columns, **kwargs)
                converted_an_index = True

            if not converted_an_index:
                kwargs['allow_nan'] = True
                kwargs['drop_na'] = False
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
                            assert orig_col not in df.columns
                            orig_cols_added.append(orig_col)
                            df[orig_col] = df[c].copy()

                        # TODO TODO test other normalize case + other fn name
                        # lookup cases in c in equivalent_col_names case, since
                        # that caused the failure here
                        norm_fn_name = 'normalize_' + ft
                        if norm_fn_name in globals():
                            if verbose:
                                print('Applying {} to input column'.format(
                                    norm_fn_name))

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
                    err_str = 'Conversion from {} to {} failed for:\n'.format(
                        [c for c in cols], to_type
                    )
                    missing = df[missing].drop_duplicates(
                        subset=cols).dropna(subset=cols, how='all')[cols]

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
                values = attempts[cols_in_order].apply(
                    lambda x: x.dropna().unique(), axis=1)
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
                if drop_na:
                    df.dropna(subset=[to_type], inplace=True)

                if conflicts.any():
                    if not allow_conflicts or report_conflicts:
                        # TODO make sure only unique conflicts are printed
                        print('Conflicting lookup results:')
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

    if from_type in equivalent_col_names:
        from_type = equivalent_col_names[from_type]

    if verbose:
        print('Trying to convert {} from {} to {}'.format(
            chem_id, from_type, to_type
        ))

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
                print('{}({}) -> {}'.format(norm_fn_name, old_chem_id, chem_id))

    # Since sometimes (just normalize_cas, for now) normalize fns can return
    # null.
    if pd.isnull(chem_id):
        return chem_id

    # TODO TODO for relationships that are *guaranteed* to be one-to-one,
    # could also populate the reverse direction in the cache when the other
    # direction is filled (CID <-> inchi? definitely inchikey <-> inchi, right?)
    # ...and do they have to be 1:1?

    if not ignore_cache:
        # TODO should this fail into elif if cached value is None?
        # (be consistent w/ all branches on try_non_normalized)
        if chem_id in cache[from_type][to_type]:
            val = cache[from_type][to_type][chem_id]
            if verbose:
                print('Returning {} from cache'.format(val))
            return val

        elif try_non_normalized and old_chem_id in cache[from_type][to_type]:
            val = cache[from_type][to_type][old_chem_id]
            if verbose:
                # TODO here and in other similar places, replace
                # 'chem_id' w/ value of `from_type`
                print('Falling back to non-normalized chem_id')
                print('Returning {} from cache'.format(val))
            return val

    if not ignore_cache and chem_id in cache[from_type]['cid']:
        cid = cache[from_type]['cid'][chem_id]
        if verbose:
            print('{} of type {} had CID {} in cache'.format(chem_id, from_type,
                cid
            ))
        if cid is None:
            if verbose:
                print('CID in cache was None!')
            conversion_fail_err()
            return None

    elif (not ignore_cache and try_non_normalized and
        old_chem_id in cache[from_type]['cid']):

        cid = cache[from_type]['cid'][old_chem_id]
        if verbose:
            print('Falling back to non-normalized chem_id')
            print('{} of type {} had CID {} in cache'.format(old_chem_id,
                from_type, cid
            ))
        if cid is None:
            if verbose:
                print('CID in cache was None!')
            conversion_fail_err()
            return None
    else:
        f2cid_fn_name = from_type + '2cid'
        if f2cid_fn_name not in globals():
            raise NotImplementedError(('define function {} to support ' +
                'conversion from type {}').format(f2cid_fn_name, from_type)
            )

        if verbose:
            print('Calling CID lookup function', f2cid_fn_name)

        cid = globals()[f2cid_fn_name](chem_id)

        if cid is None:
            if verbose:
                print('Looking up CID for {} of type {} failed!'.format(
                    chem_id, from_type))

            if try_non_normalized:
                if verbose:
                    print('Falling back to non-normalized chem_id')
                cid = globals()[f2cid_fn_name](old_chem_id)

        if cid is None:
            if try_non_normalized and verbose:
                print('Looking up CID for {} of type {} failed!'.format(
                    old_chem_id, from_type))

            for tt in chem_id_types:
                to_type_cache = cache[from_type][tt]
                # To not overwrite hardcoded values for other types.
                if chem_id not in to_type_cache:
                    to_type_cache[chem_id] = None

            conversion_fail_err()
            return None

        if verbose:
            print('CID={}'.format(cid))

        cache[from_type]['cid'][chem_id] = cid

    # TODO TODO should i also support some way of directly going from from_type
    # to to_type, w/o requiring a set of fns about a Compound intermediate?
    # somethings i want that can't be pushed through that first?

    # TODO way to go direct to Compound?
    # (__init__ takes something called a "record dict from the PubChem PUG REST
    # service", which we might already have in the results...) idk...
    compound = pcp.Compound.from_cid(cid)
    if compound is None:
        # TODO should i just assert false here or something?
        if verbose:
            print('Creating Compound from CID {} failed!'.format(cid))

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
        raise NotImplementedError(('define function {} to support ' +
            'conversion to type {}').format(compound2t_fn_name, to_type))

    if verbose:
        print('Calling function {} to get {} from Compound'.format(
            compound2t_fn_name, to_type))

    to_type_val = globals()[compound2t_fn_name](compound)
    cache[from_type][to_type][chem_id] = to_type_val

    if verbose and to_type_val is None:
        print('Conversion of {} from Compound to {} failed!'.format(
            chem_id, to_type))
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
    name = re.sub(r'\([^\(\)]{4,}\)', '', name)

    parts = name.split()
    normed_name = parts[0]
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
            print('Indicator {} was True for {}'.format(
                i.__name__, normed_name))

            normed_name = c(normed_name)

            print('Applying correction {} to get {}.'.format(
                c.__name__, normed_name))

    return normed_name


def normalize_cas(cas):
    # TODO maybe remove first null check in convert if this is the expected
    # behavior of the norm fns? 
    if pd.isnull(cas):
        return cas

    normed_cas = ''.join(cas.replace('"','').split())
    # TODO library seems to return 0-00-0 sometimes... but this is incorrect,
    # right? replace w/ NaN?
    if normed_cas == '0-00-0':
        return None
    return normed_cas


def pubchem_url(cid):
    return 'https://pubchem.ncbi.nlm.nih.gov/compound/{}'.format(cid)


def inchi2pubchem_url(inchi):
    cid = convert(inchi, from_type='inchi', to_type='cid')
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
        warnings.warn('{}\nReturning None.'.format(e))
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
        print('Got multiple results from name={}:'.format(name))
        n_inchi_parts = [len(r.inchi.split('/')) for r in results]
        fewest_inchi_parts = sorted(n_inchi_parts)[0]
        counts = Counter(n_inchi_parts)
        print('Fewest InChI parts: {}'.format(fewest_inchi_parts))
        print('# InChIs w/ that many parts: {}'.format(
            counts[fewest_inchi_parts]))
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
        warnings.warn('{}\nReturning None.'.format(e))
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
        print('Got multiple results from CAS={}'.format(cas))
        n_inchi_parts = [len(r.inchi.split('/')) for r in results]
        fewest_inchi_parts = sorted(n_inchi_parts)[0]
        counts = Counter(n_inchi_parts)
        print('Fewest InChI parts: {}'.format(fewest_inchi_parts))
        print('# InChIs w/ that many parts: {}'.format(
            counts[fewest_inchi_parts]))
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
        warnings.warn('{}\nReturning None.'.format(e))
        return None
    assert len(results) == 1
    return results[0].cid


def inchi2cid(inchi, verbose=False):
    inchi = 'InChI=' + inchi
    try:
        results = pcp.get_compounds(inchi, 'inchi')
    except urllib.error.URLError as e:
        warnings.warn('{}\nReturning None.'.format(e))
        return None
    assert len(results) == 1
    return results[0].cid


def compound2cid(compound):
    return compound.cid


# TODO TODO TODO scrape name at the top of pubchem page if there is no other way
# of getting this. most of the time this would probably be preferable to IUPAC
# name (as long as it doesn't depend on how the page is accessed...).
# e.g. 'linalool' vs IUPAC '3,7-dimethylocta-1,6-dien-3-ol'
def compound2name(compound):
    # TODO this work? 
    return compound.iupac_name


def compound2inchi(compound):
    inchi = compound.inchi
    # TODO sometimes, is it just the prefix?
    assert inchi.startswith('InChI=')
    inchi = inchi[6:]
    # TODO actually check format meets some minimum of the inchi standard
    assert len(inchi) > 0
    return inchi


def compound2cas(compound, verbose=False):
    # TODO TODO how to start w/ a compound and get CAS?
    r = results[0]
    cas_number_candidates = []
    for syn in r.get('Synonym', []):
        match = re.match('(\d{2,7}-\d\d-\d)', syn)
        if match:
            cas_number_candidates.append(match.group(1))
        # TODO so not every entry in pubchem has an associated CAS?

    if len(cas_number_candidates) == 0:
        if verbose:
            print('No CAS numbers found online for {}'.format(cas))
        cas_num = None
    else:
        cas_num = sorted(cas_number_candidates)[0]
    return cas_num


def compound2smiles(compound):
    return compound.canonical_smiles


def fmt_id_type(chem_id_type):
    if chem_id_type in equivalent_col_names:
        chem_id_type = equivalent_col_names[chem_id_type]

    if chem_id_type not in chem_id_types:
        raise KeyError('unrecognized chem_id type. options are: {}'.format(
            chem_id_types))

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


def chemspider_url(csid):
    return f'http://www.chemspider.com/Chemical-Structure.{csid}.html'


def chemspider_experimental_density(csid):
    import requests
    from bs4 import BeautifulSoup
    import time

    url = chemspider_url(csid)
    print(url)
    # TODO implement some requests rate limiting
    response = requests.get(url)
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
        curr_density = float(row.text.split()[0].strip())
        # To filter out what seem to be plainly wrong densities for e.g.
        # http://www.chemspider.com/Chemical-Structure.4444608.html
        # (20 g/mL)
        # Lowered threshold b/c 6.4 g/mL from Biosynth at this source:
        # http://www.chemspider.com/Chemical-Structure.66391.html
        # TODO maybe just filter out data from Merck Millipore?
        if curr_density < 2:
            density = curr_density

    return density


nist_url_prefix = 'http://webbook.nist.gov'
def nist_webbook_url(inchi):
    inchi_str = 'InChI=' + inchi.replace('/','%2F').replace(',','%2C').replace(
        '+', '%2B')
    return f'{nist_url_prefix}/cgi/cbook.cgi?{inchi_str}&Units=SI'


def nist_henrys_law_url(inchi):
    return nist_webbook_url(inchi) + '&cSO=on'


# TODO get temp dependence coeff if available
# TODO assert temp it's defined at is same-ish (or is it always same?)?
@cached
def inchi2k_henry(inchi):
    import requests
    from bs4 import BeautifulSoup
    import time

    url = nist_henrys_law_url(inchi)
    print(url)
    response = requests.get(url)
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
        print(len(soup.find('ol').find_all('a')), 'matches')
        # Should only be one of these <ol> tags in the html.
        for match in soup.find('ol').find_all('a'):
            href = match.attrs['href']
            url = f'{nist_url_prefix}{href}'
            print(url)

            response = requests.get(url)
            # crude rate limit
            time.sleep(0.5)
            html = response.text
            soup = BeautifulSoup(html, features='html.parser')
            inchi_matches = soup.find_all('span', string=inchi_str)
            if len(inchi_matches) == 0:
                continue

            print(len(inchi_matches))
            assert len(inchi_matches) == 1
            inchi_match = inchi_matches[0]
            # TODO test substrings dont match
            print(inchi_str)
            print(inchi_match.text)

            table = soup.find('table',
                attrs={'aria-label': "Henry's Law constant (water solution)"}
            )
            if table is not None:
                print('NEW URL:', url)
                break
            #import ipdb; ipdb.set_trace()

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
        print('MATCH WAS UNAMBIGUOUS')
        # TODO units always mol/kg*bar (probably)?
        table = soup.find('table',
            attrs={'aria-label': "Henry's Law constant (water solution)"}
        )

    if table is None:
        print('NO HENRY LAW DATA\n')
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
                #import ipdb; ipdb.set_trace()
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
    
    print('Kh:', kh)
    print()
    return kh


concentration_dimensionalities = [
    '[length]^-3 [mass]',
    # This seems to be how pint formats molarity (M).
    '[length]^-3 [substance]'
]
preferred_concentration_unit = 'g/L'
dimensionality2preferred_units = {
    '[length]^-3 [mass]': preferred_concentration_unit,
    '[length]^-3 [substance]': 'molar',
    '[length]^-1 [mass] [time]^-2': 'kPa'

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
        'concentration' or something `pint` can parse.

    `test_units` is a `pint` `Quantity` or `Unit` object.

    Equivalent in the sense that they share the same number of factors of each
    base unit type.
    """
    # Lone numbers, for instance, parse to dimensionless units.
    if expected_units is None:
        return not test_units.dimensionless

    elif expected_units == 'concentration':
        return is_concentration_unit(test_units)

    else:
        return test_units.dimensionality == ureg[expected_units].dimensionality


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

    # This happens if p == ':'
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


float_re = r'(\d+(?:,\d{3})*(?:\.\d*)?)'
# TODO also need to match case where degree sign is not there?
# case where F/C not there (maybe C assumed?)?
# TODO maybe check there were not multiple matches to this regex
temp_re = r'(?:@|at)\s?' + float_re + r'\s?°\s?([CFK])'
def parse_temperature_c(string):
    """Returns None or scalar temparture in celsius and part of string not used.
    """
    match = re.search(temp_re, string)
    if match is None:
        return None, string
    else:
        groups = match.groups()
        temp_unit_str = groups[-1]
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
def parse_pubchemprops_string(string, expected_units=None, target_temp_c=None):
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
                expected_units=expected_units, target_temp_c=target_temp_c
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
    range_re = r'(?:^|\s)' + float_re + r'\s?(?:-|to)\s?' + float_re
    match = re.search(range_re, string)
    if match is not None:
        parsed_dict['from_range'] = True

        groups = match.groups()
        range_min = float(groups[0])
        range_max = float(groups[1])
        assert range_min < range_max, 'range not in expected order'
        quantity = range_min + (range_max - range_min) / 2

        # Assuming we won't need to check for inequalities in this range case.
        parts_before_quantity = []

        # So the components of the range are not also matched in loop over
        # parts.
        string = re.sub(range_re, '', string)

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
                curr_quantity = float(groups[0].replace(',','') +
                    'e' + groups[-1]
                )
            else:
                remaining_parts.append(p)
                continue

        if quantity is None:
            parts_before_quantity = list(remaining_parts)
            quantity = curr_quantity
        else:
            raise ValueError('multiple parseable quantities')

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

        # TODO TODO check dimensions of units match one of expected things
        # (or are at least not dimensionless, if no expected units provided)
        # TODO maybe store cache in chemutils so units can be checked consistent
        # against other lookups of same property, in case expected units not 
        # specified

        if units is None:
            units = curr_units
        else:
            raise ValueError('multiple unit definitions')

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


def pubchemprops_lookup(inchi, prop_name, expected_units=None,
    target_temp_c=None):
    """
    Requires my fork of pubchemprops which adds the `cid` kwarg to
    `get_second_layer_props`.
    """
    from pubchemprops.pubchemprops import get_second_layer_props

    # TODO also pass ignore_cache added by `cached` through so it's usable here?
    # TODO TODO also handle case where this is null
    cid = convert(inchi, from_type='inchi', to_type='cid')
    # TODO delete after handling null case
    assert type(cid) is int
    #

    # TODO delete
    print(f'LOOKING UP {prop_name} OF {inchi}')
    print(pubchem_url(cid))
    #

    # TODO probably only do these lookups if expected_units == 'concentration'
    compound = pcp.Compound.from_cid(cid)
    # TODO TODO cache these in my cache? i assume the from_cid + this may
    # make at least 1-2 requests? from inchi or from CID? probably inchi?
    molecular_weight = ureg.Quantity(compound.molecular_weight, 'g/mol')
    #

    # TODO where will i have to check for / return None in case of failed
    # lookup?
    ret = get_second_layer_props(cid, [prop_name], cid=True)

    parsed_dicts = []
    # TODO also aggregate and (optionally?) return sources of the data
    if len(ret) > 0:
        ret = ret[prop_name]
        # TODO TODO maybe take closest to room temp / average across all
        for result in ret:
            value_dict = result['Value']

            if 'Number' in value_dict:
                # TODO use result['Reference'] for source in this case

                # I've only seen it have 'Unit' and 'Number' keys.
                assert len(value_dict) <= 2

                # TODO maybe just share from parsing fn above?
                keys = ['quantity', 'units', 'temperature_c', 'source',
                    'from_range'
                ]
                parsed_dict = {k: None for k in keys}

                # TODO TODO are there cases where i'd also want to share all of
                # the unit parsing logic above here?
                unit_str = value_dict['Unit']

                temp_c, unit_str = parse_temperature_c(unit_str)
                # `temp_c` can be `None`, but that will just overwrite w/
                # existing value.
                parsed_dict['temperature_c'] = temp_c

                unit_str = unit_str.replace('()','').strip()
                units = ureg[unit_str]
                assert in_expected_units(expected_units, units)
                parsed_dict['units'] = units.u

                nums = value_dict['Number']
                assert len(nums) == 1
                parsed_dict['quantity'] = nums[0]

                parsed_dict['from_range'] = False

                parsed_dicts.append(parsed_dict)

            elif 'StringWithMarkup' in value_dict:
                assert len(value_dict) == 1
                pstr = value_dict['StringWithMarkup'][0]['String']

                parsed_dict = parse_pubchemprops_string(pstr,
                    expected_units=expected_units,
                    target_temp_c=target_temp_c
                )
                if parsed_dict['quantity'] is None:
                    continue
                else:
                    parsed_dicts.append(parsed_dict)

                assert len(value_dict['StringWithMarkup']) == 1
                # If there is an extra Key under the first element, it seems
                # to have always just been a reference to water...

            else:
                raise ValueError('unexpected value dict')

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
        for dims, preferred_unit in dimensionality2preferred_units.items():
            if uq.check(dims):
                found_dims = True

                parsed_dict['quantity'] = uq.to(preferred_unit).magnitude
                # Changing type to a string here.
                # Consider keeping as a pint obj (though that may be tricky,
                # which is why I opted to do it this way).
                parsed_dict['units'] = preferred_unit

                # Assuming check will only be True once (as it should).
                break

        assert found_dims, 'no matching dimensionality with preferred units'

    for parsed_dict in parsed_dicts:
        m = parsed_dict['quantity']
        u = parsed_dict['units']
        uq = ureg.Quantity(m, u)
        # It must be a molar concentration then (assuming it's just between that
        # option and mass/vol)
        if '[substance]' in uq.dimensionality.keys():
            parsed_dict['quantity'] = (uq * molecular_weight
                ).to(preferred_concentration_unit).magnitude
            parsed_dict['units'] = preferred_concentration_unit

    if len(parsed_dicts) == 0:
        # TODO or do i want to return null_dict? as long as this works w/ pandas
        # fns...
        return None

    # TODO maybe considate w/ below...
    elif len(parsed_dicts) == 1:
        return parsed_dicts[0]['quantity']

    elif len(parsed_dicts) > 1:
        # TODO TODO outlier / within-order-of-magnitude detection here
        # TODO TODO + averaging / selection based on temperatures

        # TODO TODO maybe also prefer stuff that has a temperature to stuff that
        # doesn't... (assuming it's in a reasonable range about temp we want)

        quantities = [d['quantity'] for d in parsed_dicts]
        return np.mean(quantities)


# TODO hardcode these values:
# "Soluble in about 720 parts water, in many organic solvents" (l117 in test
# data)
# "0.43% (by wt) in water" (l71)
# "In water, 6,700 ppm at 25 °C" (l24)
# "Soluble in 120 parts water at 25 °C" (l25)
# maybe this and all the other NIOSH stuff
# "5 % (NIOSH, 2016)" (l39)
@cached
def inchi2water_solubility(inchi):
    # TODO TODO maybe convert to str repr before caching so that i don't
    # have to worry about that pint app registry thing?
    # (see notes in docs about comparing stuff across ureg instances)
    # (or just test stuff as i'm doing it now... ?)

    # TODO also pass ignore_cache added by `cached` through so it's usable here?
    return pubchemprops_lookup(inchi, 'Solubility',
        expected_units='concentration'
    )
    # TODO TODO is this something that can be retrieved via chemspipy??
    # fallback to that if so.
    # TODO TODO NIST?


@cached
def inchi2vapor_pressure(inchi):
    return pubchemprops_lookup(inchi, 'Vapor Pressure', expected_units='kPa')
    # TODO TODO is this something that can be retrieved via chemspipy??
    # fallback to that if so.
    # TODO TODO NIST?


@cached
def inchi2density(inchi, chemspipy_api_key=None):
    # TODO in the future, only warn if any are missing
    from chemspipy import ChemSpider

    density = pubchemprops_lookup(inchi, 'Density')
    if density is not None:
        return density

    # TODO maybe cache these separately if we can't tell where None came from?
    if chemspipy_api_key is None:
        return None

    cs = ChemSpider(chemspipy_api_key)
    # Docs seem to claim inchi could be passed just as if it were name,
    # but my testing did not seem to bear that out.
    # TODO probably also cache these (if using for anything other than density)
    results = cs.search('InChI=' + inchi, domain='inchi')
    # TODO handle lookup failure here / above
    if len(results) == 1:
        csid = results[0].csid
        # Since I couldn't find a way to get arbitrary properties like
        # "Experimental Density" through the API.
        density = chemspider_experimental_density(csid)
        print('CHEMSPIDER DENSITY:', density)
    #else:
    #    print(f'{len(results)} CHEMSPIDER RESULTS')

    #print()
    return density


_odor_inventory_gsheet = None
def odor_inventory_gsheet(use_cache=False, verbose=False):
    '''Returns a DataFrame with data odor inventory data from Google sheet.
    '''
    global _odor_inventory_gsheet
    if _odor_inventory_gsheet is not None:
        return _odor_inventory_gsheet

    gsheet_cache_file = '.odor_inventory_gsheet_cache.p'
    if use_cache and exists(gsheet_cache_file):
        print('Loading odor inventory sheet data from cache at {}'.format(
            gsheet_cache_file))
        # TODO use pandas interface (if not factoring out whole gsheet reading
        # thing)
        with open(gsheet_cache_file, 'rb') as f:
            df = pickle.load(f)
    else:
        pkg_data_dir = split(__file__)[0]
        with open(
            join(pkg_data_dir, 'odor_inventory_sheet_link.txt'), 'r') as f:

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

    df = convert(df, to_type='inchi', allow_nan=False,
        verbose=verbose, report_conflicts=False
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
        print('could not find inchi for odor {}!'.format(odor_name))
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
                    warnings.warn('no abbreviation in odor inventory '
                        + 'gsheet for odor "{}"!'.format(o))
                    # consistent case be damned here.
                    abbrev = chr(ord('A') + i)
                    i += 1

                o2a[o] = abbrev
    return o2a

