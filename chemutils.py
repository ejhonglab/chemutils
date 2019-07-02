
import re
import os
import pickle
import atexit
import urllib.error
import warnings
from collections import Counter

import pubchempy as pcp
import pandas as pd


# Procedure for entering CAS numbers in here (follow this to maintain
# consistency of name2cas fn):
# 1) Find pubchem page for compound of interest.
# 2) Get one of the names under "Chemical names"
# 3) Run name2cas on that chemical name.
# 4) Use the output cas number here.
# TODO TODO TODO deal w/ linalool
manual_name2cas = {
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
    'E2-hexenyl acetate': '10094-40-3',
    'spontaneous firing rate': None,
    'odor': None
}
# TODO TODO delete all old-style cache stuff

oldcache_file = os.path.expanduser('~/.chemutils_oldcache.p')
cache_file = os.path.expanduser('~/.chemutils_cache.p')

def delete_oldcache():
    """Deletes and clears the cache at oldcache_file
    """
    global to_cas_cache
    global to_name_cache
    global to_inchi_cache
    global to_smiles_cache

    to_cas_cache = dict()
    to_name_cache = dict()
    to_inchi_cache = dict()
    to_smiles_cache = dict()

    if os.path.exists(oldcache_file):
        os.remove(oldcache_file)


if os.path.exists(oldcache_file):
    try:
        with open(oldcache_file, 'rb') as f:
            to_cas_cache, to_name_cache, to_inchi_cache, to_smiles_cache = \
                pickle.load(f)
    except ValueError:
        print('Old-style cache was in unreadable format. Deleting.')
        delete_oldcache()
else:
    to_cas_cache = dict()
    to_name_cache = dict()
    to_inchi_cache = dict()
    to_smiles_cache = dict()


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

# TODO TODO TODO flag to ignore manual overrides, to the extent that they are
# causing problems. maybe phase them out altogether?

# TODO TODO try deleting / ignoring this hardcoded stuff and see if it still
# works. my changes to normalize_name probably fixed a few of these cases.
# (and if it can be ignored, delete it)
manual = {
    'name': {
        'cas': {
            # TODO TODO TODO deal w/ linalool (or maybe at cid level?)
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
        }
    }
}


# TODO TODO move inventory loading here? (maybe w/o specific link?)

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
    _cache = {ft: {tt: dict() for tt in chem_id_types} for ft in chem_id_types}
    add_manual_overrides(_cache)
    return _cache


def clear_cache():
    global cache
    cache = init_cache()


def delete_cache():
    """Deletes and clears the cache at cache_file
    """
    clear_cache()
    if os.path.exists(cache_file):
        os.remove(cache_file)


def load_cache():
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'rb') as f:
                _cache = pickle.load(f)
            return _cache
        except ValueError:
            print('Cache was in unreadable format. Deleting.')
            delete_cache()
            # Will have already been init'd in delete_cache call.
            return cache
    else:
        return init_cache()


def save_oldcache():
    with open(oldcache_file, 'wb') as f:
        pickle.dump((to_cas_cache, to_name_cache, to_inchi_cache,
            to_smiles_cache), f)


def save_cache():
    """
    Load old cache first, so that clear_cache doesn't have to worry about
    screwing up on-disk cache when atexit call to save_cache triggers.
    """
    # TODO need "global cache", as long as cache def is below fn def, or what?
    # if that's all, move "cache = load_cache()" just above this, and declare
    # as global in load_cache, cause circularity
    old_cache = load_cache()
    # So any current values overwrite old_cache values.
    old_cache.update(cache)
    
    with open(cache_file, 'wb') as f:
        pickle.dump(old_cache, f)


# TODO test
cache = load_cache()
atexit.register(save_oldcache)
atexit.register(save_cache)


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


def basic_inchi(inchi, include_h=True):
    parts = inchi.split('/')
    # TODO do i want h? can this information not be mostly / entirely inferred?
    # TODO actually, when is h layer NOT included in my data?
    if include_h:
        keep = {'c','h','b'}
    else:
        keep = {'c','b'}
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


# TODO better name?
# TODO maybe just delete the include_no_h option and code for that case
def inchi_diff_in_details(df, include_no_h=False):
    assert 'inchi' in df.columns
    other_keys = [x for x in df.columns if x in chem_id_types]
    other_keys += [x for x in df.columns if x in equivalent_col_names]
    cols_to_show = other_keys + ['inchi']
    # TODO also count name/cas/(name,cas) here, as w/ inchi_counts?
    ncdf = df.drop_duplicates(subset=cols_to_show)
    print('InChI with multiple combinations of {}:'.format(other_keys))
    for gn, gdf in ncdf.groupby('inchi'):
        if len(gdf) <= 1:
            continue
        print_full_df(gdf[cols_to_show])
        print('')
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
    if include_no_h:
        df['basic_inchi_no_h'] = df.inchi.apply(
            lambda x: basic_inchi(x, include_h=False))

    cols_to_show += ['inchi_counts', 'basic_inchi']

    # TODO maybe just print part after common prefix for each of these?
    # (part after basic inchi / basic inchi w/o h)

    print('Multiple standard InChI that map to InChI w/o chirality info:')
    for gn, gdf in df.groupby('basic_inchi'):
        if len(gdf) <= 1:
            continue
        print_full_df(gdf[cols_to_show])
        for c in convert(gdf.inchi, to_type='cid'):
            print_pubchem_link(c)
        print('')
    print('')

    if include_no_h:
        # TODO these two sections might be redundant. check + delete one if so.
        print('Multiple std InChI that map to InChI w/o chirality OR hydrogen' +
            ' info:')
        cols_to_show = cols_to_show + ['basic_inchi_no_h']
        for gn, gdf in df.groupby('basic_inchi_no_h'):
            if len(gdf) <= 1:
                continue
            print_full_df(gdf[cols_to_show])
            for c in convert(gdf.inchi, to_type='cid'):
                print_pubchem_link(c)
            print('')
        print('')

        hdf = df.drop_duplicates(subset=['basic_inchi','basic_inchi_no_h']
            ).copy()
        print('Multiple no-chirality InChI that map to InChI w/o hydrogen ' +
            'info:')
        for gn, gdf in hdf.groupby('basic_inchi_no_h'):
            if len(gdf) <= 1:
                continue
            print_full_df(gdf[cols_to_show])
            for c in convert(gdf.inchi, to_type='cid'):
                print_pubchem_link(c)
            print('')
        print('')


# TODO why does "tetrahedral stereochemistry of atoms and allenes" layer
# show up more than "double bonds..."?
# TODO maybe draw structures with / without this info (if possible) to see what
# information it contains?
def count_inchi_layers(inchis):
    counts = Counter(inchis.apply(inchi_layers).agg('sum'))
    return counts


# TODO TODO TODO some fn to strip chirality information from inchi for further
# normalization
# 1) possible?
# 2) always what we want (i.e. does the biology only make one sometimes / can we
# actually separate some on our column?)?
# 3) usually what we want?


# TODO rename to can_convert? is_chemid_type?
def convertable(chem_id_type):
    if chem_id_type in chem_id_types or chem_id_type in equivalent_col_names:
        return True
    return False


# TODO in ignore_cache case, still make a interpreter run / call specific
# cache, or like de-dupe and re-dupe, so as to still test new behavior, but also
# not waste time. (not as much of a priority if clear_cache approach works)
def convert(chem_id, from_type=None, to_type='inchi', verbose=False,
    allow_nan=False, allow_conflicts=True, ignore_cache=False, exclude_cols=[],
    already_normalized=False, drop_na=True, report_missing=True):

    def conversion_fail_errmsg():
        return 'Conversion from {} to {} failed for'.format(from_type, to_type)

    def conversion_fail_err():
        msg = conversion_fail_errmsg() + ' {}'.format(chem_id)
        if not allow_nan:
            raise ValueError(msg)
        elif report_missing:
            print(msg)

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
                already_normalized=already_normalized, report_missing=False)

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
                # TODO check this unique / dropna is what i want
                for m in chem_id[missing].dropna().unique():
                    err_str += str(m) + '\n'

                if not allow_nan:
                    raise ValueError(err_str)
                elif report_missing:
                    print(err_str)

            if drop_na:
                converted = converted.dropna()

            return converted

        # This covers DataFrames
        elif len(chem_id.shape) == 2:
            # TODO also support checking names of single rows/columns and
            # converting those (adding a column)?
            if from_type is not None:
                raise NotImplementedError('only conversion of named indices' +
                    'supported for DataFrames')

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
                converted_an_index = True
                df.index = convert(df.index, **kwargs)

            if convertable(df.columns.name):
                converted_an_index = True
                df.columns = convert(df.columns, **kwargs)

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
                    cols = attempts.columns
                    err_str = \
                        'Conversion from {} to {} failed for:\n'.format(
                            [c for c in cols], to_type)

                    with pd.option_context('display.max_colwidth', -1):
                        err_str += df[missing].drop_duplicates(
                            subset=cols).dropna(subset=cols, how='all'
                            ).to_string()

                    if not allow_nan:
                        raise ValueError(err_str)
                    elif report_missing:
                        print(err_str)

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
                    lambda x: None if len(x) == 0 else x[0])
                if drop_na:
                    df.dropna(subset=[to_type], inplace=True)

                if conflicts.any():
                    # Assuming we always want to know if there are conflicts,
                    # verbose or not.
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
        raise ValueError('specify from_type if not passing pandas object w/ ' +
            'named axes')

    if from_type in equivalent_col_names:
        from_type = equivalent_col_names[from_type]

    if verbose:
        print('Trying to convert {} from {} to {}'.format(
            chem_id, from_type, to_type))

    # To short-circuit normalization in the more-common case where null exists
    # before normalization.
    if pd.isnull(chem_id):
        return chem_id

    # TODO unit test each of these cases, including None handling

    if not already_normalized:
        norm_fn_name = 'normalize_' + from_type
        if norm_fn_name in globals():
            old_chem_id = chem_id
            chem_id = globals()[norm_fn_name](chem_id)
            if verbose and old_chem_id != chem_id:
                print('{}({}) -> {}'.format(norm_fn_name, old_chem_id, chem_id))

    # Since sometimes (just normalize_cas, for now) normalize fns can return
    # null.
    if pd.isnull(chem_id):
        return chem_id

    # TODO TODO for relationships that are *guaranteed* to be one-to-one,
    # could also populate the reserve direction in the cache when the other
    # direction is filled (CID <-> inchi? definitely inchikey <-> inchi, right?)

    if not ignore_cache:
        if chem_id in cache[from_type][to_type]:
            val = cache[from_type][to_type][chem_id]
            if verbose:
                print('Returning {} from cache'.format(val))
            return val

    if not ignore_cache and chem_id in cache[from_type]['cid']:
        cid = cache[from_type]['cid'][chem_id]
        if verbose:
            print('{} of type {} had CID {} in cache'.format(chem_id, from_type,
                cid))
        if cid is None:
            if verbose:
                print('CID in cache was None!')
            conversion_fail_err()
            return None
    else:
        f2cid_fn_name = from_type + '2cid'
        if f2cid_fn_name not in globals():
            raise NotImplementedError(('define function {} to support ' +
                'conversion from type {}').format(f2cid_fn_name, from_type))

        if verbose:
            print('Calling CID lookup function', f2cid_fn_name)

        cid = globals()[f2cid_fn_name](chem_id)
        if cid is None:
            if verbose:
                print('Looking up CID for {} of type {} failed!'.format(
                    chem_id, from_type))

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

    parts = name.split()
    normed_name = parts[0]
    for a, b in zip(parts, parts[1:]):
        # We don't want to join adjacent words.
        if a[-1].isalpha() and b[0].isalpha():
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
        'e3-hexenol': '(E)-hex-3-en-1-ol'
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


def print_pubchem_link(cid):
    print('https://pubchem.ncbi.nlm.nih.gov/compound/{}'.format(cid))


# TODO maybe also take fns <type>2results or something? which always gets
# CID as here? probably not worth it...
def name2cid(name, verbose=False):
    try:
        # TODO TODO should this have been get_compounds, as in cas2name above?
        # what's the difference?
        #results = pcp.get_synonyms(name, 'name')
        results = pcp.get_compounds(name, 'name')
    except urllib.error.URLError as e:
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
            print_pubchem_link(r.cid)
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
            print_pubchem_link(r.cid)
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


#def compound2inchikey(compound):


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


def name2cas(name, verbose=False):
    """Returns the CAS number for the chemical with the given name.

    If there seem to be multiple CAS numbers for this chemical, the first
    (as sorted by Python string comparison) is returned, so the answer should at
    least be consistent.
    """
    if pd.isnull(name):
        return name

    if name in manual_name2cas:
        return manual_name2cas[name]

    if name in to_cas_cache:
        if verbose:
            print('{} in to_cas_cache'.format(name))
        return to_cas_cache[name]

    try:
        results = pcp.get_synonyms(name, 'name')
    except urllib.error.URLError as e:
        warnings.warn('{}\nReturning None.'.format(e))
        return None

    if len(results) == 0:
        cas_num = None
        to_cas_cache[name] = cas_num
        return cas_num

    # TODO TODO if results if len > 1, maybe err / warn
    r = results[0]
    cas_number_candidates = []
    for syn in r.get('Synonym', []):
        match = re.match('(\d{2,7}-\d\d-\d)', syn)
        if match:
            cas_number_candidates.append(match.group(1))
        # TODO so not every entry in pubchem has an associated CAS?

    if len(cas_number_candidates) == 0:
        if verbose:
            print('No CAS numbers for {}'.format(name))
        cas_num = None

    else:
        # TODO maybe taking the lowest number CAS makes the most sense?
        # (moreso than just sorting at least?) (key on first part? matter?)
        cas_num = sorted(cas_number_candidates)[0]

    # TODO delete
    name_from_inverse = cas2name(cas_num)
    if pd.isnull(name_from_inverse):
        print('original name:', name)
        print('CAS from lookup:', cas_num)
        print('inverting cas number back to name failed!')
        import ipdb; ipdb.set_trace()
    '''
    elif name_from_inverse != name:
        print('inverse name did not match original name!')
        print('original name:', name)
        print('name from inverse:', name_from_inverse)
    '''
    #

    to_cas_cache[name] = cas_num
    return cas_num


def cas2cas(cas, verbose=False):
    """
    """
    if pd.isnull(cas):
        return cas

    if cas in to_cas_cache:
        if verbose:
            print('{} in to_cas_cache'.format(cas))
        return to_cas_cache[cas]

    try:
        results = pcp.get_synonyms(cas, 'name')
    except urllib.error.URLError as e:
        warnings.warn('{}\nReturning None.'.format(e))
        return None

    if len(results) == 0:
        cas_num = None
        to_cas_cache[cas] = cas_num
        return cas_num

    # TODO TODO if results if len > 1, maybe err
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

    to_cas_cache[cas] = cas_num
    return cas_num


def cas2name(cas, verbose=False):
    """
    """
    if pd.isnull(cas):
        return cas

    if cas in to_name_cache:
        if verbose:
            print('{} in to_name_cache'.format(cas))
        return to_name_cache[cas]

    try:
        # TODO TODO why get_compounds rather than get_synonyms here? what's the
        # difference?
        results = pcp.get_compounds(cas, 'name')
    except urllib.error.URLError as e:
        warnings.warn('{}\nReturning None.'.format(e))
        return None

    if len(results) == 0:
        name = None
        to_name_cache[cas] = name
        return name

    # TODO TODO if results len > 1, maybe err
    r = results[0]
    # TODO was this a pcp.Compound?
    name = r.iupac_name

    to_name_cache[cas] = name
    return name


def name2compound(name, verbose=False):
    cid = name2cid(name, verbose=verbose)
    compound = pcp.Compound.from_cid(cid)
    return compound


def name2inchi(name, verbose=True):
    """
    """
    if pd.isnull(name):
        return name

    if name in to_inchi_cache:
        if verbose:
            print('{} in to_inchi_cache'.format(name))
        return to_inchi_cache[name]

    compound = name2compound(name)

    if compound is None:
        inchi = None
        to_inchi_cache[name] = inchi
        return inchi

    inchi = compound.inchi
    # TODO sometimes, is it just the prefix?
    assert inchi.startswith('InChI=')
    inchi = inchi[6:]
    # TODO actually check format meets some minimum of the inchi standard
    assert len(inchi) > 0

    to_inchi_cache[name] = inchi
    return inchi


def name2smiles(name, verbose=False):
    if pd.isnull(name):
        return name

    if name in to_smiles_cache:
        if verbose:
            print('{} in to_smiles_cache'.format(name))
        return to_smiles_cache[name]

    compound = name2compound(name)

    if compound is None:
        smiles = None
        to_smiles_cache[name] = smiles
        return smiles

    # TODO diff between this and isomeric_smiles ?
    smiles = compound.canonical_smiles

    to_smiles_cache[name] = smiles
    return smiles


# TODO TODO cache (provide consist interface for this...)
def inchikey2inchi(inchikey):
    """
    """
    if pd.isnull(inchikey):
        return inchikey

    if len(inchikey) != 27:
        # TODO or is it 14? wikipedia seems to have conflicting information
        # 14 from hash of connectivity + hyphen + 8 from hash of other layers
        # + single character indicating kind of inchikey + single char
        # indicating inchi version + hyphen + char for protonation?
        raise ValueError('InChI Keys are 27 characters long')

    results = pcp.get_compounds(inchikey, 'inchikey')
    assert len(results) == 1

    result = results[0]
    raise NotImplementedError
    # TODO TODO fix error (not defined)
    assert r.inchi.startswith('InChI=')
    inchi = result[6:]

    return inchi


if __name__ == '__main__':
    import pandas as pd

    # test
    print(name2cas('Aspirin'))

    # would need to download this file to run this example
    data = pd.read_csv('Odor_inventory _Sheet1.csv')
    for name in data['Chemical']:
        cas = name2cas(name)
        print(name, cas)

    # TODO seems rdkit wants the 'InChI=' prefix
    inchi = '1S/C3H8O/c1-2-3-4/h4H,2-3H2,1H3'

