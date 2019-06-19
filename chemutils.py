
import re
import os
import pickle
import atexit
import urllib.error
import warnings

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


chem_id_types = [
    'name',
    'cas',
    'inchi',
    'inchikey',
    'smiles',
    'cid'
]
# No matter the type these are to be converted to, they will always return a
# null value. (maybe handle some other way?)
manual_type2null_keys = {
    'name': {
        'spontaneous firing rate',
        'odor'
    }
}
# TODO TODO TODO just have manual keys override values saved in cache when cache
# is loaded, to never have to explicitly refer to the separate hardcoded stuff
manual = {
    'name': {
        'cas': {
            # TODO TODO TODO deal w/ linalool (or maybe at cid level?)
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


def init_cache():
    global cache
    cache = {ft: {tt: dict() for tt in chem_id_types} for ft in chem_id_types}


def add_manual_overrides():
    global cache

    for ft, null_keys in manual_type2null_keys.items():
        for nk in null_keys:
            for tt in chem_id_types:
                cache[ft][tt][nk] = None

    for ft, tts in manual.items():
        for tt in tts:
            cache[ft][tt].update(tts[tt])


def delete_cache():
    """Deletes and clears the cache at cache_file
    """
    global cache
    init_cache()
    if os.path.exists(cache_file):
        os.remove(cache_file)

#delete_oldcache()
#delete_cache()
if os.path.exists(cache_file):
    try:
        with open(cache_file, 'rb') as f:
            cache = pickle.load(f)
        # TODO TODO overwrite w/ manual stuff here?
    except ValueError:
        print('Cache was in unreadable format. Deleting.')
        delete_cache()
else:
    init_cache()
add_manual_overrides()


def save_oldcache():
    with open(oldcache_file, 'wb') as f:
        pickle.dump((to_cas_cache, to_name_cache, to_inchi_cache,
            to_smiles_cache), f)


def save_cache():
    with open(cache_file, 'wb') as f:
        pickle.dump(cache, f)


atexit.register(save_oldcache)
atexit.register(save_cache)


# TODO why did i have this? delete...
def clear_inchi_cache():
    """Just clears the InChI cache.
    """
    global to_inchi_cache
    to_inchi_cache = dict()
    save_oldcache()
#


# TODO TODO TODO some fn to strip chirality information from inchi for further
# normalization
# 1) possible?
# 2) always what we want (i.e. does the biology only make one sometimes / can we
# actually separate some on our column?)?
# 3) usually what we want?


def convert(chem_id, from_type=None, to_type='inchi', verbose=False,
    allow_nan=False, ignore_cache=False):

    equivalent_col_names = {
        'odor': 'name',
        'cas_number': 'cas'
    }
    def conversion_fail_errmsg():
        return 'Conversion from {} to {} failed for'.format(from_type, to_type)

    def conversion_fail_err():
        msg = conversion_fail_errmsg() + ' {}'.format(chem_id)
        raise ValueError(msg)

    # TODO test w/ index/series/df (w/ index / columns of matching name)
    if hasattr(chem_id, 'shape') and len(chem_id.shape) > 0:
        if to_type in equivalent_col_names:
            to_name = to_type
            to_type = equivalent_col_names[to_type]
        else:
            to_name = to_type

        # This covers both Index and Series.
        if len(chem_id.shape) == 1:
            if from_type is None:
                if chem_id.name is None:
                    raise ValueError('either pass from_type or name axes')

                from_type = chem_id.name
                if from_type in equivalent_col_names:
                    from_type = equivalent_col_names[from_type]

            # Mostly for recursion from DataFrame case.
            if from_type not in chem_id_types:
                if verbose:
                    print('from_type {} not recognized. not converting.'.format(
                        from_type))
                return chem_id

            # Forcing allow_nan to True so we can report each of the failing
            # lookups.
            fn = lambda x: convert(x, from_type=from_type, to_type=to_type,
                verbose=verbose, ignore_cache=ignore_cache, allow_nan=True)

            # TODO test each of these cases
            # TODO also test w/ passing / failing allow_nan case for each
            # This means it was a Series, not an Index.
            # TODO TODO which case do i actually want to support:
            # series w/ values to convert in index, or w/ them as the value of
            # the series???? either?
            if hasattr(chem_id, 'index'):
                raise NotImplementedError
                converted = chem_id.rename(fn).rename(to_name)
                to_check = chem_id.index
            else:
                converted = chem_id.map(fn)
                converted.name = to_type
                to_check = converted

            missing = to_check.isnull()
            if not allow_nan and missing.any():
                err_str = conversion_fail_errmsg() + ':\n'
                err_str = 'Conversion from {} to {} failed for:\n'.format(
                    from_type, to_type)

                for m in chem_id[missing]:
                    err_str += str(m) + '\n'

                raise ValueError(err_str)

            return converted

        # This covers DataFrames
        elif len(chem_id.shape) == 2:
            # TODO also support checking names of single rows/columns and
            # converting those (adding a column)?
            if from_type is not None:
                raise NotImplementedError('only conversion of named indices' +
                    'supported for DataFrames')

            df = chem_id.copy()
            # TODO maybe pass all kwargs at once some way
            df.index = convert(df.index, to_type=to_type, verbose=verbose,
                allow_nan=allow_nan, ignore_cache=ignore_cache)
            df.columns = convert(df.columns, to_type=to_type, verbose=verbose,
                allow_nan=allow_nan, ignore_cache=ignore_cache)
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

    if pd.isnull(chem_id):
        return chem_id

    # TODO unit test each of these cases, including None handling

    norm_fn_name = 'normalize_' + from_type
    if norm_fn_name in globals():
        old_chem_id = chem_id
        chem_id = globals()[norm_fn_name](chem_id)
        if verbose and old_chem_id != chem_id:
            print('{}({}) -> {}'.format(norm_fn_name, old_chem_id, chem_id))

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
            if not allow_nan:
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

            if not allow_nan:
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

        if not allow_nan:
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

        if not allow_nan:
            conversion_fail_err()

    return to_type_val


def normalize_name(name):
    parts = name.split()
    normed_name = parts[0]
    for a, b in zip(parts, parts[1:]):
        # We don't want to join adjacent words.
        if a[-1].isalpha() and b[0].isalpha():
            b = ' ' + b
        normed_name += b

    corrections = {
        'linalool oxide': 'trans-linalool oxide',
        '4-ethyl guaiacol': '4-ethylguaiacol'
    }
    if normed_name in corrections:
        return corrections[normed_name]

    prefix_corrections = {
        'a-': 'alpha-',
        'b-': 'beta-',
        'g-': 'gamma-',
        'E2-': '(E)-2-',
        'E3-': '(E)-3-',
        'Z2-': '(Z)-2-',
        'Z3-': '(Z)-3-',
    }
    for p, c in prefix_corrections.items():
        if normed_name.startswith(p):
            normed_name = c + normed_name[len(p):]

    return normed_name


# TODO maybe also take fns <type>2results or something? which always gets
# CID as here? probably not worth it...
def name2cid(name, verbose=False):
    try:
        # TODO TODO should this have been get_compounds, as in cas2name above?
        # what's the difference?
        results = pcp.get_synonyms(name, 'name')
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

    # TODO TODO if results len > 1, maybe err
    return results[0]['CID']


def cas2cid(cas, verbose=False):
    try:
        # TODO TODO why get_compounds rather than get_synonyms here? what's the
        # difference?
        # TODO TODO TODO if no difference, call name2cid here, as otherwise they
        # are the same
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

    # TODO TODO if results len > 1, maybe err

    # TODO why did this seem to work in one other case?
    # is this one of / the only synonyms / compounds difference?
    #return results[0]['CID']
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

