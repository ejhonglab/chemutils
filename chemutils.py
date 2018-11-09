
import re
import os
import pickle
import atexit

import pubchempy as pcp


# Procedure for entering CAS numbers in here (follow this to maintain
# consistency of name2cas fn):
# 1) Find pubchem page for compound of interest.
# 2) Get one of the names under "Chemical names"
# 3) Run name2cas on that chemical name.
# 4) Use the output cas number here.
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


# TODO allow configuration s.t. to_cas_cache is disabled (env var?)?
cache_file = os.path.expanduser('~/.chemutils_cache.p')
if os.path.exists(cache_file):
    with open(cache_file, 'rb') as f:
        to_cas_cache, to_name_cache = pickle.load(f)
else:
    to_cas_cache = dict()
    to_name_cache = dict()

def save_cache():
    with open(cache_file, 'wb') as f:
        pickle.dump((to_cas_cache, to_name_cache), f)

atexit.register(save_cache)


def delete_to_cas_cache():
    """Deletes the to_cas_cache pickle at ~/.chemutils_to_cas_cache.p
    """
    if os.path.exists(to_cas_cache_file):
        os.remove(to_cas_cache_file)


def name2cas(name, verbose=False):
    """Returns the CAS number for the chemical with the given name.

    If there seem to be multiple CAS numbers for this chemical, the first
    (as sorted by Python string comparison) is returned, so the answer should at
    least be consistent.
    """
    if name in manual_name2cas:
        return manual_name2cas[name]

    if name in to_cas_cache:
        if verbose:
            print('{} in to_cas_cache'.format(name))
        return to_cas_cache[name]

    results = pcp.get_synonyms(name, 'name')

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
        cas_num = sorted(cas_number_candidates)[0]

    to_cas_cache[name] = cas_num
    return cas_num


def cas2cas(cas, verbose=False):
    """
    """
    if cas in to_cas_cache:
        if verbose:
            print('{} in to_cas_cache'.format(cas))
        return to_cas_cache[cas]

    results = pcp.get_synonyms(cas, 'name')

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
    if cas in to_name_cache:
        if verbose:
            print('{} in to_name_cache'.format(cas))
        return to_name_cache[cas]

    results = pcp.get_compounds(cas, 'name')

    if len(results) == 0:
        name = None
        to_name_cache[cas] = name
        return name

    # TODO TODO if results if len > 1, maybe err
    r = results[0]
    name = r.iupac_name

    to_name_cache[cas] = name
    return name


if __name__ == '__main__':
    import pandas as pd

    # test
    print(name2cas('Aspirin'))

    # would need to download this file to run this example
    data = pd.read_csv('Odor_inventory _Sheet1.csv')
    for name in data['Chemical']:
        cas = name2cas(name)
        print(name, cas)

