#!/usr/bin/env python3

"""
Takes chemical name as the lone command line argument.

Prints PubChem URL for this chemical and the volume soluble in 2mL and 20mL of
water.
"""

import argparse
from os import linesep
from os.path import isfile

import chemutils as cu


# TODO load chemspider api key from some default places / env var for lookups
# that can use that (make fn in chemutils that can do this)

def ul_soluble_in_ml_water(ml_water, water_solubility, density=None,
    molar_mass=None):
    """
    Args:
    ml_water (float): milliliters of water chemical will be dissolved in

    water_solubility (pint.Quantity): mass / volume water solubility of target
        solute

    density (pint.Quantity): mass / volume density of target solute
    """

    if density is None or molar_mass is not None:
        raise NotImplementedError

    # pint docs claim you shouldn't try to combine quantities constructed from
    # different UnitRegistry objects
    ureg = cu.ureg

    # vv = volume/volume.
    vv_water_sol = water_solubility / density
    # (mass/vol) / (mass/vol) => dimensionless
    assert vv_water_sol.dimensionless

    water_vol = ml_water * ureg('mL')
    return (water_vol * vv_water_sol).to('uL').m


def print_vol_soluble_in_water_vols(name_or_names, water_vols_ml=(2.0, 20.0),
    ignore_cache=False, inchi=True, pubchem_url=True, properties=True,
    verbose=False, line_prefix=None, between_chemicals='\n\n',
    failed_at_end=True, _recursing=False):
    """
    Looks up properties of chemical(s) and prints the maximal volume that is
    soluble in each volume in `water_vols_ml`.

    Returns `bool` indicating whether (all) lookup(s) were successful.
    """

    if type(name_or_names) is not str:
        if _recursing:
            raise ValueError('name_or_names must be str or iterable of str')

        if line_prefix is None:
            line_prefix = ' '

        all_success = True
        failed_names = []
        for i, name in enumerate(name_or_names):
            # TODO maybe print indicator next to name / colorize (if available)
            # if any lookups failed (would need to return str rather than
            # printing or capture prints somehow...)
            success, report = print_vol_soluble_in_water_vols(name,
                water_vols_ml=water_vols_ml, ignore_cache=ignore_cache,
                inchi=inchi, pubchem_url=pubchem_url, properties=properties,
                verbose=verbose, line_prefix=line_prefix, _recursing=True
            )
            if success:
                print(name)
            else:
                print(name + ' (FAILED!)')
                failed_names.append(name)

            print(report)
            if i < (len(name_or_names) - 1):
                print(between_chemicals, end='')

            all_success = all_success and success

        if not all_success and failed_at_end:
            print('\nLookups failed for:')
            for fn in failed_names:
                print(' ' + fn)
            print()

        return all_success
    else:
        name = name_or_names
        if line_prefix is None:
            line_prefix = ''

    lines_out = []
    def lprint(*args):
        strs = tuple(str(a) for a in args)
        lines_out.append(line_prefix + ' '.join(strs))

    def retvals(success):
        report = linesep.join(lines_out)
        if _recursing:
            return success, report
        else:
            print(report, end='')
            return success

    # TODO once i get a fn to lookup name at the top of the pubchem page
    # working, print that name if it differs from input name

    # TODO also warn if lookup is non-unique? does that already?
    # (or print answers for everything inchi lookup returns)
    # (should this have happened when i entered 'hexenal' instead of
    # 'hexanal', by accident?)
    _inchi = cu.convert(name, from_type='name', to_type='inchi',
        verbose=verbose
    )
    if _inchi is None:
        lprint('PubChem based InChI lookup failed!')
        return retvals(False)

    if inchi:
        lprint(_inchi)

    if pubchem_url:
        url = cu.inchi2pubchem_url(_inchi)
        assert type(url) is str
        lprint(url)

    if inchi or pubchem_url:
        lprint()

    density = cu.inchi2density(_inchi, ignore_cache=ignore_cache,
        verbose=verbose, pint=True
    )
    failed = False
    if density is None:
        lprint('Density lookup failed!')
        failed = True
    elif properties:
        lprint(f'Density: {density:.3f}')

    # TODO TODO does function that looks this up also have temperature
    # information at some point? try to preserve that in the output, so it can
    # be printed here
    water_solubility = cu.inchi2water_solubility(_inchi,
        ignore_cache=ignore_cache, verbose=verbose, pint=True
    )
    if water_solubility is None:
        lprint('Water solubility lookup failed!')
        failed = True
    elif properties:
        lprint(f'Water solubility: {water_solubility:.3f}')
        lprint()

    if failed:
        return retvals(False)

    # would only be necessary if solubility was reported in molar units
    '''
    molar_mass = cu.inchi2molar_mass(_inchi, ignore_cache=ignore_cache,
        verbose=verbose, pint=True
    )
    print(f'Molar mass: {molar_mass:.1f}')
    '''

    # (or handle case where it's of different dimensionality)
    ds = '[length]^-3 [mass]'
    # TODO test this assert actually prints out correctly when triggered
    assert water_solubility.check(ds), (f'expected dimensionality of {ds} for '
        f'water solubility, but got {water_solubility.dimensionality}'
    )

    for water_vol_ml in (2.0, 20.0):
        ul_soluble = ul_soluble_in_ml_water(water_vol_ml, water_solubility,
            density
        )
        # TODO maybe scientific notion, or some other guarantee all digits won't
        # be zeros? / limit sigfigs if lots of non-zero digits?
        lprint(f'Up to {ul_soluble:.2f} uL of {name} soluble in '
            f'{water_vol_ml:.1f} mL of water'
        )
    return retvals(True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',
        help='Chemical name or a file with one chemical name per line.'
    )
    name_cols = ['name'] + [k for k, v in cu.equivalent_col_names.items()
        if v == 'name'
    ]
    parser.add_argument('-f', '--file', default=False, action='store_true',
        help='Indicates that input is a filename with one chemical name per '
        f'line, or a CSV with exactly one of {name_cols} among its columns.'
    )
    parser.add_argument('-i', '--ignore-cache', default=False,
        action='store_true', help='All lookups are done from scratch, ignoring '
        'local cache of lookup results.'
    )
    parser.add_argument('-v', '--verbose', default=False,
        action='store_true', help='Prints some information during lookups.'
    )
    args = parser.parse_args()

    if not args.file:
        name_or_names = args.input
    else:
        fname = args.input
        if not isfile(fname):
            raise IOError('-f/--file indicated input is a file, but '
                f'{fname} does not exist'
            )

        if fname.lower().endswith('.csv'):
            # Just in case dependence on pandas is weakened in chemutils.py in
            # the future...
            import pandas as pd
            
            data = pd.read_csv(fname)
            nc = [c for c in data.columns if c in name_cols]
            if len(nc) != 1:
                print(f'Exactly one of CSV columns must be in {name_cols}')
                return
            name_col = nc[0]
            name_or_names = data[name_col]
        else:
            with open(fname, 'w') as f:
                lines = [x.strip() for x in f.readlines()]
                name_or_names = [x for x in lines if len(x) > 0]

    print_vol_soluble_in_water_vols(name_or_names,
        ignore_cache=args.ignore_cache, verbose=args.verbose
    )


if __name__ == '__main__':
    main()

