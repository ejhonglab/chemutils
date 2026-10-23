#!/usr/bin/env python3
"""
Compares `chemutils.inchi2est_vapor_pressure` (SIMPOL.1) estimates at 25 C to the vapor
pressures of all compounds with non-null cached `inchi2vapor_pressure` outputs.

Cached `inchi2vapor_pressure` values are the mean of PubChem values near 25 C (or with no
stated temperature, if there are none), unless cached by older versions, which used
values at any temperature (see scripts/fix_vapor_pressure_cache.py). With --pubchem-25c,
PubChem values are downloaded again for the same compounds, and only values stated to be
within 0.5 C of 25 C (and not estimated / extrapolated) are used. The median is taken
across sources, ignoring NTP values (which are often far off) unless they are the only
ones. Downloaded values are saved to (and re-used from) the JSON file passed to that
option.
"""

import argparse
import json
from pathlib import Path
import re
import time
from typing import Optional

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

import chemutils as cu


estimate_re = re.compile(r'\best|extrapolat|calc|predict', re.IGNORECASE)

def fetch_pubchem_vapor_pressure_strings(cid: int) -> list:
    """Returns strings in the 'Vapor Pressure' section of the PubChem page for `cid`.
    """
    url = (f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON'
        '?heading=Vapor+Pressure'
    )
    response = cu.session.get(url, timeout=30)
    # (no 'Vapor Pressure' section)
    if response.status_code == 404:
        return []
    response.raise_for_status()

    strings = []
    def add_strings(section):
        for info in section.get('Information', []):
            for markup in info.get('Value', {}).get('StringWithMarkup', []):
                strings.append(markup['String'])

        for subsection in section.get('Section', []):
            add_strings(subsection)

    add_strings(response.json()['Record'])
    return strings


def pubchem_25c_vapor_pressure(strings) -> Optional[float]:
    """Returns median vapor pressure [kPa] of values stated at 25 C in `strings`.
    """
    kpa = {True: [], False: []}
    for string in strings:
        if estimate_re.search(string):
            continue

        is_ntp = 'NTP' in string
        # NTP strings can have values at several temperatures, e.g.
        # '1 mm Hg at 41 °F ; 5 mm Hg at 86° F; 100 mm Hg at 196.5° F (NTP, 1992)'
        for part in string.split(';'):
            try:
                parsed = cu.parse_pubchemprops_string(part, expected_units='kPa')
            # (parsing can fail many ways, and these values are just skipped)
            except Exception:
                continue

            temp_c = parsed['temperature_c']
            if parsed['quantity'] is None or temp_c is None or abs(temp_c - 25) > 0.5:
                continue

            value = cu.ureg.Quantity(parsed['quantity'], parsed['units']).to('kPa')
            if value.magnitude > 0:
                kpa[is_ntp].append(value.magnitude)

    values = kpa[False] if len(kpa[False]) > 0 else kpa[True]
    if len(values) == 0:
        return None

    return 10**np.median(np.log10(values))


def add_pubchem_25c(df: pd.DataFrame, json_path: Path) -> pd.DataFrame:
    """Adds 'pubchem_25c' column, downloading strings not already in `json_path`.
    """
    cid2strings = dict()
    if json_path.exists():
        cid2strings = {int(k): v for k, v in json.loads(json_path.read_text()).items()}

    to_fetch = [int(c) for c in df.cid.dropna().unique() if int(c) not in cid2strings]
    for i, cid in enumerate(to_fetch):
        if i % 100 == 0:
            print(f'downloading PubChem vapor pressures ({i}/{len(to_fetch)})')

        cid2strings[cid] = fetch_pubchem_vapor_pressure_strings(cid)
        # PubChem asks for at most 5 requests per second.
        time.sleep(0.25)

        if i % 50 == 0 or i == len(to_fetch) - 1:
            json_path.write_text(json.dumps(cid2strings))

    df['pubchem_25c'] = [
        np.nan if pd.isnull(c) else pubchem_25c_vapor_pressure(cid2strings[int(c)])
        for c in df.cid
    ]
    return df


# Structural classes (besides SIMPOL.1 groups) to also summarize errors for.
class_smarts = {
    'furan (aromatic O)': '[o]',
    'tertiary alcohol': '[OX2H][CX4]([#6])([#6])[#6]',
    'phenol with ortho O': '[OX2H]c:c[$([OX2]),$([#6]=O)]',
    'alkyl ether': '[OX2;!R]([CX4])[CX4]',
}

def print_error_table(title: str, err: pd.Series, masks: dict, min_n: int = 1
    ) -> None:

    rows = []
    for name, mask in masks.items():
        if mask.sum() >= min_n:
            rows.append({title: name, 'n': mask.sum(),
                'median_err': err[mask].median(),
                'median_abs_err': err[mask].abs().median(),
            })

    print()
    print(pd.DataFrame(rows).to_string(index=False, float_format='%.2f'))


def print_summary(df: pd.DataFrame, ref_col: str, n_outliers: int) -> None:
    df = df.dropna(subset=[ref_col]).copy()
    ref = np.log10(df[ref_col].astype(float))

    # Compounds inchi2est_vapor_pressure excludes, with SIMPOL.1 estimate errors.
    excluded = df.excluded.notna()
    excluded_err = np.log10(df.simpol[excluded]) - ref[excluded]
    excluded_df = df[excluded]

    df = df[~excluded].copy()
    err = np.log10(df.est) - ref[~excluded]
    df['log10_err'] = err

    print(f'\nestimate vs {ref_col} (errors in log10 units, positive if estimate '
        f'higher):\n{len(df)} compounds, median error {err.median():+.2f}, median '
        f'absolute error {err.abs().median():.2f}, ' + ', '.join(
            f'{(err.abs() <= np.log10(f)).mean():.0%} within {f}x' for f in (2, 3, 10)
        )
    )

    print_error_table('SIMPOL.1 group (>=5 compounds)', err,
        {g: df[g] > 0 for g in cu._simpol_coefs if g not in ('zeroeth', 'carbon')},
        min_n=5
    )

    mols = [Chem.MolFromSmiles(s) for s in df.smiles]
    print_error_table('other class', err, {name: np.array(
            [m.HasSubstructMatch(Chem.MolFromSmarts(smarts)) for m in mols]
        ) for name, smarts in class_smarts.items()
    })

    print_error_table('excluded (SIMPOL.1 errors)', excluded_err,
        {r: (excluded_df.excluded == r).values for r in excluded_df.excluded.unique()}
    )

    print(f'\n{n_outliers} largest errors:')
    outliers = df.loc[err.abs().sort_values(ascending=False).index[:n_outliers]]
    print(outliers[['cid', 'smiles', ref_col, 'est', 'log10_err']].to_string(
        index=False, float_format='%.3g'
    ))


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--pubchem-25c', type=Path, metavar='JSON',
        help='also compare to PubChem values stated at 25 C (downloading those not '
        'already in this file)'
    )
    parser.add_argument('-o', '--output', help='write per-compound results to this CSV')
    parser.add_argument('-n', '--n-outliers', type=int, default=20,
        help='number of largest errors to print (default: %(default)s)'
    )
    args = parser.parse_args()

    RDLogger.DisableLog('rdApp.*')

    rows = []
    n_unsupported = 0
    for inchi, value in cu.cache['inchi']['vapor_pressure'].items():
        if cu.isnull(value) or cu.isnull(value['vapor_pressure']):
            continue

        assert value['vapor_pressure_units'] == 'kPa'
        mol = Chem.MolFromInchi(cu._inchi_with_prefix(inchi))
        groups = None if mol is None else cu._simpol_groups(mol)
        if groups is None:
            n_unsupported += 1
            continue

        excluded = cu._simpol_exclusion_reason(mol, groups)
        rows.append({'inchi': inchi,
            'cid': cu.cache['inchi']['cid'].get(inchi),
            'smiles': Chem.MolToSmiles(mol),
            'cached': value['vapor_pressure'],
            # (so errors for excluded compounds can also be summarized)
            'simpol': 10**cu._simpol_log10_p_atm(groups, 298.15) * 101.325,
            'excluded': excluded,
            'est': np.nan if excluded is not None else
                cu.inchi2est_vapor_pressure(inchi)['vapor_pressure'],
            **groups
        })

    df = pd.DataFrame(rows)
    print(f'{len(df) + n_unsupported} compounds with cached inchi2vapor_pressure '
        f'values ({n_unsupported} not supported by SIMPOL.1, '
        f'{df.excluded.notna().sum()} excluded by inchi2est_vapor_pressure)'
    )

    print_summary(df, 'cached', args.n_outliers)

    if args.pubchem_25c is not None:
        df = add_pubchem_25c(df, args.pubchem_25c)
        print_summary(df, 'pubchem_25c', args.n_outliers)

    if args.output is not None:
        df.to_csv(args.output, index=False)


if __name__ == '__main__':
    main()
