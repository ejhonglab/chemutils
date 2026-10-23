#!/usr/bin/env python3
"""
Re-computes cached `inchi2vapor_pressure` values with the current lookup, which only
uses PubChem values at or near 25 C (older lookups averaged values at any temperature).

PubChem values are downloaded again for the entries selected with --check (by default,
only non-null, non-manual values). Entries whose value changes are updated, or deleted
if there is no longer any PubChem value to use (including compounds PubChem no longer
has any vapor pressure data for). Entries with the same value, but without (or with a
different) 'vapor_pressure_measured', are also updated. Manually entered values are only
changed if 'manual' is in --check (then saved as non-manual). Null values (checked if
'nulls' is in --check) are updated if there is now a value.

Without --write, only reports what would change. With --write, the cache is first backed
up to <cache path>.<time>.bak
"""

import argparse
import sqlite3
import time

import numpy as np
import pandas as pd
import requests

import chemutils as cu


# Times to try each lookup PubChem responds to with 429 (Too Many Requests).
max_tries = 4

def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-c', '--check', nargs='+', default=['values'],
        choices=['values', 'manual', 'nulls'], help='which entries to check: non-null '
        'non-manual values, non-null manual values, and/or null values (default: '
        '%(default)s)'
    )
    parser.add_argument('-w', '--write', action='store_true',
        help='update / delete changed entries (otherwise only reports them)'
    )
    parser.add_argument('-r', '--report', help='also write report to this CSV')
    parser.add_argument('-i', '--inchis', metavar='CSV', help="only check entries for "
        "InChIs in the 'inchi' column of this CSV (e.g. failed lookups from a report)"
    )
    args = parser.parse_args()

    # (warnings, e.g. from strings failing to parse, still go to cu.log_path)
    cu.logger.removeHandler(cu.sh)

    def entry_kind(value, manual):
        if cu.isnull(value) or cu.isnull(value['vapor_pressure']):
            return 'nulls'
        return 'manual' if manual else 'values'

    inchis = None
    if args.inchis is not None:
        inchis = set(pd.read_csv(args.inchis).inchi)

    conn = cu._get_conn()
    entries = [(inchi, value, manual) for _, _, inchi, value, manual in
        cu._iter_entries(conn, 'inchi', 'vapor_pressure')
        if entry_kind(value, manual) in args.check and (inchis is None or inchi in inchis)
    ]
    print(f'checking {len(entries)} cached vapor pressures ({", ".join(args.check)})')

    if args.write:
        backup_path = f'{cu.cache_path()}.{time.strftime("%Y%m%d_%H%M%S")}.bak'
        backup = sqlite3.connect(backup_path)
        conn.backup(backup)
        backup.close()
        print(f'backed up cache to {backup_path}')

    # (undecorated, so new values are only saved below)
    lookup = cu.inchi2vapor_pressure.__wrapped__
    sub_cache = cu.cache['inchi']['vapor_pressure']
    inchi2cid = cu.cache['inchi']['cid']
    rows = []
    for i, (inchi, old, manual) in enumerate(entries):
        if i % 100 == 0:
            print(f'{i}/{len(entries)}')

        cid = inchi2cid.get(inchi)
        start = time.time()
        new = None
        for n_tries in range(1, max_tries + 1):
            error = None
            try:
                new = lookup(inchi)
                break

            except requests.HTTPError as err:
                error = f'{type(err).__name__}: {err}'
                # PubChem throttling requests (e.g. if others are also using it)
                if err.response is None or err.response.status_code != 429:
                    break
                if n_tries < max_tries:
                    time.sleep(10 * 2**n_tries)

            # (so one failure, e.g. from new PubChem data the parsing can't handle,
            # doesn't stop the rest)
            except Exception as err:
                error = f'{type(err).__name__}: {err}'
                break

        # PubChem asks for at most 5 requests per second. (lookups without a CID don't
        # make any requests)
        if not cu.isnull(cid):
            time.sleep(max(0, 0.2 - (time.time() - start)))

        old_kpa = None if cu.isnull(old) else old['vapor_pressure']
        new_kpa = None if new is None else new['vapor_pressure']
        same_value = cu.isnull(old_kpa) and cu.isnull(new_kpa) or (
            not cu.isnull(old_kpa) and not cu.isnull(new_kpa) and
            np.isclose(old_kpa, new_kpa, rtol=1e-6)
        )
        old_measured = None if cu.isnull(old) else old.get('vapor_pressure_measured')
        # (null values are left without this key)
        same_measured = cu.isnull(new_kpa) or (
            old_measured == new['vapor_pressure_measured']
        )
        if error is None and same_value and same_measured:
            continue

        if error is not None:
            action = f'none (lookup failed: {error})'
        elif manual and 'manual' not in args.check:
            action = 'none (manual)'
        elif cu.isnull(new_kpa):
            action = 'delete'
        elif same_value:
            action = 'update (vapor_pressure_measured only)'
        else:
            action = 'update'

        if args.write and action == 'delete':
            del sub_cache[inchi]
        elif args.write and action.startswith('update'):
            sub_cache[inchi] = new

        rows.append({'inchi': inchi, 'cid': cid, 'manual': bool(manual),
            'old_kpa': old_kpa, 'new_kpa': new_kpa,
            'new_measured': None if new is None else new['vapor_pressure_measured'],
            'action': action, 'time': time.strftime('%Y-%m-%d %H:%M:%S')
        })

    df = pd.DataFrame(rows, columns=['inchi', 'cid', 'manual', 'old_kpa', 'new_kpa',
        'new_measured', 'action', 'time'
    ])
    verb = '' if args.write else 'would be '
    print(f'\nof {len(entries)} checked cached vapor pressures:')
    for action, n in df.action.value_counts().items():
        print(f'{n} {verb}changed (action: {action})')

    changed = df.dropna(subset=['old_kpa', 'new_kpa'])
    log10_change = np.log10(changed.new_kpa.astype(float)) - np.log10(
        changed.old_kpa.astype(float)
    )
    if len(log10_change) > 0:
        print('\nlog10(new / old) for changed non-null values:')
        print(log10_change.describe().to_string(float_format='%.3g'))

    if args.report is not None:
        df.to_csv(args.report, index=False)


if __name__ == '__main__':
    main()
