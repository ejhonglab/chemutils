#!/usr/bin/env python3
"""
Migrates the old pickle chemutils cache (~/.chemutils_cache.p) to the SQLite cache
(~/.chemutils_cache.sqlite, or the path in the CHEMUTILS_CACHE environment variable).

The pickle is only read (never modified or deleted).

Changes made to the cache contents while migrating:
- Keys are re-normalized with the current `normalize_<type>` fns, and entries whose keys
  now normalize to the same key are merged. When merging, manually entered values take
  precedence, then non-null values, then the value under the key that was already in
  normalized form. Differing non-null values are reported as conflicts.

- Name keys produced by the old normalize_name joining words that it should not have
  (e.g. 'ethyl3-hydroxyhexanoate', from 'ethyl 3-hydroxyhexanoate') are dropped
  entirely if their CID differs from the CID of the name they came from (all of their
  values are for the wrong compound). Otherwise, for each value that is null for only
  one of the two names, the null entry is dropped (so it gets looked up again).

- Null values for other ID types (cached because a CID lookup failed) are dropped for
  keys that have a non-null CID after merging.

- CID values are converted to int, and InChI / InChIKey / CAS values normalized.

- Entries with null keys, entries identical to `chemutils.hardcoded` values (which are
  now applied separately, and not stored in the cache), and property values with
  inconsistent units are dropped.
"""

import argparse
from collections import defaultdict
import contextlib
import csv
from dataclasses import dataclass
import io
import os
from pathlib import Path
import pickle
import re
from typing import Optional

import pandas as pd

# Also needed before unpickling, so pint's application registry is chemutils.ureg
import chemutils as cu


# Types whose values (not just keys) are normalized while migrating.
normalized_value_types = ('cid', 'inchi', 'inchikey', 'cas')


@dataclass
class Entry:
    from_type: str
    to_type: str
    raw_key: object
    value: object
    manual: bool = False


def _quiet(fn, *args):
    # (normalize_name prints when applying some corrections)
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*args)


def old_normalize_name(name) -> Optional[str]:
    """`chemutils.normalize_name` as of commit 08e7dc4 (before word-join fix).

    Only used to identify cache keys created by its incorrect joining of words.
    """
    if pd.isnull(name):
        return name

    name = name.replace('(CAS)', '').strip().lower()
    name = re.sub(r'\([^\(\)]{4,}\)$', '', name)

    parts = name.split()
    if len(parts) == 0:
        return None

    normed_name = parts[0]
    for a, b in zip(parts, parts[1:]):
        if a[-1] == ',' or (a[-1].isalpha() and b[0].isalpha()):
            b = ' ' + b
        normed_name += b

    corrections = {
        'linalool oxide': 'trans-linalool oxide',
        '4-ethyl guaiacol': '4-ethylguaiacol',
        'dihydrazide ethanediimidic acid': 'ethanediimidic acid, dihydrazide',
        'e3-hexenol': '(E)-hex-3-en-1-ol',
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
    for p, c in prefix_corrections.items():
        if normed_name.startswith(p):
            normed_name = c + normed_name[len(p):]

    if ' ester ' in normed_name and ' acid' in normed_name:
        parts = normed_name.split()
        normed_name = ' '.join([parts[2], parts[3] + ',', parts[0], parts[1]])

    return normed_name


def load_pickle_cache(path):
    """Returns (cache, is_manual) dicts from old pickle cache at `path`.
    """
    with open(path, 'rb') as f:
        data = pickle.load(f)

    is_manual = data.get('is_manual')
    if is_manual is None:
        is_manual = dict()

    return data['cache'], is_manual


def _hardcoded_raw_entries():
    """Returns {(from_type, to_type, raw key): value} the old cache added from hardcoded.
    """
    entries = dict()
    for ft, null_keys in cu.hardcoded_type2null_keys.items():
        for k in null_keys:
            for tt in cu.chem_id_types:
                entries[(ft, tt, k)] = None

    for ft, tts in cu.hardcoded.items():
        for tt, f2t_dict in tts.items():
            for k, v in f2t_dict.items():
                entries[(ft, tt, k)] = v

    return entries


def _drop_old_name_normalization_artifacts(entries, report):
    name2entries = defaultdict(dict)
    for e in entries:
        if e.from_type == 'name' and isinstance(e.raw_key, str):
            name2entries[e.raw_key][e.to_type] = e

    # (to_type, raw name)
    to_drop = set()
    for raw, raw_entries in name2entries.items():
        try:
            old_normed = _quiet(old_normalize_name, raw)
        # (old fn raised this for some inputs, so they can't have been added by it)
        except IndexError:
            continue

        if old_normed is None or old_normed == raw or old_normed not in name2entries:
            continue

        # Keys that now normalize the same are merged later instead.
        if (_quiet(cu.normalize_name, raw) ==
            _quiet(cu.normalize_name, old_normed)):
            continue

        report['old_normalization_artifacts'].append((raw, old_normed))
        old_normed_entries = name2entries[old_normed]

        raw_cid = getattr(raw_entries.get('cid'), 'value', None)
        old_normed_cid = getattr(old_normed_entries.get('cid'), 'value', None)
        if (not cu.isnull(raw_cid) and not cu.isnull(old_normed_cid) and
            raw_cid != old_normed_cid):

            report['wrong_cid_dropped'].append(
                (raw, raw_cid, old_normed, old_normed_cid)
            )
            to_drop.update((tt, old_normed) for tt in old_normed_entries)
            continue

        for tt in set(raw_entries) & set(old_normed_entries):
            raw_value = raw_entries[tt].value
            old_normed_value = old_normed_entries[tt].value
            if cu.isnull(old_normed_value) and not cu.isnull(raw_value):
                report['one_null_dropped'].append(('name', tt, old_normed))
                to_drop.add((tt, old_normed))

            elif cu.isnull(raw_value) and not cu.isnull(old_normed_value):
                report['one_null_dropped'].append(('name', tt, raw))
                to_drop.add((tt, raw))

    return [e for e in entries
        if not (e.from_type == 'name' and (e.to_type, e.raw_key) in to_drop)
    ]


def _merge(key, members, report) -> Entry:
    """Returns which of `members` (old entries normalizing to `key`) to keep.
    """
    candidates = [m for m in members if m.manual]
    if len(candidates) == 0:
        candidates = members

    nonnull = [m for m in candidates if not cu.isnull(m.value)]
    if len(nonnull) == 0:
        return candidates[0]

    if len(nonnull) < len(members):
        report['null_merged'].append(key)

    distinct = []
    for m in nonnull:
        if not any(cu._values_equal(m.value, d.value) for d in distinct):
            distinct.append(m)

    if len(distinct) == 1:
        return nonnull[0]

    # Preferring value that was looked up under the key itself.
    already_normed = [m for m in nonnull if m.raw_key == key[2]]
    chosen = already_normed[0] if len(already_normed) > 0 else nonnull[0]
    report['conflicts'].append(
        (key, [(m.raw_key, m.value) for m in distinct], chosen.raw_key)
    )
    return chosen


def plan_migration(old_cache, old_is_manual):
    """Returns ({(from_type, to_type, key): Entry} to write, report).

    `report` is a dict of category -> list of details (see `print_report`).
    """
    report = defaultdict(list)
    hardcoded = _hardcoded_raw_entries()

    entries = []
    for ft, tts in old_cache.items():
        for tt, f2t_dict in tts.items():
            for k, v in f2t_dict.items():
                if cu.isnull(k):
                    report['null_key_dropped'].append((ft, tt, k))
                    continue

                if (ft, tt, k) in hardcoded and cu._values_equal(
                    hardcoded[(ft, tt, k)], v):

                    report['hardcoded_dropped'].append((ft, tt, k))
                    continue

                manual = bool(old_is_manual.get((ft, tt, k), False))
                entries.append(Entry(ft, tt, k, v, manual))

    report['n_input'] = len(entries) + len(report['null_key_dropped']) + len(
        report['hardcoded_dropped']
    )

    entries = _drop_old_name_normalization_artifacts(entries, report)

    groups = defaultdict(list)
    for e in entries:
        key = _quiet(cu.normalize_id, e.from_type, e.raw_key)
        if cu.isnull(key):
            report['null_normalized_key_dropped'].append(
                (e.from_type, e.to_type, e.raw_key)
            )
            continue

        if _quiet(cu.normalize_id, e.from_type, key) != key:
            report['non_idempotent'].append((e.from_type, e.raw_key, key))

        if key != e.raw_key:
            report['rekeyed'].append((e.from_type, e.to_type, e.raw_key, key))

        if e.to_type in normalized_value_types and not cu.isnull(e.value):
            e.value = cu.normalize_id(e.to_type, e.value)

        groups[(e.from_type, e.to_type, key)].append(e)

    resolved = {k: _merge(k, members, report) for k, members in groups.items()}

    # Dropping null values cached because a CID lookup failed, if we now have a CID.
    for (ft, tt, key), e in list(resolved.items()):
        if (ft == 'cid' or tt == 'cid' or tt not in cu.chem_id_types or
            not cu.isnull(e.value)):
            continue

        cid_entry = resolved.get((ft, 'cid', key))
        if cid_entry is None or cu.isnull(cid_entry.value):
            continue

        old_f2cid = old_cache.get(ft, dict()).get('cid', dict())
        if e.raw_key in old_f2cid and cu.isnull(old_f2cid[e.raw_key]):
            report['stale_null_dropped'].append((ft, tt, key))
            del resolved[(ft, tt, key)]

    for (ft, tt, key), e in list(resolved.items()):
        try:
            cu._check_property_value(tt, e.value)
        except ValueError as err:
            report['invalid_property_dropped'].append(((ft, tt, key), str(err)))
            del resolved[(ft, tt, key)]

    report['n_output'] = len(resolved)
    return resolved, report


def print_report(report, max_examples=10) -> None:
    descriptions = {
        'null_key_dropped': 'entries with null keys dropped',
        'hardcoded_dropped': 'entries identical to chemutils.hardcoded values dropped',
        'old_normalization_artifacts': ('(raw name, old normalized name) pairs, where '
            'old normalize_name had joined words'
        ),
        'wrong_cid_dropped': ('old normalized names dropped (with all their values) '
            'because CID differed from raw name CID'
        ),
        'one_null_dropped': ('entries dropped because null for only one of a raw / '
            'old normalized name pair'
        ),
        'null_normalized_key_dropped': 'entries whose keys normalize to null dropped',
        'non_idempotent': 'keys where normalization is not idempotent (should be 0)',
        'rekeyed': 'entries whose keys changed with normalization',
        'null_merged': 'merged keys where null values were dropped for non-null ones',
        'conflicts': 'merged keys with conflicting non-null values',
        'stale_null_dropped': ('null values (from failed CID lookup) dropped for '
            'keys with a CID'
        ),
        'invalid_property_dropped': 'property values with unexpected units dropped',
    }
    print(f'{report["n_input"]} entries in old cache')
    for category, description in descriptions.items():
        details = report.get(category, [])
        print(f'{len(details)} {description}')

        # Printing all of the ones that may need manual attention.
        n = len(details) if category in ('wrong_cid_dropped', 'conflicts',
            'invalid_property_dropped', 'non_idempotent') else max_examples

        for d in details[:n]:
            print(f'    {d!r}')
        if len(details) > n:
            print('    ...')

    print(f'{report["n_output"]} entries to write')


def write_report_csv(report, path) -> None:
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['category', 'details'])
        for category, details in report.items():
            if not isinstance(details, list):
                continue
            for d in details:
                writer.writerow([category, repr(d)])


def write_db(db_path, resolved, *, updated: float, force: bool = False) -> None:
    """Writes `resolved` entries (from `plan_migration`) to a new cache at `db_path`.
    """
    conn = cu._connect(db_path)
    try:
        n_existing = conn.execute('SELECT COUNT(*) FROM entries').fetchone()[0]
        if n_existing > 0 and not force:
            raise RuntimeError(f'{db_path} already has {n_existing} entries. pass '
                '--force to replace them'
            )

        rows = [
            (ft, tt, key, cu._encode_value(e.value), int(e.manual), updated)
            for (ft, tt, key), e in resolved.items()
        ]
        conn.execute('BEGIN')
        conn.execute('DELETE FROM entries')
        conn.executemany('INSERT INTO entries VALUES (?, ?, ?, ?, ?, ?)', rows)
        conn.execute('COMMIT')
    finally:
        conn.close()


def verify_db(db_path, resolved) -> None:
    """Raises ValueError if database contents do not match `resolved` entries.
    """
    conn = cu._connect(db_path)
    try:
        rows = conn.execute('SELECT from_type, to_type, key, value, manual FROM entries'
            ).fetchall()
    finally:
        conn.close()

    if len(rows) != len(resolved):
        raise ValueError(f'{db_path} has {len(rows)} entries. expected '
            f'{len(resolved)}'
        )

    for ft, tt, key, value, manual in rows:
        e = resolved[(ft, tt, key)]
        if not cu._values_equal(cu._decode_value(value), e.value):
            raise ValueError(f'{(ft, tt, key)} value in {db_path} does not match '
                f'expected {e.value!r}'
            )

        if manual != int(e.manual):
            raise ValueError(f'{(ft, tt, key)} manual flag in {db_path} != {e.manual}')


def migrate(pickle_path, db_path, *, force: bool = False, dry_run: bool = False,
    report_path=None, verbose: bool = True):
    """Migrates old cache at `pickle_path` to new cache at `db_path`.

    Returns (resolved entries, report). See module docstring for details.
    """
    old_cache, old_is_manual = load_pickle_cache(pickle_path)
    resolved, report = plan_migration(old_cache, old_is_manual)

    if verbose:
        print_report(report)

    if report_path is not None:
        write_report_csv(report, report_path)

    if dry_run:
        return resolved, report

    # So entries reflect when they were (last possibly) saved in the old cache.
    write_db(db_path, resolved, updated=os.path.getmtime(pickle_path), force=force)
    verify_db(db_path, resolved)

    if verbose:
        print(f'wrote {len(resolved)} entries to {db_path}')

    return resolved, report


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--pickle', default=cu.old_pickle_cache_file,
        help='old cache to read (default: %(default)s)'
    )
    parser.add_argument('--db', default=None, help='new cache to write (default: '
        '$CHEMUTILS_CACHE, if set, otherwise %s)' % cu.default_cache_file
    )
    parser.add_argument('-n', '--dry-run', action='store_true',
        help='only print what would be migrated'
    )
    parser.add_argument('-f', '--force', action='store_true',
        help='replace any existing entries in new cache'
    )
    parser.add_argument('-r', '--report', help='also write report to this CSV')
    args = parser.parse_args()

    db_path = Path(args.db) if args.db is not None else cu.cache_path()
    try:
        migrate(args.pickle, db_path, force=args.force, dry_run=args.dry_run,
            report_path=args.report
        )
    except RuntimeError as err:
        raise SystemExit(str(err))


if __name__ == '__main__':
    main()
