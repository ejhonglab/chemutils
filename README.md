
### chemutils

Convenience functions wrapping `pubchempy` library, supporting common operations
we may want to do with the chemical components of odor mixtures.

### Installation
```
pip install git+https://github.com/ejhonglab/chemutils
```

### Examples
```
import chemutils as cu

cas_number_str = cu.convert('ethyl acetate', from_type='name', to_type='cas')
```

Warnings are printed (with the stack of the code that issued them) and appended to
`~/.chemutils.log` (or the path in the `CHEMUTILS_LOG` environment variable). To stop
printing them: `cu.logger.removeHandler(cu.sh)`.

### Vapor pressure estimates
For compounds where `cu.inchi2vapor_pressure` finds no PubChem value,
`cu.inchi2est_vapor_pressure` estimates it from structure alone (SIMPOL.1 group
contribution method, offline). See its docstring for accuracy and which compounds it
supports. To compare its estimates to PubChem values:
```
python scripts/validate_est_vapor_pressure.py --pubchem-25c pubchem_vp_strings.json
```

### Henry's law constants in paraffin oil
`cu.inchi2k_henry_pfo` estimates Henry's law solubility constants in paraffin oil (pfo),
in the same units as `cu.inchi2k_henry` values from henrys-law.org (mol / m^3 / Pa), with
modified UNIFAC in an n-alkane surrogate for the oil (as `oil_headspace.py` does). No
measured solubilities in oil or alkanes are used unless passed in, and the
`k_henry_pfo_method` and `k_henry_pfo_measured_partitioning` outputs say which was used.

### Cache
Lookup results are cached in an SQLite database at `~/.chemutils_cache.sqlite` (set the
`CHEMUTILS_CACHE` environment variable to use another path). Input IDs are normalized
before lookup (e.g. `cu.normalize_name`, `cu.normalize_inchi`), and values are cached
under the normalized ID, so raw and normalized inputs share the same cache entries.

To migrate the old pickle cache (`~/.chemutils_cache.p`, which is left unchanged):
```
python scripts/migrate_cache.py --dry-run   # to see what will change
python scripts/migrate_cache.py
```

Vapor pressures cached by older versions averaged PubChem values at any temperature. To
re-compute them with only values at or near 25 C (the cache is backed up first):
```
python scripts/fix_vapor_pressure_cache.py --check values manual nulls          # to see what will change
python scripts/fix_vapor_pressure_cache.py --check values manual nulls --write
```
Changes this made to the cache it was first run on (2026-10-23) are listed in
[cache_changes/2026-10-23_vapor_pressure.md](cache_changes/2026-10-23_vapor_pressure.md).

### Tests
Tests use a temporary cache (never the one above). Tests that need internet access are
marked `network`. To skip those: `pytest -m "not network"`.
