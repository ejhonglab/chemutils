
### chemutils

Convenience functions wrapping `pubchempy` library, supporting common operations
we may want to do with the chemical components of odor mixtures.

### Installation
```
pip install git+https://github.com/ejhonglab/chemutils
```

Functions using `pubchemprops` (anything calling `cu.pubchemprops_lookup`) will require
that you have a copy of [my fork](https://github.com/tom-f-oconnell/pubchemprops) of
`pubchemprops` installed. You can install it with:
```
pip install git+https://github.com/tom-f-oconnell/pubchemprops
```

### Examples
```
import chemutils as cu

cas_number_str = cu.name2cas('ethyl acetate')
```

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

### Tests
Tests use a temporary cache (never the one above). Tests that need internet access are
marked `network`. To skip those: `pytest -m "not network"`.
