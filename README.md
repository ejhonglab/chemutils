
### chemutils

Convenience functions wrapping `pubchempy` library, supporting common operations
we may want to do with the chemical components of odor mixtures.

### Installation
```
pip install --user git+https://github.com/ejhonglab/chemutils.git
```
Replace `pip` with `pip3` if you are on a computer with both Python 2 and 3, and
you want to use this with Python 3.

### Examples
```
import chemutils

cas_number_str = chemutils.name2cas('ethyl acetate')
```
