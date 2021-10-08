
### chemutils

Convenience functions wrapping `pubchempy` library, supporting common operations
we may want to do with the chemical components of odor mixtures.

### Installation
```
pip install --user git+https://github.com/ejhonglab/chemutils.git
```
Replace `pip` with `pip3` if you are on a computer with both Python 2 and 3, and
you want to use this with Python 3.

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
