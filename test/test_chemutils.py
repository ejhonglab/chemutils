#!/usr/bin/env python3

import pandas as pd
import pubchempy as pcp

from chemutils import is_one2one, compound2cas, allowed_kwarg
import chemutils as cu


def _is_one2one(xs, ys, **kwargs):
    cols = ['x', 'y']
    df = pd.DataFrame(dict(zip(cols, [xs, ys])))
    r1 = is_one2one(df, *cols, **kwargs)
    r2 = is_one2one(df, *cols[::-1], **kwargs)
    assert r1 == r2
    return r1


def test_is_one2one():
    assert _is_one2one([1,2], [2,4])

    # if either col has >1 unique value for any 1 value of the other col, is_one2one
    # should return False
    assert not _is_one2one([1, 1], [2, 3])

    # duplicate combos of [x,y] should not prevent is_one2one from returning True
    assert _is_one2one([1, 1, 2], [2, 2, 4])


def test_is_one2one_with_null():
    xs1 = [1,None,None]
    xs2 = [1,3,None]
    ys1 = [2,4,8]
    ys2 = [2,0,0]
    for dropna in (True, False):
        assert _is_one2one(xs1, ys2, dropna=dropna)
        assert _is_one2one(xs2, ys1, dropna=dropna)

        assert dropna == _is_one2one(xs1, ys1, dropna=dropna)
        assert dropna == _is_one2one(xs2, ys2, dropna=dropna)


def test_compound2cas():
    # alpha-pinene
    # https://pubchem.ncbi.nlm.nih.gov/compound/6654
    cid = 6654
    compound = pcp.Compound.from_cid(cid)
    # first one i can see under "Synonyms" section on page above is "80-56-8", and this
    # one also seems to list the greatest number of sources under "Other Identifiers" >
    # "CAS" section, though there are also:
    # - 25766-18-1 (CAS Common Chemistry; ChemIDplus)
    # - 2437-95-8 (ECHA)
    # - 7785-70-8 (EPA)
    cas = compound2cas(compound)
    assert cas == '80-56-8'


def test_allowed_kwarg():
    # compound2cas DOES have a verbose=True kwarg (specified as keyword only, via
    # preceding `*, ` in def
    assert allowed_kwarg(compound2cas, 'verbose')

    # compound2name takes no arguments other than a single required positional arg
    assert not allowed_kwarg(cu.compound2name, 'verbose')

    def f1(x, verbose=True):
        return x

    def f2(x, **kwargs):
        return x

    # despite this being as positional argument, can still pass as f3(verbose=True)
    def f3(verbose):
        pass

    # here verbose is positional only, so can not specify it as kwarg
    def f4(verbose, /):
        pass

    assert allowed_kwarg(f1, 'verbose')
    assert allowed_kwarg(f2, 'verbose')
    assert allowed_kwarg(f3, 'verbose')
    assert not allowed_kwarg(f4, 'verbose')

