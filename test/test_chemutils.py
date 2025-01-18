#!/usr/bin/env python3

import pandas as pd

from chemutils import is_one2one


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


