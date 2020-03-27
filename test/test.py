#!/usr/bin/env python3

import chemutils as cu


def test_name2smiles():
    cu.name2smiles('ethanol')


if __name__ == '__main__':
    cu.delete_cache()
    test_name2smiles()
