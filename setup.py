#!/usr/bin/env python3

from setuptools import setup


setup(
    name='chemutils',
    version='0.1.0',
    py_modules=['chemutils'],
    # TODO specify my fork of pubchemprops
    install_requires=[
        'pandas', 'pubchempy', 'pint', 'openpyxl', 'XlsxWriter',
        # TODO add current verison of requests
        # TODO add bs4
        # This version seems to have working InChI support installed directly from PyPi.
        # Earlier versions probably do as well but untested.
        'rdkit>=2023.3.1',
    ]
)
