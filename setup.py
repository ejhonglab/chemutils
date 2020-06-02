#!/usr/bin/env python3

from setuptools import setup


setup(
    name='chemutils',
    version='0.1.0',
    py_modules=['chemutils'],
    # xlrd is so that pandas can actually write Excel files.
    install_requires=[
        'pandas', 'pubchempy', 'pint', 'xlrd'
    ]
)
