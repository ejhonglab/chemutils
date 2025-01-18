#!/usr/bin/env python3
# NOTE: currently ignored via ../pytest.ini

"""
Using some data from lookups on strawberry mean sample
(see `natural_odors/strawberry.py`) to test parsing some strings `pubchemprops`
returns into quantities, units, conditions they were measured at, and sources.
"""

from pprint import pprint

import chemutils as cu


# TODO TODO actually turn these into tests

# "In water, 47,000 mg/L at 20 °C; 38,000 mg/L at 100 °C" (l42)
# (harcode if not going to support)
def test_solubility_parsing():
    with open('pubchemprops_solubility_strings.txt', 'r') as f:
        lines = f.readlines()

    expected_units = 'concentration'

    start_at = 19
    for i, pstr in enumerate(lines):
        if (i + 1) < start_at:
            continue

        print(f'Line {i+1}')
        print(f'"{pstr.strip()}"')
        parsed_dict = cu.parse_pubchemprops_string(pstr,
            expected_units=expected_units
        )
        if all([x is None for x in parsed_dict.values()]):
            print('<no data>')
        else:
            pprint(parsed_dict)

        if start_at > 0:
            import ipdb; ipdb.set_trace()

        print()


# TODO refactor to share more w/ above
def test_vapor_pressure_parsing():
    with open('pubchemprops_vapor_pressure_strings.txt', 'r') as f:
        lines = f.readlines()

    expected_units = 'kPa'

    start_at = 0
    for i, pstr in enumerate(lines):
        if (i + 1) < start_at:
            continue

        print(f'Line {i+1}')
        print(f'"{pstr.strip()}"')
        parsed_dict = cu.parse_pubchemprops_string(pstr,
            expected_units=expected_units
        )
        if all([x is None for x in parsed_dict.values()]):
            print('<no data>')
        else:
            pprint(parsed_dict)

        if start_at > 0:
            import ipdb; ipdb.set_trace()

        print()


if __name__ == '__main__':
    #test_solubility_parsing()
    test_vapor_pressure_parsing()

