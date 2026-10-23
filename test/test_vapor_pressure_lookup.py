#!/usr/bin/env python3

import pubchemprops.pubchemprops as pp
import pytest

import chemutils as cu


@pytest.mark.parametrize('string,expected', [
    ('93.2 mm Hg at 25 °C', [25.0]),
    ('1.5X10-4 mm Hg at 25 °C', [25.0]),
    ('75 [mm Hg] @20.399999999999999 °C', [20.4]),
    ('(77°F): 35 mmHg', [25.0]),
    ('40 mm Hg at 85.3° F', [29.61]),
    ('Vapor pressure 4,040 kPa (-1.5 °C)', [-1.5]),
    ('10kPa at -83.8 °C', [-83.8]),
    ('Vapor pressure range 1-760 mm Hg in the temperature range 63.8-233.0 °C',
        [63.8, 233.0]
    ),
    ('5.45 Pa @ 294 deg K', [20.85]),
    ('VP: 18.7 kPa at 300 K', [26.85]),
    ('0.28 [mmHg]', []),
    ('Vapor pressure: 6.3X10-7 mm Hg', []),
    # (boiling point at reduced pressure)
    ('103.2 °C at 10kPa', None),
    ('0.0135 mm Hg at 25 /extrapolated/', None),
    ('VP: 5 mm Hg at 80 to 81 mm Hg', None),
])
def test_stated_temps_c(string, expected):
    temps_c = cu._stated_temps_c(string)
    if expected is None:
        assert temps_c is None
    else:
        assert temps_c == pytest.approx(expected, abs=0.01)


def _mmhg2kpa(mmhg):
    return cu.ureg.Quantity(mmhg, 'mmHg').to('kPa').magnitude


@pytest.mark.parametrize('strings,expected_kpa', [
    # (some of the values PubChem has for ethyl acetate)
    (['73 mmHg at 68 °F ; 100 mmHg at 81 °F (NTP, 1992)', '93.2 [mmHg]',
        '93.2 mm Hg at 25 °C', 'Vapor pressure, kPa at 20 °C: 10', '73 mmHg',
        '75 [mm Hg] @20.4 °C'],
        _mmhg2kpa(93.2)
    ),
    (['Vapor pressure, kPa at 20 °C: 10', '73 mmHg'], 10),
    (['0.28 [mmHg]', '0.3 mmHg'], _mmhg2kpa(0.29)),
    (['1 mmHg at 146.8 °F (NTP, 1992)', '0.06 [mmHg]'], _mmhg2kpa(0.06)),
    (['VP: 0.664 mm Hg at 70 °C', '1 mmHg at 146.8 °F (NTP, 1992)'], None),
    (['10.34 mm Hg (USCG, 1999)'], None),
    (['15.5 mmHg at 77 °F ; 31.3 mmHg at 108 °F (NTP, 1992)',
        'Vapor pressure, Pa at 20 °C: 5'],
        0.005
    ),
    (['15.5 mmHg at 77 °F ; 31.3 mmHg at 108 °F (NTP, 1992)'], _mmhg2kpa(15.5)),
    (['11 mmHg (@ 77 °F)', '7.5 [mm Hg] @18 °C'], _mmhg2kpa(11)),
    (['0.0012 torr (2% aqueous solution)'], None),
])
def test_inchi2vapor_pressure_temps(monkeypatch, strings, expected_kpa):
    results = [{'Value': {'StringWithMarkup': [{'String': s}]}} for s in strings]
    monkeypatch.setattr(pp, 'get_second_layer_props',
        lambda name, props, cid=False: {'Vapor Pressure': results}
    )
    monkeypatch.setattr(cu, 'convert', lambda *args, **kwargs: 1)

    ret = cu.inchi2vapor_pressure('1S/CH4/h1H4')
    if expected_kpa is None:
        assert ret.isnull().all()
    else:
        assert ret['vapor_pressure'] == pytest.approx(expected_kpa)
        assert ret['vapor_pressure_units'] == 'kPa'
