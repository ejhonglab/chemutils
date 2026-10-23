#!/usr/bin/env python3

import pytest
from rdkit import Chem

import chemutils as cu


def _inchi(smiles: str) -> str:
    return Chem.MolToInchi(Chem.MolFromSmiles(smiles))


# Examples from the footnotes of Table 5 in Pankow & Asher (2008), with the expected
# counts of all groups except 'zeroeth' (always 1) and 'carbon'.
paper_examples = {
    # acetamide
    'CC(N)=O': {'carbon_acid_side_amide': 2, 'amide_primary': 1},
    # methyl acetamide
    'CNC(C)=O': {'carbon_acid_side_amide': 2, 'amide_secondary': 1},
    # methyl ethyl acetamide
    'CCN(C)C(C)=O': {'carbon_acid_side_amide': 2, 'amide_tertiary': 1},
    # N-propyl-butyramide
    'CCCNC(=O)CCC': {'carbon_acid_side_amide': 4, 'amide_secondary': 1},
    # biphenyl
    'c1ccc(-c2ccccc2)cc1': {'aromatic_ring': 2},
    # anthracene
    'c1ccc2cc3ccccc3cc2c1': {'aromatic_ring': 3},
    # cyclohexane
    'C1CCCCC1': {'nonaromatic_ring': 1},
    # 1-hexene
    'C=CCCCC': {'c_c': 1},
    # cyclohex-2-enone
    'O=C1C=CCCC1': {'nonaromatic_ring': 1, 'c_c': 1, 'c_c_c_o_in_ring': 1, 'ketone': 1},
    # methyl acetate
    'COC(C)=O': {'ester': 1},
    # methyl nitroacetate
    'COC(=O)C[N+](=O)[O-]': {'nitro': 1, 'nitroester': 1},
    # dimethyl ether
    'COC': {'ether': 1},
    # tetrahydrofuran
    'C1CCOC1': {'nonaromatic_ring': 1, 'ether_alicyclic': 1},
    # methyl phenyl ether
    'COc1ccccc1': {'aromatic_ring': 1, 'ether_aromatic': 1},
    # 2-nitropropane
    'CC(C)[N+](=O)[O-]': {'nitro': 1},
    # nitrobenzene
    'O=[N+]([O-])c1ccccc1': {'aromatic_ring': 1, 'nitro': 1},
    # 2-methylphenol
    'Cc1ccccc1O': {'aromatic_ring': 1, 'aromatic_hydroxyl': 1},
    # 2-nitrophenol
    'O=[N+]([O-])c1ccccc1O': {'aromatic_ring': 1, 'nitro': 1, 'nitrophenol': 1},
    # methylamine
    'CN': {'amine_primary': 1},
    # dimethylamine
    'CNC': {'amine_secondary': 1},
    # trimethylamine
    'CN(C)C': {'amine_tertiary': 1},
    # N-benzyl-dimethylamine
    'CN(C)Cc1ccccc1': {'aromatic_ring': 1, 'amine_tertiary': 1},
    # phenylamine
    'Nc1ccccc1': {'aromatic_ring': 1, 'amine_aromatic': 1},
    # N-methyl-N-phenylamine
    'CNc1ccccc1': {'aromatic_ring': 1, 'amine_aromatic': 1},
    # N,N-dimethyl-N-phenylamine
    'CN(C)c1ccccc1': {'aromatic_ring': 1, 'amine_aromatic': 1},
    # peroxy propanyl nitrate
    'CCC(=O)OO[N+](=O)[O-]': {'carbonylperoxynitrate': 1},
    # N-propyl-N-butyl peroxide
    'CCCCOOCCC': {'peroxide': 1},
    # N-butyl hydroperoxide
    'CCCCOO': {'hydroperoxide': 1},
    # peroxyacetic acid
    'CC(=O)OO': {'carbonylperoxyacid': 1},
}

@pytest.mark.parametrize('smiles,expected', paper_examples.items())
def test_simpol_groups_paper_examples(smiles, expected):
    mol = Chem.MolFromSmiles(smiles)
    counts = cu._simpol_groups(mol)
    assert counts['zeroeth'] == 1
    assert counts['carbon'] == sum(a.GetSymbol() == 'C' for a in mol.GetAtoms())

    nonzero = {k: v for k, v in counts.items()
        if v != 0 and k not in ('zeroeth', 'carbon')
    }
    assert nonzero == expected


def test_simpol_paper_log10_p():
    # Examples from section 3.3 of Pankow & Asher (2008), computed with rounded b_k
    # values from their Table 6.
    for smiles, expected in [('C1=CCCCC1', -0.94), ('OC(=O)CCCCC(O)=O', -7.99)]:
        counts = cu._simpol_groups(Chem.MolFromSmiles(smiles))
        assert cu._simpol_log10_p_atm(counts, 293.15) == pytest.approx(expected,
            abs=0.02
        )


# Vapor pressures at 25 C [kPa] of some common odorants (not in classes with known
# biases, listed in `inchi2est_vapor_pressure` doc), from PubChem values stated at 25 C
# (see scripts/validate_est_vapor_pressure.py for comparison on all compounds).
odorant_vapor_pressures_25c = {
    # ethyl acetate (CID 8857)
    'CCOC(C)=O': 12.4,
    # methyl acetate (CID 6584)
    'COC(C)=O': 28.8,
    # butyl acetate (CID 31272)
    'CCCCOC(C)=O': 1.53,
    # isoamyl acetate (CID 31276)
    'CC(=O)OCCC(C)C': 0.747,
    # hexyl acetate (CID 8908)
    'CCCCCCOC(C)=O': 0.176,
    # ethyl propionate (CID 7749)
    'CCOC(=O)CC': 4.77,
    # hexanal (CID 6184)
    'CCCCCC=O': 1.51,
    # decanal (CID 8175)
    'CCCCCCCCCC=O': 0.0137,
    # benzaldehyde (CID 240)
    'O=Cc1ccccc1': 0.169,
    # acetone (CID 180)
    'CC(C)=O': 30.8,
    # 2-butanone (CID 6569)
    'CCC(C)=O': 12.1,
    # 2-heptanone (CID 8051)
    'CCCCCC(C)=O': 0.32,
    # acetophenone (CID 7410)
    'CC(=O)c1ccccc1': 0.0529,
    # ethanol (CID 702)
    'CCO': 6.67,
    # 1-propanol (CID 1031)
    'CCCO': 2.8,
    # 1-pentanol (CID 6276)
    'CCCCCO': 0.293,
    # 1-hexanol (CID 8103)
    'CCCCCCO': 0.124,
    # 2-phenylethanol (CID 6054)
    'OCCc1ccccc1': 0.0116,
    # citronellol (CID 8842)
    'CC(C)=CCCC(C)CCO': 0.00267,
    # p-cresol (CID 2879)
    'Cc1ccc(O)cc1': 0.0148,
    # toluene (CID 1140)
    'Cc1ccccc1': 3.79,
    # alpha-pinene (CID 6654)
    'CC1=CCC2CC1C2(C)C': 0.633,
    # myrcene (CID 31253)
    'C=CC(=C)CCC=C(C)C': 0.279,
    # camphor (CID 2537)
    'CC12CCC(CC1=O)C2(C)C': 0.0867,
}

@pytest.mark.parametrize('smiles,ref_kpa', odorant_vapor_pressures_25c.items())
def test_est_vapor_pressure_odorants(smiles, ref_kpa):
    est = cu.inchi2est_vapor_pressure(_inchi(smiles))
    assert est['vapor_pressure_units'] == 'kPa'
    assert ref_kpa / 3 < est['vapor_pressure'] < ref_kpa * 3


def test_est_vapor_pressure_inputs():
    inchi = _inchi('CCCCCCO')
    est = cu.inchi2est_vapor_pressure(inchi)
    assert est.equals(cu.inchi2est_vapor_pressure(cu.normalize_inchi(inchi)))
    assert est['vapor_pressure_sources'].startswith('SIMPOL.1')

    p20 = cu.inchi2est_vapor_pressure(inchi, temp_c=20)['vapor_pressure']
    assert p20 < est['vapor_pressure']

    assert cu.inchi2est_vapor_pressure(None).isnull().all()


@pytest.mark.parametrize('smiles', [
    # 3-methylthiopropanol (S not supported)
    'CSCCCO',
    # 2-chloroethanol
    'OCCCl',
    # pyridine
    'c1ccncc1',
    # benzonitrile
    'N#Cc1ccccc1',
    # sodium acetate
    'CC(=O)[O-].[Na+]',
    # (excluded classes below)
    # gamma-butyrolactone
    'O=C1CCCO1',
    # gamma-decalactone
    'CCCCCCC1CCC(=O)O1',
    # alanine
    'CC(N)C(=O)O',
])
def test_est_vapor_pressure_unsupported(smiles):
    assert cu.inchi2est_vapor_pressure(_inchi(smiles)).isnull().all()
