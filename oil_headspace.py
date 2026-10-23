#!/usr/bin/env python3
r"""
oil_headspace.py - predicted equilibrium headspace concentrations of dilute
odorants dissolved in paraffin (mineral) oil, in a closed vial.

Method
------
1. Resolve each odorant from a PubChem CID or an InChI. The `chemicals`
   database bundled with thermo is used offline; PubChem is queried only for
   compounds that database lacks.
2. Assign modified-UNIFAC (Dortmund) groups: DDBST's published assignments
   (bundled with thermo, ~29,000 compounds) first, then ugropy + RDKit, or
   your own via the `groups` input column.
3. Compute the infinite-dilution activity coefficient gamma_inf in an n-alkane
   surrogate for the oil (n-hexadecane by default; squalane is also reported
   as a sensitivity check).
4. Convert to the oil/air partition coefficient
       K = C_oil / C_air = c_solv * R * T / (gamma_inf * Psat)
   with c_solv = rho/M of the surrogate [mol/m3]. Because K is volume-based,
   most of the solvent-size effect cancels; hexadecane and squalane typically
   differ by ~0.15-0.2 log units (factor ~1.5), and both are reported so you
   can see that spread for your compounds.
5. Solve the closed-vial mass balance for each dilution:
       C_air = C0 * V_liq / (K * V_liq + V_gas)

Limitations
-----------
* gamma_inf describes isolated solute molecules. Carboxylic acids (which
  dimerize) and, more weakly, alcohols self-associate in oil, so at higher
  concentrations the real headspace falls below the prediction. UNIFAC also
  treats carboxylic acids crudely: regard their K as uncertain.
* For simple functional groups, modified UNIFAC is typically within about a
  factor of 2 in K. If you have an Abraham L (e.g. from the UFZ-LSER database)
  or a measured K, supply it: by default it takes precedence over UNIFAC.
* Equilibrium only: no wall adsorption, no losses through the cap, no
  dynamics of a flowing or repeatedly sampled headspace.

Requirements
------------
    pip install thermo            # required (brings chemicals, numpy, scipy)
    pip install rdkit ugropy      # optional: groups for compounds DDBST lacks

Examples
--------
    python oil_headspace.py 702 8857 240 7362 7991
    python oil_headspace.py "InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3" --T 22
    python oil_headspace.py --input odorants.tsv --out results.csv

Input file
----------
Tab-separated is easiest, because InChIs contain commas (a CSV works if those
fields are quoted, as Excel does). Columns, only `id` required:
    id            PubChem CID or InChI
    name          label used in the output
    psat_pa       vapor pressure at the run temperature [Pa]; overrides thermo.
                  For a solid odorant give the subcooled-liquid value.
    density_g_ml  density of the neat liquid odorant [g/mL] (for v/v)
    abraham_L     Abraham L descriptor = log10 K(hexadecane/air) at 25 C
    logK          measured or other log10 K(oil/air) at the run temperature
    groups        Dortmund subgroups "id:count id:count", e.g. "1:1 2:1 14:1"
A file without an `id` header row is read as one identifier per line.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
import urllib.parse
import urllib.request
import warnings
from dataclasses import dataclass, field

try:
    from chemicals.identifiers import search_chemical
    from thermo import Chemical
    from thermo import unifac as tu
    from thermo.unifac import DOUFSG, UNIFAC
except ImportError:  # pragma: no cover
    sys.exit("This script needs the 'thermo' package:  pip install thermo")

__version__ = "1.0"

R = 8.314462618  # J/(mol K)
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
SOLVENTS = {"hexadecane": "544-76-3", "squalane": "111-01-3"}  # CAS numbers
DEFAULT_DILUTIONS = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
UNITS = {"vv": "v/v", "wv": "g/mL", "ww": "w/w", "molL": "mol/L"}
NONDILUTE_ACTIVITY = 0.1   # flag when gamma*x exceeds this
DEPLETION_FLAG = 0.05      # flag when >5% of the solute sits in the headspace
WALDEN_DSFUS = 56.5        # J/(mol K), Walden's rule for a missing heat of fusion
ESTIMATED_PSAT = {"BOILING_CRITICAL", "LEE_KESLER_PSAT", "AMBROSE_WALTON",
                  "SANJARI", "EDALAT", "EOS"}  # thermo methods that are estimates
ESTIMATED_TM = {"JOBACK"}


# --------------------------------------------------------------------------
# Data containers
# --------------------------------------------------------------------------
@dataclass
class Solvent:
    name: str
    groups: dict
    c: float      # molar density [mol/m3]
    Tm: float | None


@dataclass
class Odorant:
    query: str
    name: str = ""
    cid: int | None = None
    cas: str | None = None
    inchikey: str | None = None
    smiles: str | None = None
    mw: float | None = None            # g/mol
    groups: dict | None = None
    groups_source: str = ""
    params: str = ""
    psat: float | None = None          # Pa
    psat_source: str = ""
    rho: float | None = None           # kg/m3, neat liquid
    rho_source: str = ""
    Tm: float | None = None
    Hfus: float | None = None          # J/mol
    abraham_L: float | None = None
    logK_input: float | None = None
    gamma: dict = field(default_factory=dict)        # solvent -> gamma_inf
    logK_unifac: dict = field(default_factory=dict)  # solvent -> log10 K
    logK_used: float | None = None
    k_source: str = ""
    gamma_eff: float | None = None     # gamma consistent with the K used
    Csat: float | None = None          # predicted solid solubility [mol/m3]
    p_solid: float | None = None       # predicted vapor pressure of the solid [Pa]
    notes: list = field(default_factory=list)
    error: str = ""


# --------------------------------------------------------------------------
# Identifier resolution
# --------------------------------------------------------------------------
def parse_id(text):
    """Return ('inchi', str) or ('cid', int)."""
    s = text.strip()
    if s.startswith("InChI="):
        return "inchi", s
    m = re.fullmatch(r"(?i)(?:pubchem|cid)?[\s:=#]*(\d+)", s)
    if m:
        return "cid", int(m.group(1))
    raise ValueError("not a PubChem CID or an InChI starting with 'InChI='")


def offline_search(query):
    try:
        return search_chemical(query)
    except Exception:
        return None


def rdkit_from_inchi(inchi):
    """(inchikey, smiles, mw) from an InChI via RDKit, or None if unavailable."""
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import Descriptors
    except ImportError:
        return None
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        return None
    return Chem.MolToInchiKey(mol), Chem.MolToSmiles(mol), Descriptors.MolWt(mol)


def pubchem_properties(kind, value, timeout=20):
    """Basic identifiers for a CID or an InChI from PubChem PUG-REST."""
    props = "InChI,InChIKey,MolecularWeight,Title"
    if kind == "cid":
        url, data = f"{PUG}/compound/cid/{value}/property/{props}/JSON", None
    else:
        url = f"{PUG}/compound/inchi/property/{props}/JSON"
        data = urllib.parse.urlencode({"inchi": value}).encode()
    req = urllib.request.Request(
        url, data=data, headers={"User-Agent": f"oil_headspace.py/{__version__}"})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        payload = json.load(resp)
    p = payload["PropertyTable"]["Properties"][0]
    mw = p.get("MolecularWeight")
    return {"cid": p.get("CID"), "inchi": p.get("InChI"),
            "inchikey": p.get("InChIKey"), "title": p.get("Title"),
            "mw": float(mw) if mw not in (None, "") else None}


def resolve(od, offline=False):
    kind, value = parse_id(od.query)
    inchi = value if kind == "inchi" else None
    if kind == "cid":
        od.cid = value
        meta = offline_search(f"pubchem={value}")
    else:
        meta = offline_search(value)

    if meta is None:
        if kind == "inchi":
            rd = rdkit_from_inchi(value)
            if rd:
                od.inchikey, od.smiles, od.mw = rd
        if not od.inchikey and not offline:
            try:
                info = pubchem_properties(kind, value)
            except Exception as exc:
                raise LookupError(f"not in the offline database and the PubChem "
                                  f"lookup failed ({exc})") from None
            od.cid = od.cid or info["cid"]
            od.inchikey = info["inchikey"]
            od.mw = od.mw or info["mw"]
            od.name = od.name or (info["title"] or "")
            inchi = inchi or info["inchi"]
            if inchi and not od.smiles:
                rd = rdkit_from_inchi(inchi)
                if rd:
                    od.smiles = rd[1]
        if not od.inchikey:
            raise LookupError("not in the offline database; run without --offline "
                              "or install RDKit")
        meta = offline_search("InChIKey=" + od.inchikey)

    exact = meta is not None
    if not exact and od.inchikey:
        # same skeleton without stereo/isotope layers (e.g. E/Z or R/S forms)
        generic = od.inchikey.split("-")[0] + "-UHFFFAOYSA-N"
        if generic != od.inchikey:
            meta = offline_search("InChIKey=" + generic)

    if meta is not None:
        od.cas = meta.CASs
        od.mw = meta.MW
        od.smiles = od.smiles or meta.smiles
        od.name = od.name or meta.common_name or ""
        if exact:
            od.inchikey = meta.InChI_key
            od.cid = od.cid or meta.pubchemid
        else:
            od.notes.append(f"properties taken from {meta.common_name or meta.CASs} "
                            f"(same skeleton; stereo/isotope not distinguished)")
    else:
        od.notes.append("not in thermo's property database: supply psat_pa "
                        "(and density_g_ml for v/v)")
    od.name = od.name or (f"CID {od.cid}" if od.cid else od.query[:30])


# --------------------------------------------------------------------------
# UNIFAC groups and activity coefficients
# --------------------------------------------------------------------------
_DDBST = None


def ddbst_table():
    global _DDBST
    if _DDBST is None:
        tu.load_group_assignments_DDBST()
        _DDBST = tu.DDBST_MODIFIED_UNIFAC_assignments
    return _DDBST


def ddbst_groups(inchikey):
    table = ddbst_table()
    if inchikey in table:
        return dict(table[inchikey]), "DDBST"
    block = inchikey.split("-")[0] + "-"
    matches = sorted(k for k in table if k.startswith(block))
    if matches:
        best = next((k for k in matches if k.endswith("-UHFFFAOYSA-N")), matches[0])
        return dict(table[best]), f"DDBST, via {best} (same skeleton)"
    return None, ""


def ugropy_groups(smiles):
    try:
        from ugropy import dortmund
    except ImportError:
        return None, "not in DDBST table; install rdkit + ugropy for automatic groups"
    try:
        result = dortmund.get_groups(smiles, "smiles")
    except Exception as exc:
        return None, f"ugropy failed ({exc})"
    groups = getattr(result, "subgroups_num", None) or {}
    if not groups:
        return None, "ugropy found no complete Dortmund fragmentation"
    return {int(k): int(v) for k, v in groups.items()}, "ugropy (automatic, check it)"


def parse_groups(text):
    pairs = re.findall(r"(\d+)\s*[:*x]\s*(\d+)", text)
    if not pairs:
        raise ValueError(f"cannot parse groups {text!r}; use 'id:count id:count'")
    return {int(k): int(n) for k, n in pairs}


def groups_text(groups, ids=False):
    return " ".join(f"{k if ids else DOUFSG[k].group}:{n}"
                    for k, n in sorted(groups.items()))


def assign_groups(od):
    if od.groups:                                   # user override
        od.groups_source = "input"
    else:
        reason = ""
        if od.inchikey:
            od.groups, od.groups_source = ddbst_groups(od.inchikey)
        if not od.groups and od.smiles:
            od.groups, reason = ugropy_groups(od.smiles)
            od.groups_source = reason if od.groups else ""
        if not od.groups:
            od.notes.append(f"no UNIFAC groups: {reason or 'no structure available'}")
            return
    unknown = [k for k in od.groups if k not in DOUFSG]
    if unknown:
        od.notes.append(f"unknown Dortmund subgroup ids {unknown}")
        od.groups = None
    elif od.groups_source != "DDBST":
        od.notes.append(f"groups {groups_text(od.groups)} from {od.groups_source}")


def parameter_sets():
    sets = {"dortmund2016": tu.DOUFIP2016}
    try:
        sets["unifac2"] = tu.DOUF2IP   # modified UNIFAC 2.0 (newer thermo only)
    except AttributeError:
        pass
    return sets


def missing_pairs(group_dicts, ip):
    mains = sorted({DOUFSG[k].main_group_id for g in group_dicts for k in g})
    return [(a, b) for a in mains for b in mains
            if a != b and b not in ip.get(a, {})]


def gamma_inf(g_solute, g_solvent, T, ip):
    model = UNIFAC.from_subgroups(T=T, xs=[0.5, 0.5],
                                  chemgroups=[g_solute, g_solvent], version=1,
                                  interaction_data=ip, subgroups=DOUFSG)
    return float(model.gammas_infinite_dilution()[0])


# --------------------------------------------------------------------------
# Pure-component properties
# --------------------------------------------------------------------------
def pure_properties(od, T, P):
    if not od.cas:
        return
    try:
        ch = Chemical(od.cas, T=T, P=P)
    except Exception as exc:
        od.notes.append(f"thermo could not load properties ({exc})")
        return
    od.Tm = ch.Tm
    if od.Tm and getattr(ch, "Tm_source", None) in ESTIMATED_TM:
        if od.Tm > T:
            od.notes.append(f"melting point only estimated ({od.Tm - 273.15:.0f} C, "
                            f"{ch.Tm_source}); treated as a liquid")
        od.Tm = None
    od.Hfus = ch.Hfusm
    if od.mw is None:
        od.mw = ch.MW
    solid = bool(od.Tm and T < od.Tm)
    if od.psat is None:
        try:
            psat = ch.VaporPressure(T)
        except Exception:
            psat = None
        if psat and psat > 0:
            method = ch.VaporPressure.method
            od.psat = psat
            od.psat_source = f"thermo {method}"
            if method in ESTIMATED_PSAT:
                od.notes.append(f"Psat is estimated ({method}), not measured; "
                                f"supply psat_pa from a reference if you can")
            if solid:
                od.notes.append("solid at T: Psat is the extrapolated "
                                "subcooled-liquid value")
    elif solid:
        od.notes.append("solid at T: psat_pa should be the subcooled-liquid value")
    if od.rho is None:
        try:
            vm = ch.VolumeLiquid(T, P)
        except Exception:
            vm = None
        if vm:
            od.rho = ch.MW / 1000.0 / vm
            od.rho_source = f"thermo {ch.VolumeLiquid.method}"


def load_solvent(name, T, P):
    meta = search_chemical(SOLVENTS[name])
    groups = dict(ddbst_table()[meta.InChI_key])
    ch = Chemical(SOLVENTS[name], T=T, P=P)
    return Solvent(name=name, groups=groups, c=1.0 / ch.VolumeLiquid(T, P), Tm=ch.Tm)


# --------------------------------------------------------------------------
# Partitioning and headspace
# --------------------------------------------------------------------------
def analyze(od, args, solvents, psets):
    T = args.T_K
    s = solvents[args.solvent]
    if od.groups:
        order = [args.params] + [p for p in psets if p != args.params]
        for pname in order:
            if not any(missing_pairs([od.groups, sv.groups], psets[pname])
                       for sv in solvents.values()):
                od.params = pname
                break
        if not od.params:
            pairs = missing_pairs([od.groups, s.groups], psets[args.params])
            od.notes.append(f"no Dortmund interaction parameters for main-group "
                            f"pairs {pairs}")
        else:
            if od.params != args.params:
                od.notes.append(f"used {od.params} parameters because "
                                f"{args.params} lacks a group pair; these differ "
                                f"noticeably, so treat this K with extra caution")
            for sname, sv in solvents.items():
                g = gamma_inf(od.groups, sv.groups, T, psets[od.params])
                od.gamma[sname] = g
                if od.psat:
                    od.logK_unifac[sname] = math.log10(sv.c * R * T / (g * od.psat))
        mains = {DOUFSG[k].main_group for k in od.groups}
        if mains & {"COOH", "HCOOH"}:
            od.notes.append("carboxylic acid: UNIFAC ignores acid dimerization; "
                            "K is uncertain, compare with Abraham L or measure")
        elif mains & {"OH", "CH3OH", "ACOH"}:
            od.notes.append("OH group: self-association lowers the headspace at "
                            "higher concentrations")

    auto = args.k_source == "auto"
    if auto and od.logK_input is not None:
        od.logK_used, od.k_source = od.logK_input, "input logK"
    elif auto and od.abraham_L is not None:
        od.logK_used, od.k_source = od.abraham_L, "Abraham L"
        if abs(T - 298.15) > 0.5:
            od.notes.append("Abraham L refers to 25 C; no temperature correction")
    elif args.solvent in od.logK_unifac:
        od.logK_used, od.k_source = od.logK_unifac[args.solvent], "UNIFAC"
    else:
        if od.gamma:
            missing, fix = "no vapor pressure", "psat_pa, abraham_L or logK"
        elif not od.groups:
            missing, fix = "no UNIFAC groups", "groups, abraham_L or logK"
        else:
            missing, fix = "no UNIFAC parameters", "abraham_L or logK"
        od.error = f"no partition coefficient ({missing}); supply {fix}"
        return

    K = 10.0 ** od.logK_used
    if od.psat:
        od.gamma_eff = s.c * R * T / (K * od.psat)
        if od.Tm and T < od.Tm:
            hfus = od.Hfus
            if not hfus:
                hfus = WALDEN_DSFUS * od.Tm
                od.notes.append("heat of fusion estimated by Walden's rule")
            a_sat = math.exp(-hfus / R * (1.0 / T - 1.0 / od.Tm))
            od.Csat = a_sat / od.gamma_eff * s.c
            od.p_solid = a_sat * od.psat     # approx. vapor pressure of the solid


def oil_concentration(value, unit, od, rho_oil):
    """Initial concentration in the oil [mol/m3], or None if not computable."""
    if unit == "molL":
        return value * 1000.0
    if od.mw is None:
        return None
    if unit == "vv":
        return value * od.rho / (od.mw / 1000.0) if od.rho else None
    if unit == "wv":                   # g per mL of solution
        return value * 1e6 / od.mw
    if unit == "ww":                   # g per g of solution
        return value * rho_oil * 1000.0 / od.mw
    raise ValueError(unit)


def headspace_rows(od, args, solvents):
    T, s = args.T_K, solvents[args.solvent]
    K = 10.0 ** od.logK_used
    v_liq, v_gas = args.v_liquid, args.v_vial - args.v_liquid
    if not args.no_depletion:
        frac = v_gas / (K * v_liq + v_gas)   # independent of concentration
        if frac > DEPLETION_FLAG:
            od.notes.append(f"{frac:.0%} of the solute moves into the headspace "
                            f"(included in the mass balance)")
    rows = []
    for d in args.dilutions:
        row = {"dilution": d, "unit": UNITS[args.unit]}
        flags = []
        c0 = oil_concentration(d, args.unit, od, args.oil_density * 1000.0)
        if c0 is None:
            row["flags"] = "cannot convert concentration (need density_g_ml / MW)"
            rows.append(row)
            continue
        if args.no_depletion:
            c_gas, frac = c0 / K, None
        else:
            c_gas = c0 * v_liq / (K * v_liq + v_gas)
            frac = c_gas * v_gas / (c0 * v_liq)
        c_oil = K * c_gas
        x = c_oil / (s.c + c_oil)
        activity = od.gamma_eff * x if od.gamma_eff else None
        if activity is not None and not od.Csat:
            if activity >= 1.0:
                flags.append("gamma*x >= 1: may exceed miscibility")
            elif activity > NONDILUTE_ACTIVITY:
                flags.append(f"gamma*x = {activity:.2f}: not fully dilute")
        if od.Csat and c0 > od.Csat:
            cap = od.p_solid / args.pressure * 1e6
            flags.append(f"above predicted solid solubility; a saturated "
                         f"headspace would be ~{cap:.2g} ppmv")
        p = c_gas * R * T
        row.update({
            "C_oil_initial_mol_L": c0 / 1000.0,
            "x_oil": x,
            "activity": activity,
            "frac_in_headspace": frac,
            "p_Pa": p,
            "ppmv": p / args.pressure * 1e6,
            "C_gas_mg_m3": c_gas * od.mw * 1000.0 if od.mw else None,
            "C_gas_mol_m3": c_gas,
            "flags": "; ".join(flags),
        })
        rows.append(row)
    return rows


# --------------------------------------------------------------------------
# Input / output
# --------------------------------------------------------------------------
def _num(text):
    text = (text or "").strip()
    return float(text) if text else None


def read_input_file(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        lines = [ln for ln in f.read().splitlines()
                 if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return []
    delim = "\t" if ("\t" in lines[0] or path.lower().endswith(".tsv")) else ","
    reader = csv.reader(lines, delimiter=delim)
    header = [h.strip().lower() for h in next(reader)]
    if "id" not in header:
        return [Odorant(query=ln.strip()) for ln in lines]
    aliases = {"l": "abraham_l", "logk_oil_air": "logk", "psat": "psat_pa",
               "density": "density_g_ml", "cid": "id", "inchi": "id"}
    header = [aliases.get(h, h) for h in header]
    odorants = []
    for n, row in enumerate(reader, start=2):
        rec = {h: (v.strip() if v else "") for h, v in zip(header, row)}
        if not rec.get("id"):
            continue
        od = Odorant(query=rec["id"], name=rec.get("name", ""))
        try:
            od.psat = _num(rec.get("psat_pa"))
            if od.psat:
                od.psat_source = "input"
            rho = _num(rec.get("density_g_ml"))
            if rho:
                od.rho, od.rho_source = rho * 1000.0, "input"
            od.abraham_L = _num(rec.get("abraham_l"))
            od.logK_input = _num(rec.get("logk"))
            if rec.get("groups"):
                od.groups = parse_groups(rec["groups"])
        except ValueError as exc:
            od.error = f"input line {n}: {exc}"
        odorants.append(od)
    return odorants


SUMMARY_FIELDS = [
    "query", "name", "cid", "cas", "inchikey", "dortmund_groups",
    "dortmund_group_ids", "groups_source", "params", "T_C", "solvent_surrogate",
    "psat_Pa", "psat_source", "density_g_ml", "density_source",
    "gamma_inf_hexadecane", "gamma_inf_squalane", "logK_unifac_hexadecane",
    "logK_unifac_squalane", "abraham_L", "logK_input", "logK_used", "k_source",
    "H_px_Pa", "Hcp_mol_m3_Pa", "Csat_solid_mol_L", "p_solid_Pa", "notes"]
ROW_FIELDS = ["dilution", "unit", "C_oil_initial_mol_L", "x_oil", "activity",
              "frac_in_headspace", "p_Pa", "ppmv", "C_gas_mg_m3", "C_gas_mol_m3",
              "flags", "error"]


def summary_record(od, args):
    T = args.T_K
    K = 10.0 ** od.logK_used if od.logK_used is not None else None
    return {
        "query": od.query, "name": od.name, "cid": od.cid, "cas": od.cas,
        "inchikey": od.inchikey,
        "dortmund_groups": groups_text(od.groups) if od.groups else "",
        "dortmund_group_ids": groups_text(od.groups, ids=True) if od.groups else "",
        "groups_source": od.groups_source, "params": od.params,
        "T_C": args.T, "solvent_surrogate": args.solvent,
        "psat_Pa": od.psat, "psat_source": od.psat_source,
        "density_g_ml": od.rho / 1000.0 if od.rho else None,
        "density_source": od.rho_source,
        "gamma_inf_hexadecane": od.gamma.get("hexadecane"),
        "gamma_inf_squalane": od.gamma.get("squalane"),
        "logK_unifac_hexadecane": od.logK_unifac.get("hexadecane"),
        "logK_unifac_squalane": od.logK_unifac.get("squalane"),
        "abraham_L": od.abraham_L, "logK_input": od.logK_input,
        "logK_used": od.logK_used, "k_source": od.k_source,
        "H_px_Pa": od.gamma_eff * od.psat if od.gamma_eff and od.psat else None,
        "Hcp_mol_m3_Pa": K / (R * T) if K else None,
        "Csat_solid_mol_L": od.Csat / 1000.0 if od.Csat else None,
        "p_solid_Pa": od.p_solid,
        "notes": "; ".join(od.notes),
    }


def fmt(value, spec):
    if value is None or value == "":
        return "-"
    return format(value, spec)


def print_report(odorants, args, rows_by_od):
    v_gas = args.v_vial - args.v_liquid
    print(f"\nOil/air partitioning at {args.T:g} C, modified UNIFAC (Dortmund, "
          f"{args.params} parameters); log K = log10(C_oil/C_air)")
    print(f"{'name':<22}{'CID':>10}{'gamma_inf':>10}{'Psat [Pa]':>10}"
          f"{'C16':>7}{'C30':>7}{'used':>7}  K from")
    for od in odorants:
        if od.error and od.logK_used is None:
            print(f"{od.name[:21]:<22}{fmt(od.cid, 'd'):>10}  ERROR: {od.error}")
            continue
        print(f"{od.name[:21]:<22}{fmt(od.cid, 'd'):>10}"
              f"{fmt(od.gamma.get(args.solvent), '.3g'):>10}"
              f"{fmt(od.psat, '.3g'):>10}"
              f"{fmt(od.logK_unifac.get('hexadecane'), '.2f'):>7}"
              f"{fmt(od.logK_unifac.get('squalane'), '.2f'):>7}"
              f"{fmt(od.logK_used, '.2f'):>7}  {od.k_source}")
    print(f"(C16 = hexadecane, C30 = squalane surrogate; gamma_inf in {args.solvent}; "
          f"'used' drives the headspace table)")

    geometry = ("no depletion (oil as infinite reservoir)" if args.no_depletion else
                f"closed vial, {args.v_liquid:g} mL oil + {v_gas:g} mL headspace")
    print(f"\nEquilibrium headspace [ppmv]; {geometry}; dilutions in {UNITS[args.unit]}")
    print(f"{'name':<22}" + "".join(f"{d:>10.0e}" for d in args.dilutions))
    for od in odorants:
        rows = rows_by_od.get(id(od))
        if not rows:
            continue
        cells = "".join(f"{fmt(r.get('ppmv'), '.2e'):>10}" for r in rows)
        marks = "*" if any(r.get("flags") for r in rows) else ""
        print(f"{od.name[:21]:<22}{cells} {marks}")

    print("\nNotes")
    for od in odorants:
        row_flags = []
        for r in rows_by_od.get(id(od), []):
            if r.get("flags"):
                row_flags.append(f"{r['dilution']:.0e}: {r['flags']}")
        items = od.notes + row_flags
        if items:
            print(f"  {od.name}:")
            for item in items:
                print(f"    - {item}")
    print(f"\n(* = see notes)  Full results: {args.out}")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Predict equilibrium headspace concentrations of dilute "
                    "odorants in paraffin oil (modified UNIFAC, alkane surrogate).",
        epilog="Example: python oil_headspace.py 702 8857 240 7362 7991",
    )
    ap.add_argument("ids", nargs="*",
                    help="PubChem CIDs and/or InChI strings (quote InChIs)")
    ap.add_argument("-i", "--input",
                    help="TSV/CSV with an 'id' column, or one identifier per line")
    ap.add_argument("--dilutions", nargs="+", type=float, default=DEFAULT_DILUTIONS,
                    help="concentrations in oil (default: 1e-3 ... 1e-8)")
    ap.add_argument("--unit", choices=list(UNITS), default="vv",
                    help="vv = volume of neat odorant per volume of solution "
                         "(default); wv = g/mL; ww = g/g; molL = mol/L")
    ap.add_argument("--T", type=float, default=25.0, help="temperature [C]")
    ap.add_argument("--v-liquid", type=float, default=2.0, help="oil volume [mL]")
    ap.add_argument("--v-vial", type=float, default=10.0,
                    help="total vial volume [mL]")
    ap.add_argument("--no-depletion", action="store_true",
                    help="ignore depletion of the oil (C_air = C0/K)")
    ap.add_argument("--solvent", choices=list(SOLVENTS), default="hexadecane",
                    help="alkane surrogate for the oil (default hexadecane)")
    ap.add_argument("--params", choices=["dortmund2016", "unifac2"],
                    default="dortmund2016",
                    help="interaction parameters; the other set is used "
                         "automatically if a group pair is missing")
    ap.add_argument("--k-source", choices=["auto", "unifac"], default="auto",
                    help="auto: input logK > abraham_L > UNIFAC; unifac: always UNIFAC")
    ap.add_argument("--oil-density", type=float, default=0.85,
                    help="paraffin oil density [g/mL], used only for --unit ww")
    ap.add_argument("--pressure", type=float, default=101325.0,
                    help="ambient pressure for ppmv [Pa]")
    ap.add_argument("--offline", action="store_true", help="never query PubChem")
    ap.add_argument("-o", "--out", default="headspace_predictions.csv",
                    help="output CSV (one row per odorant and dilution)")
    args = ap.parse_args(argv)

    if args.v_vial <= args.v_liquid:
        ap.error("--v-vial must be larger than --v-liquid")
    if any(d <= 0 for d in args.dilutions):
        ap.error("dilutions must be positive")
    args.T_K = args.T + 273.15

    odorants = [Odorant(query=q) for q in args.ids]
    if args.input:
        odorants += read_input_file(args.input)
    if not odorants:
        ap.error("give CIDs/InChIs on the command line or an --input file")

    warnings.filterwarnings("ignore", category=RuntimeWarning)
    psets = parameter_sets()
    if args.params not in psets:
        ap.error(f"parameter set {args.params} needs a newer thermo version")
    P = args.pressure
    solvents = {name: load_solvent(name, args.T_K, P) for name in SOLVENTS}
    if solvents[args.solvent].Tm and args.T_K < solvents[args.solvent].Tm:
        print(f"Note: {args.solvent} is below its melting point at {args.T:g} C; "
              f"it is used as a hypothetical liquid (consider --solvent squalane).")

    rows_by_od = {}
    for od in odorants:
        if not od.error:
            try:
                resolve(od, offline=args.offline)
                pure_properties(od, args.T_K, P)
                assign_groups(od)
                analyze(od, args, solvents, psets)
            except Exception as exc:
                od.error = str(exc)
        if not od.name:
            od.name = od.query[:30]
        if od.logK_used is not None:
            rows_by_od[id(od)] = headspace_rows(od, args, solvents)

    with open(args.out, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=SUMMARY_FIELDS + ROW_FIELDS,
                                extrasaction="ignore")
        writer.writeheader()
        for od in odorants:
            base = summary_record(od, args)
            rows = rows_by_od.get(id(od)) or [{"error": od.error}]
            for row in rows:
                writer.writerow({**base, **row, "error": od.error})

    print_report(odorants, args, rows_by_od)
    return 0


if __name__ == "__main__":
    sys.exit(main())
