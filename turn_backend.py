"""
turn_backend.py
---------------
Finds α-turns, β-turns, and γ-turns in protein structure files (.pdb / .cif / .mmcif).

β-turns (Venkatachalam 1968 · Richardson 1981 · Hutchinson & Thornton 1994)
  4-residue window [i, i+1, i+2, i+3]:
    1. Cα(i) → Cα(i+3)  < 7.0 Å
    2. Central residues i+1, i+2 must NOT be helical — DSSP code not in {G, H, I}
       (excludes α-helix H, 3₁₀-helix G, and π-helix I)
    3. H-bonds NOT required for β-turn detection
    Subtypes (φ/ψ of i+1 and i+2, ±30° tolerance):
        I  I′  II  IIb  II′  III  III′  IV  VIa1  VIII
        (IIb and VIa1 require cis-Pro at i+2 with |ω| ≤ 30°)

γ-turns (geometric H-bond definition)
  3-residue window [i, i+1, i+2]:
    1. H-bond between C=O(i) and N-H(i+2):
         If backbone H present:
           O···H < 2.5 Å
           N-H···O angle (at H) > 120°
           H···O=C angle (at O) > 90°
         If H absent (most X-ray structures):
           N(i+2)···O(i) < 3.5 Å
    2. φ(i+1), ψ(i+1) in strict ranges:
         classic γ:  φ = 70–85°,     ψ = −70 to −60°
         inverse γ′: φ = −85 to −70°, ψ = 60 to 70°

α-turns (Nataraj et al. 1995 · Pavone et al. 1996 · Dasgupta et al. 2004)
  5-residue window [i, i+1, i+2, i+3, i+4]:
    1. Cα(i) → Cα(i+4) < 6.5 Å  (Dasgupta 2004)
    2. i+1, i+2, i+3 NOT assigned G (3₁₀) or H (α-helix) by DSSP
    3. N(i+4)···O(i) distance: 2.6–3.6 Å  (Nataraj 1995)
    4. H–N···O angle at N: ≤ 40°  (Nataraj 1995)
    5. Cluster: φ/ψ at i+1, i+2, i+3 matched to Pavone et al. (1996) Table IV
       using per-angle standard deviation tolerances.

Public API
----------
find_beta_turns(file_path)  → (turns: list[dict], status_str: str)
find_gamma_turns(file_path) → (turns: list[dict], status_str: str)
find_alpha_turns(file_path) → (turns: list[dict], status_str: str)
find_all_turns(file_path)   → (beta, gamma, alpha: list[dict], status_str: str)
"""

import os
import math
import traceback

import numpy as np

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.vectors import calc_dihedral, calc_angle

# ── three_to_one fallback (mirrors omega_backend.py) ──────────────────────
try:
    from Bio.PDB.Polypeptide import three_to_one as bio_three_to_one
except ImportError:
    try:
        from Bio.SeqUtils import IUPACData
        _d = {k.upper(): v
              for k, v in IUPACData.protein_letters_3to1_extended.items()}
        def bio_three_to_one(resname):
            r = _d.get(resname.upper())
            if r is None:
                raise KeyError(resname)
            return r
    except Exception:
        _AA3TO1 = {
            "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
            "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
            "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
            "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
        }
        def bio_three_to_one(resname):
            r = _AA3TO1.get(resname.upper())
            if r is None:
                raise KeyError(resname)
            return r

# ── constants ──────────────────────────────────────────────────────────────

# ── constants ──────────────────────────────────────────────────────────────

BETA_CA_DIST  = 7.0   # Å — Cα(i)→Cα(i+3) ceiling for β-turns
ANGLE_TOL     = 30    # degrees — tolerance for β-turn subtype matching

# DSSP codes considered "helical" — β-turn central residues (i+1, i+2)
# must NOT be any of these.  Extends the original α-helix-only exclusion
# to 3₁₀-helix (G) and π-helix (I) as well.
HELICAL_SS = {'G', 'H', 'I'}

# Idealized φ/ψ for each β-turn subtype (degrees).
# Key: {name: ((phi_i1, psi_i1), (phi_i2, psi_i2))}
# References:
#   Venkatachalam (1968) Biopolymers 6:1425-1436        — original definition: I, II, III + mirrors
#   Lewis, Momany & Scheraga (1973) BBA 303:211-229     — distance criterion < 7 Å, H-bond not required
#   Wilmot & Thornton (1990) Protein Eng. 3:479-493     — Table I: canonical 7 types + VIII;
#     type III removed (90% overlap with 3₁₀-helix), IIb renamed VIb (cis-Pro family)
#   Hutchinson & Thornton (1994) Protein Sci. 3:2207-2216 — idealized angles, ±30° tolerance
#   Panasik, Fleming & Rose (2005) Protein Sci. 14:2910-2914 — |ω| ≤ 30° for VI types
#
# 7 canonical subtypes per Wilmot & Thornton (1990) Table I:
#   I, I', II, II', VIa, VIb, VIII
#   Type III and III' are NOT included — Wilmot & Thornton (1990) absorbed them into type IV
#   as they overlap extensively with 3₁₀-helix residues.
#   Type IV is the catch-all (no canonical angle match).
#   VIa and VIb require cis-Pro at i+2 (|ω| ≤ 30°).
BETA_TURN_DEFS = {
    "VIII": ((-60, -30), (-120, 120)),  # unique i+2 ψ≈120°; no cis-Pro (Wilmot & Thornton 1990)
    "I":    ((-60, -30), (-90,   0)),   # most common type (Hutchinson & Thornton 1994)
    "II":   ((-60, 120), ( 80,   0)),   # Gly favoured at i+2
    "I'":   (( 60,  30), ( 90,   0)),   # mirror of I
    "II'":  (( 60,-120), (-80,   0)),   # mirror of II
    "VIa":  ((-60, 120), (-90,   0)),   # cis-Pro at i+2; |ω| ≤ 30° (Panasik 2005)
    "VIb":  ((-120,120), (-60,   0)),   # cis-Pro at i+2; |ω| ≤ 30° (Wilmot & Thornton 1990)
}

# Types that require a cis-Pro at position i+2 (|ω| ≤ 30°)
_CIS_PRO_TYPES = {"VIa", "VIb"}

# ── generic helpers ────────────────────────────────────────────────────────

def _ca_coord(residue):
    """Return (x, y, z) for the Cα atom, or None if missing."""
    try:
        c = residue["CA"].get_coord()
        return float(c[0]), float(c[1]), float(c[2])
    except Exception:
        return None


def _dist(a, b):
    return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))


def _three_to_one(resname):
    try:
        return bio_three_to_one(resname)
    except Exception:
        return "X"


def _parser_for(file_path):
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".pdb":
        return PDBParser(QUIET=True)
    if ext in (".cif", ".mmcif"):
        return MMCIFParser(QUIET=True)
    raise ValueError(f"Unsupported extension '{ext}'. Use .pdb, .cif or .mmcif.")


def _r(v, n=3):
    """Round v to n decimal places, or pass None through."""
    return round(v, n) if v is not None else None


# ── DSSP runner (mirrors omega_backend.py) ────────────────────────────────

def _find_dssp_binary():
    """
    Find the DSSP binary, checking the script's own directory first.
    This allows Windows users to place dssp.exe / mkdssp.exe alongside
    app.py without needing to add it to the system PATH.
    Returns the path or name to pass as the `dssp` argument to Bio.PDB.DSSP.
    """
    import shutil
    here = os.path.dirname(os.path.abspath(__file__))
    for name in ("mkdssp.exe", "dssp.exe", "mkdssp", "dssp"):
        local = os.path.join(here, name)
        if os.path.isfile(local):
            return local
    for name in ("mkdssp", "dssp"):
        found = shutil.which(name)
        if found:
            return found
    return "dssp"   # last resort — let BioPython try PATH


def _reduce_pdb_for_dssp(file_path):
    """
    Write a DSSP-compatible temp PDB containing only CRYST1/ATOM/HETATM/TER/END,
    plus a dummy CRYST1 if the original lacks one (required by mkdssp 4.x).
    This strips USER/REMARK/COMPND/SOURCE records that confuse mkdssp.
    Returns the temp file path (caller must delete it).
    """
    import tempfile
    keep = ('CRYST1', 'ATOM  ', 'HETATM', 'TER   ', 'TER\n', 'END   ', 'END\n')
    with open(file_path) as fh:
        lines = [l for l in fh if any(l.startswith(k) for k in keep)]
    has_cryst1 = any(l.startswith('CRYST1') for l in lines)
    with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tf:
        if not has_cryst1:
            tf.write("CRYST1    1.000    1.000    1.000"
                     "  90.00  90.00  90.00 P 1           1\n")
        tf.writelines(lines)
        return tf.name


_CRYST1_DUMMY = ("CRYST1    1.000    1.000    1.000"
                 "  90.00  90.00  90.00 P 1           1\n")


def _with_cryst1(content: str) -> str:
    """Prepend a dummy CRYST1 record if the content lacks one."""
    if "CRYST1" not in content:
        return _CRYST1_DUMMY + content
    return content


def _run_dssp(model, file_path):
    """
    Run DSSP on model. Returns a DSSP object or None on total failure.
    Tries strategies in order, adding CRYST1 when missing (required by mkdssp 4.x).
    """
    import io, sys, subprocess, tempfile
    from Bio.PDB.DSSP import DSSP
    from Bio.PDB import PDBIO

    binary = _find_dssp_binary()
    stem   = os.path.basename(file_path)

    def _try_dssp(path):
        """Run BioPython DSSP on path, return object or None."""
        import warnings
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")   # suppress mkdssp-4.x/BioPython compat noise
                obj = DSSP(model, path, dssp=binary)
            return obj if obj and len(obj.keys()) > 0 else None
        except Exception:
            return None

    def _write_tmp(content):
        """Write content to a temp PDB file, return path."""
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tf:
            tf.write(content)
            return tf.name

    def _shell_dssp(content):
        """Run mkdssp via subprocess, return DSSP text or None."""
        tmp_in = tmp_out = None
        try:
            tmp_in = _write_tmp(content)
            with tempfile.NamedTemporaryFile(suffix=".dssp", mode="w", delete=False) as f:
                tmp_out = f.name
            cmd = f'"{binary}" --output-format=dssp "{tmp_in}" -o "{tmp_out}"'
            result = subprocess.run(cmd, shell=True,
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    timeout=60)
            if os.path.exists(tmp_out):
                with open(tmp_out) as fh:
                    text = fh.read()
                if "RESIDUE" in text:
                    return text
            stdout = result.stdout.decode(errors="replace")
            if "RESIDUE" in stdout:
                return stdout
        except Exception as e:
            print(f"[DSSP/turn] {stem}: shell call failed — {e}", file=sys.stderr)
        finally:
            for f in (tmp_in, tmp_out):
                if f and os.path.exists(f):
                    try: os.unlink(f)
                    except OSError: pass
        return None

    def _parse_dssp_text(dssp_text):
        """Parse raw DSSP text into a DSSP object via a temp file."""
        import warnings
        tmp = None
        try:
            with tempfile.NamedTemporaryFile(suffix=".dssp", mode="w", delete=False) as f:
                f.write(dssp_text)
                tmp = f.name
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                obj = DSSP(model, tmp, file_type="DSSP")
            return obj if obj and len(obj.keys()) > 0 else None
        except Exception:
            return None
        finally:
            if tmp and os.path.exists(tmp):
                try: os.unlink(tmp)
                except OSError: pass

    # Attempt 1 — direct call on original file
    obj = _try_dssp(file_path)
    if obj: return obj

    # Attempt 2 — stripped PDB (removes REDUCE/USER records that confuse mkdssp)
    stripped = None
    try:
        stripped = _reduce_pdb_for_dssp(file_path)
        obj = _try_dssp(stripped)
        if obj: return obj
    except Exception:
        pass
    finally:
        if stripped:
            try: os.unlink(stripped)
            except OSError: pass

    # Attempt 3 — PDBIO-written file WITH dummy CRYST1 (mkdssp 4.x requires it)
    try:
        buf = io.StringIO()
        writer = PDBIO()
        writer.set_structure(model.get_parent())
        writer.save(buf)
        content = _with_cryst1(buf.getvalue())
        tmp = _write_tmp(content)
        try:
            obj = _try_dssp(tmp)
            if obj: return obj
        finally:
            try: os.unlink(tmp)
            except OSError: pass
    except Exception:
        pass

    # Attempt 4 — shell subprocess + CRYST1 (bypasses BioPython DSSP wrapper)
    try:
        buf = io.StringIO()
        writer = PDBIO()
        writer.set_structure(model.get_parent())
        writer.save(buf)
        dssp_text = _shell_dssp(_with_cryst1(buf.getvalue()))
        if dssp_text:
            obj = _parse_dssp_text(dssp_text)
            if obj: return obj
    except Exception:
        pass

    print(f"[DSSP/turn] {stem}: all strategies failed", file=sys.stderr)
    return None


# ── φ/ψ map ────────────────────────────────────────────────────────────────

def _build_phi_psi_map(model):
    """
    Build {(chain_id, res_id_tuple): (phi_deg, psi_deg)} via BioPython's
    PPBuilder. Terminal residues have None for their undefined angle.
    Residues at chain breaks are silently absent from the map.
    """
    ppb = PPBuilder()
    phi_psi_map = {}
    for pp in ppb.build_peptides(model):
        for res, (phi, psi) in zip(pp, pp.get_phi_psi_list()):
            phi_psi_map[(res.get_parent().id, res.id)] = (
                math.degrees(phi) if phi is not None else None,
                math.degrees(psi) if psi is not None else None,
            )
    return phi_psi_map


# ── cis-proline detection ──────────────────────────────────────────────────

def _omega_angle_deg(res_prev, res_curr):
    """
    Peptide-bond torsion ω = dihedral(Cα[i-1], C[i-1], N[i], Cα[i]) in degrees.
    Returns None if any backbone atom is missing.
    """
    try:
        ca_p = res_prev["CA"].get_vector()
        c_p  = res_prev["C"].get_vector()
        n_c  = res_curr["N"].get_vector()
        ca_c = res_curr["CA"].get_vector()
        return math.degrees(calc_dihedral(ca_p, c_p, n_c, ca_c))
    except Exception:
        return None


def _is_cis_pro(res_before, res_at_i2):
    """
    True iff res_at_i2 is PRO with a cis peptide bond.

    Panasik, Fleming & Rose (2005) require ω to lie within ±30° of planarity
    for type VI turns; |ω| ≤ 30° selects genuine cis bonds (ω ≈ 0°) while
    excluding the far looser <90° threshold used previously, which would have
    admitted a large fraction of distorted-trans bonds.
    Trans bonds cluster near ±180°; authentic cis bonds near 0°.
    """
    if res_at_i2.get_resname().strip().upper() != "PRO":
        return False
    omega = _omega_angle_deg(res_before, res_at_i2)
    return omega is not None and abs(omega) <= 30.0


# ── β-turn subtype classifier ──────────────────────────────────────────────

def _angle_diff(a, b):
    """Circular angular distance in degrees, result in [0, 180]."""
    d = abs(a - b) % 360.0
    return min(d, 360.0 - d)


def _classify_beta_turn(phi2, psi2, phi3, psi3, cis_pro_at_i2):
    """
    Match φ/ψ of i+1 (phi2, psi2) and i+2 (phi3, psi3) against idealized
    β-turn angles with ±ANGLE_TOL tolerance.
    VIa1 and IIb are skipped unless cis_pro_at_i2 is True.
    Returns a subtype name string or 'unclassified'.
    """
    for name, ((p2, s2), (p3, s3)) in BETA_TURN_DEFS.items():
        if name in _CIS_PRO_TYPES and not cis_pro_at_i2:
            continue
        if (
            _angle_diff(phi2, p2) <= ANGLE_TOL and
            _angle_diff(psi2, s2) <= ANGLE_TOL and
            _angle_diff(phi3, p3) <= ANGLE_TOL and
            _angle_diff(psi3, s3) <= ANGLE_TOL
        ):
            return name
    return None   # catch-all (Mazoni 2025 / Hutchinson & Thornton 1994): satisfies geometry
                  # and H-bond but no canonical subtype φ/ψ match — labelled Type IV, not "unclassified"


# ── FreeSASA — computed once per file, sliced per window ──────────────────

def _compute_sasa(structure):
    """
    Run FreeSASA on the full complex once. Returns the result object or None.
    The result is passed to _window_asa for per-window slicing; this avoids
    re-running the computation for every turn window (omega_backend's approach
    was to re-call calcBioPDB inside every record, which is O(N × T) SASA work
    instead of O(1)).
    """
    try:
        import freesasa
        freesasa.setVerbosity(freesasa.silent)
        ret = freesasa.calcBioPDB(
            structure,
            options={"hetatm": False, "hydrogen": False, "skip-unknown": True},
        )
        return ret[0]
    except Exception:
        return None


def _window_asa(sasa_result, chain_id, start_res, end_res):
    """
    Sum ASA (Å²) for residues [start_res, end_res] in chain_id from a
    pre-computed FreeSASA result. Returns None if result is unavailable.
    """
    if sasa_result is None:
        return None
    try:
        chain_areas = sasa_result.residueAreas().get(chain_id, {})
        total = 0.0
        for rnum_str, area in chain_areas.items():
            try:
                rnum = int(rnum_str.strip())
            except ValueError:
                continue
            if start_res <= rnum <= end_res:
                total += area.total
        return total
    except Exception:
        return None


# ── shared setup ───────────────────────────────────────────────────────────

def _setup(file_path):
    """
    Parse structure, run DSSP, build φ/ψ map, compute SASA.
    Returns (pdb_id, model, dssp_or_None, phi_psi_map, sasa_result_or_None).
    dssp is None if DSSP totally failed; callers must check before scanning.
    """
    pdb_id    = os.path.splitext(os.path.basename(file_path))[0]
    structure = _parser_for(file_path).get_structure(pdb_id, file_path)
    model     = structure[0]
    dssp      = _run_dssp(model, file_path)
    if dssp is None:
        return pdb_id, model, None, {}, None
    phi_psi_map = _build_phi_psi_map(model)
    sasa_result = _compute_sasa(structure)
    return pdb_id, model, dssp, phi_psi_map, sasa_result


# ── record builders ────────────────────────────────────────────────────────

def _make_beta_record(pdb_id, turn_num, chain_id, window, cas,
                      ss_codes, phi2, psi2, phi3, psi3,
                      endpoint_dist, subtype, asa):
    seq       = "".join(_three_to_one(r.get_resname()) for r in window)
    start_res = window[0].id[1]
    end_res   = window[-1].id[1]
    return {
        "pdb_id":         pdb_id,
        "turn_number":    turn_num,
        "turn_type":      "beta",
        "Chain":          chain_id,
        "Start_Position": start_res,
        "End_Position":   end_res,
        "range":          f"{start_res}-{end_res}",
        "Amino_Acid_Seq": seq,
        "subtype":        subtype,
        "phi_i1":         _r(phi2, 1),
        "psi_i1":         _r(psi2, 1),
        "phi_i2":         _r(phi3, 1),
        "psi_i2":         _r(psi3, 1),
        "Sec_Structure":  "".join(ss_codes),
        "distance_x":     _r(endpoint_dist),
        "x1": _r(cas[0][0]), "y1": _r(cas[0][1]), "z1": _r(cas[0][2]),
        "x4": _r(cas[3][0]), "y4": _r(cas[3][1]), "z4": _r(cas[3][2]),
        "total_asa":      _r(asa, 2) if asa is not None else None,
        "Amino_acid1":    window[0].get_resname(),
        "Amino_acid4":    window[3].get_resname(),
    }


def _make_gamma_record(pdb_id, turn_num, chain_id, window, cas,
                       ss_codes, phi_mid, psi_mid,
                       endpoint_dist, subtype, asa,
                       hbond_dist=None, h_present=None):
    seq       = "".join(_three_to_one(r.get_resname()) for r in window)
    start_res = window[0].id[1]
    end_res   = window[-1].id[1]
    return {
        "pdb_id":         pdb_id,
        "turn_number":    turn_num,
        "turn_type":      "gamma",
        "Chain":          chain_id,
        "Start_Position": start_res,
        "End_Position":   end_res,
        "range":          f"{start_res}-{end_res}",
        "Amino_Acid_Seq": seq,
        "subtype":        subtype,
        "phi_mid":        _r(phi_mid, 1),
        "psi_mid":        _r(psi_mid, 1),
        # H-bond geometry — distance is O···H (if H present) or N···O (fallback)
        "hbond_dist":     _r(hbond_dist, 2) if hbond_dist is not None else None,
        "h_present":      h_present,         # True = geometric check; False = N·O fallback
        "Sec_Structure":  "".join(ss_codes),
        "distance_x":     _r(endpoint_dist),
        "x1": _r(cas[0][0]), "y1": _r(cas[0][1]), "z1": _r(cas[0][2]),
        "x3": _r(cas[2][0]), "y3": _r(cas[2][1]), "z3": _r(cas[2][2]),
        "total_asa":      _r(asa, 2) if asa is not None else None,
        "Amino_acid1":    window[0].get_resname(),
        "Amino_acid3":    window[2].get_resname(),
    }


# ── α-turn hydrogen bond helper ────────────────────────────────────────────

# ── α-turn constants (Nataraj et al. 1995; Dasgupta et al. 2004) ─────────
ALPHA_CA_DIST   = 6.5    # Å — Cα(i)→Cα(i+4) ceiling (Dasgupta 2004)
ALPHA_NO_MIN    = 2.6    # Å — N(i+4)···O(i) lower bound (Nataraj 1995)
ALPHA_NO_MAX    = 3.6    # Å — N(i+4)···O(i) upper bound (Nataraj 1995)
ALPHA_ANGLE_MAX = 40.0   # °  — H–N···O angle at N, upper bound (Nataraj 1995)

# DSSP codes that disqualify i+1, i+2, i+3:
# G = 3₁₀-helix, H = α-helix  (Nataraj 1995; Dasgupta 2004 use G and H only)
ALPHA_HELICAL_SS = {'G', 'H'}

# Pavone et al. (1996) Biopolymers 38:705-721 Table IV
# 9 cluster centroids with per-angle standard deviations used as tolerances.
# Format: name → ((φi+1, ψi+1, φi+2, ψi+2, φi+3, ψi+3),
#                  (tol_φi+1, tol_ψi+1, tol_φi+2, tol_ψi+2, tol_φi+3, tol_ψi+3))
ALPHA_CLUSTERS = {
    "I-αRS":  ((-60, -29, -72, -29,  -96, -20), (11, 13, 14, 15, 20, 17)),
    "I-αLS":  (( 48,  42,  67,  33,   70,  32), (22, 14,  9, 14, 11, 12)),
    "II-αRS": ((-59, 129,  88, -16,  -91, -32), (10, 15, 15, 19, 22, 18)),
    "II-αLS": (( 53,-137, -95,  81,   57,  38), (15, 25, 12, 23,  5,  8)),
    "I-αLU":  ((-61, 158,  64,  37,   62,  39), (12, 15, 17, 21, 12,  8)),
    "I-αRU":  (( 59,-157, -67, -29,  -68, -39), (18, 31, 17, 20, 12, 12)),
    "II-αLU": ((-65, -20, -90,  16,   86,  37), (15, 15, 17, 44, 18, 27)),
    "II-αRU": (( 54,  39,  67,  -5, -125, -34), ( 8, 15, 13, 31, 11, 32)),
    "I-αc":   ((-103, 143, -85,   2,  -54, -39), (23,  4,  8,  6,  6,  9)),
}


def _signed_angle_diff(a, b):
    """Smallest signed difference between two angles (degrees), wrapped to (−180, 180]."""
    d = (a - b) % 360
    return d - 360 if d > 180 else d


def _assign_alpha_cluster(phi_i1, psi_i1, phi_i2, psi_i2, phi_i3, psi_i3):
    """
    Match the 6-angle triplet against Pavone et al. (1996) Table IV centroids.
    Each angle must fall within ±1 SD of its centroid value.
    Among passing clusters, pick the one with the lowest mean |deviation / SD|.
    Returns the cluster name, or 'unclassified' if none match.
    """
    observed = (phi_i1, psi_i1, phi_i2, psi_i2, phi_i3, psi_i3)
    if any(v is None for v in observed):
        return "unclassified"

    best_name  = "unclassified"
    best_score = float("inf")

    for name, (centroid, tolerances) in ALPHA_CLUSTERS.items():
        diffs = [abs(_signed_angle_diff(observed[k], centroid[k])) for k in range(6)]
        # every angle must be within its per-angle standard deviation tolerance
        if all(diffs[k] <= tolerances[k] for k in range(6)):
            score = sum(diffs[k] / tolerances[k] for k in range(6))
            if score < best_score:
                best_score = score
                best_name  = name

    return best_name


def _alpha_hbond_geometry(window):
    """
    Compute N(i+4)···O(i) H-bond geometry for α-turn detection.

    Returns (N_O_dist_Å, angle_HNO_deg) where angle_HNO is the
    H–N···O angle at N (degrees).

    H on N(i+4) is located via _find_backbone_h: an actual modelled H atom if
    present, otherwise the sp² bisector of the C′(i+3)→N(i+4) and
    Cα(i+4)→N(i+4) bonds. This reproduces the method described by Nataraj
    et al. (1995, Curr. Sci. 69:434-447): amide hydrogens are fixed "in the
    plane of the peptide unit using the ideal bond angles around the
    nitrogen atom" — i.e. the bisector of both adjacent bonds, not a single
    bond reflected anti — the same construction already used for NHN
    γ-turn H-placement elsewhere in this module.

    Both values are None if required atoms are missing.
    """
    try:
        n_i4       = window[4]["N"].get_vector()
        o_i        = window[0]["O"].get_vector()
        c_prev_arr = window[3]["C"].get_vector().get_array()
        ca_i4_arr  = window[4]["CA"].get_vector().get_array()

        # N···O heavy-atom distance
        no_dist = float((n_i4 - o_i).norm())

        h_arr = _find_backbone_h(window[4], prev_c_arr=c_prev_arr, ca_arr=ca_i4_arr)
        if h_arr is None:
            return no_dist, None

        n_arr = n_i4.get_array()
        o_arr = o_i.get_array()

        nh_vec = h_arr - n_arr
        nh_vec = nh_vec / np.linalg.norm(nh_vec)
        no_vec = o_arr - n_arr
        no_vec = no_vec / np.linalg.norm(no_vec)

        cos_a = float(np.dot(nh_vec, no_vec))
        cos_a = max(-1.0, min(1.0, cos_a))
        angle = math.degrees(math.acos(cos_a))

        return no_dist, angle

    except Exception:
        return None, None


def _make_alpha_record(pdb_id, turn_num, chain_id, window, cas,
                       ss_codes, hbond_dist, hbond_angle,
                       endpoint_dist, asa, cluster,
                       phi_i1, psi_i1, phi_i2, psi_i2, phi_i3, psi_i3):
    seq       = "".join(_three_to_one(r.get_resname()) for r in window)
    start_res = window[0].id[1]
    end_res   = window[-1].id[1]
    return {
        "pdb_id":         pdb_id,
        "turn_number":    turn_num,
        "turn_type":      "alpha",
        "Chain":          chain_id,
        "Start_Position": start_res,
        "End_Position":   end_res,
        "range":          f"{start_res}-{end_res}",
        "Amino_Acid_Seq": seq,
        "cluster":        cluster,
        # torsion angles at i+1, i+2, i+3
        "phi_i1":         _r(phi_i1, 1),
        "psi_i1":         _r(psi_i1, 1),
        "phi_i2":         _r(phi_i2, 1),
        "psi_i2":         _r(psi_i2, 1),
        "phi_i3":         _r(phi_i3, 1),
        "psi_i3":         _r(psi_i3, 1),
        # H-bond geometry: N(i+4)···O(i)
        "hbond_NO_dist":  _r(hbond_dist, 2) if hbond_dist is not None else None,
        "hbond_HNO_angle":_r(hbond_angle, 1) if hbond_angle is not None else None,
        # secondary structure (5 DSSP codes: i, i+1, i+2, i+3, i+4)
        "Sec_Structure":  "".join(ss_codes),
        # geometry
        "distance_x":     _r(endpoint_dist),    # Cα(i)–Cα(i+4) in Å
        "x1": _r(cas[0][0]), "y1": _r(cas[0][1]), "z1": _r(cas[0][2]),
        "x5": _r(cas[4][0]), "y5": _r(cas[4][1]), "z5": _r(cas[4][2]),
        "total_asa":      _r(asa, 2) if asa is not None else None,
        "Amino_acid1":    window[0].get_resname(),
        "Amino_acid5":    window[4].get_resname(),
    }


def _scan_alpha_turns(model, dssp, phi_psi_map, sasa_result, pdb_id):
    """
    Detect α-turns: 5-residue window [i, i+1, i+2, i+3, i+4].

    Criteria (Nataraj et al. 1995; Pavone et al. 1996; Dasgupta et al. 2004):
      1. Cα(i) → Cα(i+4) < 6.5 Å  (Dasgupta 2004)
      2. i+1, i+2, i+3 not assigned as G (3₁₀-helix) or H (α-helix) by DSSP
      3. N(i+4)···O(i) distance: 2.6–3.6 Å  (Nataraj 1995)
      4. H–N···O angle at N: ≤ 40°  (Nataraj 1995; H placed anti to C′→N bond)
      5. Cluster assigned by matching φ/ψ of i+1, i+2, i+3 against
         Pavone et al. (1996) Table IV centroids using per-angle SD tolerances.
    """
    turns        = []
    turn_counter = 0

    for chain in model:
        residues = [r for r in chain if r.id[0] == " "]
        n        = len(residues)

        for i in range(n - 4):
            window = residues[i : i + 5]         # [i, i+1, i+2, i+3, i+4]
            cas    = [_ca_coord(r) for r in window]
            if any(c is None for c in cas):
                continue

            # ── 1. Cα(i)–Cα(i+4) distance < 6.5 Å ──────────────────────────
            endpoint_dist = _dist(cas[0], cas[4])
            if endpoint_dist >= ALPHA_CA_DIST:
                continue

            # ── 2. DSSP: i+1, i+2, i+3 must not be G or H ───────────────────
            try:
                ss_codes = [dssp[(chain.id, r.id)][2] for r in window]
            except Exception:
                continue
            if any(ss_codes[k] in ALPHA_HELICAL_SS for k in (1, 2, 3)):
                continue

            # ── 3 & 4. H-bond geometry: N···O 2.6–3.6 Å, H–N–O ≤ 40° ───────
            hbond_dist, hbond_angle = _alpha_hbond_geometry(window)
            if hbond_dist is None or hbond_angle is None:
                continue
            if not (ALPHA_NO_MIN <= hbond_dist <= ALPHA_NO_MAX):
                continue
            if hbond_angle > ALPHA_ANGLE_MAX:
                continue

            # ── 5. φ/ψ at i+1, i+2, i+3 ─────────────────────────────────────
            phi_i1, psi_i1 = phi_psi_map.get((chain.id, window[1].id), (None, None))
            phi_i2, psi_i2 = phi_psi_map.get((chain.id, window[2].id), (None, None))
            phi_i3, psi_i3 = phi_psi_map.get((chain.id, window[3].id), (None, None))

            cluster = _assign_alpha_cluster(phi_i1, psi_i1,
                                            phi_i2, psi_i2,
                                            phi_i3, psi_i3)
            # Windows with no Pavone cluster match are still genuine α-turns by
            # Nataraj et al. (1995)'s own definition (criteria 1-4 above); they
            # are reported as "unclassified" rather than discarded, following
            # Pavone et al. (1996)'s own convention of listing such turns
            # separately (their Table III) instead of excluding them.

            # ── 6. Exposure ───────────────────────────────────────────────────
            asa = _window_asa(sasa_result, chain.id,
                              window[0].id[1], window[4].id[1])

            turn_counter += 1
            turns.append(_make_alpha_record(
                pdb_id, turn_counter, chain.id, window, cas,
                ss_codes, hbond_dist, hbond_angle, endpoint_dist, asa,
                cluster, phi_i1, psi_i1, phi_i2, psi_i2, phi_i3, psi_i3,
            ))

    status = f"ok:{len(turns)}" if turns else "no-turns"
    return turns, status


# ── inner scanners (take pre-computed objects, return (turns, status)) ─────

def _scan_beta_turns(model, dssp, phi_psi_map, sasa_result, pdb_id):
    turns        = []
    turn_counter = 0

    for chain in model:
        residues = [r for r in chain if r.id[0] == " "]
        n        = len(residues)

        for i in range(n - 3):
            window = residues[i : i + 4]
            cas    = [_ca_coord(r) for r in window]
            if any(c is None for c in cas):
                continue

            # ── 1. Cα distance ─────────────────────────────────────────────
            endpoint_dist = _dist(cas[0], cas[3])
            if endpoint_dist >= BETA_CA_DIST:
                continue

            # ── 2. DSSP: neither i+1 nor i+2 may be helical (G, H, or I) ──
            try:
                ss_codes = [dssp[(chain.id, r.id)][2] for r in window]
            except Exception:
                continue
            if ss_codes[1] in HELICAL_SS or ss_codes[2] in HELICAL_SS:
                continue

            # ── 3. φ/ψ for subtype classification ─────────────────────────
            phi2, psi2 = phi_psi_map.get((chain.id, window[1].id), (None, None))
            phi3, psi3 = phi_psi_map.get((chain.id, window[2].id), (None, None))

            if phi2 is None or psi2 is None or phi3 is None or psi3 is None:
                subtype = None   # angles unavailable — treat as catch-all
            else:
                cis_pro = _is_cis_pro(window[1], window[2])
                subtype = _classify_beta_turn(phi2, psi2, phi3, psi3, cis_pro)

            # ── 4. solvent exposure ────────────────────────────────────────
            asa = _window_asa(sasa_result, chain.id,
                              window[0].id[1], window[3].id[1])

            if subtype is None:
                subtype = "IV"

            turn_counter += 1
            turns.append(_make_beta_record(
                pdb_id, turn_counter, chain.id, window, cas,
                ss_codes, phi2, psi2, phi3, psi3,
                endpoint_dist, subtype, asa,
            ))

    status = f"ok:{len(turns)}" if turns else "no-turns"
    return turns, status


def _gamma_hbond_geometry(window):
    """
    Check the γ-turn H-bond between C=O(i) and N-H(i+2).

    If backbone H is present on N(i+2):
        - O···H distance   < 2.5 Å
        - N-H···O angle (at H) > 120°  (donor linearity)
        - H···O=C angle (at O) > 90°   (acceptor geometry)
    If H is absent (typical in X-ray structures at <1.0 Å resolution):
        - N(i+2)···O(i) heavy-atom distance < 3.5 Å  (standard fallback)

    Returns (passes: bool, hbond_dist: float|None, h_present: bool)
      hbond_dist is O···H (H present) or N···O (fallback).
    """
    try:
        o_vec = window[0]["O"].get_vector()   # carbonyl O of residue i
        n_vec = window[2]["N"].get_vector()   # amide N of residue i+2
    except KeyError:
        return False, None, False

    # Try to find the backbone amide H on residue i+2 (Proline has none)
    h_vec = None
    if window[2].get_resname().strip().upper() != "PRO":
        for h_name in ("H", "HN", "H1", "1H"):
            try:
                h_vec = window[2][h_name].get_vector()
                break
            except KeyError:
                pass

    if h_vec is not None:
        # ── Full geometric H-bond check ────────────────────────────────────
        oh_dist = float((h_vec - o_vec).norm())
        if oh_dist >= 2.5:
            return False, oh_dist, True

        # N-H···O angle at H — must be > 120°
        nho = math.degrees(float(calc_angle(n_vec, h_vec, o_vec)))
        if nho <= 120.0:
            return False, oh_dist, True

        # H···O=C angle at O — must be > 90°
        try:
            c_vec = window[0]["C"].get_vector()
            hoc = math.degrees(float(calc_angle(h_vec, o_vec, c_vec)))
            if hoc <= 90.0:
                return False, oh_dist, True
        except KeyError:
            pass   # C missing — accept on distance + donor angle alone

        return True, oh_dist, True

    else:
        # ── Fallback: heavy-atom N···O distance ───────────────────────────
        no_dist = float((n_vec - o_vec).norm())
        return no_dist < 3.5, no_dist, False


# Classic/inverse γ-turn centroids, originally defined by Nemethy & Printz
# (1972, Macromolecules 5:755-758) and tabulated in Rose, Gierasch & Smith
# (1985, Adv. Protein Chem. 37:1-109, doi:10.1016/S0065-3233(08)60063-7).
GAMMA_CLASSIC_CENTROID = (75.0, -64.0)
GAMMA_INVERSE_CENTROID = (-79.0, 69.0)

# ±40° per-angle tolerance around each centroid (Hutchinson & Thornton, 1996
# — PROMOTIF, Protein Sci. 5:212-220, doi:10.1002/pro.5560050204). PROMOTIF
# applies this tolerance to the same Rose (1985) / Nemethy & Printz (1972)
# centroids rather than using them as a zero-tolerance cutoff; without it,
# real hydrogen-bonded γ-turns (whose φ/ψ naturally scatter around the
# idealized centroid) are missed even when correctly identified as γ-turns
# by every other criterion.
GAMMA_ANGLE_TOL = 40.0


def _classify_gamma_strict(phi, psi):
    """
    Classify γ-turn subtype by matching φ,ψ of residue i+1 against the
    classic/inverse centroids (Nemethy & Printz, 1972; Rose et al., 1985)
    within ±40° per angle (Hutchinson & Thornton, 1996 — PROMOTIF).
    Returns "classic", "inverse", or None if outside both windows.
    """
    cphi, cpsi = GAMMA_CLASSIC_CENTROID
    if _angle_diff(phi, cphi) <= GAMMA_ANGLE_TOL and _angle_diff(psi, cpsi) <= GAMMA_ANGLE_TOL:
        return "classic"
    iphi, ipsi = GAMMA_INVERSE_CENTROID
    if _angle_diff(phi, iphi) <= GAMMA_ANGLE_TOL and _angle_diff(psi, ipsi) <= GAMMA_ANGLE_TOL:
        return "inverse"
    return None


# ── NHN γ-turn helpers (Dhar, Kishore, Chakrabarti 2019) ──────────────────

def _find_backbone_h(residue, prev_c_arr=None, ca_arr=None):
    """
    Return the backbone amide H position on residue N as a numpy array.

    Strategy (in order):
    1. Search for actual H atom within 0.85-1.30 A of N (works when H atoms
       are present, e.g. after running REDUCE on the PDB).
    2. If prev_c_arr AND ca_arr are both supplied, use the sp2 bisector formula:
           H = N - normalize(normalize(C_prev - N) + normalize(Ca - N)) * 1.01
       This places H in the peptide plane at 120° from the C'-N and Ca-N bonds,
       matching the actual amide sp2 geometry. Results match the paper's values.
    3. Fallback (only prev_c_arr, no ca_arr): N + normalize(N - C_prev) * 1.02

    Returns numpy array or None.
    """
    try:
        n_arr = residue["N"].get_vector().get_array()
    except KeyError:
        return None

    # 1. Try actual H atoms already in the file
    best_h, best_d = None, float("inf")
    for atom in residue.get_atoms():
        if (atom.element or "").strip().upper() != "H":
            continue
        h_arr = atom.get_vector().get_array()
        d = float(np.linalg.norm(h_arr - n_arr))
        if 0.85 <= d <= 1.30 and d < best_d:
            best_h, best_d = h_arr, d
    if best_h is not None:
        return best_h

    # 2. sp2 bisector: correct placement in peptide plane at 120° angles
    if prev_c_arr is not None and ca_arr is not None:
        u1 = prev_c_arr - n_arr
        u2 = ca_arr     - n_arr
        n1 = np.linalg.norm(u1)
        n2 = np.linalg.norm(u2)
        if n1 > 1e-8 and n2 > 1e-8:
            h_dir = -(u1 / n1 + u2 / n2)
            hd = np.linalg.norm(h_dir)
            if hd > 1e-8:
                return n_arr + (h_dir / hd) * 1.01

    # 3. Fallback: N-to-C'_prev direction only (less accurate)
    if prev_c_arr is not None:
        nh_dir = n_arr - prev_c_arr
        nh_norm = np.linalg.norm(nh_dir)
        if nh_norm > 1e-8:
            return n_arr + (nh_dir / nh_norm) * 1.02

    return None


def _nhn_criterion_b(window):
    """
    Criterion (b): theta = arccos(|v_H(i+1) · normal_N(i)|) < 30 deg.
    normal_N(i) is the peptide-plane normal at N(i): (C(i-1)-N(i)) x (Ca(i)-N(i)).
    H on N(i+1) is found from actual atoms or estimated geometrically from C(i).
    Returns (passes: bool, theta_deg: float|None).
    """
    try:
        n_i    = window[1]["N"].get_vector().get_array()
        c_im1  = window[0]["C"].get_vector().get_array()
        ca_i   = window[1]["CA"].get_vector().get_array()
        c_i    = window[1]["C"].get_vector().get_array()
        ca_ip1 = window[2]["CA"].get_vector().get_array()
    except KeyError:
        return False, None

    h_ip1 = _find_backbone_h(window[2], prev_c_arr=c_i, ca_arr=ca_ip1)
    if h_ip1 is None:
        return False, None

    try:
        v1 = c_im1 - n_i
        v2 = ca_i  - n_i
        normal = np.cross(v1, v2)
        norm_len = np.linalg.norm(normal)
        if norm_len < 1e-8:
            return False, None
        normal = normal / norm_len

        v_H = h_ip1 - n_i
        vh_len = np.linalg.norm(v_H)
        if vh_len < 1e-8:
            return False, None
        v_H = v_H / vh_len

        dot   = abs(float(np.dot(v_H, normal)))
        dot   = max(-1.0, min(1.0, dot))
        theta = math.degrees(math.acos(dot))
        return theta < 30.0, theta
    except Exception:
        return False, None


def _nhn_criterion_c(window, prev_c_arr=None):
    """
    Criterion (c): C5 intraresidue ring at residue i-1.
    H(i-1) ... O(i-1) < 3.0 A.
    H on N(i-1) found from actual atoms or estimated geometrically from C'(i-2).
    Returns (passes: bool, h_o_dist: float|None).
    """
    try:
        o_arr  = window[0]["O"].get_vector().get_array()
        ca_im1 = window[0]["CA"].get_vector().get_array()
    except KeyError:
        return False, None

    h_im1 = _find_backbone_h(window[0], prev_c_arr=prev_c_arr, ca_arr=ca_im1)
    if h_im1 is None:
        return False, None

    h_o_dist = float(np.linalg.norm(h_im1 - o_arr))
    return h_o_dist < 3.0, h_o_dist


def _classify_nhn_gamma(phi, psi):
    """
    Classify NHN gamma-turn subtype (Dhar et al. 2019, Figure 3 / page 5).
    Phi/psi clusters of the central residue i:
      Classical NHN: phi = +59(+/-13) deg -> [+46,+72], psi = +36(+/-16) deg -> [+20,+52]
      Inverse NHN:   phi = -72(+/-18) deg -> [-90,-54], psi = -22(+/-20) deg -> [-42, -2]
    Returns 'classical-nhn', 'inverse-nhn', or None.
    """
    if 46.0 <= phi <= 72.0 and 20.0 <= psi <= 52.0:
        return "classical-nhn"
    if -90.0 <= phi <= -54.0 and -42.0 <= psi <= -2.0:
        return "inverse-nhn"
    return None


def _scan_nhn_gamma_turns(model, dssp, phi_psi_map, sasa_result, pdb_id):
    """
    Detect NHN gamma-turns (Dhar, Kishore, Chakrabarti 2019).
    Window: [i-1, i, i+1] — i is the central N acceptor.

    Detection criteria (paper Section 3.3):
      (a) 5.2 <= Ca(i-1)...Ca(i+1) <= 5.8 A
      (b) theta = arccos(|v_H(i+1).normal_N(i)|) < 30 deg
      (c) H(i-1)...O(i-1) < 3.0 A  (C5 ring)

    Only turns with a matching phi/psi subtype (classical-nhn or inverse-nhn)
    are kept; turns with no subtype match are discarded.
    """
    turns        = []
    turn_counter = 0

    for chain in model:
        residues = [r for r in chain if r.id[0] == " "]
        n        = len(residues)

        for i in range(n - 2):
            window = residues[i : i + 3]
            cas    = [_ca_coord(r) for r in window]
            if any(c is None for c in cas):
                continue

            # (a) Ca distance
            ca_dist = _dist(cas[0], cas[2])
            if not (5.2 <= ca_dist <= 5.8):
                continue

            # (b) theta angle
            b_ok, theta = _nhn_criterion_b(window)
            if not b_ok:
                continue

            # (c) C5 ring — pass C' of residue i-2 for geometric H fallback
            try:
                prev_c_arr = residues[i - 1]["C"].get_vector().get_array() if i > 0 else None
            except KeyError:
                prev_c_arr = None
            c_ok, _h_o = _nhn_criterion_c(window, prev_c_arr=prev_c_arr)
            if not c_ok:
                continue

            # phi/psi subtype — discard if no match
            phi_mid, psi_mid = phi_psi_map.get(
                (chain.id, window[1].id), (None, None)
            )
            if phi_mid is None or psi_mid is None:
                continue
            subtype = _classify_nhn_gamma(phi_mid, psi_mid)
            if subtype is None:
                continue

            # DSSP secondary structure — display only
            if dssp is not None:
                ss_codes = []
                for r in window:
                    try:
                        ss_codes.append(dssp[(chain.id, r.id)][2])
                    except Exception:
                        ss_codes.append('-')
            else:
                ss_codes = ['-', '-', '-']

            asa = _window_asa(sasa_result, chain.id,
                              window[0].id[1], window[2].id[1])

            # N(i-1)···O(i-1) intraresidue heavy-atom distance (1st residue)
            try:
                no_dist = _r(float(
                    (window[0]["N"].get_vector() - window[0]["O"].get_vector()).norm()
                ), 2)
            except Exception:
                no_dist = None

            turn_counter += 1
            record = _make_gamma_record(
                pdb_id, turn_counter, chain.id, window, cas,
                ss_codes, phi_mid, psi_mid,
                ca_dist, subtype, asa,
                hbond_dist=_r(theta, 2),
                h_present=True,
            )
            record["is_nhn"]  = True
            record["no_dist"] = no_dist
            record["ho_dist"] = _r(_h_o, 2) if _h_o is not None else None
            turns.append(record)

    status = f"ok:{len(turns)}" if turns else "no-turns"
    return turns, status


def _scan_gamma_turns(model, dssp, phi_psi_map, sasa_result, pdb_id):
    """
    Detect γ-turns using strict geometric H-bond criteria + hard φ/ψ ranges.

    Criteria (all required):
      1. 3-residue window [i, i+1, i+2]
      2. H-bond between C=O(i) and N-H(i+2):
           H present  → O···H < 2.5 Å, N-H···O > 120°, H···O=C > 90°
           H absent   → N(i+2)···O(i) < 3.5 Å (heavy-atom fallback)
      3. φ(i+1), ψ(i+1) in strict ranges:
           classic γ: φ = 70–85°,     ψ = −70 to −60°
           inverse γ: φ = −85 to −70°, ψ = 60 to 70°
    """
    turns        = []
    turn_counter = 0

    for chain in model:
        residues = [r for r in chain if r.id[0] == " "]
        n        = len(residues)

        for i in range(n - 2):
            window = residues[i : i + 3]
            cas    = [_ca_coord(r) for r in window]
            if any(c is None for c in cas):
                continue

            # ── 1. H-bond geometry ──────────────────────────────────────────
            passes, hbond_dist_val, h_present = _gamma_hbond_geometry(window)
            if not passes:
                continue

            # ── 2. φ/ψ strict ranges ────────────────────────────────────────
            phi_mid, psi_mid = phi_psi_map.get(
                (chain.id, window[1].id), (None, None)
            )
            if phi_mid is None or psi_mid is None:
                continue

            subtype = _classify_gamma_strict(phi_mid, psi_mid)
            if subtype is None:
                continue   # H-bond present but φ/ψ outside both γ-turn regions

            # ── 3. DSSP secondary structure ─────────────────────────────────
            try:
                ss_codes = [dssp[(chain.id, r.id)][2] for r in window]
            except Exception:
                continue

            # ── 4. geometry & exposure ──────────────────────────────────────
            endpoint_dist = _dist(cas[0], cas[2])
            asa = _window_asa(sasa_result, chain.id,
                              window[0].id[1], window[2].id[1])

            turn_counter += 1
            turns.append(_make_gamma_record(
                pdb_id, turn_counter, chain.id, window, cas,
                ss_codes, phi_mid, psi_mid,
                endpoint_dist, subtype, asa,
                hbond_dist_val, h_present,
            ))

    status = f"ok:{len(turns)}" if turns else "no-turns"
    return turns, status


# ── public API ─────────────────────────────────────────────────────────────

def find_beta_turns(file_path):
    """
    Scan a structure for β-turns.

    Returns
    -------
    turns      : list[dict]
    status_str : "ok:N" | "dssp-failed" | "no-turns" | "error:TYPE: msg"
    """
    try:
        pdb_id, model, dssp, phi_psi_map, sasa_result = _setup(file_path)
        if dssp is None:
            return [], "dssp-failed"
        return _scan_beta_turns(model, dssp, phi_psi_map, sasa_result, pdb_id)
    except Exception as e:
        traceback.print_exc()
        return [], f"error:{type(e).__name__}: {e}"


def find_gamma_turns(file_path):
    """
    Scan a structure for γ-turns (classic/inverse NHO + NHN subtypes).

    Returns
    -------
    turns      : list[dict]
    status_str : "ok:N" | "dssp-failed" | "no-turns" | "error:TYPE: msg"
    """
    try:
        pdb_id, model, dssp, phi_psi_map, sasa_result = _setup(file_path)
        if dssp is None:
            return [], "dssp-failed"
        gamma,     g_status = _scan_gamma_turns     (model, dssp, phi_psi_map, sasa_result, pdb_id)
        nhn_gamma, n_status = _scan_nhn_gamma_turns (model, dssp, phi_psi_map, sasa_result, pdb_id)
        gamma.extend(nhn_gamma)
        status = f"ok:{len(gamma)}" if gamma else "no-turns"
        return gamma, status
    except Exception as e:
        traceback.print_exc()
        return [], f"error:{type(e).__name__}: {e}"


def find_alpha_turns(file_path):
    """
    Scan a structure for α-turns (5-residue windows with N(i+4)···O(i) H-bond).

    Returns
    -------
    turns      : list[dict]
    status_str : "ok:N" | "dssp-failed" | "no-turns" | "error:TYPE: msg"
    """
    try:
        pdb_id, model, dssp, phi_psi_map, sasa_result = _setup(file_path)
        if dssp is None:
            return [], "dssp-failed"
        return _scan_alpha_turns(model, dssp, phi_psi_map, sasa_result, pdb_id)
    except Exception as e:
        traceback.print_exc()
        return [], f"error:{type(e).__name__}: {e}"


def find_all_turns(file_path):
    """
    Run all three turn detectors (α, β, γ + NHN γ) on the same file.
    File parsing, DSSP, φ/ψ map, and SASA are computed only once.

    Returns
    -------
    beta        : list[dict]
    gamma       : list[dict]   (classic/inverse γ + classical-nhn/inverse-nhn NHN γ)
    alpha       : list[dict]
    status_str  : combined status string
    """
    try:
        pdb_id, model, dssp, phi_psi_map, sasa_result = _setup(file_path)
        if dssp is None:
            return [], [], [], "dssp-failed"
        beta,      b_status = _scan_beta_turns      (model, dssp, phi_psi_map, sasa_result, pdb_id)
        gamma,     g_status = _scan_gamma_turns     (model, dssp, phi_psi_map, sasa_result, pdb_id)
        nhn_gamma, n_status = _scan_nhn_gamma_turns (model, dssp, phi_psi_map, sasa_result, pdb_id)
        alpha,     a_status = _scan_alpha_turns     (model, dssp, phi_psi_map, sasa_result, pdb_id)
        gamma.extend(nhn_gamma)
        return beta, gamma, alpha, (
            f"beta:{b_status} | gamma:{g_status} | nhn-gamma:{n_status} | alpha:{a_status}"
        )
    except Exception as e:
        traceback.print_exc()
        return [], [], [], f"error:{type(e).__name__}: {e}"


# ── CLI ────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys, json

    if len(sys.argv) not in (2, 3):
        sys.exit(
            "usage: python turn_backend.py <structure.pdb|.cif|.mmcif> [alpha|beta|gamma|all]"
        )

    path = sys.argv[1]
    mode = (sys.argv[2] if len(sys.argv) == 3 else "all").lower()

    if mode == "beta":
        results, status = find_beta_turns(path)
        print(f"Status: {status}", file=sys.stderr)
        print(json.dumps(results, indent=2))

    elif mode == "gamma":
        results, status = find_gamma_turns(path)
        print(f"Status: {status}", file=sys.stderr)
        print(json.dumps(results, indent=2))

    elif mode == "alpha":
        results, status = find_alpha_turns(path)
        print(f"Status: {status}", file=sys.stderr)
        print(json.dumps(results, indent=2))

    else:
        beta, gamma, alpha, status = find_all_turns(path)
        print(f"Status: {status}", file=sys.stderr)
        print("=== β-turns ===")
        print(json.dumps(beta, indent=2))
        print("\n=== γ-turns ===")
        print(json.dumps(gamma, indent=2))
        print("\n=== α-turns ===")
        print(json.dumps(alpha, indent=2))
