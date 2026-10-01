"""
omega_backend.py
----------------
Finds Ω-loops in protein structure files (.pdb / .cif / .mmcif).

An Ω-loop is a 6-residue window where:
  - At least 5 of the 6 DSSP codes are S or T (bend / turn)
  - All codes are in {S, T, E, B, G, H, I, C}
  - Cα distance between first and last residue < 10 Å  (loop "closes")
  - That endpoint distance < (2/3) × the max pairwise Cα distance in the window

DSSP is run automatically on the uploaded file — no pre-computed .dssp needed.
ASA is computed via FreeSASA (same library motif_backend already uses) so
there is no dependency on external .asa files either.

Returns plain Python lists of dicts, exactly like motif_backend.find_motifs.

Public API
----------
find_omega_loops(file_path) -> (loops, status_str)
    loops      : list of dicts (see _make_record for field names)
    status_str : "ok:N", "dssp-failed", "no-loops", or "error:..."
"""

import os
import math
import traceback
from collections import OrderedDict

from Bio.PDB import PDBParser, MMCIFParser
# BioPython ≥1.80 moved three_to_one to Bio.Data.IUPACData / Bio.SeqUtils;
# older versions had it in Bio.PDB.Polypeptide. Try both.
try:
    from Bio.PDB.Polypeptide import three_to_one as bio_three_to_one
except ImportError:
    try:
        from Bio.SeqUtils import IUPACData
        _d = {k.upper(): v for k, v in IUPACData.protein_letters_3to1_extended.items()}
        def bio_three_to_one(resname):
            r = _d.get(resname.upper())
            if r is None:
                raise KeyError(resname)
            return r
    except Exception:
        # Last-resort: inline the standard 20
        _AA3TO1 = {
            "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C",
            "GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
            "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P",
            "SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
        }
        def bio_three_to_one(resname):
            r = _AA3TO1.get(resname.upper())
            if r is None:
                raise KeyError(resname)
            return r

# ── constants ──────────────────────────────────────────────────────────────
MIN_LOOP_LENGTH = 6
MAX_LOOP_LENGTH = 16

MAX_END_DIST     = 10.0       # Å — max Cα endpoint distance
ALLOWED_SS       = set("STEBGHIC")   # valid DSSP codes inside the loop
ST_CODES         = {"S", "T"}


# ── helpers ────────────────────────────────────────────────────────────────
def _ca_coord(residue):
    """Return (x, y, z) for the Cα of a residue, or None."""
    try:
        coord = residue["CA"].get_coord()
        return (float(coord[0]), float(coord[1]), float(coord[2]))
    except Exception:
        return None


def _dist(a, b):
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def _sphere_geometry(coords):
    """
    Bounding sphere matching the original PhD script's formula exactly.
    center   = midpoint of first and last Cα
    bounding_r = max distance from that center to any Cα in the window
                 (the PhD script stores this in a variable called 'diameter' —
                  confusing name, but we preserve the same arithmetic)
    radius   = bounding_r / 2   (matches original script)
    diameter = bounding_r       (matches original script's 'diameter' column)
    """
    if not coords:
        return (None,) * 7
    center     = [(coords[0][i] + coords[-1][i]) / 2.0 for i in range(3)]
    bounding_r = max(_dist(center, c) for c in coords)   # original's "diameter"
    radius     = bounding_r / 2.0                         # original's "radius"
    diameter   = bounding_r                               # keep original column value
    area       = 4 * math.pi * radius ** 2
    volume     = (4 / 3) * math.pi * radius ** 3
    return radius, area, volume, center[0], center[1], center[2], diameter


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


# ── DSSP runner ────────────────────────────────────────────────────────────

def _find_dssp_binary():
    """Find DSSP binary — checks local folder first so Windows users can
    place dssp.exe / mkdssp.exe alongside app.py without touching PATH."""
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
    return "dssp"


_CRYST1_DUMMY = ("CRYST1    1.000    1.000    1.000"
                 "  90.00  90.00  90.00 P 1           1\n")

_DSSP_KEEP = ('CRYST1', 'ATOM  ', 'HETATM', 'TER   ', 'TER\n', 'END   ', 'END\n')


def _strip_for_dssp(file_path):
    """Return cleaned PDB text: only ATOM/HETATM/TER/END lines + dummy CRYST1 if missing.
    Avoids PDBIO which writes bad element codes that mkdssp 4.x rejects."""
    with open(file_path, errors="replace") as fh:
        lines = [l for l in fh if any(l.startswith(k) for k in _DSSP_KEEP)]
    if not any(l.startswith('CRYST1') for l in lines):
        lines = [_CRYST1_DUMMY] + lines
    # Drop alternate conformations — keep blank altloc or 'A'
    cleaned = []
    for l in lines:
        if l.startswith(('ATOM  ', 'HETATM')):
            altloc = l[16] if len(l) > 16 else ' '
            if altloc not in (' ', '', 'A'):
                continue
        cleaned.append(l)
    return "".join(cleaned)


def _run_dssp(model, file_path):
    """
    Run DSSP and return a DSSP-compatible object, or None on failure.
    Avoids PDBIO (produces bad element codes for mkdssp 4.x).
    Strategies: direct → text-strip+CRYST1 → shell subprocess on stripped text.
    """
    import sys, subprocess, tempfile, warnings as _warnings
    from Bio.PDB.DSSP import DSSP

    binary = _find_dssp_binary()
    stem   = os.path.basename(file_path)

    def _try_dssp(path):
        try:
            with _warnings.catch_warnings():
                _warnings.simplefilter("ignore")
                obj = DSSP(model, path, dssp=binary)
            if obj and len(obj.keys()) > 0:
                return obj
        except Exception:
            pass
        return None

    def _run_binary_dssp(pdb_content):
        """Run mkdssp via subprocess on pdb_content text; return DSSP text or None."""
        tmp_in = tmp_out = None
        try:
            with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
                f.write(pdb_content); tmp_in = f.name
            with tempfile.NamedTemporaryFile(suffix=".dssp", mode="w", delete=False) as f:
                tmp_out = f.name
            cmd = f'"{binary}" --output-format=dssp "{tmp_in}" -o "{tmp_out}"'
            subprocess.run(cmd, shell=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=60)
            if os.path.exists(tmp_out):
                with open(tmp_out) as fh:
                    text = fh.read()
                if "RESIDUE" in text:
                    return text
        except Exception as e:
            print(f"[DSSP/omega] {stem}: shell call failed — {e}", file=sys.stderr)
        finally:
            for f in (tmp_in, tmp_out):
                if f and os.path.exists(f):
                    try: os.unlink(f)
                    except OSError: pass
        return None

    def _dssp_from_text(dssp_text):
        tmp = None
        try:
            with tempfile.NamedTemporaryFile(suffix=".dssp", mode="w", delete=False) as f:
                f.write(dssp_text); tmp = f.name
            with _warnings.catch_warnings():
                _warnings.simplefilter("ignore")
                obj = DSSP(model, tmp, file_type="DSSP")
            return obj if obj and len(obj.keys()) > 0 else None
        except Exception as e:
            print(f"[DSSP/omega] {stem}: parse failed — {e}", file=sys.stderr)
            return None
        finally:
            if tmp and os.path.exists(tmp):
                try: os.unlink(tmp)
                except OSError: pass

    # Strategy 1: direct call on original file
    obj = _try_dssp(file_path)
    if obj:
        return obj

    # Strategy 2: text-stripped file + CRYST1 (avoids PDBIO element issues)
    try:
        clean = _strip_for_dssp(file_path)
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
            f.write(clean); tmp_clean = f.name
        try:
            obj = _try_dssp(tmp_clean)
            if obj:
                return obj
        finally:
            try: os.unlink(tmp_clean)
            except OSError: pass
    except Exception:
        pass

    # Strategy 3: shell subprocess on stripped text
    try:
        clean = _strip_for_dssp(file_path)
        dssp_text = _run_binary_dssp(clean)
        if dssp_text:
            obj = _dssp_from_text(dssp_text)
            if obj:
                return obj
    except Exception:
        pass

    print(f"[DSSP/omega] {stem}: all strategies failed — Ω-loops will be empty",
          file=sys.stderr)
    return None


# ── FreeSASA: computed ONCE per file, sliced per loop ─────────────────────
#
# Previous implementation called freesasa.calcBioPDB(structure, ...) inside
# _loop_asa, i.e. for every candidate Ω-loop window — turning a structure-
# level O(N) computation into O(N · L) work where L = number of loops. On
# real proteins (~300 residues, ~20-30 Ω-loops, ~0.5-1 s per FreeSASA pass)
# that became 15-30 s per file and made the web UI appear to hang.
# Now we compute residueAreas() exactly once and slice it per window.

def _compute_sasa(structure):
    """Run FreeSASA on the whole complex. Returns the result object or None."""
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


def _loop_asa(sasa_result, chain_id, start_res, end_res):
    """
    Return summed SASA (Å²) for residues [start_res, end_res] in chain_id
    from a pre-computed FreeSASA result. None if result is unavailable.
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


def _exposure_index(asa, radius):
    import math
    if asa is None or radius is None or radius == 0:
        return None
    sphere_sa = 4 * math.pi * radius ** 2
    return round(asa / sphere_sa, 4) if sphere_sa > 0 else None


# ── record builder ─────────────────────────────────────────────────────────

def _make_record(pdb_id, loop_number, chain_id, window, cas,
                 ss_codes, ss_before, ss_after,
                 endpoint_dist, asa, radius, area, volume,
                 cx, cy, cz, diameter):
    seq = "".join(_three_to_one(r.get_resname()) for r in window)
    start_res = window[0].id[1]
    end_res   = window[-1].id[1]
    exposure  = _exposure_index(asa, radius)

    def _r(v, n=3):
        return round(v, n) if v is not None else None

    return {
        # identity
        "pdb_id":               pdb_id,
        "loop_number":          loop_number,
        "Chain":                chain_id,
        "Start_Position":       start_res,
        "End_Position":         end_res,
        "range":                f"{start_res}-{end_res}",   # for viewer highlight
        # sequence
        "Amino_Acid_Seq":       seq,
        # secondary structure
        "Sec_Structure_before": ss_before,
        "Sec_Structure":        "".join(ss_codes),
        "Sec_Structure_After":  ss_after,
        # geometry
        "distance_x":           _r(endpoint_dist),
        "x1":                   _r(cas[0][0]),
        "y1":                   _r(cas[0][1]),
        "z1":                   _r(cas[0][2]),
        "x2":                   _r(cas[-1][0]),
        "y2":                   _r(cas[-1][1]),
        "z2":                   _r(cas[-1][2]),
        "Center_x":             _r(cx),
        "Center_y":             _r(cy),
        "Center_z":             _r(cz),
        "Diameter":             _r(diameter),
        "radius":               _r(radius),
        "Area":                 _r(area),
        "Volume":               _r(volume),
        # solvent exposure
        "total_asa":            _r(asa, 2) if asa is not None else None,
        "exposure_index":       exposure,
        # start/end residue names
        "Amino_acid1":          window[0].get_resname(),
        "Amino_acid2":          window[-1].get_resname(),
    }


# ── main scanner ───────────────────────────────────────────────────────────

def find_omega_loops(file_path):
    """
    Scan a single structure file for Ω-loops.

    Returns
    -------
    loops      : list of dicts
    status_str : "ok:N" | "dssp-failed" | "no-loops" | "error:..."
    """
    try:
        pdb_id    = os.path.splitext(os.path.basename(file_path))[0]
        parser    = _parser_for(file_path)
        structure = parser.get_structure(pdb_id, file_path)
        model     = structure[0]

        # ── DSSP ──────────────────────────────────────────────────────────
        dssp = _run_dssp(model, file_path)
        if dssp is None:
            return [], "dssp-failed"

        # ── FreeSASA: compute ONCE per file (was per-loop — major hang fix) ─
        sasa_result = _compute_sasa(structure)

        loops        = []
        loop_counter = 0

        for chain in model:
            residues = [r for r in chain if r.id[0] == " "]   # standard only
            n        = len(residues)

            for loop_length in range(MIN_LOOP_LENGTH, MAX_LOOP_LENGTH + 1):
                # S/T threshold scales with window size:
                # MIN_LOOP_LENGTH=6 requires 5 S/T (~83%). Keep same ratio.
                min_st_count = loop_length - 1

                for i in range(n - loop_length + 1):
                    window = residues[i : i + loop_length]
                    cas    = [_ca_coord(r) for r in window]
                    if any(c is None for c in cas):
                        continue

                    # Skip windows that contain UNK or non-standard residues
                    if any(r.get_resname().strip() not in {
                               'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY',
                               'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
                               'THR','TRP','TYR','VAL'} for r in window):
                        continue

                    # DSSP codes for this window
                    try:
                        ss_codes = [dssp[(chain.id, r.id)][2] for r in window]
                    except Exception:
                        continue

                    # Ω-loop filter
                    if sum(s in ST_CODES for s in ss_codes) < min_st_count:
                        continue
                    if not all(s in ALLOWED_SS for s in ss_codes):
                        continue

                    endpoint_dist = _dist(cas[0], cas[-1])
                    if endpoint_dist >= MAX_END_DIST:
                        continue

                    max_pair = max(
                        _dist(cas[a], cas[b])
                        for a in range(loop_length)
                        for b in range(a + 1, loop_length)
                    )
                    if endpoint_dist >= (2 / 3) * max_pair:
                        continue

                    # Flanking SS codes
                    ss_before = ""
                    if i > 0:
                        key = (chain.id, residues[i - 1].id)
                        if key in dssp:
                            ss_before = dssp[key][2]

                    ss_after = ""
                    if (i + loop_length) < n:
                        key = (chain.id, residues[i + loop_length].id)
                        if key in dssp:
                            ss_after = dssp[key][2]

                    # Geometry
                    radius, area, volume, cx, cy, cz, diameter = _sphere_geometry(cas)

                    # Solvent exposure — sliced from the one precomputed result
                    start_res_num = window[0].id[1]
                    end_res_num   = window[-1].id[1]
                    asa = _loop_asa(sasa_result, chain.id, start_res_num, end_res_num)

                    loop_counter += 1
                    loops.append(_make_record(
                        pdb_id, loop_counter, chain.id, window, cas,
                        ss_codes, ss_before, ss_after,
                        endpoint_dist, asa,
                        radius, area, volume, cx, cy, cz, diameter,
                    ))

        status = f"ok:{len(loops)}" if loops else "no-loops"
        return loops, status

    except Exception as e:
        traceback.print_exc()
        return [], f"error:{type(e).__name__}: {e}"


# ── CLI ────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import sys, json
    if len(sys.argv) != 2:
        sys.exit("usage: python omega_backend.py <structure.pdb|.cif|.mmcif>")
    loops, status = find_omega_loops(sys.argv[1])
    print(f"Status: {status}", file=sys.stderr)
    print(json.dumps(loops, indent=2))
