"""
motif_backend.py
----------------
Finds X***X sequence motifs in:
  - Protein structures (.pdb / .cif / .mmcif)  via Bio.PDB
  - Plain-text sequences (.txt)                 via find_motifs_in_sequence()
    Accepts raw one-letter sequences and FASTA (single or multi-record).

Changes from the original, and why:
  * Works on BOTH .pdb and .cif/.mmcif (auto-detected) instead of PDB only.
  * Returns plain Python lists instead of writing two hard-coded Excel paths.
  * The motif scan is reset PER CHAIN. The original kept one shared buffer
    across every chain in the model, so a motif could silently span a chain
    break and report a meaningless residue range / wrong chain. That was a
    bug; this version scans each chain independently. Set reset_per_chain=False
    to reproduce the old behaviour exactly.
  * The seven near-identical copy-pasted blocks (A/G/E/I/L/K/V) collapse into
    one parameterised pass over MOTIF_RESIDUES.
  * New: find_motifs_in_sequence() / parse_text_file() for .txt input.
"""

import os
import re
from collections import OrderedDict

from Bio.PDB import PDBParser, MMCIFParser, PDBIO
import io

# 3-letter -> 1-letter for the 20 standard amino acids.
AA_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Valid one-letter amino acid codes (uppercase).
VALID_AA = set(AA_3TO1.values())

# The motifs the original script hunted for: <residue>***<residue>.
MOTIF_RESIDUES = OrderedDict([
    ("ALA", "A***A"),
    ("GLY", "G***G"),
    ("GLU", "E***E"),
    ("ILE", "I***I"),
    ("LEU", "L***L"),
    ("LYS", "K***K"),
    ("VAL", "V***V"),
])

# The 7 target one-letter codes we care about
TARGET_LETTERS = {seq_type[0] for seq_type in MOTIF_RESIDUES.values()}

# How many residues make up one motif window (first + 3 wildcards + last).
WINDOW = 5

# ΔSASA threshold (Å²) for calling a residue interface-participating.
# Residues that lose more than this much SASA upon complex formation are
# considered to be at the chain–chain interface.
IFACE_DSASA_THRESHOLD = 1.0


class UnsupportedFileError(ValueError):
    pass


# ── structure parsing ──────────────────────────────────────────────────────

def aa_3to1(resname):
    """1-letter code for a residue, or None if it is not a standard amino acid."""
    return AA_3TO1.get(resname)


def _parser_for(file_path):
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".pdb":
        return PDBParser(QUIET=True)
    if ext in (".cif", ".mmcif"):
        return MMCIFParser(QUIET=True)
    raise UnsupportedFileError(
        f"Unsupported extension '{ext}'. Use .pdb, .cif or .mmcif."
    )


def parse_structure(file_path, structure_id=None):
    """Parse a PDB or mmCIF file into a Bio.PDB Structure."""
    if structure_id is None:
        structure_id = os.path.splitext(os.path.basename(file_path))[0]
    parser = _parser_for(file_path)
    return parser.get_structure(structure_id, file_path)


def get_secondary_structure(file_path, structure):
    """
    Run DSSP on a structure and return (ss_map, status).

      ss_map : {(chain_id, res_seq_int): ss_code}   (ss_code in H E B T S G I or '-')
      status : "ok:N (via <method>)"        — DSSP ran, N residues mapped
               "dssp-failed: <details>"     — every attempt failed (details = why)
               "error: <details>"           — DSSP could not even be imported

    Why this is different from the old version (and from omega_backend):
      * The old code wrapped everything in `except Exception: return {}`, so
        when DSSP failed for ANY reason the SS column just went blank with no
        explanation. omega_backend "worked" on the same files only because it
        carried a PDBIO temp-file fallback that motif_backend never had — so a
        direct DSSP call that failed on a quirky / non-PDB-HEADER / mmCIF input
        was silently fatal here but recoverable there.
      * This version tries the same ladder of strategies omega uses (direct
        call, explicit mkdssp/dssp, and a PDBIO-normalised temp .pdb), records
        the exact error from each failed attempt, and hands that string back to
        the caller so it can be surfaced in the UI warnings list.
    """
    import sys, io, os as _os, tempfile, shutil

    try:
        from Bio.PDB.DSSP import DSSP
    except Exception as e:
        return {}, f"error: cannot import Bio.PDB.DSSP ({type(e).__name__}: {e})"

    model = structure[0]

    # Find DSSP binary — local folder first so Windows users can place
    # dssp.exe / mkdssp.exe alongside app.py without adding to PATH
    _here = _os.path.dirname(_os.path.abspath(__file__))
    _local_binary = None
    for _name in ("mkdssp.exe", "dssp.exe", "mkdssp", "dssp"):
        _p = _os.path.join(_here, _name)
        if _os.path.isfile(_p):
            _local_binary = _p
            break
    if _local_binary is None:
        for _name in ("mkdssp", "dssp"):
            _p = shutil.which(_name)
            if _p:
                _local_binary = _p
                break

    def _via_temp_pdb(binary=None):
        from Bio.PDB import PDBIO
        buf = io.StringIO()
        writer = PDBIO()
        writer.set_structure(structure)
        writer.save(buf)
        with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as tf:
            tf.write(buf.getvalue())
            tf_name = tf.name
        try:
            return DSSP(model, tf_name, **({'dssp': binary} if binary else {}))
        finally:
            try:
                _os.unlink(tf_name)
            except OSError:
                pass

    attempts = []
    if _local_binary:
        attempts += [
            ("local-binary/direct",   lambda b=_local_binary: DSSP(model, file_path, dssp=b)),
            ("local-binary/temp-pdb", lambda b=_local_binary: _via_temp_pdb(b)),
        ]
    attempts += [
        ("direct/default", lambda: DSSP(model, file_path)),
        ("direct/mkdssp",  lambda: DSSP(model, file_path, dssp="mkdssp")),
        ("direct/dssp",    lambda: DSSP(model, file_path, dssp="dssp")),
        ("temp-pdb",       _via_temp_pdb),
    ]

    dssp   = None
    used   = None
    errors = []
    for label, call in attempts:
        try:
            cand = call()
            if cand is not None and len(cand.keys()) > 0:
                dssp, used = cand, label
                break
            errors.append(f"{label}: ran but produced 0 residues")
        except Exception as e:
            errors.append(f"{label}: {type(e).__name__}: {e}")

    if dssp is None:
        detail = " | ".join(errors) if errors else "unknown failure"
        print(f"[DSSP] {_os.path.basename(file_path)} — all attempts failed:\n  "
              + "\n  ".join(errors), file=sys.stderr)
        return {}, f"dssp-failed: {detail}"

    ss_map   = {}
    bad_keys = 0
    for key in dssp.keys():
        try:
            chain_id = key[0]
            res_seq  = key[1][1]          # key[1] = (hetflag, resseq, icode)
            ss       = dssp[key][2]
            ss_map[(chain_id, res_seq)] = ss if ss not in (" ", "") else "-"
        except Exception:
            bad_keys += 1
            continue

    status = f"ok:{len(ss_map)} (via {used})"
    if bad_keys:
        status += f"; {bad_keys} unparseable DSSP keys"
    return ss_map, status


def _interface_by_contact(structure, cutoff=5.0):
    """
    Dependency-free interface detector — no FreeSASA, no compiled wheels.

    A residue is "interface-participating" if any of its atoms lies within
    `cutoff` Å of an atom belonging to a DIFFERENT chain. This is a direct
    geometric definition of a chain–chain interface and uses only Bio.PDB's
    NeighborSearch (KD-tree). It is the fallback for environments where
    FreeSASA cannot be installed (e.g. Python 3.14, which has no wheel and
    whose source build fails).

    Returns (iface_set, status) with the same contract as the FreeSASA path:
      iface_set : set of (chain_id, res_seq_int)
      status    : "ok:N (...)" | "single-chain" | "error: ..."
    """
    import sys
    try:
        from Bio.PDB import NeighborSearch

        model  = structure[0]
        chains = [c for c in model.get_chains()]
        if len(chains) < 2:
            return set(), "single-chain"

        atoms = [a for c in chains for r in c if r.id[0] == " " for a in r]
        if not atoms:
            return set(), "error: no standard-residue atoms found"

        ns    = NeighborSearch(atoms)
        iface = set()
        for a, b in ns.search_all(cutoff, level="A"):
            ch_a = a.get_parent().get_parent().id
            ch_b = b.get_parent().get_parent().id
            if ch_a != ch_b:
                for at in (a, b):
                    res = at.get_parent()
                    if res.id[0] == " ":
                        iface.add((at.get_parent().get_parent().id, res.id[1]))
        chain_ids = [c.id for c in chains]
        print(f"[Interface/contact] chains {chain_ids}: "
              f"{len(iface)} residues within {cutoff} A of another chain",
              file=sys.stderr)
        return iface, f"ok:{len(iface)} (inter-chain contact < {cutoff} A)"
    except Exception as e:
        import traceback
        traceback.print_exc(file=sys.stderr)
        return set(), f"error: contact method failed ({type(e).__name__}: {e})"


def get_interface_residues(file_path, structure, threshold=None):
    """
    Identify residues participating in chain-chain interfaces.

    Primary method  : FreeSASA dSASA (isolated chain - complex > threshold).
    Fallback method : inter-chain atomic contact (Bio.PDB only) - used whenever
                      FreeSASA is not importable OR its computation fails. This
                      keeps the feature working on Python 3.14, where FreeSASA
                      has no wheel and won't build from source.

    The returned status string names which method was used so the UI can be
    honest about it.

    Strategy:
      1. Compute SASA for the full complex via freesasa.calcBioPDB.
      2. For each isolated chain compute SASA:
           .pdb  → freesasa.structureArray  (fastest; reads file directly)
           .cif  → PDBIO round-trip per chain + calcBioPDB
      3. Residue is interface-participating if ΔSASA (isolated − complex) > threshold.

    Returns (iface_set, status_str):
      iface_set  : set of (chain_id, res_seq_num_int); empty on any error/skip
      status_str : human-readable outcome, e.g. "ok:12" or "error:..." or "single-chain"
    Never raises.
    """
    import sys, traceback as _tb

    if threshold is None:
        threshold = IFACE_DSASA_THRESHOLD

    # ── FreeSASA import ───────────────────────────────────────────────────
    # If FreeSASA is unavailable (the Python-3.14 case), don't give up — fall
    # back to the dependency-free inter-chain contact detector.
    try:
        import freesasa
        freesasa.setVerbosity(freesasa.silent)
    except ImportError:
        iface, status = _interface_by_contact(structure)
        if status.startswith("ok"):
            status += " [freesasa unavailable - contact fallback]"
        return iface, status

    try:
        model  = structure[0]
        chains = list(model.get_chains())
        chain_ids = [c.id for c in chains]

        if len(chains) < 2:
            return set(), "single-chain"

        print(f"[FreeSASA] {os.path.basename(file_path)}: "
              f"{len(chains)} chains {chain_ids}", file=sys.stderr)

        # ── 1. Full-complex SASA ──────────────────────────────────────────
        try:
            ret = freesasa.calcBioPDB(
                structure,
                options={"hetatm": False, "hydrogen": False, "skip-unknown": True},
            )
            result_complex = ret[0]
            areas_complex  = result_complex.residueAreas()
            n_complex = sum(len(v) for v in areas_complex.values())
            print(f"[FreeSASA] complex: {n_complex} residues across "
                  f"chains {list(areas_complex.keys())}", file=sys.stderr)
        except Exception as e:
            _tb.print_exc(file=sys.stderr)
            print(f"[FreeSASA] complex SASA failed ({type(e).__name__}: {e}); "
                  f"falling back to contact-based interface detection", file=sys.stderr)
            iface, status = _interface_by_contact(structure)
            return iface, status + " [freesasa-failed - contact fallback]"

        # ── 2. Per-chain isolated SASA ────────────────────────────────────
        ext = os.path.splitext(file_path)[1].lower()
        areas_isolated = {}   # orig_chain_id -> {resnum_str: float}

        if ext == ".pdb":
            # structureArray reads the file directly; chainLabel(0) gives original ID
            try:
                chain_structs = freesasa.structureArray(
                    file_path,
                    options={
                        "hetatm": False, "hydrogen": False,
                        "separate-chains": True, "separate-models": False,
                    },
                )
                print(f"[FreeSASA] structureArray: {len(chain_structs)} chain structures",
                      file=sys.stderr)
                for cs in chain_structs:
                    orig_cid  = cs.chainLabel(0)
                    result_iso = freesasa.calc(cs)
                    iso_areas  = result_iso.residueAreas()
                    for rdict in iso_areas.values():
                        areas_isolated[orig_cid] = {
                            rnum: ra.total for rnum, ra in rdict.items()
                        }
                    print(f"[FreeSASA]   chain {repr(orig_cid)}: "
                          f"{len(areas_isolated.get(orig_cid, {}))} residues",
                          file=sys.stderr)
            except Exception as e:
                _tb.print_exc(file=sys.stderr)
                return set(), f"error in structureArray: {type(e).__name__}: {e}"

        else:
            # mmCIF path: PDBIO round-trip per chain
            import tempfile
            from Bio.PDB import Structure as _S, Model as _M, PDBIO as _IO
            for chain in chains:
                cid = chain.id
                tmp_s = _S.Structure("tmp")
                tmp_m = _M.Model(0)
                tmp_s.add(tmp_m)
                tmp_m.add(chain.copy())
                buf = io.StringIO()
                pdbio = _IO()
                pdbio.set_structure(tmp_s)
                pdbio.save(buf)
                pdb_str = buf.getvalue()
                with tempfile.NamedTemporaryFile(
                        suffix=".pdb", mode="w", delete=False) as tf:
                    tf.write(pdb_str)
                    tf_name = tf.name
                try:
                    chain_struct = PDBParser(QUIET=True).get_structure("c", tf_name)
                    ret_iso = freesasa.calcBioPDB(
                        chain_struct,
                        options={"hetatm": False, "hydrogen": False,
                                 "skip-unknown": True},
                    )
                    iso_areas = ret_iso[0].residueAreas()
                    for rdict in iso_areas.values():
                        areas_isolated[cid] = {
                            rnum: ra.total for rnum, ra in rdict.items()
                        }
                    print(f"[FreeSASA]   chain {repr(cid)}: "
                          f"{len(areas_isolated.get(cid, {}))} residues",
                          file=sys.stderr)
                except Exception as e:
                    print(f"[FreeSASA]   chain {repr(cid)}: FAILED — "
                          f"{type(e).__name__}: {e}", file=sys.stderr)
                    _tb.print_exc(file=sys.stderr)
                    areas_isolated[cid] = {}
                finally:
                    os.unlink(tf_name)

        # ── 3. ΔSASA → interface residues ─────────────────────────────────
        iface = set()
        unmatched_chains = []
        for orig_cid, iso_map in areas_isolated.items():
            complex_chain = areas_complex.get(orig_cid, {})
            if not complex_chain:
                unmatched_chains.append(orig_cid)
                print(f"[FreeSASA] WARNING: chain {repr(orig_cid)} in isolated "
                      f"but not in complex areas — chain ID mismatch?",
                      file=sys.stderr)
                continue
            matched = 0
            for resnum_str, iso_sasa in iso_map.items():
                comp_area = complex_chain.get(resnum_str)
                if comp_area is None:
                    continue
                matched += 1
                delta = iso_sasa - comp_area.total
                if delta > threshold:
                    try:
                        iface.add((orig_cid, int(resnum_str.strip())))
                    except ValueError:
                        pass   # insertion codes — can't map to Bio.PDB integer seq
            print(f"[FreeSASA]   chain {repr(orig_cid)}: {matched} residues matched, "
                  f"{sum(1 for c,r in iface if c==orig_cid)} interface",
                  file=sys.stderr)

        status = f"ok:{len(iface)} (freesasa dSASA > {threshold} A^2)"
        if unmatched_chains:
            status += f" (unmatched chains: {unmatched_chains})"
        print(f"[FreeSASA] done: {len(iface)} interface residues total", file=sys.stderr)
        return iface, status

    except Exception as e:
        _tb.print_exc(file=sys.stderr)
        # FreeSASA was importable but the computation blew up — rather than
        # returning a hard error, fall back to the geometric contact method.
        iface, cstatus = _interface_by_contact(structure)
        if cstatus.startswith("ok"):
            cstatus += f" [freesasa failed: {type(e).__name__}; contact fallback]"
        return iface, cstatus


def _extract_ss_for_motif(ss_map, chain_id, range_str):
    """Return 5-char SS string for a motif window, '-' for missing residues."""
    start, end = map(int, range_str.split('-'))
    return ''.join(ss_map.get((chain_id, r), '-') for r in range(start, end + 1))


def _format_motif(window):
    """Collapse a 5-residue window into one output record."""
    return {
        "pdb_id":   window[0]["pdb_id"],
        "sn":       "".join(r["sn"] for r in window),
        "range":    f'{window[0]["id"]}-{window[-1]["id"]}',
        "chain_id": window[0]["chain_id"],
        "seq_type": window[0]["seq_type"],
        "source":   "structure",
        "ss":       "",    # filled in by find_motifs after DSSP runs
        "iface":    "",    # filled in by find_motifs after FreeSASA runs
    }


def _scan(residues, target_resname, seq_type):
    """
    Walk a sequence of residues and emit every <target>***<target> window.

    A motif must cover five CONSECUTIVE residue numbers (N, N+1, N+2, N+3, N+4).
    Without this check the scanner happily steps over unmodelled gaps or
    intervening hetero/non-standard residues, producing impossible "5-residue"
    motifs that span 17 positions in the structure (e.g. GRTTG at 140-156 in
    2ALP). Reset the buffer whenever residue numbering breaks, or whenever
    we encounter a non-standard residue mid-window.
    """
    def _record(residue, one):
        return {
            "pdb_id":   residue.full_id[0],
            "id":       residue.id[1],
            "chain_id": residue.full_id[2],
            "sn":       one,
            "seq_type": seq_type,
        }

    found = []
    buf = []
    for residue in residues:
        one = aa_3to1(residue.resname)
        is_target = residue.resname == target_resname
        res_num   = residue.id[1]

        # Anything that isn't a standard amino acid breaks the window.
        # Hetero flag (residue.id[0] != ' ') also breaks it — waters, ligands,
        # modified residues etc. must not be silently skipped over.
        if one is None or residue.id[0] != ' ':
            buf = []
            continue

        # Enforce strictly consecutive residue numbering.
        if buf and res_num != buf[-1]["id"] + 1:
            buf = []

        if is_target:
            buf.append(_record(residue, one))
            if len(buf) == WINDOW:
                found.append(_format_motif(buf))
                buf = []
        else:
            if buf:
                buf.append(_record(residue, one))
                # A non-target residue can never close a window; if it lands
                # in slot 5 the window is invalid, discard and move on.
                if len(buf) == WINDOW:
                    buf = []
    return found


def _extract_iface_for_motif(iface_set, status, chain_id, range_str):
    """
    Return a 5-char interface string for a motif window.
    'I' = interface-participating, '.' = not, '' = data unavailable.
    """
    # Propagate unavailability: empty set AND non-ok status → no data
    if not iface_set and not status.startswith("ok"):
        return ""
    start, end = map(int, range_str.split('-'))
    return ''.join(
        'I' if (chain_id, r) in iface_set else '.'
        for r in range(start, end + 1)
    )


def find_motifs(file_path, structure_id=None, reset_per_chain=True):
    """
    Find every motif in MOTIF_RESIDUES for a single structure file.
    Annotates each motif with:
      ss        : 5-char secondary structure string (DSSP, best-effort)
      iface     : 5-char interface string ('I' = interface, '.' = not, '' = unavailable)

    Returns (motifs, counts, sasa_status, ss_status):
      motifs      : list of dicts
      counts      : list of dicts
      sasa_status : interface-detection outcome ("ok:N ...", "single-chain", "error:...")
      ss_status   : DSSP outcome ("ok:N (via ...)", "dssp-failed: ...", "error:...")
    """
    structure = parse_structure(file_path, structure_id)
    model = structure[0]

    motifs = []
    for resname, seq_type in MOTIF_RESIDUES.items():
        if reset_per_chain:
            for chain in model:
                motifs.extend(_scan(list(chain), resname, seq_type))
        else:
            all_residues = [res for chain in model for res in chain]
            motifs.extend(_scan(all_residues, resname, seq_type))

    # Annotate with secondary structure. ss_status carries the reason if DSSP
    # failed, so the caller can surface it instead of showing silent dashes.
    ss_map, ss_status = get_secondary_structure(file_path, structure)
    for m in motifs:
        m["ss"] = _extract_ss_for_motif(ss_map, m["chain_id"], m["range"])

    # Annotate with interface participation
    iface_set, sasa_status = get_interface_residues(file_path, structure)
    for m in motifs:
        m["iface"] = _extract_iface_for_motif(
            iface_set, sasa_status, m["chain_id"], m["range"]
        )

    return motifs, summarise_counts(motifs), sasa_status, ss_status


# ── sequence parsing ───────────────────────────────────────────────────────

def _clean_sequence(raw):
    """Strip whitespace and digits (common in numbered FASTA exports)."""
    return re.sub(r'[\s\d]', '', raw).upper()


def parse_text_file(text, source_name):
    """
    Parse a plain-text file containing one or more sequences.

    Accepts:
      - Raw sequence  : one or more lines of one-letter codes, no '>' present.
      - FASTA         : one or more records starting with '>header'.

    Returns a list of (record_name, sequence_string) tuples.
    Unknown / non-AA characters are silently dropped per-position
    (a warning list is also returned as second element).
    """
    records = []
    warnings = []

    if '>' in text:
        # FASTA
        current_name = None
        current_seq  = []
        for line in text.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name is not None:
                    records.append((current_name, _clean_sequence(''.join(current_seq))))
                # Use everything after '>' up to first space as the record id
                current_name = line[1:].split()[0] or source_name
                current_seq  = []
            else:
                current_seq.append(line)
        if current_name is not None:
            records.append((current_name, _clean_sequence(''.join(current_seq))))
    else:
        # Raw sequence — treat entire file as one record
        seq = _clean_sequence(text)
        records.append((source_name, seq))

    # Validate and strip non-AA characters, collect warnings
    clean_records = []
    for name, seq in records:
        bad = set(c for c in seq if c not in VALID_AA)
        if bad:
            warnings.append(
                f"'{name}': ignored non-standard characters: {', '.join(sorted(bad))}"
            )
        clean_seq = ''.join(c for c in seq if c in VALID_AA)
        if clean_seq:
            clean_records.append((name, clean_seq))

    return clean_records, warnings


def _scan_sequence(seq, target_letter, seq_type, seq_name):
    """
    Walk a 1-letter string and emit every <target>***<target> window.
    Positions are 1-indexed to match structure convention.
    Non-target non-standard characters have already been stripped by parse_text_file.
    """
    found = []
    buf   = []   # list of (pos, letter)

    for i, aa in enumerate(seq):
        is_target = (aa == target_letter)
        is_valid  = (aa in VALID_AA)

        if is_target:
            buf.append((i + 1, aa))
            if len(buf) == WINDOW:
                sn    = ''.join(aa for _, aa in buf)
                start = buf[0][0]
                end   = buf[-1][0]
                found.append({
                    "pdb_id":   seq_name,
                    "sn":       sn,
                    "range":    f"{start}-{end}",
                    "chain_id": "—",
                    "seq_type": seq_type,
                    "source":   "sequence",
                })
                buf = []
        else:
            if buf and is_valid:
                buf.append((i + 1, aa))
            if len(buf) == WINDOW:
                buf = []

    return found


def find_motifs_in_sequence(text, source_name):
    """
    Find every X***X motif in a plain-text sequence file (raw or FASTA).

    Returns (motifs, counts, warnings).
    """
    records, warnings = parse_text_file(text, source_name)
    motifs = []
    for rec_name, seq in records:
        for resname, seq_type in MOTIF_RESIDUES.items():
            target_letter = seq_type[0]   # e.g. 'A' from 'A***A'
            motifs.extend(_scan_sequence(seq, target_letter, seq_type, rec_name))
    return motifs, summarise_counts(motifs), warnings


# ── shared helpers ─────────────────────────────────────────────────────────

def summarise_counts(motifs):
    """Count motifs per (pdb_id, seq_type), preserving first-seen order."""
    counts = []
    index  = {}
    for m in motifs:
        key = (m["pdb_id"], m["seq_type"])
        if key in index:
            counts[index[key]]["count"] += 1
        else:
            index[key] = len(counts)
            counts.append({
                "pdb_id":   m["pdb_id"],
                "seq_type": m["seq_type"],
                "count":    1,
            })
    return counts


def export_excel(motifs, counts, motifs_path, counts_path):
    """Optional: reproduce the original two-Excel-file output."""
    import pandas as pd
    pd.DataFrame(motifs).to_excel(motifs_path, index=False)
    pd.DataFrame(counts).to_excel(counts_path, index=False)


if __name__ == "__main__":
    import sys, json
    if len(sys.argv) != 2:
        sys.exit("usage: python motif_backend.py <structure.pdb|.cif|.txt>")
    path = sys.argv[1]
    ext  = os.path.splitext(path)[1].lower()
    if ext == ".txt":
        text = open(path).read()
        m, c, w = find_motifs_in_sequence(text, os.path.splitext(os.path.basename(path))[0])
        if w:
            print("Warnings:", w, file=sys.stderr)
    else:
        m, c, sasa_status, ss_status = find_motifs(path)
        print(f"SS:   {ss_status}", file=sys.stderr)
        print(f"SASA: {sasa_status}", file=sys.stderr)
    print(json.dumps({"motifs": m, "counts": c}, indent=2))
