# Motif Scanner

Detection of sequence and structural motifs in proteins, from the Structural
Bioinformatics Lab, Department of Bioscience and Bioengineering, Indian Institute of
Technology Jodhpur.

**Motif Scanner is a web service. Use it at
<https://structbioinfo.iitj.ac.in/motif-scanner>** — nothing needs to be installed.
This repository holds the source code of the deployed service, published so that the
methods behind the results can be read and reviewed.

## What it does

Given protein structures (`.pdb`, `.cif`, `.mmcif`) or plain sequences (`.fasta`,
`.txt`), Motif Scanner reports:

* **Sequence motifs** of the X***X family (AXXXA, GXXXG, and the rest), each with its
  secondary-structure and interface context;
* **Ω-loops**, identified from DSSP secondary structure and loop geometry;
* **α-, β- and γ-turns**, classified against published geometric criteria.

Structure files receive the full analysis; sequence files are scanned for sequence
motifs alone. Results are shown in the browser and can be downloaded as CSV or Excel,
emailed, or taken as a ZIP of every table.

## Repository contents

| File | What it is |
|------|------------|
| `app.py` | Flask application: the web interface, analysis routes, exports and email |
| `motif_backend.py` | Sequence-motif detection, DSSP secondary structure, interface context |
| `omega_backend.py` | Ω-loop detection |
| `turn_backend.py` | α-, β- and γ-turn classification |
| `check_interface.py` | Helper for interface checks |
| `NHN_gamma_turn_criteria.txt` | Published geometric criteria used for γ-turns |
| `demo.pdb`, `1aap.pdb` | Example structures used by the service's demo button |
| `requirements.txt` | Python dependencies of the deployed environment |
| `.env.example`, `setup_env.sh` | Shape of the server-side mail configuration (no credentials) |
| `launch.sh`, `launch.bat` | Development launchers kept from the original tool |
| `README.txt` | The original user-facing notes |

## How the service runs

The deployed instance runs as part of the StructBio suite at IIT Jodhpur: one process
serving Motif Scanner together with the lab's other tools, behind Apache with HTTPS.
The analysis environment is Python 3.12 with the packages in `requirements.txt`, plus
two external programs:

* **DSSP** (`mkdssp`) — required for Ω-loops and the secondary-structure column;
* **FreeSASA** — used for interface accessibility, with a dependency-free fallback
  detector when it is unavailable.

Mail settings (used for the "email results" feature) are supplied to the server through
environment variables; see `.env.example` for the variable names. **No credentials are
stored in this repository or in the code.**

## Citation

If you use Motif Scanner in your work, please cite the accompanying paper
(details to be added on publication).

## Contact

structbioiitj@gmail.com — Structural Bioinformatics Lab, IIT Jodhpur.
