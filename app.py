"""
app.py -- web frontend for motif_backend.find_motifs + omega_backend.find_omega_loops.

Run:
    pip install flask biopython freesasa
    python app.py
    # open http://127.0.0.1:5000

Supports:
  .pdb / .cif / .mmcif  -- structure files, shown with 3D viewer
                           Tab 1: X***X sequence motifs
                           Tab 2: Ω-loops (auto DSSP + FreeSASA)
  .txt / .fasta ...     -- plain-text sequences, motifs only

Structure and sequence results are rendered in separate sections.
CSV export is done client-side.
"""

import csv
import io
import os
import smtplib
import tempfile
import traceback
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

try:
    from dotenv import load_dotenv
    load_dotenv()
except ImportError:
    print("WARNING: python-dotenv not installed in this environment — "
          ".env will not be loaded, SMTP_USER/SMTP_PASS will be empty.")

from flask import Flask, request, jsonify, render_template_string, Response, send_file
from werkzeug.utils import secure_filename

from motif_backend import (
    find_motifs, find_motifs_in_sequence,
    MOTIF_RESIDUES, UnsupportedFileError, summarise_counts,
)
from omega_backend import find_omega_loops
from turn_backend import find_all_turns

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 64 * 1024 * 1024  # 64 MB ceiling
STRUCT_EXT = {".pdb", ".cif", ".mmcif"}
SEQ_EXT    = {".txt", ".fasta", ".fa", ".faa"}
ALLOWED    = STRUCT_EXT | SEQ_EXT

PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>motif scanner · X***X</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:ital,opsz,wght@0,9..144,400;0,9..144,600;1,9..144,500&family=IBM+Plex+Mono:wght@400;500;600&family=IBM+Plex+Sans:wght@400;500;600&display=swap" rel="stylesheet">
<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>
<style>
  :root{
    --bg:#f7f6f2; --panel:#ffffff; --panel-2:#f0efe9; --line:#e0ddd6;
    --ink:#1a1a1a; --ink-dim:#5a5a5a; --ink-faint:#9a9a9a;
    --accent:#0e7c6b; --accent-seq:#7c3aed;
    --A:#b8860b; --G:#2d7a3a; --E:#c0392b; --I:#1a6b8a;
    --L:#6b3a9a; --K:#1a4f9a; --V:#9a1a6b;
    --mono:'IBM Plex Mono',ui-monospace,monospace;
    --sans:'IBM Plex Sans',system-ui,sans-serif;
    --serif:'Playfair Display',serif;
  }
  [data-theme="dark"]{
    --bg:#0f0f0f; --panel:#1a1a1a; --panel-2:#141414; --line:#2a2a2a;
    --ink:#f0ede6; --ink-dim:#a0a0a0; --ink-faint:#555555;
    --accent:#1ab899; --accent-seq:#9d6eff;
    --A:#e6b450; --G:#7bd88f; --E:#ff6f6f; --I:#56c8e0;
    --L:#c09cff; --K:#6aadff; --V:#ff8de0;
  }

  /* ── dark mode toggle ── */
  #theme-toggle{
    position:absolute;top:8px;right:16px;z-index:998;
    background:var(--panel);border:1px solid var(--line);border-radius:999px;
    padding:7px 14px;cursor:pointer;
    display:flex;align-items:center;gap:8px;
    font-family:var(--mono);font-size:12px;color:var(--ink-dim);
    box-shadow:0 2px 8px #00000018;
    transition:background .25s,border-color .25s,color .25s,box-shadow .25s;
  }
  #theme-toggle:hover{color:var(--ink);border-color:var(--ink-dim);
    box-shadow:0 3px 12px #00000028;}
  #theme-toggle .icon{font-size:15px;line-height:1;transition:transform .4s ease}
  #theme-toggle:hover .icon{transform:rotate(20deg)}
  *{box-sizing:border-box}
  html,body{margin:0;min-height:100%}
  body{
    background:var(--bg);
    color:var(--ink); font-family:var(--sans);
    -webkit-font-smoothing:antialiased; line-height:1.5;
  }

  /* ── haemoglobin background canvas ── */
  #hb-bg{
    position:fixed; inset:0; z-index:0;
    pointer-events:none; opacity:0;
    transition:opacity 1.2s ease;
    will-change:transform;
    transform:translateZ(0);
    backface-visibility:hidden;
  }
  #hb-bg.loaded{opacity:0.32}

  .wrap{position:relative;z-index:1;max-width:1440px;margin:0 auto;padding:48px 40px 80px}

  /* contain panels so internal changes don't trigger full-page layout ── */
  .result-section{contain:content}
  li.row,li.omega-row{contain:layout style}

  /* header */
  header{margin-bottom:20px;padding-top:36px;text-align:center;position:relative}
  /* Anchored to the viewport, so it sits under StructBio's suite bar when this
     app is mounted there. --sb-chrome-height is 0 when running standalone. */
  .header-logo{position:absolute;top:calc(8px + var(--sb-chrome-height, 0px));left:16px;height:120px;width:auto;object-fit:contain;z-index:998}
  .motif-logo{display:block;height:260px;width:auto;object-fit:contain;margin:-24px auto 8px}
  h1{font-family:var(--serif);font-weight:800;font-size:clamp(76px,12vw,104px);
     letter-spacing:-.02em;margin:0;line-height:1;color:var(--accent)}
  h1 em{font-style:italic;color:var(--accent)}
  .tag{font-family:var(--mono);font-size:12px;color:var(--ink-dim);
       border:1px solid var(--line);border-radius:999px;padding:4px 11px;
       text-transform:uppercase;letter-spacing:.12em;background:var(--panel)}
  .subtitle{max-width:900px;margin:12px auto 24px;text-align:center;font-size:18px;line-height:1.6;color:var(--ink-dim);}
  .lede{color:var(--ink-dim);max-width:68ch;margin:14px 0 34px;font-size:15px}
  .lede code{font-family:var(--mono);color:var(--ink);background:var(--panel-2);
    padding:1px 6px;border-radius:5px;font-size:13px;border:1px solid var(--line)}

  /* upload area */
  .drop{border:1.5px dashed var(--line);border-radius:16px;background:var(--panel);
    padding:40px 28px;text-align:center;cursor:pointer;transition:.18s;position:relative;
    box-shadow:0 1px 4px #0000000a}
  .drop:hover{border-color:#aaa;background:var(--panel-2)}
  .drop.hot{border-color:var(--accent);background:#e8f5f2;box-shadow:0 0 0 4px #0e7c6b18}
  .drop .glyph{font-family:var(--mono);font-size:30px;color:var(--accent);letter-spacing:.18em}
  .drop h2{font-family:var(--sans);font-weight:600;font-size:17px;margin:14px 0 4px;color:var(--ink)}
  .drop p{color:var(--ink-faint);font-size:13px;margin:0;font-family:var(--mono)}
  .files{display:flex;flex-wrap:wrap;gap:8px;margin-top:14px}
  .chip{font-family:var(--mono);font-size:12px;background:var(--panel);
    border:1px solid var(--line);border-radius:8px;padding:6px 10px;
    display:flex;align-items:center;gap:8px;box-shadow:0 1px 3px #0000000a}
  .chip.chip-seq{border-color:#7c3aed44}
  .chip b{color:var(--ink-dim);cursor:pointer}
  .chip b:hover{color:var(--E)}
  .chip-badge{font-size:10px;padding:1px 6px;border-radius:4px;letter-spacing:.06em;
    text-transform:uppercase;font-weight:600}
  .chip-badge.str{background:#0e7c6b18;color:var(--accent)}
  .chip-badge.seq{background:#7c3aed18;color:var(--accent-seq)}
  .chip-badge.inv{background:#c0392b18;color:var(--E)}
  .chip.chip-invalid{border-color:#c0392b55;background:#c0392b08}
  .chip.chip-invalid b{color:var(--E)}
  .drop-err{font-family:var(--mono);font-size:12px;color:var(--E);
    background:#c0392b08;border:1px solid #c0392b33;border-radius:10px;
    padding:8px 14px;margin-top:10px;display:none;gap:6px;flex-direction:column}
  .drop-err.visible{display:flex}
  /* PDB fetch panel */
  .pdb-fetch-panel{margin:10px 0 0;border:1.5px solid var(--line);border-radius:12px;
    padding:12px 16px 10px;background:var(--panel)}
  .pdb-fetch-label{font-family:var(--sans);font-size:11px;color:var(--ink-faint);
    margin-bottom:8px;letter-spacing:.06em;text-transform:uppercase}
  .pdb-fetch-row{display:flex;gap:8px}
  #pdb-ids-input{flex:1;font-family:var(--mono);font-size:13px;padding:8px 12px;
    border:1.5px solid var(--line);border-radius:8px;background:var(--panel-2);
    color:var(--ink);outline:none;transition:border-color .15s}
  #pdb-ids-input:focus{border-color:var(--accent)}
  #pdb-ids-input::placeholder{color:var(--ink-faint)}
  .fetch-btn{font-family:var(--sans);font-size:13px;font-weight:600;padding:8px 18px;
    border:none;border-radius:8px;background:var(--accent);color:#fff;
    cursor:pointer;transition:filter .15s;white-space:nowrap}
  .fetch-btn:hover{filter:brightness(1.08)}
  .fetch-btn:disabled{opacity:.5;cursor:not-allowed}
  .pdb-fetch-status{font-family:var(--mono);font-size:12px;margin-top:6px;
    min-height:16px;color:var(--ink-dim)}
  .pdb-fetch-status.err{color:var(--E)}
  .pdb-fetch-status.ok{color:#0e7c6b}
  /* email row */
  .email-row{display:flex;align-items:center;gap:10px;margin-top:10px;flex-wrap:wrap}
  .email-row label{font-family:var(--sans);font-size:12px;color:var(--ink-faint);
    white-space:nowrap;letter-spacing:.03em;display:flex;align-items:center;gap:5px}
  .email-sep{font-family:var(--sans);font-size:12px;color:var(--ink-faint)}
  #email-input{font-family:var(--mono);font-size:13px;padding:7px 12px;flex:1;min-width:220px;
    border:1.5px solid var(--line);border-radius:8px;background:var(--panel-2);
    color:var(--ink);outline:none;transition:border-color .15s}
  #email-input:focus{border-color:var(--accent)}
  #email-input::placeholder{color:var(--ink-faint)}
  #email-status{font-family:var(--sans);font-size:13px;color:var(--ink-dim)}
  #email-status.ok{color:#0e7c6b;font-family:var(--mono);font-size:12px;white-space:nowrap}
  #email-status.err{color:var(--E);font-family:var(--mono);font-size:12px;white-space:nowrap}
  .scan-notice{display:block;margin-top:8px;padding:14px 18px;border-radius:10px;
    background:var(--panel-2);border:1.5px solid var(--accent);
    color:var(--ink);line-height:1.7;font-size:13px}
  .scan-notice b{color:var(--accent)}

  /* action bar */
  .bar{display:flex;gap:12px;align-items:center;margin-top:22px;flex-wrap:wrap}
  .scan-row{margin-top:12px}
  button.run{font-family:var(--sans);font-weight:600;font-size:15px;color:#fff;
    background:var(--accent);border:none;border-radius:10px;padding:12px 26px;
    cursor:pointer;transition:.15s;letter-spacing:.01em;box-shadow:0 2px 8px #0e7c6b33}
  button.run:hover{filter:brightness(1.08);transform:translateY(-1px)}
  button.run:disabled{opacity:.4;cursor:not-allowed;transform:none;box-shadow:none}
  .csv{font-family:var(--mono);font-size:13px;color:var(--ink-dim);background:var(--panel);
    border:1px solid var(--line);border-radius:9px;padding:11px 16px;cursor:pointer}
  .csv:hover{color:var(--ink);border-color:#aaa;background:var(--panel-2)}

  /* ── demo button ── */
  .demo-row{margin-top:14px;display:flex;justify-content:center}
  .demo-btn{
    font-family:var(--mono);font-size:12.5px;color:var(--accent);
    background:var(--panel);border:1.5px dashed var(--accent);border-radius:10px;
    padding:10px 20px;cursor:pointer;transition:all .18s;letter-spacing:.02em;
    display:flex;align-items:center;gap:9px;
  }
  .demo-btn:hover{background:#0e7c6b10;border-style:solid;
    box-shadow:0 2px 12px #0e7c6b18;transform:translateY(-1px)}
  .demo-btn .demo-icon{font-size:16px;line-height:1}
  .demo-btn .demo-sub{font-size:10px;color:var(--ink-faint);display:block;margin-top:1px}

  /* ── tour guide card ── */
  #tour-card{
    display:none;position:fixed;bottom:28px;right:28px;z-index:8000;
    background:var(--panel);border:1px solid var(--line);border-radius:16px;
    width:340px;box-shadow:0 12px 40px #00000030;overflow:hidden;
    animation:slideUp .3s ease;
  }
  @keyframes slideUp{from{opacity:0;transform:translateY(16px)}to{opacity:1;transform:translateY(0)}}
  .tour-hdr{
    background:var(--accent);padding:12px 16px;
    display:flex;justify-content:space-between;align-items:flex-start;gap:8px;
  }
  .tour-hdr-left{display:flex;flex-direction:column;gap:2px}
  .tour-eyebrow{font-family:var(--mono);font-size:9px;color:#ffffffaa;
    letter-spacing:.14em;text-transform:uppercase}
  .tour-title{font-family:var(--serif);font-size:15px;color:#fff;
    font-weight:400;font-style:italic;line-height:1.25}
  .tour-step-ind{font-family:var(--mono);font-size:11px;color:#ffffffaa;
    white-space:nowrap;padding-top:2px;flex-shrink:0}
  .tour-body{padding:14px 16px 10px;font-size:13.5px;
    color:var(--ink-dim);line-height:1.65}
  .tour-body b{color:var(--ink);font-weight:600}
  .tour-body code{font-family:var(--mono);font-size:12px;background:var(--panel-2);
    border:1px solid var(--line);border-radius:4px;padding:1px 5px;color:var(--accent)}
  .tour-footer{
    display:flex;align-items:center;justify-content:space-between;
    padding:10px 16px 14px;border-top:1px solid var(--line);
  }
  .tour-nav{display:flex;gap:6px}
  .tour-btn{
    font-family:var(--mono);font-size:12px;cursor:pointer;border-radius:7px;
    padding:6px 13px;border:1px solid var(--line);background:var(--panel-2);
    color:var(--ink-dim);transition:all .12s;
  }
  .tour-btn:hover:not(:disabled){color:var(--ink);border-color:var(--ink-dim)}
  .tour-btn:disabled{opacity:.35;cursor:default}
  .tour-btn.primary{background:var(--accent);color:#fff;border-color:var(--accent)}
  .tour-btn.primary:hover{filter:brightness(1.08)}
  .tour-close-lnk{font-family:var(--mono);font-size:11px;color:var(--ink-faint);
    cursor:pointer;background:none;border:none;padding:0;text-decoration:underline}
  .tour-close-lnk:hover{color:var(--ink-dim)}
  /* progress dots */
  .tour-dots{display:flex;gap:5px;align-items:center}
  .tour-dot{width:6px;height:6px;border-radius:50%;background:var(--line);transition:background .2s}
  .tour-dot.active{background:var(--accent)}

  /* ── tour target highlight ── */
  .tour-highlight{
    outline:2.5px solid var(--accent) !important;outline-offset:5px;
    border-radius:12px !important;
    animation:tour-pulse 1.5s ease-in-out infinite;
  }
  @keyframes tour-pulse{
    0%,100%{box-shadow:0 0 0 0 #0e7c6b30}
    50%{box-shadow:0 0 0 8px #0e7c6b00}
  }
  .info-btn{
    font-family:var(--mono);font-size:12px;color:var(--ink-dim);
    background:var(--panel);border:1px solid var(--line);border-radius:999px;
    padding:5px 14px;cursor:pointer;letter-spacing:.06em;
    transition:color .15s,border-color .15s,background .15s;
  }
  .info-btn:hover{color:var(--accent);border-color:var(--accent);background:var(--panel-2)}

  /* ── definitions overlay ── */
  #defs-overlay{
    display:none;position:fixed;inset:0;z-index:10000;
    background:#00000099;align-items:center;justify-content:center;
    padding:20px;
  }
  .defs-panel{
    background:var(--panel);border:1px solid var(--line);border-radius:16px;
    max-width:800px;width:100%;max-height:88vh;display:flex;flex-direction:column;
    overflow:hidden;box-shadow:0 20px 60px #00000055;
  }
  .defs-hdr{
    display:flex;align-items:center;justify-content:space-between;
    padding:14px 22px;border-bottom:1px solid var(--line);
    background:var(--panel-2);flex-shrink:0;
  }
  .defs-hdr-left{display:flex;flex-direction:column;gap:2px}
  .defs-eyebrow{font-family:var(--mono);font-size:9px;color:var(--ink-faint);
    letter-spacing:.18em;text-transform:uppercase}
  .defs-title{font-family:var(--serif);font-size:19px;font-weight:400;color:var(--ink)}
  .defs-title em{font-style:italic;color:var(--accent)}
  .defs-close-btn{
    background:none;border:1px solid var(--line);border-radius:6px;
    color:var(--ink-dim);cursor:pointer;padding:5px 12px;font-size:14px;transition:all .15s;
  }
  .defs-close-btn:hover{background:#da363318;border-color:#da3633;color:#da3633}
  .defs-body{overflow-y:auto;padding:28px 24px 32px;
    scrollbar-width:thin;scrollbar-color:var(--line) transparent;}
  .defs-section{margin-bottom:40px}
  .defs-section:last-child{margin-bottom:0}
  .defs-sec-hdr{
    display:flex;align-items:center;gap:10px;margin-bottom:14px;
    padding-bottom:10px;border-bottom:1px solid var(--line);
  }
  .defs-tag{font-family:var(--mono);font-size:10px;font-weight:700;
    letter-spacing:.12em;text-transform:uppercase;padding:3px 9px;border-radius:4px;flex-shrink:0}
  .defs-sec-title{font-family:var(--serif);font-size:18px;font-weight:400;color:var(--ink)}
  .defs-prose{color:var(--ink-dim);font-size:13.5px;margin-bottom:14px;line-height:1.6}
  .defs-prose b{color:var(--ink)}
  .defs-note{
    background:#ffa65712;border:1px solid #ffa65733;border-radius:7px;
    padding:10px 14px;font-size:12.5px;color:#c47f00;margin-bottom:14px;line-height:1.5;
  }
  [data-theme="dark"] .defs-note{color:#ffa657}
  .defs-note b{font-weight:700}
  .defs-criteria{background:var(--panel-2);border:1px solid var(--line);
    border-radius:8px;overflow:hidden;margin-bottom:14px}
  .defs-criteria-lbl{font-family:var(--mono);font-size:9.5px;letter-spacing:.14em;
    text-transform:uppercase;color:var(--ink-faint);padding:7px 14px;
    border-bottom:1px solid var(--line)}
  .defs-criterion{display:flex;gap:10px;padding:7px 14px;
    border-bottom:1px solid var(--line);font-size:13px;
    line-height:1.5;color:var(--ink-dim)}
  .defs-criterion:last-child{border-bottom:none}
  .defs-cnum{font-family:var(--mono);font-size:10px;color:var(--ink-faint);
    padding-top:2px;flex-shrink:0;width:18px}
  .defs-criterion b{color:var(--ink)}
  .defs-tbl{width:100%;border-collapse:collapse;font-family:var(--mono);
    font-size:12px;margin-bottom:14px}
  .defs-tbl th{background:var(--panel-2);color:var(--ink-faint);
    font-size:9.5px;letter-spacing:.12em;text-transform:uppercase;
    padding:6px 10px;border:1px solid var(--line);text-align:left;font-weight:700}
  .defs-tbl td{padding:6px 10px;border:1px solid var(--line);
    background:var(--panel);color:var(--ink-dim)}
  .defs-tbl .tn{font-weight:700;color:var(--accent-seq)}
  .defs-tbl .nd{color:var(--ink-faint);font-size:11px}
  .defs-sources{background:var(--panel-2);border:1px solid var(--line);
    border-radius:8px;padding:10px 14px}
  .defs-src-lbl{font-family:var(--mono);font-size:9.5px;letter-spacing:.14em;
    text-transform:uppercase;color:var(--ink-faint);margin-bottom:8px}
  .defs-src{display:flex;gap:8px;padding:5px 0;border-bottom:1px solid var(--line);
    font-size:12px;color:var(--ink-dim);line-height:1.5}
  .defs-src:last-child{border-bottom:none}
  .defs-src-dot{color:var(--ink-faint);font-family:var(--mono);font-size:10px;
    padding-top:2px;flex-shrink:0}
  .defs-src b{color:var(--ink)}
  .defs-src a{color:var(--accent)}
  .err{color:var(--E);font-family:var(--mono);font-size:13px;margin-left:4px}

  /* warnings */
  .warn-box{font-family:var(--mono);font-size:12px;color:#8a6000;
    background:#fdf3d8;border:1px solid #e6c96a;border-radius:10px;
    padding:10px 14px;margin-top:12px;display:flex;flex-direction:column;gap:4px}
  .warn-box b{color:#8a6000}

  /* section separators */
  #out{margin-top:40px}
  .result-section{margin-bottom:48px}

  /* tab bar */
  .tab-bar{
    display:flex;align-items:stretch;gap:0;margin-bottom:24px;
    border-bottom:2px solid var(--line);position:relative;
  }
  .tab-btn{
    font-family:var(--mono);font-size:12px;font-weight:600;letter-spacing:.14em;
    text-transform:uppercase;cursor:pointer;
    padding:10px 20px 12px;border:none;background:none;
    color:var(--ink-faint);border-bottom:3px solid transparent;
    margin-bottom:-2px;transition:color .18s, border-color .18s, background .18s;
    border-radius:10px 10px 0 0;
  }
  .tab-btn:hover{color:var(--ink-dim);background:var(--panel-2)}
  .tab-btn.active{color:var(--accent);border-bottom-color:var(--accent);background:var(--panel)}
  .tab-btn.tab-seq.active{color:var(--accent-seq);border-bottom-color:var(--accent-seq)}
  .tab-btn .tab-count{
    display:inline-block;margin-left:7px;font-size:10px;
    background:var(--line);color:var(--ink-dim);
    border-radius:999px;padding:1px 7px;transition:background .18s,color .18s;
  }
  .tab-btn.active .tab-count{background:var(--accent);color:#fff}
  .tab-btn.tab-seq.active .tab-count{background:var(--accent-seq)}

  /* tab panels */
  .tab-panel{display:none}
  .tab-panel.active{display:block}
  /* sub-tab panels */
  .sub-tab-panel{display:none}
  .sub-tab-panel.active{display:block}
  .sub-tab-bar{margin-top:12px;border-top:none !important}
  .sub-tab-bar .tab-btn{font-size:12px;padding:6px 14px;text-transform:none;letter-spacing:.08em}

  /* keep section-banner for sequence section */
  .section-banner{
    display:flex;align-items:center;gap:12px;margin-bottom:20px;
    font-family:var(--mono);font-size:11px;letter-spacing:.18em;text-transform:uppercase;
  }
  .section-banner .pill{
    padding:4px 12px;border-radius:999px;font-weight:600;font-size:11px;
  }
  .section-banner.struct .pill{background:#0e7c6b18;color:var(--accent);border:1px solid #0e7c6b33}
  .section-banner.seq   .pill{background:#7c3aed18;color:var(--accent-seq);border:1px solid #7c3aed33}
  .section-banner::after{content:"";flex:1;height:1px;background:var(--line)}

  /* stats row */
  .summary{display:flex;flex-wrap:wrap;gap:10px;margin-bottom:26px}
  .stat{border:1px solid var(--line);border-radius:12px;padding:13px 16px;
    background:var(--panel);min-width:96px;position:relative;
    box-shadow:0 1px 4px #0000000a}
  .stat .k{font-family:var(--mono);font-size:12px;letter-spacing:.06em;
    color:var(--ink-dim);opacity:1}
  .stat .v{font-family:var(--serif);font-size:30px;font-weight:600;
    line-height:1.05;margin-top:3px;color:var(--ink)}
  .stat::before{content:"";position:absolute;left:0;top:12px;bottom:12px;
    width:3px;border-radius:3px;background:var(--c,var(--accent))}

  /* sub-section header */
  .section-h{font-family:var(--mono);font-size:12px;letter-spacing:.16em;text-transform:uppercase;
    color:var(--ink-faint);margin:30px 0 12px;display:flex;gap:10px;align-items:center}
  .section-h::after{content:"";flex:1;height:1px;background:var(--line)}

  /* motif rows */
  ol.list{list-style:none;margin:0;padding:0;display:flex;flex-direction:column;gap:8px}
  li.row{
    display:grid;
    grid-template-columns:36px 110px 1fr 230px 190px 180px;
    gap:20px;align-items:center;min-width:860px;
    background:var(--panel);border:1px solid var(--line);border-left:4px solid var(--c,var(--accent));
    border-radius:12px;padding:14px 22px;animation:rise .35s both;
    cursor:pointer;transition:all .15s;box-shadow:0 1px 4px #0000000a;
  }
  li.row:hover{background:var(--panel-2);box-shadow:0 2px 8px #00000014}
  li.row.active{background:#e8f5f2;box-shadow:0 0 0 2px var(--c,var(--accent))33}
  li.row.seq-row{cursor:default}
  @keyframes rise{from{opacity:0;transform:translateY(6px)}to{opacity:1;transform:none}}
  li.row .idx{font-family:var(--mono);font-size:12px;color:var(--ink-faint);text-align:right}
  .type{font-family:var(--mono);font-weight:600;font-size:14px;color:var(--c,var(--accent))}
  .seq-str{font-family:var(--mono);font-size:22px;letter-spacing:.36em;color:var(--ink)}
  .seq-str .end{color:var(--c,var(--accent));font-weight:600}
  .seq-str .mid{color:var(--ink-dim)}
  .meta{font-family:var(--mono);font-size:13px;color:var(--ink-dim);display:flex;flex-direction:column;gap:1px;min-width:0;overflow:hidden}
  .meta b{color:var(--ink);font-size:13px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  .meta .meta-sep{display:none}
  .meta .meta-sub{font-size:11px;color:var(--ink-faint);white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  /* secondary structure column */
  .ss-col{font-family:var(--mono);font-size:13px;letter-spacing:.18em;white-space:nowrap;
    display:flex;gap:1px;align-items:center}
  .ss-col span{display:inline-block;width:16px;text-align:center;border-radius:3px;padding:1px 0}
  .ss-H{color:#c0392b;background:#c0392b12;font-weight:600}
  .ss-E{color:#0e7c6b;background:#0e7c6b12;font-weight:600}
  .ss-B{color:#1a6b8a;background:#1a6b8a12}
  .ss-T{color:#6b3a9a;background:#6b3a9a12}
  .ss-S{color:#2d7a3a;background:#2d7a3a12}
  .ss-G{color:#b8860b;background:#b8860b12}
  .ss-I{color:#9a1a6b;background:#9a1a6b12}
  .ss--{color:var(--ink-faint);background:transparent}
  .ss-na{color:var(--ink-faint);font-size:11px;letter-spacing:.04em}

  /* interface indicator */
  .iface-col{font-family:var(--mono);font-size:12px;white-space:nowrap;
    display:flex;gap:3px;align-items:center;flex-wrap:wrap}
  .iface-num{display:inline-block;background:#e05c00;color:#fff;font-weight:700;
    font-size:10px;border-radius:4px;padding:1px 4px;
    box-shadow:0 1px 3px #e05c0033;letter-spacing:0}
  .iface-na{color:var(--ink-faint);font-size:11px;letter-spacing:.04em;font-family:var(--mono)}
  .iface-legend{font-family:var(--mono);font-size:10px;color:var(--ink-faint);
    display:flex;align-items:center;gap:4px;margin-bottom:8px}
  .iface-legend .swatch{display:inline-block;width:10px;height:10px;border-radius:2px;
    background:#e05c00;vertical-align:middle;margin-right:2px}
  .empty{font-family:var(--mono);color:var(--ink-faint);padding:24px;
    border:1px dashed var(--line);border-radius:12px;text-align:center;background:var(--panel)}

  /* structure results: list left, viewer right */
  .results-grid{display:grid;grid-template-columns:minmax(0,1fr) 370px;gap:20px;align-items:start}
  .results-grid>div:first-child{min-width:0;overflow:hidden}
  @media(max-width:900px){.results-grid{grid-template-columns:minmax(0,1fr)}}

  /* 3D viewer */
  .viewer-panel{position:sticky;top:24px;border:1px solid var(--line);
    border-radius:16px;background:var(--panel);overflow:hidden;
    box-shadow:0 2px 12px #00000012}
  .viewer-header{padding:12px 16px;border-bottom:1px solid var(--line);
    display:flex;align-items:center;gap:10px;flex-wrap:wrap;background:var(--panel-2)}
  .viewer-title{font-family:var(--mono);font-size:12px;letter-spacing:.12em;
    text-transform:uppercase;color:var(--ink-faint);flex:1}
  .viewer-select{font-family:var(--mono);font-size:12px;color:var(--ink);
    background:var(--panel);border:1px solid var(--line);border-radius:7px;
    padding:5px 10px;cursor:pointer;outline:none}
  .viewer-select:focus{border-color:var(--accent)}
  #mol-container{width:100%;height:380px;position:relative;background:#f0efe9}
  .viewer-legend{padding:10px 14px;border-top:1px solid var(--line);
    display:flex;flex-wrap:wrap;gap:8px;background:var(--panel-2)}
  .legend-pill{font-family:var(--mono);font-size:11px;padding:3px 9px;
    border-radius:999px;border:1px solid;display:flex;align-items:center;gap:5px}
  .legend-dot{width:8px;height:8px;border-radius:50%;flex-shrink:0}
  .viewer-hint{font-family:var(--mono);font-size:11px;color:var(--ink-faint);
    padding:6px 14px 10px;text-align:center;background:var(--panel-2)}

  .motif-list-wrap{overflow-x:auto;-webkit-overflow-scrolling:touch}
  /* column header row */
  .col-headers{
    display:grid;
    grid-template-columns:36px 110px 1fr 230px 190px 180px;
    gap:20px;align-items:center;min-width:860px;
    padding:7px 22px 7px 26px;
    border-bottom:1.5px solid var(--line);
    margin-bottom:6px;
  }
  .col-hdr{font-family:var(--mono);font-size:10.5px;font-weight:700;
    letter-spacing:.13em;text-transform:uppercase;color:var(--ink-dim);
    white-space:nowrap;user-select:none}
  .col-hdr:first-child{text-align:right}

  @media(max-width:700px){
    li.row,.col-headers{grid-template-columns:24px 1fr;row-gap:6px}
    li.row .type,.meta,.ss-col,.iface-col{grid-column:2}
    .col-hdr:not(:nth-child(-n+2)){display:none}
    li.row .idx{align-self:start}
  }

  /* ── omega loop rows ── */
  /* scrollable wrapper — prevents sequence column from collapsing to 0 */
  .omega-list-wrap{overflow-x:auto;-webkit-overflow-scrolling:touch}
  li.omega-row{
    display:grid;
    grid-template-columns:28px 105px minmax(100px,1fr) 155px 105px 232px;
    gap:12px;align-items:center;min-width:800px;
    background:var(--panel);border:1px solid var(--line);border-left:4px solid #e05c00;
    border-radius:12px;padding:12px 16px;animation:rise .35s both;
    cursor:pointer;transition:all .15s;box-shadow:0 1px 4px #0000000a;
  }
  li.omega-row:hover{background:var(--panel-2);box-shadow:0 2px 8px #00000014}
  li.omega-row.active{background:#fff3ec;box-shadow:0 0 0 2px #e05c0055}
  [data-theme="dark"] li.omega-row.active{background:#2a1800}
  .omega-col-headers{
    display:grid;
    grid-template-columns:28px 105px minmax(100px,1fr) 155px 105px 232px;
    gap:12px;align-items:center;min-width:800px;
    padding:7px 16px 7px 20px;
    border-bottom:1.5px solid var(--line);
    margin-bottom:6px;
  }
  /* align # header with the centered .idx numbers in rows */
  .omega-col-headers .col-hdr:first-child{text-align:center}
  .omega-id{
    display:flex;flex-direction:column;gap:2px;min-width:0;overflow:hidden;
  }
  .omega-id .oid-pdb{
    font-family:var(--mono);font-size:12px;font-weight:700;color:var(--ink);
    white-space:nowrap;overflow:hidden;text-overflow:ellipsis;
  }
  .omega-id .oid-loc{
    font-family:var(--mono);font-size:10px;color:var(--ink-faint);white-space:nowrap;
  }
  /* sequence cell — clamp but never invisible */
  li.omega-row .seq-str{min-width:0;overflow:hidden;display:block;white-space:nowrap}
  /* Radius · Area · Volume stacked */
  .geo-vals{
    font-family:var(--mono);font-size:11px;color:var(--ink-dim);
    display:flex;flex-direction:column;gap:2px;align-items:flex-start;
  }
  .geo-vals b{color:var(--ink);font-weight:600}
  .geo-lbl{display:inline-block;width:52px;color:var(--ink-faint);font-size:10px}
  /* SS loop codes — plain flex, no flanking chips */
  .omega-ss{font-family:var(--mono);font-size:11px;
    display:flex;gap:2px;align-items:center;}
  .omega-dist{font-family:var(--mono);font-size:12px;color:var(--ink-dim);white-space:nowrap}
  .omega-dist b{color:#e05c00;font-weight:700}
  @media(max-width:700px){
    li.omega-row,.omega-col-headers{grid-template-columns:24px 1fr;row-gap:6px;min-width:unset}
    li.omega-row .geo-vals,.omega-dist,.omega-ss,.omega-id{grid-column:2}
  }

  /* ── β/γ-turn rows ── */
  .turn-list-wrap{overflow-x:auto;-webkit-overflow-scrolling:touch}
  li.turn-row{
    display:grid;
    grid-template-columns:28px 105px minmax(80px,1fr) 115px 115px 115px 155px;
    gap:10px;align-items:center;min-width:820px;
    background:var(--panel);border:1px solid var(--line);
    border-left:4px solid var(--turn-c,#7c3aed);
    border-radius:12px;padding:11px 14px;animation:rise .35s both;
    cursor:pointer;transition:all .15s;box-shadow:0 1px 4px #0000000a;
  }
  li.turn-row:hover{background:var(--panel-2);box-shadow:0 2px 8px #00000014}
  li.turn-row.active{box-shadow:0 0 0 2px var(--turn-c,#7c3aed)55}
  li.turn-row.beta-turn{--turn-c:#7c3aed}
  li.turn-row.beta-turn.active{background:#f3f0ff}
  li.turn-row.gamma-turn{--turn-c:#1a6b8a}
  li.turn-row.gamma-turn.active{background:#e8f3f8}
  li.turn-row.alpha-turn{--turn-c:#b71c1c}
  li.turn-row.alpha-turn.active{background:#ffebee}
  [data-theme="dark"] li.turn-row.beta-turn.active{background:#1a0f2e}
  [data-theme="dark"] li.turn-row.gamma-turn.active{background:#001f2e}
  [data-theme="dark"] li.turn-row.alpha-turn.active{background:#2d0606}
  .turn-col-headers{
    display:grid;
    grid-template-columns:28px 105px minmax(80px,1fr) 115px 115px 115px 155px;
    gap:10px;align-items:center;min-width:820px;
    padding:7px 14px 7px 20px;
    border-bottom:1.5px solid var(--line);margin-bottom:6px;
  }
  .turn-subtype{
    font-family:var(--mono);font-size:11px;font-weight:700;
    padding:2px 7px;border-radius:5px;letter-spacing:.06em;
    display:inline-block;white-space:nowrap;
  }
  .turn-subtype.beta{background:#7c3aed18;color:#7c3aed;border:1px solid #7c3aed33}
  .turn-subtype.gamma-cls{background:#1a6b8a18;color:#1a6b8a;border:1px solid #1a6b8a33}
  .turn-subtype.gamma-inv{background:#0e7c6b18;color:#0e7c6b;border:1px solid #0e7c6b33}
  .turn-subtype.gamma-cls-nhn{background:#6d3a9a18;color:#6d3a9a;border:1px solid #6d3a9a33;
    white-space:normal;line-height:1.3;font-size:10px;text-align:center;max-width:113px}
  .turn-subtype.gamma-inv-nhn{background:#8a4a1a18;color:#8a4a1a;border:1px solid #8a4a1a33;
    white-space:normal;line-height:1.3;font-size:10px;text-align:center;max-width:113px}
  .turn-subtype.unc{background:var(--panel-2);color:var(--ink-faint);border:1px solid var(--line)}
  .turn-angles{font-family:var(--mono);font-size:11px;color:var(--ink-dim);
    display:flex;flex-direction:column;gap:2px}
  .turn-angles .angle-row{display:flex;gap:6px;align-items:center}
  .turn-angles .alabel{color:var(--ink-faint);font-size:10px;width:18px;flex-shrink:0}
  .turn-angles b{color:var(--ink)}
  .turn-dist{font-family:var(--mono);font-size:12px;color:var(--ink-dim);white-space:nowrap}
  .turn-dist b{color:var(--turn-c,#7c3aed);font-weight:700}
  @media(max-width:700px){
    li.turn-row,.turn-col-headers{grid-template-columns:24px 1fr;row-gap:6px;min-width:unset}
    li.turn-row .turn-angles,.turn-dist,.turn-subtype,.omega-id{grid-column:2}
  }
  /* ── fullscreen 3D explorer ──────────────────────────────────────────── */
  #fs-overlay{display:none;position:fixed;inset:0;z-index:9999;background:#0d1117;flex-direction:row}
  #fs-mol-container{flex:1;height:100%;min-width:0}
  #fs-panel{width:340px;height:100%;background:#161b22;display:flex;flex-direction:column;
    overflow:hidden;border-left:1px solid #30363d;flex-shrink:0}
  .fs-panel-hdr{display:flex;justify-content:space-between;align-items:center;
    padding:14px 18px;border-bottom:1px solid #30363d;background:#0d1117;flex-shrink:0}
  .fs-panel-hdr-left{display:flex;flex-direction:column;gap:2px}
  .fs-explorer-lbl{font-family:var(--mono);font-size:9px;color:#8b949e;letter-spacing:.14em;text-transform:uppercase}
  #fs-pdb-title{font-family:var(--mono);font-size:15px;font-weight:700;color:#e6edf3}
  .fs-hdr-actions{display:flex;gap:8px;align-items:center}
  .fs-reset-btn{background:none;border:1px solid #30363d;border-radius:6px;color:#8b949e;
    cursor:pointer;padding:5px 10px;font-size:11px;font-family:var(--mono);transition:all .15s}
  .fs-reset-btn:hover{background:#21262d;color:#e6edf3}
  .fs-close-btn{background:none;border:1px solid #30363d;border-radius:6px;color:#8b949e;
    cursor:pointer;padding:5px 12px;font-size:13px;transition:all .15s}
  .fs-close-btn:hover{background:#da3633;border-color:#da3633;color:#fff}
  .fs-open-btn{background:none;border:1px solid var(--line);border-radius:6px;color:var(--ink-dim);
    cursor:pointer;padding:3px 8px;font-size:14px;line-height:1;transition:all .15s;margin-left:6px}
  .fs-open-btn:hover{background:var(--panel-2);color:var(--ink);border-color:var(--accent)}
  .fs-sections{flex:1;overflow-y:auto;scrollbar-width:thin;scrollbar-color:#30363d transparent}
  .fs-section{padding:14px 18px;border-bottom:1px solid #21262d}
  .fs-sec-title{font-family:var(--mono);font-size:9px;color:#8b949e;
    letter-spacing:.12em;text-transform:uppercase;margin-bottom:10px}
  .fs-repr-grid{display:grid;grid-template-columns:1fr 1fr;gap:6px}
  .fs-repr-btn{background:#21262d;border:1px solid #30363d;border-radius:6px;color:#8b949e;
    cursor:pointer;padding:7px;font-size:12px;font-family:var(--mono);transition:all .15s;text-align:center}
  .fs-repr-btn:hover{background:#30363d;color:#e6edf3}
  .fs-repr-btn.active{background:#1f4168;border-color:#388bfd;color:#79c0ff}
  .fs-chain-row{display:flex;align-items:center;gap:10px;padding:7px 0;border-bottom:1px solid #1c2128}
  .fs-chain-row:last-child{border-bottom:none}
  .fs-color-pick{width:28px;height:28px;border:2px solid #30363d;border-radius:6px;
    cursor:pointer;padding:0;flex-shrink:0;background:none}
  .fs-color-pick::-webkit-color-swatch-wrapper{padding:2px}
  .fs-color-pick::-webkit-color-swatch{border:none;border-radius:3px}
  .fs-chain-badge{font-family:var(--mono);font-size:12px;font-weight:700;padding:3px 10px;
    border-radius:5px;color:#fff;text-shadow:0 1px 2px rgba(0,0,0,.5);flex-shrink:0}
  .fs-chain-count{font-family:var(--mono);font-size:11px;color:#8b949e;margin-left:auto}
  .fs-hint{font-family:var(--mono);font-size:11px;color:#8b949e;font-style:italic}
  .fs-info-card{display:flex;align-items:center;gap:14px;background:#21262d;
    border:1px solid #30363d;border-radius:8px;padding:12px}
  .fs-info-aa{font-family:var(--mono);font-size:30px;font-weight:700;color:#79c0ff;
    width:38px;text-align:center;flex-shrink:0}
  .fs-info-details{font-family:var(--mono);font-size:12px;color:#c9d1d9;
    display:flex;flex-direction:column;gap:3px}
  .fs-info-details b{color:#e6edf3}
  .fs-info-atom{font-size:10px;color:#8b949e;margin-top:2px}
  .fs-chain-seq{margin-bottom:10px}
  .fs-seq-label{font-family:var(--mono);font-size:10px;color:#8b949e;margin-bottom:4px}
  .fs-seq-strip{display:flex;gap:1px;overflow-x:auto;padding-bottom:4px;
    scrollbar-width:thin;scrollbar-color:#30363d transparent}
  .fs-res{font-family:var(--mono);font-size:11px;display:inline-flex;align-items:center;
    justify-content:center;width:17px;height:19px;border-radius:3px;flex-shrink:0;
    cursor:pointer;color:#8b949e;background:#21262d;transition:background .1s}
  .fs-res:hover{background:#388bfd;color:#fff}
  .fs-res.active{background:#1f6feb;color:#fff;outline:2px solid #79c0ff;outline-offset:1px}
  /* motif list replaced by fs-annot-list — see annotation tabs CSS below */
  /* ── fs annotation tabs (X***X / Ω / β / γ / α) ── */
  .fs-annot-section{padding-bottom:6px !important}
  .fs-annot-tabs{display:flex;gap:3px;margin-bottom:10px}
  .fs-annot-tab{
    flex:1;font-family:var(--mono);font-size:11px;font-weight:700;
    background:#21262d;border:1px solid #30363d;border-radius:5px;
    color:#8b949e;cursor:pointer;padding:5px 2px;text-align:center;
    transition:all .15s;letter-spacing:.01em;
  }
  .fs-annot-tab:hover{background:#30363d;color:#e6edf3}
  .fs-annot-tab[data-tab="motif"].active{background:#1f4168;border-color:#388bfd;color:#79c0ff}
  .fs-annot-tab[data-tab="omega"].active{background:#2d1500;border-color:#e05c00;color:#ff8c42}
  .fs-annot-tab[data-tab="beta"].active{background:#1a0f2e;border-color:#7c3aed;color:#a78bfa}
  .fs-annot-tab[data-tab="gamma"].active{background:#001f2e;border-color:#1a6b8a;color:#4dd0e1}
  .fs-annot-tab[data-tab="alpha"].active{background:#2d0606;border-color:#b71c1c;color:#ff6f6f}
  .fs-annot-pane{display:none}
  .fs-annot-pane.active{display:block}
  .fs-annot-list{display:flex;flex-direction:column;gap:1px;max-height:220px;
    overflow-y:auto;scrollbar-width:thin;scrollbar-color:#30363d transparent}
  .fs-annot-item{display:grid;grid-template-columns:38px 1fr auto;align-items:center;
    gap:6px;padding:6px 4px;border-radius:5px;cursor:pointer;transition:background .12s}
  .fs-annot-item:hover{background:#21262d}
  .fs-annot-item.active{outline:1px solid #388bfd;background:#1f4168}
  .fs-annot-item[data-atype="omega"].active{background:#2d1500;outline-color:#e05c00}
  .fs-annot-item[data-atype="beta"].active{background:#1a0f2e;outline-color:#7c3aed}
  .fs-annot-item[data-atype="gamma"].active{background:#001f2e;outline-color:#1a6b8a}
  .fs-annot-item[data-atype="alpha"].active{background:#2d0606;outline-color:#b71c1c}
  .fs-annot-badge{font-family:var(--mono);font-size:10px;font-weight:700;
    text-align:center;padding:2px 4px;border-radius:4px;background:#ffffff12;
    letter-spacing:.03em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}
  .fs-annot-seq{font-family:var(--mono);font-size:12px;letter-spacing:.15em;color:#e6edf3;
    min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
  .fs-annot-loc{font-family:var(--mono);font-size:10px;color:#8b949e;
    white-space:nowrap;text-align:right}
</style>
</head>
<body>
<!-- ── fullscreen 3D explorer overlay ──────────────────────────── -->
<div id="fs-overlay">
  <div id="fs-mol-container"></div>
  <div id="fs-panel">
    <div class="fs-panel-hdr">
      <div class="fs-panel-hdr-left">
        <span class="fs-explorer-lbl">3D Explorer</span>
        <span id="fs-pdb-title">—</span>
      </div>
      <div class="fs-hdr-actions">
        <button class="fs-reset-btn" onclick="fsResetView()">↺ Reset</button>
        <button class="fs-close-btn" onclick="closeFsViewer()" title="Close (Esc)">✕</button>
      </div>
    </div>
    <div class="fs-sections">
      <div class="fs-section">
        <div class="fs-sec-title">Representation</div>
        <div class="fs-repr-grid">
          <button class="fs-repr-btn active" onclick="fsSetStyle('cartoon',this)">Cartoon</button>
          <button class="fs-repr-btn" onclick="fsSetStyle('stick',this)">Stick</button>
          <button class="fs-repr-btn" onclick="fsSetStyle('sphere',this)">Sphere</button>
          <button class="fs-repr-btn" onclick="fsSetStyle('line',this)">Line</button>
        </div>
      </div>
      <div class="fs-section">
        <div class="fs-sec-title">Chains &amp; Colours</div>
        <div id="fs-chains"><span class="fs-hint">Loading…</span></div>
      </div>
      <div class="fs-section fs-annot-section">
        <div class="fs-sec-title">Annotations</div>
        <div class="fs-annot-tabs">
          <button class="fs-annot-tab active" data-tab="motif" onclick="switchFsAnnotTab('motif')" title="X***X Sequence Motifs">X***X</button>
          <button class="fs-annot-tab" data-tab="omega" onclick="switchFsAnnotTab('omega')" title="Ω-loops">Ω</button>
          <button class="fs-annot-tab" data-tab="beta"  onclick="switchFsAnnotTab('beta')"  title="β-Turns">β</button>
          <button class="fs-annot-tab" data-tab="gamma" onclick="switchFsAnnotTab('gamma')" title="γ-Turns">γ</button>
          <button class="fs-annot-tab" data-tab="alpha" onclick="switchFsAnnotTab('alpha')" title="α-Turns">α</button>
        </div>
        <div class="fs-annot-pane active" data-pane="motif" id="fs-annot-pane-motif">
          <span class="fs-hint">No motifs loaded</span>
        </div>
        <div class="fs-annot-pane" data-pane="omega" id="fs-annot-pane-omega">
          <span class="fs-hint">No loops loaded</span>
        </div>
        <div class="fs-annot-pane" data-pane="beta" id="fs-annot-pane-beta">
          <span class="fs-hint">No β-turns loaded</span>
        </div>
        <div class="fs-annot-pane" data-pane="gamma" id="fs-annot-pane-gamma">
          <span class="fs-hint">No γ-turns loaded</span>
        </div>
        <div class="fs-annot-pane" data-pane="alpha" id="fs-annot-pane-alpha">
          <span class="fs-hint">No α-turns loaded</span>
        </div>
      </div>
      <div class="fs-section">
        <div class="fs-sec-title">Selected Residue</div>
        <div id="fs-info"><span class="fs-hint">Click any atom to view details</span></div>
      </div>
      <div class="fs-section">
        <div class="fs-sec-title">Sequence</div>
        <div id="fs-seq-area"><span class="fs-hint">Loading…</span></div>
      </div>
    </div>
  </div>
</div>
<!-- ── tour guide card ──────────────────────────────────────────── -->
<div id="tour-card">
  <div class="tour-hdr">
    <div class="tour-hdr-left">
      <span class="tour-eyebrow">guided tour</span>
      <span class="tour-title" id="tour-title"></span>
    </div>
    <span class="tour-step-ind" id="tour-step-ind"></span>
  </div>
  <div class="tour-body" id="tour-body"></div>
  <div class="tour-footer">
    <button class="tour-close-lnk" onclick="endTour()">close</button>
    <div class="tour-dots" id="tour-dots"></div>
    <div class="tour-nav">
      <button class="tour-btn" id="tour-prev" onclick="moveTour(-1)">← prev</button>
      <button class="tour-btn primary" id="tour-next" onclick="moveTour(1)">next →</button>
    </div>
  </div>
</div>

<!-- ── definitions overlay ──────────────────────────────────────── -->
<div id="defs-overlay" onclick="if(event.target===this)closeDefs()">
  <div class="defs-panel">
    <div class="defs-hdr">
      <div class="defs-hdr-left">
        <span class="defs-eyebrow">motif scanner · reference</span>
        <span class="defs-title">Detection <em>Definitions</em></span>
      </div>
      <button class="defs-close-btn" onclick="closeDefs()">✕</button>
    </div>
    <div class="defs-body">

      <!-- X***X -->
      <div class="defs-section">
        <div class="defs-sec-hdr">
          <span class="defs-tag" style="background:#0e7c6b18;color:var(--accent);border:1px solid #0e7c6b33">Sequence Motifs</span>
          <span class="defs-sec-title">X***X Motifs</span>
        </div>
        <p class="defs-prose">A 5-residue motif in which the first and fifth residues are A, E, I, L, V, K, or G (e.g., AXXXA, GXXXG, LXXXL, IXXXI, KXXXK, VXXXV, EXXXE), with any three amino acids between them. These motifs can be detected from input sequence or structure files.</p>
        <div class="defs-criteria">
          <div class="defs-criteria-lbl">Detection criteria</div>
          <div class="defs-criterion"><span class="defs-cnum">1</span><span>Window: <b>5 residues</b> [i, i+1, i+2, i+3, i+4]</span></div>
          <div class="defs-criterion"><span class="defs-cnum">2</span><span>The residues at positions i and i + 4 must be identical <b> amino acids.</b></span></div>
          <div class="defs-criterion"><span class="defs-cnum">3</span><span>All five residues must belong to the <b>same chain.</b> </span></div>
          <div class="defs-criterion"><span class="defs-cnum">4</span><span>For structure files: secondary structure assignments are calculated using DSSP, and interface residues associated with each motif are identified using FreeSASA.</span></div>
        </div>
        <div class="defs-sources">
          <div class="defs-src-lbl">Sources</div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>DSSP</b> — Kabsch W, Sander C &amp; (1983). Dictionary of protein secondary structure. <i>PROTEIN SCIENCE.</i> <a href="  https://doi.org/10.1002/bip.360221211" target="_blank">doi:10.1093/nar/gkg626</a></span></div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>FreeSASA</b> — Mitternacht S (2016). FreeSASA. <i>F1000Research</i> 5:189. <a href="https://doi.org/10.12688/f1000research.7931.1" target="_blank">doi:10.12688/f1000research.7931.1</a> </span></div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>Motif research article</b> —Tajane et.al (2024).<a href="https://doi.org/10.1016/j.jsb.2024.108129" target="_blank">doi:10.1016/j.jsb.2024.108129</a></span></div>
        </div>
      </div>

      <!-- γ-turns -->
      <div class="defs-section">
        <div class="defs-sec-hdr">
          <span class="defs-tag" style="background:#e05c0018;color:#e05c00;border:1px solid #e05c0033">Structural Motifs</span>
          <span class="defs-tag" style="background:#1a6b8a18;color:#56c8e0;border:1px solid #1a6b8a33"> 1 </span>
          <span class="defs-sec-title">γ-Turns</span>
        </div>
        <p class="defs-prose"> A three-residue structural motif stabilized by a C=O(i)···H-N(i+2) Hydrogen-bond.</p>
        <div class="defs-criteria">
          <div class="defs-criteria-lbl">Detection criteria (all required)</div>
          <div class="defs-criterion"><span class="defs-cnum">1</span><span>Window: <b>3 residues</b> [i, i+1, i+2]</span></div>
          <div class="defs-criterion"><span class="defs-cnum">2</span><span><b>Hydrogen bond between C=O(i) and N–H(i+2)</b>: O(i)···H(i+2) ≤ 2.5 Å ; N(i+2)···O(i) ≤ 3.5 Å , N–H···O ≥ 120°, and H···O=C ≥ 90°.</span></div>
          <div class="defs-criterion"><span class="defs-cnum">3</span><span>Types of γ-turns depending on φ/ψ value of i+1 residue, matched within <b>±40°</b> of the classic/inverse centroid (per angle):<br>
            &nbsp;&nbsp;<b>classic γ-turn</b> centroid: φ = 75°, ψ = −64°<br>
            &nbsp;&nbsp;<b>inverse γ-turn</b> centroid: φ = −79°, ψ = 69°</span></div>
        </div>
        <table class="defs-tbl">
          <thead><tr><th>Subtype</th><th>φ(i+1) centroid</th><th>ψ(i+1) centroid</th><th>Tolerance</th></tr></thead>
          <tbody>
            <tr><td class="tn">classic γ-turn</td><td>75°</td><td>−64°</td><td>±40° per angle</td></tr>
            <tr><td class="tn">inverse γ-turn</td><td>−79°</td><td>69°</td><td>±40° per angle</td></tr>
          </tbody>
        </table>
        <div class="defs-sources">
          <div class="defs-src-lbl">Sources</div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>γ-turn definition, classic vs inverse classification, original φ/ψ centroid values</b> — Rose GD, Gierasch LM &amp; Smith JA (1985). <i>Adv. Protein Chem.</i> 37:1–109. <a href="https://doi.org/10.1016/S0065-3233(08)60063-7" target="_blank">doi:10.1016/S0065-3233(08)60063-7</a></span></div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>H-bond geometry criteria (O···H &lt; 2.5 Å, angles &gt; 90°/120°)</b> — Milner-White EJ, Ross BM, Ismail R, Belhadj-Mostefa K &amp; Poet R (1988). <i>J. Mol. Biol.</i> 204:777–782. <a href="https://doi.org/10.1016/0022-2836(88)90368-3" target="_blank">doi:10.1016/0022-2836(88)90368-3</a></span></div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>±40° per-angle tolerance around the classic/inverse centroid</b> — Hutchinson EG &amp; Thornton JM (1996). PROMOTIF. <i>Protein Sci.</i> 5:212–220. <a href="https://doi.org/10.1002/pro.5560050204" target="_blank">doi:10.1002/pro.5560050204</a></span></div>
        </div>

        <!-- NHN γ-turns subsection -->
        <div class="defs-sec-hdr" style="margin-top:22px">
          <span class="defs-tag" style="background:#6d3a9a18;color:#9c6fd4;border:1px solid #6d3a9a33">NHN γ-Turns</span>
        </div>
        <p class="defs-prose">A subclass of γ-turns in which stabilizing by <b>N–H···N</b> hydrogen bond rather than the classical C=O···H–N bond. The central residue i+1 acts as the nitrogen acceptor, receiving the N–H(i+2) donor hydrogen. Identified by Dhar &amp; et.al (2019).</p>
        <div class="defs-criteria">
          <div class="defs-criteria-lbl">Detection criteria (all required)</div>
          <div class="defs-criterion"><span class="defs-cnum">1</span><span>Window: <b>3 residues</b> [i, i+1, i+2] </span></div>
          <div class="defs-criterion"><span class="defs-cnum">2</span><span>Cα(i)···Cα(i+2) distance: <b>5.2–5.8 Å</b></span></div>
          <div class="defs-criterion"><span class="defs-cnum">3</span><span>θ value for N–H<sub>i+2</sub>···N(p<sub>z</sub>)<sub>i+1</sub> interaction: <b>θ &lt; 30°</b></span></div>
          <div class="defs-criterion"><span class="defs-cnum">4</span><span>C5 intraresidue ring at i: <b>N-H(i)···O=C(i) &lt; 3.0 Å</b></span></div>
          <div class="defs-criterion"><span class="defs-cnum">5</span><span>φ/ψ of central residue i+1 must match a classical or inverse NHN cluster</span></div>
        </div>
        <p class="defs-prose" style="margin-top:6px"><b>Subtypes</b> are assigned by the φ/ψ angles of the central residue i:</p>
        <table class="defs-tbl">
          <thead><tr><th>Subtype</th><th>φ(i+1)</th><th>ψ(i+1)</th></tr></thead>
          <tbody>
            <tr>
              <td class="tn" style="color:#6d3a9a">classical NHN γ-turn</td>
              <td>+59±13°</td><td>+36±16°</td>
            </tr>
            <tr>
              <td class="tn" style="color:#8a4a1a">inverse NHN γ-turn</td>
              <td>−72±18°</td><td>−22±20°</td>
            </tr>
          </tbody>
        </table>
        <div class="defs-sources">
          <div class="defs-src-lbl">Sources</div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>NHN γ-turn definition, detection criteria, and φ/ψ cluster ranges</b> — Dhar J, Kishore R &amp; Chakrabarti P (2019). <i>Proteins</i> 87:965–978. <a href="https://doi.org/10.1002/prot.25820" target="_blank">doi:10.1002/prot.25820</a></span></div>
        </div>
      </div>

      <!-- β-turns -->
      <div class="defs-section">
        <div class="defs-sec-hdr">
          <span class="defs-tag" style="background:#e05c0018;color:#e05c00;border:1px solid #e05c0033"> 2 </span>
          <span class="defs-sec-title">β-Turns</span>
        </div>
        <p class="defs-prose">A four-residue motif that helps reverse the direction of the polypeptide chain. β-turns are the most abundant type of motif.</p>
        <div class="defs-criteria">
          <div class="defs-criteria-lbl">Detection criteria</div>
          <div class="defs-criterion"><span class="defs-cnum">1</span><span>Window: <b>4 residues</b> [i, i+1, i+2, i+3]</span></div>
          <div class="defs-criterion"><span class="defs-cnum">2</span><span>Cα(i)→Cα(i+3) distance <b>&lt; 7.0 Å</b> (End dist (Å))(Wilmot &amp; Thornton 1990)</span></div>
          <div class="defs-criterion"><span class="defs-cnum">3</span><span>Central residues i+1 and i+2  <b>not</b> be helical — DSSP code <b>∉ {G, H, I}</b> (Wilmot &amp; Thornton 1990)</span></div>
          <div class="defs-criterion"><span class="defs-cnum">4</span><span>A hydrogen bond may or may not be present. (Wilmot &amp; Thornton 1990)</span></div>
        </div>
        <p class="defs-prose" style="margin-top:4px"><b>Subtype classification:</b> φ/ψ of i+1 and i+2 matched with <b>±30° tolerance.</b>  Types III and III′ are <b>not included</b> — they overlap extensively with 3₁₀-helix. (Wilmot &amp; Thornton 1990).</p>
        <table class="defs-tbl">
          <thead><tr><th>SUBTYPE</th><th>φ(i+1)</th><th>ψ(i+1)</th><th>φ(i+2)</th><th>ψ(i+2)</th></tr></thead>
          <tbody>
            <tr><td class="tn">I</td><td>−60°</td><td>−30°</td><td>−90°</td><td>0°</td><td class="nd"></tr>
            <tr><td class="tn">I′</td><td>+60°</td><td>+30°</td><td>+90°</td><td>0°</td><td class="nd"></tr>
            <tr><td class="tn">II</td><td>−60°</td><td>+120°</td><td>+80°</td><td>0°</td><td class="nd"></tr>
            <tr><td class="tn">II′</td><td>+60°</td><td>−120°</td><td>−80°</td><td>0°</td><td class="nd"></tr>
            <tr><td class="tn">VIa</td><td>−60°</td><td>+120°</td><td>−90°</td><td>0°</td><td class="nd"></tr>
            <tr><td class="tn">VIb</td><td>−120°</td><td>+120°</td><td>−60°</td><td>0°</td><td class="nd"></tr>
            <tr><td class="tn">VIII</td><td>−60°</td><td>−30°</td><td>−120°</td><td>+120°</td><td class="nd"></tr>
            <tr><td class="tn">IV</td><td colspan="4">catch-all — satisfies criteria 1–4 above but no angle φ/ψ matches.</td><td class="nd"></tr>
          </tbody>
        </table>
        <div class="defs-sources">
          <div class="defs-src-lbl">Sources</div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>Distance criterion &lt; 7.0 Å; Terminal residue Hydrogen-bond may or may not present.</b> — Alexandre G.de Brevern (2022). <i>International journal of molecular science</i> <a href="https://doi.org/10.3390/ijms232012314" target="_blank">doi:10.3390/ijms232012314</a></span></div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>Idealized φ/ψ angles, ±30° tolerance</b> — C.M.Wilmot &amp; Thornton JM (1990). <i>Protein Engineering</i> 3:479–493. <a href="https://doi.org/10.1093/protein/3.6.479" target="_blank">doi:10.1093/protein/3.6.479</a></span></div>
        </div>
      </div>

      <!-- α-turns -->
      <div class="defs-section">
        <div class="defs-sec-hdr">
          <span class="defs-tag" style="background:#b71c1c18;color:#b71c1c;border:1px solid #b71c1c33"> 3 </span>
          <span class="defs-sec-title">α-Turns</span>
        </div>
        <p class="defs-prose"> A five-residue structural motifs stabilized by a N-H(i+4)→C=O(i) hydrogen-bond.</p>
        <div class="defs-criteria">
          <div class="defs-criteria-lbl">Detection criteria (all required)</div>
          <div class="defs-criterion"><span class="defs-cnum">1</span><span>Window: <b>5 residues</b> [i, i+1, i+2, i+3, i+4]</span></div>
          <div class="defs-criterion"><span class="defs-cnum">2</span><span>Cα(i) → Cα(i+4) distance <b>&lt; 6.5 Å</b> (End dist (Å)) (Dasgupta et al. 2004)</span></div>
          <div class="defs-criterion"><span class="defs-cnum">3</span><span>i+1, i+2, i+3 <b>not-helical(not assigned G or H)</b>(Dasgupta et al. 2004)</span></div>
          <div class="defs-criterion"><span class="defs-cnum">4</span><span>N(i+4)···O(i) distance: <b>2.6–3.6 Å</b> (Nataraj et al. 1995)</span></div>
          <div class="defs-criterion"><span class="defs-cnum">5</span><span>H–N···O angle: <b>≤ 40°</b> (Nataraj et al. 1995)</span></div>
        </div>
        <p class="defs-prose" style="margin-top:4px"><b></b> φ/ψ of middle residue i+1, i+2, i+3 from Pavone et al. (1996) Table IV. Each angle must be within <b>± SD</b></p>
        <table class="defs-tbl">
          <thead><tr><th>Subtype</th><th>φ(i+1)</th><th>ψ(i+1)</th><th>φ(i+2)</th><th>ψ(i+2)</th><th>φ(i+3)</th><th>ψ(i+3)</th></tr></thead>
          <tbody>
            <tr><td class="tn">I-αRS</td><td>−60±11°</td><td>−29±13°</td><td>−72±14°</td><td>−29±15°</td><td>−96±20°</td><td>−20±17°</td></tr>
            <tr><td class="tn">I-αLS</td><td>+48±22°</td><td>+42±14°</td><td>+67±9°</td><td>+33±14°</td><td>+70±11°</td><td>+32±12°</td></tr>
            <tr><td class="tn">II-αRS</td><td>−59±10°</td><td>+129±15°</td><td>+88±15°</td><td>−16±19°</td><td>−91±22°</td><td>−32±18°</td></tr>
            <tr><td class="tn">II-αLS</td><td>+53±15°</td><td>−137±25°</td><td>−95±12°</td><td>+81±23°</td><td>+57±5°</td><td>+38±8°</td></tr>
            <tr><td class="tn">I-αLU</td><td>−61±12°</td><td>+158±15°</td><td>+64±17°</td><td>+37±21°</td><td>+62±12°</td><td>+39±8°</td></tr>
            <tr><td class="tn">I-αRU</td><td>+59±18°</td><td>−157±31°</td><td>−67±17°</td><td>−29±20°</td><td>−68±12°</td><td>−39±12°</td></tr>
            <tr><td class="tn">II-αLU</td><td>−65±15°</td><td>−20±15°</td><td>−90±17°</td><td>+16±44°</td><td>+86±18°</td><td>+37±27°</td></tr>
            <tr><td class="tn">II-αRU</td><td>+54±8°</td><td>+39±15°</td><td>+67±13°</td><td>−5±31°</td><td>−125±11°</td><td>−34±32°</td></tr>
            <tr><td class="tn">I-αc</td><td>−103±23°</td><td>+143±4°</td><td>−85±8°</td><td>+2±6°</td><td>−54±6°</td><td>−39±9°</td></tr>
            <tr><td class="tn">unclassified</td><td colspan="6">satisfies criteria 1–5 above but matches no cluster template — reported per Pavone et al. (1996), who list such turns separately rather than excluding them</td></tr>
          </tbody>
        </table>
        <div class="defs-sources">
          <div class="defs-src-lbl">Sources</div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>α-turn identification, H-bond criteria (N···O 2.6–3.6 Å, H–N···O ≤ 40°)</b>— Nataraj DV, Srinivasan N, Sowdhamini R &amp; Ramakrishnan C (1995).<i>Current Science</i> 69:434–447.</span></div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>Nine α-turn types based on φ and ψ angles of middle residues i+1, i+2, and i+3.</b>— Pavone V, Gaeta G, Lombardi A, Nastri F, Maglio O, Isernia C &amp; Saviano M (1996). <i>Biopolymers</i> 38(6):705–721. <a href="https://doi.org/10.1002/(SICI)1097-0282(199606)38:6%3C705::AID-BIP3%3E3.0.CO;2-V" target="_blank">doi:10.1002/(SICI)1097-0282(199606)38:6&lt;705::AID-BIP3&gt;3.0.CO;2-V</a></span></div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>Cα(i)–Cα(i+4) &lt; 6.5 Å distance criteria:</b> (End dist (Å))  — Dasgupta B, Pal L, Basu G &amp; Chakrabarti P (2004). <i>Proteins</i> 55:305–315. <a href="https://doi.org/10.1002/prot.20064" target="_blank">doi:10.1002/prot.20064</a></span></div>
        </div>
      </div>

      <!-- Ω-loops -->
      <div class="defs-section">
        <div class="defs-sec-hdr">
          <span class="defs-tag" style="background:#e05c0018;color:#e05c00;border:1px solid #e05c0033"> 4 </span>
          <span class="defs-sec-title">Ω-Loops</span>
        </div>
        <p class="defs-prose">A six to sixteen-residue nonregular secondary structural element that resembles the Greek letter Ω (omega). It was identified by Leszczynski &amp; Rose (1986) as a distinct structural motif separate from regular secondary structure elements such as helices and strands.</p>
        <div class="defs-criteria">
          <div class="defs-criteria-lbl">Detection criteria (all required)</div>
          <div class="defs-criterion"><span class="defs-cnum">1</span><span>Window: <b>6-16 residues</b></span></div>
          <div class="defs-criterion"><span class="defs-cnum">2</span><span>At least <b>5 out of 6 residue</b> secondary structure must be <b>S</b> (bend) or <b>T</b> (turn) [Loop length -1 residue should be in S or T]</span></div>
          <div class="defs-criterion"><span class="defs-cnum">4</span><span>Cα(i)→Cα(i+5) <b>&lt; 10 Å.</b> (End dist (Å))</span></div>
          <div class="defs-criterion"><span class="defs-cnum">5</span><span>Endpoint distance <b>&lt; ⅔ × max pairwise Cα distance</b> in loop.</span></div>
        </div>
        <div class="defs-sources">
          <div class="defs-src-lbl">Sources</div>
          <div class="defs-src"><span class="defs-src-dot">·</span><span><b>Definition</b> — Leszczynski JF &amp; Rose GD (1986). Loops in globular proteins. <i>Science</i> 234:849–855. <a href="https://doi.org/10.1126/science.3775366" target="_blank">doi:10.1126/science.3775366</a></span></div>
        </div>
      </div>

    </div><!-- .defs-body -->
  </div><!-- .defs-panel -->
</div><!-- #defs-overlay -->
<button id="theme-toggle" title="Toggle dark / light mode">
  <span class="icon" id="theme-icon">☀️</span>
  <span id="theme-label">light</span>
</button>
<div id="hb-bg"></div>
<div class="wrap">
  <header>
    <img src="__APP_ROOT__/motif-logo" alt="Motif Scanner Logo" class="motif-logo">
    <h1>Motif Scanner</h1>
    <p class="subtitle">
      Detects sequence motifs (AXXXA, GXXXG, LXXXL, IXXXI, VXXXV, KXXXK, EXXXE) and structural motifs including β-turns, γ-turns, α-turns, and Ω-loops.
    </p>
    <button class="info-btn" onclick="openDefs()" title="View motif definitions and sources">? Definitions</button>
  </header>

  <input type="file" id="file" accept=".pdb,.cif,.mmcif,.txt,.fasta,.fa,.faa" multiple style="display:none">
  <div class="drop" id="drop" role="button" tabindex="0">
    <div class="glyph">[ + ]</div>
    <h2>Drop files here</h2>
    <p>Upload files: .pdb, .cif, .mmcif, .txt, .fasta, .fa </p>
  </div>
  <div class="pdb-fetch-panel">
    <div class="pdb-fetch-label">or fetch from the Protein Data Bank</div>
    <div class="pdb-fetch-row">
      <input type="text" id="pdb-ids-input"
             placeholder="4-letter PDB IDs, comma-separated · e.g. 1TIM, 2PIA, 4HHB"
             autocomplete="off" autocorrect="off" spellcheck="false" />
      <button class="fetch-btn" id="fetch-pdb-btn">Fetch</button>
    </div>
    <div class="pdb-fetch-status" id="pdb-fetch-status"></div>
  </div>
  <div class="demo-row">
    <button class="demo-btn" onclick="loadDemo()">
      <span class="demo-icon">▶</span>
      <span>
        Try the demo
        <span class="demo-sub">7OM2 · 684 residues · guided tour included</span>
      </span>
    </button>
  </div>
  <div class="files" id="files"></div>
  <div class="drop-err" id="drop-err"></div>

  <div class="bar">
    <button class="csv" id="csvStructMotifs"  hidden>↓ Sequence-Motifs.csv</button>
    <button class="csv" id="csvStructCounts"  hidden>↓ Sequence-Motifs-Counts.csv</button>
    <button class="csv" id="csvOmegaLoops"    hidden>↓ Omega-Loops.csv</button>
    <button class="csv" id="csvBetaTurns"     hidden>↓ Beta-Turns.csv</button>
    <button class="csv" id="csvGammaTurns"    hidden>↓ Gamma-Turns.csv</button>
    <button class="csv" id="csvAlphaTurns"    hidden>↓ Alpha-Turns.csv</button>
    <button class="csv" id="csvSeqMotifs"     hidden>↓ Seq-Motifs.csv</button>
    <button class="csv" id="csvSeqCounts"     hidden>↓ Seq-Counts.csv</button>
    <button class="csv" id="csvZipAll"        hidden style="display:none"></button>
    <span class="err" id="err"></span>
  </div>
  <div class="email-row">
    <label><input type="checkbox" id="auto-zip-chk"> Auto-download ZIP when scan finishes</label>
    <span class="email-sep">or</span>
    <input type="email" id="email-input" placeholder="Enter email to receive results" autocomplete="email" />
    <span id="email-status"></span>
  </div>
  <div class="scan-row">
    <button class="run" id="run" disabled>Scan</button>
  </div>
  <div id="warnings"></div>
  <div id="out"></div>
</div>

<script>
// ── haemoglobin background ───────────────────────────────────────────────
(function initHbBackground(){
  const container = document.getElementById('hb-bg');
  container.style.width  = '100vw';
  container.style.height = '100vh';

  let bgViewer = null;
  let bgRafId  = null;
  let bgLast   = 0;
  const BG_FPS      = 24;                        // ~24 renders/sec
  const BG_INTERVAL = 1000 / BG_FPS;
  const BG_ANGLE    = 0.18 * (60 / BG_FPS);     // same apparent speed at lower fps

  function bgPause() { if(bgRafId){cancelAnimationFrame(bgRafId);bgRafId=null;} }
  function bgResume(){ if(!bgRafId) bgRafId=requestAnimationFrame(bgRotate); }
  window.bgPause  = bgPause;
  window.bgResume = bgResume;

  // free the GPU/compositor from the continuous rotation while the user is
  // actively scrolling, so fixed-position elements don't compete with it
  let scrollIdleTimer = null;
  window.addEventListener('scroll', function(){
    bgPause();
    clearTimeout(scrollIdleTimer);
    scrollIdleTimer = setTimeout(bgResume, 200);
  }, {passive:true});

  function bgRotate(t) {
    bgRafId = requestAnimationFrame(bgRotate);
    if (t - bgLast < BG_INTERVAL) return;   // skip frames to hit target fps
    bgLast = t;
    bgViewer.rotate(BG_ANGLE, {x:0.3, y:1, z:0.1});
    bgViewer.render();
  }

  try {
    bgViewer = $3Dmol.createViewer(container, {
      backgroundColor: document.documentElement.getAttribute('data-theme') === 'dark' ? '#0f0f0f' : 'white',
      antialias: false,
      id: 'hb-viewer',
    });
    container._viewer = bgViewer;
  } catch(e) { return; }

  fetch('https://files.rcsb.org/download/2HHB.pdb')
    .then(r => r.ok ? r.text() : Promise.reject(r.status))
    .then(pdbData => {
      bgViewer.addModel(pdbData, 'pdb');
      bgViewer.setStyle({}, {cartoon:{color:'spectrum', opacity:0.85, thickness:0.6}});
      bgViewer.zoomTo();
      bgViewer.render();
      container.classList.add('loaded');
      bgRafId = requestAnimationFrame(bgRotate);   // start throttled loop
    })
    .catch(() => { container.style.display = 'none'; });

  // Pause when tab is hidden
  document.addEventListener('visibilitychange', () => {
    document.hidden ? bgPause() : bgResume();
  });
})();

// ── dark mode ────────────────────────────────────────────────────────────
(function initTheme(){
  const root    = document.documentElement;
  const btn     = document.getElementById('theme-toggle');
  const icon    = document.getElementById('theme-icon');
  const label   = document.getElementById('theme-label');
  const DARK    = 'dark';
  const LIGHT   = 'light';

  // Restore saved preference, defaulting to light mode
  const saved   = localStorage.getItem('motif-theme');
  let current   = saved || LIGHT;

  function apply(theme){
    current = theme;
    root.setAttribute('data-theme', theme);
    icon.textContent  = theme === DARK ? '🌙' : '☀️';
    label.textContent = theme === DARK ? 'dark' : 'light';
    localStorage.setItem('motif-theme', theme);

    // Update 3Dmol background viewers if they exist
    const bgColor = theme === DARK ? '#0f0f0f' : 'white';
    const viewerBg = theme === DARK ? '#1a1a1a' : '#f0efe9';
    // Background haemoglobin viewer
    const hbEl = document.getElementById('hb-bg');
    if(hbEl && hbEl._viewer){
      try{ hbEl._viewer.setBackgroundColor(bgColor); hbEl._viewer.render(); }catch(e){}
    }
    // Main 3Dmol panel viewer
    if(window.viewer){
      try{ window.viewer.setBackgroundColor(viewerBg); window.viewer.render(); }catch(e){}
    }
  }

  btn.addEventListener('click', () => apply(current === DARK ? LIGHT : DARK));
  apply(current);
})();

// ── colour maps ──────────────────────────────────────────────────────────
const CSS_VAR   = {A:'--A',G:'--G',E:'--E',I:'--I',L:'--L',K:'--K',V:'--V'};
const HEX_COLOR = {A:'#e6b450',G:'#7bd88f',E:'#ff6f6f',I:'#56c8e0',
                   L:'#b98cff',K:'#5b9cff',V:'#ff7bd0'};

// ── state ────────────────────────────────────────────────────────────────
let queue = [];
let structMotifs=[], structCounts=[], seqMotifs=[], seqCounts=[], omegaLoops=[];
let betaTurns=[], gammaTurns=[], alphaTurns=[];
let fileContents = {};
let viewer = null;

// ── DOM ──────────────────────────────────────────────────────────────────
const $  = s => document.querySelector(s);
const drop=$('#drop'), input=$('#file'), filesBox=$('#files'),
      runBtn=$('#run'), errEl=$('#err'), out=$('#out'), warnBox=$('#warnings'),
      dropErr=$('#drop-err'),
      csvSM=$('#csvStructMotifs'), csvSC=$('#csvStructCounts'),
      csvOL=$('#csvOmegaLoops'),
      csvBT=$('#csvBetaTurns'),   csvGT=$('#csvGammaTurns'),  csvAT=$('#csvAlphaTurns'),
      csvQM=$('#csvSeqMotifs'),   csvQC=$('#csvSeqCounts'),
      csvZipAll=$('#csvZipAll');

const STRUCT_EXT = new Set(['pdb','cif','mmcif']);
const SEQ_EXT    = new Set(['txt','fasta','fa','faa']);
const ALL_EXT    = new Set([...STRUCT_EXT, ...SEQ_EXT]);
function fileKind(f){ const e=f.name.split('.').pop().toLowerCase(); return STRUCT_EXT.has(e)?'struct':SEQ_EXT.has(e)?'seq':null; }
function colorVar(t){ return CSS_VAR[t[0]]||'--accent'; }
function colorHex(t){ return HEX_COLOR[t[0]]||'#7fd3c1'; }

// ── file queue ────────────────────────────────────────────────────────────
function renderFiles(){
  filesBox.innerHTML='';
  // invalid files shown as red chips (still in queue so user can see and remove)
  queue.forEach((f,i)=>{
    const kind=fileKind(f);
    const isValid = kind !== null;
    const c=document.createElement('span');
    c.className='chip'+(kind==='seq'?' chip-seq':!isValid?' chip-invalid':'');
    const badgeClass = isValid ? kind : 'inv';
    const badgeLabel = kind==='struct'?'INPUT':kind==='seq'?'SEQ':'ERR';
    c.innerHTML=`<span class="chip-badge ${badgeClass}">${badgeLabel}</span>`+
                `${f.name} <b data-i="${i}" title="remove">✕</b>`;
    filesBox.appendChild(c);
  });
  const invalidFiles = queue.filter(f=>fileKind(f)===null);
  const validFiles   = queue.filter(f=>fileKind(f)!==null);
  // show inline drop error if any invalid files present
  if(invalidFiles.length){
    const names = invalidFiles.map(f=>`· ${f.name}`).join('<br>');
    dropErr.innerHTML=`<b>unsupported file type — remove before scanning:</b><span>${names}</span>`;
    dropErr.classList.add('visible');
  } else {
    dropErr.innerHTML='';
    dropErr.classList.remove('visible');
  }
  runBtn.disabled = queue.length===0 || invalidFiles.length>0;
}
filesBox.addEventListener('click',e=>{
  if(e.target.dataset.i!==undefined){queue.splice(+e.target.dataset.i,1);renderFiles();}
});
function addFiles(list){
  for(const f of list){
    // accept all files into the queue (invalid ones shown as red chips)
    if(!queue.some(q=>q.name===f.name&&q.size===f.size)) queue.push(f);
  }
  renderFiles();
}
input.addEventListener('change',e=>addFiles(e.target.files));
drop.addEventListener('click',()=>input.click());
drop.addEventListener('keydown',e=>{if(e.key==='Enter'||e.key===' ')input.click();});
['dragenter','dragover'].forEach(ev=>drop.addEventListener(ev,e=>{e.preventDefault();drop.classList.add('hot');}));
['dragleave','drop'].forEach(ev=>drop.addEventListener(ev,e=>{e.preventDefault();drop.classList.remove('hot');}));
drop.addEventListener('drop',e=>addFiles(e.dataTransfer.files));

// ── PDB ID fetch ──────────────────────────────────────────────────────────
(function(){
  const fetchBtn  = document.getElementById('fetch-pdb-btn');
  const pdbInput  = document.getElementById('pdb-ids-input');
  const statusEl  = document.getElementById('pdb-fetch-status');

  function setStatus(msg, cls){
    statusEl.textContent = msg;
    statusEl.className = 'pdb-fetch-status' + (cls ? ' '+cls : '');
  }

  async function fetchPdbIds(){
    const raw = pdbInput.value.trim();
    if(!raw){ setStatus('Enter at least one PDB ID.', 'err'); return; }

    const tokens = raw.split(/[\s,;]+/).map(s=>s.trim().toUpperCase()).filter(Boolean);
    const valid   = tokens.filter(s=>/^[A-Z0-9]{4}$/.test(s));
    const bad     = tokens.filter(s=>!/^[A-Z0-9]{4}$/.test(s));

    if(bad.length){
      setStatus('Invalid IDs (must be 4 characters): '+bad.join(', '), 'err'); return;
    }
    if(!valid.length){ setStatus('No valid PDB IDs found.', 'err'); return; }

    fetchBtn.disabled = true;
    setStatus(`Fetching ${valid.length} structure(s) from RCSB PDB…`);

    const results = await Promise.all(valid.map(async id => {
      try{
        const r = await fetch(`https://files.rcsb.org/download/${id}.pdb`);
        if(!r.ok) throw new Error(`HTTP ${r.status}`);
        const text = await r.text();
        if(!text || text.trim().length < 10) throw new Error('empty response');
        return {id, text, err: null};
      } catch(e){
        return {id, text: null, err: e.message};
      }
    }));

    const ok     = results.filter(r=>!r.err);
    const failed = results.filter(r=>r.err);

    for(const r of ok){
      const file = new File([r.text], `${r.id}.pdb`, {type:'text/plain'});
      addFiles([file]);
    }

    if(failed.length && ok.length){
      setStatus(`Fetched ${ok.length} — failed: ${failed.map(r=>r.id+' ('+r.err+')').join(', ')}`, 'err');
    } else if(failed.length){
      setStatus(`Failed to fetch: ${failed.map(r=>r.id+' ('+r.err+')').join(', ')}`, 'err');
    } else {
      setStatus(`Fetched ${ok.length} structure(s): ${ok.map(r=>r.id).join(', ')}`, 'ok');
      pdbInput.value = '';
    }
    fetchBtn.disabled = false;
  }

  fetchBtn.addEventListener('click', fetchPdbIds);
  pdbInput.addEventListener('keydown', e=>{ if(e.key==='Enter') fetchPdbIds(); });
})();

// ── helpers ───────────────────────────────────────────────────────────────
function readFileText(f){
  return new Promise((res,rej)=>{const r=new FileReader();r.onload=()=>res(r.result);r.onerror=()=>rej(r.error);r.readAsText(f);});
}
const AA_NAMES = {
  A:'Alanine',R:'Arginine',N:'Asparagine',D:'Aspartate',C:'Cysteine',
  E:'Glutamate',Q:'Glutamine',G:'Glycine',H:'Histidine',I:'Isoleucine',
  L:'Leucine',K:'Lysine',M:'Methionine',F:'Phenylalanine',P:'Proline',
  S:'Serine',T:'Threonine',W:'Tryptophan',Y:'Tyrosine',V:'Valine'
};
function aaSpan(c, cls){
  const name = AA_NAMES[c.toUpperCase()] || c;
  return `<span class="${cls}" title="${c} · ${name}">${c}</span>`;
}
function seqHtml(sn){
  if(sn.length<2) return aaSpan(sn,'end');
  return aaSpan(sn[0],'end') + sn.slice(1,-1).split('').map(c=>aaSpan(c,'mid')).join('') + aaSpan(sn[sn.length-1],'end');
}
function turnSeqHtml(seq){
  if(!seq) return '';
  return seq.split('').map(c => {
    const name = AA_NAMES[c.toUpperCase()] || c;
    return `<span title="${c} · ${name}">${c}</span>`;
  }).join('');
}
function hexToRgb(hex){
  return {r:parseInt(hex.slice(1,3),16)/255,g:parseInt(hex.slice(3,5),16)/255,b:parseInt(hex.slice(5,7),16)/255};
}

// ── 3Dmol ─────────────────────────────────────────────────────────────────
function initViewer(container){
  if(viewer){try{viewer.clear();}catch(e){}}
  const isDark = document.documentElement.getAttribute('data-theme') === 'dark';
  viewer=$3Dmol.createViewer(container,{backgroundColor: isDark ? '#1a1a1a' : '#f0efe9',antialias:true});
  return viewer;
}
function loadStructure(fileData, motifsForStruct, focusMotif){
  const container=$('#mol-container');
  if(!container) return;
  const v=initViewer(container);
  v.addModel(fileData.data, fileData.format);
  v.setStyle({},{cartoon:{color:'#adb5bd',opacity:0.85}});
  motifsForStruct.forEach(m=>{
    const [s,e]=m.range.split('-').map(Number);
    const hex=colorHex(m.seq_type);
    for(let resi=s;resi<=e;resi++){
      v.addStyle({chain:m.chain_id,resi},{stick:{color:hex,radius:0.18},sphere:{color:hex,radius:0.28,opacity:0.85}});
    }
  });
  if(focusMotif){
    const [s,e]=focusMotif.range.split('-').map(Number);
    v.zoomTo({chain:focusMotif.chain_id,resi:Array.from({length:e-s+1},(_,i)=>s+i)});
  } else { v.zoomTo(); }
  v.render();
}
function loadForStructure(pdbId, focusMotif){
  const fd=fileContents[pdbId]; if(!fd) return;
  loadStructure(fd, structMotifs.filter(m=>m.pdb_id===pdbId), focusMotif);
}

// ── viewer panel HTML ─────────────────────────────────────────────────────
function viewerPanelHtml(structIds, activeId){
  const opts=structIds.map(id=>`<option value="${id}"${id===activeId?' selected':''}>${id}</option>`).join('');

  let legend='';

  // Correctly detect which top-level tab is active
  const structuralTabActive =
    document.querySelector('.tab-btn[data-tab="struct-motifs"].active') !== null;

  if(structuralTabActive){

    const items=[];

    if(typeof alphaTurns!=='undefined' && alphaTurns.length)
      items.push(`<span class="legend-pill" style="border-color:#b71c1c22;color:#b71c1c"><span class="legend-dot" style="background:#b71c1c"></span>α-turn</span>`);

    if(typeof betaTurns!=='undefined' && betaTurns.length)
      items.push(`<span class="legend-pill" style="border-color:#7c3aed22;color:#7c3aed"><span class="legend-dot" style="background:#7c3aed"></span>β-turn</span>`);

    if(typeof gammaTurns!=='undefined' && gammaTurns.length)
      items.push(`<span class="legend-pill" style="border-color:#1a6b8a22;color:#1a6b8a"><span class="legend-dot" style="background:#1a6b8a"></span>γ-turn</span>`);

    if(typeof omegaLoops!=='undefined' && omegaLoops.length)
      items.push(`<span class="legend-pill" style="border-color:#e05c0022;color:#e05c00"><span class="legend-dot" style="background:#e05c00"></span>Ω-loop</span>`);

    legend = items.join('');
  }else{
    const presentTypes=[...new Set(structMotifs.map(m=>m.seq_type))];

    legend=presentTypes.map(type=>{
      const k=type[0];
      const hex=HEX_COLOR[k] || '#666';
      return `<span class="legend-pill" style="border-color:${hex}22;color:${hex}">
        <span class="legend-dot" style="background:${hex}"></span>${type}</span>`;
    }).join('');
  }

  return `<div class="viewer-panel" id="viewer-panel">
    <div class="viewer-header">
      <span class="viewer-title">3D · structure</span>
      <select class="viewer-select" id="struct-select">${opts}</select>
      <button class="fs-open-btn" onclick="openFsViewer(document.getElementById('struct-select')?.value)" title="Open fullscreen explorer (3D Explorer)">↗</button>
    </div>
    <div id="mol-container"></div>
    <div class="viewer-legend">${legend}</div>
    <div class="viewer-hint">scroll to zoom · drag to rotate · click motif row to focus</div>
  </div>`;
}

// ── stat cards ────────────────────────────────────────────────────────────
function statsHtml(motifs, counts, label, accentVar){
  const byType={};
  counts.forEach(c=>{byType[c.seq_type]=(byType[c.seq_type]||0)+c.count;});
  const typeStats=Object.keys(byType).map(t=>{
    const cv=colorVar(t);
    return `<div class="stat" style="--c:var(${cv})">
      <div class="k" style="color:var(${cv})">${t}</div>
      <div class="v">${byType[t]}</div></div>`;
  }).join('');
  const total=motifs.length;
  const srcCount=[...new Set(motifs.map(m=>m.pdb_id))].length;
  return `<div class="summary">
    <div class="stat" style="--c:var(${accentVar})"><div class="k">total motifs</div><div class="v">${total}</div></div>
    <div class="stat" style="--c:var(${accentVar})"><div class="k">${label}</div><div class="v">${srcCount}</div></div>
    ${typeStats}</div>`;
}

// ── motif list HTML ───────────────────────────────────────────────────────
const SS_LABELS = {H:'α',E:'β',B:'β',G:'3',I:'π',T:'T',S:'S','-':' '};
function ssHtml(ss){
  if(!ss) return `<span class="ss-col ss-na" title="DSSP not available">n/a</span>`;
  const spans = ss.split('').map(c=>{
    const cls = 'ss-'+(c==='-'?'-':c);
    const title = {H:'α-helix',E:'β-strand',B:'β-bridge',T:'turn',S:'bend',G:'3₁₀-helix',I:'π-helix','-':'coil'}[c]||c;
    return `<span class="${cls}" title="${title}">${c}</span>`;
  }).join('');
  return `<span class="ss-col">${spans}</span>`;
}

function ifaceHtml(iface, rangeStr){
  if(!iface) return `<span class="iface-col iface-na" title="FreeSASA unavailable or single chain">—</span>`;
  const hasAny = iface.includes('I');
  if(!hasAny) return `<span class="iface-col iface-na" title="No interface residues in this window">·····</span>`;

  // Parse start residue number from range string e.g. "200-204"
  const startRes = parseInt((rangeStr||'0').split('-')[0], 10);

  // Collect residue numbers for each 'I' position
  const ifaceNums = iface.split('').map((c,i) => c==='I' ? startRes+i : null).filter(n=>n!==null);

  const badges = ifaceNums.map(n =>
    `<span class="iface-num" title="Residue ${n}: interface">${n}</span>`
  ).join('');
  return `<span class="iface-col">${badges}</span>`;
}

function motifTypeTip(seqType){
  const c = (seqType||'')[0];
  const name = AA_NAMES[c] || c;
  return `${seqType} — ${name}-X-X-X-${name} motif`;
}

function motifListHtml(motifs, extraClass='', showChain=true){
  if(motifs.length===0) return `<div class="empty">No X***X motifs found.</div>`;

  const isStruct = extraClass === 'struct-row';

  // Determine whether any motif in this batch has iface data
  const hasIface = isStruct && motifs.some(m => typeof m.iface === 'string');
  const legendHtml = '';

  // Column header row — matches the 6-col grid exactly
  const hdrs = isStruct
    ? `<div class="col-headers">
        <span class="col-hdr">#</span>
        <span class="col-hdr">Motif type</span>
        <span class="col-hdr">Sequence</span>
        <span class="col-hdr">Structure · Chain · Range</span>
        <span class="col-hdr">Secondary structure</span>
        <span class="col-hdr">Interface residue</span>
      </div>`
    : `<div class="col-headers">
        <span class="col-hdr">#</span>
        <span class="col-hdr">Motif type</span>
        <span class="col-hdr">Sequence</span>
        <span class="col-hdr">Source · Position</span>
      </div>`;

  return `<div class="motif-list-wrap">`+legendHtml + hdrs + '<ol class="list">'+motifs.map((m,i)=>{
    const cv=colorVar(m.seq_type);
    const chainLine = (showChain && m.chain_id && m.chain_id!=='—')
      ? `chain ${m.chain_id} · 5 res · ${m.range}`
      : `5 res · ${m.range}`;
    const metaTip = isStruct
      ? `PDB: ${m.pdb_id} · Chain: ${m.chain_id||'—'} · Range: ${m.range}`
      : `Source: ${m.pdb_id} · Position: ${m.range}`;
    const ssStr = m.ss||'';
    const ssTip = ssStr ? `DSSP secondary structure (${ssStr.length} residues): ${ssStr}` : 'DSSP not available';
    const ssCol   = isStruct ? `<span title="${ssTip}">${ssHtml(ssStr)}</span>` : '';
    const ifaceCol = isStruct ? ifaceHtml(m.iface||'', m.range||'') : '';
    return `<li class="row ${extraClass}" style="--c:var(${cv});animation-delay:${Math.min(i*14,400)}ms" data-idx="${i}">
      <span class="idx" title="Motif #${i+1}">${i+1}</span>
      <span class="type" title="${motifTypeTip(m.seq_type)}">${m.seq_type}</span>
      <span class="seq-str">${seqHtml(m.sn)}</span>
      <span class="meta" title="${metaTip}"><b>${m.pdb_id}</b><span class="meta-sep">·</span><span class="meta-sub">${chainLine}</span></span>
      ${ssCol}
      ${ifaceCol}
    </li>`;
  }).join('')+'</ol></div>';
}

// ── omega loop list HTML ──────────────────────────────────────────────────
function omegaListHtml(loops){
  if(!loops || loops.length===0) return `<div class="empty">No Ω-loops found in this structure.</div>`;

  const hdrs = `<div class="omega-col-headers">
    <span class="col-hdr" style="text-align:left">#</span>
    <span class="col-hdr">PDB · Chain</span>
    <span class="col-hdr">Sequence</span>
    <span class="col-hdr">Secondary structure</span>
    <span class="col-hdr">End dist (Å)</span>
    <span class="col-hdr">Radius · Area · Volume</span>
  </div>`;

  const rows = loops.map((m, i) => {
    // Flanking SS for tooltip only
    const ssBefore = m.Sec_Structure_before || '—';
    const ssAfter  = m.Sec_Structure_After  || '—';
    // Use same ssHtml() as sequence motifs tab — identical coloured tile rendering
    const ssHtmlStr = `<span title="Before: ${ssBefore} · Loop · After: ${ssAfter}">${ssHtml(m.Sec_Structure||'')}</span>`;

    const dist   = m.distance_x !== null && m.distance_x !== undefined ? parseFloat(m.distance_x).toFixed(3) : '?';
    const radius = m.radius     !== null && m.radius     !== undefined ? parseFloat(m.radius).toFixed(3)     : '?';
    const area   = m.Area       !== null && m.Area       !== undefined ? parseFloat(m.Area).toFixed(3)       : '?';
    const vol    = m.Volume     !== null && m.Volume     !== undefined ? parseFloat(m.Volume).toFixed(3)     : '?';

    const geoVals = `<span class="geo-vals" title="Loop geometry · Radius: ${radius} Å · Area: ${area} Å² · Volume: ${vol} Å³">
      <span title="Radius: ${radius} Å"><span class="geo-lbl">Radius</span><b>${radius}</b> Å</span>
      <span title="Area: ${area} Å²"><span class="geo-lbl">Area</span><b>${area}</b> Å²</span>
      <span title="Volume: ${vol} Å³"><span class="geo-lbl">Volume</span><b>${vol}</b> Å³</span>
    </span>`;

    return `<li class="omega-row" data-idx="${i}" data-pdb="${m.pdb_id}" data-chain="${m.Chain}" data-start="${m.Start_Position}" data-end="${m.End_Position}" style="animation-delay:${Math.min(i*12,400)}ms">
      <span class="idx">${i+1}</span>
      <span class="omega-id">
        <span class="oid-pdb">${m.pdb_id}</span>
        <span class="oid-loc">ch.${m.Chain}&nbsp;·&nbsp;${m.Start_Position}–${m.End_Position}</span>
      </span>
      <span class="seq-str" style="font-size:17px;letter-spacing:.22em">${turnSeqHtml(m.Amino_Acid_Seq||'')}</span>
      ${ssHtmlStr}
      <span class="omega-dist" title="Cα–Cα distance between terminal residues i and i+5: ${dist} Å"><b>${dist}</b>&nbsp;Å</span>
      ${geoVals}
    </li>`;
  }).join('');

  return `<div class="omega-list-wrap">${hdrs}<ol class="list">${rows}</ol></div>`;
}

// ── β/γ turn list HTML ────────────────────────────────────────────────────
function turnSubtypeHtml(t, turnType){
  if(!t) return '';
  if(turnType==='alpha'){
    if(!t) return '<span class="turn-subtype" style="background:#b71c1c18;color:#b71c1c;border:1px solid #b71c1c33" title="α-turn">α-turn</span>';
    const label = t==='isolated' ? 'α isolated' : t==='helix-cap' ? 'α helix-cap' : 'α-turn';
    const tip   = t==='isolated' ? 'α-turn isolated' : t==='helix-cap' ? 'α-turn helix-cap' : 'α-turn';
    return `<span class="turn-subtype" style="background:#b71c1c18;color:#b71c1c;border:1px solid #b71c1c33" title="${tip}">${label}</span>`;
  }
  if(t==='classic')       return '<span class="turn-subtype gamma-cls"     title="Classic γ-turn">classic γ-turn</span>';
  if(t==='inverse')       return '<span class="turn-subtype gamma-inv"     title="Inverse γ-turn">inverse γ-turn</span>';
  if(t==='classical-nhn') return '<span class="turn-subtype gamma-cls-nhn" title="Classical NHN γ-turn">classical NHN γ-turn</span>';
  if(t==='inverse-nhn')   return '<span class="turn-subtype gamma-inv-nhn" title="Inverse NHN γ-turn">inverse NHN γ-turn</span>';
  if(t===null||t===undefined||t==='') return '';
  if(t==='IV')            return '<span class="turn-subtype unc"  title="β-turn type IV">β-IV</span>';
  return `<span class="turn-subtype beta" title="β-turn type ${t}">β-${t}</span>`;
}

function fmt1(v){ return (v!=null && v!==undefined) ? parseFloat(v).toFixed(1) : '—'; }

function turnAnglesHtml(turn){
  if(turn.turn_type==='beta'){
    const phi1=fmt1(turn.phi_i1), psi1=fmt1(turn.psi_i1);
    const phi2=fmt1(turn.phi_i2), psi2=fmt1(turn.psi_i2);
    return `<span class="turn-angles" title="Ramachandran torsion angles · i+1: φ ${phi1}° ψ ${psi1}° · i+2: φ ${phi2}° ψ ${psi2}°">
      <span class="angle-row" title="Residue i+1 torsion angles: φ ${phi1}° ψ ${psi1}°"><span class="alabel">i+1</span>φ <b>${phi1}</b>° ψ <b>${psi1}</b>°</span>
      <span class="angle-row" title="Residue i+2 torsion angles: φ ${phi2}° ψ ${psi2}°"><span class="alabel">i+2</span>φ <b>${phi2}</b>° ψ <b>${psi2}</b>°</span>
    </span>`;
  } else if(turn.turn_type==='gamma'){
    const phi=fmt1(turn.phi_mid), psi=fmt1(turn.psi_mid);
    const isNhn = turn.is_nhn === true;
    if(isNhn){
      const hod   = turn.ho_dist   != null ? parseFloat(turn.ho_dist).toFixed(2)   : '—';
      const theta = turn.hbond_dist != null ? parseFloat(turn.hbond_dist).toFixed(1) : '—';
      return `<span class="turn-angles" title="Middle residue torsion angles: φ ${phi}° ψ ${psi}° · H···O ${hod} Å · θ ${theta}°">
        <span class="angle-row" title="Middle residue torsion angles: φ ${phi}° ψ ${psi}°"><span class="alabel">mid</span>φ <b>${phi}</b>° ψ <b>${psi}</b>°</span>
        <span class="angle-row" title="N-H(i)···O=C(i) C5 ring distance (criterion c, must be < 3.0 Å): ${hod} Å"><span class="alabel" style="color:#4caf50">H·O</span><b style="color:#4caf50">${hod}</b><span style="color:var(--ink-faint);font-size:10px"> Å</span></span>
        <span class="angle-row" title="NHN theta: angle between N-H(i+2) and N(i+1) pz orbital (criterion b, must be < 30°): ${theta}°"><span class="alabel" style="color:#9c6fd4">θ</span><b style="color:#6d3a9a">${theta}</b><span style="color:var(--ink-faint);font-size:10px">°</span></span>
      </span>`;
    }
    const hbd = turn.hbond_dist != null ? parseFloat(turn.hbond_dist).toFixed(2) : '—';
    const hbLabel = turn.h_present ? 'O·H' : 'N·O';
    const hbTitle = turn.h_present
      ? 'O···H distance (full geometric H-bond check)'
      : 'N···O heavy-atom distance (H absent — fallback criterion)';
    return `<span class="turn-angles" title="Middle residue torsion angles: φ ${phi}° ψ ${psi}° · H-bond ${hbLabel} ${hbd} Å">
      <span class="angle-row" title="Middle residue torsion angles: φ ${phi}° ψ ${psi}°"><span class="alabel">mid</span>φ <b>${phi}</b>° ψ <b>${psi}</b>°</span>
      <span class="angle-row" title="${hbTitle}: ${hbd} Å"><span class="alabel" style="color:#4caf50">${hbLabel}</span><b style="color:#4caf50">${hbd}</b><span style="color:var(--ink-faint);font-size:10px"> Å</span></span>
    </span>`;
  } else {
    // alpha
    const hbd = turn.hbond_NO_dist!=null ? parseFloat(turn.hbond_NO_dist).toFixed(2) : '—';
    const phi2=fmt1(turn.phi_i2), psi2=fmt1(turn.psi_i2);
    return `<span class="turn-angles" title="Residue i+2 torsion angles: φ ${phi2}° ψ ${psi2}° · N(i+4)···O(i) H-bond: ${hbd} Å">
      <span class="angle-row" title="Torsion angles of central residue i+2: φ ${phi2}° ψ ${psi2}°"><span class="alabel">i+2</span>φ <b>${fmt1(turn.phi_i2)}</b>° ψ <b>${fmt1(turn.psi_i2)}</b>°</span>
      <span class="angle-row" title="N(i+4)···O(i) H-bond distance: ${hbd} Å"><span class="alabel" style="color:#e57373">N·O</span><b style="color:#b71c1c">${hbd}</b><span style="color:var(--ink-faint);font-size:10px"> Å</span></span>
    </span>`;
  }
}

function turnSectionHtml(turns, turnType, color, label, defn){
  if(!turns.length) return `<div class="empty" style="margin-bottom:16px">No ${label} found in this structure.</div>`;

  const MAX_ROWS = 300;
  const capped   = turns.length > MAX_ROWS;
  const shown    = capped ? turns.slice(0, MAX_ROWS) : turns;
  const capNote  = capped
    ? `<div class="iface-legend" style="margin-bottom:10px;background:var(--panel-2);
        padding:8px 12px;border-radius:8px;border:1px solid var(--line)">
        Showing first ${MAX_ROWS} of ${turns.length} — download CSV for the full dataset.
       </div>`
    : '';

  const hdrs = `<div class="turn-col-headers">
    <span class="col-hdr" style="text-align:center">#</span>
    <span class="col-hdr">PDB · Chain</span>
    <span class="col-hdr">Sequence</span>
    <span class="col-hdr">TYPE</span>
    <span class="col-hdr">Secondary str.</span>
    <span class="col-hdr">End dist (Å)</span>
    <span class="col-hdr">Torsion angles</span>
  </div>`;
  const rows = shown.map((t, i) => {
    const dist = (t.distance_x!=null) ? parseFloat(t.distance_x).toFixed(3) : '?';
    const asa  = (t.total_asa!=null)  ? parseFloat(t.total_asa).toFixed(1)  : null;
    const asaTip = asa ? `· ASA ${asa} Å²` : '';
    const distLabel = turnType==='beta'
      ? `Cα–Cα distance between terminal residues i and i+3 (criterion: < 7 Å): ${dist} Å`
      : turnType==='gamma'
      ? `Cα–Cα distance between terminal residues i and i+2: ${dist} Å`
      : `Cα–Cα distance between terminal residues: ${dist} Å`;
    const ssStr = t.Sec_Structure||'';
    const ssTip = ssStr ? `DSSP secondary structure (${ssStr.length} residues): ${ssStr}` : 'DSSP not available';
    return `<li class="turn-row ${turnType}-turn" data-type="${turnType}" data-idx="${i}"
         data-pdb="${t.pdb_id}" data-chain="${t.Chain}"
         data-start="${t.Start_Position}" data-end="${t.End_Position}"
         style="animation-delay:${Math.min(i*10,400)}ms"
         title="${t.pdb_id} ch.${t.Chain} ${t.Start_Position}–${t.End_Position}${asaTip}">
      <span class="idx">${i+1}</span>
      <span class="omega-id">
        <span class="oid-pdb">${t.pdb_id}</span>
        <span class="oid-loc">ch.${t.Chain}&nbsp;·&nbsp;${t.Start_Position}–${t.End_Position}</span>
      </span>
      <span class="seq-str" style="font-size:17px;letter-spacing:.26em">${turnSeqHtml(t.Amino_Acid_Seq||'')}</span>
      ${turnSubtypeHtml(t.subtype, t.turn_type)}
      <span title="${ssTip}">${ssHtml(ssStr)}</span>
      <span class="turn-dist" title="${distLabel}"><b>${dist}</b>&nbsp;Å</span>
      ${turnAnglesHtml(t)}
    </li>`;
  }).join('');
  return `<div class="turn-list-wrap">${capNote}${hdrs}<ol class="list">${rows}</ol></div>`;
}

// ── α-turn section with cluster column ───────────────────────────────────
function alphaTurnSectionHtml(alphaTurns){
  if(!alphaTurns.length) return `<div class="empty">No α-turns found in this structure.</div>`;

  const MAX_ROWS = 300;
  const shown    = alphaTurns.slice(0, MAX_ROWS);
  const capNote  = alphaTurns.length > MAX_ROWS
    ? `<div class="iface-legend" style="margin-bottom:10px;background:var(--panel-2);padding:8px 12px;border-radius:8px;border:1px solid var(--line)">
        Showing first ${MAX_ROWS} of ${alphaTurns.length} — download CSV for the full dataset.
       </div>` : '';

  // Cluster colour map
  const CLUSTER_COLORS = {
    'I-αRS':  '#c62828', 'I-αLS':  '#ad1457',
    'II-αRS': '#e65100', 'II-αLS': '#f57f17',
    'I-αLU':  '#2e7d32', 'I-αRU':  '#1565c0',
    'II-αLU': '#6a1b9a', 'II-αRU': '#00695c',
    'I-αc':   '#4e342e',
  };

  const hdrs = `<div class="turn-col-headers" >
    <span class="col-hdr" style="text-align:center">#</span>
    <span class="col-hdr">PDB · Chain</span>
    <span class="col-hdr">Sequence</span>
    <span class="col-hdr">TYPE</span>
    <span class="col-hdr">Secondary str.</span>
    <span class="col-hdr">End dist (Å)</span>
    <span class="col-hdr">Torsion angles (i+1, i+2, i+3)</span>
  </div>`;

  const rows = shown.map((t, i) => {
    const dist   = t.distance_x!=null ? parseFloat(t.distance_x).toFixed(3) : '?';
    const hbd    = t.hbond_NO_dist!=null ? parseFloat(t.hbond_NO_dist).toFixed(2) : '—';
    const cl     = (t.cluster && t.cluster !== 'unclassified') ? t.cluster : null;
    const clColor = cl ? (CLUSTER_COLORS[cl] || '#555') : '#aaa';
    const clBadge = cl
      ? `<span class="turn-subtype" style="background:${clColor}18;color:${clColor};border:1px solid ${clColor}44;font-weight:600" title="α-turn Ramachandran cluster: ${cl}">${cl}</span>`
      : `<span style="color:var(--ink-faint);font-size:11px" title="α-turn cluster: unclassified">—</span>`;
    const phi1 = t.phi_i1!=null?parseFloat(t.phi_i1).toFixed(1):'—';
    const psi1 = t.psi_i1!=null?parseFloat(t.psi_i1).toFixed(1):'—';
    const phi2 = t.phi_i2!=null?parseFloat(t.phi_i2).toFixed(1):'—';
    const psi2 = t.psi_i2!=null?parseFloat(t.psi_i2).toFixed(1):'—';
    const phi3 = t.phi_i3!=null?parseFloat(t.phi_i3).toFixed(1):'—';
    const psi3 = t.psi_i3!=null?parseFloat(t.psi_i3).toFixed(1):'—';
    const ssStr = t.Sec_Structure||'';
    const ssTip = ssStr ? `DSSP secondary structure (${ssStr.length} residues): ${ssStr}` : 'DSSP not available';
    return `<li class="turn-row alpha-turn" style="animation-delay:${Math.min(i*10,400)}ms"
         data-type="alpha" data-idx="${i}"
         data-pdb="${t.pdb_id}" data-chain="${t.Chain}"
         data-start="${t.Start_Position}" data-end="${t.End_Position}"
         title="${t.pdb_id} ch.${t.Chain} ${t.Start_Position}–${t.End_Position}">
      <span class="idx">${i+1}</span>
      <span class="omega-id">
        <span class="oid-pdb">${t.pdb_id}</span>
        <span class="oid-loc">ch.${t.Chain}&nbsp;·&nbsp;${t.Start_Position}–${t.End_Position}</span>
      </span>
      <span class="seq-str" style="font-size:17px;letter-spacing:.26em">${turnSeqHtml(t.Amino_Acid_Seq||'')}</span>
      ${clBadge}
      <span title="${ssTip}">${ssHtml(ssStr)}</span>
      <span class="turn-dist" title="Cα–Cα distance between terminal residues i and i+4: ${dist} Å · N(i+4)···O(i) H-bond: ${hbd} Å"><b>${dist}</b>&nbsp;Å&nbsp;<span style="color:var(--ink-faint);font-size:10px">N·O ${hbd} Å</span></span>
      <span class="turn-angles" title="Ramachandran torsion angles · i+1: φ ${phi1}° ψ ${psi1}° · i+2: φ ${phi2}° ψ ${psi2}° · i+3: φ ${phi3}° ψ ${psi3}°">
        <span class="angle-row" title="Residue i+1 torsion angles: φ ${phi1}° ψ ${psi1}°"><span class="alabel">i+1</span>φ <b>${phi1}</b>° ψ <b>${psi1}</b>°</span>
        <span class="angle-row" title="Residue i+2 torsion angles: φ ${phi2}° ψ ${psi2}°"><span class="alabel">i+2</span>φ <b>${phi2}</b>° ψ <b>${psi2}</b>°</span>
        <span class="angle-row" title="Residue i+3 torsion angles: φ ${phi3}° ψ ${psi3}°"><span class="alabel">i+3</span>φ <b>${phi3}</b>° ψ <b>${psi3}</b>°</span>
      </span>
    </li>`;
  }).join('');

  return `<div class="turn-list-wrap">${capNote}${hdrs}<ol class="list">${rows}</ol></div>`;
}

function turnPanelHtml(betaTurns, gammaTurns, alphaTurns){
  const aCount = alphaTurns.length, bCount = betaTurns.length, gCount = gammaTurns.length;
  const statCards = `<div class="summary">
    <div class="stat" style="--c:#1a6b8a"><div class="k">γ-turns</div><div class="v">${gCount}</div></div>
    <div class="stat" style="--c:#7c3aed"><div class="k">β-turns</div><div class="v">${bCount}</div></div>
    <div class="stat" style="--c:#b71c1c"><div class="k">α-turns</div><div class="v">${aCount}</div></div>
    <div class="stat" style="--c:#555"><div class="k">structures</div><div class="v">${[...new Set([...gammaTurns,...betaTurns,...alphaTurns].map(t=>t.pdb_id))].length}</div></div>
  </div>`;
  return `<div class="iface-legend" style="margin-bottom:16px">
      <span class="swatch" style="background:#1a6b8a"></span>
      γ: 3-res · H-bond C=O(i)···N-H(i+2): O·H &lt;2.5 Å + angles if H present, else N·O &lt;3.5 Å · classic: φ,ψ within 40° of 75°,−64° · inverse: φ,ψ within 40° of −79°,69° · NHN γ (Dhar et al. 2019): Ca–Ca 5.2–5.8 Å, θ &lt;30°, C5 H···O &lt;3.0 Å · classical NHN: φ 46–72°, ψ 20–52° · inverse NHN: φ −90 to −54°, ψ −42 to −2°.&nbsp;&nbsp;
      <span class="swatch" style="background:#7c3aed;margin-left:8px"></span>
      β: 4-res · Cα(i)→Cα(i+3) &lt;7 Å · neither i+1 nor i+2 in {G,H,I} · H-bond not required · Types I I′ II II′ VIa VIb VIII IV.&nbsp;&nbsp;
      <span class="swatch" style="background:#b71c1c;margin-left:8px"></span>
      α: 5-res · DSSP H-bond NH(i+4)→C=O(i) &lt;−0.5 kcal/mol · middle 3 not all H (Toniolo &amp; Benedetti 1980).
    </div>
    ${statCards}
    <div class="section-h" style="--c:#1a6b8a;color:#1a6b8a">γ-turns · list</div>
    ${turnSectionHtml(gammaTurns, 'gamma', '#1a6b8a', 'γ-turns', '')}
    <div class="section-h" style="--c:#7c3aed;color:#7c3aed;margin-top:22px">β-turns · list</div>
    ${turnSectionHtml(betaTurns,  'beta',  '#7c3aed', 'β-turns', '')}
    <div class="section-h" style="--c:#b71c1c;color:#b71c1c;margin-top:22px">α-turns · list</div>
    ${turnSectionHtml(alphaTurns, 'alpha', '#b71c1c', 'α-turns', '')}`;
}

// ── loadForOmegaLoop — highlight a specific loop in the 3Dmol viewer ──────
function loadForOmegaLoop(loop){
  const fd = fileContents[loop.pdb_id];
  if(!fd) return;
  const container = $('#mol-container');
  if(!container) return;
  const v = initViewer(container);
  v.addModel(fd.data, fd.format);
  v.setStyle({}, {cartoon:{color:'#adb5bd', opacity:0.7}});
  // Highlight all omega loops for this structure faintly
  omegaLoops.filter(l=>l.pdb_id===loop.pdb_id).forEach(l=>{
    for(let resi=l.Start_Position; resi<=l.End_Position; resi++){
      v.addStyle({chain:l.Chain, resi}, {stick:{color:'#e0a87a', radius:0.12}});
    }
  });
  // Highlight the selected loop brightly
  for(let resi=loop.Start_Position; resi<=loop.End_Position; resi++){
    v.addStyle({chain:loop.Chain, resi}, {
      stick:  {color:'#e05c00', radius:0.22},
      sphere: {color:'#e05c00', radius:0.32, opacity:0.9},
    });
  }
  v.zoomTo({chain:loop.Chain, resi:Array.from({length:loop.End_Position-loop.Start_Position+1},(_,k)=>loop.Start_Position+k)});
  v.render();
}

// ── loadForTurn — highlight a β/γ-turn in the 3Dmol viewer ──────────────
function loadForTurn(turn){
  const fd = fileContents[turn.pdb_id];
  if(!fd) return;
  const container = $('#mol-container');
  if(!container) return;
  let color = '#b71c1c';
  let faintColor = '#e6a5a5';
  let allTurns = alphaTurns;

  if(turn.turn_type === 'beta'){
    color = '#7c3aed';
    faintColor = '#c9b0f0';
    allTurns = betaTurns;
  }else if(turn.turn_type === 'gamma'){
    color = '#1a6b8a';
    faintColor = '#8ab8cc';
    allTurns = gammaTurns;
  }
  const v = initViewer(container);
  v.addModel(fd.data, fd.format);
  v.setStyle({}, {cartoon:{color:'#adb5bd', opacity:0.65}});
  // faintly tint all other turns of same type
  allTurns.filter(t=>t.pdb_id===turn.pdb_id).forEach(t=>{
    for(let resi=t.Start_Position; resi<=t.End_Position; resi++){
      v.addStyle({chain:t.Chain, resi}, {stick:{color:faintColor, radius:0.10}});
    }
  });
  // highlight the selected turn
  for(let resi=turn.Start_Position; resi<=turn.End_Position; resi++){
    v.addStyle({chain:turn.Chain, resi}, {
      stick:  {color, radius:0.22},
      sphere: {color, radius:0.30, opacity:0.88},
    });
  }
  v.zoomTo({chain:turn.Chain, resi:Array.from(
    {length:turn.End_Position-turn.Start_Position+1},(_,k)=>turn.Start_Position+k
  )});
  v.render();
}

// ── main render ────────────────────────────────────────────────────────────
function render(data){
  structMotifs=data.struct_motifs; structCounts=data.struct_counts;
  seqMotifs=data.seq_motifs;       seqCounts=data.seq_counts;
  omegaLoops=data.omega_loops||[];
  betaTurns=data.beta_turns||[];   gammaTurns=data.gamma_turns||[];
  alphaTurns=data.alpha_turns||[];

  csvSM.hidden=structMotifs.length===0; csvSC.hidden=structCounts.length===0;
  csvOL.hidden=omegaLoops.length===0;
  csvBT.hidden=betaTurns.length===0;    csvGT.hidden=gammaTurns.length===0;
  csvAT.hidden=alphaTurns.length===0;
  csvQM.hidden=seqMotifs.length===0;    csvQC.hidden=seqCounts.length===0;
  const _anyResults = structMotifs.length||omegaLoops.length||betaTurns.length||
                      gammaTurns.length||alphaTurns.length||seqMotifs.length;
  csvZipAll.hidden = !_anyResults;

  // warnings
  warnBox.innerHTML='';
  if(data.warnings && data.warnings.length){
    warnBox.innerHTML=`<div class="warn-box"><b>⚠ warnings</b>`+
      data.warnings.map(w=>`<span>${w}</span>`).join('')+`</div>`;
  }

  let html='';

  // ── structure section ──────────────────────────────────────────────────
  if(data.struct_files > 0){
    const structIds=[...new Set(structMotifs.map(m=>m.pdb_id))];
    const activeId=structIds[0]||null;
    const hasViewer=Object.keys(fileContents).length>0;

    const motifCount  = structMotifs.length;
    const omegaCount  = omegaLoops.length;
    const betaCount   = betaTurns.length;
    const gammaCount  = gammaTurns.length;
    const alphaCount  = alphaTurns.length;
    const structMotifCount = omegaCount + betaCount + gammaCount + alphaCount;

    // Top-level tab bar
    const tabBarHtml=`
      <div class="tab-bar" id="struct-tab-bar">
        <button class="tab-btn active" data-tab="seq-motifs">
          Sequence motifs<span class="tab-count">${motifCount}</span>
        </button>
        <button class="tab-btn" data-tab="struct-motifs" style="--tab-accent:#2a7c5f">
          Structural motifs<span class="tab-count" style="${structMotifCount?'background:#2a7c5f;color:#fff':''}">${structMotifCount}</span>
        </button>
      </div>`;

    // Tab 1 — sequence motifs (unchanged)
    const seqMotifPanel=`
      <div class="tab-panel active" id="tab-panel-seq-motifs">
        ${statsHtml(structMotifs, structCounts, 'structures', '--accent')}
        <div class="section-h">motifs · list</div>
        ${motifListHtml(structMotifs, 'struct-row', true)}
      </div>`;

    // Tab 2 — structural motifs (nested sub-tabs for β, γ, α, Ω)
    const structMotifPanel=`
      <div class="tab-panel" id="tab-panel-struct-motifs">
        <div class="tab-bar sub-tab-bar" id="sub-tab-bar" style="margin-bottom:16px;border-bottom:1px solid var(--line);padding-bottom:0">
          <button class="tab-btn sub-tab-btn active" data-subtab="gamma" style="--tab-accent:#1a6b8a">
            γ-turns<span class="tab-count" style="${gammaCount?'background:#1a6b8a;color:#fff':''}">${gammaCount}</span>
          </button>
          <button class="tab-btn sub-tab-btn" data-subtab="beta" style="--tab-accent:#7c3aed">
            β-turns<span class="tab-count" style="${betaCount?'background:#7c3aed;color:#fff':''}">${betaCount}</span>
          </button>
          <button class="tab-btn sub-tab-btn" data-subtab="alpha" style="--tab-accent:#b71c1c">
            α-turns<span class="tab-count" style="${alphaCount?'background:#b71c1c;color:#fff':''}">${alphaCount}</span>
          </button>
          <button class="tab-btn sub-tab-btn" data-subtab="omega" style="--tab-accent:#e05c00">
            Ω-loops<span class="tab-count" style="${omegaCount?'background:#e05c00;color:#fff':''}">${omegaCount}</span>
          </button>
        </div>
        <div class="sub-tab-panel active" id="sub-tab-panel-gamma">
          ${turnSectionHtml(gammaTurns, 'gamma', '#1a6b8a', 'γ-turns', '')}
        </div>
        <div class="sub-tab-panel" id="sub-tab-panel-beta">
          ${turnSectionHtml(betaTurns, 'beta', '#7c3aed', 'β-turns', '')}
        </div>
        <div class="sub-tab-panel" id="sub-tab-panel-alpha">
          ${alphaTurnSectionHtml(alphaTurns)}
        </div>
        <div class="sub-tab-panel" id="sub-tab-panel-omega">
          ${omegaListHtml(omegaLoops)}
        </div>
      </div>`;

    const leftHtml = tabBarHtml + seqMotifPanel + structMotifPanel;

    if(hasViewer && structIds.length){
      const rightHtml=viewerPanelHtml(structIds, activeId);
      html+=`<div class="result-section"><div class="results-grid"><div>${leftHtml}</div>${rightHtml}</div></div>`;
    } else {
      html+=`<div class="result-section">${leftHtml}</div>`;
    }
  }

  // ── sequence section ───────────────────────────────────────────────────
  if(data.seq_files > 0){
    let seqHtmlStr=`<div class="section-banner seq"><span class="pill">sequences</span></div>`;
    seqHtmlStr+=statsHtml(seqMotifs, seqCounts, 'sequences', '--accent-seq');
    seqHtmlStr+=`<div class="section-h">motifs · list</div>`;
    seqHtmlStr+=motifListHtml(seqMotifs, 'seq-row', false);
    html+=`<div class="result-section">${seqHtmlStr}</div>`;
  }

  if(!html) html=`<div class="empty">No files were processed.</div>`;
  out.innerHTML=html;

  // ── top-level tab switching ───────────────────────────────────────────
  const tabBar = document.getElementById('struct-tab-bar');
  if(tabBar){
    tabBar.addEventListener('click', e=>{
      const btn = e.target.closest('.tab-btn');
      if(!btn || btn.classList.contains('sub-tab-btn')) return;
      const target = btn.dataset.tab;
      tabBar.querySelectorAll('.tab-btn:not(.sub-tab-btn)').forEach(b=>b.classList.remove('active'));
      out.querySelectorAll('.tab-panel').forEach(p=>p.classList.remove('active'));
      btn.classList.add('active');
      const panel = document.getElementById(`tab-panel-${target}`);
      if(panel) panel.classList.add('active');
      // viewer sync
      updateViewerLegend();
      if(target==='seq-motifs'){
        const sel2=$('#struct-select');
        const activeStruct = sel2 ? sel2.value : structMotifs[0]&&structMotifs[0].pdb_id;
        if(activeStruct) setTimeout(()=>loadForStructure(activeStruct, null), 60);
      } else if(target==='struct-motifs'){
        // activate first sub-tab viewer sync
        const activeSubBtn = document.querySelector('#sub-tab-bar .sub-tab-btn.active');
        const subTarget = activeSubBtn ? activeSubBtn.dataset.subtab : 'gamma';
        _syncSubTabViewer(subTarget);
      }
    });
  }

  
function updateViewerLegend(){
  const legendEl = document.querySelector('.viewer-legend');
  if(!legendEl) return;

  const structuralTabActive =
    document.querySelector('.tab-btn[data-tab="struct-motifs"].active') !== null;

  if(structuralTabActive){
    const items=[];
    if(typeof betaTurns!=='undefined' && betaTurns.length)
      items.push(`<span class="legend-pill" style="border-color:#7c3aed22;color:#7c3aed"><span class="legend-dot" style="background:#7c3aed"></span>β-turn</span>`);
    if(typeof gammaTurns!=='undefined' && gammaTurns.length)
      items.push(`<span class="legend-pill" style="border-color:#1a6b8a22;color:#1a6b8a"><span class="legend-dot" style="background:#1a6b8a"></span>γ-turn</span>`);
    if(typeof alphaTurns!=='undefined' && alphaTurns.length)
      items.push(`<span class="legend-pill" style="border-color:#b71c1c22;color:#b71c1c"><span class="legend-dot" style="background:#b71c1c"></span>α-turn</span>`);
    if(typeof omegaLoops!=='undefined' && omegaLoops.length)
      items.push(`<span class="legend-pill" style="border-color:#e05c0022;color:#e05c00"><span class="legend-dot" style="background:#e05c00"></span>Ω-loop</span>`);
    legendEl.innerHTML = items.join('');
  } else {
    const presentTypes=[...new Set(structMotifs.map(m=>m.seq_type))];
    legendEl.innerHTML = presentTypes.map(type=>{
      const k=type[0];
      const hex=HEX_COLOR[k] || '#666';
      return `<span class="legend-pill" style="border-color:${hex}22;color:${hex}"><span class="legend-dot" style="background:${hex}"></span>${type}</span>`;
    }).join('');
  }
}

// ── sub-tab switching (inside Structural Motifs) ──────────────────────
  function _syncSubTabViewer(subTarget){
    const sel2=$('#struct-select');
    if(subTarget==='omega'){
      const activeStruct = sel2 ? sel2.value : (omegaLoops[0]||{}).pdb_id;
      if(activeStruct && fileContents[activeStruct]){
        const firstLoop = omegaLoops.find(l=>l.pdb_id===activeStruct);
        if(firstLoop) setTimeout(()=>loadForOmegaLoop(firstLoop), 60);
      }
    } else if(subTarget==='beta'||subTarget==='gamma'||subTarget==='alpha'){
      const pool = subTarget==='beta' ? betaTurns : subTarget==='gamma' ? gammaTurns : alphaTurns;
      const activeStruct = sel2 ? sel2.value : (pool[0]||{}).pdb_id;
      if(activeStruct && fileContents[activeStruct]){
        const first = pool.find(t=>t.pdb_id===activeStruct);
        if(first) setTimeout(()=>loadForTurn(first), 60);
      }
    }
  }

  const subTabBar = document.getElementById('sub-tab-bar');
  if(subTabBar){
    subTabBar.addEventListener('click', e=>{
      const btn = e.target.closest('.sub-tab-btn');
      if(!btn) return;
      const subTarget = btn.dataset.subtab;
      subTabBar.querySelectorAll('.sub-tab-btn').forEach(b=>b.classList.remove('active'));
      out.querySelectorAll('.sub-tab-panel').forEach(p=>p.classList.remove('active'));
      btn.classList.add('active');
      const panel = document.getElementById(`sub-tab-panel-${subTarget}`);
      if(panel) panel.classList.add('active');
      _syncSubTabViewer(subTarget);
    });
  }

  // wire viewer dropdown
  const sel=$('#struct-select');
  if(sel) sel.addEventListener('change',()=>loadForStructure(sel.value,null));

  // wire X***X motif row clicks
  out.querySelectorAll('li.struct-row').forEach(li=>{
    li.addEventListener('click',()=>{
      out.querySelectorAll('li.struct-row').forEach(r=>r.classList.remove('active'));
      li.classList.add('active');
      const m=structMotifs[+li.dataset.idx];
      const sel2=$('#struct-select'); if(sel2) sel2.value=m.pdb_id;
      loadForStructure(m.pdb_id, m);
    });
  });

  // wire omega loop row clicks
  out.querySelectorAll('li.omega-row').forEach(li=>{
    li.addEventListener('click',()=>{
      out.querySelectorAll('li.omega-row').forEach(r=>r.classList.remove('active'));
      li.classList.add('active');
      const loop = omegaLoops[+li.dataset.idx];
      const sel2=$('#struct-select'); if(sel2) sel2.value=loop.pdb_id;
      loadForOmegaLoop(loop);
    });
  });

  // wire β/γ/α turn row clicks
  out.querySelectorAll('li.turn-row').forEach(li=>{
    li.addEventListener('click',()=>{
      out.querySelectorAll('li.turn-row').forEach(r=>r.classList.remove('active'));
      li.classList.add('active');
      const idx  = +li.dataset.idx;
      const turn = li.dataset.type === 'beta'  ? betaTurns[idx]
                 : li.dataset.type === 'gamma' ? gammaTurns[idx]
                 : alphaTurns[idx];
      const sel2 = $('#struct-select'); if(sel2) sel2.value = turn.pdb_id;
      loadForTurn(turn);
    });
  });

  // initial viewer load
  if(data.struct_files>0){
    const structIds=[...new Set(structMotifs.map(m=>m.pdb_id))];
    const activeId=structIds[0];
    if(activeId && fileContents[activeId]){
      setTimeout(()=>loadForStructure(activeId,null), 60);
    }
  }

  // auto-start guided tour if this was a demo load
  if(_demoAfterRender){
    _demoAfterRender = false;
    setTimeout(startTour, 700);
  }
}

// ── scan ──────────────────────────────────────────────────────────────────
async function scan(){
  errEl.textContent=''; out.innerHTML=''; warnBox.innerHTML='';

  // hard block if any invalid files are still in the queue
  const invalid = queue.filter(f=>fileKind(f)===null);
  if(invalid.length){
    errEl.textContent = `Remove unsupported file${invalid.length>1?'s':''} before scanning: `+
                        invalid.map(f=>f.name).join(', ');
    return;
  }
  runBtn.disabled=true; runBtn.textContent='Scanning…';

  // Show "scanning started" notice if email is provided
  const _emailVal = (document.getElementById('email-input')||{}).value||'';
  const _emailStatusEl = document.getElementById('email-status');
  if(_emailVal.trim() && _emailStatusEl){
    _emailStatusEl.innerHTML =
      `<span class="scan-notice">
        Your files have been uploaded successfully.<br>
        Scanning has started.<br>
        You will receive the results at <b>${_emailVal.trim()}</b>.<br>
        You may safely close this browser.
      </span>`;
    _emailStatusEl.className = '';
    // Send job-submitted confirmation email
    fetch('__APP_ROOT__/send-confirmation', {
      method: 'POST', headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({email: _emailVal.trim()})
    }).catch(()=>{});
  }

  // Read structure file contents client-side for 3Dmol
  fileContents={};
  for(const f of queue){
    const ext=f.name.split('.').pop().toLowerCase();
    if(!['pdb','cif','mmcif'].includes(ext)) continue;
    const stem=f.name.replace(/\.[^.]+$/,'');
    const fmt=ext==='pdb'?'pdb':'cif';
    try{ fileContents[stem]={data:await readFileText(f),format:fmt}; }
    catch(e){ console.warn('Could not read',f.name,e); }
  }

  const fd=new FormData();
  queue.forEach(f=>fd.append('files',f));
  try{
    const r=await fetch('__APP_ROOT__/analyze',{method:'POST',body:fd});
    const data=await r.json();
    if(!r.ok||data.error){ errEl.textContent=data.error||('HTTP '+r.status); }
    else{
      render(data);
      // auto-download ZIP if checkbox is ticked
      const autoZip = document.getElementById('auto-zip-chk');
      if(autoZip && autoZip.checked) csvZipAll.click();
      // send CSV results by email if address entered
      const emailEl  = document.getElementById('email-input');
      const statusEl = document.getElementById('email-status');
      const email = emailEl ? emailEl.value.trim() : '';
      if(email){
        statusEl.textContent = 'Sending…';
        statusEl.className = '';
        try{
          const er  = await fetch('__APP_ROOT__/send-results', {
            method:'POST', headers:{'Content-Type':'application/json'},
            body: JSON.stringify({email, ...data})
          });
          const ed = await er.json();
          if(er.ok && ed.ok){
            statusEl.innerHTML = `<span class="scan-notice" style="border-color:#0e7c6b">
              ✓ Scan complete! Results have been sent to <b>${email}</b>.<br>
              Please check your inbox (and spam folder if not received).
            </span>`;
            statusEl.className = 'ok';
          } else {
            statusEl.textContent = `✗ ${ed.error||'Send failed'}`;
            statusEl.className = 'err';
          }
        } catch(e){
          statusEl.textContent = `✗ ${e.message}`;
          statusEl.className = 'err';
        }
      }
    }
  }catch(e){ errEl.textContent='Request failed: '+e.message; }
  finally{ runBtn.disabled=queue.length===0; runBtn.textContent='Scan'; }
}
runBtn.addEventListener('click',scan);

// ── CSV export ────────────────────────────────────────────────────────────
function ifaceToResidues(ifaceStr, rangeStr){
  if(!ifaceStr || !rangeStr) return '';
  try{
    const start = parseInt(rangeStr.split('-')[0], 10);
    const nums = [...ifaceStr].map((ch,i)=>ch==='I'?start+i:null).filter(n=>n!==null);
    return nums.length ? nums.join(', ') : 'none';
  }catch(e){ return ifaceStr; }
}
// cols: array of strings OR array of [key, header] pairs
function toCSV(rows, cols){
  const esc = v => `"${String(v==null?'':v).replace(/"/g,'""')}"`;
  const isTuple = cols.length && Array.isArray(cols[0]);
  const keys    = isTuple ? cols.map(c=>c[0]) : cols;
  const headers = isTuple ? cols.map(c=>c[1]) : cols;
  const dataRows = rows.map(r =>
    keys.map((k,i) => {
      let v = r[k] != null ? r[k] : '';
      if(k === 'iface') v = ifaceToResidues(v, r['range']||'');
      return esc(v);
    }).join(',')
  );
  return [headers.join(','), ...dataRows].join('\n');
}
function dl(name,text){
  const b=new Blob([text],{type:'text/csv'}),u=URL.createObjectURL(b),a=document.createElement('a');
  a.href=u; a.download=name; a.click(); URL.revokeObjectURL(u);
}
// [data_key, csv_header]
const MOTIF_COLS=[['pdb_id','Pdb_id'],['seq_type','Sequence_type'],['sn','Sequence'],
  ['range','Range'],['chain_id','Chain_id'],['ss','Secondary structure'],['iface','Interface residue']];
const COUNT_COLS=[['pdb_id','Pdb_id'],['seq_type','Sequence_type'],['count','Count']];
const OMEGA_COLS=[
  ['pdb_id','Pdb_id'],['loop_number','Loop_number'],['Chain','Chain'],
  ['Start_Position','Start_Position'],['End_Position','End_Position'],
  ['Amino_Acid_Seq','Sequence'],
  ['Sec_Structure_before','Secondary_Structure_before'],
  ['Sec_Structure','Secondary_Structure'],
  ['Sec_Structure_After','Secondary_Structure_After'],
  ['distance_x','Cα(first)-Cα(last) distance'],['radius','Radius'],['Area','Area'],['Volume','Volume'],
  ['Center_x','Center_x'],['Center_y','Center_y'],['Center_z','Center_z'],
  ['Diameter','Diameter'],['Amino_acid1','Amino_acid1'],['Amino_acid2','Amino_acid2']];
const BETA_COLS=[
  ['pdb_id','Pdb_id'],['turn_number','Turn_number'],['Chain','Chain'],
  ['Start_Position','Start_Position'],['End_Position','End_Position'],
  ['Amino_Acid_Seq','Sequence'],['subtype','Subtype'],
  ['phi_i1','Phi_i+1'],['psi_i1','Psi_i+1'],['phi_i2','Phi_i+2'],['psi_i2','Psi_i+2'],
  ['Sec_Structure','Secondary_Structure'],['distance_x','Cα(i)-Cα(i+3) distance'],
  ['total_asa','Total_asa'],['Amino_acid1','Amino_acid1'],['Amino_acid4','Amino_acid4']];
const ALPHA_COLS=[
  ['pdb_id','Pdb_id'],['turn_number','Turn_number'],['Chain','Chain'],
  ['Start_Position','Start_Position'],['End_Position','End_Position'],
  ['Amino_Acid_Seq','Sequence'],['cluster','Type'],
  ['phi_i1','Phi_i+1'],['psi_i1','Psi_i+1'],['phi_i2','Phi_i+2'],['psi_i2','Psi_i+2'],
  ['phi_i3','Phi_i+3'],['psi_i3','Psi_i+3'],
  ['hbond_NO_dist','N(i+4)-O(i) distance'],['hbond_HNO_angle','H-N(i+4)···O(i) angle'],
  ['Sec_Structure','Secondary_Structure'],['distance_x','Cα(i)-Cα(i+4) distance'],
  ['total_asa','Total_asa'],['Amino_acid1','Amino_acid1'],['Amino_acid5','Amino_acid5']];
// Gamma: split hbond_dist into N-O (classic/inverse) and H-O (NHN) columns
const _NHN = new Set(['classical-nhn','inverse-nhn']);
function prepGamma(rows){
  return rows.map(r=>({...r,
    '_no_dist': _NHN.has(r.subtype) ? '' : (r.hbond_dist??''),
    '_ho_dist': _NHN.has(r.subtype) ? (r.ho_dist??'') : '',
    '_nhn_angle': _NHN.has(r.subtype) ? (r.hbond_dist??'') : '',
  }));
}
const GAMMA_COLS=[
  ['pdb_id','Pdb_id'],['turn_number','Turn_number'],['Chain','Chain'],
  ['Start_Position','Start_Position'],['End_Position','End_Position'],
  ['Amino_Acid_Seq','Sequence'],['subtype','Subtype'],
  ['phi_mid','Phi_mid'],['psi_mid','Psi_mid'],
  ['_no_dist','O(i)-H(i+2)/N(i+2) distance'],['_ho_dist','H(i)-O(i) distance'],['_nhn_angle','N-H(i+2)···N(i+1) angle'],
  ['Sec_Structure','Secondary_Structure'],
  ['total_asa','Total_asa'],['Amino_acid1','Amino_acid1'],['Amino_acid3','Amino_acid3']];
csvSM.addEventListener('click',()=>dl('Sequence-Motifs.csv',toCSV(structMotifs,MOTIF_COLS)));
csvSC.addEventListener('click',()=>dl('Sequence-Motifs-Counts.csv',toCSV(structCounts,COUNT_COLS)));
csvOL.addEventListener('click',()=>dl('Omega-Loops.csv',  toCSV(omegaLoops,  OMEGA_COLS)));
csvBT.addEventListener('click',()=>dl('Beta-Turns.csv',   toCSV(betaTurns,   BETA_COLS)));
csvGT.addEventListener('click',()=>dl('Gamma-Turns.csv',  toCSV(prepGamma(gammaTurns), GAMMA_COLS)));
csvAT.addEventListener('click',()=>dl('Alpha-Turns.csv',  toCSV(alphaTurns,  ALPHA_COLS)));
csvQM.addEventListener('click',()=>dl('Seq-Motifs.csv',   toCSV(seqMotifs,   MOTIF_COLS)));
csvQC.addEventListener('click',()=>dl('Seq-Counts.csv',   toCSV(seqCounts,   COUNT_COLS)));

// ── Download all as ZIP ───────────────────────────────────────────────────
csvZipAll.addEventListener('click', async ()=>{
  const datasets = [
    ['Sequence-Motifs.csv',        structMotifs,         MOTIF_COLS],
    ['Sequence-Motifs-Counts.csv', structCounts,         COUNT_COLS],
    ['Omega-Loops.csv',            omegaLoops,           OMEGA_COLS],
    ['Gamma-Turns.csv',            prepGamma(gammaTurns),GAMMA_COLS],
    ['Beta-Turns.csv',             betaTurns,            BETA_COLS],
    ['Alpha-Turns.csv',            alphaTurns,           ALPHA_COLS],
    ['Seq-Motifs.csv',             seqMotifs,            MOTIF_COLS],
    ['Seq-Counts.csv',             seqCounts,            COUNT_COLS],
  ].filter(([,rows])=>rows && rows.length);

  // Use server-side ZIP endpoint
  const payload = {};
  datasets.forEach(([name, rows, cols])=>{
    payload[name] = {rows, cols};
  });
  csvZipAll.textContent = 'Zipping…';
  csvZipAll.disabled = true;
  try{
    const r = await fetch('__APP_ROOT__/download-zip', {
      method:'POST',
      headers:{'Content-Type':'application/json'},
      body: JSON.stringify(payload)
    });
    if(!r.ok){ alert('ZIP failed: HTTP '+r.status); return; }
    const blob = await r.blob();
    const u = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href=u; a.download='motif-scanner-results.zip'; a.click();
    URL.revokeObjectURL(u);
  } catch(e){ alert('ZIP failed: '+e.message); }
  finally{ csvZipAll.textContent='↓ Download all as ZIP'; csvZipAll.disabled=false; }
});
// ── fullscreen 3D explorer ─────────────────────────────────────────────────
let fsViewer = null;
let fsCurrentPdb = null;
let fsChainColors = {};
let fsCurrStyle = 'cartoon';
let fsAnnotTab  = 'motif';  // active annotation tab: motif | omega | beta | gamma | alpha
const FS_PALETTE = ['#4fc3f7','#81c784','#ffb74d','#f06292','#ba68c8','#4dd0e1','#ff8a65','#a5d6a7'];
const AA3TO1 = {ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',
                HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',
                THR:'T',TRP:'W',TYR:'Y',VAL:'V'};

function parsePdbChains(pdbData) {
  const chains = {};
  const seen = new Set();
  for (const line of pdbData.split('\n')) {
    if (!line.startsWith('ATOM')) continue;
    if (line.substring(12,16).trim() !== 'CA') continue;
    const chain = line[21] || 'A';
    const resi  = parseInt(line.substring(22,26).trim());
    const resn  = line.substring(17,20).trim();
    const key   = chain + ':' + resi;
    if (seen.has(key)) continue;
    seen.add(key);
    if (!chains[chain]) chains[chain] = [];
    chains[chain].push({resi, resn, aa: AA3TO1[resn] || 'X'});
  }
  return chains;
}

// ── motif overlay for fullscreen viewer ───────────────────────────────────

// Paint all X***X motifs for pdbId as coloured sticks on top of the
// current base style.  Called after every chain-colour / style change so
// motifs stay visible regardless of representation mode.
function fsPaintMotifs(pdbId) {
  if (!fsViewer || !pdbId) return;
  const motifs = structMotifs.filter(m => m.pdb_id === pdbId);
  motifs.forEach(m => {
    const [s, e] = m.range.split('-').map(Number);
    const hex = colorHex(m.seq_type);
    for (let resi = s; resi <= e; resi++) {
      fsViewer.addStyle({chain: m.chain_id, resi}, {
        stick:  {color: hex, radius: 0.20},
        sphere: {color: hex, radius: 0.28, opacity: 0.85},
      });
    }
  });
  // no render() here — caller handles it to avoid double-render
}

// ── annotation tab switcher ───────────────────────────────────────────────
function switchFsAnnotTab(tab) {
  fsAnnotTab = tab;
  document.querySelectorAll('.fs-annot-tab').forEach(b =>
    b.classList.toggle('active', b.dataset.tab === tab));
  document.querySelectorAll('.fs-annot-pane').forEach(p =>
    p.classList.toggle('active', p.dataset.pane === tab));
  document.querySelectorAll('.fs-annot-item').forEach(el => el.classList.remove('active'));
  // repaint viewer with new annotation type overlaid
  fsApplyChainColors();
}

// ── paint active annotation type onto viewer ──────────────────────────────
function fsPaintActiveAnnotations() {
  if (!fsViewer || !fsCurrentPdb) return;
  if (fsAnnotTab === 'motif')      fsPaintMotifs(fsCurrentPdb);
  else if (fsAnnotTab === 'omega') fsPaintLoops(fsCurrentPdb);
  else if (fsAnnotTab === 'beta')  fsPaintTurns(fsCurrentPdb, 'beta',  '#7c3aed');
  else if (fsAnnotTab === 'gamma') fsPaintTurns(fsCurrentPdb, 'gamma', '#1a6b8a');
  else if (fsAnnotTab === 'alpha') fsPaintTurns(fsCurrentPdb, 'alpha', '#b71c1c');
}

function fsPaintLoops(pdbId) {
  if (!fsViewer || !pdbId) return;
  const loops = (typeof omegaLoops !== 'undefined' ? omegaLoops : [])
    .filter(l => l.pdb_id === pdbId);
  loops.forEach(l => {
    for (let r = l.Start_Position; r <= l.End_Position; r++) {
      fsViewer.addStyle({chain: l.Chain, resi: r}, {
        stick:  {color: '#e05c00', radius: 0.20},
        sphere: {color: '#e05c00', radius: 0.28, opacity: 0.85},
      });
    }
  });
}

function fsPaintTurns(pdbId, type, color) {
  if (!fsViewer || !pdbId) return;
  const pool = type === 'beta' ? (typeof betaTurns  !== 'undefined' ? betaTurns  : [])
             : type === 'gamma'? (typeof gammaTurns !== 'undefined' ? gammaTurns : [])
             :                   (typeof alphaTurns !== 'undefined' ? alphaTurns : []);
  pool.filter(t => t.pdb_id === pdbId).forEach(t => {
    for (let r = t.Start_Position; r <= t.End_Position; r++) {
      fsViewer.addStyle({chain: t.Chain, resi: r}, {
        stick:  {color, radius: 0.20},
        sphere: {color, radius: 0.28, opacity: 0.85},
      });
    }
  });
}

// ── populate annotation sidebar lists ────────────────────────────────────
function renderFsAnnotations(pdbId) {
  renderFsMotifList(pdbId);
  renderFsLoopList(pdbId);
  renderFsTurnList(pdbId, 'beta');
  renderFsTurnList(pdbId, 'gamma');
  renderFsTurnList(pdbId, 'alpha');
}

// X***X motif list
function renderFsMotifList(pdbId) {
  const el = document.getElementById('fs-annot-pane-motif');
  if (!el) return;
  const motifs = structMotifs.filter(m => m.pdb_id === pdbId);
  if (!motifs.length) {
    el.innerHTML = '<span class="fs-hint">No X\u2217\u2217\u2217X motifs for this structure</span>';
    return;
  }
  const MAX_SHOW = 60;
  const shown = motifs.slice(0, MAX_SHOW);
  const note  = motifs.length > MAX_SHOW
    ? `<div style="font-family:var(--mono);font-size:10px;color:#8b949e;padding:4px 2px">Showing ${MAX_SHOW} of ${motifs.length}</div>` : '';
  el.innerHTML = `<div class="fs-annot-list">` +
    shown.map((m, i) => {
      const hex = colorHex(m.seq_type);
      return `<div class="fs-annot-item" data-atype="motif" data-aidx="${i}"
        onclick="fsFocusMotif(${i})"
        title="${m.seq_type} · chain ${m.chain_id} · pos ${m.range}">
        <span class="fs-annot-badge" style="color:${hex}">${m.seq_type}</span>
        <span class="fs-annot-seq">${m.sn||''}</span>
        <span class="fs-annot-loc">ch.${m.chain_id}&nbsp;${m.range}</span>
      </div>`;
    }).join('') + `</div>${note}`;
}

// Ω-loop list
function renderFsLoopList(pdbId) {
  const el = document.getElementById('fs-annot-pane-omega');
  if (!el) return;
  const loops = (typeof omegaLoops !== 'undefined' ? omegaLoops : [])
    .filter(l => l.pdb_id === pdbId);
  if (!loops.length) {
    el.innerHTML = '<span class="fs-hint">No \u03a9-loops for this structure</span>';
    return;
  }
  el.innerHTML = `<div class="fs-annot-list">` +
    loops.map((l, i) =>
      `<div class="fs-annot-item" data-atype="omega" data-aidx="${i}"
        onclick="fsFocusAnnot('omega',${i})"
        title="\u03a9-loop · chain ${l.Chain} · ${l.Start_Position}\u2013${l.End_Position}">
        <span class="fs-annot-badge" style="color:#ff8c42">\u03a9</span>
        <span class="fs-annot-seq">${(l.Amino_Acid_Seq||'').substring(0,8)}</span>
        <span class="fs-annot-loc">ch.${l.Chain}&nbsp;${l.Start_Position}-${l.End_Position}</span>
      </div>`
    ).join('') + `</div>`;
}

// β / γ / α turn list
function renderFsTurnList(pdbId, type) {
  const colors  = {beta:'#a78bfa', gamma:'#4dd0e1', alpha:'#ff6f6f'};
  const labels  = {beta:'\u03b2', gamma:'\u03b3', alpha:'\u03b1'};
  const el = document.getElementById(`fs-annot-pane-${type}`);
  if (!el) return;
  const pool = type === 'beta'  ? (typeof betaTurns  !== 'undefined' ? betaTurns  : [])
             : type === 'gamma' ? (typeof gammaTurns !== 'undefined' ? gammaTurns : [])
             :                    (typeof alphaTurns !== 'undefined' ? alphaTurns : []);
  const turns = pool.filter(t => t.pdb_id === pdbId);
  if (!turns.length) {
    el.innerHTML = `<span class="fs-hint">No ${labels[type]}-turns for this structure</span>`;
    return;
  }
  const color = colors[type];
  const MAX_SHOW = 60;
  const shown = turns.slice(0, MAX_SHOW);
  const note  = turns.length > MAX_SHOW
    ? `<div style="font-family:var(--mono);font-size:10px;color:#8b949e;padding:4px 2px">Showing ${MAX_SHOW} of ${turns.length}</div>` : '';
  el.innerHTML = `<div class="fs-annot-list">` +
    shown.map((t, i) => {
      const badge = t.subtype || t.cluster || labels[type];
      return `<div class="fs-annot-item" data-atype="${type}" data-aidx="${i}"
        onclick="fsFocusAnnot('${type}',${i})"
        title="${labels[type]}-turn · chain ${t.Chain} · ${t.Start_Position}\u2013${t.End_Position}">
        <span class="fs-annot-badge" style="color:${color}" title="${badge}">${badge.substring(0,6)}</span>
        <span class="fs-annot-seq">${t.Amino_Acid_Seq||''}</span>
        <span class="fs-annot-loc">ch.${t.Chain}&nbsp;${t.Start_Position}-${t.End_Position}</span>
      </div>`;
    }).join('') + `</div>${note}`;
}

// ── focus a loop or turn in the viewer ───────────────────────────────────
function fsFocusAnnot(type, idx) {
  if (!fsViewer) return;
  const colorMap = {omega:'#e05c00', beta:'#7c3aed', gamma:'#1a6b8a', alpha:'#b71c1c'};
  const color = colorMap[type];

  let item, chain, start, end;
  if (type === 'omega') {
    const pool = (typeof omegaLoops !== 'undefined' ? omegaLoops : [])
      .filter(l => l.pdb_id === fsCurrentPdb);
    item = pool[idx]; if (!item) return;
    chain = item.Chain; start = item.Start_Position; end = item.End_Position;
  } else {
    const raw = type === 'beta'  ? (typeof betaTurns  !== 'undefined' ? betaTurns  : [])
              : type === 'gamma' ? (typeof gammaTurns !== 'undefined' ? gammaTurns : [])
              :                    (typeof alphaTurns !== 'undefined' ? alphaTurns : []);
    const pool = raw.filter(t => t.pdb_id === fsCurrentPdb);
    item = pool[idx]; if (!item) return;
    chain = item.Chain; start = item.Start_Position; end = item.End_Position;
  }

  // Reset base style with gray fallback
  fsViewer.setStyle({}, {});
  const cids = Object.keys(fsChainColors);
  if (cids.length) {
    fsViewer.setStyle({}, {[fsCurrStyle]: {color:'#666666'}});
    cids.forEach(c => fsViewer.setStyle({chain:c}, {[fsCurrStyle]: {color: fsChainColors[c]}}));
  } else {
    fsViewer.setStyle({}, {cartoon:{colorscheme:'chain'}});
  }
  fsPaintActiveAnnotations();

  // Brighten the focused item
  const resiArr = [];
  for (let r = start; r <= end; r++) resiArr.push(r);
  resiArr.forEach(resi => {
    fsViewer.addStyle({chain, resi}, {
      stick:  {color, radius: 0.34},
      sphere: {color, radius: 0.46, opacity: 1.0},
    });
  });

  fsViewer.zoomTo({chain, resi: resiArr});
  fsViewer.render();

  // Highlight sidebar row
  document.querySelectorAll(`.fs-annot-item[data-atype="${type}"]`).forEach(el =>
    el.classList.toggle('active', +el.dataset.aidx === idx)
  );
}

// Zoom to a motif and briefly brighten it; highlight the sidebar row.
function fsFocusMotif(idx) {
  if (!fsViewer) return;
  const motifs = structMotifs.filter(m => m.pdb_id === fsCurrentPdb);
  const m = motifs[idx];
  if (!m) return;
  const [s, e] = m.range.split('-').map(Number);
  const hex    = colorHex(m.seq_type);

  // Reset base style, repaint all annotations at normal weight
  fsViewer.setStyle({}, {});
  const cids = Object.keys(fsChainColors);
  if (cids.length) {
    fsViewer.setStyle({}, {[fsCurrStyle]: {color:'#666666'}});
    cids.forEach(c => fsViewer.setStyle({chain:c}, {[fsCurrStyle]: {color: fsChainColors[c]}}));
  } else {
    fsViewer.setStyle({}, {cartoon:{colorscheme:'chain'}});
  }
  fsPaintActiveAnnotations();

  // Brighten the focused motif on top of the existing overlay
  const resiArr = Array.from({length: e - s + 1}, (_, k) => s + k);
  resiArr.forEach(resi => {
    fsViewer.addStyle({chain: m.chain_id, resi}, {
      stick:  {color: hex, radius: 0.34},
      sphere: {color: hex, radius: 0.46, opacity: 1.0},
    });
  });

  fsViewer.zoomTo({chain: m.chain_id, resi: resiArr});
  fsViewer.render();

  // Highlight sidebar row
  document.querySelectorAll('.fs-annot-item[data-atype="motif"]').forEach(el =>
    el.classList.toggle('active', +el.dataset.aidx === idx)
  );
}

function openFsViewer(pdbId) {
  fsCurrentPdb = pdbId || Object.keys(fileContents)[0];
  if (!fsCurrentPdb || !fileContents[fsCurrentPdb]) return;
  const fd = fileContents[fsCurrentPdb];

  document.getElementById('fs-overlay').style.display = 'flex';
  document.getElementById('theme-toggle').style.display = 'none';
  { const hl = document.querySelector('.header-logo'); if (hl) hl.style.display = 'none'; }
  if(window.bgPause) window.bgPause();   // free GPU for fs viewer
  document.getElementById('fs-pdb-title').textContent = fsCurrentPdb;
  document.getElementById('fs-info').innerHTML = '<span class="fs-hint">Click any atom to view details</span>';

  // Parse chains from PDB data
  const chains = parsePdbChains(fd.data);
  let chainIds = Object.keys(chains);

  // Initialise viewer
  const fsContainer = document.getElementById('fs-mol-container');
  if (fsViewer) { try { fsViewer.clear(); } catch(e) {} }
  fsViewer = $3Dmol.createViewer(fsContainer, {backgroundColor:'#0d1117', antialias:true});
  fsViewer.addModel(fd.data, fd.format);

  // If no chains parsed (mmCIF), extract from model
  if (!chainIds.length) {
    const atoms = fsViewer.getModel().selectedAtoms({});
    chainIds = [...new Set(atoms.map(a => a.chain).filter(Boolean))];
    chainIds.forEach(c => { chains[c] = []; });
  }

  // Assign palette colors (persist across opens)
  chainIds.forEach((c,i) => { if(!fsChainColors[c]) fsChainColors[c] = FS_PALETTE[i%FS_PALETTE.length]; });

  // Apply chain colours & reset style
  fsCurrStyle = 'cartoon';
  fsAnnotTab  = 'motif';
  document.querySelectorAll('.fs-repr-btn').forEach(b=>b.classList.remove('active'));
  document.querySelector('.fs-repr-btn').classList.add('active');
  // Reset annotation tabs to X***X
  document.querySelectorAll('.fs-annot-tab').forEach(b=>b.classList.toggle('active',b.dataset.tab==='motif'));
  document.querySelectorAll('.fs-annot-pane').forEach(p=>p.classList.toggle('active',p.dataset.pane==='motif'));
  fsApplyChainColors();

  // Clickable residues — show info on click
  fsViewer.setClickable({}, true, atom => fsOnAtomClick(atom, chains));

  // Hover labels
  fsViewer.setHoverable({}, true,
    atom => {
      fsViewer.addLabel(`${atom.resn} ${atom.resi} · Ch.${atom.chain}`, {
        position:atom, fontSize:12, fontColor:'#e6edf3',
        backgroundColor:'#161b22', backgroundOpacity:0.92,
        borderColor:'#388bfd', borderThickness:1, borderOpacity:1
      });
      fsViewer.render();
    },
    () => { fsViewer.removeAllLabels(); fsViewer.render(); }
  );

  fsViewer.zoomTo(); fsViewer.render();

  // Render control panels
  renderFsChains(chainIds, chains);
  renderFsSeq(chainIds, chains);
  renderFsAnnotations(fsCurrentPdb);
}

function closeFsViewer() {
  document.getElementById('fs-overlay').style.display = 'none';
  document.getElementById('theme-toggle').style.display = '';
  { const hl = document.querySelector('.header-logo'); if (hl) hl.style.display = ''; }
  if(window.bgResume) window.bgResume();  // resume background rotation
  if (fsViewer) { try { fsViewer.clear(); } catch(e) {} fsViewer = null; }
}

function fsApplyChainColors() {
  if (!fsViewer) return;
  const cids = Object.keys(fsChainColors);
  if (cids.length) {
    fsViewer.setStyle({}, {[fsCurrStyle]: {color:'#666666'}});
    cids.forEach(c => fsViewer.setStyle({chain:c}, {[fsCurrStyle]:{color:fsChainColors[c]}}));
  } else {
    fsViewer.setStyle({}, {cartoon:{colorscheme:'chain'}});
  }
  fsPaintActiveAnnotations();
  fsViewer.render();
}

function fsOnAtomClick(atom, chains) {
  if (!fsViewer) return;
  // Dim all chains, highlight the clicked residue
  Object.keys(fsChainColors).forEach(c =>
    fsViewer.setStyle({chain:c}, {[fsCurrStyle]:{color:fsChainColors[c], opacity:0.2}})
  );
  fsViewer.setStyle({chain:atom.chain, resi:atom.resi}, {
    cartoon:{color:'#79c0ff', opacity:1},
    stick:{colorscheme:'Jmol', radius:0.22}
  });
  fsViewer.render();

  // Info card
  const aa1 = AA3TO1[atom.resn] || '?';
  document.getElementById('fs-info').innerHTML = `
    <div class="fs-info-card">
      <div class="fs-info-aa">${aa1}</div>
      <div class="fs-info-details">
        <div><b>${atom.resn}</b> · residue <b>${atom.resi}</b></div>
        <div>Chain <b>${atom.chain}</b></div>
        <div class="fs-info-atom">atom: ${atom.atom}</div>
      </div>
    </div>`;

  // Highlight residue in sequence strip
  document.querySelectorAll('.fs-res.active').forEach(el=>el.classList.remove('active'));
  const resEl = document.querySelector(`.fs-res[data-chain="${atom.chain}"][data-resi="${atom.resi}"]`);
  if (resEl) { resEl.classList.add('active'); resEl.scrollIntoView({behavior:'smooth',block:'nearest',inline:'center'}); }
}

function renderFsChains(chainIds, chains) {
  const html = chainIds.map(c => {
    const count = (chains[c]||[]).length;
    return `<div class="fs-chain-row">
      <input type="color" class="fs-color-pick" value="${fsChainColors[c]}"
        oninput="fsUpdateColor('${c}',this.value)">
      <span class="fs-chain-badge" id="fs-badge-${c}" style="background:${fsChainColors[c]}">Ch. ${c}</span>
      <span class="fs-chain-count">${count?count+' res':''}</span>
    </div>`;
  }).join('');
  document.getElementById('fs-chains').innerHTML = html || '<span class="fs-hint">No chains detected</span>';
}

function renderFsSeq(chainIds, chains) {
  const html = chainIds.map(c => {
    const residues = chains[c] || [];
    if (!residues.length) return '';
    const resHtml = residues.map(r =>
      `<span class="fs-res" data-chain="${c}" data-resi="${r.resi}"
        onclick="fsClickSeqRes('${c}',${r.resi})"
        title="${r.resn} ${r.resi} · Ch.${c}">${r.aa}</span>`
    ).join('');
    return `<div class="fs-chain-seq">
      <div class="fs-seq-label">Chain ${c}</div>
      <div class="fs-seq-strip" id="fs-seq-${c}">${resHtml}</div>
    </div>`;
  }).join('');
  document.getElementById('fs-seq-area').innerHTML =
    html || '<span class="fs-hint">Sequence available for .pdb files</span>';
}

function fsUpdateColor(chain, color) {
  fsChainColors[chain] = color;
  const badge = document.getElementById(`fs-badge-${chain}`);
  if (badge) badge.style.background = color;
  if (!fsViewer) return;
  fsViewer.setStyle({chain}, {[fsCurrStyle]:{color}});
  fsViewer.render();
}

function fsSetStyle(style, btn) {
  if (!fsViewer) return;
  fsCurrStyle = style;
  document.querySelectorAll('.fs-repr-btn').forEach(b=>b.classList.remove('active'));
  btn.classList.add('active');
  const cids = Object.keys(fsChainColors);
  if (cids.length) {
    fsViewer.setStyle({}, {[style]: {color:'#666666'}});
    cids.forEach(c => fsViewer.setStyle({chain:c},{[style]:{color:fsChainColors[c]}}));
  } else {
    fsViewer.setStyle({},{[style]:{colorscheme:'Jmol'}});
  }
  fsPaintActiveAnnotations();
  fsViewer.render();
}

function fsClickSeqRes(chain, resi) {
  if (!fsViewer) return;
  Object.keys(fsChainColors).forEach(c =>
    fsViewer.setStyle({chain:c},{[fsCurrStyle]:{color:fsChainColors[c],opacity:0.2}})
  );
  fsViewer.setStyle({chain,resi},{cartoon:{color:'#79c0ff',opacity:1},stick:{colorscheme:'Jmol',radius:0.22}});
  fsViewer.zoomTo({chain,resi}); fsViewer.render();
  document.querySelectorAll('.fs-res.active').forEach(el=>el.classList.remove('active'));
  const el = document.querySelector(`.fs-res[data-chain="${chain}"][data-resi="${resi}"]`);
  if (el) { el.classList.add('active'); el.scrollIntoView({behavior:'smooth',block:'nearest',inline:'center'}); }
}

function fsResetView() {
  if (!fsViewer) return;
  fsViewer.setStyle({},{});
  fsApplyChainColors();
  document.querySelectorAll('.fs-res.active').forEach(el=>el.classList.remove('active'));
  document.getElementById('fs-info').innerHTML = '<span class="fs-hint">Click any atom to view details</span>';
  fsViewer.removeAllLabels(); fsViewer.zoomTo(); fsViewer.render();
}

function openDefs(){
  document.getElementById('defs-overlay').style.display='flex';
  if(window.bgPause) window.bgPause();
}
function closeDefs(){
  document.getElementById('defs-overlay').style.display='none';
  if(window.bgResume) window.bgResume();
}

// ── demo loader ────────────────────────────────────────────────────────────
let _demoAfterRender = false;

async function loadDemo(){
  errEl.textContent = '';
  try {
    const resp = await fetch('__APP_ROOT__/demo');
    if(!resp.ok) throw new Error('demo.pdb not found — place it alongside app.py');
    const blob = await resp.blob();
    const file = new File([blob], '7om2.pdb', {type: 'text/plain'});
    queue.length = 0;
    addFiles([file]);
    _demoAfterRender = true;
    setTimeout(()=>runBtn.click(), 180);
  } catch(err) {
    errEl.textContent = '⚠ ' + err.message;
  }
}

// ── guided tour ────────────────────────────────────────────────────────────
const TOUR_STEPS = [
  {
    title: 'Welcome to the Motif Scanner',
    body:  `You loaded <b>7OM2</b> — the <b>Thosea Asigna Virus RNA-dependent RNA
            Polymerase (RdRP)</b>, solved by X-ray diffraction at 2.35 Å resolution.
            This tour walks through each results panel. Use next / prev to navigate,
            or close at any time.`,
    target: null, tab: null,
  },
  {
    title: 'Sequence Motifs',
    body:  `Each row is an <b>X***X motif</b> — a 5-residue window where the first and
            last residues are the same type from <code>A G E I L K V</code>.
            The three middle residues can be anything. Click any row to see DSSP
            secondary structure, chain, and solvent-accessibility data.`,
    target: '#tab-panel-seq-motifs', tab: 'seq-motifs',
  },
  {
    title: 'Structural Motifs — Ω-Loops',
    body:  `<b>Ω-loops</b> are 6-residue loops with ≥ 5 DSSP S/T codes (bend/turn)
            and closed geometry: the endpoint Cα atoms are within 10 Å and the loop is
            compact (endpoint dist &lt; ⅔ max pairwise dist). Click any row to highlight it
            in the 3D viewer.`,
    target: '#tab-panel-omega-loops', tab: 'omega-loops',
  },
  {
    title: 'α / β / γ Turns',
    body:  `Turns are chain-reversal motifs connects regular secondary structure elements.
            <b>β-turns</b> (4-res) are the most common, classified into 10 subtypes by
            Ramachandran angles. <b>α-turns</b> (5-res) have the same i→i+4 H-bond as
            one helix unit. <b>γ-turns</b> (3-res) are the rarest.
            <b>Click any row</b> to focus it in the 3D viewer.`,
    target: '#tab-panel-turns', tab: 'turns',
  },
  {
    title: '3D Viewer & Fullscreen Explorer',
    body:  `The viewer overlays <b>X***X motifs as coloured sticks</b> on the chain
            cartoon. Click <b>↗</b> (top-right of the viewer panel) to open the
            fullscreen explorer: chain colour controls, atom-click details, a full
            sequence strip, and a Motifs sidebar where clicking any motif zooms to it.`,
    target: '#mol-container', tab: null,
  },
];
let _tourIdx = 0;

function startTour(){
  _tourIdx = 0;
  _renderTourStep();
  document.getElementById('tour-card').style.display = 'block';
}

function _renderTourStep(){
  const step  = TOUR_STEPS[_tourIdx];
  const total = TOUR_STEPS.length;
  document.getElementById('tour-title').textContent    = step.title;
  document.getElementById('tour-body').innerHTML       = step.body;
  document.getElementById('tour-step-ind').textContent = `${_tourIdx + 1} / ${total}`;
  document.getElementById('tour-prev').disabled = _tourIdx === 0;
  document.getElementById('tour-next').textContent = _tourIdx === total - 1 ? 'finish ✓' : 'next →';

  // progress dots
  document.getElementById('tour-dots').innerHTML = TOUR_STEPS.map((_,i) =>
    `<span class="tour-dot${i===_tourIdx?' active':''}"></span>`
  ).join('');

  // remove previous highlight
  document.querySelectorAll('.tour-highlight').forEach(el => el.classList.remove('tour-highlight'));

  // switch tab if needed
  if(step.tab){
    const btn = document.querySelector(`[data-tab="${step.tab}"]`);
    if(btn && !btn.classList.contains('active')) btn.click();
  }

  // highlight and scroll
  if(step.target){
    const el = document.querySelector(step.target);
    if(el){
      el.classList.add('tour-highlight');
      setTimeout(() => el.scrollIntoView({behavior:'smooth', block:'nearest'}), 250);
    }
  } else {
    window.scrollTo({top: 0, behavior: 'smooth'});
  }
}

function moveTour(dir){
  const next = _tourIdx + dir;
  if(next >= TOUR_STEPS.length){ endTour(); return; }
  if(next < 0) return;
  _tourIdx = next;
  _renderTourStep();
}

function endTour(){
  document.getElementById('tour-card').style.display = 'none';
  document.querySelectorAll('.tour-highlight').forEach(el => el.classList.remove('tour-highlight'));
}

// ESC closes fullscreen and/or definitions overlay
document.addEventListener('keydown', e => {
  if(e.key==='Escape'){ closeFsViewer(); closeDefs(); endTour(); }
});

</script>
</body>
</html>"""


@app.route("/demo")
def serve_demo():
    """Serve the bundled demo PDB file for the tutorial."""
    demo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "demo.pdb")
    if not os.path.exists(demo_path):
        return (
            "demo.pdb not found. Place the demo PDB file alongside app.py and name it demo.pdb.",
            404,
        )
    return send_file(demo_path, mimetype="text/plain", as_attachment=False,
                     download_name="7om2_rdrp_demo.pdb")


@app.route("/")
def index():
    # PAGE is a single inline document whose links and fetch() calls are written
    # against __APP_ROOT__. Flask reports where the app is mounted as script_root:
    # "" when this file is run on its own, "/motif-scanner" when StructBio mounts it.
    page = PAGE.replace("__APP_ROOT__", request.script_root)
    return Response(page, mimetype="text/html")


@app.route("/analyze", methods=["POST"])
def analyze():
    uploaded = request.files.getlist("files")
    if not uploaded:
        return jsonify(error="No files received."), 400

    struct_motifs, seq_motifs = [], []
    omega_loops               = []
    beta_turns, gamma_turns, alpha_turns = [], [], []
    struct_files, seq_files   = 0, 0
    all_warnings = []
    sasa_statuses = []   # one per structure file

    with tempfile.TemporaryDirectory() as tmp:
        for fs in uploaded:
            name = secure_filename(fs.filename or "")
            ext  = os.path.splitext(name)[1].lower()
            if ext not in ALLOWED:
                return jsonify(error=f"'{name}': unsupported type (need .pdb/.cif/.mmcif or .txt/.fasta/.fa/.faa)."), 400

            path = os.path.join(tmp, name)
            fs.save(path)

            if ext in SEQ_EXT:
                # sequence file
                try:
                    text = open(path, encoding="utf-8", errors="replace").read()
                    source = os.path.splitext(name)[0]
                    m, _, warnings = find_motifs_in_sequence(text, source)
                    seq_motifs.extend(m)
                    all_warnings.extend(warnings)
                    seq_files += 1
                except Exception as e:
                    traceback.print_exc()
                    return jsonify(error=f"Error processing sequence file '{name}': {type(e).__name__}: {e}"), 400
            else:
                # structure file — run X***X motif scan
                stem = os.path.splitext(name)[0]
                try:
                    _result    = find_motifs(path)
                    m          = _result[0]
                    sasa_status = _result[2] if len(_result) > 2 else "unknown"
                    struct_motifs.extend(m)
                    struct_files += 1
                    if sasa_status == "single-chain":
                        sasa_statuses.append(
                            f"[FreeSASA] {stem}: single-chain — no interface to compute"
                        )
                    elif sasa_status.startswith("error"):
                        sasa_statuses.append(
                            f"[FreeSASA] {stem}: {sasa_status} — check server terminal for details"
                        )
                except UnsupportedFileError as e:
                    return jsonify(error=str(e)), 400
                except Exception as e:
                    traceback.print_exc()
                    return jsonify(error=f"Error processing '{name}': {type(e).__name__}: {e}"), 400

                # structure file — run Ω-loop scan (only .pdb supported by DSSP)
                if ext == ".pdb":
                    try:
                        loops, omega_status = find_omega_loops(path)
                        omega_loops.extend(loops)
                        if omega_status == "dssp-failed":
                            all_warnings.append(
                                f"[DSSP] {stem}: DSSP could not run — "
                                f"is the 'dssp' binary installed? (apt install dssp)"
                            )
                        elif omega_status.startswith("error"):
                            all_warnings.append(f"[Ω-loop] {stem}: {omega_status}")
                        # "ok:N" and "no-loops" need no warning
                    except Exception:
                        traceback.print_exc()
                        all_warnings.append(f"[Ω-loop] {stem}: unexpected error — check server terminal")

                    # β/γ-turn scan — reuses same DSSP + SASA via find_all_turns
                    try:
                        b, g, a, turn_status = find_all_turns(path)
                        beta_turns.extend(b)
                        gamma_turns.extend(g)
                        alpha_turns.extend(a)
                        if "dssp-failed" in turn_status:
                            all_warnings.append(
                                f"[β/γ-turn] {stem}: DSSP failed — turns require DSSP binary"
                            )
                        elif "error:" in turn_status:
                            all_warnings.append(f"[β/γ-turn] {stem}: {turn_status}")
                    except Exception:
                        traceback.print_exc()
                        all_warnings.append(f"[β/γ-turn] {stem}: unexpected error — check server terminal")

                else:
                    all_warnings.append(
                        f"[Ω-loop / β/γ-turn] {stem}: scan skipped for .cif/.mmcif "
                        f"(DSSP requires .pdb format)"
                    )

    return jsonify(
        struct_motifs  = struct_motifs,
        struct_counts  = summarise_counts(struct_motifs),
        struct_files   = struct_files,
        omega_loops    = omega_loops,
        beta_turns     = beta_turns,
        gamma_turns    = gamma_turns,
        alpha_turns    = alpha_turns,
        seq_motifs     = seq_motifs,
        seq_counts     = summarise_counts(seq_motifs),
        seq_files      = seq_files,
        warnings       = all_warnings + sasa_statuses,
    )


# Column defs: (data_key_or_callable, csv_header)
# callable receives the full row dict and returns the cell value

_NHN_SUBTYPES = {'classical-nhn', 'inverse-nhn'}

_MOTIF_COLS = [
    ('pdb_id',   'Pdb_id'),
    ('seq_type', 'Sequence_type'),
    ('sn',       'Sequence'),
    ('range',    'Range'),
    ('chain_id', 'Chain_id'),
    ('ss',       'Secondary structure'),
    ('iface',    'Interface residue'),
]
_COUNT_COLS = [
    ('pdb_id',   'Pdb_id'),
    ('seq_type', 'Sequence_type'),
    ('count',    'Count'),
]
_OMEGA_COLS = [
    ('pdb_id',               'Pdb_id'),
    ('loop_number',          'Loop_number'),
    ('Chain',                'Chain'),
    ('Start_Position',       'Start_Position'),
    ('End_Position',         'End_Position'),
    ('Amino_Acid_Seq',       'Sequence'),
    ('Sec_Structure_before', 'Secondary_Structure_before'),
    ('Sec_Structure',        'Secondary_Structure'),
    ('Sec_Structure_After',  'Secondary_Structure_After'),
    ('distance_x',           'Cα(first)-Cα(last) distance'),
    ('radius',               'Radius'),
    ('Area',                 'Area'),
    ('Volume',               'Volume'),
    ('Center_x',             'Center_x'),
    ('Center_y',             'Center_y'),
    ('Center_z',             'Center_z'),
    ('Diameter',             'Diameter'),
    ('Amino_acid1',          'Amino_acid1'),
    ('Amino_acid2',          'Amino_acid2'),
]
_BETA_COLS = [
    ('pdb_id',         'Pdb_id'),
    ('turn_number',    'Turn_number'),
    ('Chain',          'Chain'),
    ('Start_Position', 'Start_Position'),
    ('End_Position',   'End_Position'),
    ('Amino_Acid_Seq', 'Sequence'),
    ('subtype',        'Subtype'),
    ('phi_i1',         'Phi_i+1'),
    ('psi_i1',         'Psi_i+1'),
    ('phi_i2',         'Phi_i+2'),
    ('psi_i2',         'Psi_i+2'),
    ('Sec_Structure',  'Secondary_Structure'),
    ('distance_x',     'Cα(i)-Cα(i+3) distance'),
    ('total_asa',      'Total_asa'),
    ('Amino_acid1',    'Amino_acid1'),
    ('Amino_acid4',    'Amino_acid4'),
]
_GAMMA_COLS = [
    ('pdb_id',         'Pdb_id'),
    ('turn_number',    'Turn_number'),
    ('Chain',          'Chain'),
    ('Start_Position', 'Start_Position'),
    ('End_Position',   'End_Position'),
    ('Amino_Acid_Seq', 'Sequence'),
    ('subtype',        'Subtype'),
    ('phi_mid',        'Phi_mid'),
    ('psi_mid',        'Psi_mid'),
    # N-O distance for classic/inverse; H-O distance + N-H...N angle for NHN types
    (lambda r: r.get('hbond_dist', '') if r.get('subtype', '') not in _NHN_SUBTYPES else '', 'O(i)-H(i+2)/N(i+2) distance'),
    (lambda r: r.get('ho_dist', '') if r.get('subtype', '') in _NHN_SUBTYPES else '',        'H(i)-O(i) distance'),
    (lambda r: r.get('hbond_dist', '') if r.get('subtype', '') in _NHN_SUBTYPES else '',     'N-H(i+2)···N(i+1) angle'),
    ('Sec_Structure',  'Secondary_Structure'),
    ('total_asa',      'Total_asa'),
    ('Amino_acid1',    'Amino_acid1'),
    ('Amino_acid3',    'Amino_acid3'),
]
_ALPHA_COLS = [
    ('pdb_id',           'Pdb_id'),
    ('turn_number',      'Turn_number'),
    ('Chain',            'Chain'),
    ('Start_Position',   'Start_Position'),
    ('End_Position',     'End_Position'),
    ('Amino_Acid_Seq',   'Sequence'),
    ('cluster',          'Type'),
    ('phi_i1',           'Phi_i+1'),
    ('psi_i1',           'Psi_i+1'),
    ('phi_i2',           'Phi_i+2'),
    ('psi_i2',           'Psi_i+2'),
    ('phi_i3',           'Phi_i+3'),
    ('psi_i3',           'Psi_i+3'),
    ('hbond_NO_dist',    'N(i+4)-O(i) distance'),
    ('hbond_HNO_angle',  'H-N(i+4)···O(i) angle'),
    ('Sec_Structure',    'Secondary_Structure'),
    ('distance_x',       'Cα(i)-Cα(i+4) distance'),
    ('total_asa',        'Total_asa'),
    ('Amino_acid1',      'Amino_acid1'),
    ('Amino_acid5',      'Amino_acid5'),
]

# Lookup used by /download-zip to avoid sending col defs over the wire
_CSV_COLS_MAP = {
    'Sequence-Motifs.csv':        _MOTIF_COLS,
    'Sequence-Motifs-Counts.csv': _COUNT_COLS,
    'Omega-Loops.csv':            _OMEGA_COLS,
    'Beta-Turns.csv':             _BETA_COLS,
    'Gamma-Turns.csv':            _GAMMA_COLS,
    'Alpha-Turns.csv':            _ALPHA_COLS,
    'Seq-Motifs.csv':             _MOTIF_COLS,
    'Seq-Counts.csv':             _COUNT_COLS,
}


def _iface_to_residues(iface_str, range_str):
    """Convert '...II' + '10-14' → '12, 13' (residue numbers at interface)."""
    if not iface_str or not range_str:
        return ''
    try:
        start = int(range_str.split('-')[0])
        nums = [str(start + i) for i, ch in enumerate(iface_str) if ch == 'I']
        return ', '.join(nums) if nums else 'none'
    except Exception:
        return iface_str


def _make_csv(rows, cols):
    buf = io.StringIO()
    w = csv.writer(buf)
    # cols: list of (key_or_callable, header) tuples, or plain strings
    if cols and isinstance(cols[0], (tuple, list)):
        keys    = [c[0] for c in cols]
        headers = [c[1] for c in cols]
    else:
        keys = headers = list(cols)
    w.writerow(headers)
    for r in rows:
        row_vals = []
        for k in keys:
            if callable(k):
                val = k(r)
            else:
                val = r.get(k, '')
                if k == 'iface':
                    val = _iface_to_residues(val, r.get('range', ''))
            row_vals.append(val)
        w.writerow(row_vals)
    return buf.getvalue()

@app.route("/send-results", methods=["POST"])
def send_results():
    data = request.get_json(silent=True)
    if not data:
        return jsonify(error="No data received."), 400

    email_to = data.get("email", "").strip()
    if not email_to or "@" not in email_to or "." not in email_to.split("@")[-1]:
        return jsonify(error="Invalid email address."), 400

    smtp_user = os.environ.get("SMTP_USER", "").strip()
    smtp_pass = os.environ.get("SMTP_PASS", "").strip().replace(" ", "")
    smtp_host = os.environ.get("SMTP_HOST", "smtp.gmail.com")
    smtp_port = int(os.environ.get("SMTP_PORT", "587"))

    if not smtp_user or not smtp_pass:
        return jsonify(error="Email not configured on this server."), 503

    datasets = [
        ("Sequence-Motifs.csv",        data.get("struct_motifs", []), _MOTIF_COLS),
        ("Sequence-Motifs-Counts.csv", data.get("struct_counts", []), _COUNT_COLS),
        ("Omega-Loops.csv",    data.get("omega_loops",   []), _OMEGA_COLS),
        ("Gamma-Turns.csv",    data.get("gamma_turns",   []), _GAMMA_COLS),
        ("Beta-Turns.csv",     data.get("beta_turns",    []), _BETA_COLS),
        ("Alpha-Turns.csv",    data.get("alpha_turns",   []), _ALPHA_COLS),
        ("Seq-Motifs.csv",     data.get("seq_motifs",    []), _MOTIF_COLS),
        ("Seq-Counts.csv",     data.get("seq_counts",    []), _COUNT_COLS),
    ]
    attachments = [(fname, _make_csv(rows, cols)) for fname, rows, cols in datasets if rows]
    if not attachments:
        return jsonify(error="No results to send — run a scan first."), 400

    # Build email with CSV attachments
    msg = MIMEMultipart()
    msg["From"]    = f"Motif Scanner <{smtp_user}>"
    msg["To"]      = email_to
    msg["Subject"] = "Motif Scanner Results"
    file_list = "\n".join(f"  • {fname}" for fname, _ in attachments)
    msg.attach(MIMEText(
        f"Hello,\n\n"
        f"Your motif scan has completed successfully.\n\n"
        f"Please find the attached CSV file(s) ({len(attachments)} file(s)):\n"
        f"{file_list}\n\n"
        f"Thank you for using Motif Scanner.\n\n"
        f"— Motif Scanner Team\n",
        "plain"
    ))
    for fname, content in attachments:
        part = MIMEBase("application", "octet-stream")
        part.set_payload(content.encode("utf-8"))
        encoders.encode_base64(part)
        part.add_header("Content-Disposition", f'attachment; filename="{fname}"')
        msg.attach(part)

    try:
        with smtplib.SMTP(smtp_host, smtp_port) as s:
            s.ehlo()
            s.starttls()
            s.ehlo()
            s.login(smtp_user, smtp_pass)
            s.sendmail(smtp_user, email_to, msg.as_string())
    except smtplib.SMTPAuthenticationError:
        return jsonify(error=(
            "Gmail login failed. Use a Gmail App Password, not your regular password. "
            "Go to myaccount.google.com/apppasswords to create one."
        )), 500
    except Exception as e:
        traceback.print_exc()
        return jsonify(error=f"Failed to send email: {e}"), 500

    return jsonify(ok=True, message=f"Results sent to {email_to}")


@app.route("/send-confirmation", methods=["POST"])
def send_confirmation():
    data = request.get_json(silent=True)
    if not data:
        return jsonify(error="No data received."), 400

    email_to = data.get("email", "").strip()
    if not email_to or "@" not in email_to or "." not in email_to.split("@")[-1]:
        return jsonify(error="Invalid email address."), 400

    smtp_user = os.environ.get("SMTP_USER", "").strip()
    smtp_pass = os.environ.get("SMTP_PASS", "").strip().replace(" ", "")
    smtp_host = os.environ.get("SMTP_HOST", "smtp.gmail.com")
    smtp_port = int(os.environ.get("SMTP_PORT", "587"))

    if not smtp_user or not smtp_pass:
        return jsonify(error="Email not configured on this server."), 503

    msg = MIMEMultipart()
    msg["From"]    = f"Motif Scanner <{smtp_user}>"
    msg["To"]      = email_to
    msg["Subject"] = "Motif Scanner — Scan Started"
    msg.attach(MIMEText(
        f"Hello,\n\n"
        f"Your files have been uploaded successfully and your scan has started.\n\n"
        f"You will receive another email with your results once the scan is complete.\n\n"
        f"Thank you for using Motif Scanner.\n\n"
        f"— Motif Scanner Team\n",
        "plain"
    ))

    try:
        with smtplib.SMTP(smtp_host, smtp_port) as s:
            s.ehlo()
            s.starttls()
            s.ehlo()
            s.login(smtp_user, smtp_pass)
            s.sendmail(smtp_user, email_to, msg.as_string())
    except smtplib.SMTPAuthenticationError:
        return jsonify(error=(
            "Gmail login failed. Use a Gmail App Password, not your regular password."
        )), 500
    except Exception as e:
        traceback.print_exc()
        return jsonify(error=f"Failed to send confirmation email: {e}"), 500

    return jsonify(ok=True, message=f"Confirmation sent to {email_to}")


@app.route("/download-zip", methods=["POST"])
def download_zip():
    import zipfile
    data = request.get_json(silent=True) or {}
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zf:
        for fname, payload in data.items():
            rows = payload.get('rows', [])
            if not rows:
                continue
            cols = _CSV_COLS_MAP.get(fname, payload.get('cols', []))
            zf.writestr(fname, _make_csv(rows, cols))
    buf.seek(0)
    return send_file(buf, mimetype='application/zip',
                     as_attachment=True,
                     download_name='motif-scanner-results.zip')


@app.route("/logo")
def serve_logo():
    logo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "IITJ COLOURED.png")
    return send_file(logo_path, mimetype="image/png")


@app.route("/motif-logo")
def serve_motif_logo():
    logo_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logo", "motif1.png")
    return send_file(logo_path, mimetype="image/png")


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
