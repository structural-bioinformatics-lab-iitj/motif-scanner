"""
check_interface.py  —  one-off diagnostic for the interface detector.
    python check_interface.py 4i71_1.pdb
Tells you (a) how many interface residues were found total, per chain,
and (b) whether ANY X***X motif window actually sits on the interface.
"""
import sys, warnings
warnings.filterwarnings("ignore")
import motif_backend as mb

f = sys.argv[1]
s = mb.parse_structure(f)

iface, status = mb.get_interface_residues(f, s)
print("interface status :", status)
print("total interface residues:", len(iface))

from collections import Counter
print("per chain        :", dict(Counter(c for c, r in iface)))

# also report the raw chain inventory + whether the two chains touch at all
model = s[0]
chains = [c.id for c in model.get_chains()]
print("chains in file   :", chains)

motifs, _, sasa, ss = mb.find_motifs(f)
hits = [m for m in motifs if "I" in m["iface"]]
print(f"\nmotifs total: {len(motifs)}  |  motifs touching interface: {len(hits)}")
for m in hits:
    print(f"   {m['seq_type']}  chain {m['chain_id']}  {m['range']}  iface={m['iface']}")
