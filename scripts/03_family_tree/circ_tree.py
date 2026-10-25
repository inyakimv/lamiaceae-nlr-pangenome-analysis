import re, sys, glob
from collections import Counter
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch
from Bio import Phylo

TREE = sys.argv[1] if len(sys.argv) > 1 else "nbarc_family_tree_plus_outgroup.treefile"
OUT = sys.argv[2] if len(sys.argv) > 2 else "Fig_nbarc_circular"
ROOT = "tnl"         # "tnl" (root on the Arabidopsis TIR clade), "asis" (as written by FastTree), "midpoint", or a leaf name
TIR_GLOBS = ["../nlr_domain_search/**/*.domtblout", "../arabidopsis_outgroup/**/*.domtblout"]
TOP_N = 8            # number of Lamiaceae species coloured (the ones with most NB-ARC leaves)
SPECIES = []         # or give your own list, e.g. ["Salvia miltiorrhiza", "Mentha aquatica"]
PALETTE = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#CC79A7", "#7B3294", "#8C510A", "#01665E"]
ALGA, ARAB, OTHER = "#d62728", "#1a237e", "#c9c9c9"

tir = set()
for pat in TIR_GLOBS:
    for f in glob.glob(pat, recursive=True):
        for line in open(f):
            if line.startswith("#"): continue
            p = line.split()
            if len(p) > 12 and p[0].startswith("TIR") and float(p[12]) < 1e-5: tir.add(p[3])
def istir(l): return l.name.split("|", 1)[-1] in tir

def species(name):
    if name.startswith("Knitens"): return "Klebsormidium nitens (alga)"
    if name.startswith("Athaliana"): return "Arabidopsis thaliana"
    m = re.match(r"^(?:GC[AF]_\d+\.\d+_)?(.+?)_genomic\|", name)
    return m.group(1).replace("_x_", " x ").replace("_", " ") if m else name

tree = Phylo.read(TREE, "newick")
if ROOT == "midpoint": tree.root_at_midpoint()
elif ROOT == "tnl":
    tree.root_with_outgroup([t for t in tree.get_terminals() if t.name.startswith("Athaliana") and istir(t)])
elif ROOT != "asis":
    tree.root_with_outgroup(next(t for t in tree.get_terminals() if t.name == ROOT))
tree.ladderize()

leaves = tree.get_terminals()
N = len(leaves)
sp = {l: species(l.name) for l in leaves}
cnt = Counter(v for v in sp.values() if "alga" not in v and "Arabidopsis" not in v)
chosen = SPECIES if SPECIES else [s for s, _ in cnt.most_common(TOP_N)]
col = dict(zip(chosen, PALETTE))
def color(l):
    s = sp[l]
    if "alga" in s: return ALGA
    if "Arabidopsis" in s: return ARAB
    return col.get(s, OTHER)

depth = tree.depths()
ang, lo, hi, parent = {}, {}, {}, {}
for c in tree.find_clades(order="level"):
    for ch in c.clades: parent[ch] = c
for i, l in enumerate(leaves): ang[l] = 2 * np.pi * i / N
for c in tree.find_clades(order="postorder"):
    if c.is_terminal(): continue
    a = [ang[ch] for ch in c.clades]
    lo[c], hi[c] = min(a), max(a)
    ang[c] = (lo[c] + hi[c]) / 2

def xy(r, a): return (r * np.cos(a), r * np.sin(a))
segs, cols, lw = [], [], []
for c in tree.find_clades():
    if c in parent:
        p = parent[c]
        segs.append([xy(depth[p], ang[c]), xy(depth[c], ang[c])])
        if c.is_terminal() and color(c) != OTHER:
            cols.append(color(c)); lw.append(1.0 if color(c) == ALGA else 0.5)
        else:
            cols.append("#8a8a8a"); lw.append(0.25)
    if not c.is_terminal():
        n = max(2, int((hi[c] - lo[c]) / (2 * np.pi) * 700))
        aa = np.linspace(lo[c], hi[c], n)
        segs.append([xy(depth[c], a) for a in aa]); cols.append("#8a8a8a"); lw.append(0.25)

rmax = max(depth[l] for l in leaves)
R1, R2 = rmax * 1.03, rmax * 1.10
ring = [[xy(R1, ang[l]), xy(R2, ang[l])] for l in leaves]
ringc = [color(l) for l in leaves]
ring2 = [[xy(R2 * 1.01, ang[l]), xy(R2 * 1.045, ang[l])] for l in leaves if istir(l)]

fig, ax = plt.subplots(figsize=(11, 11))
ax.add_collection(LineCollection(segs, colors=cols, linewidths=lw))
ax.add_collection(LineCollection(ring, colors=ringc, linewidths=0.7))
ax.add_collection(LineCollection(ring2, colors="black", linewidths=0.7))
for l in leaves:
    if "alga" in sp[l]:
        x, y = xy(depth[l], ang[l]); ax.plot(x, y, "o", color=ALGA, ms=6, zorder=5)
        xt, yt = xy(R2 * 1.12, ang[l])
        ax.annotate(l.name.split("|")[1], (x, y), (xt, yt), fontsize=8, color=ALGA,
                    ha="center", va="center", arrowprops=dict(arrowstyle="-", color=ALGA, lw=0.6))
ax.set_xlim(-R2 * 1.25, R2 * 1.25); ax.set_ylim(-R2 * 1.25, R2 * 1.25)
ax.set_aspect("equal"); ax.axis("off")
handles = [Patch(color=col[s], label="%s (%d)" % (s, cnt[s])) for s in chosen]
handles += [Patch(color=OTHER, label="other Lamiaceae (%d)" % sum(v for s, v in cnt.items() if s not in chosen)),
            Patch(color=ARAB, label="Arabidopsis thaliana (%d)" % sum("Arabidopsis" in v for v in sp.values())),
            Patch(color=ALGA, label="Klebsormidium nitens, alga (%d)" % sum("alga" in v for v in sp.values())),
            Patch(color="black", label="outer black ring: TIR domain present (%d)" % sum(istir(l) for l in leaves))]
ax.legend(handles=handles, loc="lower left", fontsize=8, frameon=False, bbox_to_anchor=(-0.02, -0.02))
rt = {"asis": "as written by FastTree, arbitrary", "tnl": "branch of the Arabidopsis TIR (TNL) clade"}.get(ROOT, ROOT)
ax.set_title("NB-ARC family tree, %d domains (root: %s; display only)" % (N, rt), fontsize=11)
fig.savefig(OUT + ".png", dpi=220, bbox_inches="tight"); fig.savefig(OUT + ".pdf", bbox_inches="tight")
print("saved", OUT + ".png", "and .pdf;", N, "leaves; coloured species:", chosen)
