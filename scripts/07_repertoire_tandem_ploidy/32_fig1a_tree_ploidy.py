#!/usr/bin/env python3
"""Fig. 1a (ANNEVO version) with species dendrogram and ploidy colours.
Usage: python3 32_fig1a_tree_ploidy.py TANDEM_ANNEVO.csv [OUT_PREFIX=Fig1a_tree_bars_ploidy] [PLOIDY_CSV] [EXCLUDE=species,substrings]
  PLOIDY_CSV: table with columns species,ploidy (diploid/polyploid), e.g. figures_data/nlr_ploidy_binary_FINAL_source_data.csv.
  If omitted, the built-in list of 7 polyploids (from the previous figure) is used.
Input: output of 29_tandem_annevo.py (columns species, n_nbarc_genes). Writes <out>.png, <out>.pdf, <out>_source_data.csv
Tree: cladogram topology of the 46-species tree of the previous Fig. 1a (recovered from the figure, no branch lengths).
Ploidy classes below are those used in the previous Fig. 1a (7 polyploids): EDIT PLOIDY if your master table differs."""
import sys, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd
NEWICK = "(((Ocimum_basilicum,Ajuga_reptans),((Ajuga_chamaepitys,(Pogostemon_cablin,Schizonepeta_tenuifolia)),((Scutellaria_minor,Thymus_quinquecostatus),(Lycopus_lucidus,(Stachys_alpina,Clinopodium_barosmum))))),((Salvia_rosmarinus,((Teucrium_scorodonia,(Coleus_hadiensis,Lavandula_angustifolia)),(Scutellaria_barbata,((Agastache_rugosa,Salvia_sclarea),(Plectranthus_ornatus,Isodon_lophanthoides))))),(Galeopsis_tetrahit,(Scutellaria_baicalensis,(((Mentha_longifolia,(Drepanocaryum_sewerzowii,Mentha_suaveolens)),(Salvia_divinorum,(Salvia_miltiorrhiza,(Phlomoides_rotata,Salvia_hispanica)))),(Dracocephalum_rupestre,(Elsholtzia_splendens,(Glechoma_hederacea,(Thymus_mandschuricus,(((Prunella_vulgaris,Salvia_splendens),(Mentha_x_piperita,Mentha_aquatica)),(Agastache_foeniculum,((Lamium_purpureum,Marmoritis_complanata,Teucrium_botrys),(Melittis_melissophyllum,(Scutellaria_galericulata,(Lamium_amplexicaule,(Callicarpa_nudiflora,(Paraleonurus_japonicus,Ballota_nigra)))))))))))))))));"
POLYPLOID = {"Ocimum basilicum", "Pogostemon cablin", "Galeopsis tetrahit", "Glechoma hederacea",
             "Salvia splendens", "Mentha x piperita", "Mentha aquatica"}
C_DIP, C_POLY = "#6a3d9a", "#e08214"
src = sys.argv[1]; out = sys.argv[2] if len(sys.argv) > 2 else "Fig1a_tree_bars_ploidy"
pl_csv = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] not in ("", "-") else None
excl = [e.replace("_", " ") for e in (sys.argv[4] if len(sys.argv) > 4 else "").split(",") if e]
if pl_csv:
    _p = pd.read_csv(pl_csv); _p["species"] = _p.species.str.strip()
    POLYPLOID = set(_p[_p.ploidy.str.lower().str.strip() == "polyploid"].species)
    _known = set(_p.species)
def parse(s):
    pos = 0
    def node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1; ch = [node()]
            while s[pos] == ",": pos += 1; ch.append(node())
            pos += 1; return ch
        j = pos
        while s[j] not in ",();": j += 1
        name = s[pos:j].replace("_", " "); pos = j; return name
    return node()
def prune(t, drop):
    if isinstance(t, str): return None if t in drop else t
    ch = [c for c in (prune(c, drop) for c in t) if c is not None]
    return None if not ch else (ch[0] if len(ch) == 1 else ch)
t = pd.read_csv(src)
t["label"] = t.species.apply(lambda s: re.sub(r"^GC[AF]_\d+\.\d+_", "", s).replace("_", " "))
t = t[~t.label.apply(lambda s: any(e in s for e in excl))]
counts = dict(zip(t.label, t.n_nbarc_genes))
tree = parse(NEWICK.rstrip(";"))
all_leaves = []
def collect(n, acc):
    if isinstance(n, str): acc.append(n)
    else: [collect(c, acc) for c in n]
collect(tree, all_leaves)
extra = [l for l in counts if l not in all_leaves]
if extra: sys.exit("species in the table but not in the tree: %s" % extra)
drop = [l for l in all_leaves if l not in counts]
if drop: print("pruned from the tree (not in the table or excluded):", drop)
tree = prune(tree, set(drop))
leaves = []
collect(tree, leaves)
ypos = {l: i for i, l in enumerate(leaves)}
def height(n): return 0 if isinstance(n, str) else 1 + max(height(c) for c in n)
H = height(tree)
plt.rcParams["font.family"] = "DejaVu Sans"
n = len(leaves)
fig = plt.figure(figsize=(8.2, 0.215 * n + 1.2))
axT = fig.add_axes([0.01, 0.07, 0.19, 0.90]); axB = fig.add_axes([0.40, 0.07, 0.58, 0.90])
def draw(nd):
    if isinstance(nd, str): return ypos[nd], 0
    cs = [draw(c) for c in nd]; h = 1 + max(c[1] for c in cs); ys = [c[0] for c in cs]
    for (y, ch) in cs: axT.plot([-h, -ch], [y, y], color="#444444", lw=1.0, solid_capstyle="butt")
    axT.plot([-h, -h], [min(ys), max(ys)], color="#444444", lw=1.0, solid_capstyle="butt")
    return (min(ys) + max(ys)) / 2, h
draw(tree)
axT.set_xlim(-H - 0.3, 0); axT.set_ylim(n - 0.4, -0.6); axT.axis("off")
vals = [counts[l] for l in leaves]
cols = [C_POLY if l in POLYPLOID else C_DIP for l in leaves]
axB.barh(range(n), vals, color=cols, height=0.72)
axB.set_ylim(n - 0.4, -0.6); axB.set_yticks(range(n)); axB.set_yticklabels([])
axB.tick_params(axis="y", length=0)
axB.set_xlabel("NB-ARC genes (ANNEVO)", fontsize=9)
for s in ("top", "right"): axB.spines[s].set_visible(False)
axB.tick_params(axis="x", labelsize=8)
# species labels, right-aligned just left of the bars axis
for i, l in enumerate(leaves):
    axB.text(-0.012, i, l, ha="right", va="center", fontsize=7, transform=axB.get_yaxis_transform())
axB.legend(handles=[Patch(color=C_DIP, label="Diploid"), Patch(color=C_POLY, label="Polyploid")], title="Ploidy",
           frameon=False, fontsize=8, title_fontsize=8, loc="upper right")
fig.patch.set_facecolor("white")
fig.savefig(out + ".png", dpi=300); fig.savefig(out + ".pdf")
pd.DataFrame({"species": leaves, "n_nbarc_genes": vals, "ploidy": ["Polyploid" if l in POLYPLOID else "Diploid" for l in leaves]}).to_csv(out + "_source_data.csv", index=False)
if pl_csv:
    unk = [l for l in leaves if l not in _known]
    if unk: print("WARNING: no ploidy in PLOIDY_CSV for:", unk, "(drawn as diploid)")
print("polyploid species drawn:", sorted(l for l in leaves if l in POLYPLOID))
print("species: %d | polyploid: %d | genes: median %d, range %d-%d" % (n, sum(l in POLYPLOID for l in leaves), sorted(vals)[n // 2], min(vals), max(vals)))
