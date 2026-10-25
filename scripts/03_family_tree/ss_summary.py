import openpyxl, statistics as st
from collections import defaultdict
wb = openpyxl.load_workbook("superasterid_Data_S2.xlsx", read_only=True, data_only=True)
rows = list(wb.worksheets[0].iter_rows(values_only=True))
h = [str(x) for x in rows[0]]; print("COLUMNAS:", h)
ix = {k: i for i, k in enumerate(h)}
def num(x):
    try: return float(x)
    except: return 0.0
cols = ["Total NLR", "CC-NLR", "CC-NLR or CCR-NLR or CCG10-N", "CCR-NLR", "TIR-NLR"]
data = [r for r in rows[1:] if r[0]]
def table(keyname, label):
    g = defaultdict(list)
    for r in data: g[r[ix[keyname]]].append(r)
    print(f"\n== Por {label}: n especies | mediana Total, CC-NLR, CCR(RNL), TIR-NLR | %TIR mediano ==")
    out = []
    for k, rs in g.items():
        med = [st.median(num(r[ix[c]]) for r in rs) for c in ("Total NLR", "CC-NLR", "CCR-NLR", "TIR-NLR")]
        pt = st.median(100 * num(r[ix["TIR-NLR"]]) / max(1, num(r[ix["Total NLR"]])) for r in rs)
        out.append((k, len(rs), med, pt))
    for k, n, m, pt in sorted(out, key=lambda x: x[3]):
        print(f"{str(k):22s} n={n:3d}  tot={m[0]:6.0f} CC={m[1]:5.0f} RNL={m[2]:4.1f} TIR={m[3]:5.0f}  %TIR={pt:5.1f}")
table("Order", "ORDEN")
print("\n== Especies de Lamiales (por familia) ==")
for r in sorted((r for r in data if r[ix["Order"]] == "Lamiales"), key=lambda r: (str(r[ix["Family"]]), str(r[ix["Organism Name"]]))):
    print(f"{r[ix['Family']]:16s} {r[ix['Organism Name']]:34s} tot={num(r[ix['Total NLR']]):4.0f} CC={num(r[ix['CC-NLR']]):4.0f} RNL={num(r[ix['CCR-NLR']]):3.0f} TIR={num(r[ix['TIR-NLR']]):3.0f}")
