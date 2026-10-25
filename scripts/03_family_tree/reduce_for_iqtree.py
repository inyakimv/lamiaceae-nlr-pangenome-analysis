import subprocess
ALN = "nbarc_hmm_refs_outgroups.fasta"
seqs, order, k = {}, [], None
for l in open(ALN):
    l = l.rstrip()
    if l.startswith(">"): k = l[1:].split()[0]; seqs[k] = []; order.append(k)
    else: seqs[k].append(l)
seqs = {a: "".join(b) for a, b in seqs.items()}
always = [n for n in order if n.startswith(("REF|", "OUT|", "Athaliana", "Knitens"))]
core = [n for n in order if not n.startswith(("REF|", "OUT|", "Athaliana", "Knitens"))]
print(f"Always kept: {len(always)}; Lamiaceae to cluster: {len(core)}")
with open("core.faa", "w") as o:
    for n in core: o.write(f">{n}\n{seqs[n].replace('-', '')}\n")
for ident, w in ((0.6, 4), (0.5, 3)):
    p = int(ident * 100)
    subprocess.run(["cd-hit", "-i", "core.faa", "-o", f"core_c{p}.faa", "-c", str(ident), "-n", str(w),
                    "-T", "16", "-M", "0", "-d", "0"], stdout=subprocess.DEVNULL, check=True)
    reps = [l[1:].split()[0] for l in open(f"core_c{p}.faa") if l.startswith(">")]
    with open(f"reduced_c{p}.fasta", "w") as o:
        for n in always + reps: o.write(f">{n}\n{seqs[n]}\n")
    print(f"CD-HIT {ident}: {len(reps)} Lamiaceae representatives -> reduced_c{p}.fasta with {len(always) + len(reps)} sequences")
