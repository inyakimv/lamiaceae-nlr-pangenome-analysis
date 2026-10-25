#!/usr/bin/env python3
"""Compare the AtADR1 used in the structure analysis with the UniProt reference (Q9FW44, ADR1_ARATH, reviewed, 787 aa; TAIR AT1G33560).
Usage: python3 26_check_atadr1.py MODEL.cif        (the AF3 model of AtADR1; sequence is read from its CA atoms)
   or: python3 26_check_atadr1.py SEQ.fasta        (any FASTA with the sequence you used)
Prints length, identity and every mismatch."""
import sys
UNIPROT_Q9FW44 = ("MASFIDLFAGDITTQLLKLLALVANTVYSCKGIAERLITMIRDVQPTIREIQYSGAELSNHHQTQLGVFYEILEKARKLCEKVLRCNRWNLKHVYHANKMKDLEKQISRFLNSQILLFVLAEVCHLRVNGDRIERNMDRLLTERNDSLSFPETMMEIETVSDPEIQTVLELGKKKVKEMMFKFTDTHLFGISGMSGSGKTTLAIELSKDDDVRGLFKNKVLFLTVSRSPNFENLESCIREFLYDGVHQRKLVILDDVWTRESLDRLMSKIRGSTTLVVSRSKLADPRTTYNVELLKKDEAMSLLCLCAFEQKSPPSPFNKYLVKQVVDECKGLPLSLKVLGASLKNKPERYWEGVVKRLLRGEAADETHESRVFAHMEESLENLDPKIRDCFLDMGAFPEDKKIPLDLLTSVWVERHDIDEETAFSFVLRLADKNLLTIVNNPRFGDVHIGYYDVFVTQHDVLRDLALHMSNRVDVNRRERLLMPKTEPVLPREWEKNKDEPFDAKIVSLHTGEMDEMNWFDMDLPKAEVLILNFSSDNYVLPPFIGKMSRLRVLVIINNGMSPARLHGFSIFANLAKLRSLWLKRVHVPELTSCTIPLKNLHKIHLIFCKVKNSFVQTSFDISKIFPSLSDLTIDHCDDLLELKSIFGITSLNSLSITNCPRILELPKNLSNVQSLERLRLYACPELISLPVEVCELPCLKYVDISQCVSLVSLPEKFGKLGSLEKIDMRECSLLGLPSSVAALVSLRHVICDEETSSMWEMVKKVVPELCIEVAKKCFTVDWLDD")
T3 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I","LEU":"L","LYS":"K",
      "MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
def from_fasta(p):
    return "".join(l.strip() for l in open(p) if not l.startswith(">")).upper().replace("*", "")
def from_cif(p):
    cols, rows, inloop = [], [], False
    lines = open(p).read().splitlines()
    for i, l in enumerate(lines):
        if l.strip() == "loop_" and i + 1 < len(lines) and lines[i + 1].startswith("_atom_site."):
            j = i + 1
            while lines[j].startswith("_atom_site."): cols.append(lines[j].strip()); j += 1
            while j < len(lines) and not lines[j].startswith(("loop_", "#", "_")):
                rows.append(lines[j].split()); j += 1
            break
    ix = {c: k for k, c in enumerate(cols)}
    need = ["_atom_site.label_atom_id", "_atom_site.label_comp_id", "_atom_site.label_seq_id", "_atom_site.label_asym_id"]
    if not all(n in ix for n in need): sys.exit("unexpected cif layout: " + str(cols[:12]))
    seq, chain = {}, None
    for r in rows:
        if r[ix["_atom_site.label_atom_id"]] != "CA": continue
        a = r[ix["_atom_site.label_asym_id"]]; chain = chain or a
        if a != chain: continue
        seq[int(r[ix["_atom_site.label_seq_id"]])] = T3.get(r[ix["_atom_site.label_comp_id"]], "X")
    return "".join(seq[k] for k in sorted(seq))
f = sys.argv[1]
mine = from_fasta(f) if f.lower().endswith((".fa", ".fasta", ".faa")) else from_cif(f)
ref = UNIPROT_Q9FW44
print("reference Q9FW44: %d aa | yours: %d aa" % (len(ref), len(mine)))
if len(ref) == len(mine):
    mm = [(i + 1, a, b) for i, (a, b) in enumerate(zip(ref, mine)) if a != b]
    print("identity %.2f%% (%d/%d)" % (100 * (len(ref) - len(mm)) / len(ref), len(ref) - len(mm), len(ref)))
    for p, a, b in mm: print("  position %d: UniProt %s, yours %s" % (p, a, b))
    print("IDENTICAL" if not mm else "DIFFERENT at %d position(s)" % len(mm))
else:
    print("different length: not comparable position by position; yours starts %s... ends ...%s" % (mine[:15], mine[-15:]))
