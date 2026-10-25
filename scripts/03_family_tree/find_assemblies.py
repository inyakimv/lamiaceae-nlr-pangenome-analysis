import json, time, urllib.parse, urllib.request
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
ORGS = ["Klebsormidium nitens", "Chara braunii", "Mesotaenium endlicherianum", "Marchantia polymorpha",
        "Physcomitrium patens", "Anthoceros agrestis", "Selaginella moellendorffii", "Amborella trichopoda"]
def get(u):
    with urllib.request.urlopen(u, timeout=30) as r: return json.load(r)
for o in ORGS:
    q = urllib.parse.quote(f'"{o}"[Organism] AND latest[filter]')
    ids = get(f"{B}esearch.fcgi?db=assembly&term={q}&retmode=json&retmax=4")["esearchresult"]["idlist"]
    time.sleep(0.4)
    print(f"== {o}: {len(ids)} assemblies")
    if not ids: continue
    s = get(f"{B}esummary.fcgi?db=assembly&id={','.join(ids)}&retmode=json")["result"]
    time.sleep(0.4)
    for i in ids:
        d = s[i]
        print("  ", d.get("assemblyaccession"), "|", d.get("assemblyname"), "|", d.get("assemblystatus"),
              "|", d.get("refseq_category"), "| RefSeq:", d.get("synonym", {}).get("refseq") or "-")
