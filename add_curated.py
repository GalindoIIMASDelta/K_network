#!/usr/bin/env python3
"""
add_curated.py — convert a literature-curated TRN table into a curated-network
layer for the viewer, and attach it to the build as data/<genome>.curated.json.

Usage:
    python3 add_curated.py CURATED.txt --genome <accession> -o <build_dir>

The curated TRN table is tab-separated with a header; the columns used are:
    regulator (symbol), relation (activatory/inhibitory/unknown/dual),
    regulated_gene (symbol), regulator Protein_ID, regulated Protein_ID,
    pmid, evidence.
Rows repeat per supporting reference; edges are aggregated by (regulator,target),
collecting the sign and the set of PMIDs. Node id = protein accession when
available, else the gene symbol. The layer keeps ALL curated nodes.
"""
import sys, os, csv, json, argparse, re

SIGN = {"activatory": "a", "inhibitory": "r", "dual": "d"}  # else -> unknown "u"

def sign_of(rel):
    return SIGN.get((rel or "").strip().lower(), "u")

def merge_sign(cur, new):
    if cur == new:
        return cur
    if "u" in (cur, new):
        return cur if new == "u" else new
    # different definite signs -> dual
    return "d"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("curated")
    ap.add_argument("--genome", required=True, help="assembly accession this subnet belongs to")
    ap.add_argument("-o", "--out", default="klebsiella-red")
    ap.add_argument("--evlen", type=int, default=280, help="max evidence chars kept per edge")
    args = ap.parse_args()

    # --- parse rows, locate columns by header name ---
    with open(args.curated, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
    # first occurrence of each needed header
    def col(name, occ=0):
        seen = 0
        for i, h in enumerate(header):
            if h == name:
                if seen == occ:
                    return i
                seen += 1
        return None
    c_pmid = col("pmid"); c_reg = col("regulator"); c_rel = col("relation")
    c_rgene = col("regulated_gene"); c_ev = col("evidence")
    c_regprot = col("Protein_ID", 0)        # regulator protein id (first Protein_ID)
    c_tgtprot = col("Protein_ID", 3)        # regulated protein id (4th Protein_ID col, after the regulated_gene block)

    # collect raw triples
    edges = {}          # (srcid,tgtid) -> {"sign","pmids":set,"ev"}
    symbols = {}        # id -> gene symbol
    regulators = set()
    targets = set()

    with open(args.curated, encoding="utf-8", errors="replace") as fh:
        rd = csv.reader(fh, delimiter="\t")
        next(rd)
        for r in rd:
            if len(r) <= max(c_tgtprot or 0, c_rgene or 0):
                continue
            reg_sym = (r[c_reg] or "").strip()
            tgt_sym = (r[c_rgene] or "").strip()
            if not reg_sym or not tgt_sym:
                continue
            reg_prot = (r[c_regprot] or "").strip() if c_regprot is not None else ""
            tgt_prot = (r[c_tgtprot] or "").strip() if c_tgtprot is not None else ""
            sid = reg_prot if reg_prot and reg_prot != "-" else reg_sym
            tid = tgt_prot if tgt_prot and tgt_prot != "-" else tgt_sym
            if sid == tid:
                pass  # keep self-loops
            symbols.setdefault(sid, reg_sym)
            symbols.setdefault(tid, tgt_sym)
            regulators.add(sid); targets.add(tid)
            s = sign_of(r[c_rel])
            pmid = (r[c_pmid] or "").strip()
            ev = (r[c_ev] or "").strip().strip('"')
            key = (sid, tid)
            e = edges.get(key)
            if e is None:
                edges[key] = {"sign": s, "pmids": set([pmid]) if pmid else set(), "ev": ev[:args.evlen]}
            else:
                e["sign"] = merge_sign(e["sign"], s)
                if pmid: e["pmids"].add(pmid)
                if not e["ev"] and ev: e["ev"] = ev[:args.evlen]

    # --- build node table with degrees and tier ---
    allnodes = sorted(set([k[0] for k in edges] + [k[1] for k in edges]))
    idx = {n: i for i, n in enumerate(allnodes)}
    outdeg = {n: 0 for n in allnodes}; indeg = {n: 0 for n in allnodes}
    for (s, t), e in edges.items():
        if s != t:
            outdeg[s] += 1; indeg[t] += 1
    tfs = set(regulators)
    def tier(n):
        if n in tfs:
            return 0 if indeg[n] == 0 else 1
        return 2
    nodes = [[n, 1 if n in tfs else 0, outdeg[n], indeg[n], tier(n)] for n in allnodes]
    names = {n: [symbols.get(n, n), ""] for n in allnodes}  # label = gene symbol, no product

    # edges: [srcIdx, tgtIdx, signCode, [pmids]]   signCode: a/r/d/u
    e_out = []
    ev_map = {}
    for (s, t), e in edges.items():
        ei = len(e_out)
        e_out.append([idx[s], idx[t], e["sign"], sorted(e["pmids"])])
        if e["ev"]:
            ev_map[str(ei)] = e["ev"]

    nsign = {}
    for _, _, s, _ in e_out:
        nsign[s] = nsign.get(s, 0) + 1

    g = {
        "genome": args.genome, "curated": True,
        "n": nodes, "e": e_out, "names": names, "ev": ev_map,
        "stats": {"nodes": len(nodes), "tfs": len(tfs),
                  "targets": len(allnodes) - len(tfs), "edges": len(e_out)},
        "signcounts": nsign,
    }
    datadir = os.path.join(args.out, "data")
    os.makedirs(datadir, exist_ok=True)
    path = os.path.join(datadir, args.genome + ".curated.json")
    with open(path, "w") as out:
        json.dump(g, out, separators=(",", ":"))

    # flag the manifest entry
    mpath = os.path.join(args.out, "manifest.json")
    if os.path.exists(mpath):
        man = json.load(open(mpath))
        for m in man["genomes"]:
            if m["id"] == args.genome:
                m["curated"] = True
                m["curated_edges"] = len(e_out)
        json.dump(man, open(mpath, "w"), separators=(",", ":"))

    print("curated layer:", g["stats"], "sign counts:", nsign)
    print("wrote", path)

if __name__ == "__main__":
    main()
