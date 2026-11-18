#!/usr/bin/env python3
"""
add_combined.py — merge the predicted and curated layers of one genome into a
single "combined" network, keeping edge provenance (predicted-only / curated-only
/ both), the curated sign and PubMed references, and the predicted confidence.

Usage:
    python3 add_combined.py --genome <accession> -o <build_dir>

Reads <build_dir>/data/<genome>.json (predicted) and
<build_dir>/data/<genome>.curated.json (curated); writes
<build_dir>/data/<genome>.combined.json and flags the manifest entry.
Node set = union of both layers (all nodes).
"""
import os, json, argparse

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", required=True)
    ap.add_argument("-o", "--out", default="klebsiella-red")
    args = ap.parse_args()
    datadir = os.path.join(args.out, "data")
    pred = json.load(open(os.path.join(datadir, args.genome + ".json")))
    cur  = json.load(open(os.path.join(datadir, args.genome + ".curated.json")))

    pids = [r[0] for r in pred["n"]]
    cids = [r[0] for r in cur["n"]]

    # merge edges by (src,tgt) -> {src: p/c/b, conf, sign, pmids, ev}
    edges = {}
    for s, t, conf in pred["e"]:
        edges[(pids[s], pids[t])] = {"src": "p", "conf": conf, "sign": None, "pmids": [], "ev": ""}
    cur_ev = cur.get("ev", {})
    for ei, e in enumerate(cur["e"]):
        k = (cids[e[0]], cids[e[1]])
        sign = e[2]; pmids = e[3] if len(e) > 3 else []
        ev = cur_ev.get(str(ei), "")
        if k in edges and edges[k]["src"] == "p":
            edges[k].update({"src": "b", "sign": sign, "pmids": pmids, "ev": ev})
        else:
            edges[k] = {"src": "c", "conf": -1, "sign": sign, "pmids": pmids, "ev": ev}

    # names: predicted first, then curated
    names = {}
    names.update(cur.get("names", {}))
    names.update(pred.get("names", {}))  # predicted labels win (richer: symbol+product)

    # node set + degrees on the merged edges
    tfs = set(k[0] for k in edges)
    allnodes = sorted(set([k[0] for k in edges] + [k[1] for k in edges]))
    idx = {n: i for i, n in enumerate(allnodes)}
    outdeg = {n: 0 for n in allnodes}; indeg = {n: 0 for n in allnodes}
    for (s, t) in edges:
        if s != t:
            outdeg[s] += 1; indeg[t] += 1
    def tier(n):
        if n in tfs:
            return 0 if indeg[n] == 0 else 1
        return 2
    nodes = [[n, 1 if n in tfs else 0, outdeg[n], indeg[n], tier(n)] for n in allnodes]

    # edge list: [srcIdx, tgtIdx, srcCode, conf, sign, pmids]  (sign None for predicted-only)
    e_out = []; ev_map = {}; srccount = {"p": 0, "c": 0, "b": 0}
    for (s, t), d in edges.items():
        ei = len(e_out)
        e_out.append([idx[s], idx[t], d["src"], d.get("conf", -1), d["sign"], d["pmids"]])
        if d["ev"]:
            ev_map[str(ei)] = d["ev"]
        srccount[d["src"]] += 1

    g = {"genome": args.genome, "combined": True,
         "n": nodes, "e": e_out, "names": names, "ev": ev_map,
         "stats": {"nodes": len(nodes), "tfs": len(tfs),
                   "targets": len(allnodes) - len(tfs), "edges": len(e_out)},
         "srccount": srccount}
    path = os.path.join(datadir, args.genome + ".combined.json")
    json.dump(g, open(path, "w"), separators=(",", ":"))

    mpath = os.path.join(args.out, "manifest.json")
    if os.path.exists(mpath):
        man = json.load(open(mpath))
        for m in man["genomes"]:
            if m["id"] == args.genome:
                m["combined"] = True
        json.dump(man, open(mpath, "w"), separators=(",", ":"))

    print("combined:", g["stats"], "provenance:", srccount)
    print("wrote", path)

if __name__ == "__main__":
    main()
