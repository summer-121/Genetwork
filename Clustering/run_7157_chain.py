import os
import sys
import argparse
import csv
import requests

# Allow local imports from this folder
sys.path.append(os.path.dirname(__file__))
from rcr_ranker import (  # type: ignore
    _esummary_gene,
    gene_id_to_official_and_aliases,
    _esearch_pubmed_text,
    merge_unique_values_lists,
)
from pubmed_mesh import fetch_mesh_terms, merge_mesh_terms  # type: ignore
from Mesh_term_weight import terms_to_pairs, count_node_pairs_in_docs  # type: ignore
from clustering import (  # type: ignore
    leiden_cluster,
    intra_cluster_mean_weight,
    find_pairs_below_cluster_mean,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Chain: _esummary_gene -> gene_id_to_official_and_aliases -> _esearch_pubmed_text -> merge_unique_values_lists for gene 7157",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--gene-ids", type=str, default="7157", help="Comma-separated NCBI Gene IDs")
    p.add_argument("--retmax", type=int, default=1000, help="Max PMIDs to retrieve per gene")
    p.add_argument("--tool", type=str, default="GENETWORK")
    p.add_argument("--email", type=str, default=os.environ.get("NCBI_EMAIL"))
    p.add_argument("--api-key", type=str, default=os.environ.get("NCBI_API_KEY"))
    return p.parse_args()


def main() -> None:
    args = parse_args()
    gene_ids = [g.strip() for g in args.gene_ids.split(",") if g.strip()]
    if not gene_ids:
        print("No gene IDs provided.")
        return

    with requests.Session() as s:
        # 1) ESummary for gene metadata
        esummary = _esummary_gene(s, gene_ids, tool=args.tool, email=args.email, api_key=args.api_key)
        print(f"[INFO] ESummary received for {len(gene_ids)} gene IDs")

        # 2) Build official + aliases mapping
        gene_map = gene_id_to_official_and_aliases(esummary)
        print(f"[INFO] gene_map built for {len(gene_map)} genes: {list(gene_map.keys())}")
        for gid, d in gene_map.items():
            print(f"  - {gid}: official='{d.get('official')}', aliases={len(d.get('aliases') or [])}")

        # 3) PubMed ESearch (Title/Abstract over official+aliases)
        pmids_by_gid = _esearch_pubmed_text(
            s,
            gene_map,
            retmax=args.retmax,
            tool=args.tool,
            email=args.email,
            api_key=args.api_key,
        )
        for gid, lst in pmids_by_gid.items():
            print(f"[INFO] PMIDs for {gid}: {len(lst)}")

    # 4) Merge unique PMIDs across genes
    unique_pmids = merge_unique_values_lists(pmids_by_gid)
    print(f"[INFO] Unique PMIDs total: {len(unique_pmids)}")
    if unique_pmids:
        preview = ", ".join(unique_pmids[:10])
        if len(unique_pmids) > 10:
            preview += ", ..."
        print(f"[INFO] Preview PMIDs: {preview}")

    # 5) Fetch MeSH terms for unique PMIDs
    if not unique_pmids:
        return
    with requests.Session() as s:
        mesh_by_pmid = fetch_mesh_terms(
            s,
            unique_pmids,
            tool=args.tool,
            email=args.email,
            api_key=args.api_key,
            include_qualifiers=False,
        )
    print(f"[INFO] MeSH fetched for {len(mesh_by_pmid)} PMIDs")
    # 6) Merge MeSH terms into a unique list
    unique_terms = merge_mesh_terms(mesh_by_pmid)
    print(f"[INFO] Unique MeSH terms: {len(unique_terms)}")
    if unique_terms:
        tprev = ", ".join(unique_terms[:10])
        if len(unique_terms) > 10:
            tprev += ", ..."
        print(f"[INFO] Preview terms: {tprev}")

    # 7) Build term pairs and count co-occurrence; write CSV outputs
    candidate_pairs = terms_to_pairs(unique_terms)
    print(f"[INFO] Candidate pairs (nC2): {len(candidate_pairs)}")

    pair_counts = count_node_pairs_in_docs(
        candidate_pairs,
        mesh_by_pmid,
        case_insensitive=False,
        dedup_within_doc=True,
    )
    nonzero = sum(1 for v in pair_counts.values() if v)
    print(f"[INFO] Nonzero pair counts: {nonzero} (of {len(pair_counts)})")

    # Output CSVs under ../data/
    out_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "data"))
    os.makedirs(out_dir, exist_ok=True)
    label = "_".join(gene_ids)
    counts_csv = os.path.join(out_dir, f"{label}_pair_counts.csv")

    with open(counts_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["u", "v", "weight"])
        for (a, b), wgt in pair_counts.items():
            w.writerow([a, b, int(wgt)])

    print(f"[INFO] Wrote: {counts_csv}")

    # 8) Run clustering with candidate pairs and pair_counts
    nodes = sorted({t for (a, b) in candidate_pairs for t in (a, b)})
    try:
        membership = leiden_cluster(
            nodes,
            pair_counts,
            resolution=1.0,
            n_iterations=-1,
            seed=42,
            min_weight=1,
        )
    except RuntimeError as e:
        print("[ERROR] Clustering skipped:", e)
        print("[HINT] pip install python-igraph leidenalg   (or conda-forge)")
        return

    n_clusters = len(set(membership.values()))
    print(f"[INFO] Clustering complete: {len(membership)} nodes, {n_clusters} clusters")

    cluster_means = intra_cluster_mean_weight(membership, pair_counts)
    pairs_below = find_pairs_below_cluster_mean(
        membership,
        pair_counts,
        inclusive=False,
        include_missing_pairs=True,
    )

    # Write clustering outputs as CSV
    membership_csv = os.path.join(out_dir, f"{label}_membership.csv")
    means_csv = os.path.join(out_dir, f"{label}_cluster_means.csv")
    below_csv = os.path.join(out_dir, f"{label}_pairs_below_mean.csv")

    with open(membership_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["node", "cluster_id"])
        for n, cid in sorted(membership.items(), key=lambda x: (x[1], x[0])):
            w.writerow([n, int(cid)])

    with open(means_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cluster_id", "mean_weight"])
        for cid in sorted(cluster_means):
            w.writerow([int(cid), float(cluster_means[cid])])

    with open(below_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["cluster_id", "u", "v"])
        for cid in sorted(pairs_below):
            for (u, v) in sorted(pairs_below[cid]):
                w.writerow([int(cid), u, v])

    print(f"[INFO] Wrote: {membership_csv}")
    print(f"[INFO] Wrote: {means_csv}")
    print(f"[INFO] Wrote: {below_csv}")


if __name__ == "__main__":
    main()
