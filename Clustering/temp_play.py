import json
import os
import sys
import requests

# 기존 파이프라인 코드 비활성화
# sys.path.append(os.path.dirname(__file__))
# from rcr_ranker import (
#     _esummary_gene,
#     gene_id_to_official_and_aliases,
#     _esearch_pubmed_text,
#     fetch_pmids_rcr_by_gene,
#     top_rcr_by_gene,
# )
# try:
#     import pandas as pd
# except ImportError:
#     pd = None
# if __name__ == "__main__":
#     ...

# 새 파이프라인: rcr_ranker → pubmed_mesh → Excel
sys.path.append(os.path.dirname(__file__))
from rcr_ranker import (
    _esummary_gene,
    gene_id_to_official_and_aliases,
    _esearch_pubmed_text,
    fetch_pmids_rcr_by_gene,
    top_rcr_by_gene,
    unique_pmids_from_top_by_gene,
)
from pubmed_mesh import fetch_mesh_terms, merge_mesh_terms
from Mesh_term_weight import count_node_pairs_in_docs
from clustering import (
    leiden_cluster,
    intra_cluster_mean_weight,
    find_pairs_below_cluster_mean,
)
from itertools import combinations

try:
    import pandas as pd
except ImportError:
    pd = None


if __name__ == "__main__":
    gene_ids = ["7157", "1956", "672", "675", "3845", "4893", "673", "5728", "5290", "207", "4609", "324", "4193", "5925", "1029", "2625", "2064", "5594", "4089", "1499"
]
    top_n = 20

    with requests.Session() as s:
        esummary = _esummary_gene(s, gene_ids, tool="GENETWORK", email=None, api_key=None)
        gene_map = gene_id_to_official_and_aliases(esummary)
        pmids_by_gid = _esearch_pubmed_text(
            s, gene_map, retmax=1000, tool="GENETWORK", email=None, api_key=None
        )
        result2 = fetch_pmids_rcr_by_gene(s, pmids_by_gid, default=0.0)
        top_by_gene = top_rcr_by_gene(result2, top_n=top_n)

        # MeSH terms for unique PMIDs from top_by_gene
        unique_pmids = unique_pmids_from_top_by_gene(top_by_gene)
        mesh_by_pmid = fetch_mesh_terms(
            s, unique_pmids, tool="GENETWORK", email=None, api_key=None, include_qualifiers=False
        )

    # ----- Mesh_term_weight 단계: 쌍 집계 -----
    # 1) 전체 고유 MeSH 용어 수집 (표시용)
    unique_terms = list(dict.fromkeys(merge_mesh_terms(mesh_by_pmid)))

    # 2) 후보 쌍(targets): 문서 내에서 함께 등장한 적 있는 쌍들의 합집합
    candidate_pairs = set()
    for terms in mesh_by_pmid.values():
        if not terms:
            continue
        # 문서 내 중복 제거 후 가능한 쌍 생성
        doc_set = sorted(set(t for t in terms if isinstance(t, str) and t))
        if len(doc_set) >= 2:
            for a, b in combinations(doc_set, 2):
                candidate_pairs.add((a, b))

    # 3) 문서 집합에서 각 쌍의 동시 등장 횟수 카운트 (대소문자 무시, 문서 내 중복 무시)
    pair_counts = count_node_pairs_in_docs(
        candidate_pairs,
        mesh_by_pmid,
        case_insensitive=True,
        dedup_within_doc=True,
    )

    # ----- clustering 단계: Leiden + 평균 계산 + 평균 이하 간선 -----
    # 노드 목록(정규화: casefold) 구성 — pair_counts의 키가 이미 정규화되어 있음
    norm_nodes = sorted({t for pair in pair_counts.keys() for t in pair})

    membership = leiden_cluster(
        norm_nodes,
        pair_counts,
        resolution=1.0,
        n_iterations=-1,
        seed=42,
        min_weight=1,
    )

    cluster_means = intra_cluster_mean_weight(membership, pair_counts)
    pairs_below = find_pairs_below_cluster_mean(
        membership,
        pair_counts,
        inclusive=False,            # 필요 시 True로 변경 가능(이하 포함)
        include_missing_pairs=True, # 존재하지 않는 간선은 0으로 간주
    )
    # pandas가 없으면 JSON 파일로만 저장
    if pd is None:
        print("pandas가 설치되어 있지 않습니다. pip install pandas openpyxl 후 다시 실행하세요.")
        base = os.path.join(os.path.dirname(__file__), "pipeline_output")
        with open(base + "_esummary.json", "w", encoding="utf-8") as f:
            json.dump(esummary, f, ensure_ascii=False, indent=2)
        with open(base + "_gene_map.json", "w", encoding="utf-8") as f:
            json.dump(gene_map, f, ensure_ascii=False, indent=2)
        with open(base + "_pmids_by_gid.json", "w", encoding="utf-8") as f:
            json.dump(pmids_by_gid, f, ensure_ascii=False, indent=2)
        with open(base + "_result2.json", "w", encoding="utf-8") as f:
            json.dump(result2, f, ensure_ascii=False, indent=2)
        with open(base + "_top_by_gene.json", "w", encoding="utf-8") as f:
            json.dump(top_by_gene, f, ensure_ascii=False, indent=2)
        with open(base + "_mesh_by_pmid.json", "w", encoding="utf-8") as f:
            json.dump(mesh_by_pmid, f, ensure_ascii=False, indent=2)
        raise SystemExit(0)

    # DataFrame 유틸
    def df_esummary(es):
        rows = []
        for gid, doc in es.items():
            if not isinstance(doc, dict) or gid == "uids":
                continue
            row = {"gene_id": gid}
            for k in [
                "name",
                "nomenclaturesymbol",
                "otheraliases",
                "description",
                "genetype",
                "chromosome",
                "maplocation",
            ]:
                if k in doc:
                    row[k] = doc.get(k)
            rows.append(row)
        return pd.DataFrame(rows)

    def df_gene_map(gm):
        rows = []
        for gid, d in gm.items():
            rows.append(
                {
                    "gene_id": gid,
                    "official": d.get("official"),
                    "aliases": "|".join(d.get("aliases") or []),
                }
            )
        return pd.DataFrame(rows)

    def df_pmids_by_gid(pb):
        rows = []
        for gid, lst in pb.items():
            for p in lst:
                rows.append({"gene_id": gid, "pmid": str(p)})
        return pd.DataFrame(rows)

    def df_result2(r2):
        rows = []
        for gid, lst in r2.items():
            for d in lst:
                if isinstance(d, dict) and d:
                    pmid, rcr = next(iter(d.items()))
                    rows.append({"gene_id": gid, "pmid": str(pmid), "rcr": float(rcr)})
        return pd.DataFrame(rows)

    def df_mesh(mesh_map):
        rows = []
        for pmid, terms in mesh_map.items():
            rows.append({"pmid": str(pmid), "mesh_terms": "|".join(terms or [])})
        return pd.DataFrame(rows)

    out_path = os.path.join(os.path.dirname(__file__), "pipeline_output.xlsx")
    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        df_esummary(esummary).to_excel(writer, index=False, sheet_name="esummary")
        df_gene_map(gene_map).to_excel(writer, index=False, sheet_name="gene_map")
        df_pmids_by_gid(pmids_by_gid).to_excel(writer, index=False, sheet_name="pmids_by_gid")
        df_result2(result2).to_excel(writer, index=False, sheet_name="pmids_rcr_by_gene")
        df_result2(top_by_gene).to_excel(writer, index=False, sheet_name="top_by_gene")
        df_mesh(mesh_by_pmid).to_excel(writer, index=False, sheet_name="mesh_terms")

        # ----- 추가 결과: Mesh_term_weight / clustering -----
        # 고유 용어 목록
        pd.DataFrame({"term": unique_terms}).to_excel(writer, index=False, sheet_name="unique_terms")

        # 쌍-가중치 표
        if pair_counts:
            rows_pairs = [
                {"u": a, "v": b, "weight": int(w)} for (a, b), w in pair_counts.items()
            ]
            pd.DataFrame(rows_pairs).to_excel(writer, index=False, sheet_name="pair_counts")
        else:
            pd.DataFrame(columns=["u", "v", "weight"]).to_excel(writer, index=False, sheet_name="pair_counts")

        # 군집 멤버십(노드→군집)
        if membership:
            rows_mem = [{"node": n, "cluster_id": cid} for n, cid in membership.items()]
            pd.DataFrame(rows_mem).to_excel(writer, index=False, sheet_name="membership")
        else:
            pd.DataFrame(columns=["node", "cluster_id"]).to_excel(writer, index=False, sheet_name="membership")

        # 군집별 평균 가중치
        if cluster_means:
            rows_means = [{"cluster_id": cid, "mean_weight": float(m)} for cid, m in cluster_means.items()]
            pd.DataFrame(rows_means).to_excel(writer, index=False, sheet_name="cluster_means")
        else:
            pd.DataFrame(columns=["cluster_id", "mean_weight"]).to_excel(writer, index=False, sheet_name="cluster_means")

        # 군집별 평균보다 낮은 간선 쌍
        def get_weight(a: str, b: str) -> float:
            w = pair_counts.get((a, b))
            if w is None:
                w = pair_counts.get((b, a))
            return float(w) if w is not None else 0.0

        rows_below = []
        for cid, pairs in pairs_below.items():
            mean_val = float(cluster_means.get(cid, 0.0))
            for a, b in pairs:
                rows_below.append({
                    "cluster_id": cid,
                    "u": a,
                    "v": b,
                    "weight": get_weight(a, b),
                    "cluster_mean": mean_val,
                })
        pd.DataFrame(rows_below or [], columns=["cluster_id", "u", "v", "weight", "cluster_mean"]).to_excel(
            writer, index=False, sheet_name="pairs_below_mean"
        )

    print(f"Saved: {out_path}")
