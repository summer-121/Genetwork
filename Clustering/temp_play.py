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
from pubmed_mesh import fetch_mesh_terms

try:
    import pandas as pd
except ImportError:
    pd = None


if __name__ == "__main__":
    gene_ids = ["7157", "1956", "672", "2597", "55818"]
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

    print(f"Saved: {out_path}")
