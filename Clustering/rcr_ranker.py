# -*- coding: utf-8 -*-
"""
유전자명으로 PubMed 논문을 모으고(i) NIH iCite에서 RCR을 조회한 뒤(ii)
RCR 높은 순으로 정렬해 반환(iii)합니다.

Gene_ids: List[str] 
-> (_esummary_gene) / (gene_id_to_official_and_aliases) 
-> Dict[Gene_id,dict[official,aliases:list[str]]]

pip install requests
"""

import time
import requests
from typing import List, Dict, Optional, Iterable, Tuple
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ICITE_API = "https://icite.od.nih.gov/api/pubs"

def _sleep(sec=0.34):
    # NCBI 권장(rate-limit 배려)
    time.sleep(sec)

def _esummary_gene(session, gene_ids: List[str], tool, email, api_key):
    """
    ncbi gene에서 유전자에 대한 요약 정보인 esummary 정보를 받아옴
    """
    r = session.get(
        f"{EUTILS}/esummary.fcgi",
        params={
            "db": "gene", "retmode": "json",
            "id": ",".join(gene_ids), "tool": tool,
            **({"email": email} if email else {}),
            **({"api_key": api_key} if api_key else {}),
        },
        timeout=30,
    )
    r.raise_for_status()
    return r.json().get("result", {})


def gene_id_to_official_and_aliases(esummary_result: Dict) -> Dict[str, dict]:
    """
    Gene ID를 키로 하여 {"official": 공식심볼, "aliases": [별칭,...]} 매핑을 생성.
    _esummary_gene(...)가 반환한 result 딕셔너리를 그대로 입력으로 사용.
    """
    out: Dict[str, dict] = {}
    for gid, doc in esummary_result.items():
        if gid == "uids" or not isinstance(doc, dict):
            continue
        official = (doc.get("name") or doc.get("nomenclaturesymbol") or "").strip()
        if not official:
            continue
        raw = doc.get("otheraliases") or ""
        aliases = [a.strip() for a in raw.replace(";", ",").split(",") if len(a.strip()) >= 2]
        # 중복 제거 및 자기 자신 제외
        uniq: List[str] = []
        seen = set()
        for a in aliases:
            if a and a != official and a not in seen:
                uniq.append(a)
                seen.add(a)
        out[gid] = {"official": official, "aliases": uniq}
    return out





def _esearch_pubmed_text(session, gene_map: Dict[str, dict], retmax: int,
                         tool, email, api_key) -> Dict[str, List[str]]:
    """
    gene_id_to_official_and_aliases(...) 결과를 입력으로 받아
    official과 aliases를 모두 포함한 Title/Abstract 검색어를 구성해 Gene ID와 그에 따른 PubMed PMIDs가 매핑된 Dict를 반환.

    입력 형식: {gene_id: {"official": str, "aliases": [str, ...]}, ...}
    출력 형식: Dict{Gene_id: str , [PMID1,PMID2,...]: List[str]}
    """
    results: Dict[str, List[str]] = {}

    for gid, item in gene_map.items():
        if not isinstance(item, dict):
            continue
        official = item.get("official")
        aliases = item.get("aliases") or []

        # No dedup: include official and aliases as-is (basic sanitization only)
        terms: List[str] = []
        if isinstance(official, str) and len(official.strip()) >= 2:
            terms.append(official.strip())
        for a in aliases:
            if isinstance(a, str) and len(a.strip()) >= 2:
                terms.append(a.strip())

        if not terms:
            results[gid] = []
            continue

        # Terms-only PubMed search (no species/date filters)
        term_or = " OR ".join([f'"{t}"[Title/Abstract]' for t in terms])
        q = f'({term_or})'
 
        # 1) Snapshot search to get count + WebEnv/QueryKey (usehistory)
        r0 = session.get(
            f"{EUTILS}/esearch.fcgi",
            params={
                "db": "pubmed",
                "retmode": "json",
                "term": q,
                "retmax": 0,
                "usehistory": "y",
                "tool": tool,
                **({"email": email} if email else {}),
                **({"api_key": api_key} if api_key else {}),
            },
            timeout=30,
        )
        r0.raise_for_status()
        esr0 = r0.json().get("esearchresult", {})
        try:
            total = int(esr0.get("count") or 0)
        except Exception:
            total = 0
        webenv = esr0.get("webenv") or esr0.get("WebEnv")
        query_key = esr0.get("querykey") or esr0.get("QueryKey")

        if not total:
            results[gid] = []
            try:
                _sleep(0.34)
            except Exception:
                pass
            continue

        # 2) Page through results with PubMed cap at first 9,999 records
        page_size = max(1, min(int(retmax) if isinstance(retmax, int) else 10000, 10000))
        limit = min(int(total or 0), 9999)  # last valid start index is 9998
        collected: List[str] = []
        for retstart in range(0, limit, page_size):
            r = session.get(
                f"{EUTILS}/esearch.fcgi",
                params={
                    "db": "pubmed",
                    "retmode": "json",
                    "term": q,
                    "retstart": retstart,
                    # ensure retstart+retmax does not exceed 9,999
                    "retmax": min(page_size, max(0, limit - retstart)),
                    "usehistory": "y",
                    **({"WebEnv": webenv} if webenv else {}),
                    **({"query_key": query_key} if query_key else {}),
                    "tool": tool,
                    **({"email": email} if email else {}),
                    **({"api_key": api_key} if api_key else {}),
                },
                timeout=30,
            )
            r.raise_for_status()
            pmids_page = r.json().get("esearchresult", {}).get("idlist", [])
            if isinstance(pmids_page, list):
                collected.extend(str(p) for p in pmids_page)

            try:
                _sleep(0.34)
            except Exception:
                pass

        results[gid] = collected

    return results


def merge_unique_values_lists(d: Dict[str, List[str]]) -> List[str]:
    """
    Dict[str, List[str]] 형태의 매핑에서 values의 리스트들을
    순서를 유지하면서 중복 없이 하나의 리스트로 병합합니다.

    예) {"a": ["1", "2"], "b": ["2", "3"]} -> ["1", "2", "3"]
    """
    out: List[str] = []
    seen = set()
    for lst in d.values():
        if not isinstance(lst, list):
            continue
        for item in lst:
            if item not in seen:
                seen.add(item)
                out.append(item)
    return out





## fetch_pmids_rcr_by_gene 함수 주석 처리
# def fetch_pmids_rcr_by_gene(
#     session,
#     pmids_by_gid: Dict[str, List[str]],
#     *,
#     default: float = 0.0,
# ) -> Dict[str, List[Dict[str, float]]]:
#     """
#     _esearch_pubmed_text 결과(Result1: {gene_id: [PMID,...]})를 받아 iCite에서 RCR을 조회하고,
#     {gene_id: [{PMID: RCR}, ...]} 형태로 반환합니다. gene별 PMID 순서를 유지합니다.
#     """
#     # 1) 전역 PMID 유니크(순서 유지, 숫자만)
#     all_pmids: List[str] = []
#     seen = set()
#     for lst in pmids_by_gid.values():
#         if not isinstance(lst, list):
#             continue
#         for p in lst:
#             sp = str(p)
#             if sp.isdigit() and sp not in seen:
#                 seen.add(sp)
#                 all_pmids.append(sp)
#
#     # 2) iCite 조회하여 PMID -> RCR 매핑 생성
#     pmid_to_rcr: Dict[str, float] = {}
#     for i in range(0, len(all_pmids), 200):
#         chunk = all_pmids[i:i+200]
#         if not chunk:
#             continue
#         r = session.get(ICITE_API, params={"pmids": ",".join(chunk)}, timeout=60)
#         r.raise_for_status()
#         js = r.json()
#         items = js if isinstance(js, list) else (
#             js.get("data") or js.get("results") or js.get("articles") or js.get("pubs") or js.get("records") or []
#         )
#         if isinstance(items, dict):
#             items = [items]
#         for it in items if isinstance(items, list) else []:
#             try:
#                 pmid = str(it.get("pmid"))
#             except Exception:
#                 pmid = None
#             if not pmid or not pmid.isdigit():
#                 continue
#             v = it.get("relative_citation_ratio")
#             try:
#                 rcr = float(v) if v is not None else default
#             except Exception:
#                 rcr = default
#             pmid_to_rcr[pmid] = rcr
#
#     # 3) 출력 구성(원래 순서 유지)
#     out: Dict[str, List[Dict[str, float]]] = {}
#     for gid, lst in pmids_by_gid.items():
#         rows: List[Dict[str, float]] = []
#         if isinstance(lst, list):
#             for p in lst:
#                 sp = str(p)
#                 rcr = pmid_to_rcr.get(sp, default)
#                 try:
#                     rcr_val = float(rcr) if rcr is not None else default
#                 except Exception:
#                     rcr_val = default
#                 rows.append({sp: rcr_val})
#         out[gid] = rows
#     return out

## top_rcr_by_gene 함수 주석 처리
# def top_rcr_by_gene(
#     result2: Dict[str, List[Dict[str, float]]],
#     top_n: int,
# ) -> Dict[str, List[Dict[str, float]]]:
#     """
#     Result2({gene_id: [{PMID: RCR}, ...]})를 입력으로 받아,
#     각 gene_id마다 RCR 상위 top_n개만 남겨 반환합니다.
#     반환 형태: {gene_id: [{PMID: RCR}, ...]}
#     """
#     out: Dict[str, List[Dict[str, float]]] = {}
#
#     for gid, rows in result2.items():
#         pairs: List[Tuple[str, float]] = []
#         if isinstance(rows, list):
#             for item in rows:
#                 if isinstance(item, dict) and item:
#                     # one-item dict expected: {pmid: rcr}
#                     (pmid, rcr_raw) = next(iter(item.items()))
#                     try:
#                         rcr_val = float(rcr_raw) if rcr_raw is not None else 0.0
#                     except Exception:
#                         rcr_val = 0.0
#                     pairs.append((str(pmid), rcr_val))
#
#         # sort by RCR desc and take top_n
#         if top_n is not None and top_n > 0:
#             pairs.sort(key=lambda t: t[1], reverse=True)
#             pairs = pairs[:top_n]
#
#         out[gid] = [{pmid: rcr} for (pmid, rcr) in pairs]
#
#     return out


## unique_pmids_from_top_by_gene 함수 주석 처리
# def unique_pmids_from_top_by_gene(
#     result_top: Dict[str, List[Dict[str, float]]]
# ) -> List[str]:
#     """
#     top_rcr_by_gene 결과({gene_id: [{PMID: RCR}, ...]})에서
#     PMID들만 순서 유지하며 중복 없이 하나의 리스트로 반환.
#     """
#     out: List[str] = []
#     seen = set()
#     for rows in result_top.values():
#         if not isinstance(rows, list):
#             continue
#         for item in rows:
#             if isinstance(item, dict) and item:
#                 pmid = next(iter(item.keys()))
#                 sp = str(pmid)
#                 if sp not in seen:
#                     seen.add(sp)
#                     out.append(sp)
#     return out
