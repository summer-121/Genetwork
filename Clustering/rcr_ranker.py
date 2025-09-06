# -*- coding: utf-8 -*-
"""
유전자명으로 PubMed 논문을 모으고(i) NIH iCite에서 RCR을 조회한 뒤(ii)
RCR 높은 순으로 정렬해 반환(iii)합니다.

pip install requests
"""

import time
import requests
from typing import List, Dict, Optional, Iterable

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


def alias_to_gene_name_from_esummary(esummary_result: Dict) -> Dict[str, str]:
    """
    _esummary_gene(...)가 반환한 result 딕셔너리에서
    {별칭 -> 공식 심볼(name/nomenclaturesymbol)} 매핑을 만들어 반환.
    """
    out: Dict[str, str] = {}
    for gid, doc in esummary_result.items():
        if gid == "uids" or not isinstance(doc, dict):
            continue
        official = (doc.get("name") or doc.get("nomenclaturesymbol") or "").strip()
        if not official:
            continue
        raw = doc.get("otheraliases") or ""
        for part in raw.replace(";", ",").split(","):
            alias = part.strip()
            if len(alias) < 2:
                continue
            out.setdefault(alias, official)  # 같은 alias가 여러 유전자에 있으면 최초 매핑 유지
    return out

def alias_to_gene_name_from_summary_doc(summary_doc: Dict) -> Dict[str, str]:
    """
    ESummary에서 특정 Gene 문서 하나(summary_doc)만 있을 때
    {별칭 -> 공식 심볼} 매핑을 만들어 반환.
    """
    official = (summary_doc.get("name") or summary_doc.get("nomenclaturesymbol") or "").strip()
    if not official:
        return {}
    raw = summary_doc.get("otheraliases") or ""
    aliases = [s.strip() for s in raw.replace(";", ",").split(",") if len(s.strip()) >= 2]
    return {a: official for a in aliases}


def _esearch_pubmed_text(session, terms: List[str], species: str,
                         date_from: Optional[str], retmax: int,
                         tool, email, api_key) -> List[str]:
    # Terms-only PubMed search (no species/date filters)
    term_or = " OR ".join([f'"{t}"[Title/Abstract]' for t in terms])
    q = f'({term_or})'
    r = session.get(
        f"{EUTILS}/esearch.fcgi",
        params={
            "db": "pubmed", "retmode": "json", "term": q,
            "retmax": retmax, "tool": tool,
            **({"email": email} if email else {}),
            **({"api_key": api_key} if api_key else {}),
        },
        timeout=30,
    )
    r.raise_for_status()
    return r.json().get("esearchresult", {}).get("idlist", [])


def _fetch_icite(session, pmids: List[str]) -> List[Dict]:
    """
    iCite API 응답이 list 또는 dict(에러/래핑)로 올 수 있으므로
    안전하게 리스트만 추출하고, dict가 아니면 버립니다.
    429/5xx에 대한 간단한 재시도도 포함합니다.
    """
    out: List[Dict] = []

    def _normalize_payload(js):
        # 응답 형태를 일관된 리스트로 변환
        if isinstance(js, list):
            return js
        if isinstance(js, dict):
            for key in ("data", "results", "articles", "pubs", "records"):
                v = js.get(key)
                if isinstance(v, list):
                    return v
            # dict이지만 리스트 래핑이 없다면, 단일 레코드일 수도 있으므로 감싸기
            # 단, 최소한 pmid 같은 필드가 있는 dict일 때만
            if any(k in js for k in ("pmid", "relative_citation_ratio", "title")):
                return [js]
            return []
        # 그 외(str 등)는 무시
        return []

    def _is_good_item(x):
        return isinstance(x, dict) and ("pmid" in x or "relative_citation_ratio" in x)

    # iCite가 한 번에 200개까지 권장
    for i in range(0, len(pmids), 200):
        chunk = pmids[i:i+200]
        # 숫자만 남기기(비정상 ID 제거)
        chunk = [str(p) for p in chunk if str(p).isdigit()]
        if not chunk:
            continue

        # 간단 재시도(최대 3회)
        backoff = 0.5
        for attempt in range(3):
            try:
                r = session.get(ICITE_API, params={"pmids": ",".join(chunk)}, timeout=60)
                # 429/5xx면 재시도
                if r.status_code >= 500 or r.status_code == 429:
                    time.sleep(backoff)
                    backoff *= 2
                    continue
                r.raise_for_status()
                js = r.json()
                items = _normalize_payload(js)
                out.extend([it for it in items if _is_good_item(it)])
                break
            except requests.RequestException:
                time.sleep(backoff)
                backoff *= 2
        time.sleep(0.2)  # 예의상 rate-limit

    return out



def sort_icite_by_rcr(items: Iterable[Dict], top_n: Optional[int] = None,
                      keep_only_known_keys: bool = True) -> List[Dict]:
    """
    iCite 결과(딕셔너리들의 반복자)를 RCR 기준으로 내림차순 정렬합니다.
    - 정렬 키: (relative_citation_ratio, citations_per_year, year)
    - top_n: 상위 N개만 반환(없으면 전체)
    - keep_only_known_keys: True면 공통 키만 추려서 반환
    """
    def _f(d: Dict, key: str, default: float = 0.0) -> float:
        v = d.get(key)
        try:
            return float(v) if v is not None else default
        except Exception:
            return default

    def _year(d: Dict) -> int:
        try:
            return int(d.get("year") or 0)
        except Exception:
            return 0

    rows = [x for x in items if isinstance(x, dict)]
    rows.sort(key=lambda d: (_f(d, "relative_citation_ratio"),
                             _f(d, "citations_per_year"),
                             _year(d)),
              reverse=True)

    if top_n:
        rows = rows[:top_n]

    if not keep_only_known_keys:
        return rows

    keep_keys = [
        "pmid", "title", "year",
        "relative_citation_ratio", "nih_percentile",
        "citation_count", "citations_per_year",
        "field_citation_rate", "provisional",
    ]
    return [{k: d.get(k) for k in keep_keys} for d in rows]


"""
# play.py 이건 실행파일임
from rcr_ranker import papers_sorted_by_rcr

if __name__ == "__main__":
    rows = papers_sorted_by_rcr(
        "KDM3B", species="Homo sapiens",
        max_papers=1000, top_n=20,
        include_synonyms=True,
        use_gene2pubmed_first=True,
        date_from="2015",
        email="you@example.com",
        api_key=None,
    )
    for i, r in enumerate(rows, 1):
        print(f"{i:>2}. PMID {r['pmid']} | RCR {r['relative_citation_ratio']} | {r['year']} | {r['title']}")

"""
