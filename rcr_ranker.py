# -*- coding: utf-8 -*-
"""
유전자명으로 PubMed 논문을 모으고(i) NIH iCite에서 RCR을 조회한 뒤(ii)
RCR 높은 순으로 정렬해 반환(iii)합니다.

pip install requests
"""

import time
import requests
from typing import List, Dict, Optional

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ICITE_API = "https://icite.od.nih.gov/api/pubs"

def _sleep(sec=0.34):
    # NCBI 권장(rate-limit 배려)
    time.sleep(sec)

def _esummary_gene(session, gene_ids: List[str], tool, email, api_key):
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

def _esearch_gene_ids(session, symbol: str, species: str, tool, email, api_key) -> List[str]:
    q = f'({symbol}[All Fields]) AND "{species}"[Organism]'
    r = session.get(
        f"{EUTILS}/esearch.fcgi",
        params={
            "db": "gene", "retmode": "json", "term": q,
            "tool": tool, **({"email": email} if email else {}),
            **({"api_key": api_key} if api_key else {}),
        },
        timeout=30,
    )
    r.raise_for_status()
    ids = r.json().get("esearchresult", {}).get("idlist", [])
    return ids

def _best_gene_ids_for_symbol(session, symbol: str, species: str, tool, email, api_key) -> List[str]:
    ids = _esearch_gene_ids(session, symbol, species, tool, email, api_key)
    if not ids:
        return []
    summ = _esummary_gene(session, ids, tool, email, api_key)
    # 정확히 같은 심볼 매칭 우선
    exact = []
    for gid in ids:
        doc = summ.get(gid, {})
        name = (doc.get("name") or doc.get("nomenclaturesymbol") or "").upper()
        if name == symbol.upper():
            exact.append(gid)
    return exact or ids

def _gene_synonyms_from_summary(summary_doc: Dict) -> List[str]:
    syns = set()
    # otheraliases: "JMJD1B, X, Y" 형태
    otheraliases = summary_doc.get("otheraliases") or ""
    for s in otheraliases.replace(";", ",").split(","):
        s = s.strip()
        if s:
            syns.add(s)
    # 공식 심볼/이전 심볼 후보도 포함
    for k in ("name", "nomenclaturesymbol"):
        v = summary_doc.get(k)
        if v:
            syns.add(str(v).strip())
    return sorted({x for x in syns if len(x) >= 2})

def _elink_gene_to_pubmed(session, gene_ids: List[str], tool, email, api_key) -> List[str]:
    """
    Gene→PubMed 큐레이션 링크에서 PMID를 수집.
    ELink JSON은 케이스/형식이 다양한 편이라(예: {"Id": "123"} 또는 "123"),
    문자열/딕셔너리 모두 안전하게 처리합니다.
    """
    pmids = set()

    def _extract_id_from_item(item):
        # item이 문자열인 경우
        if isinstance(item, str):
            return item.strip()
        # item이 숫자형인 경우
        if isinstance(item, (int, float)):
            return str(int(item))
        # item이 dict인 경우: 'id' 또는 'Id' 또는 그 외 키에서 'id' 유사 키 찾기
        if isinstance(item, dict):
            if "id" in item:
                return str(item["id"]).strip()
            if "Id" in item:
                return str(item["Id"]).strip()
            for k, v in item.items():
                if str(k).lower() == "id":
                    return str(v).strip()
        return None

    for i in range(0, len(gene_ids), 200):
        chunk = gene_ids[i:i+200]
        r = session.get(
            f"{EUTILS}/elink.fcgi",
            params={
                "dbfrom": "gene", "db": "pubmed",
                "id": ",".join(chunk),
                "linkname": "gene_pubmed",
                "retmode": "json",
                "tool": tool, **({"email": email} if email else {}),
                **({"api_key": api_key} if api_key else {}),
            },
            timeout=30,
        )
        r.raise_for_status()
        data = r.json()

        # 표준 구조: linksets → linksetdbs → (links 또는 ids)
        for ls in data.get("linksets", []):
            for ldb in ls.get("linksetdbs", []):
                # linkname 체크(혹시 구조 다르면 dbto가 'pubmed'인지로도 필터)
                if ldb.get("linkname") == "gene_pubmed" or ldb.get("dbto") == "pubmed":
                    # 1) 'links' 처리 (문자열 또는 딕셔너리 혼재 가능)
                    for item in ldb.get("links", []) or []:
                        pid = _extract_id_from_item(item)
                        if pid:
                            pmids.add(pid)
                    # 2) 일부 응답은 'ids' 키를 따로 쓸 때가 있음
                    for item in ldb.get("ids", []) or []:
                        pid = _extract_id_from_item(item)
                        if pid:
                            pmids.add(pid)

        _sleep(0.34)  # 예의상 rate-limit

    return list(pmids)

def _esearch_pubmed_text(session, terms: List[str], species: str,
                         date_from: Optional[str], retmax: int,
                         tool, email, api_key) -> List[str]:
    # ("KDM3B" OR "JMJD1B" ...) AND (Homo sapiens OR Humans) AND (날짜필터)
    term_or = " OR ".join([f'"{t}"[Title/Abstract]' for t in terms])
    species_part = f'("{species}"[Organism] OR Humans[MeSH Terms] OR human[Title/Abstract])'
    q = f'({term_or}) AND {species_part}'
    if date_from:
        q += f' AND ("{date_from}"[Date - Publication] : "3000"[Date - Publication])'
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


def papers_sorted_by_rcr(
    gene_symbol: str,
    species: str = "Homo sapiens",
    *,
    max_papers: int = 1000,          # 최대로 모을 PubMed 논문 수
    top_n: Optional[int] = None,      # 상위 N개만 반환(없으면 전부)
    include_synonyms: bool = True,    # Gene 동의어까지 텍스트 검색에 포함
    use_gene2pubmed_first: bool = True,  # 큐레이션 링크 우선
    date_from: Optional[str] = None,  # 예: "2015/01/01" 또는 "2020"
    email: Optional[str] = None,      # NCBI 권장: 연락 이메일
    api_key: Optional[str] = None,    # 선택: E-utilities API key
    tool: str = "genetwork",
) -> List[Dict]:
    """
    반환: 논문(딕셔너리) 리스트, RCR 내림차순.
    각 항목 예시 키: pmid, year, title, relative_citation_ratio, citation_count, citations_per_year, nih_percentile, provisional
    """
    session = requests.Session()

    # 1) Gene → PubMed(큐레이션) + 2) 텍스트 검색(보완)
    pmid_set = set()

    gene_ids = []
    if use_gene2pubmed_first:
        gene_ids = _best_gene_ids_for_symbol(session, gene_symbol, species, tool, email, api_key)
        if gene_ids:
            # 큐레이션된 PubMed 링크
            pmid_set.update(_elink_gene_to_pubmed(session, gene_ids, tool, email, api_key))

    # 동의어(별칭) 수집
    syn_terms = [gene_symbol]
    if include_synonyms and gene_ids:
        summ = _esummary_gene(session, gene_ids, tool, email, api_key)
        # 가장 유력한 하나만 사용
        main_doc = summ.get(gene_ids[0], {})
        syns = _gene_synonyms_from_summary(main_doc)
        # 너무 짧거나 일반단어스러운 별칭은 나중에 필요시 필터링 가능
        for s in syns:
            if s.upper() != gene_symbol.upper():
                syn_terms.append(s)

    # 텍스트 검색으로 보완
    # (큐레이션 결과가 충분해도 텍스트 검색을 추가로 합쳐 최신 누락을 줄임)
    pmids_text = _esearch_pubmed_text(
        session, terms=syn_terms, species=species, date_from=date_from,
        retmax=max_papers, tool=tool, email=email, api_key=api_key
    )
    pmid_set.update(pmids_text)

    if not pmid_set:
        return []  # 아무 것도 못 찾음

    # 최대 개수 제한
    pmids = list(pmid_set)
    if max_papers and len(pmids) > max_papers:
        pmids = pmids[:max_papers]

    # 2) iCite에서 RCR 등 메트릭 조회
    icite = _fetch_icite(session, pmids)
    icite = [d for d in icite if isinstance(d, dict)]


    # 3) 정렬(RCR None → 0으로 간주). tie-breaker로 citations_per_year/연도 사용
    def _rcr(x):
        v = x.get("relative_citation_ratio")
        try:
            return float(v) if v is not None else 0.0
        except Exception:
            return 0.0

    def _cpy(x):
        v = x.get("citations_per_year")
        try:
            return float(v) if v is not None else 0.0
        except Exception:
            return 0.0

    def _year(x):
        try:
            return int(x.get("year") or 0)
        except Exception:
            return 0

    icite_sorted = sorted(icite, key=lambda d: (_rcr(d), _cpy(d), _year(d)), reverse=True)

    # 4) 상위 N개 슬라이싱
    if top_n:
        icite_sorted = icite_sorted[:top_n]

    # 반환 필드 최소화/정돈(필요 필드만 남기고 싶으면 여기서 골라내세요)
    cleaned = []
    keep_keys = [
        "pmid", "title", "year",
        "relative_citation_ratio", "nih_percentile",
        "citation_count", "citations_per_year",
        "field_citation_rate", "provisional"
    ]
    for d in icite_sorted:
        cleaned.append({k: d.get(k) for k in keep_keys})
    return cleaned


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