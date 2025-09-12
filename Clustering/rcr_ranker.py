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
import datetime as _dt
import calendar as _cal

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ICITE_API = "https://icite.od.nih.gov/api/pubs"

def _sleep(sec=0.34):
    # NCBI 권장(rate-limit 배려)
    time.sleep(sec)


def _json_or_debug(resp, context: str):
    """Parse JSON or print diagnostics and re-raise on failure.
    Prints status code, content-type, and first 500 chars of body to help diagnose HTML/errors.
    """
    try:
        return resp.json()
    except Exception:
        try:
            ct = resp.headers.get("Content-Type")
        except Exception:
            ct = None
        try:
            body = resp.text[:500]
        except Exception:
            body = "<no text>"
        print(f"[DEBUG] JSON decode failed at {context}. status={getattr(resp, 'status_code', '?')}, content-type={ct}")
        print(f"[DEBUG] body[:500]= {body}")
        raise


def _get_json_with_retry(session, url: str, params: dict, context: str, *, retries: int = 2, timeout: int = 30):
    last_exc = None
    for attempt in range(retries + 1):
        try:
            r = session.get(url, params=params, timeout=timeout, headers={"Accept": "application/json"})
            r.raise_for_status()
            return _json_or_debug(r, context)
        except Exception as e:
            last_exc = e
            # brief backoff and retry
            try:
                _sleep(1.0 * (attempt + 1))
            except Exception:
                pass
    # Exhausted retries
    raise last_exc


def _esearch_snapshot(
    session,
    term: str,
    *,
    tool: str,
    email: Optional[str],
    api_key: Optional[str],
    extra_params: Optional[Dict] = None,
    timeout: int = 30,
):
    """ESearch snapshot (retmax=0, usehistory=y) returning esearchresult dict."""
    params = {
        "db": "pubmed",
        "retmode": "json",
        "term": term,
        "retmax": 0,
        "usehistory": "y",
        "tool": tool,
        **({"email": email} if email else {}),
        **({"api_key": api_key} if api_key else {}),
    }
    if extra_params:
        params.update(extra_params)
    esr = _get_json_with_retry(
        session,
        f"{EUTILS}/esearch.fcgi",
        params,
        context="ESearch snapshot",
        retries=2,
        timeout=timeout,
    ).get("esearchresult", {})
    try:
        _sleep(0.34)
    except Exception:
        pass
    return esr


def _page_ids_from_snapshot(
    session,
    *,
    webenv: str,
    query_key: str,
    total: int,
    batch: int,
    tool: str,
    email: Optional[str],
    api_key: Optional[str],
    timeout: int = 30,
) -> List[str]:
    """Download UIDs from an existing snapshot, respecting PubMed 9,999 limit."""
    out: List[str] = []
    # PubMed hard limit: can only retrieve first 9,999 records via ESearch
    max_retrievable = min(int(total or 0), 9999)
    page_size = max(1, min(int(batch) if isinstance(batch, int) else 10000, 9999))
    for retstart in range(0, max_retrievable, page_size):
        js_page = _get_json_with_retry(
            session,
            f"{EUTILS}/esearch.fcgi",
            {
                "db": "pubmed",
                "retmode": "json",
                "retstart": retstart,
                "retmax": page_size,
                "usehistory": "y",
                "WebEnv": webenv,
                "query_key": query_key,
                "tool": tool,
                **({"email": email} if email else {}),
                **({"api_key": api_key} if api_key else {}),
            },
            context=f"ESearch page retstart={retstart}",
            retries=2,
            timeout=timeout,
        )
        pmids_page = js_page.get("esearchresult", {}).get("idlist", [])
        if isinstance(pmids_page, list):
            out.extend(str(p) for p in pmids_page)
        try:
            _sleep(0.34)
        except Exception:
            pass
    return out


def _ymd(y: int, m: int = 1, d: int = 1) -> str:
    return f"{y:04d}/{m:02d}/{d:02d}"


def _year_month_ranges(y: int):
    for m in range(1, 13):
        last_day = _cal.monthrange(y, m)[1]
        yield (y, m, 1, y, m, last_day)


def _collect_by_date_partitions(
    session,
    term: str,
    *,
    tool: str,
    email: Optional[str],
    api_key: Optional[str],
    batch: int,
    y_start: int,
    y_end: int,
) -> List[str]:
    """Recursively split by publication date to bypass 9,999 limit."""
    collected: List[str] = []

    def count_range(y0: int, y1: int) -> int:
        esr = _esearch_snapshot(
            session,
            term,
            tool=tool,
            email=email,
            api_key=api_key,
            extra_params={"datetype": "pdat", "mindate": _ymd(y0, 1, 1), "maxdate": _ymd(y1, 12, 31)},
        )
        try:
            return int(esr.get("count") or 0)
        except Exception:
            return 0

    def fetch_range(y0: int, y1: int):
        # Snapshot for range then page
        esr = _esearch_snapshot(
            session,
            term,
            tool=tool,
            email=email,
            api_key=api_key,
            extra_params={"datetype": "pdat", "mindate": _ymd(y0, 1, 1), "maxdate": _ymd(y1, 12, 31)},
        )
        webenv = esr.get("webenv") or esr.get("WebEnv")
        qk = esr.get("querykey") or esr.get("QueryKey")
        try:
            total = int(esr.get("count") or 0)
        except Exception:
            total = 0
        if not total or not webenv or not qk:
            return
        ids = _page_ids_from_snapshot(
            session,
            webenv=webenv,
            query_key=qk,
            total=total,
            batch=batch,
            tool=tool,
            email=email,
            api_key=api_key,
        )
        collected.extend(ids)

    # Worklist of year ranges
    stack: List[Tuple[int, int]] = [(y_start, y_end)]
    while stack:
        y0, y1 = stack.pop()
        cnt = count_range(y0, y1)
        if cnt == 0:
            continue
        if cnt <= 9999:
            fetch_range(y0, y1)
            continue
        # Needs splitting
        if y0 == y1:
            # split into months
            for (sy, sm, sd, ey, em, ed) in _year_month_ranges(y0):
                esr = _esearch_snapshot(
                    session,
                    term,
                    tool=tool,
                    email=email,
                    api_key=api_key,
                    extra_params={
                        "datetype": "pdat",
                        "mindate": _ymd(sy, sm, sd),
                        "maxdate": _ymd(ey, em, ed),
                    },
                )
                try:
                    subcnt = int(esr.get("count") or 0)
                except Exception:
                    subcnt = 0
                if subcnt == 0:
                    continue
                if subcnt <= 9999:
                    webenv = esr.get("webenv") or esr.get("WebEnv")
                    qk = esr.get("querykey") or esr.get("QueryKey")
                    if webenv and qk:
                        ids = _page_ids_from_snapshot(
                            session,
                            webenv=webenv,
                            query_key=qk,
                            total=subcnt,
                            batch=batch,
                            tool=tool,
                            email=email,
                            api_key=api_key,
                        )
                        collected.extend(ids)
                else:
                    # If a month still exceeds 9,999 (very rare), split by day
                    last_day = _cal.monthrange(y0, sm)[1]
                    for day in range(1, last_day + 1):
                        esr_d = _esearch_snapshot(
                            session,
                            term,
                            tool=tool,
                            email=email,
                            api_key=api_key,
                            extra_params={
                                "datetype": "pdat",
                                "mindate": _ymd(y0, sm, day),
                                "maxdate": _ymd(y0, sm, day),
                            },
                        )
                        try:
                            dcnt = int(esr_d.get("count") or 0)
                        except Exception:
                            dcnt = 0
                        if dcnt == 0:
                            continue
                        webenv_d = esr_d.get("webenv") or esr_d.get("WebEnv")
                        qk_d = esr_d.get("querykey") or esr_d.get("QueryKey")
                        if not webenv_d or not qk_d:
                            continue
                        ids = _page_ids_from_snapshot(
                            session,
                            webenv=webenv_d,
                            query_key=qk_d,
                            total=min(dcnt, 9999),
                            batch=batch,
                            tool=tool,
                            email=email,
                            api_key=api_key,
                        )
                        collected.extend(ids)
        else:
            mid = (y0 + y1) // 2
            stack.append((y0, mid))
            stack.append((mid + 1, y1))

    return collected

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
        esr0 = _esearch_snapshot(
            session,
            q,
            tool=tool,
            email=email,
            api_key=api_key,
            extra_params=None,
            timeout=30,
        )
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
        # 2) If total <= 9,999: page directly; else split by publication date ranges
        if total <= 9999 and webenv and query_key:
            page_size = max(1, min(int(retmax) if isinstance(retmax, int) else 10000, 9999))
            collected = _page_ids_from_snapshot(
                session,
                webenv=webenv,
                query_key=query_key,
                total=total,
                batch=page_size,
                tool=tool,
                email=email,
                api_key=api_key,
            )
            results[gid] = collected
        else:
            this_year = _dt.date.today().year
            collected = _collect_by_date_partitions(
                session,
                q,
                tool=tool,
                email=email,
                api_key=api_key,
                batch=max(1, min(int(retmax) if isinstance(retmax, int) else 10000, 9999)),
                y_start=1900,
                y_end=this_year,
            )
            results[gid] = collected

    return results


def fetch_pmids_rcr_by_gene(
    session,
    pmids_by_gid: Dict[str, List[str]],
    *,
    default: float = 0.0,
) -> Dict[str, List[Dict[str, float]]]:
    """
    _esearch_pubmed_text 결과(Result1: {gene_id: [PMID,...]})를 받아 iCite에서 RCR을 조회하고,
    {gene_id: [{PMID: RCR}, ...]} 형태로 반환합니다. gene별 PMID 순서를 유지합니다.
    """
    # 1) 전역 PMID 유니크(순서 유지, 숫자만)
    all_pmids: List[str] = []
    seen = set()
    for lst in pmids_by_gid.values():
        if not isinstance(lst, list):
            continue
        for p in lst:
            sp = str(p)
            if sp.isdigit() and sp not in seen:
                seen.add(sp)
                all_pmids.append(sp)

    # 2) iCite 조회하여 PMID -> RCR 매핑 생성
    pmid_to_rcr: Dict[str, float] = {}
    for i in range(0, len(all_pmids), 200):
        chunk = all_pmids[i:i+200]
        if not chunk:
            continue
        r = session.get(ICITE_API, params={"pmids": ",".join(chunk)}, timeout=60)
        r.raise_for_status()
        js = r.json()
        items = js if isinstance(js, list) else (
            js.get("data") or js.get("results") or js.get("articles") or js.get("pubs") or js.get("records") or []
        )
        if isinstance(items, dict):
            items = [items]
        for it in items if isinstance(items, list) else []:
            try:
                pmid = str(it.get("pmid"))
            except Exception:
                pmid = None
            if not pmid or not pmid.isdigit():
                continue
            v = it.get("relative_citation_ratio")
            try:
                rcr = float(v) if v is not None else default
            except Exception:
                rcr = default
            pmid_to_rcr[pmid] = rcr

    # 3) 출력 구성(원래 순서 유지)
    out: Dict[str, List[Dict[str, float]]] = {}
    for gid, lst in pmids_by_gid.items():
        rows: List[Dict[str, float]] = []
        if isinstance(lst, list):
            for p in lst:
                sp = str(p)
                rcr = pmid_to_rcr.get(sp, default)
                try:
                    rcr_val = float(rcr) if rcr is not None else default
                except Exception:
                    rcr_val = default
                rows.append({sp: rcr_val})
        out[gid] = rows
    return out

def top_rcr_by_gene(
    result2: Dict[str, List[Dict[str, float]]],
    top_n: int,
) -> Dict[str, List[Dict[str, float]]]:
    """
    Result2({gene_id: [{PMID: RCR}, ...]})를 입력으로 받아,
    각 gene_id마다 RCR 상위 top_n개만 남겨 반환합니다.
    반환 형태: {gene_id: [{PMID: RCR}, ...]}
    """
    out: Dict[str, List[Dict[str, float]]] = {}

    for gid, rows in result2.items():
        pairs: List[Tuple[str, float]] = []
        if isinstance(rows, list):
            for item in rows:
                if isinstance(item, dict) and item:
                    # one-item dict expected: {pmid: rcr}
                    (pmid, rcr_raw) = next(iter(item.items()))
                    try:
                        rcr_val = float(rcr_raw) if rcr_raw is not None else 0.0
                    except Exception:
                        rcr_val = 0.0
                    pairs.append((str(pmid), rcr_val))

        # sort by RCR desc and take top_n
        if top_n is not None and top_n > 0:
            pairs.sort(key=lambda t: t[1], reverse=True)
            pairs = pairs[:top_n]

        out[gid] = [{pmid: rcr} for (pmid, rcr) in pairs]

    return out


def unique_pmids_from_top_by_gene(
    result_top: Dict[str, List[Dict[str, float]]]
) -> List[str]:
    """
    top_rcr_by_gene 결과({gene_id: [{PMID: RCR}, ...]})에서
    PMID들만 순서 유지하며 중복 없이 하나의 리스트로 반환.
    """
    out: List[str] = []
    seen = set()
    for rows in result_top.values():
        if not isinstance(rows, list):
            continue
        for item in rows:
            if isinstance(item, dict) and item:
                pmid = next(iter(item.keys()))
                sp = str(pmid)
                if sp not in seen:
                    seen.add(sp)
                    out.append(sp)
    return out
