#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
GeneMANIA 추천 유전자 가져오기 (동의어 정리 제외 버전)
- 전제: 입력 유전자 심볼은 이미 표준화되어 있음 (예: HGNC 심볼)
- 필요: Cytoscape 실행 + GeneMANIA 앱 설치 + CyREST(기본 1234포트)
"""

import requests
import time
from typing import List, Dict, Optional

CY_BASE = "http://127.0.0.1:1234/v1"


class CyError(RuntimeError):
    pass


def _get(url: str, **kwargs):
    r = requests.get(url, timeout=kwargs.pop("timeout", 10))
    if r.status_code >= 400:
        raise CyError(f"GET {url} -> {r.status_code}: {r.text}")
    return r.json() if r.text and "application/json" in r.headers.get("Content-Type", "") else r.text


def _post(url: str, data=None, json=None, **kwargs):
    r = requests.post(url, data=data, json=json, timeout=kwargs.pop("timeout", 30))
    if r.status_code >= 400:
        raise CyError(f"POST {url} -> {r.status_code}: {r.text}")
    try:
        return r.json()
    except Exception:
        return r.text


def check_cytoscape_alive() -> Dict:
    """Cytoscape/CyREST 동작 확인"""
    return _get(f"{CY_BASE}/")  # 버전/상태 반환


def list_network_suids() -> List[int]:
    res = _get(f"{CY_BASE}/networks")
    return res if isinstance(res, list) else []


def get_node_rows(network_suid: int) -> List[Dict]:
    # 전체 노드 테이블 행
    return _get(f"{CY_BASE}/networks/{network_suid}/tables/defaultnode/rows/")


def infer_name_column(row: Dict) -> Optional[str]:
    """노드 이름 컬럼 추정: 보통 'name' 또는 'shared name'을 사용"""
    for cand in ["name", "shared name", "shared_name", "Shared Name"]:
        if cand in row:
            return cand
    # fallback: 문자열 타입인 첫 컬럼
    for k, v in row.items():
        if isinstance(v, str):
            return k
    return None


def infer_score_column(rows: List[Dict]) -> Optional[str]:
    """점수 컬럼 추정: 이름에 'score'/'weight'/'discriminant' 포함 + 숫자형"""
    if not rows:
        return None
    keys = rows[0].keys()
    candidates = [k for k in keys if any(s in k.lower() for s in ["score", "weight", "discriminant"])]
    def _is_numeric_series(col):
        vals = [r.get(col) for r in rows]
        nums = [v for v in vals if isinstance(v, (int, float))]
        return len(nums) >= max(3, int(0.5*len(rows)))
    cands = [c for c in candidates if _is_numeric_series(c)]
    # 여러 개면 분산이 큰 것을 선택
    best, best_var = None, -1.0
    for c in cands:
        vals = [float(r[c]) for r in rows if isinstance(r.get(c), (int, float))]
        if len(vals) >= 3:
            m = sum(vals)/len(vals)
            var = sum((x-m)**2 for x in vals) / len(vals)
            if var > best_var:
                best, best_var = c, var
    return best or (cands[0] if cands else None)


def discover_genemania_command() -> Dict:
    """
    Commands Swagger에서 genemania 네임스페이스 커맨드 자동 탐색
    return: {'path': '/v1/commands/genemania/<command>', 'params': ['genes', 'organism', ...]}
    """
    swag = _get(f"{CY_BASE}/commands/swagger.json")
    paths = swag.get("paths", {})
    hits = []
    for p, meta in paths.items():
        if "/v1/commands/" in p and "genemania" in p.lower():
            # POST만 고려
            post = meta.get("post") or {}
            params = [pr.get("name") for pr in post.get("parameters", []) if pr.get("in") in ("query", "formData")]
            hits.append({"path": p, "params": params})
    if not hits:
        raise CyError("GeneMANIA 명령을 찾지 못했습니다. GeneMANIA 앱 설치/활성 여부를 확인하세요.")
    # 'gene' 관련 파라미터를 가진 엔드포인트를 우선
    def score(hit):
        params = " ".join([p.lower() for p in hit["params"]])
        s = 0
        if "gene" in params: s += 2
        if "organism" in params: s += 1
        if "max" in params or "result" in params: s += 1
        return s
    hits.sort(key=score, reverse=True)
    return hits[0]


def run_genemania_query(
    genes: List[str],
    organism: str = "Homo sapiens",
    max_resultant_genes: int = 20,
    extra_params: Optional[Dict] = None
) -> int:
    """
    GeneMANIA 'search' 명령을 직접 호출해서 네트워크 SUID 반환
    - 일부 버전에서 자동탐색이 organisms 같은 비검색 엔드포인트를 잡는 문제 회피
    - organism은 9606(인간 tax id) 또는 "Homo sapiens" 둘 다 시도
    - genes 구분자도 콤마/파이프 모두 시도
    """
    before = set(list_network_suids())
    url = "http://127.0.0.1:1234/v1/commands/genemania/search"
    genes_csv = ",".join(genes)
    genes_pipe = "|".join(genes)

    # 시도 조합(가장 흔한 조합부터)
    tries = [
        {"genes": genes_csv,  "organism": "9606",         "geneLimit": str(max_resultant_genes), "attrLimit": "10"},
        {"genes": genes_csv,  "organism": organism,       "geneLimit": str(max_resultant_genes), "attrLimit": "10"},
        {"genes": genes_pipe, "organism": "9606",         "geneLimit": str(max_resultant_genes), "attrLimit": "10"},
        {"genes": genes_pipe, "organism": organism,       "geneLimit": str(max_resultant_genes), "attrLimit": "10"},
    ]
    if extra_params:
        for t in tries:
            t.update(extra_params)

    for i, payload in enumerate(tries, start=1):
        print(f"[GeneMANIA] TRY {i}: POST {url} json={payload}")
        resp = _post(url, json=payload)  # ← JSON 바디로 보냅니다
        print(f"[GeneMANIA] RESP (first 300 chars): {str(resp)[:300]}")

        # 생성까지 대기 (최대 12초)
        for _ in range(24):
            after = set(list_network_suids())
            new_ids = list(after - before)
            if new_ids:
                print(f"[GeneMANIA] New network SUID: {new_ids[0]}")
                return new_ids[0]
            time.sleep(0.5)

    raise CyError(
        "GeneMANIA 'search' 호출 후에도 네트워크가 생성되지 않았습니다.\n"
        "- Cytoscape에서 GeneMANIA 앱/데이터셋 설치 확인 후,\n"
        "- View→Show Command Panel에서 `help genemania` 또는 CyREST Commands API에서 "
        "`genemania/search` 가 보이는지 확인하세요.\n"
        "- 그래도 실패하면 GUI로 TP53,BRCA1을 수동 검색해 데이터셋 정상 여부를 먼저 확인하세요."
    )



def recommend_genes_from_network(network_suid: int, seed_genes: List[str], top_n: int = 20) -> Dict:
    rows = get_node_rows(network_suid)
    if not rows:
        return {"recommended": [], "score_column": None}

    name_col = infer_name_column(rows[0]) or "name"
    seed_set = set(g.lower() for g in seed_genes)

    # 점수 컬럼 추정
    score_col = infer_score_column(rows)

    recs = []
    for r in rows:
        g = str(r.get(name_col, "")).strip()
        if not g:
            continue
        if g.lower() in seed_set:
            continue  # 시드는 제외
        score = r.get(score_col) if score_col else None
        recs.append({"gene": g, "score": score})

    # 점수 정렬(있으면 내림차순)
    if score_col:
        recs.sort(key=lambda x: (x["score"] is None, -(x["score"] or float("-inf"))))
    else:
        recs.sort(key=lambda x: x["gene"])

    return {
        "network_suid": network_suid,
        "score_column": score_col,
        "recommended": recs[:top_n]
    }


def genemania_recommend(
    seed_genes: List[str],
    organism: str = "Homo sapiens",
    n_result: int = 20
) -> Dict:
    """엔드투엔드: 실행 → 추천 추출"""
    check_cytoscape_alive()
    suid = run_genemania_query(seed_genes, organism=organism, max_resultant_genes=n_result)
    result = recommend_genes_from_network(suid, seed_genes, top_n=n_result)
    result.update({"seed_genes": seed_genes, "organism": organism, "n_result": n_result})
    return result


if __name__ == "__main__":
    # 사용 예시
    seeds = ["TP53", "BRCA1"]
    out = genemania_recommend(seeds, organism="Homo sapiens", n_result=20)
    print(out)
