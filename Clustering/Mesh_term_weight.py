# -*- coding: utf-8 -*-
"""
MeSH term 유틸리티

목적
Mesh-term으로 구성된 순서쌍에 weight(언급빈도)를 부여하는 dictionary를 만드는 것을 목표로 함

flatten_mesh_terms

"""

from typing import Iterable, List, Tuple, Dict
from itertools import combinations



def flatten_mesh_terms(
    mesh_by_pmid: Dict[str, List[str]],
    *,
    unique: bool = True,
    case_insensitive: bool = False,
    sort: bool = False,
) -> List[str]:
    """
    {PMID -> [MeSH, ...]} 매핑(dictionary)에서 MeSH 문자열만 쭉 모아 하나의 리스트로 반환.

    - unique=True: 중복 제거(첫 등장 순서 유지)
    - case_insensitive=True: 대소문자 무시하고 중복 판단(원문 표기는 유지)
    - sort=True: 최종 리스트를 알파벳 순으로 정렬(기본은 입력 순서 유지)
    """
    # 중복 제거 없이 단순 병합
    if not unique:
        out: List[str] = []
        for terms in mesh_by_pmid.values():
            if not terms:
                continue
            for t in terms:
                if isinstance(t, str):
                    out.append(t)
        return sorted(out, key=(lambda s: s.casefold()) if case_insensitive else None) if sort else out

    # 중복 제거(등장 순서 유지)
    seen = set()
    out: List[str] = []
    for terms in mesh_by_pmid.values():
        if not terms:
            continue
        for t in terms:
            if not isinstance(t, str):
                continue
            key = t.casefold() if case_insensitive else t
            if key in seen:
                continue
            seen.add(key)
            out.append(t)

    if sort:
        out = sorted(out, key=(lambda s: s.casefold()) if case_insensitive else None)
    return out


def terms_to_pairs(terms: List[str]) -> List[Tuple[str, str]]:
    """Mesh_term 모아놓은 list 순서쌍으로 만드는 함수"""
    return list(combinations(terms, 2))  # [(t1, t2), (t1, t3), ...]


def count_node_pairs_in_docs(
    node_pairs: Iterable[Tuple[str, str]],
    docs: Dict[str, List[str]],
    *,
    case_insensitive: bool = True,
    dedup_within_doc: bool = True,
) -> Dict[Tuple[str, str], int]:
    """
    주어진 node_pairs(튜플)만 대상으로, 각 문서(list[str])에 두 용어가 함께 등장하면 1씩 카운트.
    - 키는 정규화된(순서 무시, 필요 시 대소문자 무시) 튜플.
    - node_pairs에 없는 쌍은 집계하지 않음(0으로 초기화 후 누적).
    """

    def norm_term(t: str) -> str:
        return t.casefold() if case_insensitive else t

    def norm_pair(a: str, b: str) -> Tuple[str, str]:
        if case_insensitive:
            return tuple(sorted((a, b), key=str.casefold))  # 순서 무시
        return (a, b) if a <= b else (b, a)

    # 1) 입력 쌍 정규화 + 초기화
    targets = {
        norm_pair(a, b)
        for (a, b) in node_pairs
        if isinstance(a, str) and isinstance(b, str) and a != b
    }
    counts: Dict[Tuple[str, str], int] = {p: 0 for p in targets}
    if not targets or not docs:
        return counts

    # 2) 문서마다 존재 여부로 카운트
    for terms in docs.values():
        if not terms:
            continue
        norm_terms = [norm_term(t) for t in terms if isinstance(t, str)]
        if not norm_terms:
            continue
        # 존재 여부 판단이 목적이므로 set으로 변환(문서 내 중복 무시)
        doc_set = set(norm_terms) if dedup_within_doc else set(norm_terms)

        if len(doc_set) >= 2:
            for a, b in combinations(doc_set, 2):
                key = norm_pair(a, b)
                if key in counts:
                    counts[key] += 1

    return counts
