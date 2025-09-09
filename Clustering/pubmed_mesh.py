# -*- coding: utf-8 -*-
"""
PubMed에서 PMID로 MeSH term을 가져오는 헬퍼.

사용 예시:
    from Clustering.pubmed_mesh import get_pubmed_mesh_terms
    terms = get_pubmed_mesh_terms(["12345", "67890"], email=None, api_key=None)
    # {'12345': ['Gene Expression Regulation', ...], '67890': [...]} 
"""

import time
from typing import Dict, List, Optional
import requests
import xml.etree.ElementTree as ET

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _sleep(sec: float = 0.34) -> None:
    # NCBI 권장 레이트리밋 배려
    time.sleep(sec)


def fetch_mesh_terms(
    session: requests.Session,
    pmids: List[str],
    *,
    tool: str = "genetwork",
    email: Optional[str] = None,
    api_key: Optional[str] = None,
    include_qualifiers: bool = False,
) -> Dict[str, List[str]]:
    """
    PubMed EFetch(XML)로 PMID들의 MeSH Heading을 가져와 {PMID -> MeSH term 리스트}를 반환.

    - 기본적으로 DescriptorName(주제어)만 반환하며, include_qualifiers=True이면
      "Descriptor/Qualifier" 조합도 함께 포함합니다.
    - 주어진 PMID가 MeSH가 없거나 레코드가 없으면 빈 리스트를 매핑.
    """
    out: Dict[str, List[str]] = {}

    # PMID 정제: 숫자만, 문자열화
    

    for i in range(0, len(pmids), 200):
        chunk = pmids[i:i + 200]
        r = session.get(
            f"{EUTILS}/efetch.fcgi",
            params={
                "db": "pubmed",
                "id": ",".join(chunk),
                "retmode": "xml",
                "tool": tool,
                **({"email": email} if email else {}),
                **({"api_key": api_key} if api_key else {}),
            },
            timeout=60,
        )
        r.raise_for_status()

        try:
            root = ET.fromstring(r.text)
        except ET.ParseError:
            # 파싱 실패 시 이 청크는 건너뜀
            _sleep(0.2)
            continue

        for art in root.findall('.//PubmedArticle'):
            pmid_el = art.find('./MedlineCitation/PMID')
            if pmid_el is None or pmid_el.text is None:
                continue
            pmid = pmid_el.text.strip()

            # 순서 유지 + 중복 제거
            seen = set()
            terms: List[str] = []

            mh_list = art.find('./MedlineCitation/MeshHeadingList')
            if mh_list is not None:
                for mh in mh_list.findall('./MeshHeading'):
                    desc_el = mh.find('./DescriptorName')
                    desc = (desc_el.text.strip() if (desc_el is not None and desc_el.text) else "")
                    if desc and desc not in seen:
                        seen.add(desc)
                        terms.append(desc)

                    if include_qualifiers and desc:
                        for q_el in mh.findall('./QualifierName'):
                            if q_el is not None and q_el.text:
                                combo = f"{desc}/{q_el.text.strip()}"
                                if combo and combo not in seen:
                                    seen.add(combo)
                                    terms.append(combo)

            out[pmid] = terms

        _sleep(0.34)

    return out



def merge_mesh_terms(mesh_by_pmid: Dict[str, List[str]]) -> List[str]:
    """
    PMID -> MeSH 리스트 매핑을 하나의 리스트로 병합(중복 제거).

    - 입력: ``{"PMID": ["Term1", "Term2", ...], ...}``
    - 출력: 모든 PMID의 MeSH를 중복 없이 합친 리스트(첫 등장 순서 유지)
    """
    seen = set()
    merged: List[str] = []
    for _pmid, terms in mesh_by_pmid.items():
        for t in terms:
            if t and t not in seen:
                seen.add(t)
                merged.append(t)
    return merged
