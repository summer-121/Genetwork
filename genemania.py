# genemania_runner.py (견고판: SUID 비교 + 대기 + 폴백)
from __future__ import annotations
from typing import List, Optional
from pathlib import Path
import datetime, time, tempfile, os
import py4cytoscape as p4c

def run_genemania_and_export_xml(
    genes: List[str],
    out_dir: str,
    *,
    organism: str = "Homo sapiens",
    gene_limit: int = 50,
    attr_limit: int = 10,
    combining_method: Optional[str] = None,
    networks: Optional[List[str]] = None,
    filename_prefix: str = "genemania",
) -> str:
    gene_str = ",".join(sorted({g.strip() for g in genes if g and g.strip()}))
    if not gene_str:
        raise ValueError("유전자 리스트가 비었습니다.")

    out_dir_path = Path(out_dir); out_dir_path.mkdir(parents=True, exist_ok=True)

    # 1) Cytoscape 연결
    try:
        p4c.cytoscape_ping()
    except Exception as e:
        raise RuntimeError("Cytoscape에 연결할 수 없습니다. Cytoscape 실행 여부를 확인하세요.") from e

    # 2) GeneMANIA 앱 확인(없으면 설치 시도)
    try:
        if not any("GeneMANIA" in app for app in p4c.apps.get_installed_apps()):
            p4c.apps.install_app("GeneMANIA")
    except Exception:
        pass

    # 유틸: 현재 네트워크 SUID 리스트
    def suid_list():
        try:
            # py4cytoscape 1.9+ : networks.get_network_suid_list()
            return p4c.networks.get_network_suid_list()
        except Exception:
            # 구버전 호환: 이름 목록 → SUID 매핑
            names = p4c.get_network_list() or []
            suids = []
            for n in names:
                try:
                    p4c.set_current_network(n)
                    s = p4c.get_network_suid()
                    if s: suids.append(s)
                except Exception:
                    pass
            return suids

    before_suids = set(suid_list())

    # 3) GeneMANIA 실행 (organism 폴백 포함)
    def build_cmd(org_value: str) -> str:
        parts = [
            f'genemania search organism="{org_value}"' if " " in org_value else f"genemania search organism={org_value}",
            f"geneLimit={gene_limit}",
            f"attrLimit={attr_limit}",
            f'genes="{gene_str}"',
        ]
        if combining_method: parts.append(f"combiningMethod={combining_method}")
        if networks: parts.append(f'networks={",".join(networks)}')
        return " ".join(parts)

    try:
        p4c.commands.commands_run(build_cmd(organism))
    except Exception as e:
        if (isinstance(organism, str) and organism.isdigit()) or "valid organism" in str(e).lower():
            p4c.commands.commands_run(build_cmd("Homo sapiens"))
        else:
            raise RuntimeError("GeneMANIA 실행 실패(파라미터/데이터셋 확인 필요).") from e

    # 4) 새 네트워크 SUID 대기(최대 30초, 0.5초 간격)
    new_net_suid = None
    for _ in range(60):
        time.sleep(0.5)
        after_suids = set(suid_list())
        diff = list(after_suids - before_suids)
        if diff:
            # 가장 최근 생성된 것으로 보이는 SUID 선택
            new_net_suid = sorted(diff)[-1]
            break

    # 폴백: current network라도 잡기
    if not new_net_suid:
        try:
            new_net_suid = p4c.get_network_suid()
        except Exception:
            new_net_suid = None

    if not new_net_suid:
        raise RuntimeError("생성된 네트워크를 찾지 못했습니다. (데이터셋 설치/선택 여부와 genes/organism을 확인)")

    # 5) 뷰 보장
    try:
        try:
            views = p4c.view.get_network_views(new_net_suid)
        except Exception:
            views = p4c.get_network_views(new_net_suid)
    except Exception:
        views = []
    if not views:
        try:
            p4c.view.create_network_view(new_net_suid)
        except Exception:
            try:
                p4c.create_network_view(new_net_suid)
            except Exception:
                pass
        time.sleep(0.1)

    # 6) XGMML로 내보내기
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_file = str((out_dir_path / f"{filename_prefix}_{ts}").with_suffix(".xgmml"))

    def do_export(to_path: str) -> str:
        return p4c.networks.export_network(to_path, type="XGMML", network=new_net_suid)

    try:
        saved = do_export(out_file)
    except Exception:
        # 경로/권한 문제 대비 임시폴더 폴백
        tmp_file = os.path.join(tempfile.gettempdir(), Path(out_file).name)
        saved = do_export(tmp_file)

    return saved
