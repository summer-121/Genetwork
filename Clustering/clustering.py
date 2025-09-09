from typing import Dict, List, Tuple, Optional
from itertools import combinations


def leiden_cluster(
    nodes: List[str],
    pair_weights: Dict[Tuple[str, str], int],
    *,
    resolution: float = 1.0,
    n_iterations: int = -1,
    seed: Optional[int] = None,
    min_weight: int = 1,
) -> Dict[str, int]:
    """
    Leiden 알고리즘으로 노드 군집을 찾는다.

    - 입력 노드: ``list[str]`` (그래프의 전체 노드)
    - 입력 간선 가중치: ``{(u, v): count}`` (undirected, u!=v)
      여기서 count는 예: count_node_pairs_in_docs 의 value를 사용

    반환: ``{node: community_id}``

    의존성: python-igraph, leidenalg
    """

    # 지연 import: 필요할 때만 의존성 로드
    try:
        import igraph as ig  # type: ignore
        import leidenalg as la  # type: ignore
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "Leiden clustering requires 'python-igraph' and 'leidenalg'.\n"
            "Install with: pip install python-igraph leidenalg"
        ) from e

    # 노드 정제 및 인덱싱 (중복 제거, 순서 보존)
    unique_nodes: List[str] = []
    seen_nodes = set()
    for n in nodes:
        if not isinstance(n, str):
            continue
        if n not in seen_nodes:
            seen_nodes.add(n)
            unique_nodes.append(n)

    index = {n: i for i, n in enumerate(unique_nodes)}

    # 간선 집계: (i,j) with i<j -> weight 합산
    acc: Dict[Tuple[int, int], int] = {}
    for pair, w in pair_weights.items():
        if not isinstance(pair, tuple) or len(pair) != 2:
            continue
        a, b = pair
        if not isinstance(a, str) or not isinstance(b, str):
            continue
        if a == b:
            continue  # self-loop 무시
        if a not in index or b not in index:
            continue  # 정의된 노드에 없는 경우 무시
        try:
            ww = int(w)
        except Exception:
            continue
        if ww < min_weight:
            continue
        i, j = index[a], index[b]
        if i == j:
            continue
        if i > j:
            i, j = j, i
        acc[(i, j)] = acc.get((i, j), 0) + ww

    edges = list(acc.keys())
    weights = list(acc.values())

    # 그래프 구성 (undirected)
    g = ig.Graph(n=len(unique_nodes), edges=edges, directed=False)
    g.vs["name"] = unique_nodes
    if edges:
        g.es["weight"] = weights

    # Leiden 실행 (RBConfigurationVertexPartition; 가중치 사용)
    part = la.find_partition(
        g,
        la.RBConfigurationVertexPartition,
        weights="weight" if edges else None,
        resolution_parameter=resolution,
        seed=(0 if seed is None else seed),
        n_iterations=n_iterations,
    )

    membership = part.membership  # index -> community id
    return {unique_nodes[i]: int(membership[i]) for i in range(len(unique_nodes))}


def intra_cluster_mean_weight(
    membership: Dict[str, int],
    pair_weights: Dict[Tuple[str, str], int],
) -> Dict[int, float]:
    """
    군집별(community id) 내부 가중치의 평균을 계산.

    - 결과 형태: ``{cluster_id: mean_weight}``
    - 분모: 동일 군집 내 가능한 전체 쌍 수 (n_c * (n_c - 1) / 2)
    - 규칙:
        * 동일 군집 내 존재하지 않는 간선은 가중치 0으로 간주됨(분모에 포함).
        * 간선이 없거나 노드가 1개인 군집은 평균 0 반환.
    - 입력 간선은 무방향으로 취급(양방향 키가 있으면 합산).
    """

    # 1) 군집 크기 계산 (분모: 가능한 전체 쌍 수)
    cluster_sizes: Dict[int, int] = {}
    for node, cid in membership.items():
        try:
            icid = int(cid)
        except Exception:
            continue
        cluster_sizes[icid] = cluster_sizes.get(icid, 0) + 1

    # 초기화: 합계, 분모
    sums: Dict[int, float] = {cid: 0.0 for cid in cluster_sizes}
    denoms: Dict[int, float] = {
        cid: (size * (size - 1)) / 2.0 for cid, size in cluster_sizes.items()
    }

    # 2) 간선 가중치 무방향 정규화 (canonical key: (min(a,b), max(a,b)))
    canon: Dict[Tuple[str, str], float] = {}
    for pair, w in pair_weights.items():
        if not (isinstance(pair, tuple) and len(pair) == 2):
            continue
        a, b = pair
        if not isinstance(a, str) or not isinstance(b, str) or a == b:
            continue
        try:
            ww = float(w)
        except Exception:
            continue
        key = (a, b) if a <= b else (b, a)
        canon[key] = canon.get(key, 0.0) + ww

    # 3) 동일 군집 내부 간선 합산
    for (a, b), w in canon.items():
        ca = membership.get(a)
        cb = membership.get(b)
        if ca is None or cb is None or ca != cb:
            continue
        try:
            cid = int(ca)
        except Exception:
            continue
        if cid in sums:
            sums[cid] += w

    # 4) 평균 계산 (denom == 0 이면 0)
    means: Dict[int, float] = {}
    for cid, denom in denoms.items():
        if denom > 0:
            means[cid] = sums.get(cid, 0.0) / denom
        else:
            means[cid] = 0.0

    return means


def find_pairs_below_cluster_mean(
    membership: Dict[str, int],
    pair_weights: Dict[Tuple[str, str], float],
    *,
    inclusive: bool = False,
    include_missing_pairs: bool = True,
) -> Dict[int, List[Tuple[str, str]]]:
    """
    군집 평균 가중치보다 낮은 (노드, 노드) 쌍을 군집별로 반환.

    - 입력:
      * membership: {node -> cluster_id}
      * pair_weights: {(u,v) -> weight} (무방향, 한 방향만 존재한다고 가정)
      * inclusive: True면 평균 이하(<=), False면 미만(<)
      * include_missing_pairs: True면 군집 내 모든 가능한 쌍을 고려(없는 간선은 0)

    - 출력: {cluster_id -> [(u, v), ...]} (u < v 사전순으로 정규화)
    """

    # 군집 평균(가능한 전체 쌍 nC2 기준) 계산
    cluster_mean = intra_cluster_mean_weight(membership, pair_weights)

    # 군집별 노드 목록 구성
    clusters: Dict[int, List[str]] = {}
    for node, cid in membership.items():
        try:
            icid = int(cid)
        except Exception:
            continue
        clusters.setdefault(icid, []).append(node)

    # 빠른 조회를 위해 딕셔너리 그대로 사용하되, 키 미스 시 (b,a)도 조회
    def get_weight(a: str, b: str) -> float:
        if a == b:
            return 0.0
        w = pair_weights.get((a, b))
        if w is None:
            w = pair_weights.get((b, a))
        return float(w) if w is not None else 0.0

    def below(mean_val: float, w: float) -> bool:
        return (w <= mean_val) if inclusive else (w < mean_val)

    result: Dict[int, List[Tuple[str, str]]] = {cid: [] for cid in clusters}

    if include_missing_pairs:
        # 군집 내 가능한 모든 쌍을 생성하여 없는 간선은 0으로 취급
        for cid, nodes in clusters.items():
            if len(nodes) < 2:
                continue
            mean_val = float(cluster_mean.get(cid, 0.0))
            for a, b in combinations(sorted(nodes), 2):
                w = get_weight(a, b)
                if below(mean_val, w):
                    result[cid].append((a, b))  # 이미 a<b 정규화됨
    else:
        # 존재하는 간선만 고려 (안전하게 무방향 정규화 중복 제거)
        seen_pairs: Dict[int, set] = {cid: set() for cid in clusters}
        for (a, b), w in pair_weights.items():
            if not isinstance(a, str) or not isinstance(b, str) or a == b:
                continue
            ca = membership.get(a)
            cb = membership.get(b)
            if ca is None or cb is None or ca != cb:
                continue
            try:
                cid = int(ca)
            except Exception:
                continue
            key = (a, b) if a <= b else (b, a)
            if key in seen_pairs[cid]:
                continue
            seen_pairs[cid].add(key)
            mean_val = float(cluster_mean.get(cid, 0.0))
            if below(mean_val, float(w)):
                result[cid].append(key)

    return result
