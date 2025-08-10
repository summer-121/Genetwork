# Data_Access.py

# Biopython 라이브러리 활용(먼저 terminal에서 라이브러리 설치)
from Bio import Entrez
import py4cytoscape as p4c
import pandas as pd

email = "1018jjkk@gmail.com" # 이메일 주소(사용자 걸로)
Entrez.email = email

# 유전자, 종 정보 입력란
Gene_name = "UvrA"  # 다른 유전자로 바꾸면 동작 가능 (예: "recA", "lexA", "lacZ")
Organism = "Escherichia coli"  # 다른 유기체로 바꾸면 동작 가능 (예: "Salmonella enterica")

# NCBI에서 정보 가져오는 함수
def NCBI_Access(Gene_name: str, organism: str):

        # 1. NCBI Search Query 생성 및 검색
        search_query = f"{Gene_name}[Gene Name] AND {organism}[Organism]"
        handle = Entrez.esearch(db="gene", term=search_query)
        record = Entrez.read(handle)
        handle.close()

        # 2. 검색 결과에서 첫 번째 Gene ID 가져오기
        gene_id = record["IdList"][0]
        print(f"Gene ID: {gene_id}")

        # 3. 유전자 요약 정보 가져오기
        handle = Entrez.esummary(db="gene", id=gene_id)
        summary = Entrez.read(handle)
        handle.close()

        # 4. 요약 내용 출력
        gene_info = summary['DocumentSummarySet']['DocumentSummary'][0]
        print("Gene Symbol:", gene_info['NomenclatureSymbol'])
        print("Description:", gene_info['Description'])
        print("Chromosome:", gene_info['Chromosome'])
        print("Map Location:", gene_info['MapLocation'])
        print("Summary:", gene_info['Summary'])


def Cytoscape_Access(gene_name: str, species: str):

    # Cytoscape 연결 확인
    try:
        p4c.cytoscape_ping()
    except Exception as e:
        print("Cytoscape 연결 실패. Cytoscape가 켜져 있는지 확인하세요.")
        return None

    # STRING 데이터 불러오기
    print(f"'{gene_name}'의 STRING 네트워크를 불러오는 중...")
    p4c.commands.commands_post(
        f'string protein query query="{gene_name}" species="{species}"'
    )

    # 노드와 엣지 데이터 가져오기
    node_table = p4c.get_table_columns('node')
    edge_table = p4c.get_table_columns('edge')

    # 상호작용 factor와 관계 유형 추출
    results = []
    for _, row in edge_table.iterrows():
        source = row['name'].split(' (pp) ')[0]
        target = row['name'].split(' (pp) ')[1]
        interaction = row.get('interaction', 'unknown')

        if gene_name in (source, target):
            factor = target if source == gene_name else source
            results.append({
                "factor": factor,
                "interaction_type": interaction
            })

    df = pd.DataFrame(results).drop_duplicates()

    if df.empty:
        print(f"⚠️ '{gene_name}'에 대한 상호작용 정보를 찾을 수 없습니다.")
    else:
        print(f"\n'{gene_name}'와 상호작용하는 factor 목록:")
        print(df.to_string(index=False))




NCBI_Access(Gene_name, Organism)

Cytoscape_Access(Gene_name, Organism)