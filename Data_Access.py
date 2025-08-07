# Data_Access.py

# Biopython 라이브러리 활용(먼저 terminal에서 pip install biopython 실행)
from Bio import Entrez

# 유전자, 종 정보 입력란
Gene_name = "p53"  # 다른 유전자로 바꾸면 동작 가능 (예: "recA", "lexA", "lacZ")
Organism = "Homo sapiens"  # 다른 유기체로 바꾸면 동작 가능 (예: "Salmonella enterica")

# NCBI에서 정보 가져오는 함수
def NCBI_Access(gene_name: str, organism: str):

        # 1. NCBI Search Qeury 생성 및 검색
        search_query = f"{gene_name}[Gene Name] AND {organism}[Organism]"
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


NCBI_Access(Gene_name, Organism)


