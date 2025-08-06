import PyPDF2
import re
import requests
import json
from collections import Counter
from typing import List, Dict, Tuple
import os


class GeneAnalyzer:
    """PDF 논문에서 유전자를 추출하고 분석하는 클래스"""
    
    def __init__(self):
        self.gene_database = set()  # 모든 유전자 이름을 저장할 집합
        self.load_gene_database()
    
    def load_gene_database(self):
        """
        실제 유전자 데이터베이스를 로드
        여러 소스에서 유전자 정보를 가져옴
        """
        print("유전자 데이터베이스를 로딩 중...")
        
        # 1. 기본 유명 유전자들 (하드코딩)
        common_genes = {
            # 종양억제유전자
            'TP53', 'RB1', 'BRCA1', 'BRCA2', 'APC', 'VHL', 'NF1', 'NF2',
            'CDKN2A', 'PTEN', 'ATM', 'PALB2', 'MLH1', 'MSH2', 'MSH6', 'PMS2',
            
            # 종양유전자 (Oncogenes)
            'MYC', 'RAS', 'KRAS', 'NRAS', 'HRAS', 'EGFR', 'HER2', 'ERBB2',
            'PIK3CA', 'AKT1', 'BRAF', 'RAF1', 'JUN', 'FOS', 'BCL2', 'CCND1',
            
            # 전사인자
            'P53', 'MYC', 'FOS', 'JUN', 'ETS1', 'SP1', 'NFE2L2', 'STAT3',
            'CREB1', 'FOXO1', 'SMAD4', 'TCF7L2', 'RUNX1', 'GATA1',
            
            # 성장인자 및 수용체
            'VEGF', 'VEGFA', 'PDGF', 'FGF', 'IGF1', 'IGF2', 'TGFB1', 'EGF',
            'CTNNB1', 'WNT1', 'WNT3A', 'FZD1', 'NOTCH1', 'NOTCH2',
            
            # 세포주기 관련
            'CDK1', 'CDK2', 'CDK4', 'CDK6', 'CCNA1', 'CCNA2', 'CCNB1', 'CCND1',
            'CCNE1', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B',
            
            # DNA 복구 관련
            'XRCC1', 'XRCC3', 'RAD51', 'RAD52', 'BRIP1', 'FANCA', 'FANCD2',
            'MLH1', 'MSH2', 'MSH6', 'PMS2', 'MUTYH',
            
            # 면역 관련
            'CD4', 'CD8A', 'CD8B', 'CD19', 'CD20', 'CD3E', 'IFNG', 'IL2',
            'IL4', 'IL6', 'IL10', 'TNF', 'TNFRSF1A', 'CTLA4', 'PDCD1',
            
            # 대사 관련
            'INS', 'INSR', 'IRS1', 'PPARA', 'PPARG', 'SREBF1', 'FASN',
            'ACC1', 'SCD1', 'CPT1A', 'ACACA', 'ACACB',
            
            # 신호전달 관련
            'PI3K', 'AKT', 'MTOR', 'TSC1', 'TSC2', 'RHEB', 'RPS6KB1',
            'EIF4EBP1', 'GSK3B', 'FOXO3', 'AMPK', 'LKB1',
            
            # 후성유전학 관련
            'DNMT1', 'DNMT3A', 'DNMT3B', 'TET1', 'TET2', 'IDH1', 'IDH2',
            'HDAC1', 'HDAC2', 'SIRT1', 'SIRT3', 'EZH2', 'SUV39H1',
            
            # 줄기세포 관련
            'OCT4', 'POU5F1', 'SOX2', 'NANOG', 'KLF4', 'MYC', 'LIN28A',
            'CDX2', 'GATA4', 'GATA6', 'SOX17', 'FOXA2',
            
            # 세포사멸 관련
            'BCL2', 'BCL2L1', 'MCL1', 'BAX', 'BAK1', 'BID', 'BIM', 'PUMA',
            'NOXA', 'CASP3', 'CASP8', 'CASP9', 'APAF1', 'CYCS',
            
            # 혈관신생 관련
            'VEGFA', 'VEGFB', 'VEGFC', 'FLT1', 'KDR', 'FLT4', 'ANG',
            'ANGPT1', 'ANGPT2', 'TEK', 'PDGFA', 'PDGFB',
            
            # 호르몬 관련
            'ESR1', 'ESR2', 'AR', 'PGR', 'NR3C1', 'THRB', 'VDR',
            'RARA', 'RARB', 'RXRA', 'PPARA', 'PPARG'
        }
        
        self.gene_database.update(common_genes)
        
        # 2. 온라인 유전자 데이터베이스에서 추가 로드 시도
        try:
            self.load_online_gene_database()
        except Exception as e:
            print(f"온라인 데이터베이스 로드 실패: {e}")
            print("기본 유전자 데이터베이스만 사용합니다.")
        
        print(f"총 {len(self.gene_database)}개의 유전자가 로드되었습니다.")
    
    def load_online_gene_database(self):
        """
        온라인에서 추가 유전자 정보를 가져옴
        (HUGO Gene Nomenclature Committee 등의 공개 데이터)
        """
        # HGNC (HUGO Gene Nomenclature Committee) API 사용 예시
        try:
            print("HGNC API에서 유전자 정보를 가져오는 중...")
            
            # 실제 API 호출 예시 (간단한 버전)
            # 주의: 실제 사용시에는 API 제한을 고려해야 함
            url = "https://rest.genenames.org/fetch/status/Approved"
            headers = {'Accept': 'application/json'}
            
            # 실제로는 페이지네이션을 통해 모든 데이터를 가져와야 함
            response = requests.get(url, headers=headers, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                if 'response' in data and 'docs' in data['response']:
                    online_count = 0
                    for gene in data['response']['docs']:
                        if 'symbol' in gene:
                            self.gene_database.add(gene['symbol'].upper())
                            online_count += 1
                        if 'alias_symbol' in gene:
                            for alias in gene['alias_symbol']:
                                self.gene_database.add(alias.upper())
                                online_count += 1
                        if 'prev_symbol' in gene:
                            for prev in gene['prev_symbol']:
                                self.gene_database.add(prev.upper())
                                online_count += 1
                    
                    print(f"HGNC API에서 {online_count}개의 유전자를 추가로 로드했습니다.")
                    
            else:
                print(f"API 요청 실패: HTTP {response.status_code}")
            
        except requests.RequestException as e:
            print(f"온라인 데이터베이스 접근 실패: {e}")
        except Exception as e:
            print(f"데이터 처리 중 오류: {e}")
    
    def extract_text_from_pdf(self, pdf_path: str) -> str:
        """PDF 파일에서 텍스트를 추출"""
        try:
            with open(pdf_path, 'rb') as file:
                reader = PyPDF2.PdfReader(file)
                text = ""
                
                for page_num in range(len(reader.pages)):
                    page = reader.pages[page_num]
                    text += page.extract_text()
                
                return text
        
        except Exception as e:
            print(f"PDF 읽기 오류 ({pdf_path}): {e}")
            return ""
    
    def find_genes_in_text(self, text: str) -> List[str]:
        """
        텍스트에서 유전자 이름을 찾아서 반환
        """
        found_genes = []
        
        # 텍스트를 대문자로 변환하여 대소문자 구분 없이 검색
        text_upper = text.upper()
        
        # 단어 경계를 고려한 정규식으로 검색
        for gene in self.gene_database:
            # 유전자 이름 앞뒤로 단어 경계나 특수문자가 있는지 확인
            pattern = r'\b' + re.escape(gene) + r'\b'
            matches = re.findall(pattern, text_upper)
            
            # 발견된 횟수만큼 리스트에 추가
            found_genes.extend([gene] * len(matches))
        
        return found_genes
    
    def analyze_pdfs(self, pdf_paths: List[str]) -> Dict[str, int]:
        """
        여러 PDF 파일을 분석하여 유전자 빈도를 계산
        """
        all_genes = []
        
        for pdf_path in pdf_paths:
            if not os.path.exists(pdf_path):
                print(f"파일을 찾을 수 없습니다: {pdf_path}")
                continue
            
            print(f"분석 중: {os.path.basename(pdf_path)}")
            
            # PDF에서 텍스트 추출
            text = self.extract_text_from_pdf(pdf_path)
            
            if text:
                # 텍스트에서 유전자 찾기
                genes = self.find_genes_in_text(text)
                all_genes.extend(genes)
                print(f"  -> {len(genes)}개의 유전자 언급 발견")
            else:
                print(f"  -> 텍스트 추출 실패")
        
        # 빈도 계산
        gene_counts = Counter(all_genes)
        return dict(gene_counts)
    
    def get_top_genes(self, pdf_paths: List[str], top_n: int = 50) -> List[Tuple[str, int]]:
        """
        PDF들에서 가장 많이 언급된 유전자들을 순서대로 반환
        
        Args:
            pdf_paths: 분석할 PDF 파일 경로들
            top_n: 상위 몇 개까지 반환할지
            
        Returns:
            (유전자명, 빈도) 튜플의 리스트, 빈도 순으로 정렬됨
        """
        gene_counts = self.analyze_pdfs(pdf_paths)
        
        # 빈도순으로 정렬
        sorted_genes = sorted(gene_counts.items(), key=lambda x: x[1], reverse=True)
        
        return sorted_genes[:top_n]
    
    def save_results_to_file(self, results: List[Tuple[str, int]], filename: str = "gene_analysis_results.txt"):
        """분석 결과를 파일로 저장"""
        try:
            with open(filename, 'w', encoding='utf-8') as f:
                f.write("유전자 빈도 분석 결과\n")
                f.write("=" * 50 + "\n\n")
                f.write(f"{'순위':<5} {'유전자명':<15} {'빈도':<10}\n")
                f.write("-" * 30 + "\n")
                
                for i, (gene, count) in enumerate(results, 1):
                    f.write(f"{i:<5} {gene:<15} {count:<10}\n")
            
            print(f"결과가 {filename}에 저장되었습니다.")
            
        except Exception as e:
            print(f"파일 저장 오류: {e}")


def main():
    """테스트용 메인 함수"""
    analyzer = GeneAnalyzer()
    
    # 예시 사용법
    pdf_files = ["[W2-3] Stochastic pulsing of gene expression enables the generation of spatial patterns in Bacillus subtilis biofilms (Nadezhdin et al., 2020).pdf", "[W1-1] Methods and applications for single-cell and spatial multi-omics (Vandereyken et al., 2023).pdf", "[W1-2]Programming self-organizing multicellular structures with synthetic cell-cell signaling (Toda et al., 2018).pdf", "[W1-3] A guide to machine learning for biologists (Greener et al., 2022).pdf", "[W2-1] Stochasticity and determinism in cell fate decisions (Zechner et al., 2020).pdf", "[W2-2] A stochastic vs deterministic perspective on the timing of cellular events (Ham et al., 2024).pdf"
        ,"[W4-2] Engineering longevity – design of a synthetic gene oscillator to slow cellular aging (Zhou et al., 2023).pdf", "[W4-3] High-resolution and programmable RNA-IN and RNA-OUT genetic circuit in living mammalian cells(Supple).pdf", "[W4-3] High-resolution and programmable RNA-IN and RNA-OUT genetic circuit in living mammalian cells(Zhang, Min, et al., 2024).pdf", "[W3-1] Nonlinear dynamics in phosphoinositide metabolism (Fung et al., 2024).pdf", "[W3-2] The nonlinearity of regulation in biological networks (Manicka et.al., 2023).pdf", "[W3-3] Supplementary file.pdf", "[W3-3] Widespread biochemical reaction networks enable Turing patterns without imposed feedback (Paul et al., 2024).pdf", "[W4-1] Intracellular biosensor-based dynamic regulation to manipulate gene expression at the spatiotemporal level(S. Zhou et al., 2023).pdf"

        # 실제 PDF 파일 경로들을 여기에 입력
        # "path/to/paper1.pdf",
        # "path/to/paper2.pdf"
    ]
    
    if pdf_files:
        print("PDF 분석을 시작합니다...")
        results = analyzer.get_top_genes(pdf_files, top_n=30)
        
        print("\n=== 분석 결과 ===")
        for i, (gene, count) in enumerate(results, 1):
            print(f"{i:2d}. {gene:<12} ({count}번 언급)")
        
        # 결과를 파일로 저장
        analyzer.save_results_to_file(results)
    
    else:
        print("분석할 PDF 파일을 지정해주세요.")
        print("예시:")
        print('pdf_files = ["논문1.pdf", "논문2.pdf"]')


if __name__ == "__main__":
    main()