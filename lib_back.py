# 구동 버전: 3.10.18


# 외부 라이브러리 모음집
from Bio import Entrez                                                           # NCBI 접근용
import numpy as np                                                               #  넘버링 관련
import pandas as pd                                                              # 각종 툴이 잔뜩 들어있음. 분석 관련해서 가장 중요한 라이브러리
from pandas.core.interchange.dataframe_protocol import DataFrame                 # 데이터프레임 형성
from scipy import stats                                                          # 통계 모델
from scipy.stats import norm                                                     # 정규분포 모델
from sklearn.metrics import roc_auc_score                                        # ROC-AUC 계산
import networkx as nx                                                            # 네트워크 분석 모델
import time                                                                      # 시간별 분석에 필요한 거
from collections import defaultdict                                              # 정보 모아서 dictionary화
import matplotlib.pyplot as plt                                                  # 그래프 만드는 툴
from sklearn.linear_model import LinearRegression                                # 선형회귀 모델
from sklearn.preprocessing import PolynomialFeatures                             # 선형회귀 기반 다항회귀 모델
from statsmodels.tsa.statespace.sarimax import SARIMAX                           # 통계 툴, ARIMA 모델
import GEOparse                                                                  # GEO raw data 접근
import calendar                                                                  # 시계열 데이터 다루기
from openai import OpenAI                                                        # 오픈ai
import os



# 1. Data_Access: 데이터베이스 접근용

class Data:
    # 클래스 생성자
    def __init__(self, email):

        self.email = email
        Entrez.email = self.email

    # 유전자 정보 요약 가져오기
    def gene_summary(self, gene_name: str, organism: str):
        try:
            # 1. 유전자 검색
            search_term = f"{gene_name}[Gene Name] AND {organism}[Organism]"
            handle = Entrez.esearch(db="gene", term=search_term)
            record = Entrez.read(handle)
            handle.close()

            if not record["IdList"]:
                print(f"'{gene_name}' 유전자 검색 결과 없음.")
                return

            gene_id = record["IdList"][0]
            print(f"Found Gene ID for {gene_name}: {gene_id}")

            # 2. 유전자 요약 정보 가져오기
            handle = Entrez.esummary(db="gene", id=gene_id)
            summary = Entrez.read(handle)
            handle.close()

            gene_info = summary['DocumentSummarySet']['DocumentSummary'][0]

            print("\n유전자 요약 정보:")
            print(f"Gene Symbol: {gene_info['NomenclatureSymbol']}")
            print(f"Description: {gene_info['Description']}")
            print(f"Chromosome: {gene_info['Chromosome']}")
            print(f"Map Location: {gene_info['MapLocation']}")
            print(f"Summary: {gene_info['Summary']}")

        except Exception as e:
            print(f"오류 발생: {e}")

    # PubMed 기반 논문 수 찾아오기
    def search_pubmed(self, gene_name):

        handle = Entrez.esearch(db="pubmed", term=gene_name, retmode="xml")
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])

# 2. Importance: 중요도 계산 코드

class Importance:
    # 클래스 생성자
    def __init__(self,
                 alpha_diff: float = 0.5,                         # LFC와 p-value Z를 섞는 가중치 alpha (0~1)
                 beta_weights: tuple = (0.5, 0.2, 0.3),           # SCORE1 내 Diff, Rob, AUC의 가중치 (추후 회의를 통해 조정, 합=1)
                 gamma_weights: tuple = (0.25, 0.25, 0.25, 0.25), # SCORE2 내 degree, closeness, betweenness, eigen 가중치 (추후 회의를 통해 조정, 합=1)
                 omega_weights: tuple = (0.7, 0.3),               # 최종 SCORE의 SCORE1, SCORE2 가중치
                 eps: float = 1e-6):                              # 로그 계산시 0 방지용 작은 상수
        # 입력 가중치 인스턴스화
        self.alpha_diff = float(alpha_diff)                        # alpha를 float로 보관
        self.beta_weights = tuple(beta_weights)                    # beta 가중치 튜플로 보관
        self.gamma_weights = tuple(gamma_weights)                  # gamma 가중치 튜플로 보관
        self.omega_weights = tuple(omega_weights)                  # omega 가중치 튜플로 보관
        self.eps = float(eps)                                      # eps 값을 float로 보관

        # 데이터 관련 속성 초기화 (나중에 load_data로 채움)
        self.expr = None                                           # 유전자 발현 DataFrame (samples x genes)
        self.labels = None                                         # 샘플 라벨 Series (index = sample ids)
        self.edges = None                                          # 네트워크 엣지 DataFrame

        # 내부 계산 결과 초기값
        self.df_expr = None                                        # 발현 기반 계산 결과 DataFrame (genes index)
        self.df_net = None                                         # 네트워크 중심성 계산 결과 DataFrame (genes index)
        self.result = None                                         # 최종 병합 결과 DataFrame (genes index)

    @staticmethod
    def _safe_log2(x, eps=1e-6):
        return np.log2(np.clip(x + eps, a_min=eps, a_max=None))    # log2 분모 0 방지

    @staticmethod
    def _min_max_series(s: pd.Series) -> pd.Series:
        mn: float = float(np.nanmin(s))                            # 최솟값
        mx: float = float(np.nanmax(s))                            # 최댓값
        if np.isclose(mx, mn):                                     # 최솟값과 최댓값이 같으면 Series 전체에 0 출력
            return pd.Series(np.zeros_like(s), index=s.index)
        return (s - mn) / (mx - mn)                                # 일반적인 min-max 공식

    @staticmethod
    def _median_abs_deviation(arr: np.ndarray) -> float:
        med: float = float(np.nanmedian(arr))                      # 중앙값 계산
        return float(np.nanmedian(np.abs(arr - med)))              # 중앙값으로부터의 절대편차들 중앙값 반환

    def load_data_from_df(self, expr_df: pd.DataFrame, labels_df: pd.DataFrame, edges_df: pd.DataFrame):
        """
        DataFrame을 내부 속성(self.expr, self.labels, self.edges)에 저장

        expr_df : pd.DataFrame
            샘플 x 유전자 발현량 DataFrame. index = sample, columns = genes.
        labels_df : pd.DataFrame
            sample, label 두 컬럼을 가져야 함. label은 0(정상), 1(암).
        edges_df : pd.DataFrame
            최소 ['source','target'] 컬럼을 포함하는 네트워크 엣지 리스트.
        """
        # expr_df 저장 (샘플 x 유전자 매트릭스)
        self.expr = expr_df.apply(pd.to_numeric, errors='coerce')  # 숫자형으로 변환

        # labels_df에서 라벨 시리즈 추출
        if not {'sample', 'label'}.issubset(labels_df.columns):
            raise ValueError("labels_df의 컬럼 결점 발생")
        self.labels = labels_df.set_index('sample')['label'].reindex(self.expr.index)

        if self.labels.isnull().any():
            missing = self.labels[self.labels.isnull()].index.tolist()
            raise ValueError(f"라벨이 없는 샘플이 존재함: {missing}")

        # edges_df 저장
        if not {'source', 'target'}.issubset(edges_df.columns):
            raise ValueError("edges_df의 컬럼 결점 발생")
        self.edges = edges_df.copy()

    # 발현 기반 점수 계산
    def compute_expression_scores(self):

        # 지역 변수에 편의상 할당
        expr = self.expr                                            # 샘플 x 유전자 DataFrame
        labels = self.labels                                        # 샘플별 라벨 Series
        genes = expr.columns.tolist()                               # 유전자 목록(컬럼명 리스트)

        # 암/정상 샘플 인덱스 분리 (label 값이 1이면 암, 0이면 정상으로 가정)
        cancer_idx = labels[labels == 1].index                      # 암 샘플의 인덱스 리스트
        normal_idx = labels[labels == 0].index                      # 정상 샘플의 인덱스 리스트

        # 클래스가 한쪽뿐이면 계산 불가
        if len(cancer_idx) == 0 or len(normal_idx) == 0:
            raise ValueError("labels.csv에 암세포(1)와 정상세포(0) 데이터가 모두 필요함.")

        # 각 유전자의 암/정상 평균 발현 계산
        mean_cancer = expr.loc[cancer_idx].mean(axis=0)             # series: 유전자별 암 평균
        mean_normal = expr.loc[normal_idx].mean(axis=0)             # series: 유전자별 정상 평균

        # LFC 계산
        LFC = self._safe_log2(mean_cancer, self.eps) - self._safe_log2(mean_normal, self.eps)  # 유전자별 LFC
        absLFC = np.abs(LFC)                                        # LFC 절댓값 (발현 변화의 크기)

        # |LFC|의 표준화 -> Z_LFC (전체 유전자 분포의 평균/표준편차 사용)
        mu_absLFC = np.nanmean(absLFC)                              # |LFC|의 평균
        sigma_absLFC = np.nanstd(absLFC, ddof=0)                    # |LFC|의 표준편차
        if np.isclose(sigma_absLFC, 0):
            Z_LFC = pd.Series(np.zeros_like(absLFC), index=genes)   # 분산이 0이면 모두 0으로
        else:
            Z_LFC = (absLFC - mu_absLFC) / sigma_absLFC             # Z-표준화

        # two-sample t-test (Welch's t-test; equal_var=False)
        pvals = pd.Series(index=genes, dtype=float)                 # p-value를 저장할 Series
        tstats = pd.Series(index=genes, dtype=float)                # t-statistics 저장
        for g in genes:
            x = expr.loc[cancer_idx, g].values                      # 암 샘플에서의 g 유전자 발현 배열
            y = expr.loc[normal_idx, g].values                      # 정상 샘플에서의 g 유전자 발현 배열
            try:
                t, p = stats.ttest_ind(x, y, equal_var=False, nan_policy='omit')  # Welch's t-test
            except Exception:
                t, p = 0.0, 1.0                                    # 예외 발생 시 통계적으로 무의미하다고 처리
            if np.isnan(p):
                p = 1.0                                            # p가 NaN이면 1.0으로 대체 (유의하지 않음)
            tstats[g] = t
            pvals[g] = p

        # p-value를 Z-score로 변환
        pclip = np.clip(pvals.values.astype(float), 1e-300, 1.0)     # 너무 작은 값 방지
        Z_p_vals = -norm.ppf(pclip / 2.0)                                         # scipy.stats.norm.ppf: 역누적분포함수
        Z_p_vals = np.nan_to_num(Z_p_vals, nan=0.0, posinf=0.0, neginf=0.0)       # 무한대/NaN 등이 있으면 0으로 대체
        Z_p = pd.Series(Z_p_vals, index=genes)

        # Diff* 계산
        Diff_star = self.alpha_diff * Z_LFC + (1.0 - self.alpha_diff) * Z_p

        # Diff*를 정규화하여 Diff(0~1) 계산
        Diff = self._min_max_series(Diff_star)

        # 암 샘플에서의 유전자별 MAD 계산 -> Rob
        MAD_cancer = expr.loc[cancer_idx].apply(lambda col: self._median_abs_deviation(col.values), axis=0)
        Rob = self._min_max_series(MAD_cancer)

        # AUC 계산
        y_true = labels.loc[expr.index].values                      # 샘플 순서에 맞춘 라벨 배열
        AUCs = pd.Series(index=genes, dtype=float)                  # AUC 저장용
        for g in genes:
            vals = expr[g].values                                   # 유전자 g의 샘플별 발현값 배열
            try:
                # sklearn.metrics.roc_auc_score는 레이블이 0/1이어야 함
                if len(np.unique(vals[~np.isnan(vals)])) < 2:
                    AUCs[g] = 0.5                                   # 모든 값이 동일하면 분류 불가 -> AUC=0.5로 처리
                else:
                    AUCs[g] = roc_auc_score(y_true, vals)           # AUC 계산
            except Exception:
                AUCs[g] = 0.5                                       # 오류 발생 시 0.5로 대체

        # AUC 범위 정규화
        AUC_absprime = np.abs(2.0 * (AUCs - 0.5))

        # 계산된 모든 지표를 DataFrame으로 묶어서 self.df_expr에 저장 (index = gene names)
        self.df_expr = pd.DataFrame({
            'LFC': LFC,
            'absLFC': absLFC,
            'Z_LFC': Z_LFC,
            't_stat': tstats,
            'pval': pvals,
            'Z_p': Z_p,
            'Diff_star': Diff_star,
            'Diff': Diff,
            'MAD_cancer': MAD_cancer,
            'Rob': Rob,
            'AUC': AUCs,
            'AUC_absprime': AUC_absprime
        })

        # 결과 DataFrame 반환
        return self.df_expr

    # SCORE1 계산 (발현 기반 종합 점수)
    def compute_score1(self):

        b1, b2, b3 = self.beta_weights                                                                         # beta 가중치를 언패킹

        SCORE1_raw = b1 * self.df_expr['Diff'] + b2 * self.df_expr['Rob'] + b3 * self.df_expr['AUC_absprime']  # SCORE1 raw 계산 (정규화 전 가중합)

        self.df_expr['SCORE1_raw'] = SCORE1_raw                                                                # raw 값을 DataFrame에 저장

        self.df_expr['SCORE1'] = self._min_max_series(SCORE1_raw)                                              # 0~1 정규화한 SCORE1 저장

        return self.df_expr

    # 네트워크 중심성 계산
    def compute_network_scores(self):

        # edges DataFrame과 유전자 목록을 편의 변수에 할당
        edges = self.edges
        gene_list = list(self.expr.columns)

        # 그래프 객체 생성(무향 그래프로 가정)
        G = nx.Graph()

        # 엣지 리스트를 순회하면서 그래프에 추가 (source-target이 문자열로 취급되도록 변환)
        for _, row in edges.iterrows():
            s = str(row['source'])
            t = str(row['target'])
            if s == '' or t == '':                                 # 빈 값 예외 처리
                continue
            G.add_edge(s, t)                                       # 네트워크에 엣지 추가

        # 유전자 리스트 안의 유전자들이 그래프에 존재하지 않으면 노드로 추가(고립 노드 허용)
        for g in gene_list:
            if g not in G:
                G.add_node(g)                                     # 고립 노드로 추가

        # 노드 총수 N 계산
        N = G.number_of_nodes()

        # degree 중심성 계산
        if N <= 1:
            deg_centrality = {n: 0.0 for n in G.nodes()}
        else:
            deg_centrality = {n: int(G.degree(n)) / float(N - 1) for n in G.nodes()}

        # closeness 중심성 계산
        closeness_centrality = nx.closeness_centrality(G)

        # betweenness 중심성 계산
        betweenness_centrality = nx.betweenness_centrality(G, normalized=True)

        # eigenvector centrality 계산
        try:
            eigen_centrality = nx.eigenvector_centrality_numpy(G)                         # 행렬분해 기반 계산
        except Exception:
            try:
                eigen_centrality = nx.eigenvector_centrality(G, max_iter=200, tol=1e-06)  # power-iteration fallback
            except Exception:
                eigen_centrality = {n: 0.0 for n in G.nodes()}                            # 모두 실패 시 0으로 대체

        # DataFrame으로 모아서 gene_list 순서로 재정렬(없는 값은 0으로 채움)
        df = pd.DataFrame({
            'degree_node': pd.Series(deg_centrality),
            'closeness': pd.Series(closeness_centrality),
            'betweenness': pd.Series(betweenness_centrality),
            'eigenvector': pd.Series(eigen_centrality)
        }).reindex(gene_list).fillna(0.0)

        # 각 중심성에 대해 0~1 min-max 정규화된 컬럼을 추가
        df['degree_node_norm'] = self._min_max_series(df['degree_node'])
        df['closeness_norm'] = self._min_max_series(df['closeness'])
        df['betweenness_norm'] = self._min_max_series(df['betweenness'])
        df['eigenvector_norm'] = self._min_max_series(df['eigenvector'])

        # 내부 속성에 저장
        self.df_net = df
        return df

    # SCORE2 계산 (네트워크 기반 종합 점수)
    def compute_score2(self):

        g1, g2, g3, g4 = self.gamma_weights                         # gamma 가중치 언패킹
        s_raw = (g1 * self.df_net['degree_node_norm'] +
                 g2 * self.df_net['closeness_norm'] +
                 g3 * self.df_net['betweenness_norm'] +
                 g4 * self.df_net['eigenvector_norm'])
        self.df_net['SCORE2_raw'] = s_raw                           # 가중합 raw
        self.df_net['SCORE2'] = self._min_max_series(s_raw)         # 0~1 정규화된 SCORE2
        return self.df_net

    # 최종 통합 및 출력
    def compute_final_score(self, output_csv: str = "gene_scores_output.csv"):

        # df_expr와 df_net을 유전자 인덱스로 병합 (left join으로 expr 기반 모든 유전자를 포함)
        merged = self.df_expr.merge(
            self.df_net[['degree_node', 'degree_node_norm', 'closeness', 'closeness_norm',
                         'betweenness', 'betweenness_norm', 'eigenvector', 'eigenvector_norm',
                         'SCORE2', 'SCORE2_raw']],
            left_index=True, right_index=True, how='left').fillna(0.0)  # 네트워크에 없는 경우 0으로 채움

        # omega 가중치 적용하여 raw 최종점수 계산
        w1, w2 = self.omega_weights
        merged['SCORE_final_raw'] = w1 * merged['SCORE1_raw'] + w2 * merged['SCORE2_raw']

        # 0~1 정규화
        merged['SCORE_final'] = self._min_max_series(merged['SCORE_final_raw'])

        # SCORE_final을 기준으로 내림차순 정렬 (중요도 높은 유전자 상단)
        merged = merged.sort_values('SCORE_final', ascending=False)

        # CSV로 저장
        merged.to_csv(output_csv)

        # 내부 속성에 저장
        self.result = merged

        # 결과 반환
        return merged

    # 전체 실행용
    def run_pipeline(self, expr_df: DataFrame, labels_df: DataFrame, edges_df: DataFrame, output_csv: str = "gene_scores_output.csv"):

        # 데이터 로드
        self.load_data_from_df(expr_df, labels_df, edges_df)
        # 발현 기반 지표 계산
        self.compute_expression_scores()
        # SCORE1 계산
        self.compute_score1()
        # 네트워크 중심성 계산
        self.compute_network_scores()
        # SCORE2 계산
        self.compute_score2()
        # 최종 통합 및 CSV 출력
        return self.compute_final_score(output_csv)

    # 부분점수 확인 (구체적인 정보가 필요할 경우 이 함수를 켜서 제공)
    def all_results(self) -> pd.DataFrame:

        if self.df_expr is None or self.df_net is None or self.result is None:
            raise ValueError("먼저 파이프라인을 실행하거나 compute_* 메서드를 모두 호출해야 합니다.")

        # SCORE1 부분 점수
        score1_parts = self.df_expr[['Diff', 'Rob', 'AUC_absprime', 'SCORE1']]

        # SCORE2 부분 점수
        score2_parts = self.df_net[['degree_node_norm', 'closeness_norm',
                                    'betweenness_norm', 'eigenvector_norm', 'SCORE2']]

        # 최종 점수
        final_score = self.result[['SCORE_final']]

        # 모두 병합
        merged = score1_parts.join(score2_parts, how='outer').join(final_score, how='outer')

        return merged

# 3-1. 트렌드 분석: 펍메드 논문 수 함수화

class Pub_Analysis:
    def __init__(self, email: str, api_key: str = None):
        Entrez.email = email
        if api_key:
            Entrez.api_key = api_key

    def _search_pubmed(self, keyword: str, start_date: str, end_date: str) -> int:

        query = f"{keyword} AND ({start_date}:{end_date}[PDAT])"
        handle = Entrez.esearch(db="pubmed", term=query, datetype="pdat", rettype="count")
        record = Entrez.read(handle)
        handle.close()
        return int(record["Count"])

    def _fetch_yearly_counts(self, keyword: str, start_year: int, end_year: int) -> dict:
        """
        내부용 함수: 연도별 논문 수
        """
        yearly_counts = defaultdict(int)
        for year in range(start_year, end_year + 1):
            try:
                start_date = f"{year}/01/01"
                end_date = f"{year}/12/31"
                count = self._search_pubmed(keyword, start_date, end_date)
                yearly_counts[year] = count
                time.sleep(0.34)  # API rate 제한 준수
            except Exception as e:
                print(f"Error {year}: {e}")
        return dict(yearly_counts)

    def _fetch_monthly_counts(self, keyword: str, start_year: int, end_year: int) -> dict:
        """
        내부용 함수: 월별 논문 수
        """
        monthly_counts = defaultdict(int)
        for year in range(start_year, end_year + 1):
            for month in range(1, 13):
                start_date = f"{year}/{month:02d}/01"
                last_day = calendar.monthrange(year, month)[1]
                end_date = f"{year}/{month:02d}/{last_day:02d}"

                try:
                    count = self._search_pubmed(keyword, start_date, end_date)
                    monthly_counts[f"{year}-{month:02d}"] = count
                    time.sleep(0.34)
                except Exception as e:
                    print(f"Error {year}-{month:02d}: {e}")
        return dict(monthly_counts)

    # ---------------------------
    # 최종 사용자 함수 (외부 인터페이스)
    # ---------------------------
    def get_yearly_publications(self, keyword: str, start_year: int, end_year: int) -> dict:
        """
        최종 함수: 연도별 논문 수 반환
        """
        return self._fetch_yearly_counts(keyword, start_year, end_year)

    def get_monthly_publications(self, keyword: str, start_year: int, end_year: int) -> dict:
        """
        최종 함수: 월별 논문 수 반환
        """
        return self._fetch_monthly_counts(keyword, start_year, end_year)

# 3-2. 트렌드 분석: 향후 5년 출판 논문 수 예측

class Trend:
    def __init__(self, data: dict, degree: int = 2, freq: str = "year"):
        """
        Trend 분석 클래스 (연도별/월별 모두 지원)
        :param data: {연도: count} 또는 {"YYYY-MM": count} dict
        :param degree: 다항 회귀 차수
        :param freq: "year" 또는 "month"
        """
        self.degree = degree
        self.poly = PolynomialFeatures(degree=self.degree)
        self.model_poly = None
        self.model_arima = None
        self.freq = freq

        # dict → DataFrame 변환
        self.df = pd.DataFrame(list(data.items()), columns=["time", "count"]).sort_values("time")

        if freq == "year":
            # 연도 데이터는 그대로 숫자로 처리
            self.df["time"] = self.df["time"].astype(int)
        elif freq == "month":
            # YYYY-MM 문자열을 datetime으로 변환
            self.df["time"] = pd.to_datetime(self.df["time"], format="%Y-%m")
            # 시계열 모델 편의를 위해 int index 생성 (예: 0,1,2,...)
            self.df["t_index"] = np.arange(len(self.df))
        else:
            raise ValueError("freq는 'year' 또는 'month'만 지원합니다.")

    # ----------------------
    # 다항 회귀
    # ----------------------
    def fit_polynomial(self):
        if self.freq == "year":
            X = self.df[["time"]].values
        else:  # month
            X = self.df[["t_index"]].values

        y = self.df["count"].values
        X_poly = self.poly.fit_transform(X)

        self.model_poly = LinearRegression()
        self.model_poly.fit(X_poly, y)

    def predict_polynomial(self, future_periods: int = 12) -> pd.DataFrame:
        if self.model_poly is None:
            raise RuntimeError("Polynomial 모델이 학습되지 않았습니다.")

        if self.freq == "year":
            last_time = self.df["time"].max()
            future_times = np.arange(last_time + 1, last_time + future_periods + 1).reshape(-1, 1)
            X_poly_future = self.poly.transform(future_times)
            preds = self.model_poly.predict(X_poly_future)
            future_df = pd.DataFrame({"time": future_times.flatten(), "poly_pred": preds})

        else:  # month
            last_index = self.df["t_index"].max()
            future_idx = np.arange(last_index + 1, last_index + future_periods + 1).reshape(-1, 1)
            X_poly_future = self.poly.transform(future_idx)
            preds = self.model_poly.predict(X_poly_future)

            future_dates = pd.date_range(start=self.df["time"].max() + pd.offsets.MonthBegin(1),
                                         periods=future_periods, freq="MS")
            future_df = pd.DataFrame({"time": future_dates, "poly_pred": preds})

        future_df["poly_pred"] = future_df["poly_pred"].clip(lower=0)
        return future_df

    # ----------------------
    # ARIMA
    # ----------------------
    def fit_arima(self, order=(2, 1, 2), seasonal_order=(1, 1, 1, 12)):
        """
        SARIMA 모델 학습 (계절성 포함)
        :param order: (p,d,q)
        :param seasonal_order: (P,D,Q,s) - s=12 → 12개월 주기 계절성
        """
        y = self.df["count"].values
        self.model_arima = SARIMAX(y, order=order, seasonal_order=seasonal_order,
                                   enforce_stationarity=False, enforce_invertibility=False).fit()

    def predict_arima(self, future_periods: int = 12) -> pd.DataFrame:
        if self.model_arima is None:
            raise RuntimeError("ARIMA/SARIMA 모델이 학습되지 않았습니다.")

        preds = self.model_arima.forecast(steps=future_periods)

        if self.freq == "year":
            last_time = self.df["time"].max()
            future_times = np.arange(last_time + 1, last_time + future_periods + 1)
            future_df = pd.DataFrame({"time": future_times, "arima_pred": preds})
        else:  # month
            future_dates = pd.date_range(start=self.df["time"].max() + pd.offsets.MonthBegin(1),
                                         periods=future_periods, freq="MS")
            future_df = pd.DataFrame({"time": future_dates, "arima_pred": preds})

        future_df["arima_pred"] = future_df["arima_pred"].clip(lower=0)
        return future_df

    # ----------------------
    # Ensemble
    # ----------------------
    def ensemble_predictions(self, poly_future: pd.DataFrame, arima_future: pd.DataFrame,
                             w_poly: float = 0.3, w_arima: float = 0.7) -> pd.DataFrame:
        merged = pd.merge(poly_future, arima_future, on="time")
        merged["ensemble_pred"] = merged["poly_pred"] * w_poly + merged["arima_pred"] * w_arima
        merged["ensemble_pred"] = merged["ensemble_pred"].clip(lower=0)
        return merged

    # ----------------------
    # 시각화
    # ----------------------
    def plot_comparison(self, poly_future, arima_future, ensemble_future):
        plt.figure(figsize=(10, 6))

        plt.plot(self.df["time"], self.df["count"], marker="o", color="blue", label="Historical")

        plt.plot(poly_future["time"], poly_future["poly_pred"], "g--", marker="x", label=f"Polynomial (deg={self.degree})")
        plt.plot(arima_future["time"], arima_future["arima_pred"], "r--", marker="s", label="ARIMA")
        plt.plot(ensemble_future["time"], ensemble_future["ensemble_pred"], "k-", marker="d", linewidth=2,
                 label="Ensemble (0.3 Poly + 0.7 ARIMA)")

        plt.xlabel("Time")
        plt.ylabel("Number of Publications")
        plt.title("Publication Trend Analysis")
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

# raw data 가져오기 (미완료)
class GeoDataExtractor:
    def __init__(self, geo_id: str):
        self.geo_id = geo_id
        self.gse = None
        self.platform = None

    def load_geo(self):
        self.gse = GEOparse.get_GEO(geo=self.geo_id)
        self.platform = list(self.gse.gpls.values())[0]

    def extract_labels(self) -> pd.DataFrame:
        parsers = {
            "GSE15852": self._parse_labels_gse15852,
            "GSE65212": self._parse_labels_gse65212,
            "GSE152908": self._parse_labels_gse152908,
        }
        if self.geo_id in parsers:
            return parsers[self.geo_id]()
        else:
            raise ValueError(f"{self.geo_id}에 대한 맞춤형 라벨 파서가 정의되지 않았습니다.")

    def _parse_labels_gse15852(self):
        """GSE15852: paired tumor vs normal breast tissue."""
        metadata = []
        for gsm_name, gsm in self.gse.gsms.items():
            val = " ".join(sum(gsm.metadata.values(), [])).lower()
            if "tumor" in val or "tumour" in val:
                status = 1
            elif "normal" in val:
                status = 0
            else:
                status = None
            metadata.append({"sample": gsm_name, "cancer_status": status})
        return pd.DataFrame(metadata).set_index("sample").T

    def _parse_labels_gse65212(self):
        """GSE65212: tumor vs adjacent noncancerous vs complete normal."""
        metadata = []
        for gsm_name, gsm in self.gse.gsms.items():
            val = " ".join(sum(gsm.metadata.values(), [])).lower()
            if "tumor" in val or "cancer" in val:
                status = 1
            elif "adjacent" in val or "noncancerous" in val or "normal" in val:
                status = 0
            else:
                status = None
            metadata.append({"sample": gsm_name, "cancer_status": status})
        return pd.DataFrame(metadata).set_index("sample").T

    def _parse_labels_gse152908(self):
        """GSE152908: includes cell lines, primary tumors, adjacent normal, and healthy tissue."""
        metadata = []
        for gsm_name, gsm in self.gse.gsms.items():
            val = " ".join(sum(gsm.metadata.values(), [])).lower()
            if "tumor" in val or "cancer" in val or "tumour" in val:
                status = 1
            elif "adjacent" in val or "uninvolved" in val or "normal" in val:
                status = 0
            else:
                status = None
            metadata.append({"sample": gsm_name, "cancer_status": status})
        return pd.DataFrame(metadata).set_index("sample").T

    def extract_expression(self, gene_list: list[str]) -> pd.DataFrame:
        expression_data = []
        for gene in gene_list:
            probe_ids = (
                self.platform.table[self.platform.table["Gene Symbol"] == gene]["ID"].tolist()
            )
            if not probe_ids:
                continue
            gene_values = {}
            for gsm_name, gsm in self.gse.gsms.items():
                vals = gsm.table[gsm.table["ID_REF"].isin(probe_ids)]["VALUE"].values
                if len(vals) > 0:
                    gene_values[gsm_name] = float(vals[0])
            if gene_values:
                expression_data.append(pd.Series(gene_values, name=gene))
        return pd.DataFrame(expression_data)

# 논문 개수 기울기 분석
class Slope:
    def __init__(self, data: dict):
        """
        :param data: {연도: 논문 수} 또는 {"YYYY-MM": count} dict
        """
        self.df = pd.DataFrame(list(data.items()), columns=["time", "count"]).sort_values("time")

    def compute_rate_of_change(self) -> pd.DataFrame:
        """
        변화율(%) 계산
        """
        self.df["rate_change"] = self.df["count"].pct_change() * 100
        return self.df

    def detect_significant_changes(self, threshold: float = 2.0) -> tuple:
        """
        유의미한 변화 시점 탐지 (평균 ± threshold*표준편차)
        :param threshold: 표준편차 배수 (default=2.0 → 95% 수준)
        :return: (연도1, 연도2, ...) 형태의 튜플
        """
        if "rate_change" not in self.df.columns:
            self.compute_rate_of_change()

        mean_change = self.df["rate_change"].mean(skipna=True)
        std_change = self.df["rate_change"].std(skipna=True)

        lower_bound = mean_change - threshold * std_change
        upper_bound = mean_change + threshold * std_change

        significant_years = self.df[
            (self.df["rate_change"] < lower_bound) |
            (self.df["rate_change"] > upper_bound)
        ]["time"].tolist()

        return tuple(significant_years)

# ChatGPT를 활용해 분석
class HelpGPT:
    def __init__(self, model: str = "gpt-5"):
        """
        OpenAI ChatGPT API를 이용해 원인(reason) 추론
        :param model: 사용할 모델 이름 (기본 gpt-5)
        """
        api_key = os.getenv("OPENAI_API_KEY")  # 환경변수에서 API 키 읽기
        if not api_key:
            raise ValueError("환경변수 OPENAI_API_KEY가 설정되지 않았습니다. "
                             "터미널/환경설정에서 API 키를 등록하세요.")

        self.client = OpenAI(api_key=api_key)
        self.model = model

    def analyze_changes(self, slope_df: pd.DataFrame, significant_years: tuple, gene: str) -> pd.DataFrame:
        """
        PubMed 변화율을 바탕으로 GPT에게 원인(reason) 추론 요청
        :param slope_df: Slope 클래스 결과 DataFrame (time, count, rate_change 포함)
        :param significant_years: 유의미한 변화 연도 리스트/튜플
        :param gene: 분석 대상 유전자명
        :return: DataFrame(year, rate_change, reason)
        """
        results = []

        for year in significant_years:
            row = slope_df[slope_df["time"] == year].iloc[0]
            rate_change = row["rate_change"]

            # 전후년도 출판량
            prev_count = slope_df[slope_df["time"] == year - 1]["count"].values
            next_count = slope_df[slope_df["time"] == year + 1]["count"].values

            prev_count = int(prev_count[0]) if len(prev_count) > 0 else None
            next_count = int(next_count[0]) if len(next_count) > 0 else None

            # GPT 프롬프트
            prompt = f"""
            유전자: {gene}
            연도: {year}
            전년도 논문 수: {prev_count}
            해당년도 논문 수: {row['count']}
            다음년도 논문 수: {next_count}
            변화율: {rate_change:.2f}%

            위 데이터를 바탕으로, {year}년에 논문 수 변화가 발생한 학술적/사회적/기술적 배경을 추론해줘.
            가능한 원인을 간단히 요약해 reason으로 제시해줘.
            """

            response = self.client.chat.completions.create(
                model=self.model,
                messages=[{"role": "user", "content": prompt}],
                temperature=0.3
            )

            reason = response.choices[0].message.content.strip()
            results.append({"year": year, "rate_change": rate_change, "reason": reason})

        return pd.DataFrame(results, columns=["year", "rate_change", "reason"])