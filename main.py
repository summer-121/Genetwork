
from library import Data
from library import Importance
from library import Pub_Analysis
from library import Trend
import pandas as pd

# 출력 옵션 변경 (행, 열 모두 제한 해제)
pd.set_option("display.max_rows", None)   # 모든 행 출력
pd.set_option("display.max_columns", None) # 모든 열 출력
pd.set_option("display.width", None)      # 줄바꿈 없이 한 줄에 표시
pd.set_option("display.max_colwidth", None) # 컬럼 내 문자열 길이 제한 해제


# 이메일
email = "1018jjkk@gmail.com"
data = Data(email)

# 찾는 유전자, 종 이름
gene = "p53"
org = "Homo sapiens"

data.gene_summary(gene, org)

pubmed_count = data.search_pubmed("UvrA")
print(f"PubMed 논문 수: {pubmed_count}")


# 원본 튜플 데이터 (이제 이걸 raw data에서 가져오는 방법을 모색해야 함)
labels = (1, 0, 1, 0)
Gene1  = (6.6, 3.2, 5.6, 2.8)
Gene2  = (2.2, 3.3, 3.0, 2.9)
Gene3  = (3.3, 3.2, 3.1, 3.2)
Gene4  = (1.0, 4.0, 1.2, 3.8)
Gene5 = (0.0, 3.0, 1.0, 3.0)
edges_tuples = [
    ('Gene1', 'Gene5'),
    ('Gene5', 'Gene4'),
    ('Gene2', 'Gene5'),
    ('Gene1', 'Gene4'),
    ('Gene2', 'Gene1'),
    ('Gene2', 'Gene3')
]

# 샘플 ID 생성
samples = [f"S{i+1}" for i in range(len(labels))]

# expr_df 생성 (샘플 x 유전자)
expr_df = pd.DataFrame({
    "Gene1": Gene1,
    "Gene2": Gene2,
    "Gene3": Gene3,
    "Gene4": Gene4,
    "Gene5": Gene5
}, index=samples)

# labels_df 생성 (sample, label)
labels_df = pd.DataFrame({
    "sample": samples,
    "label": labels
})

# edges_df 생성
edges_df = pd.DataFrame(edges_tuples, columns=["source", "target"])


# 1. Importance 인스턴스 생성
imp = Importance(
    alpha_diff=0.5,                       # LFC vs p-value 비중 (Diff* 계산 시)
    beta_weights=(0.5, 0.2, 0.3),         # SCORE1 가중치 (Diff, Rob, AUC)
    gamma_weights=(0.25, 0.25, 0.25, 0.25), # SCORE2 가중치 (degree, closeness, betweenness, eigen)
    omega_weights=(0.7, 0.3),             # 최종 SCORE 가중치 (SCORE1 vs SCORE2)
    eps=1e-6                              # log 계산 시 작은 상수
)

# 2. DataFrame 직접 로딩
imp.load_data_from_df(expr_df, labels_df, edges_df)

# 3. 발현 기반 지표 계산
imp.compute_expression_scores()
imp.compute_score1()

# 4. 네트워크 중심성 계산
imp.compute_network_scores()
imp.compute_score2()

# 5. 최종 점수 계산 및 결과 DataFrame 생성
result_df = imp.compute_final_score(output_csv="gene_scores_output.csv")

# 6. 상위 10개 결과 출력
print("상위 10개 유전자 중요도:")
print(result_df[['SCORE_final', 'SCORE1', 'SCORE2']].head(10))

# 7. 구체적인 부분점수 확인 (구체적인 정보가 필요할 경우 이 함수를 켜서 제공)
all_results = imp.all_results()
print(all_results)

# Step 1: PubMed 데이터 수집 (p53 유전자)
analyzer = Pub_Analysis(email)
yearly_data = analyzer.get_yearly_publications(gene, 2000, 2025)
print("연도별 논문 수:", yearly_data)

# Step 2: 트렌드 분석
trend = Trend(yearly_data, degree=3)

# Polynomial
trend.fit_polynomial()
poly_future = trend.predict_polynomial(future_years=5)

# ARIMA
trend.fit_arima(order=(2,1,2))
arima_future = trend.predict_arima(future_years=5)

# Ensemble
ensemble_future = trend.ensemble_predictions(poly_future, arima_future, w_poly=0.3, w_arima=0.7)

# 결과 출력
print("\nPolynomial 예측 결과:")
print(poly_future)
print("\nARIMA 예측 결과:")
print(arima_future)
print("\nEnsemble 예측 결과 (0.3 Poly + 0.7 ARIMA):")
print(ensemble_future)

# 시각화
trend.plot_comparison(poly_future, arima_future, ensemble_future)