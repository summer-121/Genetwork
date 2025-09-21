#!/usr/bin/env python3
"""
run_importance_expression_only.py

Importance 클래스 (발현 기반 점수만) 실행 예제
임의의 작은 raw data를 포함하여 시험 실행 가능
"""

import pandas as pd
from lib_back import Importance   # 앞서 만든 클래스 파일 import

def main():
    # --------------------------
    # 1. 시험용 toy raw data 준비
    # --------------------------
    labels = (1, 0, 1, 0)  # 암(1), 정상(0)
    Gene1  = (6.6, 3.2, 5.6, 2.8)
    Gene2  = (2.2, 3.3, 3.0, 2.9)
    Gene3  = (1.0, 4.0, 1.2, 3.8)
    Gene4  = (3.3, 3.2, 3.1, 3.2)

    samples = [f"S{i+1}" for i in range(len(labels))]

    # expr_df 생성 (샘플 × 유전자)
    expr_df = pd.DataFrame({
        "Gene1": Gene1,
        "Gene2": Gene2,
        "Gene3": Gene3,
        "Gene4": Gene4
    }, index=samples)

    # labels_df 생성
    labels_df = pd.DataFrame({
        "sample": samples,
        "label": labels
    })

    # --------------------------
    # 2. Importance 실행
    # --------------------------
    imp = Importance(alpha_diff=0.5, beta_weights=(0.5, 0.2, 0.3))

    # 데이터 로드
    imp.load_data_from_df(expr_df, labels_df)

    # 발현 기반 점수 계산
    imp.compute_expression_scores()
    result = imp.compute_score1()

    # --------------------------
    # 3. 결과 출력
    # --------------------------
    print("=== 상위 유전자 중요도 (SCORE1 기반) ===")
    print(result[["SCORE1", "Diff", "Rob", "AUC_absprime"]].head())



if __name__ == "__main__":
    main()
