# 구동 버전: 3.10.18
from __future__ import annotations

# 외부 라이브러리 모음집
from Bio import Entrez                                                           # NCBI 접근용
import numpy as np                                                               #  넘버링 관련
import pandas as pd                                                              # 각종 툴이 잔뜩 들어있음. 분석 관련해서 가장 중요한 라이브러리
from scipy import stats                                                          # 통계 모델
from scipy.stats import norm                                                     # 정규분포 모델
from sklearn.metrics import roc_auc_score                                        # ROC-AUC 계산
import time                                                                      # 시간별 분석에 필요한 거
from collections import Counter, defaultdict                                     # 정보 모아서 dictionary화
import matplotlib.pyplot as plt                                                  # 그래프 만드는 툴
from sklearn.linear_model import LinearRegression                                # 선형회귀 모델
from sklearn.preprocessing import PolynomialFeatures                             # 선형회귀 기반 다항회귀 모델
from statsmodels.tsa.statespace.sarimax import SARIMAX                           # 통계 툴, ARIMA 모델
import GEOparse                                                                  # GEO raw data 접근
import calendar                                                                  # 시계열 데이터 다루기
from openai import OpenAI                                                        # 오픈ai
import os
import pathlib
import re
from typing import Dict, List
import fitz    #PyMuPDF의 namescape
import spacy



# 데이터베이스 접근용 (유전자명 자동완성)
class Data:
    def __init__(self, email = "1018jjkk@gmail.com"):
        """
        데이터 접근 클래스 (NCBI Gene + PubMed)
        """
        self.email = email
        Entrez.email = self.email

    def search_gene(self, gene_name):
        """
        NCBI Gene DB에서 특정 유전자를 검색해
        Gene ID와 Gene 이름을 리스트로 반환

        :param gene_name: 유전자 이름 (예: "UvrA")
        :return: [(gene_id, gene_name), ...] 형태의 리스트
        """
        handle = Entrez.esearch(db="gene", term=gene_name, retmode="xml")
        record = Entrez.read(handle)
        handle.close()

        gene_list = []
        for gid in record["IdList"]:
            # Gene 정보 가져오기
            h = Entrez.efetch(db="gene", id=gid, retmode="xml")
            gene_record = Entrez.read(h)
            h.close()

            try:
                official_name = gene_record[0]["Entrezgene_gene"]["Gene-ref"]["Gene-ref_locus"]
            except (KeyError, IndexError):
                official_name = "Unknown"

            gene_list.append((gid, official_name))

        return gene_list

# 논문 텍스트 마이닝 클래스
class Text:
    def __init__(
        self,
        model_name: str = "en_ner_jnlpba_md",
        hgnc_path: pathlib.Path | None = None,
        whitelist_path: pathlib.Path | None = None,
        ref_weight: float = 0.5,
        n_process: int = 1,
        batch_size: int = 1000,
    ):
        # Load SciSpaCy model
        try:
            self.nlp = spacy.load(model_name)
        except Exception as e:
            raise RuntimeError(
                f"Failed to load model {model_name}. Ensure SciSpaCy is installed. Reason: {e}"
            )

        # HGNC mapping and whitelist
        self.hgnc_map = self.load_hgnc(hgnc_path)
        self.whitelist = self.load_whitelist(whitelist_path)
        if not self.whitelist and self.hgnc_map:
            # auto whitelist from HGNC
            self.whitelist = set(self.hgnc_map.values())

        self.ref_weight = ref_weight
        self.n_process = max(1, n_process)
        self.batch_size = max(1, batch_size)

        # Valid labels from SciSpaCy
        self.VALID_LABELS = {"gene", "protein", "gene_or_gene_product"}

    # ----------------------------------------------------------------------------------
    # PDF → text
    # ----------------------------------------------------------------------------------

    def extract_pdf_pages(self, pdf_path: pathlib.Path) -> List[str]:
        """Return a list of plain-text strings, one per page."""
        pages: List[str] = []
        try:
            with fitz.open(pdf_path) as doc:
                for page in doc:
                    txt = page.get_text("text") or ""
                    pages.append(txt)
        except Exception as e:
            raise RuntimeError(f"{pdf_path.name}: {e}")
        return pages

    def normalize_page_text(self, text: str) -> str:
        """Light normalization for NER friendliness."""
        text = re.sub(r"(\w)-\n(\w)", r"\1\2", text)  # join hyphenation
        text = re.sub(r"\n{2,}", "\n\n", text)        # keep paragraph breaks
        text = text.replace("\n", " ")
        text = re.sub(r"\s{2,}", " ", text).strip()
        return text

    # ----------------------------------------------------------------------------------
    # Lexicon helpers (HGNC + whitelist)
    # ----------------------------------------------------------------------------------

    def load_hgnc(self, tsv: pathlib.Path | None) -> Dict[str, str]:
        """Load HGNC TSV and build alias→canonical mapping."""
        mapping: Dict[str, str] = {}
        if not tsv or not tsv.exists():
            return mapping

        import csv as _csv
        with tsv.open(encoding="utf-8") as fh:
            reader = _csv.DictReader(fh, delimiter="\t")
            required = {"symbol", "alias_symbol", "prev_symbol"}
            missing = [c for c in required if c not in reader.fieldnames]
            if missing:
                print(f"[WARN] HGNC TSV missing columns: {missing}")
            for row in reader:
                symbol = (row.get("symbol") or "").upper()
                if not symbol:
                    continue
                mapping[symbol] = symbol
                for col in ("alias_symbol", "prev_symbol"):
                    val = row.get(col) or ""
                    for token in (val.split("|") if val else []):
                        token = token.strip().upper()
                        if token:
                            mapping[token] = symbol
        return mapping

    def load_whitelist(self, path: pathlib.Path | None) -> set[str]:
        """Load whitelist of valid gene symbols (UPPERCASE)."""
        if not path or not path.exists():
            return set()
        syms = set()
        with path.open(encoding="utf-8") as fh:
            for line in fh:
                s = line.strip()
                if not s or s.startswith("#"):
                    continue
                syms.add(s.upper())
        return syms

    # ----------------------------------------------------------------------------------
    # Entity counting
    # ----------------------------------------------------------------------------------

    def count_entities_in_pages(self, pages: List[str]) -> Counter:
        """Run NER per page with optional down-weighting after 'References'."""
        cleaned: List[str] = [self.normalize_page_text(p) for p in pages]

        # per-page weights
        weights: List[float] = []
        seen_refs = False
        for p in pages:
            if not seen_refs and re.search(r"\bReferences\b", p, flags=re.IGNORECASE):
                seen_refs = True
            weights.append(self.ref_weight if seen_refs else 1.0)

        counts: Counter = Counter()
        for doc, w in zip(
            self.nlp.pipe(
                cleaned,
                n_process=self.n_process,
                batch_size=self.batch_size,
                disable=["parser", "attribute_ruler", "lemmatizer"],
            ),
            weights,
        ):
            for ent in doc.ents:
                if ent.label_.lower() in self.VALID_LABELS:
                    sym = ent.text.strip().upper()
                    sym = re.sub(r"\s+", " ", sym)
                    sym = sym.strip(";,:.()[]{}")
                    # HGNC normalization
                    sym = self.hgnc_map.get(sym, sym)
                    # whitelist filter
                    if self.whitelist and sym not in self.whitelist:
                        continue
                    counts[sym] += w
        return counts

    # ----------------------------------------------------------------------------------
    # Main API
    # ----------------------------------------------------------------------------------

    def extract_genes_from_pdf(self, pdf_path: pathlib.Path) -> List[str]:
        """Given a PDF path, return list of extracted gene/protein names."""
        pages = self.extract_pdf_pages(pdf_path)
        if not any(p.strip() for p in pages):
            raise RuntimeError(f"No text layer found: {pdf_path.name} — skipped")
        counts = self.count_entities_in_pages(pages)
        return list(counts.keys())


# TODO: GEO로부터 실제 raw data 가져오기 (미완료. 우선 가상데이터 이용하고, 다른 기능들이 완성된 후 작업 예정.)

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

# 중요도 계산 코드
class Importance:
    def __init__(self,
                 alpha_diff=0.5,
                 beta_weights=(0.5, 0.2, 0.3),
                 eps=1e-6):
        # 가중치
        self.alpha_diff = alpha_diff
        self.beta_weights = beta_weights
        self.eps = eps

        # 데이터
        self.expr = None
        self.labels = None

        # 결과
        self.df_expr = None
        self.result = None

    # ------------------- 유틸 -------------------
    @staticmethod
    def safe_log2(x, eps=1e-6):
        return np.log2(np.clip(x + eps, a_min=eps, a_max=None))

    @staticmethod
    def min_max_series(s):
        mn = np.nanmin(s)
        mx = np.nanmax(s)
        if np.isclose(mx, mn):
            return pd.Series(np.zeros_like(s), index=s.index)
        return (s - mn) / (mx - mn)

    @staticmethod
    def median_abs_deviation(arr):
        med = np.nanmedian(arr)
        return np.nanmedian(np.abs(arr - med))

    # ------------------- 데이터 로딩 -------------------
    def load_data_from_df(self, expr_df: pd.DataFrame, labels_df: pd.DataFrame):
        """
        DataFrame을 직접 입력받아 expr, labels를 세팅합니다.
        labels_df는 반드시 ['sample','label'] 컬럼을 가져야 합니다.
        """
        self.expr = expr_df.apply(pd.to_numeric, errors='coerce')  # 숫자형 변환
        if not {'sample','label'}.issubset(labels_df.columns):
            raise ValueError("labels_df는 반드시 'sample'과 'label' 컬럼을 포함해야 합니다.")
        self.labels = labels_df.set_index('sample')['label'].reindex(self.expr.index)
        if self.labels.isnull().any():
            missing = self.labels[self.labels.isnull()].index.tolist()
            raise ValueError(f"labels_df에 없는 샘플이 expr_df에 있습니다: {missing}")

    # ------------------- 발현 기반 점수 -------------------
    def compute_expression_scores(self):
        expr, labels = self.expr, self.labels
        genes = expr.columns
        cancer_idx = labels[labels == 1].index
        normal_idx = labels[labels == 0].index

        # 평균 기반 Fold Change
        mean_cancer = expr.loc[cancer_idx].mean(axis=0)
        mean_normal = expr.loc[normal_idx].mean(axis=0)
        LFC = self.safe_log2(mean_cancer, self.eps) - self.safe_log2(mean_normal, self.eps)
        absLFC = np.abs(LFC)

        # 표준화된 Z-score
        mu_absLFC = np.nanmean(absLFC)
        sigma_absLFC = np.nanstd(absLFC, ddof=0)
        Z_LFC = (absLFC - mu_absLFC) / (sigma_absLFC if sigma_absLFC != 0 else 1)

        # t-test
        pvals, tstats = {}, {}
        for g in genes:
            x = expr.loc[cancer_idx, g].values
            y = expr.loc[normal_idx, g].values
            t, p = stats.ttest_ind(x, y, equal_var=False, nan_policy="omit")
            if np.isnan(p):
                p = 1.0
            pvals[g], tstats[g] = p, t
        pvals = pd.Series(pvals)
        tstats = pd.Series(tstats)

        # Z_p 변환
        pclip = np.clip(pvals.values.astype(float), 1e-300, 1.0)
        Z_p = -norm.ppf(pclip / 2.0)
        Z_p = pd.Series(np.nan_to_num(Z_p, nan=0.0), index=genes)

        # Diff*
        Diff_star = self.alpha_diff * Z_LFC + (1 - self.alpha_diff) * Z_p
        Diff = self.min_max_series(Diff_star)

        # Robustness (MAD in cancer group)
        MAD_cancer = expr.loc[cancer_idx].apply(lambda col: self.median_abs_deviation(col.values), axis=0)
        Rob = self.min_max_series(MAD_cancer)

        # AUC 기반 판별력
        y_true = labels.loc[expr.index].values
        AUCs = {}
        for g in genes:
            vals = expr[g].values
            try:
                AUCs[g] = roc_auc_score(y_true, vals)
            except Exception:
                AUCs[g] = 0.5
        AUCs = pd.Series(AUCs)
        AUC_absprime = np.abs(2 * (AUCs - 0.5))

        self.df_expr = pd.DataFrame({
            "LFC": LFC,
            "absLFC": absLFC,
            "Z_LFC": Z_LFC,
            "t_stat": tstats,
            "pval": pvals,
            "Z_p": Z_p,
            "Diff_star": Diff_star,
            "Diff": Diff,
            "MAD_cancer": MAD_cancer,
            "Rob": Rob,
            "AUC": AUCs,
            "AUC_absprime": AUC_absprime
        })
        return self.df_expr

    def compute_score1(self):
        b1, b2, b3 = self.beta_weights
        SCORE1 = (b1 * self.df_expr["Diff"] +
                  b2 * self.df_expr["Rob"] +
                  b3 * self.df_expr["AUC_absprime"])
        self.df_expr["SCORE1_raw"] = SCORE1
        self.df_expr["SCORE1"] = self.min_max_series(SCORE1)

        # 최종 결과는 SCORE1 정렬
        self.result = self.df_expr.sort_values("SCORE1", ascending=False)
        return self.result


# 트렌드 분석: 펍메드 논문 수 함수화

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

# 트렌드 분석: 향후 5년 출판 논문 수 예측

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


# 논문 개수 기울기 분석
class Slope:
    def __init__(self, data: dict):
        """
        :param data: {연도: 논문 수} 형태의 dict
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
    
# 유의미한 논문 수 변화에 대한 인덱스 부여
class Index:
    def __init__(self, data_dict: dict):
        """
        :param data_dict: {gene: {연도: count, ...}, ...} 형태의 딕셔너리
                          예: {"p53": {2018: 451, 2019: 467}, "BRCA1": {...}}
        """
        self.data_dict = data_dict

    def assign_index(self, threshold: float = 2.0) -> pd.DataFrame:
        """
        각 유전자에 대해 변화율 기반 인덱싱:
        0 = 유의미한 변화 없음
        1 = 유의미한 증가만 존재
        2 = 유의미한 감소만 존재
        3 = 유의미한 증가와 감소 모두 존재
        :param threshold: 유의미성 기준 (표준편차 배수)
        :return: DataFrame(gene, index, significant_years)
        """
        results = []

        for gene, yearly_data in self.data_dict.items():
            slope = Slope(yearly_data)
            slope_df = slope.compute_rate_of_change()
            significant_years = slope.detect_significant_changes(threshold=threshold)

            if len(significant_years) == 0:
                idx = 0
            else:
                changes = slope_df[slope_df["time"].isin(significant_years)]["rate_change"].dropna()
                has_increase = any(changes > 0)
                has_decrease = any(changes < 0)

                if has_increase and has_decrease:
                    idx = 3
                elif has_increase:
                    idx = 1
                elif has_decrease:
                    idx = 2
                else:
                    idx = 0

            results.append({
                "gene": gene,
                "index": idx,
                "significant_years": ", ".join(map(str, significant_years)) if significant_years else "-"
            })

        return pd.DataFrame(results, columns=["gene", "index", "significant_years"])

