"""
lib_front_page3.py
GeNetwork - Frontend Page 3 (PySimpleGUI)

기능:
 1. GRN 네트워크 (좌측 상단, 큰 영역) - 임시 네트워크 시각화
 2. Clusters of gene nodes (좌측 하단, 작은 리스트)
 3. Expression-based importance score (우측 상단, 표)
 4. Hot nodes & interactions (우측 중간, 빈 칸)
 5. Trend of gene publications (우측 하단, 빈 칸)
"""

import PySimpleGUI as sg
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import pathlib

from lib_back import Importance

# ---------- Helper: Matplotlib Figure를 PySimpleGUI에 표시 ----------
def draw_figure(canvas, figure):
    fig_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    fig_canvas_agg.draw()
    fig_canvas_agg.get_tk_widget().pack(side="top", fill="both", expand=1)
    return fig_canvas_agg

# ---------- Page 3 실행 ----------
def run_page3(gene_list):
    sg.theme("LightBlue2")

    # 전체 화면 크기
    screen_w, screen_h = sg.Window.get_screen_size()
    win_w, win_h = int(screen_w * 0.95), int(screen_h * 0.95)

    # GRN of gene nodes (좌측 상단, 큰 영역)
    grn_col = [
        [sg.Text("GRN of gene nodes", font=("Helvetica", 12, "bold"))],
        [sg.Canvas(key="-GRN-CANVAS-", size=(int(win_w*0.9), int(win_h*0.7)))]
    ]

    # Clusters of gene nodes (좌측 하단)
    cluster_col = [
        [sg.Text("Clusters of gene nodes", font=("Helvetica", 12, "bold"))],
        [sg.Listbox(values=[], size=(120,8), key="-CLUSTERS-")]
    ]

    # Expression-based importance score (우측 상단)
    expr_col = [
        [sg.Text("Expression-based importance score", font=("Helvetica", 12, "bold"))],
        [sg.Table(values=[],
                  headings=["Gene", "Score1", "Diff", "Rob", "AUC"],
                  auto_size_columns=False,
                  col_widths=[15, 10, 10, 10, 10],
                  justification="right",
                  num_rows=10,
                  key="-TABLE-")]
    ]

    # Hot nodes (우측 중간)
    hot_col = [
        [sg.Text("Hot nodes and important interactions", font=("Helvetica", 12, "bold"))],
        [sg.Multiline("", size=(100,8), key="-HOT-", disabled=True)]
    ]

    # Trend of publications (우측 하단)
    trend_col = [
        [sg.Text("Trend of gene publications", font=("Helvetica", 12, "bold"))],
        [sg.Canvas(key="-TREND-CANVAS-", size=(int(win_w*0.8), int(win_h*0.5)))]
    ]

    # 레이아웃: 좌측(상단: GRN, 하단: Clusters), 우측(상단: Expr, 중간: Hot, 하단: Trend)
    layout = [
        [sg.Column([[sg.Column(grn_col)], [sg.Column(cluster_col)]], size=(int(win_w*0.65), win_h)),
         sg.VSeperator(),
         sg.Column([[sg.Column(expr_col)],
                    [sg.Column(hot_col)],
                    [sg.Column(trend_col)]],
                   size=(int(win_w*0.3), win_h))]
    ]

    window = sg.Window("GeNetwork - Page 3", layout, finalize=True, resizable=True, size=(win_w, win_h))

    # ----- 1. GRN 네트워크 표시 (Arbitrary code) -----
    try:
        fig_grn = plt.figure(figsize=(6,5))
        G = nx.Graph()
        G.add_nodes_from(gene_list)
        # [Arbitrary code] 랜덤 edge 생성
        for i in range(len(gene_list)-1):
            G.add_edge(gene_list[i], gene_list[i+1])
        pos = nx.spring_layout(G)
        nx.draw(G, pos, with_labels=True, node_size=500, font_size=8)
        draw_figure(window["-GRN-CANVAS-"].TKCanvas, fig_grn)
    except Exception as e:
        print("GRN drawing failed:", e)

    # ----- 2. Expression-based importance score -----
    try:
        data_path = pathlib.Path(r"C:\Users\SAMSUNG\PycharmProjects\GenetWork\data\imaginary_data.csv")
        if data_path.exists():
            df = pd.read_csv(data_path, index_col=0)
            df = pd.read_csv(data_path)
            expr = df.drop(columns=["Sample", "Labels"]).set_index(df["Sample"])
            labels = df[["Sample", "Labels"]].rename(columns={"Sample": "sample", "Labels": "label"})

            if labels is not None:
                imp = Importance()
                imp.load_data_from_df(expr, labels)
                imp.compute_expression_scores()
                result = imp.compute_score1()

                # gene_list 필터링
                filtered = result.loc[result.index.intersection(gene_list)]
                table_data = [[g, round(filtered.loc[g,"SCORE1"],3),
                               round(filtered.loc[g,"LFC"],3),
                               round(filtered.loc[g,"AUC"],3)] for g in filtered.index]
                window["-TABLE-"].update(values=table_data)
    except Exception as e:
        print("Importance score failed:", e)

    # ----- Clusters / Hot nodes / Trend: 현재는 빈 칸 -----

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break

    window.close()
