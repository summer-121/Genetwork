"""
lib_front_page3.py
GeNetwork - Frontend Page 3
"""

import PySimpleGUI as sg
import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import pandas as pd
import pathlib

from lib_back import Importance, Pub_Analysis, Trend

def draw_figure(canvas, figure):
    fig_canvas_agg = FigureCanvasTkAgg(figure, canvas)
    fig_canvas_agg.draw()
    fig_canvas_agg.get_tk_widget().pack(side="top", fill="both", expand=1)
    return fig_canvas_agg

def run_page3(gene_list):
    sg.theme("LightBlue2")

    screen_w, screen_h = sg.Window.get_screen_size()
    win_w, win_h = int(screen_w * 0.95), int(screen_h * 0.95)

    layout = [
        [sg.Column([
            [sg.Text("GRN of gene nodes", font=("Helvetica", 12, "bold"))],
            [sg.Canvas(key="-GRN-CANVAS-", size=(int(win_w*0.8), int(win_h*0.75)))],
            [sg.Text("Clusters of gene nodes", font=("Helvetica", 12, "bold"))],
            [sg.Listbox(values=[], size=(60,10), key="-CLUSTERS-")]
        ], size=(int(win_w*0.8), win_h)),
         sg.VSeperator(),
         sg.Column([
            [sg.Text("Expression-based importance score", font=("Helvetica", 12, "bold"))],
            [sg.Table(values=[],
                      headings=["Gene","SCORE1","Diff","Rob","AUC_absprime"],
                      auto_size_columns=False,
                      col_widths=[15,10,10,10,10],
                      justification="center",
                      alternating_row_color="lightgrey",
                      num_rows=20,
                      key="-TABLE-")],
            [sg.Text("Hot nodes and important interactions", font=("Helvetica", 12, "bold"))],
            [sg.Multiline("", size=(50,7), key="-HOT-", disabled=True)],
            [sg.Text("Trend of gene publications", font=("Helvetica", 12, "bold"))],
            [sg.Canvas(key="-TREND-CANVAS-", size=(int(win_w*0.45), int(win_h*0.35)))]
        ], size=(int(win_w*0.45), win_h))]
    ]

    window = sg.Window("GeNetwork - Page 3", layout, finalize=True, resizable=True, size=(win_w, win_h))

    # -------------------------
    # 1. GRN 네트워크
    # -------------------------
    fig_grn, ax_grn = plt.subplots(figsize=(9,7))
    G = nx.Graph()
    G.add_nodes_from(gene_list)
    for i in range(len(gene_list)-1):  # [Arbitrary code]
        G.add_edge(gene_list[i], gene_list[i+1])
    pos = nx.spring_layout(G)

    nodes = nx.draw_networkx_nodes(G, pos, node_size=600, ax=ax_grn, node_color="skyblue")
    nx.draw_networkx_edges(G, pos, ax=ax_grn)
    nx.draw_networkx_labels(G, pos, font_size=10, ax=ax_grn)
    ax_grn.set_title("Gene Regulatory Network", fontsize=14)
    fig_canvas_agg_grn = draw_figure(window["-GRN-CANVAS-"].TKCanvas, fig_grn)

    selected_gene = None

    def onclick(event):
        nonlocal selected_gene
        if event.inaxes == ax_grn:
            for node, (x,y) in pos.items():
                if (event.xdata - x)**2 + (event.ydata - y)**2 < 0.05:  # threshold
                    selected_gene = node
                    print("Clicked node:", selected_gene)

                    # 강조 표시: 노드 색상 업데이트
                    ax_grn.clear()
                    nx.draw_networkx_nodes(G, pos, node_size=600, ax=ax_grn,
                                           node_color=["red" if n==selected_gene else "skyblue" for n in G.nodes])
                    nx.draw_networkx_edges(G, pos, ax=ax_grn)
                    nx.draw_networkx_labels(G, pos, font_size=10, ax=ax_grn)
                    ax_grn.set_title(f"GRN (Selected: {selected_gene})", fontsize=14)
                    fig_canvas_agg_grn.draw()

                    # TODO: PubMed Trend 분석 호출
                    break

    def onmotion(event):
        """마우스가 노드 위에 올라가면 커서를 손모양으로 변경"""
        if event.inaxes == ax_grn:
            over_node = False
            for _, (x,y) in pos.items():
                if (event.xdata - x)**2 + (event.ydata - y)**2 < 0.05:
                    over_node = True
                    break
            fig_grn.canvas.set_cursor(1 if over_node else 0)  # 1=hand2, 0=arrow

    fig_grn.canvas.mpl_connect("button_press_event", onclick)
    fig_grn.canvas.mpl_connect("motion_notify_event", onmotion)

    # -------------------------
    # 2. Importance Score (변경 없음, 단 표 크기 확대됨)
    # -------------------------
    try:
        data_path = pathlib.Path(r"C:\Users\SAMSUNG\PycharmProjects\GenetWork\data\imaginary_data.csv")
        if data_path.exists():
            df = pd.read_csv(data_path)
            expr = df.drop(columns=["Sample","Labels"]).set_index(df["Sample"])
            labels = df[["Sample","Labels"]].rename(columns={"Sample":"sample","Labels":"label"})

            imp = Importance()
            imp.load_data_from_df(expr, labels)
            imp.compute_expression_scores()
            result = imp.compute_score1()

            valid_genes = list(set(gene_list) & set(expr.columns))
            filtered = result.loc[result.index.intersection(valid_genes)]
            table_data = [[g,
                           round(filtered.loc[g,"SCORE1"],5),
                           round(filtered.loc[g,"Diff"],5),
                           round(filtered.loc[g,"Rob"],5),
                           round(filtered.loc[g,"AUC_absprime"],5)] for g in filtered.index]

            window["-TABLE-"].update(values=table_data)
    except Exception as e:
        print("Importance error:", e)

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break

    window.close()
