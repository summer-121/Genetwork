"""
lib_front_page2.py
GeNetwork - Frontend Page 2 (PySimpleGUI)
"""

import PySimpleGUI as sg
import pathlib
from lib_back import Text
from lib_front_output import run_page3   # Page 3 연동

text_api = Text()

def run_page2(pdf_list, prev_gene_list):
    sg.theme("LightBlue2")

    layout = [
        [sg.Text("GeNetwork - Page 2", font=("Helvetica", 14, "bold"))],

        [sg.Text("Step 1. Run text mining on uploaded PDFs")],
        [sg.Button("Run Text Mining", key="-RUN-MINING-")],
        [sg.Text("Text mining results:")],
        [sg.Listbox(values=[], size=(50,6), key="-MINING-RESULTS-")],

        [sg.HorizontalSeparator()],

        [sg.Text("Step 2. Merge with previous gene list (remove duplicates)")],
        [sg.Button("Merge Gene Lists", key="-MERGE-")],
        [sg.Text("Merged Gene List (editable):")],
        [sg.Listbox(values=prev_gene_list.copy(), size=(50,10), key="-GENE-LIST-", enable_events=True)],

        [sg.Button("Remove Selected Gene", key="-GENE-REMOVE-"),
         sg.Button("Add Gene Manually", key="-GENE-ADD-"),
         sg.Button("Export Gene List", key="-GENE-EXPORT-")],

        [sg.HorizontalSeparator()],
        [sg.Button("Analyze and Create New GRN!", key="-NEXT-", size=(50,1))],
        [sg.Text("", key="-STATUS-", size=(60,1), text_color="green")],
    ]

    window = sg.Window("GeNetwork - Page 2", layout, finalize=True, resizable=True)

    mining_genes = []
    gene_list = prev_gene_list.copy()

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break

        if event == "-RUN-MINING-":
            try:
                mining_genes = []
                for pdf in pdf_list:
                    pdf_path = pathlib.Path(pdf)
                    if pdf_path.exists():
                        genes = text_api.extract_genes_from_pdf(pdf_path)
                        mining_genes.extend(genes)
                mining_genes = list(set(mining_genes))
                window["-MINING-RESULTS-"].update(mining_genes)
                window["-STATUS-"].update(f"Text mining completed, {len(mining_genes)} unique genes found.")
            except Exception as e:
                window["-STATUS-"].update(f"Text mining failed: {e}")

        if event == "-MERGE-":
            merged = set(gene_list) | set(mining_genes)
            gene_list = sorted(list(merged))  # [Arbitrary code] 보기 좋게 정렬
            window["-GENE-LIST-"].update(gene_list)
            window["-STATUS-"].update(f"Merged list contains {len(gene_list)} unique genes.")

        if event == "-GENE-REMOVE-":
            sel = values["-GENE-LIST-"]
            if sel:
                for g in sel:
                    if g in gene_list:
                        gene_list.remove(g)
                window["-GENE-LIST-"].update(gene_list)

        if event == "-GENE-ADD-":
            newg = sg.popup_get_text("Enter gene name to add:", "Add Gene")
            if newg:
                newg = newg.strip().upper()
                if newg not in gene_list:
                    gene_list.append(newg)
                    window["-GENE-LIST-"].update(gene_list)

        if event == "-GENE-EXPORT-":
            if gene_list:
                save_path = sg.popup_get_file("Save gene list", save_as=True, no_window=True,
                                              default_extension=".txt", file_types=(("Text Files","*.txt"),))
                if save_path:
                    with open(save_path, "w", encoding="utf-8") as f:
                        for g in gene_list:
                            f.write(g + "\n")

        if event == "-NEXT-":
            window.close()  # Page 2 닫기
            run_page3(gene_list)  # Page 3 실행
            return gene_list

    window.close()
    return gene_list
