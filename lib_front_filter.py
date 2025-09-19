"""
lib_front_page2.py
GeNetwork - Frontend Page 2 (PySimpleGUI)

기능:
 1. Page 1에서 전달된 PDF 파일 리스트 → Text 클래스 호출 → gene list 생성
 2. Page 1 gene list와 통합 (중복 제거)
 3. 사용자 편집 (추가/삭제/저장)
"""
# 터미널 환경: pip install https://s3-us-west-2.amazonaws.com/ai2-s2-scispacy/releases/v0.5.1/en_ner_jnlpba_md-0.5.1.tar.gz 필수적으로 해야 함
import PySimpleGUI as sg
import pathlib
from lib_back import Text

# ---------- Text API 초기화 ----------
text_api = Text()

# ---------- Page 2 실행 ----------
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
        [sg.Button("Next →", key="-NEXT-", size=(10,1))],
        [sg.Text("", key="-STATUS-", size=(60,1), text_color="green")],
    ]

    window = sg.Window("GeNetwork - Page 2", layout, finalize=True, resizable=True)

    mining_genes = []
    gene_list = prev_gene_list.copy()

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break

        # Step 1. Run Text Mining
        if event == "-RUN-MINING-":
            try:
                mining_genes = []
                for pdf in pdf_list:
                    pdf_path = pathlib.Path(pdf)
                    if not pdf_path.exists():
                        continue
                    # 실제 백엔드 호출
                    genes = text_api.extract_genes_from_pdf(pdf_path)
                    mining_genes.extend(genes)
                mining_genes = list(set(mining_genes))  # 중복 제거
                window["-MINING-RESULTS-"].update(mining_genes)
                window["-STATUS-"].update(f"Text mining completed, {len(mining_genes)} unique genes found.")
            except Exception as e:
                window["-STATUS-"].update(f"Text mining failed: {e}")

        # Step 2. Merge with previous list
        if event == "-MERGE-":
            merged = set(gene_list) | set(mining_genes)
            gene_list = sorted(list(merged))  # [Arbitrary code] 보기 좋게 정렬
            window["-GENE-LIST-"].update(gene_list)
            window["-STATUS-"].update(f"Merged list contains {len(gene_list)} unique genes.")

        # Gene 삭제
        if event == "-GENE-REMOVE-":
            sel = values["-GENE-LIST-"]
            if sel:
                for g in sel:
                    if g in gene_list:
                        gene_list.remove(g)
                window["-GENE-LIST-"].update(gene_list)
                window["-STATUS-"].update("Removed selected gene(s)")

        # Gene 수동 추가
        if event == "-GENE-ADD-":
            newg = sg.popup_get_text("Enter gene name to add:", "Add Gene")
            if newg:
                newg = newg.strip().upper()
                if newg not in gene_list:
                    gene_list.append(newg)
                    window["-GENE-LIST-"].update(gene_list)
                    window["-STATUS-"].update(f"Added gene: {newg}")
                else:
                    sg.popup("That gene already exists in the list.")

        # Gene list 내보내기
        if event == "-GENE-EXPORT-":
            if not gene_list:
                sg.popup("Gene list is empty")
            else:
                save_path = sg.popup_get_file("Save gene list", save_as=True, no_window=True,
                                              default_extension=".txt", file_types=(("Text Files","*.txt"),))
                if save_path:
                    try:
                        with open(save_path, "w", encoding="utf-8") as f:
                            for g in gene_list:
                                f.write(g + "\n")
                        sg.popup("Saved gene list to:", save_path)
                    except Exception as e:
                        sg.popup_error(f"Failed to save: {e}")

        # Next 버튼
        if event == "-NEXT-":
            sg.popup("Next page (Page 3) 기능은 아직 구현되지 않았습니다.", title="Next Page")
            window["-STATUS-"].update("Moving to Page 3... (not yet implemented)")

    window.close()
    return gene_list