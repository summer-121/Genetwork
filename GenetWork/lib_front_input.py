"""
lib_front_input.py
GeNetwork - Frontend Page 1 (PySimpleGUI + backend Data.search_gene)

기능:
 1. 논문 업로드 및 리스트 표시
 2. NCBI Gene 검색 (Data.search_gene 이용)
 3. Next 버튼 → Page 2 이동 (현재 미구현 안내)
"""

import PySimpleGUI as sg     # 최신 버전으로 해야 함. 그냥 pip install PySimpleGUI로 하면 오류 발생.
from lib_back import Data


# ---------- Data API ----------
data_api = Data()

def run_page1():
    sg.theme('LightBlue2')

    layout = [
        [sg.Text("GeNetwork 1.0.0 - Page 1", font=("Helvetica", 14, "bold"))],

        [sg.Text("Upload your article PDF file HERE")],
        [sg.Input(key='-PDF-PATH-', enable_events=True, size=(40,1)),
         sg.FilesBrowse("Upload PDF(s)", file_types=(("PDF Files","*.pdf"),), target='-PDF-PATH-')],
        [sg.Listbox(values=[], size=(50,6), key='-PDF-LIST-', enable_events=True)],
        [sg.Button("Remove Selected PDF", key='-PDF-REMOVE-'),
         sg.Button("Clear PDFs", key='-PDF-CLEAR-')],

        [sg.HorizontalSeparator()],

        [sg.Text("Search your gene name (via NCBI)")],
        [sg.Input(key='-NCBI-QUERY-', size=(30,1)), sg.Button("NCBI Search", key='-NCBI-SEARCH-')],
        [sg.Listbox(values=[], size=(50,6), key='-NCBI-RESULTS-', enable_events=True)],
        [sg.Button("Add Selected Gene", key='-NCBI-ADD-')],
        [sg.Text("Your gene list:")],
        [sg.Listbox(values=[], size=(50,6), key='-GENE-LIST-', enable_events=True)],
        [sg.Button("Remove Selected Gene", key='-GENE-REMOVE-')],

        [sg.HorizontalSeparator()],
        [sg.Button("Next →", key='-NEXT-', size=(10,1))],
        [sg.Text("", key='-STATUS-', size=(60,1), text_color='green')],
    ]

    window = sg.Window("GeNetwork - Page 1", layout, finalize=True, resizable=True)

    pdf_list = []
    gene_list = []

    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED:
            break

        # PDF 업로드
        if event == '-PDF-PATH-':
            paths = values['-PDF-PATH-']
            if paths:
                for p in paths.split(';'):
                    p = p.strip()
                    if p and p not in pdf_list:
                        pdf_list.append(p)
                window['-PDF-LIST-'].update(pdf_list)
                window['-STATUS-'].update(f"{len(pdf_list)} PDF(s) uploaded")

        if event == '-PDF-REMOVE-':
            sel = values['-PDF-LIST-']
            if sel:
                for s in sel:
                    if s in pdf_list:
                        pdf_list.remove(s)
                window['-PDF-LIST-'].update(pdf_list)
                window['-STATUS-'].update("Removed selected PDF(s)")

        if event == '-PDF-CLEAR-':
            pdf_list.clear()
            window['-PDF-LIST-'].update(pdf_list)
            window['-STATUS-'].update("Cleared PDFs")

        # Gene 검색 (백엔드 호출)
        if event == '-NCBI-SEARCH-':
            q = values['-NCBI-QUERY-'].strip()
            if not q:
                window['-STATUS-'].update("검색어를 입력하세요.")
            else:
                try:
                    result = data_api.search_gene(q)  # [(gene_id, gene_name), ...]
                    display_results = [f"{gid}: {gname}" for gid, gname in result]
                    window['-NCBI-RESULTS-'].update(display_results)
                    window['-STATUS-'].update(f"검색 결과 {len(result)}개 반환")
                except Exception as e:
                    window['-STATUS-'].update(f"검색 실패: {e}")

        # Gene 리스트에 추가
        if event == '-NCBI-ADD-':
            sel = values['-NCBI-RESULTS-']
            if sel:
                for item in sel:
                    gname = item.split(":",1)[1].strip()
                    if gname not in gene_list:
                        gene_list.append(gname)
                window['-GENE-LIST-'].update(gene_list)
                window['-STATUS-'].update(f"Added {len(sel)} gene(s)")

        # Gene 제거
        if event == '-GENE-REMOVE-':
            sel = values['-GENE-LIST-']
            if sel:
                for g in sel:
                    if g in gene_list:
                        gene_list.remove(g)
                window['-GENE-LIST-'].update(gene_list)
                window['-STATUS-'].update("Removed selected gene(s)")

        # Next 버튼
        if event == '-NEXT-':
            window.close()
            return pdf_list, gene_list

    window.close()
    return pdf_list, gene_list
