"""
main.py
GeNetwork - Frontend Launcher
"""

from lib_front_input import run_page1
from lib_front_filter import run_page2

def main():
    # Page 1 실행
    pdf_list, gene_list = run_page1()
    if not pdf_list and not gene_list:
        print("No input provided. Exiting.")
        return

    # Page 2 실행
    gene_list = run_page2(pdf_list, gene_list)

    print("최종 gene list:", gene_list)

if __name__ == "__main__":
    main()
