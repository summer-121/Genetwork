"""
main.py
GeNetwork - Frontend Launcher
"""

from lib_front_input import run_page1
from lib_front_filter import run_page2
from lib_front_output import run_page3

def main():
    # Page 1 실행
    pdf_list, gene_list = run_page1()
    if not pdf_list and not gene_list:
        print("No input provided. Exiting.")
        return

    # Page 2 실행
    gene_list = run_page2(pdf_list, gene_list)
    if not gene_list:
        print("No gene list after Page 2. Exiting.")
        return

    # Page 3 실행
    run_page3(gene_list)

if __name__ == "__main__":
    main()
