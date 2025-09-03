import sys
from pathlib import Path
from PyQt5.QtWidgets import QApplication
from library2 import Step1Window, Step2Window

# Step1에 필요한 상수
APP_ORG = "Genetwork"
APP_NAME = "GeneCounterStep1"
SCRIPT_DIR = Path(__file__).resolve().parent
HGNC_DEFAULT_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_DEFAULT_LOCAL = SCRIPT_DIR / "data" / "hgnc" / "hgnc_complete_set.txt"


def run_step1():
    app = QApplication(sys.argv)
    win = Step1Window(APP_ORG, APP_NAME, HGNC_DEFAULT_URL, HGNC_DEFAULT_LOCAL)
    win.show()
    sys.exit(app.exec_())


def run_step2():
    app = QApplication(sys.argv)
    win = Step2Window()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "step2":
        run_step2()
    else:
        run_step1()