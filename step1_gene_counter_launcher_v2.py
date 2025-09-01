
"""
Step 1 GUI (PyQt5): Gene Counter Launcher — HGNC 정규화 '항상' 적용
-----------------------------------------------------------------
변경점 요약
- ✅ HGNC 정규화를 **무조건 적용**하도록 변경 (선택/체크박스 제거)
- ✅ 기본 경로: ./data/hgnc/hgnc_complete_set.txt (없으면 **자동 다운로드**, 사용자 확인 팝업 없음)
- ✅ 실행 커맨드에 항상 --hgnc <경로> 부여 (다운로드 실패 시 실행 중단)
- ✅ pipeline_config.json에 use_hgnc=True, hgnc_tsv=<경로> 기록
- 기타 기존 기능(드래그앤드롭, 유전자 입력, Step 2 자동 실행, 폴더 열기, 로그 등) 유지
"""

import json
import os
import sys
import time
import warnings
from typing import List, Optional

# --- Windows Qt plugins 경로 보정 (필요 환경에서만 동작) -------------------------------
warnings.filterwarnings("ignore", category=DeprecationWarning)

def _ensure_qt_plugin_path():
    try:
        import PyQt5
        from PyQt5.QtCore import QLibraryInfo
        candidates = []
        try:
            p1 = QLibraryInfo.location(QLibraryInfo.PluginsPath)
            if p1:
                candidates.append(p1)
        except Exception:
            p1 = None
        p2 = os.path.join(os.path.dirname(PyQt5.__file__), "Qt", "plugins")
        candidates.append(p2)
        p3 = os.environ.get("QT_QPA_PLATFORM_PLUGIN_PATH", "")
        if p3:
            candidates.append(p3)
        use = None
        for p in candidates:
            if not p:
                continue
            plat = os.path.join(p, "platforms")
            qwin = os.path.join(plat, "qwindows.dll")
            if os.path.exists(qwin):
                use = p
                break
        if use:
            os.environ["QT_QPA_PLATFORM_PLUGIN_PATH"] = use
            os.environ["QT_PLUGIN_PATH"] = use
            qt_bin = os.path.join(os.path.dirname(use), "bin")
            if os.path.isdir(qt_bin):
                try:
                    if hasattr(os, "add_dll_directory"):
                        os.add_dll_directory(qt_bin)
                    else:
                        os.environ["PATH"] = qt_bin + os.pathsep + os.environ.get("PATH", "")
                except Exception:
                    pass
            print("[Qt] Using PluginsPath:", use)
    except Exception as e:
        print("[Qt] Plugin path setup skipped:", e)

_ensure_qt_plugin_path()

from pathlib import Path
from urllib.request import urlretrieve

from PyQt5.QtCore import Qt, QProcess, QSettings
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QFileDialog,
    QListWidget,
    QToolButton,
    QPushButton,
    QLabel,
    QVBoxLayout,
    QHBoxLayout,
    QTextEdit,
    QMessageBox,
    QCheckBox,
    QFrame,
    QSplitter,
    QStyle,
    QStatusBar,
)

APP_ORG = "Genetwork"
APP_NAME = "GeneCounterStep1"

# HGNC 기본 다운로드 경로/파일명
SCRIPT_DIR = Path(__file__).resolve().parent
HGNC_DEFAULT_URL = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
HGNC_DEFAULT_LOCAL = SCRIPT_DIR / "data" / "hgnc" / "hgnc_complete_set.txt"


class Card(QFrame):
    def __init__(self, title: str = "", parent=None):
        super().__init__(parent)
        self.setObjectName("Card")
        self.setFrameShape(QFrame.StyledPanel)
        self.setFrameShadow(QFrame.Raised)
        self.v = QVBoxLayout(self)
        self.v.setContentsMargins(14, 14, 14, 14)
        self.v.setSpacing(10)
        if title:
            lbl = QLabel(title)
            f = lbl.font(); f.setPointSize(f.pointSize() + 1); f.setBold(True)
            lbl.setFont(f); lbl.setObjectName("CardTitle")
            self.v.addWidget(lbl)


class PdfListWidget(QListWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setSelectionMode(QListWidget.ExtendedSelection)
        self.setAcceptDrops(True)
        self.viewport().setAcceptDrops(True)
        self.setDropIndicatorShown(True)
        self.setDragDropMode(QListWidget.NoDragDrop)
        self.setAlternatingRowColors(True)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            super().dragEnterEvent(event)

    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
        else:
            super().dragMoveEvent(event)

    def dropEvent(self, event):
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                local_path = url.toLocalFile()
                if not local_path:
                    continue
                if os.path.isdir(local_path):
                    self._add_dir_pdfs(local_path)
                elif local_path.lower().endswith(".pdf"):
                    self.add_unique_item(local_path)
            event.acceptProposedAction()
        else:
            super().dropEvent(event)

    def _add_dir_pdfs(self, directory: str):
        for root, _, files in os.walk(directory):
            for name in files:
                if name.lower().endswith(".pdf"):
                    self.add_unique_item(os.path.join(root, name))

    def add_unique_item(self, path: str):
        existing = [self.item(i).text() for i in range(self.count())]
        if path not in existing:
            self.addItem(path)

    def items(self) -> List[str]:
        return [self.item(i).text() for i in range(self.count())]


class Step1Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("STEP 1 — Gene Counter 준비 (HGNC 정규화 항상 적용)")
        self.resize(980, 760)
        QApplication.setStyle("Fusion")
        self.apply_light_styles()

        self.settings = QSettings(APP_ORG, APP_NAME)
        self.last_out_dir = None  # ✅ Step 2 자동 실행 시 사용

        # 프로세스
        self.process = QProcess(self)
        self.process.readyReadStandardOutput.connect(self.on_read_stdout)
        self.process.readyReadStandardError.connect(self.on_read_stderr)
        self.process.finished.connect(self.on_process_finished)

        # === 상단 헤더 ===
        header = Card()
        header_layout = header.v
        title = QLabel("Gene Counter 파이프라인 — Step 1")
        title_f = QFont(); title_f.setPointSize(15); title_f.setBold(True)
        title.setFont(title_f)
        subtitle = QLabel("PDF를 모으고, 유전자 목록(선택)을 입력한 뒤 ‘분석 시작’을 누르세요.\n※ HGNC 정규화는 기본적으로 항상 적용됩니다.")
        subtitle.setObjectName("SubTitle")
        header_layout.addWidget(title)
        header_layout.addWidget(subtitle)

        # === PDF 카드 ===
        pdf_card = Card("PDF 파일")
        pdf_v = pdf_card.v
        pdf_help = QLabel("여러 PDF 또는 폴더를 드래그해서 놓거나, ‘추가’로 선택하세요.")
        self.pdf_list = PdfListWidget()
        tools = QHBoxLayout()
        self.btn_pdf_add = QToolButton(); self.btn_pdf_add.setText("추가")
        self.btn_pdf_add.setIcon(self.style().standardIcon(QStyle.SP_DialogOpenButton))
        self.btn_pdf_add.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.btn_pdf_add_folder = QToolButton(); self.btn_pdf_add_folder.setText("폴더 추가")
        self.btn_pdf_add_folder.setIcon(self.style().standardIcon(QStyle.SP_DirOpenIcon))
        self.btn_pdf_add_folder.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.btn_pdf_remove = QToolButton(); self.btn_pdf_remove.setText("선택 제거")
        self.btn_pdf_remove.setIcon(self.style().standardIcon(QStyle.SP_TrashIcon))
        self.btn_pdf_remove.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.btn_pdf_clear = QToolButton(); self.btn_pdf_clear.setText("전체 지우기")
        self.btn_pdf_clear.setIcon(self.style().standardIcon(QStyle.SP_BrowserStop))
        self.btn_pdf_clear.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        tools.addWidget(self.btn_pdf_add)
        tools.addWidget(self.btn_pdf_add_folder)
        tools.addWidget(self.btn_pdf_remove)
        tools.addStretch(1)
        tools.addWidget(self.btn_pdf_clear)
        pdf_v.addWidget(pdf_help)
        pdf_v.addWidget(self.pdf_list)
        pdf_v.addLayout(tools)

        # === 유전자 직접 입력 카드 ===
        gene_card = Card("유전자 이름 직접 입력 (선택)")
        gene_v = gene_card.v
        gene_help = QLabel("한 줄에 하나씩 입력하세요. 예: TP53\n붙여넣기 후 ‘중복 제거’로 정리할 수 있습니다.")
        gene_row_tools = QHBoxLayout()
        self.btn_paste_clip = QToolButton(); self.btn_paste_clip.setText("클립보드 붙여넣기")
        self.btn_paste_clip.setIcon(self.style().standardIcon(QStyle.SP_DialogOpenButton))
        self.btn_dedup = QToolButton(); self.btn_dedup.setText("중복 제거")
        self.btn_dedup.setIcon(self.style().standardIcon(QStyle.SP_BrowserReload))
        gene_row_tools.addWidget(self.btn_paste_clip)
        gene_row_tools.addWidget(self.btn_dedup)
        gene_row_tools.addStretch(1)
        self.te_genes = QTextEdit(); self.te_genes.setPlaceholderText("TP53\nBRCA1\nEGFR\n…")
        self.te_genes.setAcceptRichText(False)
        gene_v.addWidget(gene_help)
        gene_v.addLayout(gene_row_tools)
        gene_v.addWidget(self.te_genes)

        # === 실행/옵션 행 ===
        actions_row = QHBoxLayout()
        self.chk_autostep2 = QCheckBox("완료 후 자동으로 Step 2 실행"); self.chk_autostep2.setChecked(True)
        self.chk_open_folder = QCheckBox("완료 후 폴더 열기"); self.chk_open_folder.setChecked(True)
        self.btn_run = QPushButton("분석 시작"); self.btn_run.setObjectName("PrimaryButton")
        actions_row.addWidget(self.chk_autostep2)
        actions_row.addWidget(self.chk_open_folder)
        actions_row.addStretch(1)
        actions_row.addWidget(self.btn_run)

        # === 상단 영역 묶기 ===
        upper = QWidget(); upper_v = QVBoxLayout(upper)
        upper_v.setContentsMargins(0, 0, 0, 0); upper_v.setSpacing(10)
        upper_v.addWidget(header)
        upper_v.addWidget(pdf_card)
        upper_v.addWidget(gene_card)
        upper_v.addLayout(actions_row)

        # === 로그 영역 ===
        log_card = Card("실행 로그")
        log_v = log_card.v
        self.log = QTextEdit(); self.log.setReadOnly(True)
        self.log.setPlaceholderText("실행 로그가 여기에 표시됩니다…")
        log_tools = QHBoxLayout()
        self.btn_clear_log = QToolButton(); self.btn_clear_log.setText("로그 지우기")
        self.btn_clear_log.setIcon(self.style().standardIcon(QStyle.SP_BrowserStop))
        self.btn_clear_log.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        log_tools.addStretch(1)
        log_tools.addWidget(self.btn_clear_log)
        log_v.addWidget(self.log)
        log_v.addLayout(log_tools)

        # === 스플리터 ===
        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(upper)
        splitter.addWidget(log_card)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)

        central = QWidget(self)
        central_v = QVBoxLayout(central)
        central_v.setContentsMargins(12, 12, 12, 12); central_v.setSpacing(10)
        central_v.addWidget(splitter)
        self.setCentralWidget(central)

        # 상태바 + 진행바
        self.status = QStatusBar(); self.setStatusBar(self.status)
        from PyQt5.QtWidgets import QProgressBar
        self.progress = QProgressBar(); self.progress.setMaximumHeight(14)
        self.progress.setTextVisible(False); self.progress.hide()
        self.status.addPermanentWidget(self.progress)

        # 시그널 연결
        self.btn_pdf_add.clicked.connect(self.on_add_pdfs)
        self.btn_pdf_add_folder.clicked.connect(self.on_add_pdf_folder)
        self.btn_pdf_remove.clicked.connect(self.on_remove_selected)
        self.btn_pdf_clear.clicked.connect(self.pdf_list.clear)
        self.btn_run.clicked.connect(self.on_run)
        self.btn_clear_log.clicked.connect(self.log.clear)
        self.btn_paste_clip.clicked.connect(self.on_paste_clipboard)
        self.btn_dedup.clicked.connect(self.on_dedup)

        # 세팅 복원
        self.restore_settings()

    # ------------------------ 라이트 스타일 ------------------------
    def apply_light_styles(self):
        self.setStyleSheet(
            """
            QMainWindow { background: #f7f7fb; }
            QLabel#SubTitle { color: #5f6b7a; }
            #Card { background: #ffffff; border: 1px solid #e2e8f0; border-radius: 10px; }
            #CardTitle { color: #0f172a; }
            QTextEdit, QListWidget {
                background: #ffffff; color: #0f172a; border: 1px solid #cbd5e1; border-radius: 8px; padding: 6px;
                selection-background-color: #cfe0ff; }
            QTextEdit:focus, QListWidget:focus { border: 1px solid #3b82f6; }
            QPushButton, QToolButton { background: #ffffff; color: #0f172a; border: 1px solid #cbd5e1; border-radius: 8px; padding: 6px 10px; }
            QPushButton:hover, QToolButton:hover { border-color: #3b82f6; }
            QPushButton#PrimaryButton { background: #2563eb; border-color: #2563eb; color: white; font-weight: 600; padding: 10px 16px; }
            QPushButton#PrimaryButton:hover { background: #3b82f6; border-color: #3b82f6; }
            QCheckBox { color: #0f172a; }
            QStatusBar { background: #ffffff; border-top: 1px solid #e2e8f0; }
            QProgressBar { background: #eef2ff; border: 1px solid #cbd5e1; border-radius: 7px; }
            QProgressBar::chunk { background-color: #3b82f6; border-radius: 7px; }
            """
        )

    # ------------------------ 설정 저장/복원 ------------------------
    def restore_settings(self):
        for p in self.settings.value("pdfs", [], type=list):
            if os.path.exists(p):
                self.pdf_list.add_unique_item(p)
        self.chk_open_folder.setChecked(self.settings.value("open_folder", True, type=bool))
        self.chk_autostep2.setChecked(self.settings.value("auto_step2", True, type=bool))
        self.te_genes.setPlainText(self.settings.value("user_genes_text", ""))

    def save_settings(self):
        self.settings.setValue("pdfs", self.pdf_list.items())
        self.settings.setValue("open_folder", self.chk_open_folder.isChecked())
        self.settings.setValue("auto_step2", self.chk_autostep2.isChecked())
        self.settings.setValue("user_genes_text", self.te_genes.toPlainText())

    # ------------------------ 유틸 ------------------------
    def append_log(self, text: str):
        self.log.append(text.rstrip())
        self.status.showMessage("진행 중…")

    def set_controls_enabled(self, enabled: bool):
        for w in [self.btn_pdf_add, self.btn_pdf_add_folder, self.btn_pdf_remove, self.btn_pdf_clear,
                  self.btn_run, self.pdf_list, self.te_genes,
                  self.btn_paste_clip, self.btn_dedup]:
            w.setEnabled(enabled)

    def parse_user_genes(self) -> List[str]:
        raw = self.te_genes.toPlainText().splitlines()
        seen = set(); genes = []
        for line in raw:
            g = line.strip()
            if not g:
                continue
            if g not in seen:
                seen.add(g); genes.append(g)
        return genes

    def ensure_hgnc_file(self, path: Path) -> Optional[Path]:
        """HGNC TSV가 없으면 자동 다운로드(무조건 수행). 성공 시 경로 반환, 실패 시 None."""
        path = Path(path)
        if path.exists() and path.is_file():
            return path

        # 폴더 생성
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            QMessageBox.critical(self, "오류", f"폴더 생성 실패: {e}")
            return None

        try:
            self.append_log(f"[HGNC] 파일이 없어 자동 다운로드를 시작합니다:\n  URL: {HGNC_DEFAULT_URL}\n  저장: {path}")
            urlretrieve(HGNC_DEFAULT_URL, str(path))
            self.append_log(f"[HGNC] 다운로드 완료: {path}")
            return path
        except Exception as e:
            QMessageBox.critical(self, "오류", f"HGNC TSV 다운로드 실패: {e}")
            return None

    # ------------------------ 입력 핸들러 ------------------------
    def on_add_pdfs(self):
        start_dir = self.guess_dir_from_pdfs()
        files, _ = QFileDialog.getOpenFileNames(self, "PDF 파일 선택", start_dir, "PDF Files (*.pdf)")
        for f in files:
            if f.lower().endswith(".pdf"):
                self.pdf_list.add_unique_item(f)
        self.save_settings()

    def on_add_pdf_folder(self):
        start_dir = self.guess_dir_from_pdfs()
        directory = QFileDialog.getExistingDirectory(self, "PDF 폴더 선택", start_dir)
        if directory:
            self.pdf_list._add_dir_pdfs(directory)
        self.save_settings()

    def on_remove_selected(self):
        for item in self.pdf_list.selectedItems():
            row = self.pdf_list.row(item)
            self.pdf_list.takeItem(row)
        self.save_settings()

    def on_paste_clipboard(self):
        cb = QApplication.clipboard().text()
        if cb:
            current = self.te_genes.toPlainText()
            if current and not current.endswith("\n"):
                current += "\n"
            self.te_genes.setPlainText(current + cb)
            self.save_settings()

    def on_dedup(self):
        genes = self.parse_user_genes()
        self.te_genes.setPlainText("\n".join(genes))
        self.save_settings()

    # ------------------------ 실행 플로우 ------------------------
    def validate_inputs(self) -> bool:
        if not self.pdf_list.items():
            QMessageBox.warning(self, "확인", "PDF 파일을 하나 이상 추가하세요.")
            return False
        if not self.parse_user_genes():
            ret = QMessageBox.question(
                self, "확인", "유전자 직접 입력이 비어 있습니다. 이대로 진행할까요?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes)
            if ret != QMessageBox.Yes:
                return False
        script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gene_counter.py")
        if not os.path.isfile(script_path):
            QMessageBox.critical(self, "오류", "현재 폴더에 gene_counter.py가 없습니다. 동일 폴더에 두고 실행하세요.")
            return False
        return True

    def build_command(self, out_dir: str) -> List[str]:
        script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gene_counter.py")
        pdfs = self.pdf_list.items()
        global_csv = os.path.join(out_dir, "global_counts.csv")
        per_file_csv = os.path.join(out_dir, "per_file_counts.csv")
        cmd = [sys.executable, script_path]
        cmd.extend(pdfs)
        cmd.extend(["--global-out", global_csv, "--per-file-out", per_file_csv])

        # ✅ HGNC 정규화 '항상' 적용
        ok_path = self.ensure_hgnc_file(HGNC_DEFAULT_LOCAL)
        if ok_path is None:
            raise RuntimeError("HGNC TSV 준비 실패로 실행을 중단합니다.")
        cmd.extend(["--hgnc", str(ok_path)])

        return cmd

    def compute_out_dir(self) -> str:
        first_pdf = self.pdf_list.items()[0]
        base = os.path.dirname(first_pdf)
        stamp = time.strftime("%Y%m%d-%H%M%S")
        out_dir = os.path.join(base, f"gene_counter_results-{stamp}")
        os.makedirs(out_dir, exist_ok=True)
        return out_dir

    def write_pipeline_config(self, out_dir: str):
        cfg_path = os.path.join(out_dir, "pipeline_config.json")
        data = {
            "global_csv": os.path.join(out_dir, "global_counts.csv"),
            "per_file_csv": os.path.join(out_dir, "per_file_counts.csv"),
            "user_genes": self.parse_user_genes(),
            # 기록: 항상 정규화 적용
            "use_hgnc": True,
            "hgnc_tsv": str(HGNC_DEFAULT_LOCAL),
        }
        try:
            with open(cfg_path, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
            self.append_log(f"pipeline_config.json 저장: {cfg_path}")
        except Exception as e:
            QMessageBox.warning(self, "경고", f"pipeline_config.json 저장 중 오류: {e}")

    def on_run(self):
        if not self.validate_inputs():
            return
        self.log.clear(); self.append_log("분석을 시작합니다…")
        self.set_controls_enabled(False)
        self.progress.show(); self.progress.setRange(0, 0)
        self.status.showMessage("실행 중…")

        out_dir = self.compute_out_dir()
        self.last_out_dir = out_dir  # ✅ 기억
        self.append_log(f"출력 폴더: {out_dir}")

        try:
            cmd = self.build_command(out_dir)
        except Exception as e:
            # HGNC 준비 실패 등
            self.progress.hide(); self.set_controls_enabled(True)
            QMessageBox.critical(self, "실행 중단", str(e))
            return

        pretty = " ".join([f'"{c}"' if " " in c else c for c in cmd])
        self.append_log(f"실행 커맨드:\n  {pretty}")

        self.process.start(cmd[0], cmd[1:])
        if not self.process.waitForStarted(3000):
            self.progress.hide(); self.set_controls_enabled(True)
            QMessageBox.critical(self, "오류", "프로세스 시작에 실패했습니다. Python 또는 스크립트 경로를 확인하세요.")
            return

    # ------------------------ QProcess 콜백 ------------------------
    def on_read_stdout(self):
        data = bytes(self.process.readAllStandardOutput()).decode(errors="ignore")
        if data:
            self.append_log(data)

    def on_read_stderr(self):
        data = bytes(self.process.readAllStandardError()).decode(errors="ignore")
        if data:
            self.append_log(data)

    def on_process_finished(self, exit_code: int, exit_status):
        self.append_log(f"\n프로세스 종료 코드: {exit_code}")
        # 로그에서 출력 폴더를 재획득
        out_dir = None
        for line in self.log.toPlainText().splitlines():
            if line.startswith("출력 폴더: "):
                out_dir = line.split(": ", 1)[1].strip()
        if not out_dir and self.last_out_dir:
            out_dir = self.last_out_dir
        if out_dir:
            self.write_pipeline_config(out_dir)

        self.set_controls_enabled(True)
        self.progress.hide()

        if exit_code == 0:
            # ✅ Step 2 자동 실행
            launched = False
            if out_dir and self.chk_autostep2.isChecked():
                step2_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "step2_grn_selector.py")
                cfg_path = os.path.join(out_dir, "pipeline_config.json")
                if os.path.isfile(step2_path) and os.path.isfile(cfg_path):
                    QProcess.startDetached(sys.executable, [step2_path, cfg_path])
                    self.append_log(f'Step 2 실행: "{step2_path}" "{cfg_path}"')
                    launched = True

            self.status.showMessage("완료", 5000)
            if out_dir:
                msg = "분석이 완료되었습니다.\n" + f"출력 폴더: {out_dir}\n\n"
                msg += ("Step 2가 자동으로 실행되었습니다." if launched else "Step 2에서 pipeline_config.json을 읽어 다음 단계로 진행하세요.")
                QMessageBox.information(self, "완료", msg)
                if self.chk_open_folder.isChecked():
                    self.open_folder(out_dir)
        else:
            self.status.showMessage("실행 오류", 5000)
            QMessageBox.warning(self, "실행 오류", "실행 중 오류가 발생했습니다. 로그를 확인하세요.")

    # ------------------------ 기타 ------------------------
    def open_folder(self, path: str):
        try:
            if sys.platform.startswith("win"):
                os.startfile(path)
            elif sys.platform == "darwin":
                os.system(f'open "{path}"')
            else:
                os.system(f'xdg-open "{path}"')
        except Exception:
            pass

    def guess_dir_from_pdfs(self) -> str:
        items = self.pdf_list.items()
        if items:
            base = os.path.dirname(items[0])
            if os.path.isdir(base):
                return base
        return os.path.dirname(os.path.abspath(__file__))

    def closeEvent(self, event):
        self.save_settings(); super().closeEvent(event)


def main():
    app = QApplication(sys.argv)
    win = Step1Window(); win.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
