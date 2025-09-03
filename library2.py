import os
import sys
import time
import json
import csv
from pathlib import Path
from urllib.request import urlretrieve
from typing import List, Dict, Tuple, Optional

from PyQt5.QtCore import Qt, QProcess, QSettings, QEvent
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QApplication, QFileDialog, QListWidget, QListWidgetItem,
    QToolButton, QPushButton, QLabel, QVBoxLayout, QHBoxLayout, QTextEdit, QLineEdit,
    QMessageBox, QCheckBox, QFrame, QSplitter, QStyle, QStatusBar, QProgressBar, QSpinBox
)

# -------- 공통 UI --------
class Card(QFrame):
    def __init__(self, title: str = "", parent=None):
        super().__init__(parent)
        self.setObjectName("Card")
        self.setFrameShape(QFrame.StyledPanel)
        self.setFrameShadow(QFrame.Raised)
        self.v = QVBoxLayout(self)
        self.v.setContentsMargins(12, 12, 12, 12)
        self.v.setSpacing(8)
        if title:
            lbl = QLabel(title)
            f = lbl.font(); f.setBold(True)
            lbl.setFont(f); lbl.setObjectName("CardTitle")
            self.v.addWidget(lbl)

class PdfListWidget(QListWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setSelectionMode(QListWidget.ExtendedSelection)
        self.setAcceptDrops(True); self.viewport().setAcceptDrops(True)
        self.setDropIndicatorShown(True)
        self.setDragDropMode(QListWidget.NoDragDrop)
        self.setAlternatingRowColors(True)

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls(): event.acceptProposedAction()
        else: super().dragEnterEvent(event)

    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls(): event.acceptProposedAction()
        else: super().dragMoveEvent(event)

    def dropEvent(self, event):
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                p = url.toLocalFile()
                if not p: continue
                if os.path.isdir(p): self._add_dir_pdfs(p)
                elif p.lower().endswith(".pdf"): self.add_unique_item(p)
            event.acceptProposedAction()
        else: super().dropEvent(event)

    def _add_dir_pdfs(self, directory: str):
        for root, _, files in os.walk(directory):
            for name in files:
                if name.lower().endswith(".pdf"):
                    self.add_unique_item(os.path.join(root, name))  # :contentReference[oaicite:9]{index=9}

    def add_unique_item(self, path: str):
        if path not in [self.item(i).text() for i in range(self.count())]:
            self.addItem(path)  # :contentReference[oaicite:10]{index=10}

    def items(self) -> List[str]:
        return [self.item(i).text() for i in range(self.count())]  # :contentReference[oaicite:11]{index=11}

class CheckList(QListWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAlternatingRowColors(True)
        self._checked: Dict[str, bool] = {}
        self._items: List[Tuple[str, str]] = []
        self._filter = ""

    def clear_all(self):
        self._checked.clear(); self._items.clear(); super().clear()

    def set_data(self, pairs: List[Tuple[str, str]]):
        self._items = pairs
        for g, _ in pairs:
            self._checked.setdefault(g, False)
        self._rebuild()  # :contentReference[oaicite:12]{index=12}

    def _rebuild(self):
        super().clear()
        f = self._filter.lower().strip()
        for gene, text in self._items:
            if f and (f not in gene.lower()) and (f not in text.lower()):
                continue
            it = QListWidgetItem(text)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if self._checked.get(gene, False) else Qt.Unchecked)
            it.setData(Qt.UserRole, gene)
            self.addItem(it)  # :contentReference[oaicite:13]{index=13}

    def set_filter(self, text: str):
        self._filter = text; self._rebuild()

    def select_all(self, checked: bool):
        for i in range(self.count()):
            it = self.item(i)
            it.setCheckState(Qt.Checked if checked else Qt.Unchecked)
            self._checked[it.data(Qt.UserRole)] = checked  # :contentReference[oaicite:14]{index=14}

    def set_checked_genes(self, genes: List[str], checked: bool = True):
        for g in genes: self._checked[g] = checked
        self._rebuild()

    def checked_genes(self) -> List[str]:
        return [g for g, v in self._checked.items() if v]

    def all_genes(self) -> List[str]:
        return [g for g, _ in self._items]

    def on_item_changed(self, item: QListWidgetItem):
        self._checked[item.data(Qt.UserRole)] = (item.checkState() == Qt.Checked)  # :contentReference[oaicite:15]{index=15}

# -------- Step1 --------
class Step1Window(QMainWindow):
    def __init__(self, app_org: str, app_name: str, hgnc_url: str, hgnc_local: Path):
        super().__init__()
        self.app_org = app_org
        self.app_name = app_name
        self.hgnc_url = hgnc_url
        self.hgnc_local = Path(hgnc_local)

        self.setWindowTitle("STEP 1 — Gene Counter 준비 (HGNC 정규화 항상 적용)")
        self.resize(980, 760)
        QApplication.setStyle("Fusion")

        self.settings = QSettings(self.app_org, self.app_name)
        self.last_out_dir: Optional[str] = None

        self.process = QProcess(self)
        self.process.readyReadStandardOutput.connect(self.on_read_stdout)
        self.process.readyReadStandardError.connect(self.on_read_stderr)
        self.process.finished.connect(self.on_process_finished)

        self._build_ui()
        self.apply_light_styles()
        self.restore_settings()  # :contentReference[oaicite:16]{index=16}

    def _build_ui(self):
        header = Card()
        title = QLabel("Gene Counter 파이프라인 — Step 1")
        f = QFont(); f.setPointSize(15); f.setBold(True); title.setFont(f)
        subtitle = QLabel("PDF를 모으고, 유전자 목록(선택)을 입력한 뒤 ‘분석 시작’을 누르세요.\n※ HGNC 정규화는 항상 적용됩니다.")
        subtitle.setObjectName("SubTitle")
        header.v.addWidget(title); header.v.addWidget(subtitle)

        pdf_card = Card("PDF 파일")
        self.pdf_list = PdfListWidget()
        tools = QHBoxLayout()
        self.btn_pdf_add = QToolButton(); self.btn_pdf_add.setText("추가")
        self.btn_pdf_add.setIcon(self.style().standardIcon(QStyle.SP_DialogOpenButton))
        self.btn_pdf_add_folder = QToolButton(); self.btn_pdf_add_folder.setText("폴더 추가")
        self.btn_pdf_add_folder.setIcon(self.style().standardIcon(QStyle.SP_DirOpenIcon))
        self.btn_pdf_remove = QToolButton(); self.btn_pdf_remove.setText("선택 제거")
        self.btn_pdf_remove.setIcon(self.style().standardIcon(QStyle.SP_TrashIcon))
        self.btn_pdf_clear = QToolButton(); self.btn_pdf_clear.setText("전체 지우기")
        self.btn_pdf_clear.setIcon(self.style().standardIcon(QStyle.SP_TrashIcon))
        tools.addWidget(self.btn_pdf_add); tools.addWidget(self.btn_pdf_add_folder)
        tools.addWidget(self.btn_pdf_remove); tools.addStretch(1); tools.addWidget(self.btn_pdf_clear)
        pdf_card.v.addWidget(self.pdf_list); pdf_card.v.addLayout(tools)

        gene_card = Card("유전자 입력 (선택)")
        self.btn_paste_clip = QToolButton(); self.btn_paste_clip.setText("클립보드 붙여넣기")
        self.btn_dedup = QToolButton(); self.btn_dedup.setText("중복 제거")
        row2 = QHBoxLayout(); row2.addWidget(self.btn_paste_clip); row2.addWidget(self.btn_dedup); row2.addStretch(1)
        self.te_genes = QTextEdit()
        gene_card.v.addLayout(row2); gene_card.v.addWidget(self.te_genes)

        actions = QHBoxLayout()
        self.chk_autostep2 = QCheckBox("완료 후 자동으로 Step 2 실행"); self.chk_autostep2.setChecked(True)
        self.chk_open_folder = QCheckBox("완료 후 폴더 열기"); self.chk_open_folder.setChecked(True)
        self.btn_run = QPushButton("분석 시작"); self.btn_run.setObjectName("PrimaryButton")
        actions.addWidget(self.chk_autostep2); actions.addWidget(self.chk_open_folder)
        actions.addStretch(1); actions.addWidget(self.btn_run)  # :contentReference[oaicite:17]{index=17}

        log_card = Card("실행 로그")
        self.log = QTextEdit(); self.log.setReadOnly(True)
        log_tools = QHBoxLayout(); self.btn_clear_log = QToolButton(); self.btn_clear_log.setText("로그 지우기")
        log_tools.addStretch(1); log_tools.addWidget(self.btn_clear_log)
        log_card.v.addWidget(self.log); log_card.v.addLayout(log_tools)

        splitter = QSplitter(Qt.Vertical)
        upper = QWidget(); uv = QVBoxLayout(upper)
        uv.setContentsMargins(0,0,0,0); uv.setSpacing(10)
        uv.addWidget(header); uv.addWidget(pdf_card); uv.addWidget(gene_card); uv.addLayout(actions)
        splitter.addWidget(upper); splitter.addWidget(log_card)
        splitter.setStretchFactor(0, 3); splitter.setStretchFactor(1, 2)

        central = QWidget(self); cv = QVBoxLayout(central)
        cv.setContentsMargins(12,12,12,12); cv.setSpacing(10); cv.addWidget(splitter)
        self.setCentralWidget(central)

        self.status = QStatusBar(); self.setStatusBar(self.status)
        self.progress = QProgressBar(); self.progress.setMaximumHeight(14)
        self.progress.setTextVisible(False); self.progress.hide()
        self.status.addPermanentWidget(self.progress)  # :contentReference[oaicite:18]{index=18}

        self.btn_pdf_add.clicked.connect(self.on_add_pdfs)
        self.btn_pdf_add_folder.clicked.connect(self.on_add_pdf_folder)
        self.btn_pdf_remove.clicked.connect(self.on_remove_selected)
        self.btn_pdf_clear.clicked.connect(self.pdf_list.clear)
        self.btn_run.clicked.connect(self.on_run)
        self.btn_clear_log.clicked.connect(self.log.clear)
        self.btn_paste_clip.clicked.connect(self.on_paste_clipboard)
        self.btn_dedup.clicked.connect(self.on_dedup)

    # --- 스타일 ---
    def apply_light_styles(self):
        self.setStyleSheet("""
            QMainWindow { background: #f7f7fb; }
            QLabel#SubTitle { color: #5f6b7a; }
            #Card { background: #ffffff; border: 1px solid #e2e8f0; border-radius: 10px; }
            #CardTitle { color: #0f172a; }
            QTextEdit, QListWidget, QLineEdit { background: #ffffff; }
        """)  # :contentReference[oaicite:19]{index=19}

    # --- 세팅/유틸 ---
    def restore_settings(self):
        for p in self.settings.value("pdfs", [], type=list):
            if os.path.exists(p): self.pdf_list.add_unique_item(p)
        self.chk_open_folder.setChecked(self.settings.value("open_folder", True, type=bool))
        self.chk_autostep2.setChecked(self.settings.value("auto_step2", True, type=bool))
        self.te_genes.setPlainText(self.settings.value("user_genes_text", ""))  # :contentReference[oaicite:20]{index=20}

    def save_settings(self):
        self.settings.setValue("pdfs", self.pdf_list.items())
        self.settings.setValue("open_folder", self.chk_open_folder.isChecked())
        self.settings.setValue("auto_step2", self.chk_autostep2.isChecked())
        self.settings.setValue("user_genes_text", self.te_genes.toPlainText())  # :contentReference[oaicite:21]{index=21}

    def append_log(self, text: str):
        self.log.append(text.rstrip()); self.status.showMessage("진행 중…")  # :contentReference[oaicite:22]{index=22}

    def set_controls_enabled(self, enabled: bool):
        for w in [self.btn_pdf_add, self.btn_pdf_add_folder, self.btn_pdf_remove, self.btn_pdf_clear,
                  self.btn_run, self.pdf_list, self.te_genes, self.btn_paste_clip, self.btn_dedup]:
            w.setEnabled(enabled)  # :contentReference[oaicite:23]{index=23}

    def parse_user_genes(self) -> List[str]:
        seen=set(); out=[]
        for line in self.te_genes.toPlainText().splitlines():
            g=line.strip()
            if g and g not in seen:
                seen.add(g); out.append(g)
        return out  # :contentReference[oaicite:24]{index=24}

    def ensure_hgnc_file(self, path: Path) -> Optional[Path]:
        path = Path(path)
        if path.exists() and path.is_file(): return path
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
        except Exception as e:
            QMessageBox.critical(self, "오류", f"폴더 생성 실패: {e}"); return None
        try:
            self.append_log(f"[HGNC] 다운로드: {self.hgnc_url}\n  -> {path}")
            urlretrieve(self.hgnc_url, str(path))
            self.append_log(f"[HGNC] 완료: {path}")
            return path
        except Exception as e:
            QMessageBox.critical(self, "오류", f"HGNC TSV 다운로드 실패: {e}")
            return None  # :contentReference[oaicite:25]{index=25}

    def guess_dir_from_pdfs(self) -> str:
        return os.path.dirname(self.pdf_list.items()[0]) if self.pdf_list.count() else os.path.expanduser("~")

    # --- 입력 핸들러 ---
    def on_add_pdfs(self):
        start = self.guess_dir_from_pdfs()
        files, _ = QFileDialog.getOpenFileNames(self, "PDF 파일 선택", start, "PDF Files (*.pdf)")
        for f in files:
            if f.lower().endswith(".pdf"): self.pdf_list.add_unique_item(f)
        self.save_settings()  # :contentReference[oaicite:26]{index=26}

    def on_add_pdf_folder(self):
        start = self.guess_dir_from_pdfs()
        d = QFileDialog.getExistingDirectory(self, "PDF 폴더 선택", start)
        if d: self.pdf_list._add_dir_pdfs(d)
        self.save_settings()  # :contentReference[oaicite:27]{index=27}

    def on_remove_selected(self):
        for it in self.pdf_list.selectedItems():
            self.pdf_list.takeItem(self.pdf_list.row(it))
        self.save_settings()  # :contentReference[oaicite:28]{index=28}

    def on_paste_clipboard(self):
        cb = QApplication.clipboard().text()
        if cb:
            cur = self.te_genes.toPlainText()
            if cur and not cur.endswith("\n"): cur += "\n"
            self.te_genes.setPlainText(cur + cb)
            self.save_settings()  # :contentReference[oaicite:29]{index=29}

    def on_dedup(self):
        self.te_genes.setPlainText("\n".join(self.parse_user_genes()))
        self.save_settings()  # :contentReference[oaicite:30]{index=30}

    # --- 실행 플로우 ---
    def validate_inputs(self) -> bool:
        if not self.pdf_list.items():
            QMessageBox.warning(self, "확인", "PDF 파일을 하나 이상 추가하세요."); return False
        script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gene_counter.py")
        if not os.path.isfile(script_path):
            QMessageBox.critical(self, "오류", "현재 폴더에 gene_counter.py가 없습니다."); return False
        return True  # :contentReference[oaicite:31]{index=31}

    def build_command(self, out_dir: str) -> List[str]:
        script_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "gene_counter.py")
        pdfs = self.pdf_list.items()
        global_csv = os.path.join(out_dir, "global_counts.csv")
        per_file_csv = os.path.join(out_dir, "per_file_counts.csv")
        cmd = [sys.executable, script_path, *pdfs, "--global-out", global_csv, "--per-file-out", per_file_csv]
        ok_path = self.ensure_hgnc_file(self.hgnc_local)
        if ok_path is None: raise RuntimeError("HGNC TSV 준비 실패")
        cmd.extend(["--hgnc", str(ok_path)])
        return cmd  # :contentReference[oaicite:32]{index=32}

    def compute_out_dir(self) -> str:
        base = os.path.dirname(self.pdf_list.items()[0])
        stamp = time.strftime("%Y%m%d-%H%M%S")
        out_dir = os.path.join(base, f"gene_counter_results-{stamp}")
        os.makedirs(out_dir, exist_ok=True)
        return out_dir  # :contentReference[oaicite:33]{index=33}

    def write_pipeline_config(self, out_dir: str):
        cfg_path = os.path.join(out_dir, "pipeline_config.json")
        data = {
            "global_csv": os.path.join(out_dir, "global_counts.csv"),
            "per_file_csv": os.path.join(out_dir, "per_file_counts.csv"),
            "user_genes": self.parse_user_genes(),
            "use_hgnc": True,
            "hgnc_tsv": str(self.hgnc_local),
        }
        try:
            with open(cfg_path, "w", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, indent=2)
            self.append_log(f"pipeline_config.json 저장: {cfg_path}")
        except Exception as e:
            QMessageBox.warning(self, "경고", f"pipeline_config.json 저장 중 오류: {e}")  # :contentReference[oaicite:34]{index=34}

    def on_run(self):
        if not self.validate_inputs(): return
        self.log.clear(); self.append_log("분석을 시작합니다…")
        self.set_controls_enabled(False)
        self.progress.show(); self.progress.setRange(0, 0)
        self.status.showMessage("실행 중…")
        out_dir = self.compute_out_dir(); self.last_out_dir = out_dir
        self.append_log(f"출력 폴더: {out_dir}")
        try:
            cmd = self.build_command(out_dir)
        except Exception as e:
            self.progress.hide(); self.set_controls_enabled(True)
            QMessageBox.critical(self, "실행 중단", str(e)); return
        pretty = " ".join([f'"{c}"' if " " in c else c for c in cmd])
        self.append_log(f"실행 커맨드:\n  {pretty}")
        self.process.start(cmd[0], cmd[1:])
        if not self.process.waitForStarted(3000):
            self.progress.hide(); self.set_controls_enabled(True)
            QMessageBox.critical(self, "오류", "프로세스 시작 실패"); return  # :contentReference[oaicite:35]{index=35}

    def on_read_stdout(self):
        s = bytes(self.process.readAllStandardOutput()).decode(errors="ignore")
        if s: self.append_log(s)  # :contentReference[oaicite:36]{index=36}

    def on_read_stderr(self):
        s = bytes(self.process.readAllStandardError()).decode(errors="ignore")
        if s: self.append_log(s)  # :contentReference[oaicite:37]{index=37}

    def on_process_finished(self, exit_code: int, exit_status):
        self.append_log(f"\n프로세스 종료 코드: {exit_code}")
        out_dir = self.last_out_dir
        if out_dir: self.write_pipeline_config(out_dir)
        self.set_controls_enabled(True); self.progress.hide()
        if exit_code == 0:
            # 필요 시 Step2 외부 스크립트 자동 실행 (동일 폴더에 있을 때)
            if out_dir and self.chk_autostep2.isChecked():
                step2_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "step2_grn_selector.py")
                cfg_path = os.path.join(out_dir, "pipeline_config.json")
                if os.path.isfile(step2_path) and os.path.isfile(cfg_path):
                    QProcess.startDetached(sys.executable, [step2_path, cfg_path])
                    self.append_log(f'Step 2 실행: "{step2_path}" "{cfg_path}"')
            if self.chk_open_folder.isChecked() and out_dir:
                self.open_folder(out_dir)
            QMessageBox.information(self, "완료", f"분석 완료\n출력 폴더: {out_dir}")
        else:
            QMessageBox.warning(self, "실행 오류", "실행 중 오류")  # :contentReference[oaicite:38]{index=38}

    def open_folder(self, path: str):
        try:
            if sys.platform.startswith("win"): os.startfile(path)
            elif sys.platform == "darwin": os.system(f'open "{path}"')
            else: os.system(f'xdg-open "{path}"')
        except Exception: pass

    def closeEvent(self, e):
        self.save_settings(); super().closeEvent(e)

# -------- Step2 --------
class Step2Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("STEP 2 — GRN 유전자 선택")
        self.resize(1000, 720)
        QApplication.setStyle("Fusion")

        self.config_path = ""
        self.config_dir = ""
        self.global_csv = ""
        self.per_file_csv = ""
        self.user_genes: List[str] = []
        self.extracted_counts: Dict[str, int] = {}  # :contentReference[oaicite:39]{index=39}

        header = Card()
        t = QLabel("GRN 구성 — Step 2")
        tf = QFont(); tf.setPointSize(14); tf.setBold(True); t.setFont(tf)
        s = QLabel("Step 1의 pipeline_config.json을 불러오고, 포함할 유전자를 체크한 뒤 'GRN 생성'을 누르세요.")
        header.v.addWidget(t); header.v.addWidget(s)

        self.load_card = Card("설정 불러오기 (pipeline_config.json)")
        row = QHBoxLayout()
        self.le_config = QLineEdit(); self.le_config.setPlaceholderText("pipeline_config.json 경로… (여기로 드래그 가능)")
        self.le_config.setDragEnabled(True)
        self.btn_browse = QToolButton(); self.btn_browse.setText("찾기")
        self.btn_browse.setIcon(self.style().standardIcon(QStyle.SP_DialogOpenButton))
        self.btn_load = QPushButton("불러오기")
        row.addWidget(self.le_config); row.addWidget(self.btn_browse); row.addWidget(self.btn_load)
        self.load_card.v.addLayout(row)
        self.lb_summary = QLabel("-"); self.load_card.v.addWidget(self.lb_summary)

        left = Card("추출된 유전자 (global_counts.csv)")
        left_tools = QHBoxLayout()
        self.ed_filter_left = QLineEdit(); self.ed_filter_left.setPlaceholderText("검색…")
        self.sp_view_topn = QSpinBox(); self.sp_view_topn.setRange(1, 100000); self.sp_view_topn.setValue(200)
        self.sp_topn = QSpinBox(); self.sp_topn.setRange(1, 100000); self.sp_topn.setValue(100)
        self.btn_topn = QToolButton(); self.btn_topn.setText("상위 N 선택")
        self.btn_all_left = QToolButton(); self.btn_all_left.setText("모두 선택")
        self.btn_none_left = QToolButton(); self.btn_none_left.setText("모두 해제")
        left_tools.addWidget(self.ed_filter_left)
        left_tools.addWidget(QLabel("보기 N:")); left_tools.addWidget(self.sp_view_topn)
        left_tools.addWidget(QLabel("N:")); left_tools.addWidget(self.sp_topn)
        left_tools.addWidget(self.btn_topn); left_tools.addStretch(1)
        left_tools.addWidget(self.btn_all_left); left_tools.addWidget(self.btn_none_left)
        self.list_left = CheckList()
        left.v.addLayout(left_tools); left.v.addWidget(self.list_left)  # :contentReference[oaicite:40]{index=40}

        right = Card("사용자 입력 유전자 (Step 1)")
        right_tools = QHBoxLayout()
        self.ed_filter_right = QLineEdit(); self.ed_filter_right.setPlaceholderText("검색…")
        self.btn_all_right = QToolButton(); self.btn_all_right.setText("모두 선택")
        self.btn_none_right = QToolButton(); self.btn_none_right.setText("모두 해제")
        right_tools.addWidget(self.ed_filter_right); right_tools.addStretch(1)
        right_tools.addWidget(self.btn_all_right); right_tools.addWidget(self.btn_none_right)
        self.list_right = CheckList()
        right.v.addLayout(right_tools); right.v.addWidget(self.list_right)

        mid = QSplitter(Qt.Horizontal); mid.addWidget(left); mid.addWidget(right)
        mid.setStretchFactor(0,3); mid.setStretchFactor(1,2)

        bottom = Card("선택 미리보기")
        self.preview = QTextEdit(); self.preview.setReadOnly(True)
        bottom.v.addWidget(self.preview)

        actions = QHBoxLayout()
        self.chk_open = QCheckBox("완료 후 폴더 열기"); self.chk_open.setChecked(True)
        self.btn_generate = QPushButton("GRN 생성")
        actions.addWidget(self.chk_open); actions.addStretch(1); actions.addWidget(self.btn_generate)

        upper = QWidget(); uv = QVBoxLayout(upper)
        uv.setContentsMargins(0,0,0,0); uv.setSpacing(8)
        uv.addWidget(header); uv.addWidget(self.load_card); uv.addWidget(mid); uv.addWidget(bottom); uv.addLayout(actions)

        self.log_card = Card("실행 로그")
        self.log = QTextEdit(); self.log.setReadOnly(True)
        self.log_card.v.addWidget(self.log)

        splitter = QSplitter(Qt.Vertical); splitter.addWidget(upper); splitter.addWidget(self.log_card)
        splitter.setStretchFactor(0, 3); splitter.setStretchFactor(1, 2)

        central = QWidget(); cv = QVBoxLayout(central); cv.addWidget(splitter)
        self.setCentralWidget(central)

        self.status = QStatusBar(); self.setStatusBar(self.status)

        self.btn_browse.clicked.connect(self.on_browse)
        self.btn_load.clicked.connect(self.on_load)
        self.ed_filter_left.textChanged.connect(self.list_left.set_filter)
        self.ed_filter_right.textChanged.connect(self.list_right.set_filter)
        self.btn_all_left.clicked.connect(lambda: self.list_left.select_all(True))
        self.btn_none_left.clicked.connect(lambda: self.list_left.select_all(False))
        self.btn_topn.clicked.connect(self.on_select_topn)
        self.sp_view_topn.valueChanged.connect(lambda _: self._populate_lists())
        self.btn_all_right.clicked.connect(lambda: self.list_right.select_all(True))
        self.btn_none_right.clicked.connect(lambda: self.list_right.select_all(False))
        self.btn_generate.clicked.connect(self.on_generate)
        self.list_left.itemChanged.connect(self.list_left.on_item_changed)
        self.list_right.itemChanged.connect(self.list_right.on_item_changed)

        self.le_config.setAcceptDrops(True); self.le_config.installEventFilter(self)

        if len(sys.argv) > 1:
            arg = sys.argv[1]
            if os.path.isfile(arg):
                self.le_config.setText(arg)
                self.load_config(arg)
                self._hide_loader()  # :contentReference[oaicite:41]{index=41}

    # --- 이벤트/도우미 ---
    def _hide_loader(self):
        try:
            self.le_config.setEnabled(False); self.btn_browse.setEnabled(False)
            self.btn_load.setEnabled(False); self.load_card.setVisible(False)
        except Exception: pass  # :contentReference[oaicite:42]{index=42}

    def eventFilter(self, obj, event):
        if obj is self.le_config:
            if event.type() == QEvent.DragEnter and event.mimeData().hasUrls():
                event.acceptProposedAction(); return True
            elif event.type() == QEvent.Drop:
                urls = event.mimeData().urls()
                if urls:
                    p = urls[0].toLocalFile()
                    if p and os.path.isfile(p): self.le_config.setText(p)
                return True
        return super().eventFilter(obj, event)  # :contentReference[oaicite:43]{index=43}

    def log_append(self, s: str):
        self.log.append(s.rstrip())  # :contentReference[oaicite:44]{index=44}

    def on_browse(self):
        start = os.path.dirname(self.le_config.text()) if self.le_config.text() else os.path.expanduser("~")
        p, _ = QFileDialog.getOpenFileName(self, "pipeline_config.json 선택", start, "JSON (*.json)")
        if p: self.le_config.setText(p)  # :contentReference[oaicite:45]{index=45}

    def on_load(self):
        p = self.le_config.text().strip()
        if not p: return QMessageBox.warning(self, "확인", "pipeline_config.json 경로를 입력하세요.")
        self.load_config(p)  # :contentReference[oaicite:46]{index=46}

    def load_config(self, path: str):
        if not os.path.isfile(path): return QMessageBox.critical(self, "오류", "파일을 찾을 수 없습니다.")
        try:
            with open(path, "r", encoding="utf-8") as f:
                cfg = json.load(f)
        except Exception as e:
            return QMessageBox.critical(self, "오류", f"JSON 읽기 실패: {e}")
        self.config_path = path; self.config_dir = os.path.dirname(path)
        self.global_csv = cfg.get("global_csv", ""); self.per_file_csv = cfg.get("per_file_csv", "")
        self.user_genes = list(dict.fromkeys(cfg.get("user_genes", [])))
        self.lb_summary.setText(f"global_csv: {self.global_csv}\nper_file_csv: {self.per_file_csv}\nuser_genes: {len(self.user_genes)}개")
        self.extracted_counts = self._read_global_counts(self.global_csv)
        self._populate_lists(); self._update_preview()
        self.log_append(f"불러오기 완료: {os.path.basename(path)}")  # :contentReference[oaicite:47]{index=47} :contentReference[oaicite:48]{index=48}

    def _read_global_counts(self, csv_path: str) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        if not csv_path or not os.path.isfile(csv_path):
            self.log_append("global_counts.csv 경로가 유효하지 않습니다."); return counts
        try:
            with open(csv_path, "r", encoding="utf-8") as f:
                rows = list(csv.reader(f))
        except Exception as e:
            self.log_append(f"CSV 읽기 실패: {e}"); return counts
        if not rows: return counts
        header = rows[0]
        data_rows = rows[1:] if any(not c.isdigit() for c in header[1:]) else rows
        for r in data_rows:
            if not r: continue
            gene = r[0].strip() if len(r) > 0 else ""
            cnt_str = (r[1].strip() if len(r) > 1 else "0")
            try: cnt = int(cnt_str)
            except ValueError: cnt = 0
            if gene: counts[gene] = counts.get(gene, 0) + cnt
        return counts  # :contentReference[oaicite:49]{index=49}

    def _populate_lists(self):
        pairs_all = sorted(self.extracted_counts.items(), key=lambda x: (-x[1], x[0]))
        topn = self.sp_view_topn.value() if hasattr(self, "sp_view_topn") else None
        if topn: pairs_all = pairs_all[:topn]
        left_pairs = [(g, f"{g}  ({c})") for g, c in pairs_all]
        self.list_left.clear_all(); self.list_left.set_data(left_pairs)
        right_pairs = [(g, g) for g in self.user_genes]
        self.list_right.clear_all(); self.list_right.set_data(right_pairs)  # :contentReference[oaicite:50]{index=50}

    def _update_preview(self):
        left_sel = set(self.list_left.checked_genes())
        right_sel = set(self.list_right.checked_genes())
        merged = list(dict.fromkeys(list(left_sel) + list(right_sel)))
        self.preview.setPlainText("\n".join(merged))
        self.status.showMessage(f"선택: 추출 {len(left_sel)}개, 사용자 {len(right_sel)}개, 합계 {len(merged)}개")  # :contentReference[oaicite:51]{index=51}

    def on_select_topn(self):
        n = self.sp_topn.value()
        genes = [g for g, _ in self.list_left._items][:n]
        self.list_left.set_checked_genes(genes, True)
        self._update_preview()  # :contentReference[oaicite:52]{index=52}

    def on_generate(self):
        left = set(self.list_left.checked_genes())
        right = set(self.list_right.checked_genes())
        selected = list(dict.fromkeys(list(left) + list(right)))
        if not selected: return QMessageBox.warning(self, "확인", "선택된 유전자가 없습니다.")
        if not self.config_dir: return QMessageBox.critical(self, "오류", "저장 경로를 결정할 수 없습니다.")
        out_txt = os.path.join(self.config_dir, "selected_genes.txt")
        out_json = os.path.join(self.config_dir, "grn_input.json")
        try:
            with open(out_txt, "w", encoding="utf-8") as f: f.write("\n".join(selected))
            with open(out_json, "w", encoding="utf-8") as f:
                json.dump({
                    "selected_genes": selected,
                    "source_config": self.config_path,
                    "global_csv": self.global_csv,
                    "per_file_csv": self.per_file_csv,
                }, f, ensure_ascii=False, indent=2)
            self.log_append(f"선택 저장: {out_txt}")
            self.log_append(f"구성 저장: {out_json}")
        except Exception as e:
            return QMessageBox.critical(self, "오류", f"저장 실패: {e}")
        self.preview.setPlainText("\n".join(selected))
        builder = os.path.join(os.path.dirname(self.config_path), "build_grn.py")
        if os.path.isfile(builder):
            run = QMessageBox.question(self, "확인", "build_grn.py가 발견되었습니다. 지금 실행할까요?",
                                       QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if run == QMessageBox.Yes:
                self._run_builder(builder, out_json); return
        QMessageBox.information(self, "완료", "GRN 입력 파일이 저장되었습니다.")
        if self.chk_open.isChecked(): self._open_folder(self.config_dir)  # :contentReference[oaicite:53]{index=53} :contentReference[oaicite:54]{index=54}

    def _run_builder(self, builder_script: str, json_path: str):
        self.log_append(f"실행: {builder_script}")
        self.proc = QProcess(self)
        self.proc.readyReadStandardOutput.connect(lambda: self._on_proc_out(self.proc))
        self.proc.readyReadStandardError.connect(lambda: self._on_proc_err(self.proc))
        self.proc.finished.connect(lambda code, status: self._on_proc_done(code))
        cmd = [sys.executable, builder_script, json_path]
        pretty = " ".join([f'"{c}"' if " " in c else c for c in cmd])
        self.log_append(f"커맨드:\n  {pretty}")
        self.proc.start(cmd[0], cmd[1:])  # :contentReference[oaicite:55]{index=55}

    def _on_proc_out(self, p: QProcess):
        s = bytes(p.readAllStandardOutput()).decode(errors="ignore")
        if s: self.log_append(s)  # :contentReference[oaicite:56]{index=56}

    def _on_proc_err(self, p: QProcess):
        s = bytes(p.readAllStandardError()).decode(errors="ignore")
        if s: self.log_append(s)  # :contentReference[oaicite:57]{index=57}

    def _on_proc_done(self, code: int):
        self.log_append(f"빌더 종료 코드: {code}")
        if code == 0:
            QMessageBox.information(self, "완료", "GRN 생성이 완료되었습니다.")
            if self.chk_open.isChecked() and self.config_dir:
                self._open_folder(self.config_dir)
        else:
            QMessageBox.warning(self, "오류", "GRN 생성 중 오류가 발생했습니다. 로그를 확인하세요.")  # :contentReference[oaicite:58]{index=58}

    def _open_folder(self, path: str):
        try:
            if sys.platform.startswith('win'): os.startfile(path)
            elif sys.platform == 'darwin': os.system(f'open "{path}"')
            else: os.system(f'xdg-open "{path}"')
        except Exception: pass  # :contentReference[oaicite:59]{index=59}