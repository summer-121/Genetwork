"""
Step 2 GUI (PyQt5): GRN Selector — 체크박스로 유전자 선택 후 GRN 생성
-----------------------------------------------------------------
기능 요약
- Step1에서 생성된 pipeline_config.json을 불러옵니다(버튼/드래그앤드롭/명령줄 인자 지원).
- global_counts.csv 에서 **추출된 유전자 리스트** (카운트 포함)를 표시하고 체크박스로 선택.
- pipeline_config.json 의 **user_genes** (사용자 직접 입력) 리스트도 표시하고 체크박스로 선택.
- 검색 필터, 모두 선택/해제, (추출 리스트는) 상위 N개 빠른 선택 지원.
- 하단 "GRN 생성" 버튼: 선택 유전자를 합쳐서 중복 제거 후
  - 같은 폴더에 selected_genes.txt, grn_input.json 저장
  - (옵션) build_grn.py 가 같은 폴더에 있으면 QProcess로 실행 시도

요구 사항 반영
- PyQt5, 라이트 테마, 디자이너로 수정 가능한 위젯만 사용(커스텀 위젯 최소화)
- 출력 폴더를 따로 묻지 않고 **config가 있는 폴더**에 저장
- 버튼 문구: "GRN 생성"

실행 팁
- 명령줄 인자로 config 경로 전달 가능:
  python step2_grn_selector.py "C:/.../gene_counter_results-YYYYMMDD-HHMMSS/pipeline_config.json"
"""

import csv
import json
import os
import sys
from typing import Dict, List, Tuple

from PyQt5.QtCore import Qt, QProcess
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QFileDialog,
    QListWidget,
    QListWidgetItem,
    QToolButton,
    QPushButton,
    QLabel,
    QVBoxLayout,
    QHBoxLayout,
    QLineEdit,
    QTextEdit,
    QMessageBox,
    QCheckBox,
    QFrame,
    QSplitter,
    QStatusBar,
    QSpinBox,
    QStyle,
)


class Card(QFrame):
    def __init__(self, title: str = "", parent=None):
        super().__init__(parent)
        self.setFrameShape(QFrame.StyledPanel)
        self.setFrameShadow(QFrame.Raised)
        self.v = QVBoxLayout(self)
        self.v.setContentsMargins(12, 12, 12, 12)
        self.v.setSpacing(8)
        if title:
            lbl = QLabel(title)
            f = lbl.font(); f.setBold(True)
            lbl.setFont(f)
            self.v.addWidget(lbl)


class CheckList(QListWidget):
    """체크 가능한 리스트(필터 시에도 체크 상태 유지)."""
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAlternatingRowColors(True)
        # 내부 상태: gene -> checked(bool)
        self._checked = {}
        # 표시용 데이터: [(gene, label_text)]
        self._items: List[Tuple[str, str]] = []
        self._filter = ""

    def clear_all(self):
        self._checked.clear()
        self._items.clear()
        super().clear()

    def set_data(self, pairs: List[Tuple[str, str]]):
        self._items = pairs
        # 기존 체크 상태 유지
        for g, _ in pairs:
            self._checked.setdefault(g, False)
        self._rebuild()

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
            self.addItem(it)

    def set_filter(self, text: str):
        self._filter = text
        self._rebuild()

    def select_all(self, checked: bool):
        # 현재 표시 중인 아이템만 일괄 변경
        for i in range(self.count()):
            it = self.item(i)
            it.setCheckState(Qt.Checked if checked else Qt.Unchecked)
            gene = it.data(Qt.UserRole)
            self._checked[gene] = checked

    def set_checked_genes(self, genes: List[str], checked: bool = True):
        for g in genes:
            self._checked[g] = checked
        self._rebuild()

    def checked_genes(self) -> List[str]:
        return [g for g, v in self._checked.items() if v]

    def all_genes(self) -> List[str]:
        return [g for g, _ in self._items]

    def on_item_changed(self, item: QListWidgetItem):
        gene = item.data(Qt.UserRole)
        self._checked[gene] = (item.checkState() == Qt.Checked)


class Step2Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("STEP 2 — GRN 유전자 선택")
        self.resize(1000, 720)
        QApplication.setStyle("Fusion")

        # 경로 상태
        self.config_path = ""
        self.config_dir = ""
        self.global_csv = ""
        self.per_file_csv = ""
        self.user_genes: List[str] = []
        self.extracted_counts: Dict[str, int] = {}

        # 헤더
        header = Card()
        t = QLabel("GRN 구성 — Step 2")
        tf = QFont(); tf.setPointSize(14); tf.setBold(True); t.setFont(tf)
        s = QLabel("Step 1의 pipeline_config.json을 불러오고, 포함할 유전자를 체크한 뒤 'GRN 생성'을 누르세요.")
        header.v.addWidget(t)
        header.v.addWidget(s)

        # Config 로더
        self.load_card = Card("설정 불러오기 (pipeline_config.json)")
        row = QHBoxLayout()
        self.le_config = QLineEdit(); self.le_config.setPlaceholderText("pipeline_config.json 경로… (여기로 드래그해도 됩니다)")
        self.le_config.setDragEnabled(True)
        self.btn_browse = QToolButton(); self.btn_browse.setText("찾기")
        self.btn_browse.setIcon(self.style().standardIcon(QStyle.SP_DialogOpenButton))
        self.btn_load = QPushButton("불러오기")
        row.addWidget(self.le_config); row.addWidget(self.btn_browse); row.addWidget(self.btn_load)
        self.load_card.v.addLayout(row)
        self.lb_summary = QLabel("-")
        self.load_card.v.addWidget(self.lb_summary)

        # 좌: 추출 유전자, 우: 사용자 유전자
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
        left_tools.addWidget(self.btn_topn)
        left_tools.addStretch(1)
        left_tools.addWidget(self.btn_all_left); left_tools.addWidget(self.btn_none_left)
        self.list_left = CheckList()
        left.v.addLayout(left_tools)
        left.v.addWidget(self.list_left)

        right = Card("사용자 입력 유전자 (Step 1)")
        right_tools = QHBoxLayout()
        self.ed_filter_right = QLineEdit(); self.ed_filter_right.setPlaceholderText("검색…")
        self.btn_all_right = QToolButton(); self.btn_all_right.setText("모두 선택")
        self.btn_none_right = QToolButton(); self.btn_none_right.setText("모두 해제")
        right_tools.addWidget(self.ed_filter_right)
        right_tools.addStretch(1)
        right_tools.addWidget(self.btn_all_right); right_tools.addWidget(self.btn_none_right)
        self.list_right = CheckList()
        right.v.addLayout(right_tools)
        right.v.addWidget(self.list_right)

        # 중앙 스플리터
        mid = QSplitter(Qt.Horizontal)
        mid.addWidget(left)
        mid.addWidget(right)
        mid.setStretchFactor(0, 3)
        mid.setStretchFactor(1, 2)

        # 하단: 미리보기 + 실행
        bottom = Card("선택 미리보기")
        self.preview = QTextEdit(); self.preview.setReadOnly(True)
        bottom.v.addWidget(self.preview)

        actions = QHBoxLayout()
        self.chk_open = QCheckBox("완료 후 폴더 열기"); self.chk_open.setChecked(True)
        self.btn_generate = QPushButton("GRN 생성")
        actions.addWidget(self.chk_open); actions.addStretch(1); actions.addWidget(self.btn_generate)

        # 전체 배치
        upper = QWidget(); uv = QVBoxLayout(upper)
        uv.setContentsMargins(0,0,0,0); uv.setSpacing(8)
        uv.addWidget(header)
        uv.addWidget(self.load_card)
        uv.addWidget(mid)
        uv.addWidget(bottom)
        uv.addLayout(actions)

        self.log_card = Card("실행 로그")
        self.log = QTextEdit(); self.log.setReadOnly(True)
        self.log_card.v.addWidget(self.log)

        splitter = QSplitter(Qt.Vertical)
        splitter.addWidget(upper)
        splitter.addWidget(self.log_card)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)

        central = QWidget(); cv = QVBoxLayout(central)
        cv.setContentsMargins(10,10,10,10); cv.setSpacing(8)
        cv.addWidget(splitter)
        self.setCentralWidget(central)

        self.status = QStatusBar(); self.setStatusBar(self.status)

        # 시그널 연결
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

        # 드래그앤드롭으로 config 경로 입력 지원
        self.le_config.setAcceptDrops(True)
        self.le_config.installEventFilter(self)

        # 시작시 인자 처리
        if len(sys.argv) > 1:
            arg = sys.argv[1]
            if os.path.isfile(arg):
                self.le_config.setText(arg)
                self.load_config(arg)
                self._hide_loader()

    # -------- 이벤트 핸들러 --------
    def _hide_loader(self):
        try:
            self.le_config.setEnabled(False)
            self.btn_browse.setEnabled(False)
            self.btn_load.setEnabled(False)
            self.load_card.setVisible(False)
        except Exception:
            pass

    def eventFilter(self, obj, event):  # 드래그앤드롭 for QLineEdit
        from PyQt5.QtCore import QEvent
        if obj is self.le_config:
            if event.type() == QEvent.DragEnter:
                if event.mimeData().hasUrls():
                    event.acceptProposedAction(); return True
            elif event.type() == QEvent.Drop:
                urls = event.mimeData().urls()
                if urls:
                    p = urls[0].toLocalFile()
                    if p and os.path.isfile(p):
                        self.le_config.setText(p)
                return True
        return super().eventFilter(obj, event)

    def log_append(self, s: str):
        self.log.append(s.rstrip())

    def on_browse(self):
        start = os.path.dirname(self.le_config.text()) if self.le_config.text() else os.path.expanduser("~")
        p, _ = QFileDialog.getOpenFileName(self, "pipeline_config.json 선택", start, "JSON (*.json)")
        if p:
            self.le_config.setText(p)

    def on_load(self):
        p = self.le_config.text().strip()
        if not p:
            return QMessageBox.warning(self, "확인", "pipeline_config.json 경로를 입력하세요.")
        self.load_config(p)

    def load_config(self, path: str):
        if not os.path.isfile(path):
            return QMessageBox.critical(self, "오류", "파일을 찾을 수 없습니다.")
        try:
            with open(path, "r", encoding="utf-8") as f:
                cfg = json.load(f)
        except Exception as e:
            return QMessageBox.critical(self, "오류", f"JSON 읽기 실패: {e}")

        self.config_path = path
        self.config_dir = os.path.dirname(path)
        self.global_csv = cfg.get("global_csv", "")
        self.per_file_csv = cfg.get("per_file_csv", "")
        self.user_genes = list(dict.fromkeys(cfg.get("user_genes", [])))  # 중복 제거, 순서 유지

        # 요약 표시
        self.lb_summary.setText(
            f"global_csv: {self.global_csv}\n"
            f"per_file_csv: {self.per_file_csv}\n"
            f"user_genes: {len(self.user_genes)}개"
        )

        # CSV 로드
        self.extracted_counts = self._read_global_counts(self.global_csv)
        self._populate_lists()
        self._update_preview()
        self.log_append(f"불러오기 완료: {os.path.basename(path)}")

    def _read_global_counts(self, csv_path: str) -> Dict[str, int]:
        counts: Dict[str, int] = {}
        if not csv_path or not os.path.isfile(csv_path):
            self.log_append("global_counts.csv 경로가 유효하지 않습니다.")
            return counts
        try:
            with open(csv_path, "r", encoding="utf-8") as f:
                reader = csv.reader(f)
                rows = list(reader)
        except Exception as e:
            self.log_append(f"CSV 읽기 실패: {e}")
            return counts

        if not rows:
            return counts

        # 헤더 추정: 첫 행에 문자 포함 시 헤더로 간주
        header = rows[0]
        data_rows = rows[1:] if any(not c.isdigit() for c in header[1:]) else rows

        # (gene, count) 형태로 파싱 시도
        for r in data_rows:
            if not r:
                continue
            gene = r[0].strip() if len(r) > 0 else ""
            cnt_str = (r[1].strip() if len(r) > 1 else "0")
            try:
                cnt = int(cnt_str)
            except ValueError:
                # 숫자 변환 실패 시 0 처리
                cnt = 0
            if gene:
                counts[gene] = counts.get(gene, 0) + cnt
        return counts

    def _populate_lists(self):
        # 왼쪽(추출): count 기준 내림차순 + 보기용 Top-N
        pairs_all = sorted(self.extracted_counts.items(), key=lambda x: (-x[1], x[0]))
        topn = self.sp_view_topn.value() if hasattr(self, "sp_view_topn") else None
        if topn:
            pairs_all = pairs_all[:topn]
        left_pairs = [(g, f"{g}  ({c})") for g, c in pairs_all]
        self.list_left.clear_all(); self.list_left.set_data(left_pairs)
        # 오른쪽(사용자)
        right_pairs = [(g, g) for g in self.user_genes]
        self.list_right.clear_all(); self.list_right.set_data(right_pairs)

    def _update_preview(self):
        left_sel = set(self.list_left.checked_genes())
        right_sel = set(self.list_right.checked_genes())
        merged = list(dict.fromkeys(list(left_sel) + list(right_sel)))
        self.preview.setPlainText("\n".join(merged))
        self.status.showMessage(f"선택: 추출 {len(left_sel)}개, 사용자 {len(right_sel)}개, 합계 {len(merged)}개")

    def on_select_topn(self):
        n = self.sp_topn.value()
        genes = [g for g, _ in self.list_left._items][:n]
        self.list_left.set_checked_genes(genes, True)
        self._update_preview()

    def on_generate(self):
        # 선택 취합
        left = set(self.list_left.checked_genes())
        right = set(self.list_right.checked_genes())
        selected = list(dict.fromkeys(list(left) + list(right)))
        if not selected:
            return QMessageBox.warning(self, "확인", "선택된 유전자가 없습니다.")

        if not self.config_dir:
            return QMessageBox.critical(self, "오류", "저장 경로를 결정할 수 없습니다. config를 먼저 불러오세요.")

        # 저장
        out_txt = os.path.join(self.config_dir, "selected_genes.txt")
        out_json = os.path.join(self.config_dir, "grn_input.json")
        try:
            with open(out_txt, "w", encoding="utf-8") as f:
                f.write("\n".join(selected))
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

        # 선택 미리보기 업데이트
        self.preview.setPlainText("\n".join(selected))

        # 선택: build_grn.py 자동 실행(있을 때만)
        builder = os.path.join(os.path.dirname(self.config_path), "build_grn.py")
        if os.path.isfile(builder):
            run = QMessageBox.question(self, "확인", "build_grn.py가 발견되었습니다. 지금 실행할까요?",
                                       QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
            if run == QMessageBox.Yes:
                self._run_builder(builder, out_json)
                return

        QMessageBox.information(self, "완료", "GRN 입력 파일이 저장되었습니다.")
        if self.chk_open.isChecked():
            self._open_folder(self.config_dir)

    def _run_builder(self, builder_script: str, json_path: str):
        self.log_append(f"실행: {builder_script}")
        self.proc = QProcess(self)
        self.proc.readyReadStandardOutput.connect(lambda: self._on_proc_out(self.proc))
        self.proc.readyReadStandardError.connect(lambda: self._on_proc_err(self.proc))
        self.proc.finished.connect(lambda code, status: self._on_proc_done(code))
        cmd = [sys.executable, builder_script, json_path]
        pretty = " ".join([f'"{c}"' if " " in c else c for c in cmd])
        self.log_append(f"커맨드:\n  {pretty}")
        self.proc.start(cmd[0], cmd[1:])

    def _on_proc_out(self, p: QProcess):
        s = bytes(p.readAllStandardOutput()).decode(errors="ignore")
        if s: self.log_append(s)

    def _on_proc_err(self, p: QProcess):
        s = bytes(p.readAllStandardError()).decode(errors="ignore")
        if s: self.log_append(s)

    def _on_proc_done(self, code: int):
        self.log_append(f"빌더 종료 코드: {code}")
        if code == 0:
            QMessageBox.information(self, "완료", "GRN 생성이 완료되었습니다.")
            if self.chk_open.isChecked() and self.config_dir:
                self._open_folder(self.config_dir)
        else:
            QMessageBox.warning(self, "오류", "GRN 생성 중 오류가 발생했습니다. 로그를 확인하세요.")

    # -------- 기타 --------
    def _open_folder(self, path: str):
        try:
            if sys.platform.startswith('win'):
                os.startfile(path)
            elif sys.platform == 'darwin':
                os.system(f'open "{path}"')
            else:
                os.system(f'xdg-open "{path}"')
        except Exception:
            pass


def main():
    app = QApplication(sys.argv)
    w = Step2Window(); w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
