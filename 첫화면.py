import sys
import os
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, 
                             QHBoxLayout, QPushButton, QLabel, QLineEdit, 
                             QTextEdit, QListWidget, QFileDialog, QMessageBox,
                             QGroupBox, QScrollArea, QSplitter)
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QFont, QPixmap, QDragEnterEvent, QDropEvent


class PdfDropArea(QWidget):
    """PDF 파일을 드래그 앤 드롭할 수 있는 영역"""
    files_dropped = pyqtSignal(list)
    
    def __init__(self):
        super().__init__()
        self.setAcceptDrops(True)
        self.setMinimumHeight(200)
        self.setStyleSheet("""
            QWidget {
                border: 2px dashed #aaa;
                border-radius: 10px;
                background-color: #f9f9f9;
                color: #666;
            }
            QWidget:hover {
                border-color: #007acc;
                background-color: #f0f8ff;
            }
        """)
        
        layout = QVBoxLayout()
        self.label = QLabel("PDF 파일을 여기에 드래그하거나\n아래 버튼을 클릭하세요")
        self.label.setAlignment(Qt.AlignCenter)
        self.label.setFont(QFont("Arial", 12))
        layout.addWidget(self.label)
        self.setLayout(layout)
    
    def dragEnterEvent(self, event: QDragEnterEvent):
        if event.mimeData().hasUrls():
            event.accept()
        else:
            event.ignore()
    
    def dropEvent(self, event: QDropEvent):
        files = []
        for url in event.mimeData().urls():
            file_path = url.toLocalFile()
            if file_path.lower().endswith('.pdf'):
                files.append(file_path)
        
        if files:
            self.files_dropped.emit(files)
        else:
            QMessageBox.warning(self, "경고", "PDF 파일만 업로드할 수 있습니다.")


class ResearchToolWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.uploaded_files = []
        self.gene_list = []
        self.init_ui()
        
    def init_ui(self):
        # 메인 윈도우 설정
        self.setWindowTitle('연구자용 논문 분석 도구')
        self.setGeometry(100, 100, 1000, 700)
        
        # 중앙 위젯 설정
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        
        # 메인 레이아웃 (수직 분할)
        main_layout = QHBoxLayout(central_widget)
        
        # 스플리터로 좌우 분할
        splitter = QSplitter(Qt.Horizontal)
        main_layout.addWidget(splitter)
        
        # 왼쪽 패널 - PDF 업로드
        left_panel = self.create_pdf_upload_panel()
        splitter.addWidget(left_panel)
        
        # 오른쪽 패널 - 유전자 입력
        right_panel = self.create_gene_input_panel()
        splitter.addWidget(right_panel)
        
        # 스플리터 비율 설정
        splitter.setSizes([500, 500])
        
        # 상태바
        self.statusBar().showMessage('준비됨')
        
    def create_pdf_upload_panel(self):
        """PDF 업로드 패널 생성"""
        panel = QGroupBox("논문 업로드 (PDF)")
        panel.setFont(QFont("Arial", 12, QFont.Bold))
        layout = QVBoxLayout(panel)
        
        # 설명 레이블
        desc_label = QLabel("연구하신 논문을 업로드해주세요.")
        desc_label.setFont(QFont("Arial", 10))
        desc_label.setStyleSheet("color: #555; margin-bottom: 10px;")
        layout.addWidget(desc_label)
        
        # 드래그 앤 드롭 영역
        self.drop_area = PdfDropArea()
        self.drop_area.files_dropped.connect(self.handle_dropped_files)
        layout.addWidget(self.drop_area)
        
        # 파일 선택 버튼
        select_btn = QPushButton("파일 선택")
        select_btn.setFont(QFont("Arial", 11))
        select_btn.setStyleSheet("""
            QPushButton {
                background-color: #007acc;
                color: white;
                padding: 10px 20px;
                border: none;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #005a9f;
            }
            QPushButton:pressed {
                background-color: #004080;
            }
        """)
        select_btn.clicked.connect(self.select_pdf_files)
        layout.addWidget(select_btn)
        
        # 업로드된 파일 목록
        list_label = QLabel("업로드된 논문:")
        list_label.setFont(QFont("Arial", 11, QFont.Bold))
        list_label.setStyleSheet("margin-top: 15px; color: #333;")
        layout.addWidget(list_label)
        
        self.file_list = QListWidget()
        self.file_list.setStyleSheet("""
            QListWidget {
                border: 1px solid #ddd;
                border-radius: 5px;
                padding: 5px;
                background-color: white;
            }
            QListWidget::item {
                padding: 8px;
                border-bottom: 1px solid #eee;
            }
            QListWidget::item:selected {
                background-color: #e3f2fd;
                color: #1976d2;
            }
        """)
        layout.addWidget(self.file_list)
        
        # 파일 삭제 버튼
        remove_btn = QPushButton("선택된 파일 삭제")
        remove_btn.setFont(QFont("Arial", 10))
        remove_btn.setStyleSheet("""
            QPushButton {
                background-color: #dc3545;
                color: white;
                padding: 8px 16px;
                border: none;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #c82333;
            }
        """)
        remove_btn.clicked.connect(self.remove_selected_file)
        layout.addWidget(remove_btn)
        
        return panel
    
    def create_gene_input_panel(self):
        """유전자 입력 패널 생성"""
        panel = QGroupBox("관심 유전자")
        panel.setFont(QFont("Arial", 12, QFont.Bold))
        layout = QVBoxLayout(panel)
        
        # 설명 레이블
        desc_label = QLabel("분석하고자 하는 유전자명을 입력해주세요.")
        desc_label.setFont(QFont("Arial", 10))
        desc_label.setStyleSheet("color: #555; margin-bottom: 10px;")
        layout.addWidget(desc_label)
        
        # 유전자 입력 영역
        input_layout = QHBoxLayout()
        
        self.gene_input = QLineEdit()
        self.gene_input.setPlaceholderText("유전자명 입력 (예: BRCA1, TP53, EGFR)")
        self.gene_input.setFont(QFont("Arial", 11))
        self.gene_input.setStyleSheet("""
            QLineEdit {
                padding: 10px;
                border: 2px solid #ddd;
                border-radius: 5px;
                font-size: 11px;
            }
            QLineEdit:focus {
                border-color: #007acc;
            }
        """)
        self.gene_input.returnPressed.connect(self.add_gene)
        input_layout.addWidget(self.gene_input)
        
        add_btn = QPushButton("추가")
        add_btn.setFont(QFont("Arial", 11))
        add_btn.setStyleSheet("""
            QPushButton {
                background-color: #28a745;
                color: white;
                padding: 10px 20px;
                border: none;
                border-radius: 5px;
                font-weight: bold;
            }
            QPushButton:hover {
                background-color: #218838;
            }
        """)
        add_btn.clicked.connect(self.add_gene)
        input_layout.addWidget(add_btn)
        
        layout.addLayout(input_layout)
        
        # 유전자 목록
        list_label = QLabel("관심 유전자 목록:")
        list_label.setFont(QFont("Arial", 11, QFont.Bold))
        list_label.setStyleSheet("margin-top: 15px; color: #333;")
        layout.addWidget(list_label)
        
        self.gene_list_widget = QListWidget()
        self.gene_list_widget.setStyleSheet("""
            QListWidget {
                border: 1px solid #ddd;
                border-radius: 5px;
                padding: 5px;
                background-color: white;
            }
            QListWidget::item {
                padding: 8px;
                border-bottom: 1px solid #eee;
                background-color: #e8f5e8;
                margin: 2px;
                border-radius: 3px;
            }
            QListWidget::item:selected {
                background-color: #c8e6c9;
                color: #2e7d32;
            }
        """)
        layout.addWidget(self.gene_list_widget)
        
        # 유전자 삭제 버튼
        remove_gene_btn = QPushButton("선택된 유전자 삭제")
        remove_gene_btn.setFont(QFont("Arial", 10))
        remove_gene_btn.setStyleSheet("""
            QPushButton {
                background-color: #dc3545;
                color: white;
                padding: 8px 16px;
                border: none;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #c82333;
            }
        """)
        remove_gene_btn.clicked.connect(self.remove_selected_gene)
        layout.addWidget(remove_gene_btn)
        
        # 분석 시작 버튼
        analyze_btn = QPushButton("분석 시작")
        analyze_btn.setFont(QFont("Arial", 12, QFont.Bold))
        analyze_btn.setStyleSheet("""
            QPushButton {
                background-color: #6f42c1;
                color: white;
                padding: 15px 30px;
                border: none;
                border-radius: 8px;
                margin-top: 20px;
            }
            QPushButton:hover {
                background-color: #5a2d91;
            }
            QPushButton:pressed {
                background-color: #4c1e7a;
            }
        """)
        analyze_btn.clicked.connect(self.start_analysis)
        layout.addWidget(analyze_btn)
        
        return panel
    
    def select_pdf_files(self):
        """파일 선택 다이얼로그"""
        files, _ = QFileDialog.getOpenFileNames(
            self, 
            'PDF 파일 선택', 
            '', 
            'PDF 파일 (*.pdf)'
        )
        if files:
            self.handle_dropped_files(files)
    
    def handle_dropped_files(self, files):
        """드롭된 파일 처리"""
        added_count = 0
        for file_path in files:
            if file_path not in self.uploaded_files:
                self.uploaded_files.append(file_path)
                file_name = os.path.basename(file_path)
                self.file_list.addItem(file_name)
                added_count += 1
        
        if added_count > 0:
            self.statusBar().showMessage(f'{added_count}개의 PDF 파일이 추가되었습니다.')
            QMessageBox.information(self, "성공", f'{added_count}개의 PDF 파일이 성공적으로 추가되었습니다.')
    
    def remove_selected_file(self):
        """선택된 파일 삭제"""
        current_row = self.file_list.currentRow()
        if current_row >= 0:
            removed_file = self.uploaded_files.pop(current_row)
            self.file_list.takeItem(current_row)
            self.statusBar().showMessage(f'파일이 삭제되었습니다: {os.path.basename(removed_file)}')
        else:
            QMessageBox.warning(self, "경고", "삭제할 파일을 선택해주세요.")
    
    def add_gene(self):
        """유전자 추가"""
        gene_name = self.gene_input.text().strip().upper()
        if gene_name and gene_name not in self.gene_list:
            self.gene_list.append(gene_name)
            self.gene_list_widget.addItem(gene_name)
            self.gene_input.clear()
            self.statusBar().showMessage(f'유전자가 추가되었습니다: {gene_name}')
        elif gene_name in self.gene_list:
            QMessageBox.warning(self, "경고", "이미 추가된 유전자입니다.")
        else:
            QMessageBox.warning(self, "경고", "유전자명을 입력해주세요.")
    
    def remove_selected_gene(self):
        """선택된 유전자 삭제"""
        current_row = self.gene_list_widget.currentRow()
        if current_row >= 0:
            removed_gene = self.gene_list.pop(current_row)
            self.gene_list_widget.takeItem(current_row)
            self.statusBar().showMessage(f'유전자가 삭제되었습니다: {removed_gene}')
        else:
            QMessageBox.warning(self, "경고", "삭제할 유전자를 선택해주세요.")
    
    def start_analysis(self):
        """분석 시작"""
        if not self.uploaded_files:
            QMessageBox.warning(self, "경고", "분석할 PDF 파일을 먼저 업로드해주세요.")
            return
        
        if not self.gene_list:
            QMessageBox.warning(self, "경고", "관심 유전자를 먼저 추가해주세요.")
            return
        
        # 분석 시작 확인
        msg = f"업로드된 논문 {len(self.uploaded_files)}개와 관심 유전자 {len(self.gene_list)}개로 분석을 시작하시겠습니까?"
        reply = QMessageBox.question(self, "분석 시작", msg, 
                                   QMessageBox.Yes | QMessageBox.No)
        
        if reply == QMessageBox.Yes:
            self.statusBar().showMessage("분석을 시작합니다...")
            # 여기에 실제 분석 로직을 구현
            self.run_analysis()
    
    def run_analysis(self):
        """실제 분석 실행 (구현 예정)"""
        analysis_info = f"""
        분석 정보:
        - 업로드된 PDF: {len(self.uploaded_files)}개
        - 관심 유전자: {', '.join(self.gene_list)}
        
        분석이 시작되었습니다!
        """
        QMessageBox.information(self, "분석 시작", analysis_info)
        self.statusBar().showMessage("분석이 진행 중입니다...")


def main():
    app = QApplication(sys.argv)
    
    # 애플리케이션 스타일 설정
    app.setStyle('Fusion')
    
    window = ResearchToolWindow()
    window.show()
    
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()