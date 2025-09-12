# -*- coding: utf-8 -*-
# Copyright (C) 2025 Guangzhou Molcalx Ltd.
# Released under CC-BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/).
# Originally downloaded from https://gitlab.com/gkxiao
"""
Calculate and display a 2D drug-like properties report for the selected ligand.
The report can be copied to the clipboard.

Ribbon Controls:
    Extensions -> Ligand -> Property Report
"""
from cresset import flare
from PySide2.QtWidgets import (QDialog, QVBoxLayout, QTextEdit, QPushButton,
                               QApplication, QMessageBox)
from PySide2.QtCore import Qt

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski


# ------------------------------------------------------------------------------ #
#  计算函数
# ------------------------------------------------------------------------------ #
def _generate_report(ligand):
    """Return a plain-text property report for a Flare ligand."""
    smi = ligand.smiles()
    mol = ligand.to_rdmol()

    mol_weight = Descriptors.MolWt(mol)
    logp       = Descriptors.MolLogP(mol)
    tpsa       = Descriptors.TPSA(mol)
    hbd        = Lipinski.NumHDonors(mol)
    hba        = Lipinski.NumHAcceptors(mol)
    heavy      = Descriptors.HeavyAtomCount(mol)
    rotb       = Descriptors.NumRotatableBonds(mol)

    report = f"DRUG SMILES: {smi}\n"
    report += f"{'Property':<20} {'Value':<10}\n"
    report += "-" * 30 + "\n"
    report += f"{'Molecular Weight':<20} {mol_weight:.2f}\n"
    report += f"{'logP':<20} {logp:.2f}\n"
    report += f"{'TPSA':<20} {tpsa:.2f}\n"
    report += f"{'Lipinski HBD':<20} {hbd}\n"
    report += f"{'Lipinski HBA':<20} {hba}\n"
    report += f"{'Heavy Atoms':<20} {heavy}\n"
    report += f"{'Rotatable Bonds':<20} {rotb}\n"
    return report


# ------------------------------------------------------------------------------ #
#  GUI 对话框
# ------------------------------------------------------------------------------ #
class _ReportDialog(QDialog):
    def __init__(self, text, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Ligand Property Report")
        self.resize(400, 300)

        layout = QVBoxLayout(self)

        self.editor = QTextEdit()
        self.editor.setPlainText(text)
        self.editor.setReadOnly(True)
        layout.addWidget(self.editor)

        copy_btn = QPushButton("Copy to Clipboard")
        copy_btn.clicked.connect(self._copy)
        layout.addWidget(copy_btn)

    def _copy(self):
        QApplication.clipboard().setText(self.editor.toPlainText())
        QApplication.instance().beep()


# ------------------------------------------------------------------------------ #
#  Flare 扩展入口
# ------------------------------------------------------------------------------ #
@flare.extension
class PropertyReportExtension:
    def load(self):
        """Called by Flare when the extension is loaded."""
        tab = flare.main_window().ribbon["Extensions"]
        tab.tip_key = "X"
        group = tab["Ligand"]

        control = group.add_button("Property Report", self._run)
        control.icon = QIcon("property_report")   # 自动找 property_report.png
        control.tooltip = (
            "Calculate and display 2D drug-like properties for the selected ligand. "
            "The report can be copied to the clipboard."
        )
        control.tip_key = "P"

        print(f"Loaded {self.__class__.__module__}.{self.__class__.__name__}")

    # -------------------------------------------------------------------------- #
    def _run(self):
        """Ribbon button callback."""
        main_window = flare.main_window()
        selected = main_window.selected_ligands

        if not selected:
            QMessageBox.critical(main_window.widget(),
                                 "Error",
                                 "Please select a single ligand first.")
            return
        if len(selected) > 1:
            QMessageBox.warning(main_window.widget(),
                                "Information",
                                "More than one ligand selected; "
                                "report will be generated for the first one only.")

        report_txt = _generate_report(selected[0])
        dlg = _ReportDialog(report_txt, parent=main_window.widget())
        dlg.setAttribute(Qt.WA_DeleteOnClose)
        dlg.show()          # 非模态，不阻塞主界面