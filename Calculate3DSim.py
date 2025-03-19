# -*- coding: utf-8 -*-
# Copyright (C) 2025 Cresset Biomolecular Discovery Ltd.
# Released under CC-BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/).
# Originally downloaded from https://gitlab.com/cresset
"""Calculates the 3D shape similarity of all ligands to the selected ligands.

Ribbon Controls:
    Extensions -> Ligand -> Calculate 3D Sim
        Calculates the 3D shape similarity between the selected ligands and all other ligands.
        A "3D_ShapeSim" column will be added to the Ligands table with the results.
"""
import statistics

from PySide2 import QtWidgets

from cresset import flare

from rdkit import Chem
from rdkit.Chem import AllChem


@flare.extension
class Calculate3DSimExtension:
    """Add a button to the ribbon which calculates the 3D shape similarity."""

    def __init__(self):
        parent_widget = flare.main_window().widget()
        self._dialog = QtWidgets.QDialog(parent_widget)

        button_box = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel, self._dialog
        )

        button_box.accepted.connect(self._dialog.accept)
        button_box.rejected.connect(self._dialog.reject)

        # 3D similarity calculation methods
        self._shape_tanimoto = QtWidgets.QRadioButton("Shape Tanimoto (Default)", self._dialog)
        self._shape_protrude = QtWidgets.QRadioButton("Shape Protrude", self._dialog)

        self._shape_tanimoto.setChecked(True)

        self._message_label = QtWidgets.QLabel("", self._dialog)
        self._message_label.setWordWrap(True)

        layout = QtWidgets.QVBoxLayout(self._dialog)
        layout.addWidget(self._message_label)
        layout.addWidget(self._shape_tanimoto)
        layout.addWidget(self._shape_protrude)
        layout.addStretch(1)
        layout.addWidget(button_box)

    def load(self):
        """Load the extension."""
        tab = flare.main_window().ribbon["Extensions"]
        tab.tip_key = "X"
        group = tab["Ligand"]

        control = group.add_button("Calculate 3D Sim", self._show_dialog)
        control.tooltip = (
            "Calculate the 3D shape similarity of all ligands to the selected ligands.\n"
            + "When selecting multiple ligands, the similarity will be calculated as the average\n"
            + "of the similarity values towards each selected ligand.\n"
            + "A '3_DShapeSim' column will be added to the ligand table with the results."
        )
        control.tip_key = "3"

        print(f"Loaded {self.__class__.__module__}.{self.__class__.__name__}")

    def _show_dialog(self):
        """Show a dialog asking for the 3D similarity method to use."""
        main_window = flare.main_window()
        references = main_window.selected_ligands

        if len(references) == 0:
            QtWidgets.QMessageBox.critical(
                main_window.widget(), "Error", "One or more ligands must be selected."
            )
            return

        project = main_window.project
        message = (
            f"Calculate the 3D shape similarity of all {len(project.ligands)} ligands to the "
            + f"selected {len(references)} reference ligand(s). "
            + 'The results will be placed in the "3D_ShapeSim" column.'
        )
        self._message_label.setText(message)

        if self._dialog.exec_():
            sim_func = None
            if self._shape_tanimoto.isChecked():
                sim_func = self._shape_tanimoto_sim
            elif self._shape_protrude.isChecked():
                sim_func = self._shape_protrude_sim

            self._calculate_3d_sim(references, project.ligands, sim_func)

    @staticmethod
    def _shape_tanimoto_sim(rd_mol1, rd_mol2):
        """Calculate shape similarity using Tanimoto metric.
        The metric ranges from 0 (no overlap) to 1 (perfect overlap)"""
        return 1 - AllChem.ShapeTanimotoDist(rd_mol1, rd_mol2)
    @staticmethod
    def _shape_protrude_sim(rd_mol1, rd_mol2):
        """Calculate shape similarity using protrusion metric.
        Measures volume mismatch with more weight on protruding features"""
        return 1 - AllChem.ShapeProtrudeDist(rd_mol1, rd_mol2)
    @staticmethod
    def _calculate_3d_sim(references, ligands, sim_func):
        """Calculate the 3D shape similarity between the selected ligands and all other ligands."""
        rd_mol_references = [mol.to_rdmol() for mol in references]

        # Validate 3D coordinates presence
        if any(mol.GetConformer().Is3D() == False for mol in rd_mol_references):
            QtWidgets.QMessageBox.critical(
                flare.main_window().widget(),
                "Error",
                "Reference ligands must have valid 3D coordinates"
            )
            return

        for ligand in ligands:
            ligand_mol = ligand.to_rdmol()
            if not ligand_mol.GetConformer().Is3D():
                print(f"Skipping ligand {ligand.name} - missing 3D coordinates")
                continue
                
            sims = [sim_func(ref, ligand_mol) for ref in rd_mol_references]
            sim = statistics.mean(sims)
            ligand.properties["3D_ShapeSim"].value = round(sim, 3)
