# -*- coding: utf-8 -*-
# Copyright (C) 2025 Guangxzhou Molcalx Ltd.
# Released under CC-BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/).

from cresset import flare
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFreeSASA
from PySide2 import QtWidgets

def compute_sasa(mol):
    ptable = Chem.GetPeriodicTable()
    radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()]
    sasa = rdFreeSASA.CalcSASA(mol, radii)
    return sasa

def compute_ligand_sasa_pocket(prot, lig):
    complex_mol = Chem.CombineMols(prot, lig)
    complex_h = Chem.AddHs(complex_mol, addCoords=True)
    # Calculate van der Waals radii for each atom in complex_h
    ptable = Chem.GetPeriodicTable()
    radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in complex_h.GetAtoms()]
    # Pass radii parameter when calling CalcSASA
    rdFreeSASA.CalcSASA(complex_h, radii)
    # Assume the ligand is the last fragment
    comp_lig = Chem.GetMolFrags(complex_h, asMols=True, sanitizeFrags=False)[-1]
    lig_sasa_bound = sum([float(a.GetProp("SASA")) for a in comp_lig.GetAtoms()])
    return lig_sasa_bound

def compute_ligand_sasa_pocket_cutoff(prot, lig, cutoff=8):
    lig_conf = lig.GetConformer()
    lig_xyz = lig_conf.GetPositions()

    prot_conf = prot.GetConformer()
    prot_xyz = prot_conf.GetPositions()

    # Calculate minimum distance between each protein atom and all ligand atoms
    r = np.min(np.linalg.norm(prot_xyz[:, np.newaxis, :] - lig_xyz[np.newaxis, :, :], axis=2), axis=1)
    indices = np.argwhere(r > cutoff).flatten()
    
    prot_cut = Chem.RWMol(prot)
    for idx in sorted(indices, reverse=True):
        prot_cut.RemoveAtom(int(idx))
        
    return compute_ligand_sasa_pocket(prot_cut, lig)

@flare.extension
class CalculateSASAExtension:
    """Add a button to calculate SASA metrics for all ligands: 
       SASA_Free, SASA_Bound, and SASA_Buried"""
    
    def __init__(self):
        parent_widget = flare.main_window().widget()
        self._dialog = QtWidgets.QDialog(parent_widget)
        self._dialog.setWindowTitle("Ligand SASA Calculation")
        
        self._message_label = QtWidgets.QLabel("Click OK to calculate SASA metrics for all ligands.", self._dialog)
        button_box = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel, self._dialog
        )
        button_box.accepted.connect(self._dialog.accept)
        button_box.rejected.connect(self._dialog.reject)
        
        layout = QtWidgets.QVBoxLayout(self._dialog)
        layout.addWidget(self._message_label)
        layout.addWidget(button_box)
    
    def load(self):
        tab = flare.main_window().ribbon["Extensions"]
        group = tab["Ligand"]
        control = group.add_button("Ligand SASA", self._show_dialog)
        control.tooltip = "Calculate three SASA values for the ligand based on the ligand and its associating protein: free SASA, bound SASA and buried SASA.\n Three columns will be added into ligand table:\n SASA_free, SASA_bound and SASA_buried."
        print(f"Loaded {self.__class__.__module__}.{self.__class__.__name__}")
    
    def _show_dialog(self):
        if self._dialog.exec_():
            project = flare.main_window().project
            for ligand_obj in project.ligands:
                lig = ligand_obj.to_rdmol()
                prot = ligand_obj.protein.to_rdmol()
                
                sasa_free = round(compute_sasa(lig), 2)
                sasa_bound = round(compute_ligand_sasa_pocket_cutoff(prot, lig, cutoff=8), 2)
                sasa_buried = round(sasa_free - sasa_bound, 2)
                
                ligand_obj.properties['SASA_Free'].value = str(sasa_free)
                ligand_obj.properties['SASA_Bound'].value = str(sasa_bound)
                ligand_obj.properties['SASA_Buried'].value = str(sasa_buried)
            
            # Display message using QMessageBox
            QtWidgets.QMessageBox.information(
                flare.main_window().widget(), "Info", "Ligand SASA calculation completed"
            )
