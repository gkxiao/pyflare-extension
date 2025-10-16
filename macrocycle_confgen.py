# -*- coding: utf-8 -*-
# Copyright (C) 2025 Guangzhou Molcalx Ltd.
# Released under CC-BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/).
# Based on RDKit ETKDGv3 + MMFF94s + symmetry-aware RMSD deduplication.
"""
Macrocycle Conformer Search Extension for Flare.
Generates a non-redundant, low-energy 3D conformer ensemble from the selected ligand.
Input 3D coordinates are ignored; only molecular connectivity (via SMILES) is used.
Results are added as a SINGLE ligand with multiple conformers.
"""

import numpy as np
from PySide2 import QtWidgets

from cresset import flare
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom


class MacrocycleConformerSearchDialog(QtWidgets.QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Macrocycle Conformer Search")

        info_label = QtWidgets.QLabel(
            "Generate conformers from the selected ligand's structure.\n"
            "Only molecular connectivity is used; input 3D coordinates are ignored."
        )
        info_label.setWordWrap(True)

        self.num_confs_spin = QtWidgets.QSpinBox()
        self.num_confs_spin.setRange(100, 10000)
        self.num_confs_spin.setValue(2000)

        self.rmsd_spin = QtWidgets.QDoubleSpinBox()
        self.rmsd_spin.setRange(0.1, 2.0)
        self.rmsd_spin.setSingleStep(0.1)
        self.rmsd_spin.setValue(0.5)

        self.energy_window_spin = QtWidgets.QDoubleSpinBox()
        self.energy_window_spin.setRange(1.0, 20.0)
        self.energy_window_spin.setValue(10.0)

        self.energy_tol_spin = QtWidgets.QDoubleSpinBox()
        self.energy_tol_spin.setRange(0.01, 1.0)
        self.energy_tol_spin.setSingleStep(0.01)
        self.energy_tol_spin.setValue(0.05)

        self.max_confs_spin = QtWidgets.QSpinBox()
        self.max_confs_spin.setRange(1, 1000)
        self.max_confs_spin.setValue(1000)

        self.seed_edit = QtWidgets.QLineEdit()
        self.seed_edit.setPlaceholderText("Optional, e.g., 42")

        button_box = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel
        )
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)

        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(info_label)
        form = QtWidgets.QFormLayout()
        form.addRow("Initial Conformers:", self.num_confs_spin)
        form.addRow("IsoRMSD Threshold (Å):", self.rmsd_spin)
        form.addRow("Energy Window (kcal/mol):", self.energy_window_spin)
        form.addRow("Energy Tolerance (kcal/mol):", self.energy_tol_spin)
        form.addRow("Max Output Conformers:", self.max_confs_spin)
        form.addRow("Random Seed (optional):", self.seed_edit)
        layout.addLayout(form)
        layout.addWidget(button_box)
        self.setLayout(layout)

    def get_params(self):
        seed_text = self.seed_edit.text().strip()
        seed = int(seed_text) if seed_text else None
        return {
            "num_confs": self.num_confs_spin.value(),
            "rmsd_threshold": self.rmsd_spin.value(),
            "energy_window": self.energy_window_spin.value(),
            "energy_tolerance": self.energy_tol_spin.value(),
            "max_output_confs": self.max_confs_spin.value(),
            "random_seed": seed,
        }


def calc_isormsd_fast(mol_ref, mol_prb, automorphisms):
    """Fast symmetry-aware RMSD for identical molecules (heavy atoms only)."""
    lig1 = Chem.RemoveHs(mol_ref, sanitize=False)
    lig2 = Chem.RemoveHs(mol_prb, sanitize=False)
    coords1 = lig1.GetConformer().GetPositions()
    coords2 = lig2.GetConformer().GetPositions()
    if coords1.shape[0] == 0 or coords1.shape != coords2.shape:
        return 0.0 if coords1.shape == coords2.shape else 999.0
    min_rmsd = float('inf')
    for match in automorphisms:
        if len(match) != coords1.shape[0]:
            continue
        coords2_mapped = coords2[list(match)]
        diff = coords2_mapped - coords1
        rmsd = np.sqrt(np.mean(np.sum(diff * diff, axis=1)))
        if rmsd < min_rmsd:
            min_rmsd = rmsd
            if min_rmsd < 1e-3:
                break
    return min_rmsd if min_rmsd != float('inf') else 0.0


def generate_conformers_internal(
    rd_mol,
    num_confs=2000,
    rmsd_threshold=0.5,
    energy_window=10.0,
    energy_tolerance=0.05,
    max_output_confs=1000,
    random_seed=None,
):
    mol = Chem.Mol(rd_mol)
    mol.RemoveAllConformers()
    mol = Chem.AddHs(mol)

    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError("MMFF94s parameters not available.")

    mol_heavy = Chem.RemoveHs(mol, sanitize=False)
    automorphisms = list(mol_heavy.GetSubstructMatches(mol_heavy, uniquify=False))
    if not automorphisms:
        automorphisms = [tuple(range(mol_heavy.GetNumAtoms()))]

    params = rdDistGeom.ETKDGv3()
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True
    if random_seed is not None:
        params.randomSeed = random_seed

    cids = rdDistGeom.EmbedMultipleConfs(mol, numConfs=num_confs, params=params)
    if not cids:
        raise RuntimeError("Conformer embedding failed.")

    mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
    energies = []
    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=cid)
        if ff:
            ff.Minimize(maxIts=500)
            energies.append((cid, ff.CalcEnergy()))

    if not energies:
        raise RuntimeError("No valid conformers after optimization.")

    energies.sort(key=lambda x: x[1])
    minE = energies[0][1]

    within = [(cid, e) for cid, e in energies if (e - minE) <= energy_window]

    groups = []
    current = []
    for cid, e in within:
        if not current:
            current.append((cid, e))
        else:
            if abs(e - current[0][1]) <= energy_tolerance:
                current.append((cid, e))
            else:
                groups.append(current)
                current = [(cid, e)]
    if current:
        groups.append(current)

    kept_cids = []
    energy_map = {}
    for group in groups:
        if len(kept_cids) >= max_output_confs:
            break
        rep_cid, rep_eng = group[0]
        kept_cids.append(rep_cid)
        energy_map[rep_cid] = rep_eng
        for cid, e in group[1:]:
            if len(kept_cids) >= max_output_confs:
                break
            duplicate = False
            for ref in kept_cids:
                m1 = Chem.Mol(mol)
                m2 = Chem.Mol(mol)
                m1.RemoveAllConformers()
                m2.RemoveAllConformers()
                m1.AddConformer(mol.GetConformer(ref), assignId=True)
                m2.AddConformer(mol.GetConformer(cid), assignId=True)
                rms = calc_isormsd_fast(m1, m2, automorphisms)
                if rms <= rmsd_threshold:
                    duplicate = True
                    break
            if not duplicate:
                kept_cids.append(cid)
                energy_map[cid] = e

    return mol, kept_cids, energy_map


@flare.extension
class MacrocycleConformerSearchExtension:
    def load(self):
        tab = flare.main_window().ribbon["Extensions"]
        group = tab["Ligand"]
        control = group.add_button("Macrocycle Conf Search", self._show_dialog)
        control.tooltip = (
            "Generate a non-redundant, low-energy conformer ensemble for the selected ligand.\n"
            "Only molecular connectivity is used (via SMILES). Results added as one ligand with multiple conformers."
        )

    def _show_dialog(self):
        main_window = flare.main_window()
        selected = main_window.selected_ligands

        if len(selected) != 1:
            QtWidgets.QMessageBox.critical(
                main_window.widget(),
                "Error",
                "Please select exactly ONE ligand."
            )
            return

        ligand = selected[0]
        try:
            smi = ligand.smiles()
            if not smi.strip():
                raise ValueError("Empty SMILES")
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                main_window.widget(), "Error", f"Failed to get SMILES: {e}"
            )
            return

        dialog = MacrocycleConformerSearchDialog(main_window.widget())
        if dialog.exec_():
            try:
                params = dialog.get_params()
                print("Running macrocycle conformer search...")

                rd_mol = Chem.MolFromSmiles(smi)
                if rd_mol is None:
                    raise ValueError("Invalid SMILES")

                mol_with_h, kept_cids, energy_map = generate_conformers_internal(rd_mol, **params)

                # Build single molecule with multiple conformers (WITH HYDROGENS)
                ensemble_mol = Chem.Mol(mol_with_h)
                ensemble_mol.RemoveAllConformers()

                minE = min(energy_map.values())
                for i, cid in enumerate(kept_cids):
                    conf = mol_with_h.GetConformer(cid)
                    new_conf = Chem.Conformer(conf)
                    ensemble_mol.AddConformer(new_conf, assignId=True)
                    relE = energy_map[cid] - minE
                    ensemble_mol.SetProp(f"Conf{i}_RelEnergy_kcal", f"{relE:.4f}")

                # Use .title for naming
                base_name = ligand.title or "Ligand"
                ensemble_mol.SetProp("_Name", f"{base_name}_MacrocycleEnsemble")

                # Add as ONE ligand with multiple conformers
                project = main_window.project
                project.ligands.append(ensemble_mol)

                num_confs = ensemble_mol.GetNumConformers()
                print(f"✅ Added 1 ligand with {num_confs} non-redundant conformers.")

            except Exception as e:
                QtWidgets.QMessageBox.critical(
                    main_window.widget(), "Error", f"Failed: {e}"
                )
                raise
