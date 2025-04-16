# -*- coding: utf-8 -*-
# Copyright (C) 2025 Guangxzhou Molcalx Ltd.
# Released under CC-BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/).

from cresset import flare
from rdkit import Chem
from PySide2 import QtWidgets
from langchain_ollama import OllamaLLM
import re

def predict_herg(smiles):
    """
    Predicts hERG toxicity of compounds using the txgemma model from Ollama.
    Args:
      smiles: SMILES string of the compound
    Returns:
      '0' or '1' indicating doesn't block or blocks hERG respectively,
      returns "N/A" if parsing fails
    """
    host = "192.168.31.197"
    port = "11434"
    llm = OllamaLLM(
        base_url=f"http://{host}:{port}",
        model="txgemma-predict:27b",
        temperature=0
    )
    task = (
        "\n Instructions: Answer the following question about drug properties. "
        "\n Context: Human ether-à-go-go related gene (hERG) is crucial for "
        "the coordination of the heart's beating. Thus, if a drug blocks "
        "the hERG, it could lead to severe adverse effects. Therefore, "
        "reliable prediction of hERG liability in the early stages of "
        "drug design is quite important to reduce the risk of "
        "cardiotoxicity-related attritions in the later development stages."
        "\n Question:  Given a drug SMILES string, predict whether it \n (A) does not inihibit hERG  (B) inhibitr hERG"
        "\n Drug SMILES:"
   )    
    response = llm.invoke(task + " " + smiles + "\n Answer:")
    # print(f"Model response: {response}")
    # Using regex to extract prediction result (0 or 1)
    """
    match = re.search(r'\b([01])\b', response)
    if match:
        return match.group(1)
    else:
        return "N/A"
    """

    # Extract prediction based on actual response format
    if 'A' in response:
        return '0'
    elif 'B' in response:
        return '0'
    else:
        return 'N/A'

@flare.extension
class CalculateHERGExtension:
    """
    Adds a button to the ligand table for predicting hERG toxicity using the txgemma model.
    Prediction results (0/1) will be inserted into a new column "Txgemma_hERG".
    """
    
    def __init__(self):
        parent_widget = flare.main_window().widget()
        self._dialog = QtWidgets.QDialog(parent_widget)
        self._dialog.setWindowTitle("Ligand hERG Prediction")
        
        self._message_label = QtWidgets.QLabel(
            "Click OK to predict hERG toxicity\n for all ligands using txgemma model.", 
            self._dialog
        )
        button_box = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel, 
            self._dialog
        )
        button_box.accepted.connect(self._dialog.accept)
        button_box.rejected.connect(self._dialog.reject)
        
        layout = QtWidgets.QVBoxLayout(self._dialog)
        layout.addWidget(self._message_label)
        layout.addWidget(button_box)
    
    def load(self):
        tab = flare.main_window().ribbon["Extensions"]
        group = tab["Ligand"]
        control = group.add_button("TXGemma hERG", self._show_dialog)
        control.tooltip = (
            "Predict hERG toxicity for all ligands using TXGemma model.\n"
            "The results will be inserted as a new column 'Txgemma_hERG' in the ligand table.\n"
            "0: doesn't block hERG\n"
            "1: blocks hERG"
        )
        print(f"Loaded {self.__class__.__module__}.{self.__class__.__name__}")
    
    def _show_dialog(self):
        if self._dialog.exec_():
            project = flare.main_window().project
            for ligand_obj in project.ligands:
                try:
                    # Convert Flare Ligand object to RDKit molecule
                    lig = flare.Ligand.to_rdmol(ligand_obj)
                except RuntimeError as e:
                    print(f"Error converting ligand: {e}")
                    continue  # Skip ligand on error

                # Generate SMILES string from molecule
                smiles = Chem.MolToSmiles(lig)
                cur_mol = Chem.MolFromSmiles(smiles)
                smiles = Chem.MolToSmiles(cur_mol)
                # Predict hERG toxicity using txgemma model
                herg_prediction = predict_herg(smiles)
                
                # Update molecular properties
                ligand_obj.properties['Txgemma_hERG'].value = herg_prediction
                
            # Show completion message
            QtWidgets.QMessageBox.information(
                flare.main_window().widget(), "Info", "Ligand hERG prediction completed"
            )
