# -*- coding: utf-8 -*-
from cresset import flare
import numpy as np
from PySide2 import QtWidgets, QtGui

@flare.extension
class PlaneDefinitionExtension:
    """Extension to add a button for generating atomic plane definition"""
    
    def __init__(self):
        self._dialog = None  # Lazy initialization of the dialog
        
    def load(self):
        """Create the button when the extension is loaded"""
        tab = flare.main_window().ribbon["Extensions"]
        tab.tip_key = "X"
        group = tab["Ligand"]
        
        control = group.add_button("Create AromaticFeat", self._generate_plane_definition)
        control.tooltip = (
            "Generate an aromactic feature based on the selected atoms.\n"
            "At least 3 atoms must be selected.\n"
            "The generated aromatic feature will show up in a dialog box"
        )
        control.tip_key = "P"

    def _create_dialog(self, xml_content):
        """Create a dialog to display the XML result"""
        dialog = QtWidgets.QDialog(flare.main_window().widget())
        dialog.setWindowTitle("Define Aromatic Feature")
        
        # Create the text edit widget
        text_edit = QtWidgets.QTextEdit()
        text_edit.setFont(QtGui.QFont("Courier New"))
        text_edit.setText(xml_content)
        text_edit.setReadOnly(True)
        
        # Create buttons
        copy_btn = QtWidgets.QPushButton("Copy to clipboard")
        copy_btn.clicked.connect(lambda: self._copy_to_clipboard(xml_content))
        close_btn = QtWidgets.QPushButton("Close")
        close_btn.clicked.connect(dialog.close)
        
        # Layout setup
        layout = QtWidgets.QVBoxLayout(dialog)
        layout.addWidget(text_edit)
        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addWidget(copy_btn)
        btn_layout.addWidget(close_btn)
        layout.addLayout(btn_layout)
        
        return dialog

    def _generate_plane_definition(self):
        """Handle button click event"""
        main_window = flare.main_window()
        atoms = main_window.picked_atoms
        
        # Validate atom selection
        if len(atoms) < 3:
            QtWidgets.QMessageBox.critical(
                main_window.widget(),
                "Error",
                f"At least 3 atoms must be selected. The current number of selections is: {len(atoms)}"
            )
            return
            
        # Generate the plane definition XML
        xml = self.generate_plane_xml([a.pos for a in atoms])
        
        # Show the result dialog
        self._dialog = self._create_dialog(xml)
        self._dialog.exec_()

    def _copy_to_clipboard(self, text):
        """Copy text to clipboard"""
        clipboard = QtWidgets.QApplication.clipboard()
        clipboard.setText(text)
        QtWidgets.QMessageBox.information(
            self._dialog,
            "Copyed",
            "Plane has been copied to the clipboard"
        )

    @staticmethod
    def generate_plane_xml(atoms_pos):
        """Core logic to generate plane definition XML (adapted from original code)"""
        points = np.array(atoms_pos)
        centroid = np.mean(points, axis=0)
        centered_points = points - centroid
        
        # Calculate normal vector
        _, _, vh = np.linalg.svd(centered_points)
        normal = vh[2] if vh[2].any() else vh[1]  # Handle degenerate cases
        normal /= np.linalg.norm(normal)
        
        return f'''    <plane name="AR" featureId="1" optional="false" disabled="false" weight="1.0" id="feature1">
      <position x3="{centroid[0]:.3f}" y3="{centroid[1]:.3f}" z3="{centroid[2]:.3f}" tolerance="0.9"/>
      <normal x3="{normal[0]:.7f}" y3="{normal[1]:.7f}" z3="{normal[2]:.7f}" tolerance="0.45"/>
    </plane>'''
