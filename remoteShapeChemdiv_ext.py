# -*- coding: utf-8 -*-
# Copyright (C) 2025 Guangxzhou Molcalx Ltd.
# Released under CC-BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/).

import os
import sys
import tempfile
import hashlib
import paramiko
import statistics
from scp import SCPClient
from datetime import datetime
from PySide2 import QtWidgets
from rdkit import Chem
from cresset import flare

# ======================== SSH Integration ========================
class SSHManager:
    """SSH connection manager"""
    def __init__(self, hostname='192.168.31.197', port=22, 
                 username='gkxiao', key_path='~/.ssh/id_rsa'):
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.ssh.connect(hostname, port, username, 
                         key_filename=os.path.expanduser(key_path))
        self.sftp = self.ssh.open_sftp()
        
    def __enter__(self):
        return (self.ssh, self.sftp)
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.sftp.close()
        self.ssh.close()

class ClusterJobManager:
    """Core job management"""
    @staticmethod
    def resolve_directory_conflict(sftp, path):
        """Handle directory conflicts by appending timestamp"""
        try:
            sftp.stat(path)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M")
            new_path = f"{path}_conflict_{timestamp}"
            sftp.mkdir(new_path)
            return new_path
        except FileNotFoundError:
            sftp.mkdir(path)
            return path

    @staticmethod
    def secure_file_transfer(local_path, ssh, remote_path):
        """
        Secure file transfer with MD5 verification:
        1. Calculate local MD5
        2. Transfer via SCP
        3. Verify remote MD5
        """
        def _generate_md5_fileobj(f_obj):
            hash_md5 = hashlib.md5()
            for chunk in iter(lambda: f_obj.read(4096), b""):
                hash_md5.update(chunk)
            return hash_md5.hexdigest()
        
        # Calculate local MD5
        with open(local_path, "rb") as f:
            local_hash = _generate_md5_fileobj(f)
        
        # SCP transfer
        with SCPClient(ssh.get_transport(), socket_timeout=60) as scp:
            scp.put(local_path, remote_path)
        
        # Remote MD5 verification
        sftp_remote = ssh.open_sftp()
        try:
            with sftp_remote.open(remote_path, "rb") as remote_file:
                remote_hash = _generate_md5_fileobj(remote_file)
        finally:
            sftp_remote.close()
        
        if local_hash != remote_hash:
            raise IOError("File hash verification failed")

    @staticmethod
    def generate_sbatch_script(query_path, work_dir):
        """Generate SLURM submission script"""
        file_name = os.path.basename(query_path)
        file_prefix = os.path.splitext(file_name)[0]
        return f"""#!/bin/bash
#SBATCH --job-name=shapescreen_{file_prefix}
#SBATCH --output={work_dir}/shapescreen_%j.out
#SBATCH --error={work_dir}/shapescreen_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=72:00:00

cd {work_dir} || exit 1
/public/gkxiao/bin/chemdiv_shapescreen.sh {file_name} {file_prefix}_chemdiv_vs.log
"""

# ======================== Flare Extension Implementation ========================
@flare.extension
class RemoteShapeScreenExtension:
    """Flare extension for remote shapescreen submission"""
    def __init__(self):
        parent_widget = flare.main_window().widget()
        self._dialog = QtWidgets.QDialog(parent_widget)
        
        # UI elements
        self._remote_dir_edit = QtWidgets.QLineEdit()
        self._message_label = QtWidgets.QLabel()
        self._status_label = QtWidgets.QLabel("Waiting for submission...")
        
        # Button layout
        button_box = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        button_box.accepted.connect(self._dialog.accept)
        button_box.rejected.connect(self._dialog.reject)

        # Layout setup
        layout = QtWidgets.QVBoxLayout(self._dialog)
        layout.addWidget(QtWidgets.QLabel("Remote working directory:"))
        layout.addWidget(self._remote_dir_edit)
        layout.addWidget(self._message_label)
        layout.addWidget(self._status_label)
        layout.addWidget(button_box)

    def load(self):
        """Initialize extension in Flare's ribbon interface"""
        tab = flare.main_window().ribbon["Extensions"]
        group = tab["Ligand"]
        control = group.add_button("Remote Shapescreen", self._show_dialog)
        control.tooltip = ("Submit selected ligands for remote shapescreen\n"
                           "Requirements:\n"
                           "- At least 1 reference ligand selected\n"
                           "- Ligands must contain 3D coordinates\n"
                           "- Valid SSH key configuration")
        print(f"Loaded {self.__class__.__name__}")

    def _show_dialog(self):
        """Display submission dialog with validation"""
        main_window = flare.main_window()
        ref_ligands = main_window.selected_ligands
        
        # Input validation
        if not ref_ligands:
            QtWidgets.QMessageBox.critical(
                main_window.widget(), "Error", "Please select at least one reference ligand")
            return
            
        if any(not lig.to_rdmol().GetConformer().Is3D() for lig in ref_ligands):
            QtWidgets.QMessageBox.critical(
                main_window.widget(), "Error", "All reference ligands must have 3D coordinates")
            return

        # Update UI information
        self._message_label.setText(
            f"Preparing to submit {len(ref_ligands)} reference ligands\n"
            f"Remote directory will be created with result files")
            
        if self._dialog.exec_():
            remote_dir = self._remote_dir_edit.text().strip()
            if not remote_dir:
                QtWidgets.QMessageBox.critical(
                    main_window.widget(), "Error", "Remote directory path is required")
                return
                
            self._submit_job(ref_ligands, remote_dir)

    def _submit_job(self, ref_ligands, remote_dir):
        """Core job submission logic"""
        try:
            # Generate temporary SDF file
            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tmp:
                sdf_path = tmp.name
                with Chem.SDWriter(sdf_path) as writer:
                    for idx, lig in enumerate(ref_ligands):
                        rd_mol = lig.to_rdmol()
                        rd_mol.SetProp("_Name", f"reference_{idx+1}")
                        writer.write(rd_mol)

            # SSH connection and file transfer
            with SSHManager() as (ssh, sftp):
                # Handle remote directory
                try:
                    sftp.stat(remote_dir)
                    remote_dir = ClusterJobManager.resolve_directory_conflict(sftp, remote_dir)
                except FileNotFoundError:
                    sftp.mkdir(remote_dir)

                # Transfer file
                remote_file = f"{remote_dir}/reference.sdf"
                ClusterJobManager.secure_file_transfer(sdf_path, ssh, remote_file)

                # Generate submission script
                script_content = ClusterJobManager.generate_sbatch_script(remote_file, remote_dir)
                script_path = f"{remote_dir}/submit_job.sh"
                with sftp.file(script_path, 'w') as f:
                    f.write(script_content)
                ssh.exec_command(f"chmod +x {script_path}")

                # Submit job
                _, stdout, _ = ssh.exec_command(f"sbatch {script_path}")
                job_info = stdout.read().decode().strip()
                if not job_info.startswith("Submitted batch job"):
                    raise RuntimeError(f"Job submission failed: {job_info}")

                # Display results
                QtWidgets.QMessageBox.information(
                    flare.main_window().widget(),
                    "Submission Successful",
                    f"Job submitted to: {remote_dir}\nJob ID: {job_info.split()[-1]}",
                    QtWidgets.QMessageBox.Ok
                )

        except paramiko.AuthenticationException:
            QtWidgets.QMessageBox.critical(
                flare.main_window().widget(),
                "SSH Authentication Failed",
                "Please check SSH key configuration\nPath: ~/.ssh/id_rsa"
            )
        except Exception as e:
            QtWidgets.QMessageBox.critical(
                flare.main_window().widget(),
                "Submission Error",
                f"Error details: {str(e)}"
            )
        finally:
            if os.path.exists(sdf_path):
                os.remove(sdf_path)

