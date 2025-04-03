<h2>pyflare python extension</h2>
<ol>
   <li>Calculate3DSim.py</li>
   <p>Calculates the 3D shape similarity between the selected ligands and all other ligands.</p>
   <li>remoteShapeChEMBL_ext.py</li>
   <p>Perform virtual screening of the ChEMBL35 database on a remote SLURM server, using the selected ligand as the query.</p>
   <li>remoteShapeChemdiv_ext.py</li>
   <p>Perform virtual screening of the ChemDiv database on a remote SLURM server, using the selected ligand as the query.</p>
   <li>remoteShapeRings_ext.py</li>
   <p>Perform shape-based sreening against Peter Ertl's 4M ring systems on a remote SLURM server, using the selected ligand as the query.</p>
   <p>Peter Ertl. (2024) Database of 4 Million Medicinal Chemistry-Relevant Ring Systems. Available at: https://pubs.acs.org/doi/10.1021/acs.jcim.3c01812. </p>
   <li>SASA_ext.py</li>
   <p>Calculate three SASA values for the ligand based on the ligand and its associating protein: free SASA, bound SASA and buried SASA.Three columns will be added into ligand table: SASA_free, SASA_bound and SASA_buried.</p>
   <li>txgemma_hERG_ext.py</li>
   <p>Predicts hERG toxicity of compounds using the txgemma model via Ollama server</p>
</ol>
