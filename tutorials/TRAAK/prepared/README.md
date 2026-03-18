## Tutorial

In this section, we will utilize TRAAK as an example to illustrate the capability of the Multiple-basin Gō-Martini method in sampling conformational transitions of membrane proteins. In brief, we first generate the multiple-basin Gō-Martini model of up- and down-state TRAAK. Subsequently, we perform a coarse-grained MD simulation with OpenMM. For convenience, all the necessary files used within this tutorial are supplied.

(1) Preprocess PDB files (PDB code: 7LJB and 4WFE) for the up- and down-state TRAAK. Mutations and missing residues have been repaired by Modeler and residues of proteins should be trimmed to have the same length. In this case, we reserve residues 28-286 for each chain of TRAAK.

(2) Generate the Multiple-basin Gō-Martini model for TRAAK. Firstly, we need to upload the atomistic pdb files of TRAAK to [web-server](http://info.ifpan.edu.pl/~rcsu/rcsu/index.html) and generate the contact maps with default settings. Download and uncompress the generated .tgz files. Then, we utilize `ctgomartinize` to generate the Multiple-basin Martini topology of TRAAK (TRAAK.itp and TRAAK_params.itp). If the web server is not available, we can also use the script [OVrCSU.py](https://github.com/ComputBiophys/OVrCSU) to generate the corresponding contact maps of TRAAK.

```bash
ctgomartinize -s 7LJB_clean.pdb 4WFE_clean.pdb -m 7LJB_clean.map 4WFE_clean.map -mol Up Down -mbmol TRAAK -dssp path/to/dssp -ff martini3001 -method exp -other_params "-nt -cys 0.25"
```

Then we should add the suitable mixing parameters for the Multiple-basin Gō-Martini model of TRAAK.
```bash
vi TRAAK.itp # Change beta, C1, and C2 to 1/215, -30, and 0, respectively.
```

(3) Insert the CG protein of the down-state TRAAK into a POPC bilayer (Lipid counts: Upper leaflet = 146, Lower leaflet = 152) and solvate the system with water and ions by using the script [insane.py](https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/tools/proteins_and_bilayers/insane.py). The lipid counts were derived from CHARMM-GUI based on the all-atom model. Be aware that there are issues with ion names and counts for the added ions that require manual correction.

```bash
python2 insane.py -f Down/Down_cg.pdb -box 11,11,12 -o ions.gro -salt 0.15 -l POPC -charge auto -center -sol W -dm 0.5 2>system.top

# Repair the wrong ion names and ion counts.
sed -i "s/NA+    NA+/NA      NA/g" ions.gro
sed -i "s/CL-    CL-/CL      CL/g" ions.gro
sed -i "s/NA+/NA /g" system.top
sed -i "s/CL-/CL /g" system.top
vi system.top # Delete 12 NA in system.top and ions.gro. The error arises because the virtual atoms without charges are also recognized as regular charged residues by insane.py.
vi system.top # Modify the upper and lower lipid counts to 146 and 152 in system.top.
vi ions.gro # Delete 12 NA and modify the lipid counts
```

(4) Prepare the following files: system.top and the MD parameter files (npt.inp and md.inp). 

```bash
vi system.top # Add the martini 3.0 force field and the TRAAK topology files
```

If everything goes well, our next step is to run the simulations as follows.

```bash
# Minimization and Equilibration
run_ctgomartini -i npt.inp

# Production
run_ctgomartini -i md.inp
```

(5) Finally, analyze the simulations we just obtained.

```bash
python -m ctgomartini.analysis.drms_analysis -s npt.pdb -f md.xtc -r Up/Up_cg.pdb Down/Down_cg.pdb -sel "name BB" -prefix dRMStrj
```
