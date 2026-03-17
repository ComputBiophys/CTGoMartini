## Tutorial

In this section, we will utilize GlnBP as an example to illustrate the capability of the Multiple-basin Gō-Martini method in sampling conformational transitions. In brief, we first generate the multiple-basin Gō-Martini model of open- and closed-state GlnBP. Subsequently, we perform a coarse-grained MD simulation with OpenMM. For convenience, all the necessary files used within this tutorial are supplied.

(1) Preprocess PDB files (PDB code: 1GGG and 1WDN) for the open- and closed-state GlnBP. Water around the proteins should be removed and residues of proteins should be trimmed to have the same length. In this case, we reserve residues 5-224 of GlnBP.

(2) Generate the Multiple-basin Gō-Martini model for GlnBP. Firstly, we need to upload the atomistic pdb files of GlnBP to [web-server](http://info.ifpan.edu.pl/~rcsu/rcsu/index.html) and generate the contact maps with default settings. Download and uncompress the generated .tgz files. Then, we utilize ctgomartinize.py to generate the Multiple-basin Martini topology of GlnBP (GlnBP.itp and GlnBP_params.itp).

```
python ctgomartinize.py -s 1GGG_clean.pdb 1WDN_clean.pdb -m 1GGG_clean.map 1WDN_clean.map -mol Open Closed -mbmol GlnBP -dssp dssp -ff martini3001 -method exp

```

Then we should add the suitable mixing parameters for the Multiple-basin Gō-Martini model of GlnBP.
```
vi GlnBP.itp # Change beta, C1, and C2 to 1/300, -300, and 0, respectively.

```

(3) Insert the CG protein of the open-state GlnBP into a cubic box and solvate the system with water and ions by using the script [insane.py](http://www.cgmartini.nl/images/tools/insane/insane.py). Be aware that there are issues with ion names and counts for the added ions that require manual correction.

```
python2 insane.py -f Open/Open_cg.pdb -box 9,9,9 -o ions.gro -salt 0.15 -charge auto -center -sol W 2>system.top


# Repair the wrong ion names and ion counts.
sed -i "s/NA+    NA+/NA      NA/g" ions.gro
sed -i "s/CL-    CL-/CL      CL/g" ions.gro
sed -i "s/NA+/NA /g" system.top
sed -i "s/CL-/CL /g" system.top
vi system.top # Delete 3 NA in system.top and ions.gro. The error arises because the virtual atoms without charges are also recognized as regular charged residues by insane.py.
vi ions.gro # Delete 3 NA in ions.gro.

```

(4) Prepare the following files: system.top and the MD parameter files (npt.inp and md.inp). 

```
vi system.top # Add the martini 3.0 force field and the GlnBP topology files
```

If everything goes well, our next step is to run the simulations as follows.

```
# Minimization and Equilibration
python run_ctgomartini.py -i npt.inp

# Production
python run_ctgomartini.py -i md.inp

```

(5) Finally, analyze the simulations we just obtained.