## Tutorial: Replica Exchange MD with Multiple-basin Go-Martini

In this tutorial, we will demonstrate how to perform Replica Exchange Molecular Dynamics (REMD) simulations using the Multiple-basin Gō-Martini method. REMD enhances sampling of conformational transitions by running multiple replicas with different Hamiltonians and allowing exchanges between them.

### Overview

This tutorial uses GlnBP (Glutamine-binding protein) as an example, transitioning between open (1GGG) and closed (1WDN) states. Unlike standard MD, REMD uses:
- Multiple replicas with different energy surfaces (controlled by C1 parameter)
- Exchange attempts between replicas to enhance sampling
- Single-state topologies as unsampled states for analysis

### Prerequisites

- CTGoMartini installed with all dependencies
- OpenMM with CUDA support (recommended)
- Python 2 (for `insane.py`)
- Basic understanding of Multiple-basin Go-Martini (see GlnBP standard tutorial)
- `openmmtools` for REMD (install with: `pip install openmmtools`)

---

### (1) Preprocess PDB Files

The PDB files have been preprocessed for you:
- `1GGG_clean.pdb`: Open-state GlnBP (residues 5-224)
- `1WDN_clean.pdb`: Closed-state GlnBP (residues 5-224)

Water has been removed and residues trimmed to the same length.

---

### (2) Generate Multiple-basin Topology with Single-state Extraction

Generate the multiple-basin topology along with single-state topologies for REMD analysis:

```bash
ctgomartinize -s 1GGG_clean.pdb 1WDN_clean.pdb \
              -m auto \
              -mol Open Closed \
              -mbmol GlnBP \
              -dssp \
              -ff martini3001 \
              -method exp \
              -extract-states
```

**Key differences from standard MD:**
- `-extract-states` flag generates single-state ITP files (`GlnBP_stateA.itp`, `GlnBP_stateB.itp`)
- These are used as unsampled states for computing state populations

**Note:** The topology files have been pre-generated with the following parameters:
- Multiple-basin method: `exp`
- Beta (coupling): `1/300`
- C1: `-300`
- C2: `0`

For REMD, these parameters will be overridden by the replica configuration in `remd.inp`.

---

### (3) Solvate the System

Insert the CG protein into a cubic box and solvate:

```bash
python2 insane.py -f Open/Open_cg.pdb -box 9,9,9 -o ions.gro -salt 0.15 -charge auto -center -sol W 2>system.top
```

Fix ion names and counts:

```bash
# Fix ion names in GRO file
sed -i "s/NA+    NA+/NA      NA/g" ions.gro
sed -i "s/CL-    CL-/CL      CL/g" ions.gro

# Fix ion names in TOP file
sed -i "s/NA+/NA /g" system.top
sed -i "s/CL-/CL /g" system.top
```

**Important:** The `system.top` created by `insane.py` contains 70 NA ions, but 3 of them are excess due to virtual atoms being counted as charged residues. You need to:

1. Edit `system.top` to change NA count from 70 to 67:
```
[ molecules ]
GlnBP     1
W      6109
NA        67    # Changed from 70
CL        64
```

2. Remove the last 3 NA atoms from `ions.gro` (total atoms should be 6972, not 6975).

---

### (4) Prepare System Topology Files

For REMD, you need three system topology files:

**A. Main system.top (multiple-basin)** - Edit to add includes:

Run equilibration:
```bash
run_ctgomartini -i npt.inp
```

This generates `npt.gro` as the starting structure for REMD.

**B. REMD production (remd.inp)**



**Key REMD parameters explained:**
- `replica_count`: Number of replicas (11 in this example)
- `replica_c1`: C1 parameter for each replica (kJ/mol), ranging from -480 to -80
- `replica_c2`: C2 parameter (set to 0 for all replicas)
- `replica_temp`: Temperature (constant 310K for all replicas - Hamiltonian REMD)
- `replica_coupling`: Beta = 1/300
- `remd_unsampled_topfiles`: Single-state topologies for predicting free energy profiles with different mixing parameters

---

### (6) Run REMD Simulation

```bash
# Continue from equilibration
run_ctgomartini -i remd.inp
```

Output files:
- `output.nc` - NetCDF file containing all replica trajectories and exchange information
- `remd.log` - Log file with exchange statistics
- `output_checkpoint.nc` - Checkpoint file for restarting

**Note:** REMD simulations are computationally intensive. The example uses 11 replicas, effectively running 11 parallel simulations with exchanges. For a 1 µs simulation, this is equivalent to 11 µs of aggregate sampling.

---

### (7) Analysis

REMD analysis differs from standard MD. Key analyses include:

**A. Extract replica trajectories:**

```bash
python -m ctgomartini.analysis.remd_replica_state -f output.nc
```

This generates `replica_*.xtc` files for each replica.

**B. Compute state populations using unsampled states (MBAR):**

```bash
python -m ctgomartini.analysis.remd_mbar_analysis \
    -f output.nc \
    -u system_stateA.top system_stateB.top \
    -o populations.dat
```

**C. dRMS analysis (similar to standard MD):**

```bash
python -m ctgomartini.analysis.drms_analysis \
    -s npt.pdb \
    -f replica_0.xtc \
    -r Open/Open_cg.pdb Closed/Closed_cg.pdb \
    -sel "name BB" \
    -prefix dRMStrj
```

**D. Exchange ratio analysis:**

```bash
python -m ctgomartini.analysis.remd_exchange_ratio -f output.nc
```

