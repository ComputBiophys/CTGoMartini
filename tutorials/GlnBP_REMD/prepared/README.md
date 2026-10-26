## Tutorial: Replica Exchange MD with Multiple-basin Gō-Martini

In this tutorial, we will demonstrate how to perform Replica Exchange Molecular Dynamics (REMD) simulations using the Multiple-basin Gō-Martini method. REMD enhances sampling of conformational transitions by running multiple replicas with different Hamiltonians and allowing exchanges between them.

### Overview

This tutorial uses GlnBP (Glutamine-binding protein) as an example, transitioning between open (1GGG) and closed (1WDN) states. Unlike standard MD, REMD uses:
- Multiple replicas with different energy surfaces (controlled by C1 parameter)
- Exchange attempts between replicas to enhance sampling

### Prerequisites

- CTGōMartini installed with all dependencies
- OpenMM >= 8.6 (provides the built-in REMD backend; CUDA support recommended)
- Python 2 (for `insane.py`)
- Basic understanding of Multiple-basin Go-Martini (see GlnBP standard tutorial)

---

### (1) Preprocess PDB Files

The PDB files have been preprocessed for you:
- `1GGG_clean.pdb`: Open-state GlnBP (residues 5-224)
- `1WDN_clean.pdb`: Closed-state GlnBP (residues 5-224)

Water has been removed and residues trimmed to the same length.

---

### (2) Generate Multiple-basin Topology

Generate the multiple-basin topology:

```bash
ctgomartinize -s 1GGG_clean.pdb 1WDN_clean.pdb \
              -m auto \
              -mol Open Closed \
              -mbmol GlnBP \
              -dssp \
              -ff martini3001 \
              -method exp
```

This generates the multiple-basin ITP files (`GlnBP.itp`, `GlnBP_params.itp`) used by all replicas.

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

### (4) Prepare System Topology File

Edit `system.top` (multiple-basin) to add includes:

```
#include "FF/martini_v3.0.0.itp"
#include "GlnBP_params.itp"
#include "GlnBP.itp"
#include "FF/martini_v3.0.0_solvents_v1.itp"
#include "FF/martini_v3.0.0_ions_v1.itp"

[ system ]
Title of the system

[ molecules ]
GlnBP     1
W      6109
NA        67
CL        64
```

---

### (5) Equilibration

Run equilibration:
```bash
vi GlnBP.itp # Change beta, C1, and C2 to 1/300, -300, and 0, respectively.
run_ctgomartini -i npt.inp
```

This generates `npt.gro` as the starting structure for REMD.

---

### (6) Run REMD Simulation

REMD simulations can be run in two modes:

#### A. Single-GPU Mode (for testing or small systems)

```bash
# Continue from equilibration
run_ctgomartini -i remd.inp
```

#### B. Multi-GPU Mode (recommended for production)

The built-in REMD backend can spread replicas across multiple GPUs with the
`remd_devices` input key (one worker Context per replica, integrated
concurrently; replicas are assigned to the listed devices round-robin). For
example, to run 11 replicas on 4 GPUs, add to `remd.inp`:

```
remd_devices = 0 1 2 3
```

Then run as usual:

```bash
run_ctgomartini -i remd.inp
```

Without `remd_devices`, all replicas share a single Context on one device.

---

### (7) REMD Output Files

REMD output is a directory (set by `remd_output`, default `output_remd`):

- `log.csv` - per-iteration state assignment of each replica
- `energy.csv` - reduced energy of every replica in every state (for MBAR)
- `cv_energies.csv` - per-basin energies (E1, E2) of every replica
- `exchange.csv` - cumulative per-pair exchange counters
- `volume.csv` - per-replica box volumes (pressure coupling on)
- `metadata.json` - states, timestep and reporting intervals for the analysis tools
- `checkpoint_<i>.xml` - per-replica checkpoints (used by `--append` restarts)
- `replica_<i>.xtc` - continuous per-replica trajectories (`remd_xtc_output = replica`)
- `state_<i>.xtc` - per-thermodynamic-state trajectories (`remd_xtc_output = state`)

**Note:** REMD simulations are computationally intensive. 

---

### (8) Analysis

REMD simulations generate rich trajectory data that requires specialized analysis tools. CTGōMartini provides several analysis modules for calculating observables and computing free energy surfaces.

#### A. Calculate dRMS (Distance Root Mean Square)

dRMS measures the structural difference from reference states. This is essential for defining the collective variable (CV) used in free energy analysis. It reads the per-replica XTC trajectories (`remd_xtc_output = replica` in `remd.inp`), which carry no topology, so a topology file (`--top`) is required:

```bash
python -m ctgomartini.analysis.remd_drms_analysis \
    -d output_remd \
    --top npt.pdb \
    --num-workers 10 \
    -ref Open/Open_cg.pdb Closed/Closed_cg.pdb \
    -o dRMSTraj
```

#### B. Check Exchange Efficiency

Monitor the exchange acceptance ratio between replicas:

```bash
python -m ctgomartini.analysis.remd_exchange_ratio -f output_remd
```

#### C. Monitor Replica State Population

Track which conformational state each replica visits over time:

```bash
python -m ctgomartini.analysis.remd_replica_state -f output_remd
```

#### D. Free Energy Analysis (MBAR)

The **MBAR (Multistate Bennett Acceptance Ratio)** method computes the free energy surface (FES) along the dRMS collective variable.

**Using the Jupyter Notebook (Recommended):**

Open `REMD_MBAR_Analysis.ipynb` for an interactive analysis workflow:

```bash
jupyter notebook REMD_MBAR_Analysis.ipynb
```

**Notebook Contents:**
1. **Setup** - Import libraries, configure parameters
2. **Initialize** - Load dRMS data, setup FES analyzer
3. **Single State Analysis** - Analyze each replica independently
4. **Parameter Sweep** - Convergence test with different parameters
5. **Mixed State Analysis** - Predict the free energy profile at arbitrary mixing parameters


#### E. Additional Analysis Tools

The per-replica trajectories written with `remd_xtc_output = replica` (or `state`)
are plain XTC files, so they can be fed directly into the standard dRMS
analysis (or any other trajectory tool) without an extraction step:

```bash
# Example: Analyze replica 0
python -m ctgomartini.analysis.drms_analysis \
    -s npt.pdb \
    -f output_remd/replica_0.xtc \
    -r Open/Open_cg.pdb Closed/Closed_cg.pdb \
    -sel "name BB" \
    -n 10 \
    -prefix replica_0_drms
```




