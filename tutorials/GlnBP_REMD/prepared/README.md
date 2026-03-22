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
- `openmmtools` for REMD

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
- These are used as basic states for predicting free energies profiles at defined mixing parameters

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

**B. system_stateA.top (unsampled state A)**:

```
#include "FF/martini_v3.0.0.itp"
#include "GlnBP_params.itp"
#include "GlnBP_stateA.itp"
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

**C. system_stateB.top (unsampled state B)**:

```
#include "FF/martini_v3.0.0.itp"
#include "GlnBP_params.itp"
#include "GlnBP_stateB.itp"
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

For running REMD across multiple GPUs (e.g., 11 replicas on 4 GPUs):

**Install MPI dependencies:**

```bash
conda install -c conda-forge mpi4py mpich=3 -y
```

**Generate MPI configuration files:**

The `build_mpirun_configfile.py` script automatically detects your job scheduler (PBS, LSF, or SLURM) and generates the necessary configuration files for MPICH3:

```bash
python build_mpirun_configfile.py "run_ctgomartini -i remd.inp"
```

This creates two files:
- `configfile`: MPI configuration with CUDA device assignments
- `hostfile`: List of hosts for the MPI job

**Run the REMD simulation:**

using `mpiexec.hydra` (MPI launcher):
```bash
mpiexec.hydra -f hostfile -configfile configfile
```

**Note:** The script automatically handles:
- CUDA_VISIBLE_DEVICES assignment for each replica
- Host detection from PBS_GPUFILE, LSB_HOSTS, or SLURM_JOB_NODELIST
- MPICH3-compatible configuration format

---

### (7) REMD Output Files

After the simulation completes, you will have:

- `output.nc` - NetCDF file containing all replica trajectories and exchange information
- `output_real_time_analysis.yaml` - Real-time analysis data 
- `output_checkpoint.nc` - Checkpoint file for restarting simulations

**Note:** REMD simulations are computationally intensive. 

---

### (8) Analysis

REMD simulations generate rich trajectory data that requires specialized analysis tools. CTGoMartini provides several analysis modules for extracting replica trajectories, calculating observables, and computing free energy surfaces.

#### Analysis Workflow Overview

```
┌─────────────────┐    ┌────────────────────┐    ┌──────────────────┐
│  Extract        │───▶│  Calculate dRMS    │───▶│  MBAR Free Energy│
│  Trajectories   │    │  per Replica       │    │  Analysis        │
└─────────────────┘    └────────────────────┘    └──────────────────┘
         │                       │                        │
         ▼                       ▼                        ▼
┌─────────────────┐    ┌────────────────────┐    ┌──────────────────┐
│  Check Exchange │    │  Monitor Replica   │    │  Parameter Sweep │
│  Efficiency     │    │  State Population  │    │  (Optional)      │
└─────────────────┘    └────────────────────┘    └──────────────────┘
```

#### A. Extract Replica Trajectories

Extract individual replica trajectories from the NetCDF output:

```bash
python -m ctgomartini.analysis.remd_trajectory_extractor \
    --mode replica \
    -nc output.nc \
    -c output_checkpoint.nc \
    -p npt.pdb
```

**Output:** `replica_0.xtc`, `replica_1.xtc`, ... (one per replica)

**Alternative - Extract by state:**
```bash
python -m ctgomartini.analysis.remd_trajectory_extractor \
    --mode state \
    -nc output.nc \
    -c output_checkpoint.nc \
    -p npt.pdb
```

#### B. Calculate dRMS (Distance Root Mean Square)

dRMS measures the structural difference from reference states. This is essential for defining the collective variable (CV) used in free energy analysis.

```bash
python -m ctgomartini.analysis.remd_drms_analysis \
    -nc output.nc \
    -c output_checkpoint.nc \
    --num-workers 20 \
    -ref Open/Open_cg.pdb Closed/Closed_cg.pdb \
    -prefix drms_analysis
```

**Key Parameters:**
| Parameter | Description |
|-----------|-------------|
| `-nc` | NetCDF trajectory file from REMD |
| `-c` | Checkpoint file (contains topology info) |
| `-ref` | Reference structures for state A and B |
| `--num-workers` | Parallel processing threads |
| `-prefix` | Output file prefix |

**Output Files:**
- `drms_analysis_stateA.dat` - dRMS to state A (Open) for each frame
- `drms_analysis_stateB.dat` - dRMS to state B (Closed) for each frame
- Format: `Time | Replica_ID | dRMS_A | dRMS_B`

#### C. Check Exchange Efficiency

Monitor the exchange acceptance ratio between replicas:

```bash
python -m ctgomartini.analysis.remd_exchange_ratio -f output.nc
```

**Output:**
- Exchange matrix showing acceptance rates between all replica pairs
- Optimal exchange ratio: 20-40% between adjacent replicas

#### D. Monitor Replica State Population

Track which conformational state each replica visits over time:

```bash
python -m ctgomartini.analysis.remd_replica_state -f output.nc
```

**Output:**
- Time series of state assignments for each replica
- Useful for checking sampling efficiency and transition events

#### E. Free Energy Analysis (MBAR)

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
5. **Mixed State Analysis** - Combine EXP and HAM mixing methods
6. **Results Summary** - Barrier heights, equilibrium constants, FES plots

**Key Configuration in Notebook:**
```python
# CV range for histogram (adjust based on your system)
FES_RANGES = [1.5, 11]  # nm

# Barrier search boundaries within CV range
LEFT_BOUND = 3    # nm (basin 1 region end)
RIGHT_BOUND = 8   # nm (basin 2 region start)
```

**Outputs:**
- Free Energy Surface (FES) plots
- Energy barrier height (kJ/mol)
- Equilibrium constant (Keq)
- Basin free energies

#### F. Additional Analysis Tools

**Standard dRMS analysis (for single trajectory):**
```bash
python -m ctgomartini.analysis.drms_analysis \
    -s npt.gro \
    -f trajectory.xtc \
    -r Open/Open_cg.pdb Closed/Closed_cg.pdb \
    -sel "name BB" \
    -n 10 \
    -prefix analysis
```

**Extract specific frames:**
```bash
python -m ctgomartini.analysis.remd_trajectory_extractor \
    --mode replica \
    -nc output.nc \
    -c output_checkpoint.nc \
    -p npt.pdb \
    -b 1000 -e 5000  # Extract frames 1000-5000
```

---

### Summary of Analysis Files

| File | Purpose | Tool |
|------|---------|------|
| `replica_*.xtc` | Individual replica trajectories | `remd_trajectory_extractor` |
| `*_stateA.dat`, `*_stateB.dat` | dRMS time series | `remd_drms_analysis` |
| `exchange_ratio.dat` | Exchange efficiency matrix | `remd_exchange_ratio` |
| `replica_state.dat` | State population over time | `remd_replica_state` |
| `FreeEnergy_EXP.pkl` | FES data (EXP mixing) | `REMD_MBAR_Analysis.ipynb` |
| `FreeEnergy_HAM.pkl` | FES data (HAM mixing) | `REMD_MBAR_Analysis.ipynb` |

---

### Key REMD Parameters in remd.inp

| Parameter | Description | Example Value |
|-----------|-------------|---------------|
| `replica_count` | Number of replicas | 11 |
| `replica_c1` | C1 parameter for each replica (kJ/mol) | -480 -440 ... -80 |
| `replica_c2` | C2 parameter | 0 |
| `replica_temp` | Temperature(s) in K | 310.0 |
| `replica_coupling` | Beta = 1/300 | 1/300 |
| `exc_freq` | Exchange attempt frequency (steps) | 250 |
| `remd_unsampled_topfiles` | Single-state topologies for analysis | system_stateA.top system_stateB.top |






