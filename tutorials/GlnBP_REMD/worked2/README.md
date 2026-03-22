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

REMD analysis differs from standard MD. Key analyses include:

#### A. Extract replica trajectories:

```bash
python -m ctgomartini.analysis.remd_replica_state -f output.nc
```

This generates `replica_*.xtc` files for each replica.


#### B. Compute state populations using unsampled states (MBAR):

```bash
python -m ctgomartini.analysis.remd_free_energy -f output.nc --single -o free_energy_profile.pdf
```

For comparing multiple REMD runs:
```bash
python -m ctgomartini.analysis.remd_free_energy \
    -f output1.nc output2.nc output3.nc \
    -l "Run 1" "Run 2" "Run 3" \
    -o free_energy_comparison.pdf
```

Options:
- `-f FILES`: One or more NetCDF files
- `-l LABELS`: Labels for each run
- `-o OUTPUT`: Output plot file
- `-u FILE`: Unbiased sampling data for comparison
- `--single`: Plot single profile instead of comparison

#### C. Plot replica state trajectories:

```bash
python -m ctgomartini.analysis.remd_replica_state -f output.nc -o replica_states.pdf
```

To check state occupancies:
```bash
python -m ctgomartini.analysis.remd_replica_state -f output.nc --occupancies
```

Options:
- `-f FILE`: Input NetCDF file
- `-o OUTPUT`: Output plot file
- `--dt DT`: Time step in microseconds (default: 0.005)
- `--skip N`: Frame skip for plotting (default: 100)

#### D. Exchange ratio analysis:

```bash
python -m ctgomartini.analysis.remd_exchange_ratio -f output.nc
```

Options:
- `-f FILE`: Input NetCDF file
- `-o OUTPUT`: Output file for exchange statistics

#### E. dRMS analysis on replica trajectories:

```bash
python -m ctgomartini.analysis.drms_analysis \
    -s npt.pdb \
    -f replica_0.xtc \
    -r Open/Open_cg.pdb Closed/Closed_cg.pdb \
    -sel "name BB" \
    -prefix dRMStrj_replica0
```

For all replicas:
```bash
for i in {0..10}; do
    python -m ctgomartini.analysis.drms_analysis \
        -s npt.pdb \
        -f replica_${i}.xtc \
        -r Open/Open_cg.pdb Closed/Closed_cg.pdb \
        -sel "name BB" \
        -prefix dRMStrj_replica${i}
done
```

#### F. MBAR parameter optimization analysis (if using remd_mbar_analysis):

```bash
# Analyze selected states convergence
python -m ctgomartini.analysis.remd_mbar_analysis \
    -t selected_states -f selected_states_results.pkl

# Analyze start ratio effects
python -m ctgomartini.analysis.remd_mbar_analysis \
    -t start_ratio -f start_ratio_results.pkl

# Compute equilibrium constant
python -m ctgomartini.analysis.remd_mbar_analysis \
    -t selected_states -f results.pkl --keq --temp 310
```

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






