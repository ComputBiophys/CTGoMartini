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
- Basic understanding of Multiple-basin Go-Martini (see GlnBP standard tutorial)

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

Generated files:
- `GlnBP.itp` - Multiple-basin topology (21 replicas)
- `GlnBP_params.itp` - Parameters
- `GlnBP_stateA.itp` - State A (Open) single-state topology
- `GlnBP_stateB.itp` - State B (Closed) single-state topology

Modify the mixing parameters:

```bash
vi GlnBP.itp  # Change beta to 1/300, C1 to appropriate values, C2 to 0
```

---

### (3) Solvate the System

Insert the CG protein into a cubic box and solvate:

```bash
python2 insane.py -f Open/Open_cg.pdb -box 9,9,9 -o ions.gro -salt 0.15 -charge auto -center -sol W 2>system.top
```

Fix ion names and counts:

```bash
sed -i "s/NA+    NA+/NA      NA/g" ions.gro
sed -i "s/CL-    CL-/CL      CL/g" ions.gro
sed -i "s/NA+/NA /g" system.top
sed -i "s/CL-/CL /g" system.top
```

Edit `system.top` to remove excess NA ions (3 ions) that were added due to virtual atoms being counted as charged residues.

---

### (4) Prepare System Topology Files

For REMD, you need three system topology files:

**A. Main system.top (multiple-basin)** - Already created by insane.py, edit to add includes:

```
#include "FF/martini_v3.0.0.itp"
#include "GlnBP_params.itp"
#include "GlnBP.itp"
#include "FF/martini_v3.0.0_solvents_v1.itp"
#include "FF/martini_v3.0.0_ions_v1.itp"
```

**B. system_stateA.top** (State A - Open):

```
#include "FF/martini_v3.0.0.itp"
#include "GlnBP_params.itp"
#include "GlnBP_stateA.itp"
#include "FF/martini_v3.0.0_solvents_v1.itp"
#include "FF/martini_v3.0.0_ions_v1.itp"

[ molecules ]
GlnBP_stateA       1
W             6109
NA              67
CL              64
```

**C. system_stateB.top** (State B - Closed):

```
#include "FF/martini_v3.0.0.itp"
#include "GlnBP_params.itp"
#include "GlnBP_stateB.itp"
#include "FF/martini_v3.0.0_solvents_v1.itp"
#include "FF/martini_v3.0.0_ions_v1.itp"

[ molecules ]
GlnBP_stateB       1
W             6109
NA              67
CL              64
```

---

### (5) Prepare Input Files

**A. NPT equilibration (npt.inp)**

Use the provided `npt.inp` for initial equilibration with position restraints:

```bash
run_ctgomartini -i npt.inp
```

This generates `npt.gro` as the starting structure for REMD.

**B. REMD production (remd.inp)**

Key REMD parameters in `remd.inp`:

```ini
remd        = yes
exc_freq    = 250               ; Exchange every 250 steps
replica_count = 21              ; 21 replicas
replica_c1  = -2000 ... 2000    ; C1 gradient across replicas
remd_unsampled_topfiles = system_stateA.top system_stateB.top
```

The `remd_unsampled_topfiles` parameter specifies the single-state topologies used for computing state populations during analysis.

---

### (6) Run REMD Simulation

```bash
# Continue from equilibration
run_ctgomartini -i remd.inp
```

Output files:
- `output.nc` - NetCDF file containing all replica trajectories and exchange information
- `remd.log` - Log file with exchange statistics

**Note:** REMD simulations are computationally intensive. The example uses 21 replicas, effectively running 21 parallel simulations with exchanges.

---

### (7) Analysis

REMD analysis differs from standard MD. Key analyses include:

**A. Extract replica trajectories:**

```bash
python -m ctgomartini.analysis.remd_replica_state -f output.nc
```

**B. Compute state populations using unsampled states:**

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

---

### Key Concepts

**1. Hamiltonian Replica Exchange**

Unlike temperature REMD, we use Hamiltonian REMD where replicas differ in the energy surface (controlled by C1 parameter). Each replica samples a different bias toward open/closed states.

**2. Unsampled States**

The single-state topologies (`GlnBP_stateA.itp`, `GlnBP_stateB.itp`) are not directly simulated. Instead, they are used as reference states to compute the relative free energy of open vs closed conformations using MBAR analysis.

**3. Exchange Frequency**

`exc_freq = 250` means exchange attempts occur every 250 steps (5 ps with dt=0.02). This should be tuned based on the acceptance ratio (target ~20-30%).

**4. C1 Gradient**

The `replica_c1` values create a linear gradient from -2000 to 2000. Replica 0 strongly favors state A (Open), replica 20 strongly favors state B (Closed), and middle replicas sample both.

---

### Troubleshooting

**Low exchange acceptance:**
- Increase `exc_freq` (less frequent exchanges)
- Adjust C1 gradient range (make it narrower)

**Poor sampling:**
- Increase number of replicas
- Extend simulation time (increase `nstep`)

**Analysis issues:**
- Ensure `remd_unsampled_topfiles` point to valid topology files
- Check that molecule names match in topologies

---

### References

1. CTGoMartini documentation
2. OpenMM documentation for REMD: http://docs.openmm.org
3. MBAR method: Shirts & Chodera, J. Chem. Phys. 129, 124105 (2008)
