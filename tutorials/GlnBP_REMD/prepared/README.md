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

```
#include "FF/martini_v3.0.0.itp"
#include "GlnBP_params.itp"
#include "GlnBP.itp"
#include "FF/martini_v3.0.0_solvents_v1.itp"
#include "FF/martini_v3.0.0_ions_v1.itp"

[ system ]
GlnBP in solvent

[ molecules ]
GlnBP     1
W      6109
NA        67
CL        64
```

**B. system_stateA.top (State A - Open) and system_stateB.top (State B - Closed)**

```bash
cp system.top system_stateA.top
sed -i 's/GlnBP.itp/GlnBP_stateA.itp/g' system_stateA.top

cp system.top system_stateB.top
sed -i 's/GlnBP.itp/GlnBP_stateB.itp/g' system_stateB.top
```

---

### (5) Prepare Input Files

**A. NPT equilibration (npt.inp)**

```ini
; NPT equilibration for GlnBP REMD

mini_nstep  = 10000                             ; Number of steps for minimization
mini_Tol    = 1000.0                            ; Minimization energy tolerance

gen_vel     = yes                               ; Generate initial velocities
gen_temp    = 310                               ; Temperature for generating initial velocities (K)

nstep       = 100000                            ; 100000 * 0.02 ps = 2 ns equilibration
dt          = 0.020                             ; Time step (ps)
b_step      = 0                                 ; Beginning step count
append      = no                                ; Append new dcd to the old one

input       = ions.gro
topol       = system.top
ichk        =                                   ; Load check file
nstout      = 1000                              ; Writing output frequency (steps)
nstdcd      = 1000                              ; Writing coordinates trajectory frequency (steps)
output      = npt.gro
output_pdb  = npt.pdb
odcd        = 
oxtc        = npt.xtc
ochk        = npt.chk

defines     = 

; Restraints - auto-generate
rest        = yes                               ; Turn on/off restraints
rest_ref    = ions.gro                          ; Reference structure file
rest_file   = restraints.txt
gen_rest    = yes                               ; Generate restraint file
atomname    = BB                                ; Select atom name
fc          = 1000.0                            ; Positional restraint force constant
gen_rest_file = restraints.txt                  ; Generated restraint file name

plumed      = no                                ; Turn on/off plumed
plumed_file =                                   ; Default: plumed.dat

platform    = CUDA                              ; CPU, CUDA, OpenCL, Reference
precision   = single                            ; single, mixed, double for CUDA and OpenCL
temp        = 310                               ; Temperature (K)
fric_coeff  = 0.1                               ; Friction coefficient for Langevin dynamics
nonbonded_cutoff = 1.1                          ; Cutoff of the L-J and Coulombic interactions (nm)
epsilon_r   = 15.0                              ; epsilon_r
const_tol   =                                   ; Set the distance tolerance for constraints

pcouple     = yes                               ; Turn on/off pressure coupling
p_ref       = 1.0                               ; Pressure (bar)
p_type      = isotropic                         ; MonteCarloBarostat type: isotropic, membrane
p_freq      = 100                               ; Pressure coupling frequency (steps)
```

Run equilibration:
```bash
run_ctgomartini -i npt.inp
```

This generates `npt.gro` as the starting structure for REMD.

**B. REMD production (remd.inp)**

```ini
; REMD production for GlnBP

nstep       = 50000000      ; 50M * 0.02 ps = 1 µs
dt          = 0.020
temp        = 310
platform    = CUDA
precision   = single

input       = npt.gro
topol       = system.top
output      = npt.gro
output_pdb  = npt.pdb
oxtc        = npt.xtc

; REMD specific parameters
remd        = yes
exc_freq    = 250           ; Exchange every 250 steps (5 ps)

; Replica configuration (11 replicas with C1 gradient)
replica_count = 11
replica_molname = GlnBP
replica_method = exp
replica_temp = 310.0
replica_coupling = 1/300
replica_c1 = -480 -440 -400 -360 -320 -280 -240 -200 -160 -120 -80
replica_c2 = 0

; Unsampled states for MBAR analysis
remd_unsampled_topfiles = system_stateA.top system_stateB.top

; Output
remd_output = output.nc
remd_checkpoint_interval = 5
nstlog      = 1000
nstdcd      = 5000
```

**Key REMD parameters explained:**
- `replica_count`: Number of replicas (11 in this example)
- `replica_c1`: C1 parameter for each replica (kJ/mol), ranging from -480 to -80
- `replica_c2`: C2 parameter (set to 0 for all replicas)
- `replica_temp`: Temperature (constant 310K for all replicas - Hamiltonian REMD)
- `replica_coupling`: Beta = 1/(kB*T), calculated as 1/300
- `remd_unsampled_topfiles`: Single-state topologies for computing state populations

---

### (6) Run REMD Simulation

```bash
# Continue from equilibration
run_ctgomartini -i remd.inp
```

Output files:
- `output.nc` - NetCDF file containing all replica trajectories and exchange information
- `remd.log` - Log file with exchange statistics
- `output.chk` - Checkpoint file for restarting

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

---

### Summary of Files

| File | Purpose |
|------|---------|
| `1GGG_clean.pdb`, `1WDN_clean.pdb` | Input PDB structures |
| `Open/`, `Closed/` | Martinize2 output directories |
| `GlnBP.itp` | Multiple-basin topology (main) |
| `GlnBP_params.itp` | Force field parameters |
| `GlnBP_stateA.itp`, `GlnBP_stateB.itp` | Single-state topologies |
| `ions.gro` | Solvated system structure |
| `system.top`, `system_stateA.top`, `system_stateB.top` | System topologies |
| `npt.inp` | NPT equilibration input |
| `remd.inp` | REMD production input |
| `npt.gro`, `npt.pdb` | Equilibrated structure |
| `output.nc` | REMD trajectory output |

---

### Tips and Troubleshooting

1. **Ion count mismatch**: If you get "wrong number of positions" error, check that `ions.gro` and `system.top` have matching atom counts (6972 atoms).

2. **Missing restraints.txt**: Ensure `gen_rest = yes` is set in `npt.inp` to auto-generate the restraints file.

3. **CUDA out of memory**: Reduce `replica_count` or use fewer replicas per GPU.

4. **Exchange ratio**: A good exchange ratio is 20-40%. If too low, reduce the C1 gradient; if too high, increase it.

5. **Restarting simulations**: Use the checkpoint file to restart:
   ```ini
   ichk = output.chk
   append = yes
   ```
