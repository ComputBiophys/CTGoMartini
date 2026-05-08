## Tutorial: Switching Gō-Martini for Conformational Transitions

In this tutorial, we will utilize GlnBP as an example to illustrate the Switching Gō-Martini method. Unlike Multiple-basin (EXP/HAM) where the protein spontaneously transitions between states during a single simulation, **Switching** runs each state with its own single-basin topology and manually switches between them. This is particularly useful for studying the conformational response to a sudden change in the energy landscape (e.g., ligand binding/dissociation).

### Overview

The switching workflow alternates between three types of MD runs:

| Step | Topology | dt | Purpose |
|------|----------|-----|---------|
| **StateA production** | `system_Open.top` | 0.020 ps | Sample the Open state |
| **Relax** | `system_Closed.top` | **0.002 ps** | Gently adapt to new energy surface |
| **StateB production** | `system_Closed.top` | 0.020 ps | Sample the Closed state |

The **relax step** is critical: when switching topology, old Go contacts are removed and new ones appear instantly. The relax step uses a small timestep from the previous checkpoint to let the system gently release stress, preventing simulation crashes.

---

### (1) Preprocess PDB Files

We use the same GlnBP structures as the Multiple-basin tutorial (PDB codes: 1GGG open, 1WDN closed, residues 5-224). The preprocessed files are provided:

- `1GGG_clean.pdb`: Open-state GlnBP
- `1WDN_clean.pdb`: Closed-state GlnBP

---

### (2) Generate Switching Topologies

Generate separate single-basin topologies for each state:

```bash
ctgomartinize -s 1GGG_clean.pdb 1WDN_clean.pdb \
              -m auto \
              -mol Open Closed \
              -ff martini3001 \
              -dssp \
              -method switching
```

**Output files:**
- `Open/Open_cg.pdb`, `Open/Open.itp`, `Open/Open_params.itp`
- `Closed/Closed_cg.pdb`, `Closed/Closed.itp`, `Closed/Closed_params.itp`
- `Open/` and `Closed/` directories containing martinize2 intermediate files

**Key difference from EXP/HAM:** No combined `GlnBP.itp` is created. Each state has its own independent single-basin Gō-Martini topology. The `-m auto` flag auto-generates contact maps using the OVrCSU algorithm.

---

### (3) Solvate the System

Insert the CG protein into a cubic box and solvate with water and ions using [insane.py](https://cgmartini-library.s3.ca-central-1.amazonaws.com/1_Downloads/tools/proteins_and_bilayers/insane.py). We solvate the Open state as the starting structure:

```bash
python2 insane.py -f Open/Open_cg.pdb -box 9,9,9 -o ions.gro \
    -salt 0.15 -charge auto -center -sol W 2>molecules.txt
```

`insane.py` writes the `[ molecules ]` section to stderr, which we capture to `molecules.txt` for reference.

Fix ion names in the GRO file (insane.py generates incorrect ion names for CG models):

```bash
sed -i "s/NA+    NA+/NA      NA/g" ions.gro
sed -i "s/CL-    CL-/CL      CL/g" ions.gro
```

**Important:** The virtual atoms (Go virtual sites) without charges are recognized as charged residues by insane.py, leading to excess counter-ions. Check `molecules.txt` — for GlnBP, `insane.py` typically reports 70 NA but 3 are excess. Delete them from `ions.gro`:

```bash
# Expected: W 6109, NA 67, CL 64 (total 6972 atoms)
vi ions.gro      # Delete the last 3 NA atoms, update the atom count in the header
```

---

### (4) Prepare System Topology Files

Two system topology files are provided — one for each state. They contain the complete topology (force field includes plus molecule list). After step (3), verify that the molecule counts match your `molecules.txt` (after removing the excess NA):

**system_Open.top** (provided):

```
#include "FF/martini_v3.0.0.itp"
#include "Open_params.itp"
#include "Open.itp"
#include "FF/martini_v3.0.0_solvents_v1.itp"
#include "FF/martini_v3.0.0_ions_v1.itp"

[ system ]
Title of the system

[ molecules ]
Open             1
W             6109
NA              67
CL              64
```

**system_Closed.top** (provided):

```
#include "FF/martini_v3.0.0.itp"
#include "Closed_params.itp"
#include "Closed.itp"
#include "FF/martini_v3.0.0_solvents_v1.itp"
#include "FF/martini_v3.0.0_ions_v1.itp"

[ molecules ]
Closed           1
W             6109
NA              67
CL              64
```

Ensure the molecule counts (W, NA, CL) are **identical** in both files — the only difference is the protein ITP.

---

### (5) Run Equilibration

```bash
run_ctgomartini -i npt.inp
```

This runs minimization + NPT equilibration using the Open-state topology. Generates `npt.gro`, `npt.pdb`, `npt.chk`, `npt.xtc`.

> **Note:** `npt.inp` uses `system_Open.top` as the topology. You only need to equilibrate once — the resulting structure in `npt.chk` will serve as the starting point for the production runs.

---

### (6) Run Switching Simulation

The switching simulation consists of three sequential stages. Each stage loads the checkpoint from the previous stage.

#### Stage 1: StateA Production (Open)

```bash
run_ctgomartini -i md_stateA.inp
```

Runs 500 ns of standard MD with the Open-state Gō-Martini topology. The protein samples conformations around the Open state.

#### Stage 2: Relax (Switch to Closed)

```bash
run_ctgomartini -i md_relax.inp
```

This is the critical switching step:
- Loads coordinates and velocities from `md_stateA.chk`
- Switches to `system_Closed.top` — Go contacts for the Open state are **removed**, Closed-state contacts are **added**
- Uses a small timestep (`dt = 0.002 ps`) for 1000 steps to gently relax the system
- The small dt prevents large forces from the newly-introduced contacts from crashing the simulation

#### Stage 3: StateB Production (Closed)

```bash
run_ctgomartini -i md_stateB.inp
```

Continues from the relaxed checkpoint (`md_relax.chk`) with normal timestep (`dt = 0.020 ps`) for 500 ns. The protein now samples conformations around the Closed state under the Closed Gō-Martini topology.

---

### (7) Analysis

Track the structural response after switching:

```bash
python -m ctgomartini.analysis.drms_analysis \
    -s npt.pdb \
    -f md_stateA.xtc md_relax.xtc md_stateB.xtc \
    -r Open/Open_cg.pdb Closed/Closed_cg.pdb \
    -sel "name BB" \
    -prefix dRMStraj
```

This calculates dRMS to both reference states over time. You should observe:
- During StateA: dRMS to Open is low, dRMS to Closed is high
- At the switch (Relax step): a sudden jump as the energy landscape changes
- During StateB: dRMS to Closed decreases as the protein relaxes into the new basin
