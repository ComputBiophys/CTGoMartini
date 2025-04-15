# CTGoMartini: A Python package for protein conformational transitions with Go-Martini models

CTGoMartini is a Python package for single-basin Go-Martini, Switching Go-Martini, and Multiple-basin Go-Martini simulations. The main aim of this package is to provide a user-friendly way to simulate conformational transitions of protiens, with particular focus on membrane proteins, using Go-Martini models.


## Installation
### Prerequisites
- conda
- Python >= 3.8
- dssp < 4.0
- openmm >=8.1
```
conda create -n ctgomartini python=3.8
conda activate ctgomartini
conda install -c conda-forge openmm
pip install -e .

# if use plumed function, please install openmm-plumed
conda install -c conda-forge openmm-plumed
```

Verify your installation by typing the following command
```
python -m ctgomartini.tests.tests
```

## Features
- Implement a new interaction type named "Contacts" to replace the LJ-type contact interactions in classic GoMartini3 models, which can (1) eliminate incorrect interactions between multiple protein copies in classic GoMartini3 models and (2) facilitate construction of multiple basin potentials in Multiple-basin GoMartini simulations.
- Automatically construct GoMartini3 model topologies for different simulation methods.
- Simplify the process of running coarse-grained simulations.

## Tutorial
- Please see the directory Tutorial





## Project Roadmap
- [x] Force Class Refinement: Transition the force class to utilize objective classification criteria to enhance clarity and maintainability.
- [x] Multi-basin Force Groups Overhaul: Rebuild the multiple-basin force groups by employing a more efficient approach using two collective variables (CVs) instead of analyzing across all CVs. This will improve performance.
- [ ] Vermouth Martinize Functions Update: Update the Vermouth Martinize functions to ensure compatibility with the latest standards and to improve overall functionality. Will be done until Vermouth becomes stable 
- [x] Testing Completion: Develop a comprehensive test to validate the integrity and reliability of the project, ensuring that all components function as intended.
- [x] Charge and atom number checkment: Add charge check and atom number check
- [x] XTCReporter: Use the XTCReporter to replace the DCDReporter, which can reduce the size of output.
- [x] Add Checkpoint at the ending of the simulations
- [x] Add the constraint tolerance parameter in .inp file to better manage the constraint
- [x] Position restraints: x,y,z restraint instead of the distance restraint
- [x] Reconstruct the tests module
- [ ] Analysis module should be added.
- [x] Minimization output module should be added. This function can be achieved simplely by setting the nstep in npt.inp as zero.
- [x] Plumed module should be added.
- [x] Add the single-basin contact topology generation
- [ ] Refine Create_goVirt_for_multimer by fixing the extracting CA/BB bugs
- [ ] Seperate the tests and codes to different respositories in order to make the project more clear.
- [ ] Reconstruct the codes with AI tools to make the codes more readable.
- [ ] More tutorials should be added, especially for the membrane systems.
- [ ] Modify the suitable friction coefficient parameter in the inp file. 