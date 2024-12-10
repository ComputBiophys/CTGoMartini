# CTGoMartini: A Python package for protein conformational transitions with Go-Martini models

CTGoMartini is a Python package for single-basin Go-Martini, Switching Go-Martini, and Multiple-basin Go-Martini simulations. The main aim of this package is to provide a user-friendly way to simulate conformational transitions of protiens, with particular focus on membrane proteins, using Go-Martini models.


## Installation
### Prerequisites
- conda
- Python 3.8
- dssp

```
conda creaate -n ctgomartini python=3.8
conda activate ctgomartini
pip install -e .
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
- [ ] Vermouth Martinize Functions Update: Update the Vermouth Martinize functions to ensure compatibility with the latest standards and to improve overall functionality.
- [x] Testing Completion: Develop a comprehensive test to validate the integrity and reliability of the project, ensuring that all components function as intended.
- [x] Charge and atom number checkment: Add charge check and atom number check
- [ ] XTCReporter: Use the XTCReporter to replace the DCDReporter, which can reduce the size of output.
- [x] Add Checkpoint at the ending of the simulations
- [ ] Position restraints: x,y,z restraint instead of the distance restraint
- [ ] Analysis module should be added.
- [ ] Minimization output module should be added.
- [ ] Plumed module should be added.
- [x] Add the single-basin contact topology generation