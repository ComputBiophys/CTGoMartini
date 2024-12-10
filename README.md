# CTGoMartini: A Python package for protein conformational transitions with Go-Martini models

CTGoMartini is a Python package for single-basin Go-Martini, Switching Go-Martini, and Multiple-basin Go-Martini simulations. The main aim of this package is to provide a user-friendly scritps for conformational transition simulations of protiens, especially membrane proteins with Go-Martini models.


## Installation
### Prerequisites
conda
Python 3.8
dssp


'conda creaate -n ctgomartini python=3.8
conda activate ctgomartini
pip install -e .'

## Features
- "Contacts" intereaction type is defined in this code to replace the LJ-type contact interactions in traditional GoMartini3 models, which can overcome the mistaken interactions between multiple copies of proteins in traditional GoMartini3 models and provide a more convient way to construct the multiple basin potentials in Multiple-basin GoMartini simulations.
- This code can construct the topology of the GoMartini3 models for different methods automatically.
- This can run the coarse-grained simulations easily.

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