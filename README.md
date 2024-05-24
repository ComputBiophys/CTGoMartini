CTGoMartini

# Install
pip install .

# api: the main class
### 1.1 Load topology files 
### 1.2 Generate new topology
### 1.3 Gereate the running file: Switching GoMartini and Multiple-baisn Go-Martini
### 1.4 Provide the analysis scripts
# func: 
### 2.1 Generate the new topology files
### 2.2 Generate the running file
### 2.3 Provide the analysis scripts: Pathway searching and state clustering

# util: 
### 3.1 Functions which are used to load topology files and generate the new topology
### 3.2 Convert to all-atom

# Questions:
CHOL for martini2

Vsite self.nb_force.addExclusion(index + offset, atoms[0]) # ????
Warninging: bonded types in [ atomtypes ]

# Project Roadmap
- [x] Force Class Refinement: Transition the force class to utilize objective classification criteria to enhance clarity and maintainability.
- [x] Multi-basin Force Groups Overhaul: Rebuild the multiple-basin force groups by employing a more efficient approach using two collective variables (CVs) instead of analyzing across all CVs. This will improve performance.
- [ ] Vermouth Martinize Functions Update: Update the Vermouth Martinize functions to ensure compatibility with the latest standards and to improve overall functionality.
- [x] Testing Completion: Develop a comprehensive test to validate the integrity and reliability of the project, ensuring that all components function as intended.
- [x] Charge and atom number checkment: Add charge check and atom number check
- [ ] XTCReporter: Use the XTCReporter to replace the DCDReporter, which can reduce the size of output.
- [x] Add Checkpoint at the ending of the simulations
- [ ] Position restraints: x,y,z restraint instead of the distance restraint
- [ ] Analysis module should be added.
