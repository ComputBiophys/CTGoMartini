"""
Authors: Song Yang
Last update: 20230510
"""

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

class _OpenMMReadInputs():

    def __init__(self):
        self.mini_nstep       = 0                                         # Number of steps for minimization
        self.mini_Tol         = 1000.0                                    # Minimization energy tolerance
        self.gen_vel          = 'no'                                      # Generate initial velocities
        self.gen_temp         = 300.0                                     # Temperature for generating initial velocities (K)
        self.gen_seed         = None                                      # Seed for generating initial velocities
        self.nstep            = 0                                         # Number of steps to run
        self.dt               = 0.020                                     # Time-step (ps)
        self.b_step           = 0                                         # Set begininng step count. Default: continual
        self.append           = 'no'                                      # Append the dcd to the older one                              

        self.input            = 'input.gro'
        self.topol            = 'system.top'
        self.ichk             = None
        # self.irst             = None
        self.nstout           = 100                                       # Writing output frequency (steps)
        self.nstdcd           = 0                                         # Writing coordinates trajectory frequency (steps)
        self.output           = 'output.gro'
        self.output_pdb       = 'output.pdb'
        self.odcd             = None
        self.oxtc             = None
        self.ochk             = 'output.chk'
        # self.orst             = 'output.rst'

        self.defines          = {}                                        # Defines: unsupport position restraints

        self.rest             = 'no'                                      # Turn on/off restraints
        self.rest_ref         = 'input.gro'                               # Reference structure file
        self.rest_file        = 'restraints.txt'
        self.gen_rest         = 'no'
        self.atomname         = None                                      # Atom name, such as "BB"
        self.fc               = 1000.0                                    # Positional restraint force constant for restraint generation (kJ/mol/nm^2)
        self.gen_rest_file    = 'restraints.txt'                          # Generated restraint file name
        
        self.plumed           = 'no'                                      # Turn on/off plumed
        self.plumed_file      = 'plumed.dat'                              # Plumed input file name

        self.platform         = 'CUDA'                                    # CPU, CUDA, OpenCL, Reference
        self.precision        = 'single'                                  # single, mixed, double for CUDA and OpenCL. For CPU and Reference platforms, default precision is used.
        self.GPU_id           = None                                      # None is default
        self.temp             = 310.0                                     # Temperature (K)
        self.fric_coeff       = 1                                         # Friction coefficient for Langevin dynamics
        self.nonbonded_cutoff = 1.1                                       # nm
        self.const_tol        = None                                      # Set the distance tolerance within which constraints are maintained. 

        self.pcouple          = 'no'                                      # Turn on/off pressure coupling
        self.p_ref            = 1.0                                       # Pressure (Pref or Pxx, Pyy, Pzz; bar)
        self.p_type           = 'membrane'                                # MonteCarloBarotat type: isotropic or membrane
        self.p_scale          = True, True, True                          # For MonteCarloAnisotropicBarostat
        self.p_XYMode         = MonteCarloMembraneBarostat.XYIsotropic    # For MonteCarloMembraneBarostat
        self.p_ZMode          = MonteCarloMembraneBarostat.ZFree          # For MonteCarloMembraneBarostat
        self.p_tens           = 0.0                                       # Sulface tension for MonteCarloMembraneBarostat (dyne/cm)
        self.p_freq           = 15                                        # Pressure coupling frequency (steps)

        self.epsilon_r        = 15.0

    def read(self, inputFile):
        for line in open(inputFile, 'r'):
            if line.find(';') >= 0: line = line.split(';')[0]
            line = line.strip()
            if len(line) > 0:
                segments = line.split('=')
                input_param = segments[0].upper().strip()
                try:    input_value = segments[1].strip()
                except: input_value = None
                if input_value:
                    if input_param == 'MINI_NSTEP':                     self.mini_nstep       = int(input_value)
                    if input_param == 'MINI_TOL':                       self.mini_Tol         = float(input_value)
                    if input_param == 'GEN_VEL':
                        if input_value.upper() == 'YES':                self.gen_vel          = 'yes'
                        if input_value.upper() == 'NO':                 self.gen_vel          = 'no'
                    if input_param == 'GEN_TEMP':                       self.gen_temp         = float(input_value)
                    if input_param == 'GEN_SEED':                       self.gen_seed         = int(input_value)
                    if input_param == 'NSTEP':                          self.nstep            = int(input_value)
                    if input_param == 'DT':                             self.dt               = float(input_value)
                    if input_param == 'B_STEP':                         self.b_step           = int(input_value)
                    if input_param == 'APPEND':                         self.append           = 'yes' if input_value.upper() == 'YES' else 'no'

                    if input_param == 'INPUT':                          self.input            = str(input_value)
                    if input_param == 'TOPOL':                          self.topol            = str(input_value)
                    if input_param == 'ICHK':                           self.ichk             = str(input_value)
                    # if input_param == 'IRST':                           self.irst             = str(input_value)
                    if input_param == 'NSTOUT':                         self.nstout           = int(input_value)
                    if input_param == 'NSTDCD':                         self.nstdcd           = int(input_value)
                    if input_param == 'OUTPUT':                         self.output           = str(input_value)
                    if input_param == 'OUTPUT_PDB':                     self.output_pdb       = str(input_value)
                    if input_param == 'ODCD':                           self.odcd             = str(input_value)
                    if input_param == 'OXTC':                           self.oxtc             = str(input_value)
                    if input_param == 'OCHK':                           self.ochk             = str(input_value)
                    # if input_param == 'ORST':                           self.orst             = str(input_value)

                    if input_param == 'DEFINES':
                        if len(input_value)==0:
                            self.defines = {}
                        else:
                            self.defines = {item.strip():True for item in input_value.split(',')}

                    if input_param == 'REST':
                        if input_value.upper() == 'YES':                self.rest             = 'yes'
                        if input_value.upper() == 'NO':                 self.rest             = 'no'
                    if input_param == 'REST_REF':                       self.rest_ref            = str(input_value)   
                    if input_param == 'REST_FILE':                      self.rest_file        = str(input_value)

                    if input_param == 'GEN_REST':
                        if input_value.upper() == 'YES':                self.gen_rest         = 'yes'
                        if input_value.upper() == 'NO':                 self.gen_rest         = 'no'  
                    if input_param == 'ATOMNAME':                       self.atomname         = str(input_value)
                    if input_param == 'FC':                             self.fc               = float(input_value)    
                    if input_param == 'GEN_REST_FILE':                  self.gen_rest_file    = str(input_value)
                    
                    if input_param == 'PLUMED':
                        if input_value.upper() == 'YES':                self.plumed          = 'yes'
                        if input_value.upper() == 'NO':                 self.plumed          = 'no'
                    if input_param == 'PLUMED_FILE':                    self.plumed_file      = str(input_value)

                    if input_param == 'PLATFORM':                       self.platform         = str(input_value)
                    if input_param == 'PRECISION':                      self.precision        = str(input_value)
                    if input_param == 'GPU_ID':                         self.GPU_id           = str(input_value)
                    if input_param == 'TEMP':                           self.temp             = float(input_value)
                    if input_param == 'FRIC_COEFF':                     self.fric_coeff       = float(input_value)
                    if input_param == 'NONBONDED_CUTOFF':               self.nonbonded_cutoff = float(input_value)
                    if input_param == 'CONST_TOL':                      self.const_tol        = float(input_value)

                    if input_param == 'PCOUPLE':
                        if input_value.upper() == 'YES':                self.pcouple          = 'yes'
                        if input_value.upper() == 'NO':                 self.pcouple          = 'no'
                    if input_param == 'P_REF':
                        if input_value.find(',') < 0:
                            self.p_ref = float(input_value)
                        else:
                            Pxx = float(input_value.split(',')[0])
                            Pyy = float(input_value.split(',')[1])
                            Pzz = float(input_value.split(',')[2])
                            self.p_ref = Pxx, Pyy, Pzz
                    if input_param == 'P_TYPE':
                        if input_value.upper() == 'ISOTROPIC':          self.p_type           = 'isotropic'
                        if input_value.upper() == 'MEMBRANE':           self.p_type           = 'membrane'
                        if input_value.upper() == 'ANISOTROPIC':        self.p_type           = 'anisotropic'
                    if input_param == 'P_SCALE':
                        scaleX = True
                        scaleY = True
                        scaleZ = True
                        if input_value.upper().find('X') < 0: scaleX = False
                        if input_value.upper().find('Y') < 0: scaleY = False
                        if input_value.upper().find('Z') < 0: scaleZ = False
                        self.p_scale = scaleX, scaleY, scaleZ
                    if input_param == 'P_XYMODE':
                        if input_value.upper() == 'XYISOTROPIC':        self.p_XYMode         = MonteCarloMembraneBarostat.XYIsotropic
                        if input_value.upper() == 'XYANISOTROPIC':      self.p_XYMode         = MonteCarloMembraneBarostat.XYAnisotropic
                    if input_param == 'P_ZMODE':
                        if input_value.upper() == 'ZFREE':              self.p_ZMode          = MonteCarloMembraneBarostat.ZFree
                        if input_value.upper() == 'ZFIXED':             self.p_ZMode          = MonteCarloMembraneBarostat.ZFixed
                        if input_value.upper() == 'CONSTANTVOLUME':     self.p_ZMode          = MonteCarloMembraneBarostat.ConstantVolume
                    if input_param == 'P_TENS':                         self.p_tens           = float(input_value)
                    if input_param == 'P_FREQ':                         self.p_freq           = int(input_value)

                    if input_param == 'EPSILON_R':                      self.epsilon_r        = float(input_value)

        return self

def read_inputs(inputFile):
    return _OpenMMReadInputs().read(inputFile)

