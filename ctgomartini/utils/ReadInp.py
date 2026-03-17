"""
Input file reading module for CTGoMartini.

Parses simulation parameter files in GROMACS-style format.
"""

from __future__ import annotations

from openmm.unit import *
from openmm import *
from openmm.app import *


class _OpenMMReadInputs:
    """
    Internal class for storing and parsing simulation input parameters.
    
    Attributes are organized by category:
        - Minimization: mini_nstep, mini_Tol
        - Velocity generation: gen_vel, gen_temp, gen_seed
        - Simulation steps: nstep, dt, b_step, append
        - File I/O: input, topol, ichk, nstout, nstdcd, output, etc.
        - Constraints: defines, const_tol
        - Restraints: rest, rest_ref, rest_file, gen_rest, atomname, fc
        - Enhanced sampling: plumed, plumed_file
        - Platform: platform, precision, GPU_id
        - Thermodynamics: temp, fric_coeff, nonbonded_cutoff
        - Pressure coupling: pcouple, p_ref, p_type, p_scale, p_XYMode, etc.
    """

    def __init__(self) -> None:
        # Minimization parameters
        self.mini_nstep: int = 0
        self.mini_Tol: float = 1000.0
        
        # Velocity generation
        self.gen_vel: str = 'no'
        self.gen_temp: float = 300.0
        self.gen_seed: int | None = None
        
        # Simulation steps
        self.nstep: int = 0
        self.dt: float = 0.020
        self.b_step: int = 0
        self.append: str = 'no'

        # Input/Output files
        self.input: str = 'input.gro'
        self.topol: str = 'system.top'
        self.ichk: str | None = None
        self.nstout: int = 100
        self.nstdcd: int = 0
        self.output: str = 'output.gro'
        self.output_pdb: str = 'output.pdb'
        self.odcd: str | None = None
        self.oxtc: str | None = None
        self.ochk: str = 'output.chk'

        # Defines and constraints
        self.defines: dict[str, bool] = {}
        self.const_tol: float | None = None

        # Restraint parameters
        self.rest: str = 'no'
        self.rest_ref: str = 'input.gro'
        self.rest_file: str = 'restraints.txt'
        self.gen_rest: str = 'no'
        self.atomname: str | None = None
        self.fc: float = 1000.0
        self.gen_rest_file: str = 'restraints.txt'
        
        # Plumed parameters
        self.plumed: str = 'no'
        self.plumed_file: str = 'plumed.dat'

        # Platform parameters
        self.platform: str = 'CUDA'
        self.precision: str = 'single'
        self.GPU_id: str | None = None
        
        # Thermodynamic parameters
        self.temp: float = 310.0
        self.fric_coeff: float = 1.0
        self.nonbonded_cutoff: float = 1.1
        
        # Pressure coupling parameters
        self.pcouple: str = 'no'
        self.p_ref: float | tuple[float, float, float] = 1.0
        self.p_type: str = 'membrane'
        self.p_scale: tuple[bool, bool, bool] = (True, True, True)
        self.p_XYMode = MonteCarloMembraneBarostat.XYIsotropic
        self.p_ZMode = MonteCarloMembraneBarostat.ZFree
        self.p_tens: float = 0.0
        self.p_freq: int = 15

        # Electrostatics
        self.epsilon_r: float = 15.0

        # REMD parameters
        self.remd: str = 'no'           # REMD mode (yes/no)
        self.exc_freq: int | None = None  # Exchange frequency for REMD

    def read(self, input_file: str) -> _OpenMMReadInputs:
        """
        Parse input parameters from a file.
        
        Args:
            input_file: Path to the input parameter file.
            
        Returns:
            Self with parsed parameters set as attributes.
        """
        with open(input_file, 'r') as f:
            for line in f:
                # Remove comments
                if ';' in line:
                    line = line.split(';')[0]
                line = line.strip()
                
                if not line:
                    continue
                    
                segments = line.split('=')
                input_param = segments[0].upper().strip()
                
                try:
                    input_value = segments[1].strip()
                except IndexError:
                    input_value = None
                    
                if not input_value:
                    continue
                    
                self._parse_parameter(input_param, input_value)

        return self

    def _parse_parameter(self, param: str, value: str) -> None:
        """
        Parse a single parameter value.
        
        Args:
            param: Uppercase parameter name.
            value: Parameter value as string.
        """
        # Minimization
        if param == 'MINI_NSTEP':
            self.mini_nstep = int(value)
        elif param == 'MINI_TOL':
            self.mini_Tol = float(value)
            
        # Velocity generation
        elif param == 'GEN_VEL':
            self.gen_vel = 'yes' if value.upper() == 'YES' else 'no'
        elif param == 'GEN_TEMP':
            self.gen_temp = float(value)
        elif param == 'GEN_SEED':
            self.gen_seed = int(value)
            
        # Simulation steps
        elif param == 'NSTEP':
            self.nstep = int(value)
        elif param == 'DT':
            self.dt = float(value)
        elif param == 'B_STEP':
            self.b_step = int(value)
        elif param == 'APPEND':
            self.append = 'yes' if value.upper() == 'YES' else 'no'

        # File I/O
        elif param == 'INPUT':
            self.input = str(value)
        elif param == 'TOPOL':
            self.topol = str(value)
        elif param == 'ICHK':
            self.ichk = str(value)
        elif param == 'NSTOUT':
            self.nstout = int(value)
        elif param == 'NSTDCD':
            self.nstdcd = int(value)
        elif param == 'OUTPUT':
            self.output = str(value)
        elif param == 'OUTPUT_PDB':
            self.output_pdb = str(value)
        elif param == 'ODCD':
            self.odcd = str(value)
        elif param == 'OXTC':
            self.oxtc = str(value)
        elif param == 'OCHK':
            self.ochk = str(value)

        # Defines
        elif param == 'DEFINES':
            if not value:
                self.defines = {}
            else:
                self.defines = {item.strip(): True for item in value.split(',')}

        # Restraints
        elif param == 'REST':
            self.rest = 'yes' if value.upper() == 'YES' else 'no'
        elif param == 'REST_REF':
            self.rest_ref = str(value)
        elif param == 'REST_FILE':
            self.rest_file = str(value)
        elif param == 'GEN_REST':
            self.gen_rest = 'yes' if value.upper() == 'YES' else 'no'
        elif param == 'ATOMNAME':
            self.atomname = str(value)
        elif param == 'FC':
            self.fc = float(value)
        elif param == 'GEN_REST_FILE':
            self.gen_rest_file = str(value)
        
        # Plumed
        elif param == 'PLUMED':
            self.plumed = 'yes' if value.upper() == 'YES' else 'no'
        elif param == 'PLUMED_FILE':
            self.plumed_file = str(value)

        # Platform
        elif param == 'PLATFORM':
            self.platform = str(value)
        elif param == 'PRECISION':
            self.precision = str(value)
        elif param == 'GPU_ID':
            self.GPU_id = str(value)
            
        # Thermodynamics
        elif param == 'TEMP':
            self.temp = float(value)
        elif param == 'FRIC_COEFF':
            self.fric_coeff = float(value)
        elif param == 'NONBONDED_CUTOFF':
            self.nonbonded_cutoff = float(value)
        elif param == 'CONST_TOL':
            self.const_tol = float(value)

        # Pressure coupling
        elif param == 'PCOUPLE':
            self.pcouple = 'yes' if value.upper() == 'YES' else 'no'
        elif param == 'P_REF':
            if ',' not in value:
                self.p_ref = float(value)
            else:
                parts = value.split(',')
                self.p_ref = (float(parts[0]), float(parts[1]), float(parts[2]))
        elif param == 'P_TYPE':
            if value.upper() == 'ISOTROPIC':
                self.p_type = 'isotropic'
            elif value.upper() == 'MEMBRANE':
                self.p_type = 'membrane'
            elif value.upper() == 'ANISOTROPIC':
                self.p_type = 'anisotropic'
        elif param == 'P_SCALE':
            scale_x = 'X' in value.upper()
            scale_y = 'Y' in value.upper()
            scale_z = 'Z' in value.upper()
            self.p_scale = (scale_x, scale_y, scale_z)
        elif param == 'P_XYMODE':
            if value.upper() == 'XYISOTROPIC':
                self.p_XYMode = MonteCarloMembraneBarostat.XYIsotropic
            elif value.upper() == 'XYANISOTROPIC':
                self.p_XYMode = MonteCarloMembraneBarostat.XYAnisotropic
        elif param == 'P_ZMODE':
            if value.upper() == 'ZFREE':
                self.p_ZMode = MonteCarloMembraneBarostat.ZFree
            elif value.upper() == 'ZFIXED':
                self.p_ZMode = MonteCarloMembraneBarostat.ZFixed
            elif value.upper() == 'CONSTANTVOLUME':
                self.p_ZMode = MonteCarloMembraneBarostat.ConstantVolume
        elif param == 'P_TENS':
            self.p_tens = float(value)
        elif param == 'P_FREQ':
            self.p_freq = int(value)

        # Electrostatics
        elif param == 'EPSILON_R':
            self.epsilon_r = float(value)

        # REMD parameters
        elif param == 'REMD':
            self.remd = 'yes' if value.upper() == 'YES' else 'no'
        elif param == 'EXC_FREQ':
            self.exc_freq = int(value)


def read_inputs(input_file: str) -> _OpenMMReadInputs:
    """
    Read simulation parameters from an input file.
    
    Args:
        input_file: Path to the input parameter file.
        
    Returns:
        Parsed input parameters object with all settings as attributes.
        
    Raises:
        ValueError: If REMD is enabled but exc_freq is not specified.
    """
    inputs = _OpenMMReadInputs().read(input_file)
    
    # Validate REMD configuration
    if inputs.remd == 'yes' and inputs.exc_freq is None:
        raise ValueError(
            "REMD is enabled (REMD = yes) but exc_freq is not specified. "
            "Please add 'exc_freq = <number>' to your input file."
        )
    
    return inputs
