import openmm as mm
import openmm.unit as unit
from simtk.openmm.app import Topology as mm_Topology
from collections import OrderedDict, defaultdict
from ..util import *

class MartiniTopFile(Topology):
    def __init__(
        self,
        file,
        periodicBoxVectors=None,
        unitCellDimensions=None,
        includeDir=None,
        defines=None,
    ):
        """Load a top file.

        Parameters
        ----------
        file : str
            the name of the file to load
        periodicBoxVectors : tuple of Vec3=None
            the vectors defining the periodic box
        unitCellDimensions : Vec3=None
            the dimensions of the crystallographic unit cell.  For
            non-rectangular unit cells, specify periodicBoxVectors instead.
        includeDir : string=None
            A directory in which to look for other files included from the
            top file. If not specified, we will attempt to locate a gromacs
            installation on your system. When gromacs is installed in
            /usr/local, this will resolve to /usr/local/gromacs/share/gromacs/top
        defines : dict={}
            preprocessor definitions that should be predefined when parsing the file
        """
        super().__init__(file=file, 
                         includeDir=includeDir,
                         defines=defines)

        top = mm_Topology()
        self.topology = top

        if periodicBoxVectors is not None:
            if unitCellDimensions is not None:
                raise ValueError(
                    "specify either periodicBoxVectors or unitCellDimensions, but not both"
                )
            top.setPeriodicBoxVectors(periodicBoxVectors)
        else:
            top.setUnitCellDimensions(unitCellDimensions)

        # Add atoms, residues, and bonds
        for molecule_name, molecule_count in self.molecules:
            if molecule_name not in self.moleculeTypes:
                raise ValueError(f"Unknown molecule type: {molecule_name}")
            molecule_type = self.moleculeTypes[molecule_name]

            # Create the specified number of molecules of this type.
            for i in range(molecule_count):
                atoms = []
                lastResidue = None
                c = top.addChain()
                for index, fields in enumerate(molecule_type.atoms):
                    resNumber = fields[2]
                    if resNumber != lastResidue:
                        lastResidue = resNumber
                        resName = fields[3]
                        r = top.addResidue(resName, c)
                    atomName = fields[4]
                    atoms.append(top.addAtom(atomName, None, r))

                # Add bonds to the topology
                if hasattr(molecule_type, 'bonds'):
                    for fields in molecule_type.bonds:
                        top.addBond(atoms[int(fields[0]) - 1], atoms[int(fields[1]) - 1])
            
    def create_system(
            self,
            nonbonded_cutoff=1.1 * unit.nanometer,
            epsilon_r=15.0,
            remove_com_motion=True,
    ):
        """Construct an OpenMM System representing the topology described by this
        top file.

        Parameters
        ----------
        nonbonded_cutoff : distance=1.1 * nanometer
            The cutoff distance to use for nonbonded interactions
        epsilon_r: 15.0
            The espilon for Electrostatic Interction
        remove_com_motion : boolean=True
            If true, a CMMotionRemover will be added to the System

        Returns
        -------
        System
            the newly created System
        """
        self.use_sigma_eps = True if self.forcefield.defaults[0][1] == '2' else False
        self.nonbonded_cutoff = nonbonded_cutoff
        self.epsilon_r = epsilon_r
        self.vsites = VSiteManager()

        sys = mm.System()
        box_vectors = self.topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            sys.setDefaultPeriodicBoxVectors(*box_vectors)
        else:
            raise ValueError("periodicBoxVectors must be set")   
                 
        all_exceptions = []
        all_charges = []

        # build a lookup table mapping atom types into integer indices
        # that will later be used to lookup LJ combinations
        used_atom_types = set()
        for molecule_name, _ in self.molecules:
            molecule_type = self.moleculeTypes[molecule_name]
            for atom in molecule_type.atoms:
                used_atom_types.add(atom[1])
        atom_type_map = OrderedDict()
        for i, k in enumerate(sorted(used_atom_types)):
            atom_type_map[k] = i
        # atom_type_map = {k: i for i, k in enumerate(sorted(used_atom_types))}
        
        # now we need to setup our tables of C6 and C12
        n_types = len(atom_type_map)
        C6, C12 = self._get_LJ_parameters(atom_type_map)

        # Nonbonded interaction initialization
        # custom non-bonded force
        nb_interaction = Nonbonded_interaction(self.epsilon_r, self.nonbonded_cutoff)
        sys.addForce(nb_interaction.mm_force)

        nb_interaction.mm_force.addTabulatedFunction(
            "C6", mm.Discrete2DFunction(n_types, n_types, C6)
        )
        nb_interaction.mm_force.addTabulatedFunction(
            "C12", mm.Discrete2DFunction(n_types, n_types, C12)
        )

        # custom non-bonded force to add in the electrostatic terms for self and excluded interactions
        es_self_excl_interaction = ES_self_excl_interaction(self.epsilon_r, self.nonbonded_cutoff)
        sys.addForce(es_self_excl_interaction.mm_force)   

        # custom non-bonded force for the ES exceptions
        es_except_interaction = ES_except_interaction(self.epsilon_r, self.nonbonded_cutoff)
        sys.addForce(es_except_interaction.mm_force)

        # custom bonded force to handle exceptions
        lj_except_interaction = LJ_except_interaction(self.epsilon_r, self.nonbonded_cutoff)
        sys.addForce(lj_except_interaction.mm_force)

        # Nonlocal Bonded interaction initialization
        nonlocal_bonded_force_used = set()
        nonlocal_bonded_interaction_dict = self._Interaction_dict_initialization(sys, NonLocal_BondedInteraction_dict)

        # Loop over molecules and create the specified number of each type.
        for molecule_name, molecule_count in self.molecules:
            molecule_type = self.moleculeTypes[molecule_name]
            assert molecule_type.moleculetype[0][1] == "1", f"Only support moleculetype with one exclusion bond length, {molecule_type.moleculetype}"

            for _ in range(molecule_count):
                base_atom_index = sys.getNumParticles()

                # add atoms parameters
                charges = []
                for i, fields in enumerate(molecule_type.atoms):
                    assert len(fields) == 8, f"Too few fields in [ atoms ] lines: {fields}"
                    sys.addParticle(mass=float(fields[7]))

                    index = base_atom_index + i
                    q = float(fields[6])
                    charges.append(q)
                    atomType = atom_type_map[fields[1]]
                    nb_interaction.mm_force.addParticle([atomType, q])
                    # Add in the self term for the reaction field correction
                    es_self_excl_interaction.mm_force.addBond(index, index, [0.5 * q * q])
                all_charges.extend(charges)
                
                # add bonds, angles, dihedrals, cmaps, virutal_sitesn, contacts, pairs, exclusions, ...
                for category in nonlocal_bonded_interaction_dict.keys():
                    if not hasattr(molecule_type, category):
                        continue
                    for fields in getattr(molecule_type, category):
                        fields_used = False
                        for interaction in nonlocal_bonded_interaction_dict[category]:
                            try:
                                interaction.add_interaction(fields, base_atom_index, offset=-1)
                                exceptions = interaction.get_exception(molecule_type.atoms, fields, base_atom_index, offset=-1)
                                all_exceptions.extend(exceptions)
                                if fields_used:
                                    raise ValueError(f"{fields} is used twice!")
                                else:
                                    fields_used = True
                                if interaction.mm_force is not None:
                                    nonlocal_bonded_force_used.add(interaction.mm_force)
                            except:
                                continue
                        if not fields_used:
                            raise ValueError(f"Cannot recoginze the fields: {fields}, {category}")

                # add MBP
                if hasattr(molecule_type, "multiple_basin"):
                    if molecule_type.multiple_basin[0][0].upper() == "TRUE":
                        n_states = int(molecule_type.multiple_basin[0][2])
                        method = molecule_type.multiple_basin[0][1].upper()
                        coupling_constant = eval(molecule_type.multiple_basin[0][3])
                        basin_energy_list = list(map(float, molecule_type.multiple_basin[0][4:]))
                        assert len(basin_energy_list) == n_states, "The number of enerngy basin should equal to the number of states"
                        print(f'{molecule_name} uses the multiple basin potential.\n{molecule_type.multiple_basin[0]}')

                        # elect the suitable mixing method 
                        if method == "EXP":
                            mbp_interaction = EXP_Interaction()
                        elif method == "HAM":
                            if n_states == 2:
                                mbp_interaction = HAM_Interaction()
                            else:
                                raise ValueError("Unsupport HAM multiple basin potential for more than 2 states")
                        
                        # MBP Interaction Initialization
                        mbp_force_dict = {str(i+1): set() for i in range(n_states)}
                        mbp_bonded_interaction_dict_list = [self._Interaction_dict_initialization(sys, Local_BondedInteraction_dict) for i in range(n_states)]

                        for category in mbp_bonded_interaction_dict_list[0].keys():
                            if not hasattr(molecule_type, category):
                                continue
                            for fields in getattr(molecule_type, category):
                                fields_used = False
                                for i in range(n_states):
                                    for interaction in mbp_bonded_interaction_dict_list[i]['multi_allbonds']:    # Use the multi_allbonds instead of the individual category
                                        try:
                                            interaction.add_interaction(str(i+1), category, fields, base_atom_index, offset=-1)
                                            exceptions = interaction.get_exception(molecule_type.atoms, category, fields, base_atom_index, offset=-1)
                                            all_exceptions.extend(exceptions)
                                            if fields_used:
                                                raise ValueError(f"{fields} is used twice!")
                                            else:
                                                fields_used = True
                                            mbp_force_dict[str(i+1)].add(interaction.mm_force)
                                        except:
                                            continue
                                if not fields_used:
                                    raise ValueError(f"Cannot recoginze the fields: {fields}, {category}")

                        mbp_force = mbp_interaction.addForce(mbp_force_dict, coupling_constant, basin_energy_list)
                        sys.addForce(mbp_force)
        
        # add nonlocal bonded forces
        for force in nonlocal_bonded_force_used:
            sys.addForce(force)
        
        # add virtual sites to systems
        self.vsites.convert_com_to_linear(sys, offset=0)
        for index, site in self.vsites.iter():
            if isinstance(site, OutOfPlane):
                self._add_out_of_plane_vsite(sys, index, site, offset=0)
            elif isinstance(site, LinearSite):
                self._add_linear_vsite(sys, index, site, offset=0)
            elif isinstance(site, NormalizedInPlaneSite):
                self._add_normalized_in_plane_vsite(sys, index, site, offset=0)
            elif isinstance(site, NormalizedInPlaneTwoParticleSite):
                self._add_normalized_in_plane_two_particle_vsite(sys, index, site, offset=0)
            else:
                raise RuntimeError(f"Unknown site type {type(site)}.")            

        self.all_exceptions = all_exceptions
        # Exceptions
        if all_exceptions:
            # build a map of unique exceptions
            # process in order, so that later entries trump earlier ones
            except_map = defaultdict(list)
            for exception in all_exceptions:
                i, j, q, c6, c12 = exception
                if i < j:
                    except_map[(i, j)] = [q, c6, c12]
                else:
                    except_map[(j, i)] = [q, c6, c12]

            # add in all of the exclusions
            for i, j in except_map:
                # Remove i,j from nonbonded interactions for all exceptions / exclusions
                nb_interaction.mm_force.addExclusion(i, j)
                
                q, c6, c12 = except_map[(i,j)] # ys modified
                # Handle electrostatic exceptions / exclusions.
                # We're going to assume that q==0 means that this was an
                # exclusion.
                if q == 0:
                    # In this case, we still need to add in the reaction field correction
                    # term.
                    qprod = all_charges[i] * all_charges[j]
                    # Don't bother adding the interaction if either particle has zero charge.
                    if qprod != 0.0:
                        es_self_excl_interaction.mm_force.addBond(i, j, [qprod])
                # If q !=0, then this is an exception.
                else:
                    # We don't bother adding the interaction if either the charge product is zero
                    es_except_interaction.mm_force.addBond(i, j, [q])
                # Now we'll add the LJ exceptions. We don't bother
                # adding the interaction if the combined LJ parameters are zero.
                if c6 != 0 and c12 != 0:
                    lj_except_interaction.mm_force.addBond(i, j, [c12, c6])
   
        if remove_com_motion:
            sys.addForce(mm.CMMotionRemover())
        return sys
                             

    def _get_LJ_parameters(self, atom_type_map):
        """
        Get LJ parameters for a particular atom type.

        Parameters
        ----------
        atom_type_map : dict of {str: int}
            Dictionary mapping atom types to atom indices.

        Returns
        -------
        C6 : list
            C6 paramerets based on the order of atom_type_map
        C12: list
            C12 paramerets based on the order of atom_type_map
        """
        C6 = []
        C12 = []
        for type_i in atom_type_map:
            for type_j in atom_type_map:
                type_i_sorted, type_j_sorted = sorted([type_i, type_j])
                if (type_i_sorted, type_j_sorted) in self.forcefield.nonbond_params:
                    # Parameters were specified for this pair of types,
                    # so use them.
                    params = self.forcefield.nonbond_params[(type_i_sorted, type_j_sorted)]
                    if self.use_sigma_eps:
                        sigma = float(params[3])
                        eps = float(params[4])
                        c6 = 4 * eps * sigma ** 6
                        c12 = 4 * eps * sigma ** 12
                    else:
                        c6 = float(params[3])
                        c12 = float(params[4])
                else:
                    # Parameters were not specified for this pair type,
                    # so calculate using combination rules.
                    params_i = self.forcefield.atomtypes[type_i]
                    v_i = float(params_i[6])
                    w_i = float(params_i[7])
                    params_j = self.forcefield.atomtypes[type_j]
                    v_j = float(params_j[6])
                    w_j = float(params_j[7])

                    if self.use_sigma_eps:
                        sigma = 0.5 * (v_i + v_j)
                        eps = math.sqrt(w_i * w_j)
                        c6 = 4 * eps * sigma ** 6
                        c12 = 4 * eps * sigma ** 12
                    else:
                        c6 = math.sqrt(v_i * v_j)
                        c12 = math.sqrt(w_i * w_j)

                C6.append(c6)
                C12.append(c12)
        return C6, C12
    
    def _Interaction_dict_initialization(self, sys, Interaction_dict):
        """
        Initialize the Interaction class and assign the system, self.nonbonded_cutoff, self.use_sigma_eps

        Parameters
        -----
        sys: openmm.system
        Intearction_dict: dict
            NonLocal_BondedInteraction_dict or Local_BondedInteraction_dict
            {"category": [Interaction1_class, Interaction2_class, ...]}


        Return:
        -----
        interaction_dict: dict
            {"category": [Interaction1_class(), Interaction2_class(), ...]}
        """
        interaction_dict = {}
        for catetogy_name, Interaction_list in Interaction_dict.items():
            interaction_dict[catetogy_name] = []
            for Interaction in Interaction_list:
                interaction = Interaction()
                if hasattr(interaction, 'sys'):
                    interaction.sys = sys
                if hasattr(interaction, 'vsites'):
                    interaction.vsites = self.vsites
                if hasattr(interaction, 'nonbonded_cutoff'):
                    interaction.nonbonded_cutoff = self.nonbonded_cutoff
                if hasattr(interaction, 'use_sigma_eps'):
                    interaction.use_sigma_eps = self.use_sigma_eps
                    
                interaction_dict[catetogy_name].append(interaction)
        return interaction_dict
    

    def _add_normalized_in_plane_two_particle_vsite(self, sys, index, site, offset):
        vsite = mm.LocalCoordinatesSite(
            site.atom1 + offset,
            site.atom2 + offset,
            site.atom1 + offset,
            [1.0, 0.0, 0.0],
            [-0.5, 0.5, 0.0],
            [0.0, 0.0, 0.0],
            [site.a, 0.0, 0.0],
        )
        sys.setVirtualSite(index + offset, vsite)

    def _add_normalized_in_plane_vsite(self, sys, index, site, offset):
        vsite = mm.LocalCoordinatesSite(
            site.atom1 + offset,
            site.atom2 + offset,
            site.atom3 + offset,
            [1.0, 0.0, 0.0],
            [-1.0, 1.0 - site.a, site.a],
            [0.0, 0.0, 0.0],
            [site.d, 0.0, 0.0],
        )
        sys.setVirtualSite(index + offset, vsite)

    def _add_out_of_plane_vsite(self, sys, index, site, offset):
        vsite = mm.OutOfPlaneSite(
            site.atom1 + offset,
            site.atom2 + offset,
            site.atom3 + offset,
            site.a,
            site.b,
            site.c,
        )
        sys.setVirtualSite(index + offset, vsite)

    def _add_linear_vsite(self, sys, index, site, offset):
        n = len(site.atom_weights)
        if n == 1:
            self._add_one_particle_vsite(sys, index, site, offset)
        elif n == 2:
            self._add_two_particle_vsite(sys, index, site, offset)
        elif n == 3:
            self._add_three_particle_vsite(sys, index, site, offset)
        else:
            self._add_n_particle_vsite(sys, index, site, offset)

    def _add_one_particle_vsite(self, sys, index, site, offset):
        atoms = []
        weights = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        assert len(atoms) == 1

        # We need at least two atoms, so we'll use atom 0 as
        # dummy with a weight of 0.
        atoms.append(0)
        weights.append(0.0)

        # These are dummy weights that specify a local coordinate system
        # that openmm supports, but we don't use.
        x_weights = [0.0, 0.0]
        y_weights = [0.0, 0.0]

        vsite = mm.LocalCoordinatesSite(
            atoms, weights, x_weights, y_weights, [0.0, 0.0, 0.0]
        )
        sys.setVirtualSite(index + offset, vsite)

        # Add exclusions
        # self.nb_force.addExclusion(index + offset, atoms[0]) ???

    def _add_two_particle_vsite(self, sys, index, site, offset):
        atoms = []
        weights = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        assert len(atoms) == 2

        vsite = mm.TwoParticleAverageSite(atoms[0], atoms[1], weights[0], weights[1])
        sys.setVirtualSite(index + offset, vsite)

    def _add_three_particle_vsite(self, sys, index, site, offset):
        atoms = []
        weights = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        assert len(atoms) == 3

        vsite = mm.ThreeParticleAverageSite(
            atoms[0],
            atoms[1],
            atoms[2],
            weights[0],
            weights[1],
            weights[2],
        )
        sys.setVirtualSite(index + offset, vsite)

    def _add_n_particle_vsite(self, sys, index, site, offset):
        atoms = []
        weights = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        n_atoms = len(atoms)

        # These are dummy weights that specify a local coordinate system
        # that openmm supports, but we don't use.
        x_weights = [0.0] * n_atoms
        y_weights = [0.0] * n_atoms
        vsite = mm.LocalCoordinatesSite(
            atoms, weights, x_weights, y_weights, [0.0, 0.0, 0.0]
        )
        sys.setVirtualSite(index + offset, vsite)