"""
MartiniTopology module for CTGoMartini

This module provides the MartiniTopFile class for loading and processing Martini topology files.
It supports creating OpenMM systems from Martini topologies, including handling of
Multiple Basin potentials, virtual sites, and various interaction types.
"""

from __future__ import annotations

import math
from typing import Any, OrderedDict as TypingOrderedDict
from collections import OrderedDict, defaultdict

import numpy as np
import openmm as mm
import openmm.unit as unit
from openmm.app import Topology as mm_Topology

from ..core import *


class MartiniTopFile(Topology):
    """
    A class for loading and processing Martini topology files with extended functionality.
    
    This class extends the Topology class to provide functionality for creating
    OpenMM systems from Martini topology files, including support for Multiple
    Basin potentials, virtual sites, and various interaction types.
    """
    
    def __init__(
        self,
        file: str,
        periodicBoxVectors: tuple[unit.Quantity, unit.Quantity, unit.Quantity] | None = None,
        unitCellDimensions: unit.Quantity | None = None,
        includeDir: str | None = None,
        defines: dict[str, Any] | None = None,
    ) -> None:
        """
        Initialize and load a Martini topology file.

        Args:
            file: Path to the topology file to load.
            periodicBoxVectors: Three vectors defining the periodic box.
                For non-rectangular unit cells, use this instead of unitCellDimensions.
            unitCellDimensions: Dimensions of the crystallographic unit cell.
                For non-rectangular cells, use periodicBoxVectors instead.
            includeDir: Directory to search for #include files.
                If not specified, attempts to locate GROMACS installation.
            defines: Preprocessor definitions for parsing the topology file.
                Example: {"FLEXIBLE": None, "PositionRestraint": "1000"}
        """
        super().__init__(file=file, includeDir=includeDir, defines=defines)

        top = mm_Topology()
        self.topology: mm_Topology = top

        if periodicBoxVectors is not None:
            if unitCellDimensions is not None:
                raise ValueError(
                    "Specify either periodicBoxVectors or unitCellDimensions, but not both"
                )
            top.setPeriodicBoxVectors(periodicBoxVectors)
        else:
            top.setUnitCellDimensions(unitCellDimensions)

        # Add atoms, residues, and bonds
        for molecule_name, molecule_count in self.molecules:
            if molecule_name not in self.moleculeTypes:
                raise ValueError(f"Unknown molecule type: {molecule_name}")
            molecule_type = self.moleculeTypes[molecule_name]

            # Create the specified number of molecules of this type
            for _ in range(molecule_count):
                atoms: list[mm.topology.Atom] = []
                last_residue: str | None = None
                chain = top.addChain()
                for _, fields in enumerate(molecule_type.atoms):
                    res_number = fields[2]
                    if res_number != last_residue:
                        last_residue = res_number
                        res_name = fields[3]
                        residue = top.addResidue(res_name, chain)
                    atom_name = fields[4]
                    atoms.append(top.addAtom(atom_name, None, residue))

                # Add bonds to the topology
                if hasattr(molecule_type, 'bonds'):
                    for bond_fields in molecule_type.bonds:
                        atom1 = atoms[int(bond_fields[0]) - 1]
                        atom2 = atoms[int(bond_fields[1]) - 1]
                        top.addBond(atom1, atom2)
    
    def create_system(
        self,
        nonbonded_cutoff: unit.Quantity = 1.1 * unit.nanometer,
        epsilon_r: float = 15.0,
        remove_com_motion: bool = True,
    ) -> mm.System:
        """
        Construct an OpenMM System from this topology.

        Creates a complete OpenMM System including non-bonded and bonded
        interactions, exception handling, Multiple Basin Potentials, 
        virtual sites, and system initialization.

        Args:
            nonbonded_cutoff: Cutoff distance for non-bonded interactions.
                Default 1.1 nm is typical for Martini simulations.
            epsilon_r: Relative dielectric constant for electrostatics.
                Default 15.0 is standard for Martini in aqueous environments.
            remove_com_motion: If True, add CMMotionRemover to prevent COM drift.

        Returns:
            The configured OpenMM System with all particles and forces.
        """
        # Determine parameter format: True for sigma/epsilon, False for C6/C12
        self.use_sigma_eps: bool = self.forcefield.defaults[0][1] == '2'
        self.nonbonded_cutoff: unit.Quantity = nonbonded_cutoff
        self.epsilon_r: float = epsilon_r
        self.vsites: VSiteManager = VSiteManager()
        self.charges: float | None = None

        # Create OpenMM System and set periodic boundary conditions
        system = mm.System()
        box_vectors = self.topology.getPeriodicBoxVectors()
        if box_vectors is not None:
            system.setDefaultPeriodicBoxVectors(*box_vectors)
        else:
            raise ValueError("periodicBoxVectors must be set")
                  
        all_exceptions: list[tuple[int, int, float, float, float]] = []
        all_charges: list[float] = []

        # Build lookup table mapping atom types to integer indices
        used_atom_types: set[str] = set()
        for molecule_name, _ in self.molecules:
            molecule_type = self.moleculeTypes[molecule_name]
            for atom in molecule_type.atoms:
                used_atom_types.add(atom[1])
        
        atom_type_map: TypingOrderedDict[str, int] = OrderedDict()
        for i, atom_type in enumerate(sorted(used_atom_types)):
            atom_type_map[atom_type] = i
        
        # Calculate Lennard-Jones parameters
        n_types = len(atom_type_map)
        C6, C12 = self._get_LJ_parameters(atom_type_map)

        # Initialize non-bonded interaction forces
        nb_interaction = Nonbonded_interaction(self.epsilon_r, self.nonbonded_cutoff)
        system.addForce(nb_interaction.mm_force)

        # Add tabulated functions for LJ parameters
        nb_interaction.mm_force.addTabulatedFunction(
            "C6", mm.Discrete2DFunction(n_types, n_types, C6)
        )
        nb_interaction.mm_force.addTabulatedFunction(
            "C12", mm.Discrete2DFunction(n_types, n_types, C12)
        )

        # Initialize electrostatic self and exclusion interaction force
        es_self_excl_interaction = ES_self_excl_interaction(self.epsilon_r, self.nonbonded_cutoff)
        system.addForce(es_self_excl_interaction.mm_force)

        # Initialize electrostatic exception interaction force
        es_except_interaction = ES_except_interaction(self.epsilon_r, self.nonbonded_cutoff)
        system.addForce(es_except_interaction.mm_force)

        # Initialize Lennard-Jones exception interaction force
        lj_except_interaction = LJ_except_interaction(self.epsilon_r, self.nonbonded_cutoff)
        system.addForce(lj_except_interaction.mm_force)

        # Initialize non-local bonded interactions
        nonlocal_bonded_force_used: set[Any] = set()
        nonlocal_bonded_interaction_dict = self._Interaction_dict_initialization(
            system, NonLocal_BondedInteraction_dict
        )

        # Process each molecule type
        for molecule_name, molecule_count in self.molecules:
            molecule_type = self.moleculeTypes[molecule_name]
            assert molecule_type.moleculetype[0][1] == "1", (
                f"Only support moleculetype with one exclusion bond length, "
                f"{molecule_type.moleculetype}"
            )

            for _ in range(molecule_count):
                base_atom_index = system.getNumParticles()

                # Add atoms to system
                charges: list[float] = []
                for i, fields in enumerate(molecule_type.atoms):
                    assert len(fields) == 8, f"Too few fields in [ atoms ] lines: {fields}"
                    system.addParticle(mass=float(fields[7]))

                    index = base_atom_index + i
                    q = float(fields[6])
                    charges.append(q)
                    atom_type = atom_type_map[fields[1]]
                    nb_interaction.mm_force.addParticle([atom_type, q])
                    # Add self term for reaction field correction
                    es_self_excl_interaction.mm_force.addBond(index, index, [0.5 * q * q])
                all_charges.extend(charges)
                
                # Process all bonded interactions
                for category in nonlocal_bonded_interaction_dict.keys():
                    if not hasattr(molecule_type, category):
                        continue
                    for fields in getattr(molecule_type, category):
                        fields_used = False
                        for interaction in nonlocal_bonded_interaction_dict[category]:
                            try:
                                interaction.add_interaction(fields, base_atom_index, offset=-1)
                                exceptions = interaction.get_exception(
                                    molecule_type.atoms, fields, base_atom_index, offset=-1
                                )
                                all_exceptions.extend(exceptions)
                                if fields_used:
                                    raise ValueError(f"{fields} is used twice!")
                                fields_used = True
                                if interaction.mm_force is not None:
                                    nonlocal_bonded_force_used.add(interaction.mm_force)
                            except Exception:
                                continue
                        if not fields_used:
                            raise ValueError(f"Cannot recognize the fields: {fields}, {category}")

                # Handle Multiple Basin Potential
                if hasattr(molecule_type, "multiple_basin"):
                    if molecule_type.multiple_basin[0][0].upper() == "TRUE":
                        n_states = int(molecule_type.multiple_basin[0][2])
                        method = molecule_type.multiple_basin[0][1].upper()
                        coupling_constant = eval(molecule_type.multiple_basin[0][3])
                        basin_energy_list = list(map(float, molecule_type.multiple_basin[0][4:]))
                        assert len(basin_energy_list) == n_states, (
                            "The number of energy basins should equal the number of states"
                        )
                        print(
                            f'{molecule_name} uses the multiple basin potential.\n'
                            f'{molecule_type.multiple_basin[0]}'
                        )

                        # Select mixing method
                        if method == "EXP":
                            mbp_interaction = EXP_Interaction()
                        elif method == "HAM":
                            if n_states == 2:
                                mbp_interaction = HAM_Interaction()
                            else:
                                raise ValueError("HAM method not supported for more than 2 states")
                        
                        # MBP Interaction Initialization
                        mbp_force_dict: dict[str, set[Any]] = {str(i+1): set() for i in range(n_states)}
                        mbp_bonded_interaction_dict_list = [
                            self._Interaction_dict_initialization(system, Local_BondedInteraction_dict)
                            for _ in range(n_states)
                        ]

                        # Process interactions for each basin state
                        # Note: mbp_bonded_interaction_dict_list contains 'multi_allbonds' 
                        # which handles all multi-basin interactions (angles, dihedrals, contacts)
                        for category in ['multi_angles', 'multi_dihedrals', 'multi_contacts']:
                            if not hasattr(molecule_type, category):
                                continue
                            for fields in getattr(molecule_type, category):
                                fields_used = False
                                for i in range(n_states):
                                    for interaction in mbp_bonded_interaction_dict_list[i]['multi_allbonds']:
                                        try:
                                            interaction.add_interaction(
                                                str(i+1), category, fields, base_atom_index, offset=-1
                                            )
                                            exceptions = interaction.get_exception(
                                                molecule_type.atoms, category, fields, 
                                                base_atom_index, offset=-1
                                            )
                                            all_exceptions.extend(exceptions)
                                            if fields_used:
                                                raise ValueError(f"{fields} is used twice!")
                                            fields_used = True
                                            mbp_force_dict[str(i+1)].add(interaction.mm_force)
                                        except Exception:
                                            continue
                                if not fields_used:
                                    raise ValueError(f"Cannot recognize the fields: {fields}, {category}")

                        # Create and add MBP force
                        mbp_force = mbp_interaction.addForce(
                            mbp_force_dict, coupling_constant, basin_energy_list
                        )
                        system.addForce(mbp_force)
        
        # Add non-local bonded forces
        for force in nonlocal_bonded_force_used:
            system.addForce(force)
        
        # Add virtual sites
        self.vsites.convert_com_to_linear(system, offset=0)
        for index, site in self.vsites.iter():
            if isinstance(site, OutOfPlane):
                self._add_out_of_plane_vsite(system, index, site, offset=0)
            elif isinstance(site, LinearSite):
                self._add_linear_vsite(system, index, site, offset=0)
            elif isinstance(site, NormalizedInPlaneSite):
                self._add_normalized_in_plane_vsite(system, index, site, offset=0)
            elif isinstance(site, NormalizedInPlaneTwoParticleSite):
                self._add_normalized_in_plane_two_particle_vsite(system, index, site, offset=0)
            else:
                raise RuntimeError(f"Unknown site type {type(site)}.")

        self.all_exceptions = all_exceptions
        
        # Process exceptions
        if all_exceptions:
            except_map: dict[tuple[int, int], list[float]] = defaultdict(list)
            for exception in all_exceptions:
                i, j, q, c6, c12 = exception
                if i < j:
                    except_map[(i, j)] = [q, c6, c12]
                else:
                    except_map[(j, i)] = [q, c6, c12]

            for i, j in except_map:
                nb_interaction.mm_force.addExclusion(i, j)
                
                q, c6, c12 = except_map[(i, j)]
                # Handle electrostatic exceptions and exclusions
                if q == 0:
                    qprod = all_charges[i] * all_charges[j]
                    if qprod != 0.0:
                        es_self_excl_interaction.mm_force.addBond(i, j, [qprod])
                else:
                    es_except_interaction.mm_force.addBond(i, j, [q])
                
                # Handle Lennard-Jones exceptions
                if c6 != 0 and c12 != 0:
                    lj_except_interaction.mm_force.addBond(i, j, [c12, c6])
   
        if remove_com_motion:
            system.addForce(mm.CMMotionRemover())

        self.charges = sum(all_charges)

        return system
    
    def _get_LJ_parameters(
        self, 
        atom_type_map: TypingOrderedDict[str, int]
    ) -> tuple[list[float], list[float]]:
        """
        Calculate Lennard-Jones parameters for all atom type combinations.

        Computes C6 and C12 parameters following Martini force field rules.
        Handles both explicit parameter definitions and combination rules.

        Args:
            atom_type_map: Dictionary mapping atom type names to integer indices.

        Returns:
            Tuple of (C6_list, C12_list), each containing n_types × n_types elements.
        """
        C6: list[float] = []
        C12: list[float] = []
        
        for type_i in atom_type_map:
            for type_j in atom_type_map:
                type_i_sorted, type_j_sorted = sorted([type_i, type_j])
                if (type_i_sorted, type_j_sorted) in self.forcefield.nonbond_params:
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
    
    def _Interaction_dict_initialization(
        self, 
        system: mm.System, 
        Interaction_dict: dict[str, list[Any]]
    ) -> dict[str, list[Any]]:
        """
        Initialize interaction classes with system parameters.

        Args:
            system: The OpenMM System to contain the interactions.
            Interaction_dict: Dictionary mapping categories to interaction classes.

        Returns:
            Dictionary mapping categories to initialized interaction instances.
        """
        interaction_dict: dict[str, list[Any]] = {}
        for category_name, Interaction_list in Interaction_dict.items():
            interaction_dict[category_name] = []
            for Interaction in Interaction_list:
                interaction = Interaction()
                if hasattr(interaction, 'sys'):
                    interaction.sys = system
                if hasattr(interaction, 'vsites'):
                    interaction.vsites = self.vsites
                if hasattr(interaction, 'nonbonded_cutoff'):
                    interaction.nonbonded_cutoff = self.nonbonded_cutoff
                if hasattr(interaction, 'use_sigma_eps'):
                    interaction.use_sigma_eps = self.use_sigma_eps
                    
                interaction_dict[category_name].append(interaction)
        return interaction_dict
    
    def _add_normalized_in_plane_two_particle_vsite(
        self, 
        system: mm.System, 
        index: int, 
        site: NormalizedInPlaneTwoParticleSite, 
        offset: int
    ) -> None:
        """
        Add a normalized in-plane virtual site with two reference particles.
        
        Creates a virtual site in the plane defined by two reference particles.
        
        Args:
            system: The system to add the virtual site to.
            index: Index of the virtual site.
            site: Virtual site object with atom indices and parameters.
            offset: Offset to apply to atom indices.
        """
        vsite = mm.LocalCoordinatesSite(
            site.atom1 + offset,
            site.atom2 + offset,
            site.atom1 + offset,
            [1.0, 0.0, 0.0],
            [-0.5, 0.5, 0.0],
            [0.0, 0.0, 0.0],
            [site.a, 0.0, 0.0],
        )
        system.setVirtualSite(index + offset, vsite)

    def _add_normalized_in_plane_vsite(
        self, 
        system: mm.System, 
        index: int, 
        site: NormalizedInPlaneSite, 
        offset: int
    ) -> None:
        """
        Add a normalized in-plane virtual site with three reference particles.
        
        Args:
            system: The system to add the virtual site to.
            index: Index of the virtual site.
            site: Virtual site object with atom indices and parameters.
            offset: Offset to apply to atom indices.
        """
        vsite = mm.LocalCoordinatesSite(
            site.atom1 + offset,
            site.atom2 + offset,
            site.atom3 + offset,
            [1.0, 0.0, 0.0],
            [-1.0, 1.0 - site.a, site.a],
            [0.0, 0.0, 0.0],
            [site.d, 0.0, 0.0],
        )
        system.setVirtualSite(index + offset, vsite)

    def _add_out_of_plane_vsite(
        self, 
        system: mm.System, 
        index: int, 
        site: OutOfPlane, 
        offset: int
    ) -> None:
        """
        Add an out-of-plane virtual site.
        
        Args:
            system: The system to add the virtual site to.
            index: Index of the virtual site.
            site: Virtual site object with atom indices and parameters.
            offset: Offset to apply to atom indices.
        """
        vsite = mm.OutOfPlaneSite(
            site.atom1 + offset,
            site.atom2 + offset,
            site.atom3 + offset,
            site.a,
            site.b,
            site.c,
        )
        system.setVirtualSite(index + offset, vsite)

    def _add_linear_vsite(
        self, 
        system: mm.System, 
        index: int, 
        site: LinearSite, 
        offset: int
    ) -> None:
        """
        Add a linear virtual site with variable number of reference particles.
        
        Routes to appropriate helper based on number of reference particles.
        
        Args:
            system: The system to add the virtual site to.
            index: Index of the virtual site.
            site: Virtual site object with atom indices and weights.
            offset: Offset to apply to atom indices.
        """
        n = len(site.atom_weights)
        if n == 1:
            self._add_one_particle_vsite(system, index, site, offset)
        elif n == 2:
            self._add_two_particle_vsite(system, index, site, offset)
        elif n == 3:
            self._add_three_particle_vsite(system, index, site, offset)
        else:
            self._add_n_particle_vsite(system, index, site, offset)

    def _add_one_particle_vsite(
        self, 
        system: mm.System, 
        index: int, 
        site: LinearSite, 
        offset: int
    ) -> None:
        """
        Add a one-particle virtual site.
        
        OpenMM requires at least two reference particles, so a dummy
        second particle with zero weight is used.
        
        Args:
            system: The system to add the virtual site to.
            index: Index of the virtual site.
            site: Virtual site object with atom indices and weights.
            offset: Offset to apply to atom indices.
        """
        atoms: list[int] = []
        weights: list[float] = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        assert len(atoms) == 1

        # OpenMM requires at least two atoms
        atoms.append(0)
        weights.append(0.0)

        x_weights = [0.0, 0.0]
        y_weights = [0.0, 0.0]

        vsite = mm.LocalCoordinatesSite(
            atoms, weights, x_weights, y_weights, [0.0, 0.0, 0.0]
        )
        system.setVirtualSite(index + offset, vsite)

    def _add_two_particle_vsite(
        self, 
        system: mm.System, 
        index: int, 
        site: LinearSite, 
        offset: int
    ) -> None:
        """
        Add a two-particle virtual site using weighted average.
        
        Args:
            system: The system to add the virtual site to.
            index: Index of the virtual site.
            site: Virtual site object with atom indices and weights.
            offset: Offset to apply to atom indices.
        """
        atoms: list[int] = []
        weights: list[float] = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        assert len(atoms) == 2

        vsite = mm.TwoParticleAverageSite(atoms[0], atoms[1], weights[0], weights[1])
        system.setVirtualSite(index + offset, vsite)

    def _add_three_particle_vsite(
        self, 
        system: mm.System, 
        index: int, 
        site: LinearSite, 
        offset: int
    ) -> None:
        """
        Add a three-particle virtual site using weighted average.
        
        Args:
            system: The system to add the virtual site to.
            index: Index of the virtual site.
            site: Virtual site object with atom indices and weights.
            offset: Offset to apply to atom indices.
        """
        atoms: list[int] = []
        weights: list[float] = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        assert len(atoms) == 3

        vsite = mm.ThreeParticleAverageSite(
            atoms[0], atoms[1], atoms[2],
            weights[0], weights[1], weights[2],
        )
        system.setVirtualSite(index + offset, vsite)

    def _add_n_particle_vsite(
        self, 
        system: mm.System, 
        index: int, 
        site: LinearSite, 
        offset: int
    ) -> None:
        """
        Add an N-particle virtual site using local coordinates.
        
        Args:
            system: The system to add the virtual site to.
            index: Index of the virtual site.
            site: Virtual site object with atom indices and weights.
            offset: Offset to apply to atom indices.
        """
        atoms: list[int] = []
        weights: list[float] = []
        for atom, weight in site.atom_weights.items():
            atoms.append(atom + offset)
            weights.append(weight)
        n_atoms = len(atoms)

        x_weights = [0.0] * n_atoms
        y_weights = [0.0] * n_atoms
        vsite = mm.LocalCoordinatesSite(
            atoms, weights, x_weights, y_weights, [0.0, 0.0, 0.0]
        )
        system.setVirtualSite(index + offset, vsite)
