from __future__ import annotations

import openmm as mm
import openmm.unit as unit
import math


class Interaction:
    """Base class for molecular interactions.

    This class serves as a foundation for all interaction types in the
    molecular dynamics simulation framework. It provides common attributes
    and methods for managing interaction forces.

    Attributes:
        name: Name of the interaction.
        description: Description of the interaction type.
        mm_force: The OpenMM force object associated with this interaction.
        contents: List to store interaction parameters.
    """

    def __init__(
        self,
        name: str,
        description: str,
        mm_force: mm.Force | None,
    ) -> None:
        """Initialize the Interaction base class.

        Args:
            name: Name identifier for the interaction.
            description: Human-readable description of the interaction.
            mm_force: OpenMM force object or None if not yet created.
        """
        self.name: str = name
        self.description: str = description
        self.mm_force: mm.Force | None = mm_force

        # default parameters
        self.contents: list = []

    def __str__(self) -> str:
        """Return a string representation of the interaction.

        Returns:
            A formatted string with the interaction name and description.
        """
        return f"{self.name}: {self.description}"


class Nonbonded_interaction(Interaction):
    """Non-bonded interaction combining Lennard-Jones and electrostatic potentials.

    This class implements a shifted Lennard-Jones potential combined with
    a reaction-field modified electrostatic interaction. The potential
    includes both short-range repulsion/dispersion (LJ) and long-range
    electrostatic interactions with proper cutoffs and corrections.

    The force expression includes:
    - Lennard-Jones potential: C12/r^12 - C6/r^6
    - LJ shift correction at cutoff
    - Reaction-field electrostatics with proper boundary conditions

    Attributes:
        mm_force: CustomNonbondedForce with the combined potential.
    """

    def __init__(
        self,
        epsilon_r: float = 15,
        nonbonded_cutoff: unit.Quantity = 1.1 * unit.nanometers,
    ) -> None:
        """Initialize the non-bonded interaction.

        Args:
            epsilon_r: Relative permittivity (dielectric constant).
            nonbonded_cutoff: Cutoff distance for non-bonded interactions.
        """
        super().__init__(
            name="nonbonded_interaction",
            description=(
                "Shifted Lenard-Jones potential and "
                "Reaction-field modified electrostatic interaction"
            ),
            mm_force=None,
        )
        self.mm_force = mm.CustomNonbondedForce(
            "step(rcut-r)*(LJ - corr + ES);"
            # "(LJ - corr + ES);"
            "LJ = (C12(type1, type2) / r^12 - C6(type1, type2) / r^6);"
            "corr = (C12(type1, type2) / rcut^12 - C6(type1, type2) / rcut^6);"
            "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
            "crf = 1 / rcut + krf * rcut^2;"
            "krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {epsilon_r};"
            "f = 138.935458;"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.mm_force.addPerParticleParameter("type")
        self.mm_force.addPerParticleParameter("q")
        self.mm_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        self.mm_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))


class ES_self_excl_interaction(Interaction):
    """Electrostatic self and excluded interactions correction.

    Custom bonded force to add electrostatic terms for self and excluded
    atom pairs. This handles the correction terms needed when atoms are
    excluded from the standard non-bonded interaction calculation.

    Attributes:
        mm_force: CustomBondForce for the ES correction.
    """

    def __init__(
        self,
        epsilon_r: float = 15,
        nonbonded_cutoff: unit.Quantity = 1.1 * unit.nanometers,
    ) -> None:
        """Initialize the ES self/excluded interaction.

        Args:
            epsilon_r: Relative permittivity (dielectric constant).
            nonbonded_cutoff: Cutoff distance for the reaction field.
        """
        super().__init__(
            name="es_self_excl_interaction",
            description=(
                "custom non-bonded force to add in the electrostatic terms "
                "for self and excluded interactions"
            ),
            mm_force=None,
        )

        self.mm_force = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*qprod/epsilon_r * (krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {epsilon_r};"
            f"f = 138.935458;"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.mm_force.addPerBondParameter("qprod")


class ES_except_interaction(Interaction):
    """Electrostatic exception interaction handler.

    Custom bonded force for electrostatic exceptions. This handles
    specific atom pairs that need different electrostatic treatment
    than the standard non-bonded interaction.

    Attributes:
        mm_force: CustomBondForce for ES exceptions.
    """

    def __init__(
        self,
        epsilon_r: float = 15,
        nonbonded_cutoff: unit.Quantity = 1.1 * unit.nanometers,
    ) -> None:
        """Initialize the ES exception interaction.

        Args:
            epsilon_r: Relative permittivity (dielectric constant).
            nonbonded_cutoff: Cutoff distance for the reaction field.
        """
        super().__init__(
            name="es_except_interaction",
            description="custom non-bonded force for the ES exceptions",
            mm_force=None,
        )

        self.mm_force = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*qprod/epsilon_r * (1/r + krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {epsilon_r};"
            f"f = 138.935458;"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.mm_force.addPerBondParameter("qprod")


class LJ_except_interaction(Interaction):
    """Lennard-Jones exception interaction handler.

    Custom bonded force to handle Lennard-Jones exceptions. This allows
    specific atom pairs to have custom LJ parameters that override
    the standard combination rules.

    Attributes:
        mm_force: CustomBondForce for LJ exceptions.
    """

    def __init__(
        self,
        epsilon_r: float = 15,
        nonbonded_cutoff: unit.Quantity = 1.1 * unit.nanometers,
    ) -> None:
        """Initialize the LJ exception interaction.

        Args:
            epsilon_r: Relative permittivity (unused, kept for API consistency).
            nonbonded_cutoff: Cutoff distance for the LJ potential.
        """
        super().__init__(
            name="lj_except_interaction",
            description="custom bonded force to handle exceptions",
            mm_force=None,
        )

        self.mm_force = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C6")


