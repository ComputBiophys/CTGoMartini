"""Multiple basin mixing interaction classes."""

from __future__ import annotations

import openmm as mm

from .base import InteractionError


class MixingError(InteractionError):
    """Raised when multiple basin mixing fails."""
    
    def __init__(self, method: str, n_states: int, reason: str = "") -> None:
        msg = f"Mixing error for method '{method}' with {n_states} states"
        if reason:
            msg = f"{msg}: {reason}"
        super().__init__(msg)


class EXPInteraction:
    """Exponential mixing scheme for multiple basin potential.
    
    Implements the exponential mixing formula for combining
    multiple potential energy basins.
    
    The energy is computed as:
        E = -1/β * log(sum_i exp(-β * (E_i + C_i)))
    
    where β is the coupling constant, E_i is the energy of basin i,
    and C_i is the energy offset for basin i.
    """

    def __init__(self) -> None:
        """Initialize exponential mixing interaction."""
        self.name: str = "exponential mixing scheme"
        self.description: str = 'Exponential mixing scheme for multiple basin potential'

    def addForce(
        self,
        mbp_force_dict: dict[int, list[mm.Force]],
        coupling_constant: float,
        basin_energy_list: list[float],
    ) -> mm.CustomCVForce:
        """Add exponential mixing force.
        
        Args:
            mbp_force_dict: Dictionary mapping state to list of forces.
            coupling_constant: Beta parameter for exponential mixing.
            basin_energy_list: List of basin energy offsets.
            
        Returns:
            The CustomCVForce implementing exponential mixing.
            
        Raises:
            MixingError: If force dictionary is empty or mismatched with basin energies.
        """
        if not mbp_force_dict:
            raise MixingError("EXP", 0, "Empty force dictionary")
        
        if len(mbp_force_dict) != len(basin_energy_list):
            raise MixingError(
                "EXP", len(mbp_force_dict),
                f"Force count ({len(mbp_force_dict)}) doesn't match energy count ({len(basin_energy_list)})"
            )
        
        beta = coupling_constant

        part1_list = []
        part2_list = []
        for state, force_set in mbp_force_dict.items():
            part1_list.append(f"exp(-beta * (energy{state} + C{state}))")
            energy_combined = ' + '.join([f'state{state}_force{j+1}' for j in range(len(force_set))])
            part2_list.append(f"energy{state} = {energy_combined};")

        part1 = ' + '.join(part1_list)
        part2 = '\n'.join(part2_list)
        energy = f"""-1/beta * log({part1});\n{part2}"""
        print(energy)

        self.mm_force = mm.CustomCVForce(energy)
        for state, force_set in mbp_force_dict.items():
            for j, force in enumerate(force_set):
                self.mm_force.addCollectiveVariable(f"state{state}_force{j+1}", force)
        self.mm_force.addGlobalParameter("beta", float(beta))
        for i, basin_energy in enumerate(basin_energy_list):
            self.mm_force.addGlobalParameter(f'C{i+1}', float(basin_energy))
        return self.mm_force


class HAMInteraction:
    """Hamiltonian mixing scheme for multiple basin potential.
    
    Implements the Hamiltonian mixing formula for combining
    two potential energy basins.
    
    The energy is computed as:
        E = (E1 + E2 + ΔV)/2 - sqrt(((E1 - E2 - ΔV)/2)^2 + Δ^2)
    
    where ΔV = C2 - C1 is the energy offset difference and Δ is the coupling constant.
    This method is only valid for exactly 2 states.
    """

    def __init__(self) -> None:
        """Initialize Hamiltonian mixing interaction."""
        self.name: str = "Hamiltonian mixing scheme"
        self.description: str = 'Hamiltonian mixing scheme for multiple basin potential'

    def addForce(
        self,
        mbp_force_dict: dict[int, list[mm.Force]],
        coupling_constant: float,
        basin_energy_list: list[float],
    ) -> mm.CustomCVForce:
        """Add Hamiltonian mixing force.
        
        Args:
            mbp_force_dict: Dictionary mapping state to list of forces.
            coupling_constant: Delta parameter for Hamiltonian mixing.
            basin_energy_list: List of basin energy offsets.
            
        Returns:
            The CustomCVForce implementing Hamiltonian mixing.
            
        Raises:
            MixingError: If the number of states is not exactly 2.
        """
        n_states = len(mbp_force_dict)
        if n_states != 2:
            raise MixingError(
                "HAM", n_states,
                f"HAM method requires exactly 2 states, got {n_states}"
            )
        
        if len(basin_energy_list) != 2:
            raise MixingError(
                "HAM", n_states,
                f"HAM method requires exactly 2 basin energies, got {len(basin_energy_list)}"
            )
        
        delta = coupling_constant

        part2_list = []
        for state, force_set in mbp_force_dict.items():
            energy_combined = ' + '.join([f'state{state}_force{j+1}' for j in range(len(force_set))])
            part2_list.append(f"energy{state} = {energy_combined};")

        part1 = '(energy1+energy2+deltaV)/2 - sqrt(((energy1-energy2-deltaV)/2)^2+delta^2);deltaV=mbp_energy2-mbp_energy1;'
        part2 = '\n'.join(part2_list)
        energy = f"""{part1};\n{part2}"""
        print(energy)

        self.mm_force = mm.CustomCVForce(energy)
        for state, force_set in mbp_force_dict.items():
            for j, force in enumerate(force_set):
                self.mm_force.addCollectiveVariable(f"state{state}_force{j+1}", force)
        self.mm_force.addGlobalParameter("delta", float(delta))
        for i, basin_energy in enumerate(basin_energy_list):
            self.mm_force.addGlobalParameter(f'mbp_energy{i+1}', float(basin_energy))
        return self.mm_force


__all__ = [
    "MixingError",
    "EXPInteraction",
    "HAMInteraction",
]
