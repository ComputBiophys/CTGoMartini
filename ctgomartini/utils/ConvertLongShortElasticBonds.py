"""
convert_long_short_elastic_bonds.py

This module converts long and short elastic bonds in Martini coarse-grained molecular dynamics simulations.
It processes topology files to transform elastic bonds between backbone (BB) beads into Go-model interactions
between alpha-carbon (CA) beads, which are more suitable for certain types of simulations.

The main function performs the following steps:
1. Validates input files
2. Maps BB beads to their corresponding CA beads
3. Identifies long and short elastic bonds based on force constants and distances
4. Converts these bonds into exclusions and Go-model interactions
5. Updates the relevant topology files

"""

from __future__ import annotations

import os
import subprocess
from typing import TYPE_CHECKING

import MDAnalysis as mda
import numpy as np

from ctgomartini.api import MartiniTopFile

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe


def bb_distance(atomid1: str, atomid2: str, bb_list: list[str]) -> int:
    """Calculate the sequence distance between two backbone beads in a polymer chain.

    Args:
        atomid1: Atom ID of the first backbone bead.
        atomid2: Atom ID of the second backbone bead.
        bb_list: Ordered list of backbone bead IDs in the polymer chain.

    Returns:
        The absolute sequence distance between the two backbone beads.
    """
    distance = abs(bb_list.index(atomid1) - bb_list.index(atomid2))
    return distance


def convert_long_short_elastic_bonds(
    prefix: str,
    ref_pdb: str,
    convert_long_elastic_bonds: bool = True,
    convert_short_elastic_bonds: bool = False,
    lj_epsilon: float = 12.00,
) -> None:
    """Convert long and short elastic bonds in Martini topology to Go-model interactions.

    This function processes a Martini topology file to identify elastic bonds between
    backbone beads and converts them to Go-model interactions between alpha-carbon beads.
    It updates the topology, exclusions, and Go-pair files accordingly.

    Args:
        prefix: Prefix for the input/output files (e.g., protein name).
        ref_pdb: Path to the reference PDB file for distance calculations.
        convert_long_elastic_bonds: Whether to convert long elastic bonds. Defaults to True.
        convert_short_elastic_bonds: Whether to convert short elastic bonds. Defaults to False.
        lj_epsilon: Lennard-Jones epsilon parameter for Go-model interactions. Defaults to 12.00.

    Raises:
        ValueError: If any of the required input files are not found.
    """
    # Check files
    prot_itp = f"{prefix}.itp"
    go_pair_file = f"{prefix}_go-table_VirtGoSites.itp"
    exclusion_file = f"{prefix}_exclusions_VirtGoSites.itp"
    if not os.path.exists(prot_itp):
        raise ValueError(f"Error: cannot find {prot_itp}!")
    if not os.path.exists(go_pair_file):
        raise ValueError(f"Error: cannot find {go_pair_file}!")
    if not os.path.exists(exclusion_file):
        raise ValueError(f"Error: cannot find {exclusion_file}!")
    if not os.path.exists(ref_pdb):
        raise ValueError(f"Error: cannot find {ref_pdb}!")

    # Load the topology file and extract molecule information
    top = MartiniTopFile(prot_itp)
    mol = top._moleculeTypes[prefix]

    # Create a dictionary mapping backbone (BB) beads to their corresponding alpha-carbon (CA) beads
    # This is necessary because elastic bonds are defined between BB beads but we want to convert
    # them to interactions between CA beads for the Go-model
    bb2ca_dict: dict[str, str] = {}
    atoms: list[list[str]] = mol._topology["atoms"]
    for fields in mol._topology["virtual_sitesn"]:
        # Check if this is a virtual site definition for a CA bead (type 1)
        if len(fields) == 3 and fields[1] == "1":
            # Verify that the atoms are correctly identified as CA and BB
            assert (
                atoms[int(fields[0]) - 1][0] == fields[0]
                and atoms[int(fields[0]) - 1][4] == "CA"
            )
            assert (
                atoms[int(fields[2]) - 1][0] == fields[2]
                and atoms[int(fields[2]) - 1][4] == "BB"
            )

            # Map the BB bead to its corresponding CA bead
            if fields[2] not in bb2ca_dict:
                bb2ca_dict[fields[2]] = fields[0]
            else:
                print("Warning: duplicated BB to CA for ", fields)

    # Create a list of all backbone beads in order
    bb_list: list[str] = []
    for fields in atoms:
        if fields[4] == "BB":
            bb_list.append(fields[0])

    # Identify long and short elastic bonds based on force constants and sequence distances
    # Long elastic bonds: BB distance = 3, force constant (k) = 0.970
    # Short elastic bonds: BB distance = 2, force constant (k) = 0.640
    long_elastic_bonds: list[tuple[str, str]] = []  # [(atomid1, atomid2)]
    short_elastic_bonds: list[tuple[str, str]] = []
    bonds: list[list[str]] = mol._topology["bonds"]
    for fields in bonds:
        # Check for long elastic bonds
        if float(fields[3]) == 0.970 and bb_distance(fields[0], fields[1], bb_list) == 3:
            atomid1, atomid2 = fields[0], fields[1]
            long_elastic_bonds.append((atomid1, atomid2))
        # Check for short elastic bonds
        if float(fields[3]) == 0.640 and bb_distance(fields[0], fields[1], bb_list) == 2:
            atomid1, atomid2 = fields[0], fields[1]
            short_elastic_bonds.append((atomid1, atomid2))

    # Convert the BB bead pairs to CA bead pairs using the mapping dictionary
    ca_long_elastic_bonds: list[tuple[str, str]] = []
    for item in long_elastic_bonds:
        ca_long_elastic_bonds.append((bb2ca_dict[item[0]], bb2ca_dict[item[1]]))
    ca_short_elastic_bonds: list[tuple[str, str]] = []
    for item in short_elastic_bonds:
        ca_short_elastic_bonds.append((bb2ca_dict[item[0]], bb2ca_dict[item[1]]))

    # Generate exclusion entries for the converted bonds
    # Exclusions prevent non-bonded interactions between atoms that should not interact
    line_style = " {}  {}           ;  {}  {}\n"
    exclusions_extend: list[str] = []
    if convert_long_elastic_bonds:
        for item in long_elastic_bonds:
            atomid1, atomid2 = item
            resid1, resid2 = (
                atoms[int(atomid1) - 1][2],
                atoms[int(atomid2) - 1][2],
            )
            line = line_style.format(atomid1, atomid2, resid1, resid2)
            exclusions_extend.append(line)

    if convert_short_elastic_bonds:
        for item in short_elastic_bonds:
            atomid1, atomid2 = item
            resid1, resid2 = (
                atoms[int(atomid1) - 1][2],
                atoms[int(atomid2) - 1][2],
            )
            line = line_style.format(atomid1, atomid2, resid1, resid2)
            exclusions_extend.append(line)

    # Create a backup of the exclusion file and append the new exclusions
    subprocess.run(f'cp {exclusion_file} {exclusion_file + ".bk"}', shell=True)
    with open(exclusion_file, "a+") as f:
        f.writelines(exclusions_extend)

    # Generate Go-model pair interactions for the converted bonds
    # Go-model interactions are attractive interactions between specific atom pairs
    # based on their native distances in the reference structure
    go_pairs_extend: list[str] = []
    u: Universe = mda.Universe(ref_pdb)
    line_style = " {}  {}    1  {:.10f}  {:.10f}  ;  {}  {}  {:.3f}\n"
    if convert_long_elastic_bonds:
        for item, bb_item in zip(ca_long_elastic_bonds, long_elastic_bonds):
            atomid1, atomid2 = int(item[0]), int(item[1])
            atomname1, atomname2 = atoms[atomid1 - 1][1], atoms[atomid2 - 1][1]
            # Calculate distance between CA beads in the reference structure (convert from Å to nm)
            distance = (
                np.linalg.norm(
                    u.atoms[atomid1 - 1].position - u.atoms[atomid2 - 1].position
                )
                / 10
            )
            # Calculate sigma parameter for Lennard-Jones potential
            sigma = distance / 1.12246204830
            line = line_style.format(
                atomname1, atomname2, sigma, lj_epsilon, bb_item[0], bb_item[1], distance
            )
            go_pairs_extend.append(line)

    if convert_short_elastic_bonds:
        for item, bb_item in zip(ca_long_elastic_bonds, long_elastic_bonds):
            atomid1, atomid2 = int(item[0]), int(item[1])
            atomname1, atomname2 = atoms[atomid1 - 1][1], atoms[atomid2 - 1][1]
            # Calculate distance between CA beads in the reference structure (convert from Å to nm)
            distance = (
                np.linalg.norm(
                    u.atoms[atomid1 - 1].position - u.atoms[atomid2 - 1].position
                )
                / 10
            )
            # Calculate sigma parameter for Lennard-Jones potential
            sigma = distance / 1.12246204830
            line = line_style.format(
                atomname1, atomname2, sigma, lj_epsilon, bb_item[0], bb_item[1], distance
            )
            go_pairs_extend.append(line)

    # Create a backup of the Go-pair file and append the new interactions
    subprocess.run(f'cp {go_pair_file} {go_pair_file + ".bk"}', shell=True)
    with open(go_pair_file, "a+") as f:
        f.writelines(go_pairs_extend)

    # Remove the original long and short elastic bonds from the topology file
    # This is necessary because we've converted them to Go-model interactions
    with open(prot_itp, "r") as f:
        lines = f.readlines()

    newlines: list[str] = []
    for line in lines:
        sline = line.split(";")[0].strip().split()
        # Skip long elastic bonds (k=0.970, equilibrium distance=1, force constant=2500)
        if (
            len(sline) == 5
            and (sline[0], sline[1]) in long_elastic_bonds
            and float(sline[2]) == 1
            and float(sline[3]) == 0.970
            and float(sline[4]) == 2500
        ):
            # print(line)
            continue
        # Skip short elastic bonds (k=0.640, equilibrium distance=1, force constant=2500)
        if (
            len(sline) == 5
            and (sline[0], sline[1]) in short_elastic_bonds
            and float(sline[2]) == 1
            and float(sline[3]) == 0.640
            and float(sline[4]) == 2500
        ):
            # print(line)
            continue
        newlines.append(line)

    # Create a backup of the topology file and write the updated version
    subprocess.run(f'cp {prot_itp} {prot_itp + ".bk"}', shell=True)
    with open(prot_itp, "w") as f:
        f.writelines(newlines)

