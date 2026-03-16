"""OVrCSU - Protein Contact Surface Unit Analysis Tool

This script calculates atom-atom and residue-residue contacts in protein structures.

Version: 20250323
Author: Song Yang
Email: yangsong2015@pku.edu.cn

References:
    [1] Wołek, K.; Gómez-Sicilia, À.; Cieplak, M. Determination of contact maps in proteins: 
        A combination of structural and chemical approaches. The Journal of Chemical Physics 2015, 143.
    [2] Yang, S.; Song, C. Multiple-basin Go-Martini for investigating conformational transitions 
        and environmental interactions of proteins. bioRxiv 2024.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from typing import Any

import MDAnalysis as mda
import numpy as np
from MDAnalysis.core.groups import Atom, Residue
from MDAnalysis.core.universe import Universe
from MDAnalysis.lib.distances import distance_array

# Global parameters
default_dtype = np.float64  # Default floating-point precision for calculations
R_WATER = 1.42  # Radius of water molecule in Angstroms (Å)
FIB_NUMBER = 14  # Number of Fibonacci iterations for sphere point sampling
ENLARGEMENT_FACTOR = 1.244455060259808  # Scaling factor for van der Waals radii (26/7)**(1/6)


parser = argparse.ArgumentParser(
    description="""
An example:
python OVrCSU.py -f strfile.pdb -o new_OV.map

""",
    formatter_class=argparse.RawTextHelpFormatter,
)

parser.add_argument("-f", "--strfile", required=True, help="The structure file (pdb)")
parser.add_argument("-o", "--fileout", default="new_OV.map", help="The OV contact filename")

args = parser.parse_args()
pdb_file: str = args.strfile
outfile: str = args.fileout


def CalFibonacci(fib: int = 14, fiba: int = 0, fibb: int = 1) -> tuple[int, int]:
    """Calculate Fibonacci numbers up to specified index.

    Args:
        fib: Number of Fibonacci iterations (default 14).
        fiba: Starting first number (default 0).
        fibb: Starting second number (default 1).

    Returns:
        A tuple containing Fibonacci numbers at positions fib and fib+1.
    """
    for _ in range(fib):
        fibc = fiba + fibb
        fiba = fibb
        fibb = fibc
    return fiba, fibb


def fibonacci_sphere(fiba: int, fibb: int, radius: float) -> np.ndarray:
    """Generate evenly distributed points on a sphere using Fibonacci sampling.

    Args:
        fiba: Fibonacci number F(n).
        fibb: Fibonacci number F(n+1).
        radius: Sphere radius.

    Returns:
        Array of 3D coordinates (x, y, z) for sampling points.
    """
    indices = np.arange(1, fibb + 1, dtype=default_dtype)
    theta = np.arccos(1 - 2 * (indices) / fibb)
    phi = np.pi * 2 * indices * fiba / fibb
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    return np.stack([x, y, z], axis=1)


def RemoveH(PDBFile: str) -> None:
    """Process PDB file by removing hydrogen atoms and OXT (terminal oxygen).

    Args:
        PDBFile: Path to input PDB file.

    Writes:
        _protein_noH.pdb: Processed PDB file without hydrogens.
    """
    u = mda.Universe(PDBFile)
    sel = u.select_atoms(
        "protein and not name H* OXT"
    )  # Select protein atoms excluding hydrogens and OXT
    sel.atoms.write("_protein_noH.pdb")


RemoveH(pdb_file)
strfile = "_protein_noH.pdb"
u: Universe = mda.Universe(strfile)
atoms = u.atoms
n_atoms: int = len(atoms)


vdw_dict: dict[str, dict[str, float]] = {
    "SER": {"N": 1.64, "CA": 1.88, "C": 1.61, "O": 1.42, "CB": 1.88, "OG": 1.46},
    "PRO": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.88,
        "CD": 1.88,
    },
    "TYR": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.61,
        "CD1": 1.76,
        "CD2": 1.76,
        "CE1": 1.76,
        "CE2": 1.76,
        "CZ": 1.61,
        "OH": 1.46,
    },
    "VAL": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG1": 1.88,
        "CG2": 1.88,
    },
    "TRP": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.61,
        "CD1": 1.76,
        "CD2": 1.61,
        "CE2": 1.61,
        "CE3": 1.76,
        "NE1": 1.64,
        "CZ2": 1.76,
        "CZ3": 1.76,
        "CH2": 1.76,
    },
    "GLN": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.88,
        "CD": 1.61,
        "NE2": 1.64,
        "OE1": 1.42,
    },
    "HIS": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.61,
        "CD2": 1.76,
        "ND1": 1.64,
        "CE1": 1.76,
        "NE2": 1.64,
    },
    "GLU": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.88,
        "CD": 1.61,
        "OE1": 1.46,
        "OE2": 1.42,
    },
    "LYS": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.88,
        "CD": 1.88,
        "CE": 1.88,
        "NZ": 1.64,
    },
    "THR": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG2": 1.88,
        "OG1": 1.46,
    },
    "MET": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.88,
        "SD": 1.77,
        "CE": 1.88,
    },
    "ALA": {"N": 1.64, "CA": 1.88, "C": 1.61, "O": 1.42, "CB": 1.88},  # 'OXT': 1.46
    "CYS": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "SG": 1.77,
    },
    "ASP": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.61,
        "OD1": 1.46,
        "OD2": 1.42,
    },
    "ASN": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.61,
        "ND2": 1.64,
        "OD1": 1.42,
    },
    "ARG": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.88,
        "CD": 1.88,
        "NE": 1.64,
        "CZ": 1.61,
        "NH1": 1.64,
        "NH2": 1.64,
    },
    "PHE": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.88,
        "CD1": 1.61,
        "CD2": 1.76,
        "CE1": 1.76,
        "CE2": 1.76,
        "CZ": 1.76,
    },
    "ILE": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG1": 1.88,
        "CG2": 1.88,
        "CD1": 1.88,
    },
    "GLY": {"N": 1.64, "CA": 1.88, "C": 1.61, "O": 1.42},
    "LEU": {
        "N": 1.64,
        "CA": 1.88,
        "C": 1.61,
        "O": 1.42,
        "CB": 1.88,
        "CG": 1.88,
        "CD1": 1.88,
        "CD2": 1.88,
    },
}


atom_type_dict: dict[str, dict[str, int]] = {
    "SER": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 6, "OG": 1},
    "PRO": {"N": 6, "CA": 4, "C": 6, "O": 2, "CB": 4, "CG": 4, "CD": 4},
    "TYR": {
        "N": 3,
        "CA": 7,
        "C": 6,
        "O": 2,
        "CB": 4,
        "CG": 5,
        "CD1": 5,
        "CD2": 5,
        "CE1": 5,
        "CE2": 5,
        "CZ": 5,
        "OH": 1,
    },
    "VAL": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 4, "CG1": 4, "CG2": 4},
    "TRP": {
        "N": 3,
        "CA": 7,
        "C": 6,
        "O": 2,
        "CB": 4,
        "CG": 5,
        "CD1": 5,
        "CD2": 5,
        "CE2": 5,
        "CE3": 5,
        "NE1": 3,
        "CZ2": 5,
        "CZ3": 5,
        "CH2": 5,
    },
    "GLN": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 4, "CG": 4, "CD": 6, "NE2": 3, "OE1": 2},
    "HIS": {
        "N": 3,
        "CA": 7,
        "C": 6,
        "O": 2,
        "CB": 4,
        "CG": 5,
        "CD2": 5,
        "ND1": 1,
        "CE1": 5,
        "NE2": 1,
    },
    "GLU": {
        "N": 3,
        "CA": 7,
        "C": 6,
        "O": 2,
        "CB": 4,
        "CG": 4,
        "CD": 6,
        "OE1": 10,
        "OE2": 10,
    },
    "LYS": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 4, "CG": 4, "CD": 4, "CE": 7, "NZ": 9},
    "THR": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 6, "CG2": 4, "OG1": 1},
    "MET": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 4, "CG": 4, "SD": 8, "CE": 4},
    "ALA": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 4},  # 'OXT': 1
    "CYS": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 4, "SG": 6},
    "ASP": {
        "N": 3,
        "CA": 7,
        "C": 6,
        "O": 2,
        "CB": 4,
        "CG": 6,
        "OD1": 10,
        "OD2": 10,
    },
    "ASN": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 4, "CG": 6, "ND2": 3, "OD1": 2},
    "ARG": {
        "N": 3,
        "CA": 7,
        "C": 6,
        "O": 2,
        "CB": 4,
        "CG": 4,
        "CD": 7,
        "NE": 3,
        "CZ": 6,
        "NH1": 9,
        "NH2": 9,
    },
    "PHE": {
        "N": 3,
        "CA": 7,
        "C": 6,
        "O": 2,
        "CB": 4,
        "CG": 5,
        "CD1": 5,
        "CD2": 5,
        "CE1": 5,
        "CE2": 5,
        "CZ": 5,
    },
    "ILE": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 4, "CG1": 4, "CG2": 4, "CD1": 4},
    "GLY": {"N": 3, "CA": 6, "C": 6, "O": 2},
    "LEU": {"N": 3, "CA": 7, "C": 6, "O": 2, "CB": 4, "CG": 4, "CD1": 4, "CD2": 4},
}


interaction_matrix: np.ndarray = np.array([
    [1, 1, 1, 5, 5, 6, 6, 6, 1, 1],  # Class I
    [1, 5, 1, 5, 5, 6, 6, 6, 1, 5],  # Class II
    [1, 1, 5, 5, 5, 6, 6, 6, 5, 1],  # Class III
    [5, 5, 5, 2, 2, 6, 6, 6, 5, 5],  # Class IV
    [5, 5, 5, 2, 3, 6, 6, 6, 5, 5],  # Class V
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],  # Class VI
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],  # Class VII
    [6, 6, 6, 6, 6, 6, 6, 6, 6, 6],  # Class VIII
    [1, 1, 5, 5, 5, 6, 6, 6, 5, 4],  # Class IX
    [1, 5, 1, 5, 5, 6, 6, 6, 4, 5],  # Class X
])

interaction_type_value_dict: dict[int, int] = {
    1: 1,
    2: 1,
    3: 1,
    4: 1,
    5: -1,
    6: 0,
}

# Preprocess atom types and van der Waals radii
atom_types: list[int] = [atom_type_dict[atom.resname][atom.name] for atom in atoms]
vdw_radii: list[float] = [vdw_dict[atom.resname][atom.name] for atom in atoms]

# Step 1: Calculate candidate atom pairs
positions = atoms.positions
positions = np.array(positions, dtype=default_dtype)  # Convert to numpy array
vdw_radii_arr: np.ndarray = np.array(vdw_radii, dtype=default_dtype)  # Convert to numpy array
cutoffs_matrix: np.ndarray = vdw_radii_arr[:, None] + vdw_radii_arr[None, :] + 2 * R_WATER
dist_matrix: np.ndarray = distance_array(positions, positions)
dist_matrix = np.array(dist_matrix, dtype=default_dtype)  # Convert to numpy array
candidate_pairs: np.ndarray = np.argwhere(
    (dist_matrix <= cutoffs_matrix) & (np.triu(np.ones_like(dist_matrix), k=1).astype(bool))
)


# Build atom neighbor dictionary
atom_neighbors: defaultdict[int, list[int]] = defaultdict(list)
for i, j in candidate_pairs:
    atom_neighbors[int(i)].append(int(j))
    atom_neighbors[int(j)].append(int(i))

# Step 2: Vectorized calculation of grid interactions
interaction_pairs: defaultdict[tuple[int, int], dict[str, int | float]] = defaultdict(
    lambda: {"DISTANCE": 0, "Surf": 0, "S0": 0, "Cont": 0}
)

fiba, fibb = CalFibonacci(fib=FIB_NUMBER, fiba=0, fibb=1)
for i in atom_neighbors:
    neighbors = atom_neighbors[i]
    if not neighbors:
        continue

    pos_i = positions[i]
    radius_i = vdw_radii_arr[i] + R_WATER
    fib_points = fibonacci_sphere(fiba, fibb, radius_i)
    global_points = pos_i + fib_points

    # Process neighbor data
    neighbor_pos = positions[neighbors]
    neighbor_vdw = vdw_radii_arr[neighbors] + R_WATER
    dist_atoms = np.linalg.norm(pos_i - neighbor_pos, axis=1)

    # Vectorized distance calculation
    diffs = global_points[:, np.newaxis, :] - neighbor_pos[np.newaxis, :, :]
    dists = np.linalg.norm(diffs, axis=2)
    mask = dists <= neighbor_vdw[np.newaxis, :]

    # Process each sampling point
    for s in range(fibb):
        valid_js = np.where(mask[s])[0]
        if valid_js.size == 0:
            continue

        # Find the index of the nearest atom
        min_j_idx = valid_js[np.argmin(dist_atoms[valid_js])]
        j = neighbors[min_j_idx]

        # Count the number of interactions
        interaction_pairs[(i, j)]["Cont"] += 1  # type: ignore[operator]

# Complete the interaction pairs
for (i, j) in list(interaction_pairs.keys()):
    interaction_pairs[(i, j)]["DISTANCE"] = dist_matrix[i, j]  # type: ignore[index]

    radius_i = vdw_radii_arr[i] + R_WATER
    interaction_pairs[(i, j)]["Surf"] = np.array(  # type: ignore[index]
        int(interaction_pairs[(i, j)]["Cont"]) * 4 * np.pi * radius_i**2 / fibb,  # type: ignore[arg-type]
        dtype=default_dtype,
    )
    # interaction_pairs[(i,j)]['S0'] = 4 * np.pi * radius_i**2


lines: list[str] = []
lines.append(
    """
Atom - atom contacts


I1,I2 - id of residues in contact

A1,A2 - names of atoms in contact

C - chain

T - id of atom class

I - id of interaction type

Surf - surface of interaction in CSU           algorithm

S0 - whole surface of overlap in CSU (Not supported)

Cont - number of contacts between atoms
"""
)
lines.append("    I1  AA    C A1   T       I2  AA    C A2   T     I   DISTANCE       Surf      S0    Cont\n")
lines.append("============================================================================================\n")

sorted_pairs = sorted(interaction_pairs.keys())
for i, pair in enumerate(sorted_pairs):
    value = interaction_pairs[pair]

    atom1: Atom = atoms[pair[0]]
    atom2: Atom = atoms[pair[1]]
    residue1: Residue = atom1.residue
    residue2: Residue = atom2.residue

    I1 = residue1.resid
    I2 = residue2.resid

    if residue1.ix == residue2.ix:
        continue

    AA1 = residue1.resname
    AA2 = residue2.resname

    C1 = residue1.segid
    C2 = residue2.segid

    A1 = atom1.name
    A2 = atom2.name

    T1 = atom_type_dict[residue1.resname][atom1.name]
    T2 = atom_type_dict[residue2.resname][atom2.name]

    I = interaction_matrix[T1 - 1][T2 - 1]

    DISTANCE = float(interaction_pairs[pair]["DISTANCE"])
    Surf = float(interaction_pairs[pair]["Surf"])
    S0 = float(interaction_pairs[pair]["S0"])
    Cont = int(interaction_pairs[pair]["Cont"])

    line = f"A {I1:>4} {AA1:>4}   {C1} {A1:<4}  {T1}     {I2:>4} {AA2:>4}   {C2} {A2:<4}  {T2:>2}    {I:>4} {DISTANCE:>12.4f} {Surf:>12.4f} {S0:>12.4f} {Cont:>8}\n"
    # line=f"A  {I1:>5} {AA1:>4} {C1} {I1_pdb:>4} {I2:>8} {AA2:>4} {C2} {I2_pdb:>4} {dist:>12.4f}     0 0 0 {is_rCSU}     {n_contats:>4}   0.0000   0.0000   0.0000\n"
    lines.append(line)
    # print(line)


# Write Contact
outfile = outfile
with open(outfile, "w") as fp:
    fp.writelines(lines)

# Residue pair interactions
residue_interactions: defaultdict[tuple[Residue, Residue], dict[str, float]] = defaultdict(
    lambda: {"rCSU": 0, "aSurf": 0, "rSurf": 0, "nSurf": 0, "OV": 0}
)

for (i, j), counts in interaction_pairs.items():
    # Don't include the asymmetric pairs
    if (j, i) not in interaction_pairs:
        # continue
        pass

    res_i: Residue = atoms[i].residue
    res_j: Residue = atoms[j].residue

    # Skip interactions between the same residue
    if res_i.ix == res_j.ix:
        continue

    # Interaction types
    type_i = atom_types[i]
    type_j = atom_types[j]
    interaction_type = int(interaction_matrix[type_i - 1][type_j - 1])
    interaction = interaction_type_value_dict[interaction_type]

    key = (res_i, res_j)
    residue_interactions[key]["rCSU"] += interaction * int(interaction_pairs[(i, j)]["Cont"])
    if interaction == 1:
        residue_interactions[key]["aSurf"] += float(interaction_pairs[(i, j)]["Surf"])
    elif interaction == -1:
        residue_interactions[key]["rSurf"] += float(interaction_pairs[(i, j)]["Surf"])
    elif interaction == 0:
        residue_interactions[key]["nSurf"] += float(interaction_pairs[(i, j)]["Surf"])


# OV Calculations
# positions = atoms.positions
OV_cutoffs_matrix: np.ndarray = (vdw_radii_arr[:, None] + vdw_radii_arr[None, :]) * ENLARGEMENT_FACTOR
# dist_matrix = distance_array(positions, positions)
# dist_matrix = np.array(dist_matrix, dtype=default_dtype)  # Convert to numpy array
OV_candidate_pairs: np.ndarray = np.argwhere(
    (dist_matrix <= OV_cutoffs_matrix) & (np.triu(np.ones_like(dist_matrix), k=1).astype(bool))
)
for (i, j) in OV_candidate_pairs:
    atom1: Atom = atoms[i]
    atom2: Atom = atoms[j]
    residue1: Residue = atom1.residue
    residue2: Atom = atom2.residue
    if residue1.ix == residue2.ix:
        continue
    key = (residue1, residue2)
    residue_interactions[key]["OV"] += 1
    key = (residue2, residue1)
    residue_interactions[key]["OV"] += 1


# Output
def ExtactContacts(contact_pair: list[int], u: Universe, sel: str = "CA") -> float:
    """Extract distance between specified atoms of two residues.

    Args:
        contact_pair: List containing indices of two residues [resindex1, resindex2].
        u: MDAnalysis Universe object containing the molecular structure.
        sel: Atom name to use for distance calculation (default "CA" for alpha carbon).

    Returns:
        Distance between the selected atoms of the two residues in Angstroms.
    """
    resindex1 = contact_pair[0]
    resindex2 = contact_pair[1]

    position1 = u.residues[resindex1].atoms[
        u.residues[resindex1].atoms.names == sel
    ].positions
    position2 = u.residues[resindex2].atoms[
        u.residues[resindex2].atoms.names == sel
    ].positions

    dist: float = distance_array(position1, position2)[0][0]
    return dist


lines = []
lines.append(
    """
Resideue residue contacts

I1,I1 - residue id
I(PDB) - residue number in PDB file
C - chain
AA - 3-letter code of aminoacid
DISTANCE - distance between CA
CMs - contacts in particular contactmaps, OV, CSU, oCSU, rCSU respectively, (CSU do not take into
account chemical properties of atoms) (Not supported CSU, oCSU)
aSurf - surface of attractive connections
rSurf - surface of repulsive connections
nSurf - surface of neutral connections

"""
)
lines.append(
    "            I1  AA  C I(PDB)    I2  AA  C I(PDB)    DISTANCE       CMs    rCSU    aSurf    rSurf    nSurf\n"
)
lines.append(
    "==========================================================================================================\n"
)
sorted_pairs = sorted(residue_interactions.keys())
for i, pair in enumerate(sorted_pairs):
    value = residue_interactions[pair]

    residue1 = pair[0]
    residue2 = pair[1]

    I1 = residue1.ix + 1
    I2 = residue2.ix + 1

    AA1 = residue1.resname
    AA2 = residue2.resname

    C1 = residue1.segid
    C2 = residue2.segid

    I1_pdb = residue1.resid
    I2_pdb = residue2.resid

    dist = ExtactContacts([residue1.ix, residue2.ix], u, sel="CA")

    rCSU = value["rCSU"]
    aSurf = value["aSurf"]
    rSurf = value["rSurf"]
    nSurf = value["nSurf"]
    OV = value["OV"]
    is_rCSU = 1 if rCSU > 0 else 0
    is_OV = 1 if OV > 0 else 0
    if is_OV == 0 and is_rCSU == 0:
        continue

    line = f"R {i+1:>6} {I1:>5} {AA1:>4} {C1} {I1_pdb:>4} {I2:>8} {AA2:>4} {C2} {I2_pdb:>4} {dist:>12.4f}     {is_OV} 0 0 {is_rCSU}     {rCSU:>4}   {aSurf:>8.4f}   {rSurf:>8.4f}   {nSurf:>8.4f}\n"
    lines.append(line)


# Write Contact
outfile = outfile
with open(outfile, "a") as fp:
    fp.writelines(lines)
