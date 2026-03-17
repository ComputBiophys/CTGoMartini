# Update: 20250610
# Author: Song Yang

import MDAnalysis as mda
import os
import numpy as np
from optparse import OptionParser
import re
import warnings
warnings.filterwarnings('ignore')

def parse_chain_spec(chain_spec):
    """Parse chain specification string, return a list of chains with residue lists"""
    chains_list = []
    chain_parts = chain_spec.split(';')
    
    for part in chain_parts:
        if not part.strip():
            continue
        
        chain_id, resid_ranges = part.split(':', 1)
        chain_id = chain_id.strip()
        
        chains = {}
        chains[chain_id] = {"residlist": [], "segid": chain_id}
        
        ranges = resid_ranges.split(',')
        for r in ranges:
            r = r.strip()
            if '-' in r:
                start, end = map(int, r.split('-'))
                chains[chain_id]["residlist"].extend(range(start, end + 1))
            else:
                chains[chain_id]["residlist"].append(int(r))
        chains_list.append(chains[chain_id])
    
    return chains_list

# Set up command line parser with improved help messages
parser = OptionParser(usage="usage: %prog [options]",
                     description="Change residue IDs and segment IDs in a PDB file.",
                     epilog="Example: python ChangeResidSegid_input.py -f input.pdb -s \"A:1-20,23-50;B:45-90,92-100\" -o output.pdb")
parser.add_option("-f", "--file", dest="input_file", 
                 help="Input PDB file to be processed")
parser.add_option("-s", "--structure", dest="structure", 
                 help="Chain and residue specification in format 'ChainID:ResidRange[,ResidRange];ChainID:ResidRange'. "
                      "Example: 'A:1-20,23-50;B:45-90,92-100'")
parser.add_option("-o", "--output", dest="output_file", 
                 help="Output PDB file with modified residue and segment IDs")
(options, args) = parser.parse_args()

if not options.input_file or not options.output_file:
    parser.error("Input and output files are required")

str_file = options.input_file
outfile = options.output_file

# Check MDAnalysis version
if mda.__version__ == "2.0.0":
    print("Error: Incompatible MDAnalysis version (2.0.0). Please use a different version.")
    exit()

# Load structure and create a clean PDB
print("\n" + "="*60)
print(f"Processing file: {str_file}")
print("="*60)

u = mda.Universe(str_file)
u.atoms.write('_tmp.pdb')
u = mda.Universe('_tmp.pdb')
os.remove('_tmp.pdb')

# If structure specification is provided, parse it
if options.structure:
    chainlist = parse_chain_spec(options.structure)
    print(f"Using provided chain specification: {options.structure}")
else:
    # Use default chain structure
    print("Error: Chain structure specification (-s) is required")
    raise ValueError("Missing chain structure specification")

# selection
sel = 'protein'

# replace resid
residlist = []
for chain in chainlist:
    residlist += chain["residlist"]

print("\n" + "-"*60)
print("RESIDUE ID INFORMATION:")
print("-"*60)
print(f"New residue count: {len(residlist)}")
print(f"Original protein residue count: {u.select_atoms(sel).residues.n_residues}")

if u.select_atoms(sel).residues.n_residues != len(residlist):
    raise Exception("Error: Mismatch between specified residues and protein residues in structure!")
u.residues[0:len(residlist)].resids = residlist

# add segment
addsegment = []
print("\n" + "-"*60)
print("SEGMENT ID INFORMATION:")
print("-"*60)
print(f"Original segment IDs: {u.segments.segids}")

for chain in chainlist:
    chain_segid = chain["segid"]
    if chain_segid not in u.segments.segids:
        u.add_Segment(segid=chain_segid)
        addsegment.append(chain_segid)

print(f"Added segment IDs: {addsegment if addsegment else 'None'}")
print(f"Final segment IDs: {u.segments.segids}")

# replace segid
index_start = 0
index_end = 0
for chain in chainlist:
    length = len(chain["residlist"])
    index_end += length
    
    # print(chain['segid'])
    segid_index = u.segments.segids == chain['segid']
    u.residues[index_start:index_end].segments = u.segments[segid_index][0]
    u.residues[index_start:index_end].atoms.chainIDs = chain["segid"]

    index_start = index_end

# write pdb
u.atoms.write(outfile)
print("\n" + "-"*60)
print(f"Successfully wrote modified structure to: {outfile}")
print("="*60 + "\n")

