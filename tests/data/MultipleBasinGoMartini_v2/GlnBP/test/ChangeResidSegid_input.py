# Update: 20220608: can be used in version 2.1.0

import MDAnalysis as mda
import os
import numpy as np
from optparse import OptionParser

parser = OptionParser()
(options, args) = parser.parse_args()
print(args)
name0=args[0]
name1=args[1]

if mda.__version__=="2.0.0":
    print("Something wrong with the version of MDAnalysis")
    exit()

str_file=name0
outfile=name1
u=mda.Universe(str_file)

chainA={"residlist":[i for i in range(5,225)],
       "segid":"A"}
others={"residlist":[i+1 for i in range(u.select_atoms('not protein').residues.n_residues)],
        "segid":"O"}

chainlist=[chainA, chainA] #, others ]

# selection
sel='protein'

# replace resid
residlist=[]
for chain in chainlist:
    residlist+=chain["residlist"]
print("New resid length is {}".format(len(residlist)))
print("Old resid length is {}".format(u.select_atoms(sel).residues.n_residues))
if u.select_atoms(sel).residues.n_residues != len(residlist):
    raise Exception("Error: There is something wrong with resids!")
    #print("There is something wrong with resids!")
u.residues[0:len(residlist)].resids=residlist
print("\n")


# add segment
addsegment=[]
print("Old segment segids: {}".format(u.segments.segids))
for chain in chainlist:
    chain_segid=chain["segid"]
    if chain_segid not in u.segments.segids:
        u.add_Segment(segid=chain_segid)
        addsegment.append(chain_segid)
print("Add segment segids: {}".format(addsegment))
print("New segment segids: {}".format(u.segments.segids))
print('\n')

# replace segid
index_start=0
index_end=0
for chain in chainlist:
    length=len(chain["residlist"])
    index_end+=length
    
    segid_index=u.segments.segids==chain['segid']
    u.residues[index_start:index_end].segments=u.segments[segid_index][0]
    u.residues[index_start:index_end].atoms.chainIDs=chain["segid"]
    index_start=index_end
    
    

# write pdb
u.atoms.write(outfile)
