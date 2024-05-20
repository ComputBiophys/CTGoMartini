import os
import MDAnalysis as mda
import numpy as np
from ctgomartini.api import MartiniTopFile
import subprocess


def BB_Distance(atomid1, atomid2, BB_list):
    distance = abs(BB_list.index(atomid1) - BB_list.index(atomid2))
    return distance

def ConvertLongShortElasticBonds(prefix, ref_pdb, convertLongElasticBonds=True, convertShortElasticBonds=False, LJ_epsilon=12.00):
    # Check files
    prot_itp = f"{prefix}.itp"
    go_pair_file=f'{prefix}_go-table_VirtGoSites.itp'
    exclusion_file=f'{prefix}_exclusions_VirtGoSites.itp'
    if not os.path.exists(prot_itp):
        raise ValueError(f'Error: cannot find {prot_itp}!')
    if not os.path.exists(go_pair_file):
        raise ValueError(f'Error: cannot find {go_pair_file}!')
    if not os.path.exists(exclusion_file):
        raise ValueError(f'Error: cannot find {exclusion_file}!')
    if not os.path.exists(ref_pdb):
        raise ValueError(f'Error: cannot find {ref_pdb}!')    
        
    top =MartiniTopFile(prot_itp)
    mol = top._moleculeTypes[prefix]

    # Dict of BB2CA
    BB2CA_dict = {}
    atoms = mol._topology['atoms']
    for fields in mol._topology['virtual_sitesn']:
        if len(fields) == 3 and fields[1] == "1":
            assert atoms[int(fields[0])-1][0] == fields[0] and atoms[int(fields[0])-1][4] == 'CA'
            assert atoms[int(fields[2])-1][0] == fields[2] and atoms[int(fields[2])-1][4] == 'BB'

            if fields[2] not in BB2CA_dict:
                BB2CA_dict[fields[2]] = fields[0]
            else:
                print('Warning: duplicated BB to CA for ',fields)

    # BB list
    BB_list = []
    for fields in atoms:
        if fields[4] == 'BB':
            BB_list.append(fields[0])

    # Find the long bond pairs: k value match
    LongElasticBonds = [] # [(atomid1, atomid2)]
    ShortElasticBonds = []
    bonds =  mol._topology['bonds']
    for fields in bonds:
        # BB distance is 3 and k value is 0.970
        if float(fields[3]) == 0.970 and BB_Distance(fields[0], fields[1], BB_list) == 3:
            atomid1, atomid2 = fields[0], fields[1]
            LongElasticBonds.append((atomid1, atomid2))
        # BB distance is 2 and k value is 0.640
        if float(fields[3]) == 0.640 and BB_Distance(fields[0], fields[1], BB_list) == 2:
            atomid1, atomid2 = fields[0], fields[1]
            ShortElasticBonds.append((atomid1, atomid2))    

    # Convert the BB to the CA
    CA_LongElasticBonds=[]
    for item in LongElasticBonds:
        CA_LongElasticBonds.append((BB2CA_dict[item[0]], BB2CA_dict[item[1]]))
    CA_ShortElasticBonds=[]
    for item in ShortElasticBonds:
        CA_ShortElasticBonds.append((BB2CA_dict[item[0]], BB2CA_dict[item[1]]))

    # Gen the exclusions
    line_style=' {}  {}           ;  {}  {}\n'
    exclusions_extend=[]
    if convertLongElasticBonds:
        for item in LongElasticBonds:
            atomid1, atomid2=item
            resid1, resid2=atoms[int(atomid1)-1][2], atoms[int(atomid2)-1][2]
            line=line_style.format(atomid1, atomid2, resid1, resid2)
            exclusions_extend.append(line)

    if convertShortElasticBonds:
        for item in ShortElasticBonds:
            atomid1, atomid2=item
            resid1, resid2=atoms[int(atomid1)-1][2], atoms[int(atomid2)-1][2]
            line=line_style.format(atomid1, atomid2, resid1, resid2)
            exclusions_extend.append(line)

    # Append exclusion file
    subprocess.run(f'cp {exclusion_file} {exclusion_file+".bk"}', shell=True)
    with open(exclusion_file,'a+') as f:
        f.writelines(exclusions_extend)

    # Gen go pairs
    go_pairs_extend = []
    u = mda.Universe(ref_pdb)
    line_style = ' {}  {}    1  {:.10f}  {:.10f}  ;  {}  {}  {:.3f}\n'
    if convertLongElasticBonds:
        for item, BB_item in zip(CA_LongElasticBonds, LongElasticBonds):
            atomid1, atomid2 = int(item[0]), int(item[1])
            atomname1, atomname2 = atoms[atomid1-1][1], atoms[atomid2-1][1]
            distance = np.linalg.norm(u.atoms[atomid1-1].position-u.atoms[atomid2-1].position)/10
            sigma = distance/1.12246204830
            line = line_style.format(atomname1, atomname2, sigma, LJ_epsilon, BB_item[0], BB_item[1], distance)
            go_pairs_extend.append(line)

    if convertShortElasticBonds:
        for item, BB_item in zip(CA_LongElasticBonds, LongElasticBonds):
            atomid1, atomid2 = int(item[0]), int(item[1])
            atomname1, atomname2 = atoms[atomid1-1][1], atoms[atomid2-1][1]
            distance = np.linalg.norm(u.atoms[atomid1-1].position-u.atoms[atomid2-1].position)/10
            sigma = distance/1.12246204830
            line = line_style.format(atomname1, atomname2, sigma, LJ_epsilon, BB_item[0], BB_item[1], distance)
            go_pairs_extend.append(line)

    # Append go_pair file
    subprocess.run(f'cp {go_pair_file} {go_pair_file+".bk"}', shell=True)
    with open(go_pair_file,'a+') as f:
        f.writelines(go_pairs_extend)

    # Delte the Long/Short Elastic Bonds 
    with open(prot_itp, 'r') as f:
        lines =  f.readlines()

    newlines  = []
    for line in lines:
        sline= line.split(';')[0].strip().split()
        if len(sline) == 5 and (sline[0], sline[1]) in LongElasticBonds \
            and float(sline[2]) == 1 and float(sline[3]) == 0.970 and float(sline[4]) == 2500:
            # print(line)
            continue
        if len(sline) == 5 and (sline[0], sline[1]) in ShortElasticBonds \
            and float(sline[2]) == 1 and float(sline[3]) == 0.640 and float(sline[4]) == 2500:
            # print(line)
            continue
        newlines.append(line)

    subprocess.run(f'cp {prot_itp} {prot_itp+".bk"}', shell=True)
    with open(prot_itp, 'w') as f:
        f.writelines(newlines)

