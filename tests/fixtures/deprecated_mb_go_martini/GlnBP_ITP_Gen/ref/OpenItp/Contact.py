import subprocess
import MDAnalysis as mda
import numpy as np
import re

# Input parameters
prefix='gbp_open'
strfile='../GBP_open_cg.pdb'


def SplitLines(lines,pattern='\[ (.+) \]'):
    startlist=[]
    endlist=[]
    for i,line in enumerate(lines):
        if len(re.findall(pattern, line))!=0:
            startlist.append(i)
            endlist.append(i-1)
    endlist.append(i)
    endlist=endlist[1:]

    Data_Dict={}
    for pair in zip(startlist,endlist):
        key=re.findall(pattern, lines[pair[0]])[0]

        value=lines[pair[0]+1:pair[1]+1]
        value=ProcessLines(value)
        if key in Data_Dict:
            Data_Dict[key]+=value
        else:
            Data_Dict[key]=value
    return Data_Dict


def ReadItp(file):
    with open(file,'r') as f:
        lines=f.readlines()
    newlines=[]
    for line in lines:
        line=line[:line.find(';')] if line.find(';') != -1 else line
        sline=line.strip()
        if sline=='': continue

        newlines.append(sline)

    Data_Dict=SplitLines(newlines,pattern='\[ (.+) \]')
    return Data_Dict
def ChangeTypes(alist):
    newlist=[]
    for item in alist:
        if re.match('\-?\d+\.\d',item):
            item=float(item)
        elif re.match('^\-?\d+$',item):
            item=int(item)
        else:
            item=item
        newlist.append(item)
    return newlist
def ProcessLines(lines):
    newlines=[]
    for line in lines:
        sline=line.strip().split()
        sline=ChangeTypes(sline)
        newlines.append(sline)
    return newlines


def WriteItp(itp_dict, outfile):
    lines=['\n']
    for key, values in itp_dict.items():
        lines+=[f'[ {key} ]\n']
        for value in values:
            line=' '.join(map(str, value))+'\n'
            lines+=[line]

        lines+=['\n']

    with open(outfile, 'w') as fp:
        fp.writelines(lines)
    return lines


# Load ITP dict
Data_Dict=ReadItp(f'{prefix}.itp')
go_pair_file=f'{prefix}_go-table_VirtGoSites.itp'
exclusion_file=f'{prefix}_exclusions_VirtGoSites.itp'


# Dict of the BB to the CA
BB2CA_dict={}
for item in Data_Dict['virtual_sitesn']:
    if len(item)==3 and item[1]==1:
        if item[2] not in BB2CA_dict:
            BB2CA_dict[item[2]]=item[0]
        else:
            print('Warning: duplicated BB to CA for ',item)

# Check
for key, value in BB2CA_dict.items():
    if Data_Dict['atoms'][key-1][2]!=Data_Dict['atoms'][value-1][2]:
        print('Error: BB2CA projection error!', key, value)


# Find the long bond pairs
LongBond_record=[] # [(atomid1, atomid2)]
for item in Data_Dict['bonds']:
    if len(item)!=5: continue
    if item[3]== 0.970:
        atomid1, atomid2=item[0], item[1]
        LongBond_record.append((atomid1, atomid2))

# Convert the BB to the CA
CA_LongBond_record=[]
for item in LongBond_record:
    CA_LongBond_record.append((BB2CA_dict[item[0]], BB2CA_dict[item[1]]))


# Gen the exclusions
line_style=' {}  {}           ;  {}  {}\n'
exclusions=[]
for item in LongBond_record:
    atomid1, atomid2=item
    resid1, resid2=Data_Dict['atoms'][atomid1-1][2], Data_Dict['atoms'][atomid2-1][2]
    line=line_style.format(atomid1, atomid2, resid1, resid2)
    exclusions.append(line)
    
# Append exclusion file
subprocess.run(f'cp {exclusion_file} {exclusion_file+".bk"}', shell=True)
with open(exclusion_file,'a+') as f:
    f.writelines(exclusions)


# Gen go pairs
go_pairs=[]
u=mda.Universe(strfile)
line_style=' {}  {}    1  {:.10f}  12.0000000000  ;  {}  {}  {:.3f}\n'
for item, BB_item in zip(CA_LongBond_record, LongBond_record):
    atomid1, atomid2=item
    atomname1, atomname2=Data_Dict['atoms'][atomid1-1][1], Data_Dict['atoms'][atomid2-1][1]
    distance=np.linalg.norm(u.atoms[atomid1-1].position-u.atoms[atomid2-1].position)/10
    sigma=distance/1.12246204830
    line=line_style.format(atomname1, atomname2, sigma, BB_item[0], BB_item[1], distance)
    go_pairs.append(line)

# Append go_pair file
subprocess.run(f'cp {go_pair_file} {go_pair_file+".bk"}', shell=True)
with open(go_pair_file,'a+') as f:
    f.writelines(go_pairs)
