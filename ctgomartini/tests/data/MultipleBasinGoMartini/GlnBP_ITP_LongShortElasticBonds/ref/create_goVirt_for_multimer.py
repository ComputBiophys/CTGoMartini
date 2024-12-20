#!/usr/bin/env python
# Date: 20221011 
# Author: Song Yang

#: Bug: extract is not good for CA


import argparse
import subprocess
import numpy as np

def user_input():
    parser = argparse.ArgumentParser('''
    Example:
    python create_goVir_for_multimer.py -r protein.pdb -s protein_cg.pdb -f protein.map --moltype pro --go_eps 12 --Natoms 100
    ''')
    parser.add_argument('-s', help='File containing the coarse-grained structure of the protein in pdb format.')
    parser.add_argument('-r', help='File containing the all-atomed structure of the protein in pdb format. The pdb file is the one used to generate OV+rCSU map file.')
    parser.add_argument('-f', help='File containing the contact analysis of the (atomistic) protein structure obtained from the webserver http://info.ifpan.edu.pl/~rcsu/rcsu/index.html.')
    parser.add_argument('--moltype', default='molecule_0', help='Molecule name used as prefix in your output file names and the virtual bead names (default: molecule_0). If you will combine your Go-like model with a coarse-grained protein generated with martinize2, you must use the same name as specified with the --govs-moltype flag of martinize2!')
    parser.add_argument('--go_eps', type=float, default=9.414, help='Dissociation energy [kJ/mol] of the Lennard-Jones potential used in the Go-like model (default: 9.414).')
    parser.add_argument('--cutoff_short', type=float, default=0.3, help='Lower cutoff distance [nm]: contacts with a shorter distance than cutoff_short are not included in the Go-like interactions (default: 0.3).')
    parser.add_argument('--cutoff_long', type=float, default=1.1, help='Upper cutoff distance [nm]: contacts with a longer distance than cutoff_long are not included in the Go-like interactions (default: 1.1).')
    parser.add_argument('--Natoms', type=int, help='Number of coarse-grained beads in the protein excluding the virtual Go beads.')
    args = parser.parse_args()
    return args

def get_settings():
    # names for temporary files:
    file_CA = 'CA.pdb'
    file_BB = 'BB.pdb'
    file_OV = 'OV.map'
    file_rCSU = 'rCSU.map'

    # some other rudimentary variables
    seqDist = 4               # minimal distance in the sequence to add a elastic bond (ElNedyn=3 [Perriole2009]; Go=4 [Poma2017])
    cols = [4, 5, 8, 9, 10]   # colums of interest in the OV and rCSU contact map file. 
                              # chainid A, resid A, chainid B, resid B, contact distance
    c6c12 = 0                 # if set to 1, the C6 and C12 term are expected in the .itp file; if set to 0, sigma and go_eps are used
    return file_CA, file_BB, file_OV, file_rCSU, seqDist, cols, c6c12


def read_data(struct_pdb, ref_pdb, file_contacts, file_CA, file_BB, file_OV, file_rCSU, cols):
    # preparation of temporary files for reading
    subprocess.call("grep -a '1 [01] [01] [01]' " + file_contacts + " > " + file_OV, shell=True)
    subprocess.call("echo '' >> " + file_OV, shell=True)
    subprocess.call("grep -a '0 [01] [01] 1' " + file_contacts + " > " + file_rCSU, shell=True)
    subprocess.call("echo '' >> " + file_rCSU, shell=True)
    subprocess.call("grep -a 'BB' " + struct_pdb + " > " + file_BB, shell=True)
    subprocess.call("echo '' >> " + file_BB, shell=True)
    subprocess.call("grep -a 'CA ' " + ref_pdb + " > " + file_CA, shell=True)
    subprocess.call("echo '' >> " + file_CA, shell=True)

    # read coarse-grained BB bead positions
    with open(file_BB,'r') as fid:
        dat = fid.readlines()
    dat = dat[:-1]
    print('Number of coarse-grained BB beads in your protein: ' + str(len(dat)))

    indBB = [] # newresid, newatomid, position x,y,z, 
    nameAA = [] # residue name

    for k in range(0, len(dat)):
        tmp = dat[k]
        tmp=tmp[:22]+' '+tmp[22:38]+' '+tmp[38:46]+' '+tmp[46:54]+' '+tmp[54:] # add one gap between chain id, new resid, positions, and others.
        tmp = tmp.split()
        indBB.append([int(tmp[5]), int(tmp[1]), float(tmp[6]), float(tmp[7]), float(tmp[8])]) # newresid, newatomid, px, py, pz
        nameAA.append(tmp[3]) # resname

    indBB = np.array(indBB) # newresid, newatomid, px, py, pz

    # read file_CA 
    with open(file_CA,'r') as fid:
        dat=fid.readlines()
    dat = dat[:-1]
    order = {} # (chainid, resid):[index_resid, newresid]     
    
    ## test the equality of Ref pdb and CG pdb
    if len(dat) != len(indBB):
        print("Error: Ref pdb does not equal to CG pdb")
        return 
    
    for k in range(0, len(dat)):
        tmp = dat[k]
        chainid = tmp[21]
        resid = int(tmp[22:26])
        key=(chainid, resid)
        value = [k, indBB[k,0]] # considered that the orders of CG file and ref pdbfile are same
        order[key] = value
        
    # read OV contact map
    with open(file_OV,'r') as fid:
        dat = fid.readlines()
    dat = dat[:-1]
    print('Number of contacts read from your OV contact map file: ' + str(len(dat)))

    map_OVrCSU = []
    for k in range(0, len(dat)):
        row = []
        tmp = dat[k]
        tmp = tmp.replace('\t',' ')
        tmp = tmp.split()
        for l in cols:
            row.append(tmp[l]) # chainid A, resid A, chainid B, resid B, contact distance
        row[1]=int(row[1])
        row[3]=int(row[3])
        row[4]=float(row[4])         
        map_OVrCSU.append(row)

    # read rCSU contact map
    with open(file_rCSU,'r') as fid:
        dat = fid.readlines()
    dat = dat[:-1]
    print('Number of contacts read from your rCSU contact map file: ' + str(len(dat)))

    for k in range(0, len(dat)):
        row = []
        tmp = dat[k]
        tmp = tmp.replace('\t',' ')
        tmp = tmp.split()
        for l in cols:
            row.append(tmp[l]) # chainid A, resid A, chainid B, resid B, contact distance
        row[1]=int(row[1])
        row[3]=int(row[3])
        row[4]=float(row[4]) 
        map_OVrCSU.append(row)

    return indBB, order, nameAA, map_OVrCSU


def get_go(indBB, order, nameAA, map_OVrCSU, cutoff_short, cutoff_long, go_eps, seqDist):
    # calculate the distances based on the coordinates of the CG BB bead
    for k in range(0, len(map_OVrCSU)):
        chainid_resid_1=(map_OVrCSU[k][0],map_OVrCSU[k][1])
        chainid_resid_2=(map_OVrCSU[k][2],map_OVrCSU[k][3])
        
        dist_vec = indBB[order[chainid_resid_1][0] ,2:5] - indBB[order[chainid_resid_2][0] ,2:5]
        map_OVrCSU[k][4] = np.linalg.norm(dist_vec) /10     # [Ang] to [nm]

    pairs = []
    for k in range(0, len(map_OVrCSU)):
        if (map_OVrCSU[k][4] > cutoff_short) and (map_OVrCSU[k][4] < cutoff_long) and \
           ((map_OVrCSU[k][0] != map_OVrCSU[k][2]) or (abs(map_OVrCSU[k][1]-map_OVrCSU[k][3]) >= seqDist)):
            # parameters for LJ potential
            sigma = map_OVrCSU[k][4] / 1.12246204830 # calc sigma for the LJ potential in [nm]  
            Vii = 4.0 * pow(sigma,6) * go_eps
            Wii = 4.0 * pow(sigma,12) * go_eps

            chainid_resid_1=(map_OVrCSU[k][0],map_OVrCSU[k][1])
            chainid_resid_2=(map_OVrCSU[k][2],map_OVrCSU[k][3])
            pairs.append([indBB[order[chainid_resid_1][0], 0], indBB[order[chainid_resid_1][0] ,1],
                          indBB[order[chainid_resid_2][0], 0], indBB[order[chainid_resid_2][0] ,1],
                          Vii, Wii, map_OVrCSU[k][4], sigma])
            #paris: newresid 1, newatomid 1, newresid 2, newatomid 2, Vii, Wii, contact distance, sigma
        
        elif map_OVrCSU[k][4] > cutoff_long:
            print('This contanct is excluded due to distance > cutoff_long: ' + str(map_OVrCSU[k]))
        elif map_OVrCSU[k][4] < cutoff_short:
            print('This contanct is excluded due to distance < cutoff_short: ' + str(map_OVrCSU[k]))
        elif (abs(map_OVrCSU[k][1]-map_OVrCSU[k][3]) < seqDist) and (map_OVrCSU[k][0]==map_OVrCSU[k][1]) :
            print('This contanct is excluded because the AA have less than ' + str(seqDist) + ' other AA between each other: ' + str(map_OVrCSU[k]))

    sym_pairs = []
    # count contacts only once; exclude asymmetric rCSU contacts (cf. doi 10.1063/1.4929599)
    for k in range(0, len(pairs)):
        if pairs[k][1] < pairs[k][3]:
            for l in range(k+1, len(pairs)):
                if (pairs[l][1] == pairs[k][3]) and (pairs[l][3] == pairs[k][1]):
                    sym_pairs.append(pairs[k])

    print ('- - - -')
    print ('These results exclude the contacts with distances higher than the cutoff_long (' + str(cutoff_long) + ' nm), shorter than the cutoff_short (' + 
            str(cutoff_short) + ' nm), or where the AA have less than ' + str(seqDist) + ' other AA between each other:')
    print ('Sum of symmetric (doubly counted) and asymmetric OV + rCSU contacts: ' + str(len(pairs)))
    print ('Only symmetric OV + rCSU contacts (singly counted):' + str(len(sym_pairs)))

    return sym_pairs 

def write_files(file_pref, sym_pairs, order, indBB, Natoms, nameAA, go_eps, c6c12):
    # write the interaction table for the Go-like bonds
    with open(file_pref + '_go-table_VirtGoSites.itp','w') as f:
        f.write('; OV + symmetric rCSU contacts \n')
        if (c6c12 == 1):
            for k in range(0, len(sym_pairs)):
                # to write the LJ potential itp:
                s2print = " %s_%s  %s_%s    1  %.10f  %.10f  ;  %s  %s  %.3f \n" % (file_pref, str(int(sym_pairs[k][0])), file_pref, str(int(sym_pairs[k][2])), 
                                                                          sym_pairs[k][4], sym_pairs[k][5], str(int(sym_pairs[k][0])), 
                                                                          str(int(sym_pairs[k][2])), sym_pairs[k][6]) 
                                                            
                f.write(s2print)
        else:
            for k in range(0, len(sym_pairs)):
                # to write the LJ potential itp:
                s2print = " %s_%s  %s_%s    1  %.10f  %.10f  ;  %s  %s  %.3f \n" % (file_pref, str(int(sym_pairs[k][0])), file_pref, str(int(sym_pairs[k][2])), 
                                                                          sym_pairs[k][7], go_eps, str(int(sym_pairs[k][0])), 
                                                                          str(int(sym_pairs[k][2])), sym_pairs[k][6]) 
                                                            
                f.write(s2print)

    subprocess.call("echo '#include \"" + file_pref + "_go-table_VirtGoSites.itp\"' > go-table_VirtGoSites.itp ", shell=True)

    # write supplementary file: BB virtual particle definitions for martini.itp
    with open(file_pref + '_BB-part-def_VirtGoSites.itp','w') as f:
        f.write('; protein BB virtual particles \n')
        for k in range(0, len(indBB)):
            s2print = "%s_%s 0.0 0.000 A 0.0 0.0 \n" % (file_pref, int(indBB[k][0])) 
            f.write(s2print)
    subprocess.call("echo '#include \"" + file_pref + "_BB-part-def_VirtGoSites.itp\"' > BB-part-def_VirtGoSites.itp ", shell=True)

    # write supplementary file: exclusions for protein.itp
    with open(file_pref + '_exclusions_VirtGoSites.itp','w') as f:
        f.write(';[ exclusions ] \n')
        f.write('; OV + symmetric rCSU contacts \n')
        for k in range(0, len(sym_pairs)):
            s2print = " %s  %s  \t ;  %s  %s \n" % (str(int(sym_pairs[k][1])), str(int(sym_pairs[k][3])), 
                                                    str(int(sym_pairs[k][0])), str(int(sym_pairs[k][2]))) 
                                        # atom index and residue index adapted due to missing residues
            f.write(s2print)

                
    # write supplementary file: Go-like bonds as harmonic bonds for visulization of the protein
    newresid2index={item[1]:item[0] for item in order.values()}
    index2newresid={item[0]:item[1] for item in order.values()}
    
    with open(file_pref + '_go4view_harm.itp','w') as f:
        f.write('; Go bonds as harmonic bonds between the virtual particles: \n')
        f.write('; OV + symmetric rCSU contacts \n')
        for k in range(0, len(sym_pairs)):
            # to write the harmonic bonds itp:
            s2print = " %s  %s  1  %s  1250  ; %s_%s  %s_%s \n" % (str(int(newresid2index[sym_pairs[k][0]]+1+Natoms)), str(int(newresid2index[sym_pairs[k][2]]+1+Natoms)), str(round(sym_pairs[k][6],3)), 
                                                                 file_pref, str(int(sym_pairs[k][0])), file_pref, str(int(sym_pairs[k][2])))
                                                               # the bonds are added between the virtual particles
            f.write(s2print)
        for k in range(1, len(indBB)):
            if (np.sum(np.array(sym_pairs)[:,0]==index2newresid[k]) + np.sum(np.array(sym_pairs)[:,2]==index2newresid[k])) == 0:
                s2print = " %s  %s  1  1.  1     ; %s_%s  %s_%s --> added for vmd \n" % (str(int(k+1+Natoms)), str(int(k+Natoms)), file_pref, int(index2newresid[k]), file_pref, int(index2newresid[k-1]))
                f.write(s2print)

                

def main():
    args = user_input()

    file_CA, file_BB, file_OV, file_rCSU, seqDist, cols, c6c12 = get_settings()

    indBB, order, nameAA, map_OVrCSU = read_data(args.s, args.r, args.f, file_CA, file_BB, file_OV, file_rCSU, cols)

    sym_pairs = get_go(indBB, order, nameAA, map_OVrCSU, args.cutoff_short, args.cutoff_long, args.go_eps, seqDist)

    write_files(args.moltype, sym_pairs, order, indBB, args.Natoms, nameAA, args.go_eps, c6c12)
    print('All symmetric OV and rCSU contacts written! Have fun!')

if __name__ == '__main__':
    main()
