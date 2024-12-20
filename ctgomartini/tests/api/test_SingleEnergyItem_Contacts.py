from ctgomartini.api import MartiniTopFile
from function import *

def Distance(p1,p2):
    return np.linalg.norm(p1-p2)

def Energy_LJ(r,sigma,epsilon, rcut=1.1):
    C12=4*epsilon*sigma**12
    C6=4*epsilon*sigma**6
    
    #step=1 if rcut >= r else 0
    step=np.where(rcut >= r, 1, 0)
    #energy=step*((C12/r**12-C6/r**6))
    energy=step*((C12/r**12-C6/r**6)-(C12/rcut**12-C6/rcut**6))
    return energy

def Cal_Dist_from_atoms(atom1, atom2):
    p1=atom1.position/10
    p2=atom2.position/10
    r=Distance(p1,p2)
    return r

def Cal_LJ_energy_from_atoms(atom1, atom2, sigma, epsilon, rcut=1.1):
    r=Cal_Dist_from_atoms(atom1, atom2)
    return Energy_LJ(r, sigma, epsilon, rcut=rcut)


def ContactEnergyComparison(working_dir, contact_str):
    os.chdir(os.path.join(working_dir, "openmm1"))
    strfile = "minimized.gro"
    topfile = "system.top"

    simulation = OMM_setSimulation(strfile, topfile, epsilon_r=15, temperature=310.15, double_precision=True)
    OMM_calStrfile(strfile, simulation, set_vsite=True)

    omm_energy1=Load_energy(clean=True)[0][1]

    os.chdir(os.path.join(working_dir, "openmm2"))
    strfile = "minimized.gro"
    topfile = "system.top"

    simulation = OMM_setSimulation(strfile, topfile, epsilon_r=15, temperature=310.15, double_precision=True)
    OMM_calStrfile(strfile, simulation, set_vsite=True)

    omm_energy2=Load_energy(clean=True)[0][1]

    #-----------------------------
    contact = contact_str.split()
    atomid1, atomid2 = int(contact[0]), int(contact[1])
    sigma, epsilon = float(contact[-2]), float(contact[-1])
    u = mda.Universe(strfile)
    atom1, atom2 = u.atoms[atomid1-1], u.atoms[atomid2-1]
    diff_energy = Cal_LJ_energy_from_atoms(atom1, atom2, sigma, epsilon, rcut=1.1)

    print("diff_energy:", omm_energy2 - omm_energy1)
    print("diff_energy_cal:", diff_energy)
    print('abs:', np.abs((omm_energy2 - omm_energy1)-diff_energy))
    assert np.abs((omm_energy2 - omm_energy1)-diff_energy) < 1e-5

class TestMartiniTopology:
    """
    Test the Category [ contacts ]
    """
    path = os.path.dirname(__file__)

    def test_GlnBP(self):
        working_dir = os.path.join(self.path, "../data/Contacts/GlnBP_go_m3_contacts")
        contact_str = "513 549 1 0.7097616382 12.0"
        ContactEnergyComparison(working_dir, contact_str)
    

    

    