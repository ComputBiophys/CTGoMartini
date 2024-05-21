import openmm as mm
import math
from .vsites import (
    LinearSite,
    OutOfPlane,
    VSiteManager,
    COMLinearSite,
    NormalizedInPlaneSite,
    NormalizedInPlaneTwoParticleSite,
)

class Interaction:
    def __init__(self, name, description, category, mm_force, type_label):
        self.name = name
        self.description = description

        self.category = category
        self.mm_force = mm_force
        self.type_label = type_label

        # default parameters
        self.contents = []
        self.intermolecule_sharing = True

    def __str__(self):
        return f"{self.name}: {self.description}"

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        """Add one interaction to mm_force"""
        raise NotImplementedError
    
    def get_exception(self, atoms, fields, base_atom_index=0, offset=-1):
        """Get the exception from the fields"""
        return []

# Register interaction types
BondedInteraction_types = []


# Bonds
class Harmonic_bonds(Interaction):
    def __init__(self):
        super().__init__(name='harmonic_bonds',
                         description='Harmonic bond potential 1/2*k*(r-r0)^2: atomid1, atomid2, functype, length, and k',
                         category='bonds',
                         mm_force=mm.HarmonicBondForce(),
                         type_label=[2, "1"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 5 and fields[self.type_label[0]] == self.type_label[1], "harmonic_bonds requires 5 items and the functype is 1"
        self.mm_force.addBond(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            float(fields[3]), float(fields[4]))

    def get_exception(self, atoms, fields, base_atom_index=0, offset=-1):
        return [[base_atom_index + int(fields[0]) + offset, base_atom_index + int(fields[1]) + offset, 0, 0, 0]]

BondedInteraction_types.append(Harmonic_bonds)


# Angles
class Harmonic_angles(Interaction):
    def __init__(self):
        super().__init__(name='harmonic_angles',
                         description='Harmonic angle potential 1/2*k*(theta-theta0)^2: atomid1, atomid2, atomid3, functype, angle, and k',
                         category='angles',
                         mm_force=mm.HarmonicAngleForce(),
                         type_label=[3, "1"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 6 and fields[self.type_label[0]] == self.type_label[1], "harmonic_angles requires 6 items and the functype is 1"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            base_atom_index + int(fields[2]) + offset, 
            float(fields[4])*degToRad, float(fields[5]))

BondedInteraction_types.append(Harmonic_angles)


class G96_angles(Interaction):
    def __init__(self):
        super().__init__(name='g96_angles',
                         description='g96 angles 0.5*k*(cos(theta)-cos(theta0))^2: atomid1, atomid2, atomid3, functype, angle, and k',
                         category='angles',
                         mm_force=mm.CustomAngleForce(
                             "0.5 * k * (cos(theta) - cos(theta0))^2"),
                         type_label=[3, "2"])
        self.mm_force.addPerAngleParameter('theta0')
        self.mm_force.addPerAngleParameter('k')

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 6 and fields[self.type_label[0]] == self.type_label[1], "g96_angles requires 6 items and the functype is 2"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            base_atom_index + int(fields[2]) + offset, 
            [float(fields[4])*degToRad, float(fields[5])])

BondedInteraction_types.append(G96_angles)


class Restricted_angles(Interaction):
    def __init__(self):
        super().__init__(name='restricted_angles',
                         description='Restricted angles 0.5*k*(cos(theta)-cos(theta0))^2 /sin(theta)^2: atomid1, atomid2, atomid3, functype, angle, and k',
                         category='angles',
                         mm_force=mm.CustomAngleForce(
                             "0.5 * k * (cos(theta) - cos(theta0))^2 / sin(theta)^2"),
                         type_label=[3, "10"])
        self.mm_force.addPerAngleParameter('theta0')
        self.mm_force.addPerAngleParameter('k')

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 6 and fields[self.type_label[0]] == self.type_label[1], "restricted_angles requires 6 items and the functype is 10"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            base_atom_index + int(fields[2]) + offset, 
            [float(fields[4])*degToRad, float(fields[5])])

BondedInteraction_types.append(Restricted_angles)


# Dihedrals
class Periodic_dihedrals(Interaction):
    def __init__(self):
        super().__init__(name='periodic_dihedrals',
                         description='Periodic dihedrals k*(1+cos(n*theta-theta0)): atomid1, atomid2, atomid3, atomid4, functype, theta, k, and phase',
                         category='dihedrals',
                         mm_force=mm.PeriodicTorsionForce(),
                         type_label=[4, "1", "4", "9"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 8 and fields[self.type_label[0]] in self.type_label[1:], "periodic_dihedrals requires 8 items and the functype is 1 or 4 or 9"
        degToRad = math.pi / 180
        self.mm_force.addTorsion(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            base_atom_index + int(fields[3]) + offset,
            float(fields[7]), float(fields[5])*degToRad, float(fields[6]))

BondedInteraction_types.append(Periodic_dihedrals)

class Harmonic_dihedrals(Interaction):
    def __init__(self):
        super().__init__(name='harmonic_dihedrals',
                         description='Harmonic dihedrals k*(1+cos(n*theta-theta0)): atomid1, atomid2, atomid3, atomid4, functype, theta, and k',
                         category='dihedrals',
                         mm_force=mm.CustomTorsionForce(
                                     "0.5*k*(thetap-theta0)^2; thetap = step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi); pi = %.15g"
                                     % math.pi
                         ),
                         type_label=[4, "2"])
        self.mm_force.addPerTorsionParameter("theta0")
        self.mm_force.addPerTorsionParameter("k")

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 7 and fields[self.type_label[0]] == self.type_label[1], "harmonic_dihedrals requires 7 items and the functype is 2"
        degToRad = math.pi / 180
        self.mm_force.addTorsion(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            base_atom_index + int(fields[3]) + offset,
            [float(fields[5])*degToRad, float(fields[6])])

BondedInteraction_types.append(Harmonic_dihedrals)

class Combined_bending_torsion_potentials(Interaction):
    def __init__(self):
        super().__init__(name='combined_bending_torsion_potentials',
                         description='Combined bending-torsion potentials: atomid1, atomid2, atomid3, atomid4, functype, k, a0, a1, a2, a3, and a4',
                         category='dihedrals',
                         mm_force=None,
                         type_label=[4, "11"])
        self.mm_force = mm.CustomCompoundBondForce(
                            4,
                            "k*sintheta0^3*sintheta1^3*(a0 + a1*cosphi + a2*cosphi^2 + a3*cosphi^3 + a4*cosphi^4); "
                            "sintheta0 = sin(angle(p1, p2, p3));"
                            "sintheta1 = sin(angle(p2, p3, p4));"
                            "cosphi = cos(dihedral(p1, p2, p3, p4));",
                        )
        self.mm_force.addPerBondParameter("k")
        self.mm_force.addPerBondParameter("a0")
        self.mm_force.addPerBondParameter("a1")
        self.mm_force.addPerBondParameter("a2")
        self.mm_force.addPerBondParameter("a3")
        self.mm_force.addPerBondParameter("a4")

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 11 and fields[self.type_label[0]] == self.type_label[1], "combined_bending_torsion_potentials requires 11 items and the functype is 11"
        k = float(fields[5])
        a0 = float(fields[6])
        a1 = float(fields[7])
        a2 = float(fields[8])
        a3 = float(fields[9])
        a4 = float(fields[10])
        self.mm_force.addBond(
            [
                base_atom_index + int(fields[0]) + offset,
                base_atom_index + int(fields[1]) + offset,
                base_atom_index + int(fields[2]) + offset,
                base_atom_index + int(fields[3]) + offset
            ],
            [k, a0, a1, a2, a3, a4],)

BondedInteraction_types.append(Combined_bending_torsion_potentials)

class Ryckaert_Bellemans_dihedrals(Interaction):
    def __init__(self):
        super().__init__(name='Ryckaert_Bellemans_dihedrals',
                         description='Ryckaert_Bellemans_dihedrals: atomid1, atomid2, atomid3, atomid4, functype, C0, C1, C2, C3, C4, and C5',
                         category='dihedrals',
                         mm_force=mm.RBTorsionForce(),
                         type_label=[4, "3", "5"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 11 and fields[self.type_label[0]] in self.type_label[1:], "Ryckaert_Bellemans_dihedrals requires 11 items and the functype is 3 or 5"
        c = [float(x) for x in fields[5:11]]
        if fields[self.type_label[0]] == "5":
            # Convert Fourier coefficients to RB coefficients.
            c = [
                c[1] + 0.5 * (c[0] + c[2]),
                0.5 * (-c[0] + 3 * c[2]),
                -c[1] + 4 * c[3],
                -2 * c[2],
                -4 * c[3],
                0,
            ]            
        self.mm_force.addTorsion(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            base_atom_index + int(fields[3]) + offset,
            c[0], c[1], c[2], c[3], c[4], c[5],)

BondedInteraction_types.append(Ryckaert_Bellemans_dihedrals)




# Constraints
class Constraints(Interaction):
    def __init__(self, sys=None):
        super().__init__(name='constraints',
                         description='Constraints (sys.addConstraint): atomid1, atomid2, functype, and length',
                         category='constraints',
                         mm_force=None,
                         type_label=[2, "1"])
        self.sys = sys

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        # Sometimes the fields in the constraints category have 5 items, which the 5th item perhaps has no meaning according to the Martini-openmm project.
        assert len(fields) >= 4 and fields[self.type_label[0]] == self.type_label[1], "constraints requires at least 4 items and the functype is 1"
        self.sys.addConstraint(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            float(fields[3]))

    def get_exception(self, atoms, fields, base_atom_index=0, offset=-1):
        return [[base_atom_index + int(fields[0]) + offset, base_atom_index + int(fields[1]) + offset, 0, 0, 0]]


BondedInteraction_types.append(Constraints)

# Pairs
class Pairs(Interaction):
    def __init__(self, epsilon_r=15, use_sigma_eps=True):
        super().__init__(name='pairs',
                         description='Pairs: atomid1, atomid2, functype, (C6, C12)',
                         category='pairs',
                         mm_force=None,
                         type_label=[2, "1"])
        self.epsilon_r = epsilon_r
        self.use_sigma_eps = use_sigma_eps

        self.mm_force = mm.CustomBondForce(
            "LJ + ES;"
            "LJ = C12/r^12 - C6/r^6;"
            "ES = f/epsilon_r*q1*q2 * 1/r;"
            f"epsilon_r = {self.epsilon_r};"
            "f = 138.935458;"
        )       
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C6")
        self.mm_force.addPerBondParameter("q1")
        self.mm_force.addPerBondParameter("q2")

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 5 and fields[self.type_label[0]] == self.type_label[1], "pairs requires 5 items and the functype is 1"
        
    def get_exception(self, atoms, fields, base_atom_index=0, offset=-1):
        if self.use_sigma_eps:
            sigma = float(fields[3])
            eps = float(fields[4])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[3])
            C12 = float(fields[4])
        q1 = float(atoms[int(fields[0])-1][6])
        q2 = float(atoms[int(fields[1])-1][6])
        self.mm_force.addBond(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            [C12, C6, q1, q2])
        return []

BondedInteraction_types.append(Pairs)

# Contacts
class Contacts(Interaction):
    def __init__(self, nonbonded_cutoff=1.1 * mm.unit.nanometer, use_sigma_eps=True):
        super().__init__(name='contacts',
                         description='Contacts Lenard-Jones Potential(r - cutoff): atomid1, atomid2, functype, C6/sigma, C12/epsilon',
                         category='contacts',
                         mm_force=None,
                         type_label=[2, "1"])
        self.nonbonded_cutoff = nonbonded_cutoff
        self.use_sigma_eps = use_sigma_eps

        self.mm_force = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
            f"rcut={self.nonbonded_cutoff.value_in_unit(mm.unit.nanometers)};"
        )
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C6")

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 5 and fields[self.type_label[0]] == self.type_label[1], "contacts requires 5 items and the functype is 1"
        if self.use_sigma_eps:
            sigma = float(fields[3])
            eps = float(fields[4])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[3])
            C12 = float(fields[4])

        self.mm_force.addBond(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            [C12, C6])

    def get_exception(self, atoms, fields, base_atom_index=0, offset=-1):
        q1 = float(atoms[int(fields[0])-1][6])
        q2 = float(atoms[int(fields[1])-1][6])
        return [[base_atom_index + int(fields[0]) + offset, base_atom_index + int(fields[1]) + offset, q1 * q2, 0, 0]]

BondedInteraction_types.append(Contacts)

# Cmaps

# # position_restraints
# class Position_restraints(Interaction):
#     def __init__(self, sys, name, description, category, mm_force, type_label):
#         super().__init__(name='position_restraints', 
#                          description='position_restraints are not currently supported',
#                          category='position_restraints',
#                          mm_force=None,
#                          type_label=None)



# virtual_sitesn
class VirtualSite(Interaction):
    def __init__(self, vsites, name, description, category, mm_force, type_label):
        super().__init__(name=name, 
                         description=description,
                         category=category,
                         mm_force=mm_force,
                         type_label=type_label)
        self.vsites = vsites


class Virtual_sitesn_COG(VirtualSite):
    def __init__(self, vsites=None):
        super().__init__(vsites=vsites,
                         name='virtual_sitesn_COG',
                         description="N-body virutal site (COG): atomid, functype, atomid1, atomid2, ..., atomidn",
                         category='virtual_sitesn',
                         mm_force=None,
                         type_label=[1, "1"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) >= 3 and fields[self.type_label[0]] == self.type_label[1], "virtual_sitesn_COG requires at least 3 items and the functype is 1"
        index =  int(fields[0]) + base_atom_index + offset
        from_atoms = [int(field) + base_atom_index + offset for field in fields[2:]]

        w = 1 / len(from_atoms)
        site_dict = {atom: w for atom in from_atoms}
        site = LinearSite(site_dict)
        self.vsites.add(index, site)

    def get_exception(self, atoms, fields, base_atom_index=0, offset=-1):
        if len(fields[2:]) == 1:
            return [[base_atom_index + int(fields[0]) + offset, base_atom_index + int(fields[2]) + offset, 0, 0, 0]]
        else:
            return []

BondedInteraction_types.append(Virtual_sitesn_COG)


class Virtual_sitesn_COM(VirtualSite):
    def __init__(self, vsites=None):
        super().__init__(vsites=vsites,
                         name='virtual_sitesn_COM',
                         description="N-body virutal site (COM): atomid, functype, atomid1, atomid2, ..., atomidn",
                         category='virtual_sitesn',
                         mm_force=None,
                         type_label=[1, "2"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) >= 3 and fields[self.type_label[0]] == self.type_label[1], "virtual_sitesn_COM requires at least 3 items and the functype is 2"
        index =  int(fields[0]) + base_atom_index + offset
        from_atoms = [int(field) + base_atom_index + offset for field in fields[2:]]        
        
        if len(from_atoms) == 1:
            site_dict = {from_atoms[0]: 1.0}
            site = LinearSite(site_dict)
        else:
            site = COMLinearSite(from_atoms)
            # site = site.to_linear(self.sys, offset)

        self.vsites.add(index, site)
        
    def get_exception(self, atoms, fields, base_atom_index=0, offset=-1):
        if len(fields[2:]) == 1:
            return [[base_atom_index + int(fields[0]) + offset, base_atom_index + int(fields[2]) + offset, 0, 0, 0]]
        else:
            return []
BondedInteraction_types.append(Virtual_sitesn_COM)


class Virtual_sites2(VirtualSite):
    def __init__(self, vsites=None):
        super().__init__(vsites=vsites,
                         name='virtual_sites2',
                         description="2-body virutal site: atomid, atomid1, atomid2, functype, weight",
                         category='virtual_sites2',
                         mm_force=None,
                         type_label=[3, "1"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 5 and fields[self.type_label[0]] == self.type_label[1], "virtual_sites2 requires 5 items and the functype is 1"
        index = int(fields[0]) + base_atom_index + offset
        atom1 = int(fields[1]) + base_atom_index + offset
        atom2 = int(fields[2]) + base_atom_index + offset
        w = float(fields[4])

        site_dict = {atom1: 1 - w, atom2: w}
        site = LinearSite(site_dict)
        self.vsites.add(index, site)

BondedInteraction_types.append(Virtual_sites2)

class Virtual_sites2fd(VirtualSite):
    def __init__(self, vsites=None):
        super().__init__(vsites=vsites,
                         name='virtual_sites2 (fd)',
                         description="2-body virutal site: atomid, atomid1, atomid2, functype, distance",
                         category='virtual_sites2',
                         mm_force=None,
                         type_label=[3, "2"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 5 and fields[self.type_label[0]] == self.type_label[1], "virtual_sites2 (fd) requires 5 items and the functype is 2"
        index = int(fields[0]) + base_atom_index + offset
        atom1 = int(fields[1]) + base_atom_index + offset
        atom2 = int(fields[2]) + base_atom_index + offset
        a = float(fields[4])

        site = NormalizedInPlaneTwoParticleSite(atom1, atom2, a)
        self.vsites.add(index, site)

BondedInteraction_types.append(Virtual_sites2fd)

class Virtual_sites3(VirtualSite):
    def __init__(self, vsites=None):
        super().__init__(vsites=vsites,
                         name='virtual_sites3',
                         description="3-body virutal site: atomid, atomid1, atomid2, atomid3, functype, a, b",
                         category='virtual_sites3',
                         mm_force=None,
                         type_label=[4, "1"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 7 and fields[self.type_label[0]] == self.type_label[1], "virtual_site3 requires 7 items and the functype is 1"
        index = int(fields[0]) + base_atom_index + offset
        atom1 = int(fields[1]) + base_atom_index + offset
        atom2 = int(fields[2]) + base_atom_index + offset
        atom3 = int(fields[3]) + base_atom_index + offset
        a = float(fields[5])
        b = float(fields[6])
        w1 = 1 - a - b
        w2 = a
        w3 = b
        site_dict = {
            atom1: w1,
            atom2: w2,
            atom3: w3,
        }
        site = LinearSite(site_dict)
        self.vsites.add(index, site)

BondedInteraction_types.append(Virtual_sites3)

class Virtual_sites3fd(VirtualSite):
    def __init__(self, vsites=None):
        super().__init__(vsites=vsites,
                         name='virtual_sites3 (fd)',
                         description="3-body virutal site: atomid, atomid1, atomid2, atomid3, functype, a, d",
                         category='virtual_sites3',
                         mm_force=None,
                         type_label=[4, "2"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 7 and fields[self.type_label[0]] == self.type_label[1], "virtual_site3 (fd) requires 7 items and the functype is 2"
        index = int(fields[0]) + base_atom_index + offset
        atom1 = int(fields[1]) + base_atom_index + offset
        atom2 = int(fields[2]) + base_atom_index + offset
        atom3 = int(fields[3]) + base_atom_index + offset
        a = float(fields[5])
        d = float(fields[6])
        site = NormalizedInPlaneSite(atom1, atom2, atom3, a, d)
        self.vsites.add(index, site)

BondedInteraction_types.append(Virtual_sites3fd)

class Virtual_sites3out(VirtualSite):
    def __init__(self, vsites=None):
        super().__init__(vsites=vsites,
                         name='virtual_sites3 (out)',
                         description="3-body virutal site: atomid, atomid1, atomid2, atomid3, functype, a, b, c",
                         category='virtual_sites3',
                         mm_force=None,
                         type_label=[4, "4"])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 8 and fields[self.type_label[0]] == self.type_label[1], "virtual_site3 (out) requires 8 items and the functype is 4"
        index = int(fields[0]) + base_atom_index + offset
        atom1 = int(fields[1]) + base_atom_index + offset
        atom2 = int(fields[2]) + base_atom_index + offset
        atom3 = int(fields[3]) + base_atom_index + offset
        a = float(fields[5])
        b = float(fields[6])
        c = float(fields[7])
        site = OutOfPlane(atom1, atom2, atom3, a, b, c)
        self.vsites.add(index, site)

BondedInteraction_types.append(Virtual_sites3out)



# Exclusions
class Exclusions(Interaction):
    def __init__(self):
        super().__init__(name='exclusions',
                         description='exclusions: atomid1, atomid2, ..., atomidn',
                         category='exclusions',
                         mm_force=None,
                         type_label=None)
                        #  fieldtypes=[int, '+'])

    def add_interaction(self, fields, base_atom_index=0, offset=-1):
        assert len(fields) >=2, "exclusions requires at least 2 items"
        pass

    def get_exception(self, atoms, fields, base_atom_index=0, offset=-1):
        exclusions = []
        index1 = base_atom_index + int(fields[0]) + offset
        from_atoms = fields[1:]
        for atom in from_atoms:
            index2 = base_atom_index + int(atom) + offset
            exclusions.append([index1, index2, 0, 0, 0])
        return exclusions

BondedInteraction_types.append(Exclusions)

# Multi_angles
class Multi_harmonic_angles(Interaction):
    def __init__(self):
        super().__init__(name='multi_harmonic_angles',
                         description='Multi_harmonic angle potential 1/2*k*(theta-theta0)^2: atomid1, atomid2, atomid3, total number of states, stateid, functype, angle, and k',
                         category='multi_angles',
                         mm_force=mm.HarmonicAngleForce(),
                         type_label=[5, "1"])
        self.intermolecule_sharing = False
        
    def add_interaction(self, state, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 8 and fields[self.type_label[0]] == self.type_label[1] and fields[4] == state, f"multi_harmonic_angles requires 8 items and the functype is 1 and state is {state}"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            base_atom_index + int(fields[2]) + offset, 
            float(fields[6])*degToRad, float(fields[7]))

BondedInteraction_types.append(Multi_harmonic_angles)


class Multi_g96_angles(Interaction):
    def __init__(self):
        super().__init__(name='multi_g96_angles',
                         description='multi_g96 angles 0.5*k*(cos(theta)-cos(theta0))^2: atomid1, atomid2, atomid3, total number of states, stateid, functype, angle, and k',
                         category='multi_angles',
                         mm_force=mm.CustomAngleForce(
                             "0.5 * k * (cos(theta) - cos(theta0))^2"),
                         type_label=[5, "2"])
        self.mm_force.addPerAngleParameter('theta0')
        self.mm_force.addPerAngleParameter('k')
        self.intermolecule_sharing = False

    def add_interaction(self, state, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 8 and fields[self.type_label[0]] == self.type_label[1] and fields[4] == state, f"g96_angles requires 8 items and the functype is 2 and state is {state}"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            base_atom_index + int(fields[2]) + offset, 
            [float(fields[6])*degToRad, float(fields[7])])

BondedInteraction_types.append(Multi_g96_angles)


class Multi_restricted_angles(Interaction):
    def __init__(self):
        super().__init__(name='multi_restricted_angles',
                         description='Multi_restricted angles 0.5*k*(cos(theta)-cos(theta0))^2 /sin(theta)^2: atomid1, atomid2, atomid3, total number of states, stateid, functype, angle, and k',
                         category='multi_angles',
                         mm_force=mm.CustomAngleForce(
                             "0.5 * k * (cos(theta) - cos(theta0))^2 / sin(theta)^2"),
                         type_label=[5, "10"])
        self.mm_force.addPerAngleParameter('theta0')
        self.mm_force.addPerAngleParameter('k')
        self.intermolecule_sharing = False

    def add_interaction(self, state, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 8 and fields[self.type_label[0]] == self.type_label[1] and fields[4] == state, f"restricted_angles requires 8 items and the functype is 10 and state is {state}"
        degToRad = math.pi / 180
        self.mm_force.addAngle(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            base_atom_index + int(fields[2]) + offset, 
            [float(fields[6])*degToRad, float(fields[7])])

BondedInteraction_types.append(Multi_restricted_angles)


# Multi_dihedrals
class Multi_periodic_dihedrals(Interaction):
    def __init__(self):
        super().__init__(name='multi_periodic_dihedrals',
                         description='Multi_periodic dihedrals k*(1+cos(n*theta-theta0)): atomid1, atomid2, atomid3, atomid4, total number of states, stateid, functype, theta, k, and phase',
                         category='multi_dihedrals',
                         mm_force=mm.PeriodicTorsionForce(),
                         type_label=[6, "1"])
        self.intermolecule_sharing = False

    def add_interaction(self, state, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 10 and fields[self.type_label[0]] == self.type_label[1] and fields[5] == state, f"periodic_dihedrals requires 10 items and the functype is 1 and state is {state}"
        degToRad = math.pi / 180
        self.mm_force.addTorsion(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            base_atom_index + int(fields[3]) + offset,
            float(fields[9]), float(fields[7])*degToRad, float(fields[8]))
        
BondedInteraction_types.append(Multi_periodic_dihedrals)


class Multi_harmonic_dihedrals(Interaction):
    def __init__(self):
        super().__init__(name='multi_harmonic_dihedrals',
                         description='Multi_harmonic dihedrals k*(1+cos(n*theta-theta0)): atomid1, atomid2, atomid3, atomid4, total number of states, stateid, functype, theta, and k',
                         category='multi_dihedrals',
                         mm_force=mm.CustomTorsionForce(
                                     "0.5*k*(thetap-theta0)^2; thetap = step(-(theta-theta0+pi))*2*pi+theta+step(theta-theta0-pi)*(-2*pi); pi = %.15g"
                                     % math.pi
                         ),
                         type_label=[6, "2"])
        self.mm_force.addPerTorsionParameter("theta0")
        self.mm_force.addPerTorsionParameter("k")
        self.intermolecule_sharing = False

    def add_interaction(self, state, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 9 and fields[self.type_label[0]] == self.type_label[1] and fields[5] == state, f"harmonic_dihedrals requires 9 items and the functype is 2 and state is {state}"
        degToRad = math.pi / 180
        self.mm_force.addTorsion(
            base_atom_index + int(fields[0]) + offset,
            base_atom_index + int(fields[1]) + offset,
            base_atom_index + int(fields[2]) + offset,
            base_atom_index + int(fields[3]) + offset,
            [float(fields[7])*degToRad, float(fields[8])])

BondedInteraction_types.append(Multi_harmonic_dihedrals)

# Multi_contacts
class Multi_contacts(Interaction):
    def __init__(self, nonbonded_cutoff=1.1 * mm.unit.nanometer, use_sigma_eps=True):
        super().__init__(name='multi_contacts',
                         description='Multi_contacts Lenard-Jones Potential(r - cutoff): atomid1, atomid2, total number of states, stateid, functype, C6/sigma, C12/epsilon',
                         category='multi_contacts',
                         mm_force=None,
                         type_label=[4, "1"])
        self.nonbonded_cutoff = nonbonded_cutoff
        self.use_sigma_eps = use_sigma_eps

        self.mm_force = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
            f"rcut={self.nonbonded_cutoff.value_in_unit(mm.unit.nanometers)};"
        )
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C6")

        self.intermolecule_sharing = False

    def add_interaction(self, state, fields, base_atom_index=0, offset=-1):
        assert len(fields) == 7 and fields[self.type_label[0]] == self.type_label[1] and fields[3] == state, f"contacts requires 7 items and the functype is 1 and state is {state}"
        if self.use_sigma_eps:
            sigma = float(fields[5])
            eps = float(fields[6])
            C6 = 4 * eps * sigma ** 6
            C12 = 4 * eps * sigma ** 12
        else:
            C6 = float(fields[5])
            C12 = float(fields[6])

        self.mm_force.addBond(
            base_atom_index + int(fields[0]) + offset, 
            base_atom_index + int(fields[1]) + offset, 
            [C12, C6])

    def get_exception(self, atoms, fields, base_atom_index=0, offset=-1):
        q1 = float(atoms[int(fields[0])-1][6])
        q2 = float(atoms[int(fields[1])-1][6])
        return [[base_atom_index + int(fields[0]) + offset, base_atom_index + int(fields[1]) + offset, q1 * q2, 0, 0]]

BondedInteraction_types.append(Multi_contacts)




# Summary
NonLocal_BondedInteraction_dict = {}
Local_BondedInteraction_dict = {}
for _Interaction in BondedInteraction_types:
    interaction = _Interaction()
    category = interaction.category
    if interaction.intermolecule_sharing:
        if category not in NonLocal_BondedInteraction_dict:
            NonLocal_BondedInteraction_dict[category] = [_Interaction]
        else:
            NonLocal_BondedInteraction_dict[category].append(_Interaction)
    else:
        if category not in Local_BondedInteraction_dict:
            Local_BondedInteraction_dict[category] = [_Interaction]
        else:
            Local_BondedInteraction_dict[category].append(_Interaction)


# Multiple basin
class EXP_Interaction():
    def __init__(self):
        self.name="exponential mixing scheme",
        self.description='exponential mixing scheme for multiple baisn popential',
        
        
    def addForce(self, mbp_force_dict, coupling_constant, basin_energy_list):
        beta = coupling_constant
        # print("beta", beta)
        basin_energy_list = basin_energy_list

        part1_list = []
        part2_list = []
        for state, force_set in mbp_force_dict.items():
            part1_list.append(f"exp(-beta * (energy{state} + C{state}))")
            energy_combined = ' + '.join([f'state{state}_force{j+1}' for j in range(len(force_set))])
            part2_list.append(f"energy{state} = {energy_combined};")
        
        part1 = ' + '.join(part1_list)
        part2 = '\n'.join(part2_list)
        energy= f"""-1/beta * log({part1});\n{part2}"""
        print(energy)

        self.mm_force = mm.CustomCVForce(energy)
        for state, force_set in mbp_force_dict.items():
            for j, force in enumerate(force_set):
                self.mm_force.addCollectiveVariable(f"state{state}_force{j+1}", force)
        self.mm_force.addGlobalParameter("beta", float(beta))
        for i, basin_energy in enumerate(basin_energy_list):
            self.mm_force.addGlobalParameter(f'C{i+1}', float(basin_energy))
        return self.mm_force


class HAM_Interaction():
    def __init__(self):
        self.name="hamiltonian mixing scheme",
        self.description='hamiltonian mixing scheme for multiple baisn popential',
        
    def addForce(self, mbp_force_dict, coupling_constant, basin_energy_list):
        assert len(mbp_force_dict.keys()) == 2
        delta = coupling_constant
        basin_energy_list = basin_energy_list

        part2_list = []
        for state, force_set in mbp_force_dict.items():
            energy_combined = ' + '.join([f'state{state}_force{j+1}' for j in range(len(force_set))])
            part2_list.append(f"energy{state} = {energy_combined};")
        
        part1 = '(energy1+energy2+deltaV)/2 - sqrt(((energy1-energy2-deltaV)/2)^2+delta^2);deltaV=mbp_energy2-mbp_energy1;'
        part2 = '\n'.join(part2_list)
        energy= f"""{part1};\n{part2}"""
        print(energy)

        self.mm_force = mm.CustomCVForce(energy)
        for state, force_set in mbp_force_dict.items():
            for j, force in enumerate(force_set):
                self.mm_force.addCollectiveVariable(f"state{state}_force{j+1}", force)
        self.mm_force.addGlobalParameter("delta", float(delta))
        for i, basin_energy in enumerate(basin_energy_list):
            self.mm_force.addGlobalParameter(f'mbp_energy{i+1}', float(basin_energy))
        return self.mm_force
    