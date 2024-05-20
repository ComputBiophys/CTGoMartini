import openmm as mm
import openmm.unit as unit
import math

class Interaction:
    def __init__(self, name, description, mm_force):
        self.name = name
        self.description = description
        self.mm_force = mm_force

        # default parameters
        self.contents = []

    def __str__(self):
        return f"{self.name}: {self.description}"

    # def add_interactions(self, base_atom_index=0, offset=-1):
    #     """Add all interactions to mm_force"""
    #     for fields in self.contents:
    #         self.add_interaction(fields, base_atom_index=base_atom_index, offset=offset)

    # def add_interaction(self, fields, base_atom_index=0, offset=-1):
    #     """Add one interaction to mm_force"""
    #     raise NotImplementedError

# # Lenard-Jones Interaction
# class LJInteraction(Interaction):
#     def __init__(self, nonbonded_cutoff = 1.1 * unit.nanometers):
#         super().__init__(name="LJ_interaction",
#                          description="Shifted Lenard-Jones potential",
#                          mm_force=None)
#         self.mm_force =  mm.CustomNonbondedForce(
#             "step(rcut-r)*(LJ - corr);"
#             "LJ = (C12(type1, type2) / r^12 - C6(type1, type2) / r^6);"
#             "corr = (C12(type1, type2) / rcut^12 - C6(type1, type2) / rcut^6);"
#             f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
#         )
#         self.mm_force.addPerParticleParameter("type")
#         self.mm_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
#         self.mm_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))          

# # Coulombic Interaction
# class CoulInteraction(Interaction):
#     def __init__(self, epsilon_r = 15, nonbonded_cutoff = 1.1 * unit.nanometers):
#         super().__init__(name="Coulomb_interaction",
#                          description="Coulombic interaction",
#                          mm_force=None)
#         self.mm_force = mm.CustomNonbondedForce(
#             "step(rcut-r)*ES;"
#             "ES = f/epsilon_r*q1*q2 * 1/r;"
#             f"epsilon_r = {epsilon_r};"
#             f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
#         )
#         self.mm_force.addPerParticleParameter("q")
#         self.mm_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
#         self.mm_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))        

# # Reaction-Field Interaction
# class RFInteraction(Interaction):
#     def __init__(self, epsilon_r = 15, nonbonded_cutoff = 1.1 * unit.nanometers):
#         super().__init__(name="Reaction_Field_interaction",
#                          description="Reaction Field interaction",
#                          mm_force=None)
#         self.mm_force = mm.CustomNonbondedForce(
#             "step(rcut-r)*ES;"
#             "ES = f/epsilon_r*q1*q2 * (krf * r^2 - crf);"
#             "crf = 1 / rcut + krf * rcut^2;"
#             "krf = 1 / (2 * rcut^3);"
#             f"epsilon_r = {epsilon_r};"
#             "f = 138.935458;"
#             f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
#         )
#         self.mm_force.addPerParticleParameter("q")
#         self.mm_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
#         self.mm_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))        

# nonbonded_interaction
class Nonbonded_interaction(Interaction):
    def __init__(self, epsilon_r = 15, nonbonded_cutoff = 1.1 * unit.nanometers):
        super().__init__(name="nonbonded_interaction",
                         description="Shifted Lenard-Jones potential and Reaction-field modified electrostatic interaction",
                         mm_force = None)
        self.mm_force = mm.CustomNonbondedForce(
            "step(rcut-r)*(LJ - corr + ES);"
            "LJ = (C12(type1, type2) / r^12 - C6(type1, type2) / r^6);"
            "corr = (C12(type1, type2) / rcut^12 - C6(type1, type2) / rcut^6);"
            "ES = f/epsilon_r*q1*q2 * (1/r + krf * r^2 - crf);"
            "crf = 1 / rcut + krf * rcut^2;"
            "krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {epsilon_r};"
            "f = 138.935458;"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.mm_force.addPerParticleParameter("type")
        self.mm_force.addPerParticleParameter("q")
        self.mm_force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
        self.mm_force.setCutoffDistance(nonbonded_cutoff.value_in_unit(unit.nanometer))        


# custom non-bonded force to add in the electrostatic terms for self and excluded interactions
class ES_self_excl_interaction(Interaction):
    def __init__(self, epsilon_r = 15, nonbonded_cutoff = 1.1 * unit.nanometers):
        super().__init__(name='es_self_excl_interaction',
                         description='custom non-bonded force to add in the electrostatic terms for self and excluded interactions',
                         mm_force = None)
        
        self.mm_force = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*qprod/epsilon_r * (krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {epsilon_r};"
            f"f = 138.935458;"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.mm_force.addPerBondParameter("qprod")

# custom non-bonded force for the ES exceptions
class ES_except_interaction(Interaction):
    def __init__(self, epsilon_r = 15, nonbonded_cutoff = 1.1 * unit.nanometers):
        super().__init__(name='es_except_interaction',
                         description='custom non-bonded force for the ES exceptions',
                         mm_force = None)

        self.mm_force = mm.CustomBondForce(
            f"step(rcut-r) * ES;"
            f"ES = f*qprod/epsilon_r * (1/r + krf * r^2 - crf);"
            f"crf = 1 / rcut + krf * rcut^2;"
            f"krf = 1 / (2 * rcut^3);"
            f"epsilon_r = {epsilon_r};"
            f"f = 138.935458;"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.mm_force.addPerBondParameter("qprod")

# custom bonded force to handle exceptions
class LJ_except_interaction(Interaction):
    def __init__(self, epsilon_r = 15, nonbonded_cutoff = 1.1 * unit.nanometers):
        super().__init__(name='lj_except_interaction',
                         description='custom bonded force to handle exceptions',
                         mm_force = None)

        self.mm_force = mm.CustomBondForce(
            "step(rcut-r) * (energy - corr);"
            "energy = (C12/r^12 - C6/r^6);"
            "corr = (C12/rcut^12 - C6/rcut^6);"
            f"rcut={nonbonded_cutoff.value_in_unit(unit.nanometers)};"
        )
        self.mm_force.addPerBondParameter("C12")
        self.mm_force.addPerBondParameter("C6")