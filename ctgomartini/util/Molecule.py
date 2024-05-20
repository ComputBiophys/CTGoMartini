from .Bonded_interaction import BondedInteraction_types

class Molecule:
    """
    Molecule._topology: dict
        {"category_name":contents}
    Molecule.atoms: list(list)
        contents of atoms
    Molecule._atoms: Molecule_Category class
        details of atoms
    """
    def __init__(self, name=None):
        self.name = name
        self._topology = {}

    def readLine(self, line, currentCategory):
        """read and classify every line into the suitable category"""
        line = line.strip().split()
        if currentCategory not in self._topology:
            self._topology[currentCategory] = []
        self._topology[currentCategory].append(line)

    def initialize(self, ff_parameters):
        for category_name, contents in self._topology.items():
            if category_name in Molecule_Categories:
                category = Molecule_Categories[category_name](contents)            
            elif category_name in Bonded_Categories:
                category = Bonded_Categories[category_name]
                category.contents = contents
            else:
                raise ValueError(f"Unknown category {category_name}")

            if category.parameter_missing:
                category.complement(ff_parameters)  
            self.__dict__[category_name] = category.contents
            self.__dict__[f'_{category_name}'] = category

class Molecule_Category:
    def __init__(self, name, description, category, contents):
        self.name = name
        self.description = description
        self.category = category
        self.contents = contents
        self.parameter_missing = False

    def __str__(self):
        return f"{self.name}: {self.description}"
    
    def initialize(self):
        pass

    def complement(self, ff_parameters):
        raise NotImplementedError(f"Do not implement complement for category ({self.name}) ")


Molecule_Categories = {}
# moleculetype
class Moleculetype(Molecule_Category):
    def __init__(self, contents=[]):
        super().__init__(name="moleculetype", 
                         description="defines the name of your molecule in this top and nrexcl = 3 stands for excluding non-bonded interactions between atoms that are no further than 3 bonds away.",
                         category="moleculetype",
                         contents=contents)

Molecule_Categories['moleculetype'] = Moleculetype

# atoms
class Atoms(Molecule_Category):
    def __init__(self, contents=[]):
        super().__init__(name="atoms", 
                         description="atom number; atom type; residue number; residue name; atom name; charge group number; (e); (u).",
                         category="atoms",
                         contents=contents)
        self.parameter_missing = True
        
    def complement(self, ff_parameters):
        """Complement the fields with the parameters of the force field"""
        atomtypes = ff_parameters['atomtypes']
        dict_atomtypes = {fields[0]: fields for fields in atomtypes}
        for fields in self.contents:
            if len(fields) < 5:
                raise ValueError(f"Too few fields in [ atoms ] line: {fields}")
            elif len(fields) > 8:
                raise ValueError(f"Too many fields in [ atoms ] line: {fields}")
            else:
                fields += [None] * (8 - len(fields))

            # charge missing
            if fields[6] is None:
                fields[6] = dict_atomtypes[fields[1]][4]
            # mass missing
            if fields[7] is None:
                fields[7] = dict_atomtypes[fields[1]][3]
            if None in fields:
                raise ValueError(f'None in the fields {fields}, {self.category}')

Molecule_Categories['atoms'] = Atoms

# multiple_basin
class Multiple_basin(Molecule_Category):
    def __init__(self, contents=[]):
        super().__init__(name="multiple_basin", 
                         description="bool, method, n_states, coupling_constant, energy1, energy2.",
                         category="multiple_basin",
                         contents=contents)

Molecule_Categories['multiple_basin'] = Multiple_basin

# Bonded Categories
Bonded_Categories = {}
for _Interaction in BondedInteraction_types:
    interaction = _Interaction()
    category_name = interaction.category

    if category_name not in Bonded_Categories:
        molecule_category = Molecule_Category(
                                    name=category_name,
                                    description=f"{category_name}: automated generation",
                                    category=category_name,
                                    contents=[])
        Bonded_Categories[category_name] = molecule_category
