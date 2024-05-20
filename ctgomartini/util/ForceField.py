class ForceField:
    """
    ForceField._parameters: dict
        {"category_name":contents}
    ForceField.atomtypes: list(list)
        contents of atomtypes
    ForceField._atomtypes: FF_Category class
        details of atomtypes
    """
    def __init__(self, name=None):
        self.name = name
        self._parameters = {}

    def readLine(self, line, currentCategory):
        line = line.strip().split()
        if currentCategory not in self._parameters:
            self._parameters[currentCategory] = []
        self._parameters[currentCategory].append(line)

    def initialize(self):
        for category_name, contents in self._parameters.items():
            if category_name in FF_Categories:
                category = FF_Categories[category_name](contents)
            else:
                try:
                    category = FF_Category(name = category_name,
                                           description="Automated generation",
                                           category=category_name,
                                           contents=contents)
                except:
                    raise ValueError(
                        f"[ {category_name} ] is not currently supported.\n")

            category.initialize()
            self.__dict__[category_name] = category.dict_contents
            self.__dict__[f'_{category_name}'] = category


class FF_Category:
    def __init__(self, name, description, category, contents):
        self.name = name
        self.description = description
        self.category = category
        self.contents = contents
        self.dict_contents = {}

    def __str__(self):
        return f"{self.name}: {self.description}"
    
    def initialize(self):
        pass

FF_Categories = {}

# defaults
class Defaults(FF_Category):
    def __init__(self, contents=[]):
        super().__init__(name="defaults", 
                         description="defaults for Gromacs. More information can be seen in the mannual https://manual.gromacs.",
                         category="defaults",
                         contents=contents)
        self.dict_contents = self.contents # force the dict_contents is contents
FF_Categories["defaults"]=Defaults

# atomtypes
class Atomtypes(FF_Category):
    def __init__(self, contents=[]):
        super().__init__(name="atomtypes", 
                         description="atomtypes for Gromacs (atom type; (bonded type; atomic number;) m (u); q (e); particle type; V ; W (bonded type and atomic number are optional)).",
                         category="atomtypes",
                         contents=contents)

    def initialize(self):
        """Initialize the fields for [ atomtypes ]"""
        for fields in self.contents:
            if len(fields) < 6:
                raise ValueError(f"Too few fields in [ atomtypes ] line: {fields}")
            if len(fields[3]) == 1:
                # Bonded type and atomic number are both missing.
                fields.insert(1, None)
                fields.insert(1, None)
            elif len(fields[4]) == 1 and fields[4].isalpha():
                if fields[1][0].isalpha():
                    # Atomic number is missing.
                    fields.insert(2, None)
                else:
                    # Bonded type is missing.
                    fields.insert(1, None)    

            if len(fields) != 8:
                raise ValueError(f"Something wrong in [ atomtypes ] lines: {fields}")    
            
            assert fields[1] is None, "Unsupport bonded types: {fields}"
            self.dict_contents[fields[0]] = fields

FF_Categories["atomtypes"]=Atomtypes

# nonbond_params
class Nonbond_params(FF_Category):
    def __init__(self, contents=[]):
        super().__init__(name="nonbond_params", 
                         description="nonbond_params for Gromacs (V ; W).",
                         category="nonbond_params",
                         contents=contents)
    
    def initialize(self):
        """Initialize the fields for [ nonbond_params ]"""
        for fields in self.contents:
            if len(fields) < 5:
                raise ValueError(f"Too few fields in [ nonbond_params ] line: {fields}")
            if fields[2] != "1":
                raise ValueError(
                    f"Unsupported function type in [ nonbond_params ] line: {fields}"
                )
            self.dict_contents[tuple(sorted(fields[:2]))] = fields

FF_Categories["nonbond_params"]=Nonbond_params
