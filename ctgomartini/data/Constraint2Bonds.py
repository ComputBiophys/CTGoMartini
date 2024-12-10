from ctgomartini.api import MartiniTopFile
from ctgomartini.func import WriteItp



top = MartiniTopFile('system.top')
mol = top._moleculeTypes['TREK1']

for item in mol._topology['constraints']:
    item.append('10000')
    mol._topology['bonds'].append(item)
    
mol._topology['constraints']  = []
del mol._topology['constraints']

WriteItp(mol)