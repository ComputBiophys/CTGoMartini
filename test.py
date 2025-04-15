import re
Extract=lambda atomtype: re.findall(r'^(\w+)_(\d+)$',atomtype)[0][0]

print(Extract('Cadfa_123'))
