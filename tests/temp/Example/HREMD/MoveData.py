import os
import subprocess

dir=os.path.basename(os.getcwd())
subprocess.run(f'mkdir ../../Data/{dir}', shell=True)
subprocess.run(f'cp -a dRMSTraj_nc* ../../Data/{dir}', shell=True)
subprocess.run(f'cp -a output.nc ../../Data/{dir}', shell=True)
subprocess.run(f'cp -a output.nc ../../Data/{dir}', shell=True)

