import os
from functools import partial
import ctgomartini
import subprocess

def RunScript(script_string):
    output = subprocess.run(script_string,
    shell=True, capture_output=True, encoding='utf-8')
    # print(output.args)
    if output.returncode != 0:
        stdout = output.stdout
        stderr = output.stderr
        raise Exception(f'Error! {output.args} {stdout} {stderr}')

class TestMBMartini:
    """
    Test running the multiple baisn GoMartini
    """
    path = os.path.dirname(__file__)

    def test_GenMBRun_EXP(self):
        working_dir = os.path.join(self.path, "../data/MDRun/EXP")
        os.chdir(working_dir)

        os.system('rm -r test')
        os.system('cp -a template test')
        
        os.chdir(os.path.join(working_dir, 'test'))
        # Fetch run_ctgomartini.py
        os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'data/run_ctgomartini.py')} .")
        # Generate Itp
        RunScript("python run_ctgomartini.py -i npt.inp > npt.log")
        RunScript("python run_ctgomartini.py -i md.inp > md.log")
        RunScript("python run_ctgomartini.py -i md_cpi.inp >> md.log")
    
    def test_GenMBRun_HAM(self):
        working_dir = os.path.join(self.path, "../data/MDRun/HAM")
        os.chdir(working_dir)

        os.system('rm -r test')
        os.system('cp -a template test')
        
        os.chdir(os.path.join(working_dir, 'test'))
        # Fetch run_ctgomartini.py
        os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'data/run_ctgomartini.py')} .")
        # Generate Itp
        RunScript("python run_ctgomartini.py -i npt.inp > npt.log")
        RunScript("python run_ctgomartini.py -i md.inp > md.log")
        RunScript("python run_ctgomartini.py -i md_cpi.inp >> md.log")
