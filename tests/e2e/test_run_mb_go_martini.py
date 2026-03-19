"""
End-to-end tests for molecular dynamics simulation runs.

Tests the complete run_ctgomartini.py workflow which executes:
    1. Energy minimization
    2. NPT equilibration
    3. Production MD run
    4. Checkpoint restart

Note:
    These tests may take several minutes to complete.
"""

import os
from functools import partial
import ctgomartini
import subprocess

from tests.conftest import WorkingDirectoryContext


def run_script(script_string):
    """Run a shell script and check for errors."""
    output = subprocess.run(script_string,
    shell=True, capture_output=True, encoding='utf-8')
    if output.returncode != 0:
        stdout = output.stdout
        stderr = output.stderr
        raise Exception(f'Error! {output.args} {stdout} {stderr}')

class TestRunMBGoMartini:
    """
    Test running the multiple basin GoMartini simulation.
    """
    path = os.path.dirname(__file__)

    def test_run_mb_exp(self):
        working_dir = os.path.join(self.path, "../fixtures/md_run/EXP")

        with WorkingDirectoryContext(working_dir):
            os.system('rm -r test 2>/dev/null; cp -a template test')

            with WorkingDirectoryContext(os.path.join(working_dir, 'test')):
                # Fetch run_ctgomartini.py
                os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'cli/run_ctgomartini.py')} .")
                # Generate Itp
                run_script("python run_ctgomartini.py -i npt.inp > npt.log")
                run_script("python run_ctgomartini.py -i md.inp > md.log")
                run_script("python run_ctgomartini.py -i md_cpi.inp >> md.log")

    def test_run_mb_ham(self):
        working_dir = os.path.join(self.path, "../fixtures/md_run/HAM")

        with WorkingDirectoryContext(working_dir):
            os.system('rm -r test 2>/dev/null; cp -a template test')

            with WorkingDirectoryContext(os.path.join(working_dir, 'test')):
                # Fetch run_ctgomartini.py
                os.system(f"cp {os.path.join(ctgomartini.__path__[0], 'cli/run_ctgomartini.py')} .")
                # Generate Itp
                run_script("python run_ctgomartini.py -i npt.inp > npt.log")
                run_script("python run_ctgomartini.py -i md.inp > md.log")
                run_script("python run_ctgomartini.py -i md_cpi.inp >> md.log")
