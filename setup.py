# setup.py
from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

# read version info
import re
VERSIONFILE="ctgomartini/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

# setup
setup(
    name='ctgomartini',
    version=verstr,
    description='CTGoMartini - A Python Library For Protein Conformational Transitions',
    long_description=long_description,
    # long_description_content_type='text/x-rst',
    # url='https://github.com/wlsong/PyLipID',
    author='Yang Song',
    author_email='yangsong2015@pku.edu.cn',
    classifiers=[
        # 'Development Status :: 5 - Production/Stable',
        'Development Status :: 0.3 - Production/Dev',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        # 'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.8',
    ],
    # keywords='simulation tools, network community, binding site',
    python_requires='>=3.8',
    packages=find_packages(),
    install_requires=[
        "MDAnalysis>=2.4",
        "numpy",
        "pandas",
        "matplotlib>=3.3.3",
        "vermouth==0.9.6",
        #"openmm>=8.1",
#        "dssp"
    ],
    # package_data={'ctgomartini': [
    # 'data/*',
    # 'data/ForceFields/Martini300/*',
    # 'data/Membrane/*',
    # 'data/Soluble/*',
    # # Add other specific file patterns as needed
    # ]},
    include_package_data=True,
)
