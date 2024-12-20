
from simtk import unit as u
from simtk import openmm as mm
from simtk.openmm import app
from openmm.app import *
import os
import sys
import subprocess
import MDAnalysis as mda
import numpy as np

from ctgomartini.api import MartiniTopFile

def OMM_setSimulation(
        strfile,
        topfile,
        epsilon_r=15.0,
        temperature=310.15,
        double_precision=False):

    # Box_vectors Initiation
    if strfile.split('.')[-1] == 'gro':
        conf = GromacsGroFile(strfile)
        box_vectors = conf.getPeriodicBoxVectors()
    elif strfile.split('.')[-1] == 'pdb':
        conf = PDBFile(strfile)
        box_vectors = conf.getTopology().getPeriodicBoxVectors()
    else:
        raise ValueError(f"Cannot find {strfile}")

    # Set Platform
    if double_precision:
        platform = mm.Platform.getPlatformByName("Reference")
    else:
        platform = mm.Platform.getPlatformByName("CPU")

    # get any defines
    defines = {}
    try:
        with open("defines.txt") as def_file:
            for line in def_file:
                line = line.strip()
                defines[line] = True
    except FileNotFoundError:
        pass

    # Get system and top
    top = MartiniTopFile(
        topfile,
        periodicBoxVectors=box_vectors,
        defines=defines,
    )
    system = top.create_system(
        nonbonded_cutoff=1.1 * u.nanometer,
        epsilon_r=epsilon_r,)
    integrator = mm.LangevinIntegrator(
        temperature * u.kelvin, 1.0 / u.picosecond, 20 * u.femtosecond
    )

    simulation = mm.app.Simulation(top.topology, system, integrator, platform)
    simulation.__dict__['top'] = top
    return simulation


def OMM_calStrfile(strfile, simulation, prefix='conf', isEnergy=True, isForces=True, set_vsite=True):
    """
    Calculate and save the energy and forces from the simulation.

    Parameters:
    strfile (str): Path to the input .gro or .pdb file.
    simulation (simulation): OpenMM simulation object.
    prefix (str): Prefix for the output files, default is 'conf'.
    isEnergy (bool): Whether to output the energy file, default is True.
    isForces (bool): Whether to output the force file, default is True.
    set_vsite (bool): Whether to set the force of virtual sites (vsite) to zero, default is True.
    """

    if strfile.split('.')[-1] == 'gro':
        conf = GromacsGroFile(strfile)
        box_vectors = conf.getPeriodicBoxVectors()
    elif strfile.split('.')[-1] == 'pdb':
        conf = PDBFile(strfile)
        box_vectors = conf.getTopology().getPeriodicBoxVectors()

    simulation.context.setPeriodicBoxVectors(*box_vectors)
    simulation.context.setPositions(conf.getPositions())

    # Compute energies and forces
    state = simulation.context.getState(getEnergy=isEnergy, getForces=isForces)

    if isEnergy:
        energy = np.array(state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole), ndmin=2)
        energy_array = np.zeros(shape=(energy.shape[0], 2))
        energy_array[:, 0] = 0  # frame : 0
        energy_array[:, 1] = energy
        np.savetxt(f"energy_{prefix}.dat", energy_array)

    if isForces:
        forces = np.array(
            state.getForces().value_in_unit(u.kilojoule / u.nanometer / u.mole), ndmin=2
        )
        if set_vsite:
            for i, atom in enumerate(simulation.top.topology.atoms()):
                if simulation.system.isVirtualSite(i):
                    forces[i, :] = 0.0

        forces_array = np.zeros(shape=(forces.shape[0], 4))
        forces_array[:, 0] = 0  # frame: 0
        forces_array[:, 1:] = forces
        np.savetxt(f"forces_{prefix}.dat", forces_array)

# generate run_gmx.sh
def GMX_set(
    strfile='minimized.gro',
    trjfile='minimized.gro',
    topfile='system.top',
    mdpfile='martini_md.mdp',
    indexfile='system.ndx',
    reffile='default',
    prefix='conf',
    run_file='run_gmx.sh',
    CreateMDP=True,
    double_precision=False
):
    if reffile == 'default':
        reffile = strfile
    if double_precision:
        if indexfile is not None:
            line = f"""
# source double-precision gmx_d
source /usr/local/gromacs-2020.7-d/bin/GMXRC
set -e
# use double precision for everything
gmx_d grompp -f {mdpfile} -c {strfile} -p {topfile} -n {indexfile} -o {prefix}.tpr -r {reffile} -maxwarn 1 -v
gmx_d mdrun -deffnm {prefix} -rerun {trjfile} -nt 1 -v
echo pot | gmx_d energy -f {prefix} -o energy_{prefix}
echo 0 | gmx_d traj -f {prefix} -s {prefix} -of forces_{prefix}
rm mdout.mdp
rm {prefix}.edr {prefix}.log {prefix}.tpr {prefix}.trr
"""
        else:
            line = f"""
# source double-precision gmx_d
source /usr/local/gromacs-2020.7-d/bin/GMXRC
set -e
# use double precision for everything
gmx_d grompp -f {mdpfile} -c {strfile} -p {topfile} -o {prefix}.tpr -r {reffile} -maxwarn 1 -v
gmx_d mdrun -deffnm {prefix} -rerun {trjfile} -nt 1 -v
echo pot | gmx_d energy -f {prefix} -o energy_{prefix}
echo 0 | gmx_d traj -f {prefix} -s {prefix} -of forces_{prefix}
rm mdout.mdp
rm {prefix}.edr {prefix}.log {prefix}.tpr {prefix}.trr
"""
    else:
        line = f"""
# use single precision for everything
gmx grompp -f {mdpfile} -c {strfile} -p {topfile} -n {indexfile} -o {prefix}.tpr -r {reffile} -maxwarn 1 -v
gmx mdrun -deffnm {prefix} -rerun {trjfile} -nt 1 -v
echo pot | gmx energy -f {prefix} -o energy_{prefix}
echo 0 | gmx traj -f {prefix} -s {prefix} -of forces_{prefix}
rm mdout.mdp
rm {prefix}.edr {prefix}.log {prefix}.tpr {prefix}.trr
"""        

    with open(run_file, 'w') as fp:
        fp.write(line)

    line = """
integrator               = md
tinit                    = 0.0
dt                       = 0.020
nsteps                   = 1
nstxout                  = 1
nstvout                  = 1
nstfout                  = 1
nstlog                   = 1
nstenergy                = 1
nstxout-compressed       = 5000
compressed-x-precision   = 1000
compressed-x-grps        = System
energygrps               =

cutoff-scheme            = Verlet
nstlist                  = 20
ns_type                  = grid
pbc                      = xyz
verlet-buffer-tolerance  = 0.005

epsilon_r                = 15
coulombtype              = reaction-field
rcoulomb                 = 1.1
vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1

tcoupl                   = v-rescale
tc-grps                  = system
tau_t                    = 1.0
ref_t                    = 310.15

; Pressure coupling:
Pcoupl                   = Parrinello-rahman
Pcoupltype               = isotropic
tau_p                    = 12.0
compressibility          = 3e-4
ref_p                    = 1.0

; GENERATE VELOCITIES FOR STARTUP RUN:
gen_vel                  = no

    """
    if CreateMDP:
        with open(mdpfile, 'w') as fp:
            fp.write(line)


def GMX_run(run_file='run_gmx.sh', prefix='conf'):
    result = subprocess.run(
        f'bash {run_file}', shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print("Gromacs failed to run\n", os.getcwd())
        print("Standard Output:")
        print(result.stdout)
        print("Standard Error")
        print(result.stderr)
        raise ValueError("\n")

    # remove *#
    subprocess.run("rm *#", shell=True, capture_output=True, text=True)

    # process .xvg file to .dat file
    # energy
    energy = np.loadtxt(f'energy_{prefix}.xvg', comments=[
                        '#', '@'], ndmin=2)[:, 1:2]
    os.remove(f'energy_{prefix}.xvg')

    n_frames = energy.shape[0]
    energy_array = np.zeros(shape=(n_frames, 2))
    energy_array[:, 0] = np.arange(n_frames)
    energy_array[:, 1] = energy[:, 0]
    np.savetxt(f'energy_{prefix}.dat', energy_array)

    # forces
    forces = np.loadtxt(f'forces_{prefix}.xvg', comments=[
                        '#', '@'], ndmin=2)[:, 1:]
    n_frames = forces.shape[0]
    n_atoms = forces.size//n_frames//3

    forces_array = np.zeros(shape=(n_frames, n_atoms, 4))
    for i in range(n_frames):
        forces_array[i, :, 0] = i
        forces_array[i, :, 1:] = forces[i].reshape((n_atoms, 3))

    forces_array = forces_array.reshape((n_frames * n_atoms, 4))
    np.savetxt(f'forces_{prefix}.dat', forces_array)


def Clean(prefix='conf'):
    try:
        os.remove(f'energy_{prefix}.dat')
        os.remove(f'forces_{prefix}.dat')
    except:
        pass

def Load_energy(prefix='conf', clean=False):
    energy = np.loadtxt(f"energy_{prefix}.dat", ndmin=2)
    if clean:
        os.remove(f"energy_{prefix}.dat")
    return energy


def Load_forces(prefix='conf', clean=False):
    forces = np.loadtxt(f"forces_{prefix}.dat", ndmin=2)
    if clean:
        os.remove(f"forces_{prefix}.dat")
    return forces


def Compare_energy(energy1, energy2, isPrint=True):
    energy1 = float(energy1)
    energy2 = float(energy2)

    relative_energy_error = abs(energy1 - energy2) * \
        2 / (abs(energy1)+abs(energy2))
    abs_energy_error = abs(energy1 - energy2)

    if isPrint:
        print('###Energy Compare###')
        print(f'Absolute error: {abs_energy_error:.5e}')
        print(f'Relative error: {relative_energy_error:.5e}')

    if relative_energy_error <= 1e-5 or abs_energy_error <= 1e-3:
        if isPrint:
            print("Energies match!")
        return True
    else:
        if isPrint:
            print("Error: Energies do not match!")
        return False

def Compare_energy_array(energy_array1, energy_array2, isPrint=True):
    # shape check:
    if energy_array1.shape != energy_array2.shape:
        raise ValueError(
            f"The shape of energy arrays is not the same. {energy_array1.shape} vs {energy_array2.shape}")

    n_frames = energy_array1.shape[0]
    result_list = []
    for i in range(n_frames):
        frame = energy_array1[i, 0]
        energy1 = energy_array1[i, 1]
        energy2 = energy_array2[i, 1]
        if isPrint:
            print(f"###Frame {frame}###")
        result = Compare_energy(energy1, energy2, isPrint=isPrint)
        if isPrint:
            print()

        result_list.append(result)
    result_array = np.array(result_list)
    n_True = result_array.sum()
    n_False = n_frames-n_True

    if n_False != 0:
        print("Warning: Energy array does't match!!!")
        False_frame_list = energy_array1[:, 0][result_array == False].tolist()
        print(f"Mismatch frame list: ", *False_frame_list)
        return False
    else:
        print("Energy array all match!!!")
        return True


# def Compare_forces(forces1, forces2, isPrint=True):
#     '''forces1 and forces2 should be the calculated force and standard forces, respectively.'''
#     average = 0.5*np.linalg.norm(forces1, axis=1) + 0.5*np.linalg.norm(forces2, axis=1)

#     relative_force_error = np.linalg.norm(forces1 - forces2, axis=1) / average
#     relative_force_error = np.nan_to_num(relative_force_error, nan=0)
#     max_relative_force_error, max_relative_force_error_index = relative_force_error.max(
#     ), relative_force_error.argmax()

#     abs_force_error = np.linalg.norm(forces1 - forces2, axis=1)
#     max_abs_force_error, max_abs_force_error_index = abs_force_error.max(
#     ), abs_force_error.argmax()

#     atol = 1e-5
#     rtol = 1e-4
#     allclose = abs_force_error-(atol+rtol*average)
#     max_allclose = allclose.max()
#     if isPrint:
#         print('###Forces Compare###')
#         print(f'Max absolute error: {max_abs_force_error:.5f}')
#         print(f'Max relative error: {max_relative_force_error:.5f}')
#         print(f'      Max allclose: {max_allclose:.5f}')

#     if max_allclose <= 0:
#         if isPrint:
#             print("Forces match!")
#         return True
#     else:
#         if isPrint:
#             print("Error: Forces do not match!")
#         return False 

def Compare_forces(forces1, forces2, isPrint=True):
    '''forces1 and forces2 should be the calculated force and standard forces, respectively.'''
    # Don't consider zero forces (virtual site atoms)
    nonzero_index = (np.linalg.norm(forces1, axis=1) != 0) & (np.linalg.norm(forces2, axis=1) != 0)
    forces1 = forces1[nonzero_index]
    forces2 = forces2[nonzero_index]

    ref = np.linalg.norm(forces2, axis=1)

    relative_force_errors = np.linalg.norm(forces1 - forces2, axis=1) / ref
    relative_force_errors = np.nan_to_num(relative_force_errors, nan=0)
    abs_force_errors = np.linalg.norm(forces1 - forces2, axis=1)


    atol = 1e-4
    rtol = 1e-5
    allclose = abs_force_errors-(atol+rtol*ref)
    max_allclose = allclose.max()
    max_allclose_index = allclose.argmax()
    abs_force_error = abs_force_errors[max_allclose_index]
    relative_force_error = relative_force_errors[max_allclose_index]

    if isPrint:
        print('###Forces Compare###')
        print(f'   Max allclose: {max_allclose:.5e}')
        print(f' Absolute error: {abs_force_error:.5e}')
        print(f' Relative error: {relative_force_error:.5e}')


    if max_allclose <= 0:
        if isPrint:
            print("Forces match!")
        return True
    else:
        if isPrint:
            print("Error: Forces do not match!")
        return False 

def Compare_forces_array(forces_array1, forces_array2, isPrint=True):
    # shape check:
    if forces_array1.shape != forces_array2.shape:
        raise ValueError(
            f"The shape of energy arrays is not the same. {forces_array1.shape} vs {forces_array2.shape}")

    frame_array = np.unique(forces_array1[:, 0])
    n_frames = len(frame_array)
    result_list = []
    for frame in frame_array:
        forces1 = forces_array1[forces_array1[:, 0] == frame, 1:]
        forces2 = forces_array2[forces_array2[:, 0] == frame, 1:]

        if isPrint:
            print(f"###Frame {frame}###")
        result = Compare_forces(forces1, forces2, isPrint=isPrint)
        if isPrint:
            print()
        result_list.append(result)

    result_array = np.array(result_list)
    n_True = result_array.sum()
    n_False = n_frames-n_True

    if n_False != 0:
        print("Warning: Forces array does't match!!!")
        False_frame_list = frame_array[result_array == False].astype(np.int64).tolist()
        print(f"Mismatch frame list: ", *False_frame_list)
        return False
    else:
        print("Forces array all match!!!")
        return True