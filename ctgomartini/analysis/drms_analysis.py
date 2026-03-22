

import numpy as np
import argparse
import MDAnalysis as mda
from multiprocessing import Pool
from functools import partial
import time


class DRMSAnalyzer:
    def __init__(self):
        self.drms_trajectories = None
        self.reference_distances = None
        self.start_time = None

    def _report_time(self):
        """Report elapsed time since start."""
        elapsed = time.time() - self.start_time
        print(f"Elapsed Time: {elapsed:.2f} seconds")

    def _calculate_distances(self, points1, points2):
        """Calculate distances between pairs of points."""
        points1 = np.array(points1)
        points2 = np.array(points2)
        return np.linalg.norm(points1 - points2, axis=1)

    def _check_same_residue(self, atom1, atom2):
        """Check if two atoms belong to the same residue."""
        return (atom1.segid == atom2.segid) and (atom1.resid == atom2.resid)

    def _calculate_distance_differences(self, distances):
        """Calculate absolute differences between all pairs of distances."""
        return [abs(distances[i] - distances[j]) 
                for i in range(len(distances)) 
                for j in range(i+1, len(distances))]

    def _validate_references(self, ref_atoms):
        """Validate reference structures consistency."""
        if len(ref_atoms) < 2:
            raise ValueError("At least two reference states required!")
            
        base_count = len(ref_atoms[0])
        for i, atoms in enumerate(ref_atoms[1:]):
            if len(atoms) != base_count:
                raise ValueError(f"Reference structure {i+2} has different atom count!")
                
            for j in range(base_count):
                if not self._check_same_residue(ref_atoms[0][j], atoms[j]):
                    raise ValueError(f"Atom {j} segid/resid mismatch in references")

    def calculate_reference_distances(self, ref_files, min_distance=6.0, max_distance=50.0, 
                                   min_difference=5.0, excluded_residues=4, atom_selection='name BB'):
        """Calculate reference distances between atom pairs."""
        # Load reference structures
        self.start_time = time.time()
        ref_universes = [mda.Universe(file) for file in ref_files]
        ref_atoms = [u.select_atoms(atom_selection) for u in ref_universes]
        
        self._validate_references(ref_atoms)
        self.reference_distances = self._process_reference_distances(
            ref_atoms, min_distance, max_distance, min_difference, excluded_residues)
        print("Reference distance calculation complete!")
        self._report_time()
        # return self.reference_distances

    def _process_reference_distances(self, ref_atoms, min_distance, max_distance, 
                                   min_difference, excluded_residues):
        """Process and filter reference distances."""
        num_refs = len(ref_atoms)
        num_atoms = len(ref_atoms[0])
        distances = [[] for _ in range(num_refs)]

        for i in range(num_atoms):
            for j in range(i+1, num_atoms):
                atom_i, atom_j = ref_atoms[0][i], ref_atoms[0][j]
                
                # Skip residues too close in sequence
                if (atom_i.segid == atom_j.segid and 
                    atom_i.resid + excluded_residues > atom_j.resid):
                    continue
                
                # Calculate distances across all references
                ref_dists = [self._calculate_distances([ref[i].position], [ref[j].position])[0]
                            for ref in ref_atoms]
                min_diff = min(self._calculate_distance_differences(ref_dists))

                # Apply distance filters
                if min_diff >= min_difference:
                    for k, dist in enumerate(ref_dists):
                        if min_distance <= dist <= max_distance:
                            distances[k].append([i, j, dist])

        # Convert to numpy arrays
        return [np.array(dist) for dist in distances]

    def analyze_trajectory(self, structure, trajectory, atom_selection='name BB', cores=10, skip=1):
        """Calculate dRMS for a trajectory against references."""
        if not self.reference_distances:
            raise ValueError("Reference distances not calculated!")

        # Load trajectory
        universe = mda.Universe(structure, trajectory)
        atoms = universe.select_atoms(atom_selection)

        self.start_time = time.time()
        self.drms_trajectories = []
        
        for ref_dist in self.reference_distances:
            results = self._process_trajectory(atoms, ref_dist, cores, skip)
            self.drms_trajectories.append(results)
        
        print("Trajectory analysis complete!")
        self._report_time()
        return self.drms_trajectories

    def _calculate_frame_drms(self, frame, atoms, ref_dist):
        """Calculate dRMS for a single frame."""
        atoms.universe.trajectory[frame]
        points1 = atoms[ref_dist[:, 0].astype(int)].positions
        points2 = atoms[ref_dist[:, 1].astype(int)].positions
        dist = self._calculate_distances(points1, points2)
        return [atoms.universe.trajectory.time, np.sqrt(np.mean((dist - ref_dist[:, 2])**2))]

    def _process_trajectory(self, atoms, ref_dist, cores=10, skip=1):
        """Process trajectory frames in parallel."""
        frame_func = partial(self._calculate_frame_drms, atoms=atoms, ref_dist=ref_dist)
        frames = range(0, atoms.universe.trajectory.n_frames, skip)
        
        with Pool(cores) as pool:
            results = pool.map(frame_func, frames)
        return np.array(results)

    def save_results(self, prefix='drms'):
        """Save analysis results to files."""
        if not self.drms_trajectories:
            raise ValueError("No data to save!")

        for i, trajectory in enumerate(self.drms_trajectories):
            filename = f"{prefix}_state{chr(65 + i)}.dat"
            np.savetxt(filename, trajectory)


def parse_arguments():
    """Parse and return command line arguments."""
    parser = argparse.ArgumentParser(
        description="Calculate dRMS between trajectory and references",
        epilog='Example: python -m ctgomartini.analysis.drms_analysis -s structure.pdb -f trajectory.xtc -r ref1.pdb ref2.pdb -sel "name BB" -n 10 -prefix dRMSTraj',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        )
    
    parser.add_argument('-s', '--structure', required=True,
                      help='Input structure file (gro/pdb)')
    parser.add_argument('-f', '--trajectory', required=True,
                      help='Input trajectory file (xtc/dcd)')
    parser.add_argument('-r', '--references', required=True, nargs='+',
                      help='Reference structure files (gro/pdb)')
    parser.add_argument('--skip', type=int, default=1,
                      help='Frame skipping interval')
    parser.add_argument('-sel', '--selection', default='name BB',
                      help='Atom selection criteria')
    parser.add_argument('-prefix', '--output_prefix', default='drms',
                      help='Output filename prefix')
    parser.add_argument('-n', '--cores', type=int, default=10,
                      help='Number of CPU cores to use')

    return parser.parse_args()


def main():
    args = parse_arguments()
    analyzer = DRMSAnalyzer()
    
    # Default parameters moved to main function
    default_params = {
        'min_distance': 6.0,
        'max_distance': 50.0,
        'min_difference': 5.0,
        'excluded_residues': 4,
        'atom_selection': args.selection,
        'cores': args.cores,
        'skip': args.skip
    }
    
    print("Processing reference distances...")
    analyzer.calculate_reference_distances(
        args.references,
        min_distance=default_params['min_distance'],
        max_distance=default_params['max_distance'],
        min_difference=default_params['min_difference'],
        excluded_residues=default_params['excluded_residues'],
        atom_selection=default_params['atom_selection']
    )
    
    print("Analyzing trajectory...")
    analyzer.analyze_trajectory(
        args.structure,
        args.trajectory,
        atom_selection=default_params['atom_selection'],
        cores=default_params['cores'],
        skip=default_params['skip']
    )
    
    print("Saving results...")
    analyzer.save_results(prefix=args.output_prefix)
    print("Analysis completed successfully.")


if __name__ == '__main__':
    main()