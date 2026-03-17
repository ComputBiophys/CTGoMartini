"""
Author: Song Yang
Date: 2025-05-24
Version: 0.1.0
"""
import os
import numpy as np
import argparse
import multiprocessing
import MDAnalysis as mda
import time
from typing import List, Dict, Optional, Tuple
from PositionExtractor import PositionExtractor



class ReferenceDistanceCalculator:
    def __init__(self, params: Dict):
        self.params = params
        # Pre-calculate parameters for faster access
        self.min_dist = params['minimum_distance']
        self.max_dist = params['maximum_distance']
        self.min_diff = params['minimum_difference']
        self.excl_res = params['excluded_residues']
    
    @staticmethod
    def calculate_distance(p1: np.ndarray, p2: np.ndarray) -> np.ndarray:
        """Calculate distance between two points or arrays of points."""
        return np.sqrt(np.sum((p1 - p2)**2, axis=-1))
    
    @staticmethod
    def get_difference_distances(distances: np.ndarray) -> np.ndarray:
        """Calculate pairwise differences between distances."""
        return np.abs(distances[:, None] - distances)
    
    def validate_references(self, ref_atoms_list: List[mda.AtomGroup]) -> None:
        """Validate reference structures have consistent atoms."""
        if len(ref_atoms_list) < 2:
            raise ValueError("Reference structures should have at least two states!")
        
        lengths = [len(atoms) for atoms in ref_atoms_list]
        if len(set(lengths)) > 1:
            raise ValueError("Reference structures have different atom counts!")
    
    def calculate_reference_distances(self, ref_atoms_list: List[mda.AtomGroup]) -> List[np.ndarray]:
        """
        Calculate reference distances between atom pairs.
        Returns list of numpy arrays with columns: atomid1, atomid2, ref_d
        """
        self.validate_references(ref_atoms_list)
        n_atoms = len(ref_atoms_list[0])
        n_refs = len(ref_atoms_list)
        
        # Pre-extract all positions for faster access
        all_positions = np.array([ref_atoms.positions for ref_atoms in ref_atoms_list])
        
        # Get atom properties for exclusion checks (only need from first reference)
        segids = ref_atoms_list[0].segids
        resids = ref_atoms_list[0].resids
        atoms = ref_atoms_list[0]

        # Prepare storage for results as list of lists
        results = [[] for _ in range(n_refs)]
        
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                # Fast exclusion checks using pre-loaded properties
                if segids[i] == segids[j] and resids[j] - resids[i] < self.excl_res:
                    continue
                
                # Vectorized distance calculation across all references
                distances = self.calculate_distance(all_positions[:, i], all_positions[:, j])
                
                # Check distance range validity
                valid_mask = (distances >= self.min_dist) & (distances <= self.max_dist)
                
                if not np.any(valid_mask):
                    continue
                
                # Check minimum difference requirement
                diffs = self.get_difference_distances(distances)
                if diffs[np.triu_indices_from(diffs, k=1)].min() < self.min_diff:
                    continue
                
                # Store valid distances
                for ref_idx in range(n_refs):
                    if valid_mask[ref_idx]:
                        results[ref_idx].append([atoms[i].ix, atoms[j].ix, distances[ref_idx]])
        
        # save results
        # for i, result in enumerate(results):
        #     np.savetxt(f"reference_distances_{i}.dat", np.array(result), fmt="%d %d %.3f")

        # Convert each reference's results to numpy array
        return [np.array(res) for res in results]

def test_reference_distance_calculation(ref_files: List[str], params: Dict) -> Tuple[float, List[np.ndarray]]:
    """Test function for calculating reference distances with optimized loading."""
    start_time = time.time()
    
    # Load reference structures
    ref_atoms = [mda.Universe(f).select_atoms(params['selected_atom']) for f in ref_files]
    
    # Calculate reference distances
    ref_calculator = ReferenceDistanceCalculator(params)
    ref_distances_list = ref_calculator.calculate_reference_distances(ref_atoms)
    
    elapsed = time.time() - start_time
    print(f"Reference distance calculation took {elapsed:.2f} seconds")
    
    # Print some statistics
    for i, arr in enumerate(ref_distances_list):
        print(f"State {i}: {len(arr)} reference distances calculated")
        if len(arr) > 0:
            print(f"  Distance stats - min: {arr[:, 2].min():.2f}, "
                  f"max: {arr[:, 2].max():.2f}, "
                  f"mean: {arr[:, 2].mean():.2f}")
    
    return elapsed, ref_distances_list

# if __name__ == '__main__':
#     # Example usage with timing
#     elapsed, ref_distances_list = test_reference_distance_calculation(
#         ["Up/Up_cg.pdb", "Down/Down_cg.pdb"],
#         params={
#             'minimum_distance': 6.0,
#             'maximum_distance': 50.0,
#             'minimum_difference': 5,
#             'excluded_residues': 4,
#             'selected_atom': 'name BB'
#         }
#     )


class DRMSCalculator:
    def __init__(self, position_extractor: PositionExtractor):
        self.position_extractor = position_extractor
        self._results = []
        self.replica_indices = None
        self.skip = None

    @staticmethod
    def _worker_process_chunk(
        netcdf_path: str,
        checkpoint_path: str,
        ref_distances_list: List[np.ndarray],
        frame_chunk_indices: List[int],
        replica_indices: np.ndarray,
        dt: float
    ) -> List[Tuple[np.ndarray, np.ndarray]]:
        """
        Worker function to process a chunk of frames and calculate dRMS.
        Creates its own PositionExtractor instance to avoid pickling issues.
        
        Args:
            netcdf_path: Path to the NetCDF trajectory file
            checkpoint_path: Path to the checkpoint file
            ref_distances_list: List of reference distance matrices
            frame_chunk_indices: List of frame indices to process
            replica_indices: Array of replica indices to analyze
            dt: Time step between frames (ps)
            
        Returns:
            List of tuples (times, drms) for each reference state
        """
        chunk_start_time = time.time()

        extractor = PositionExtractor(netcdf_path, checkpoint_path)
        n_states = len(ref_distances_list)
        results = []

        try:
            # Precompute all atom groups needed
            atom_groups = []
            for ref_distances in ref_distances_list:
                atom_groups.extend([np.array(ref_distances[:, 0],dtype=int), np.array(ref_distances[:, 1], dtype=int)])

            # Get all positions in one call
            positions_groups = extractor.get_positions_for_atom_groups(
                atom_groups, frame_chunk_indices, replica_indices
            )

            # Process each state
            for i in range(n_states):
                pos1, pos2 = positions_groups[2*i], positions_groups[2*i+1]
                ref_dist = ref_distances_list[i][:, 2]

                # Calculate distances and dRMS
                distances = np.linalg.norm(pos1 - pos2, axis=-1)
                drms = np.sqrt(
                    np.mean((distances - ref_dist)**2, axis=-1)
                )
                        
                # Convert frame indices to time in ps
                times = np.array(frame_chunk_indices) * dt
                
                results.append((times, drms))
        finally:
            extractor.close()

        chunk_time = time.time() - chunk_start_time
        print(f"Processed frames {frame_chunk_indices[0]}-{frame_chunk_indices[-1]} in {chunk_time:.2f} seconds")
        return results
    
    def calculate_trajectory_dRMS(
        self,
        ref_distances_list: List[np.ndarray],
        replica_indices: Optional[np.ndarray] = None,
        skip: int = 1,
        chunk_size: int = 500,
        num_workers: Optional[int] = None
    ) -> List[Tuple[np.ndarray, np.ndarray]]:
        """
        Calculate dRMS for entire trajectory using multiprocessing.
        
        Args:
            ref_distances_list: List of reference distance matrices
            replica_indices: Array of replica indices to analyze (None for all)
            skip: Process every skip-th frame
            chunk_size: Number of frames per worker task
            num_workers: Number of parallel workers (default: cpu_count)
            
        Returns:
            List of (time, drms) tuples for each reference state
        """
        self.replica_indices = replica_indices
        self.skip = skip

        # Prepare frame indices to process
        total_frames = self.position_extractor.n_frames
        frames_to_process = np.arange(0, total_frames, skip, dtype=int)
        n_frames = len(frames_to_process)

        if not n_frames:
            print("Warning: No frames to process")
            return []

        # Prepare tasks for multiprocessing
        tasks = []
        for i in range(0, n_frames, chunk_size):
            current_chunk_indices = range(i, min(i + chunk_size, n_frames))
            chunk_indices = [frames_to_process[k] for k in current_chunk_indices]

            tasks.append((
                self.position_extractor.main_reporter_nc_path,
                self.position_extractor.checkpoint_nc_path,
                ref_distances_list,
                chunk_indices,
                replica_indices,
                self.position_extractor.dt
            ))

        print(f"Processing {n_frames} frames in {len(tasks)} chunks...")
        
        # Process tasks with multiprocessing
        if num_workers is None:
            num_workers = multiprocessing.cpu_count()

        with multiprocessing.Pool(num_workers) as pool:
            chunk_results = pool.starmap(self._worker_process_chunk, tasks)

        # Combine results from all chunks
        self._results = []
        n_states = len(ref_distances_list)
        
        for state_idx in range(n_states):
            all_times = np.concatenate([r[state_idx][0] for r in chunk_results])
            all_drms = np.vstack([r[state_idx][1] for r in chunk_results])
            
            self._results.append(
                (all_times, all_drms)
            )


    def save_results(self, output_prefix: str='dRMSTraj'):
        """
        Save distance calculation results to file.
        
        Args:
            output_prefix: Prefix for output file
        """

        if self._results is None:
            raise RuntimeError("No results to save. Run calculate_trajectory_dRMS() first.")
        
        n_states = len(self._results)
        states_str = 'ABCDEFGHIJLKMN'

        for i in range(n_states):
            all_times, all_distances = self._results[i]
            
            # Prepare header information
            header_replica_cols = " ".join([f"Replica_{i}" for i in self.replica_indices])
            header_text = (f'time(ps), {header_replica_cols}\n'
                        f'NetCDF: {self.position_extractor.main_reporter_nc_path}\n'
                        f'Checkpoint: {self.position_extractor.checkpoint_nc_path}\n'
                        f'Replicas processed: {",".join(map(str, self.replica_indices))}\n'
                        f'Skip: {self.skip}, dt: {self.position_extractor.dt:.4f} ps')
            
            # Format specifications
            fmt_list = ['%.1f'] + ['%.4f'] * len(self.replica_indices)
            
            # Combine time and distance data
            data = np.hstack((all_times.reshape(-1, 1), all_distances))
            
            # Create directory if it doesn't exist
            output_file = f"{output_prefix}_State{states_str[i]}.dat"
            if os.path.dirname(output_file): 
                os.makedirs(os.path.dirname(output_file))
            
            # Save results
            np.savetxt(output_file, data, header=header_text, comments='# ', fmt=fmt_list)
            print(f"Results of State {states_str[i]} saved to {output_file}")

if __name__ == '__main__':
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description='Calculate dRMS between trajectory and reference structures.',
        epilog='Example: python dRMS_analysis_nc.py -s structure.pdb -nc output.nc --ref_str ref1.pdb ref2.pdb --sel "name BB" --num_workers 20 -o dRMSTraj_nc',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('-s', '--structure', required=True,
                      help='The structure file (pdb, gro)')
    parser.add_argument('-nc', '--netcdf', required=True,
                      help='The NetCDF trajectory file')
    parser.add_argument('-ref', '--ref_str', required=True, nargs='+',
                      help='Reference structure files (at least 2 required)')
    
    # Optional arguments with defaults
    parser.add_argument('-c', '--checkpoint', default='output_checkpoint.nc',
                      help='Checkpoint NetCDF file')
    parser.add_argument('-o', '--output_prefix', default='dRMSTraj_nc',
                      help='Prefix for output files')
    parser.add_argument('--sel', '--selected_atom', default='name BB',
                      help='Atom selection string for analysis')
    parser.add_argument('--min_dist', type=float, default=6.0,
                      help='Minimum distance between atom pairs (Å)')
    parser.add_argument('--max_dist', type=float, default=50.0,
                      help='Maximum distance between atom pairs (Å)')
    parser.add_argument('--min_diff', type=float, default=5.0,
                      help='Minimum difference between reference distances (Å)')
    parser.add_argument('--excl_res', type=int, default=4,
                      help='Exclude residues within this sequence separation')
    parser.add_argument('--skip', type=int, default=1,
                      help='Process every N-th frame')
    parser.add_argument('--chunk', type=int, default=10,
                      help='Frames per chunk for parallel processing')
    parser.add_argument('--num_workers', type=int, default=None,
                      help='Number of worker processes (default: all CPUs)')
    parser.add_argument('--replicas', type=str, default="all",
                      help='Comma-separated replica indices or "all"')
    
    args = parser.parse_args()
    
    # Create parameter dictionary
    params = {
        'minimum_distance': args.min_dist,
        'maximum_distance': args.max_dist,
        'minimum_difference': args.min_diff,
        'excluded_residues': args.excl_res,
        'selected_atom': args.sel
    }
    
    # 1. Calculate reference distances
    print("\nCalculating reference distances...")
    start_time = time.time()
    # Load reference structures
    ref_atoms = [mda.Universe(f).select_atoms(params['selected_atom']) for f in args.ref_str]
    # Calculate reference distances
    ref_calculator = ReferenceDistanceCalculator(params)
    ref_distances_list = ref_calculator.calculate_reference_distances(ref_atoms)
    elapsed = time.time() - start_time
    print(f"Reference distance calculation took {elapsed:.2f} seconds")
    
    # 2. dRMS analysis
    # Initialize position extractor and dRMS calculator
    position_extractor = PositionExtractor(args.netcdf, args.checkpoint)
    drms_calculator = DRMSCalculator(position_extractor)
    
    # Parse replica indices
    replica_indices = np.arange(position_extractor.n_replicas, dtype=int) if args.replicas == "all" else np.array([int(i) for i in args.replicas.split(',')])
    
    print("\nTrajectory analysis setup:")
    print(f"- Structure: {args.structure}")
    print(f"- NetCDF: {args.netcdf} ({position_extractor.n_frames} frames, {position_extractor.n_replicas} replicas)")
    print(f"- Reference structures: {len(args.ref_str)} states")
    print(f"- Replicas to analyze: {replica_indices}")
    print(f"- Atom selection: {args.sel}")
    print(f"- Distance range: {args.min_dist}-{args.max_dist} Å")
    print(f"- Minimum difference: {args.min_diff} Å")
    print(f"- Residue exclusion: {args.excl_res}")
    print(f"- Processing every {args.skip} frame(s)")
    print(f"- Chunk size: {args.chunk} frames")
    print(f"- Workers: {args.num_workers if args.num_workers else 'all CPUs'}")
    position_extractor.close()    
    
    # Calculate dRMS for trajectory
    print("\nCalculating trajectory dRMS...")
    drms_calculator.calculate_trajectory_dRMS(
        ref_distances_list=ref_distances_list,
        replica_indices=replica_indices,
        skip=args.skip,
        chunk_size=args.chunk,
        num_workers=args.num_workers
    )
    
    # Save results
    drms_calculator.save_results(output_prefix=args.output_prefix)
    
    print("\nAnalysis complete!")

