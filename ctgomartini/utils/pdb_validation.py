"""
PDB structure validation module for CTGoMartini.

Provides comprehensive validation of PDB file compatibility for multiple-basin
topology generation. Uses MDTraj for reliable structure parsing.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import sys

try:
    import mdtraj as md
except ImportError:
    raise ImportError(
        "MDTraj is required for PDB validation. "
        "Install with: conda install -c conda-forge mdtraj"
    )


@dataclass
class ResidueInfo:
    """Detailed information about a residue in a PDB structure."""
    index: int  # 0-based index in the structure
    res_seq: int  # Residue sequence number from PDB
    res_name: str  # 3-letter residue code
    chain_id: str  # Chain identifier
    atoms: list[str] = field(default_factory=list)  # List of atom names
    
    @property
    def atom_count(self) -> int:
        return len(self.atoms)
    
    def __repr__(self) -> str:
        return f"{self.res_name}{self.res_seq}({self.chain_id}):{self.atom_count}atoms"


@dataclass  
class StructureSummary:
    """Complete summary of a PDB structure."""
    filepath: Path
    name: str  # State name for display
    residues: list[ResidueInfo]
    total_atoms: int
    chains: list[str]
    n_residues: int
    
    @property
    def residue_sequence(self) -> list[tuple[str, int, str]]:
        """Return list of (res_name, res_seq, chain_id) tuples."""
        return [(r.res_name, r.res_seq, r.chain_id) for r in self.residues]
    
    @property
    def chain_residue_counts(self) -> dict[str, int]:
        """Count residues per chain."""
        counts: dict[str, int] = {}
        for res in self.residues:
            counts[res.chain_id] = counts.get(res.chain_id, 0) + 1
        return counts


@dataclass
class MismatchDetail:
    """Detailed information about a specific mismatch."""
    position: int  # 1-based position in sequence
    ref_residue: ResidueInfo | None
    other_residue: ResidueInfo | None
    mismatch_type: str  # 'missing', 'extra', 'name_mismatch', 'chain_mismatch'
    description: str


@dataclass
class ValidationReport:
    """Comprehensive validation report for PDB compatibility."""
    is_valid: bool
    ref_name: str
    other_name: str
    ref_file: Path
    other_file: Path
    
    # Detailed findings
    residue_count_match: bool = True
    sequence_match: bool = True
    chain_structure_match: bool = True
    
    # Specific issues
    mismatches: list[MismatchDetail] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    
    # Statistics
    ref_residue_count: int = 0
    other_residue_count: int = 0
    matching_residues: int = 0
    
    def format_error(self) -> str:
        """Format a detailed error message."""
        if self.is_valid:
            return ""
        
        lines = [
            "\n" + "=" * 75,
            "  PDB COMPATIBILITY ERROR",
            "=" * 75,
            f"\nReference: {self.ref_name}",
            f"  File: {self.ref_file}",
            f"  Residues: {self.ref_residue_count}",
            f"\nComparison: {self.other_name}",
            f"  File: {self.other_file}",
            f"  Residues: {self.other_residue_count}",
        ]
        
        # Main issue
        if not self.residue_count_match:
            lines.extend(self._format_residue_count_error())
        elif not self.chain_structure_match:
            lines.extend(self._format_chain_error())
        elif not self.sequence_match:
            lines.extend(self._format_sequence_error())
        
        # Warnings
        if self.warnings:
            lines.extend([
                "\n⚠️  WARNINGS:",
            ])
            for warning in self.warnings:
                lines.append(f"  - {warning}")
        
        lines.append("=" * 75)
        return "\n".join(lines)
    
    def _format_residue_count_error(self) -> list[str]:
        """Format residue count mismatch details."""
        diff = abs(self.ref_residue_count - self.other_residue_count)
        lines = [
            f"\n❌ RESIDUE COUNT MISMATCH (difference: {diff})",
            "\nDetailed Analysis:",
        ]
        
        # Analyze specific mismatches
        missing_in_other = []
        extra_in_other = []
        
        for mismatch in self.mismatches:
            if mismatch.mismatch_type == 'missing':
                missing_in_other.append(mismatch)
            elif mismatch.mismatch_type == 'extra':
                extra_in_other.append(mismatch)
        
        if missing_in_other:
            lines.append(f"\n  Residues present in '{self.ref_name}' but missing in '{self.other_name}':")
            for m in missing_in_other[:10]:  # Show first 10
                r = m.ref_residue
                lines.append(f"    - Position {m.position}: {r.res_name}{r.res_seq} (chain {r.chain_id})")
            if len(missing_in_other) > 10:
                lines.append(f"    ... and {len(missing_in_other) - 10} more")
        
        if extra_in_other:
            lines.append(f"\n  Residues present in '{self.other_name}' but not in '{self.ref_name}':")
            for m in extra_in_other[:10]:
                r = m.other_residue
                lines.append(f"    - Position {m.position}: {r.res_name}{r.res_seq} (chain {r.chain_id})")
            if len(extra_in_other) > 10:
                lines.append(f"    ... and {len(extra_in_other) - 10} more")
        
        return lines
    
    def _format_chain_error(self) -> list[str]:
        """Format chain structure mismatch details."""
        return [
            "\n❌ CHAIN STRUCTURE MISMATCH",
            "\nChain IDs in reference: " + str(getattr(self, '_ref_chains', 'N/A')),
            "Chain IDs in comparison: " + str(getattr(self, '_other_chains', 'N/A')),
        ]
    
    def _format_sequence_error(self) -> list[str]:
        """Format sequence mismatch details."""
        lines = [
            "\n❌ RESIDUE SEQUENCE MISMATCH",
            f"\nFirst {min(len(self.mismatches), 10)} differences:",
        ]
        
        for mismatch in self.mismatches[:10]:
            pos = mismatch.position
            ref = mismatch.ref_residue
            other = mismatch.other_residue
            
            ref_str = f"{ref.res_name}{ref.res_seq}" if ref else "<missing>"
            other_str = f"{other.res_name}{other.res_seq}" if other else "<missing>"
            
            lines.append(f"  Position {pos}:")
            lines.append(f"    {self.ref_name}:  {ref_str} (chain {ref.chain_id if ref else '?'})")
            lines.append(f"    {self.other_name}: {other_str} (chain {other.chain_id if other else '?'})")
            
            # Provide specific diagnosis
            if ref and other:
                if ref.res_name != other.res_name:
                    lines.append(f"    → Different residue types: {ref.res_name} ≠ {other.res_name}")
                    if self._is_similar_residue(ref.res_name, other.res_name):
                        lines.append(f"    → Note: These are similar residues (possible mutation)")
                elif ref.res_seq != other.res_seq:
                    lines.append(f"    → Different numbering: {ref.res_seq} ≠ {other.res_seq}")
                    if abs(ref.res_seq - other.res_seq) == 1:
                        lines.append(f"    → Off-by-one error in residue numbering")
        
        if len(self.mismatches) > 10:
            lines.append(f"  ... and {len(self.mismatches) - 10} more differences")
        
        return lines
    
    @staticmethod
    def _is_similar_residue(res1: str, res2: str) -> bool:
        """Check if two residues are similar (e.g., mutations)."""
        # Similar amino acid groups
        similar_groups = [
            {'ALA', 'VAL', 'LEU', 'ILE', 'MET'},  # Hydrophobic aliphatic
            {'PHE', 'TYR', 'TRP'},  # Aromatic
            {'SER', 'THR'},  # Hydroxyl
            {'ASP', 'GLU'},  # Acidic
            {'LYS', 'ARG', 'HIS'},  # Basic
            {'ASN', 'GLN'},  # Amide
        ]
        
        for group in similar_groups:
            if res1 in group and res2 in group:
                return True
        return False


class PDBCompatibilityValidator:
    """
    Validator for checking PDB file compatibility.
    
    Use this class to validate that multiple PDB structures are compatible
    for generating multiple-basin potential topologies.
    """
    
    # Residues that commonly differ due to protonation
    PROTONATION_VARIANTS = {
        'HIS', 'HSD', 'HSE', 'HSP', 'HID', 'HIE', 'HIP',
        'ASP', 'ASH',  # Aspartic acid vs protonated
        'GLU', 'GLH',  # Glutamic acid vs protonated
        'LYS', 'LYN',  # Lysine vs neutral
        'CYS', 'CYX', 'CYM',  # Cysteine variants
    }
    
    def __init__(self, verbose: bool = True):
        self.verbose = verbose
    
    def validate(
        self,
        pdb_files: list[str | Path],
        state_names: list[str] | None = None
    ) -> list[ValidationReport]:
        """
        Validate compatibility of multiple PDB files.
        
        Args:
            pdb_files: List of PDB file paths to validate
            state_names: Optional list of names for each PDB (for reporting)
            
        Returns:
            List of ValidationReport objects, one for each comparison
            against the first (reference) structure.
            
        Raises:
            FileNotFoundError: If any PDB file doesn't exist
            ValueError: If structures are incompatible (with detailed message)
        """
        if len(pdb_files) < 2:
            return []
        
        if state_names is None:
            state_names = [f"Structure{i+1}" for i in range(len(pdb_files))]
        
        # Parse all structures
        summaries = []
        for pdb_file, name in zip(pdb_files, state_names):
            path = Path(pdb_file)
            if not path.exists():
                raise FileNotFoundError(f"PDB file not found: {path.absolute()}")
            summaries.append(self._parse_structure(path, name))
        
        # Compare each structure against the reference
        reports = []
        ref_summary = summaries[0]
        
        for other_summary in summaries[1:]:
            report = self._compare_structures(ref_summary, other_summary)
            reports.append(report)
            
            if not report.is_valid:
                raise ValueError(report.format_error())
        
        return reports
    
    def _parse_structure(self, filepath: Path, name: str) -> StructureSummary:
        """Parse a PDB file using MDTraj."""
        try:
            traj = md.load_pdb(str(filepath))
        except Exception as e:
            raise ValueError(f"Failed to load PDB file {filepath}: {e}")
        
        topology = traj.topology
        residues = []
        
        for idx, residue in enumerate(topology.residues):
            # Get atom names
            atoms = [atom.name for atom in residue.atoms]
            
            res_info = ResidueInfo(
                index=idx,
                res_seq=residue.resSeq,
                res_name=residue.name,
                chain_id=residue.chain.chain_id if residue.chain else ' ',
                atoms=atoms
            )
            residues.append(res_info)
        
        chains = list(dict.fromkeys(r.chain_id for r in residues))  # Preserve order
        
        return StructureSummary(
            filepath=filepath,
            name=name,
            residues=residues,
            total_atoms=topology.n_atoms,
            chains=chains,
            n_residues=topology.n_residues
        )
    
    def _compare_structures(
        self,
        ref: StructureSummary,
        other: StructureSummary
    ) -> ValidationReport:
        """Compare two structures and generate a detailed report."""
        report = ValidationReport(
            is_valid=True,
            ref_name=ref.name,
            other_name=other.name,
            ref_file=ref.filepath,
            other_file=other.filepath,
            ref_residue_count=ref.n_residues,
            other_residue_count=other.n_residues
        )
        
        # Check residue count
        if ref.n_residues != other.n_residues:
            report.is_valid = False
            report.residue_count_match = False
        
        # Check chain structure
        if ref.chains != other.chains:
            report.chain_structure_match = False
            if not report.is_valid:
                report.is_valid = False
        
        # Detailed residue-by-residue comparison
        mismatches = self._find_residue_mismatches(ref, other)
        report.mismatches = mismatches
        
        if mismatches:
            report.sequence_match = False
            report.is_valid = False
        
        # Calculate matching residues
        report.matching_residues = min(ref.n_residues, other.n_residues) - len(mismatches)
        
        # Check for atom count differences (warnings only)
        self._check_atom_differences(ref, other, report)
        
        # Check for protonation differences
        self._check_protonation_differences(ref, other, report)
        
        return report
    
    def _find_residue_mismatches(
        self,
        ref: StructureSummary,
        other: StructureSummary
    ) -> list[MismatchDetail]:
        """Find detailed mismatches between two residue sequences."""
        mismatches = []
        
        # Create dictionaries for lookup by (chain, res_seq)
        ref_dict = {(r.chain_id, r.res_seq): r for r in ref.residues}
        other_dict = {(r.chain_id, r.res_seq): r for r in other.residues}
        
        # Find all unique keys
        all_keys = set(ref_dict.keys()) | set(other_dict.keys())
        
        for key in sorted(all_keys):
            chain_id, res_seq = key
            ref_res = ref_dict.get(key)
            other_res = other_dict.get(key)
            
            # Determine position (use index from reference if available)
            if ref_res:
                position = ref_res.index + 1
            elif other_res:
                position = other_res.index + 1
            else:
                position = 0
            
            if ref_res is None:
                # Extra residue in other
                mismatches.append(MismatchDetail(
                    position=position,
                    ref_residue=None,
                    other_residue=other_res,
                    mismatch_type='extra',
                    description=f"Extra residue {other_res.res_name}{other_res.res_seq} in {other.name}"
                ))
            elif other_res is None:
                # Missing residue in other
                mismatches.append(MismatchDetail(
                    position=position,
                    ref_residue=ref_res,
                    other_residue=None,
                    mismatch_type='missing',
                    description=f"Missing residue {ref_res.res_name}{ref_res.res_seq} in {other.name}"
                ))
            elif ref_res.res_name != other_res.res_name:
                # Name mismatch
                mismatches.append(MismatchDetail(
                    position=position,
                    ref_residue=ref_res,
                    other_residue=other_res,
                    mismatch_type='name_mismatch',
                    description=f"Residue mismatch: {ref_res.res_name} vs {other_res.res_name}"
                ))
        
        return mismatches
    
    def _check_atom_differences(
        self,
        ref: StructureSummary,
        other: StructureSummary,
        report: ValidationReport
    ) -> None:
        """Check for significant atom count differences (warnings only)."""
        ref_dict = {(r.chain_id, r.res_seq): r for r in ref.residues}
        other_dict = {(r.chain_id, r.res_seq): r for r in other.residues}
        
        for key in ref_dict:
            if key in other_dict:
                ref_res = ref_dict[key]
                other_res = other_dict[key]
                
                if ref_res.res_name != other_res.res_name:
                    continue  # Already reported as mismatch
                
                atom_diff = abs(ref_res.atom_count - other_res.atom_count)
                
                if atom_diff > 5:
                    report.warnings.append(
                        f"Residue {ref_res.res_name}{ref_res.res_seq} (chain {ref_res.chain_id}): "
                        f"{ref.name} has {ref_res.atom_count} atoms, "
                        f"{other.name} has {other_res.atom_count} atoms "
                        f"(difference: {atom_diff})"
                    )
    
    def _check_protonation_differences(
        self,
        ref: StructureSummary,
        other: StructureSummary,
        report: ValidationReport
    ) -> None:
        """Check for protonation state differences."""
        ref_dict = {(r.chain_id, r.res_seq): r for r in ref.residues}
        other_dict = {(r.chain_id, r.res_seq): r for r in other.residues}
        
        protonation_diffs = []
        
        for key in ref_dict:
            if key in other_dict:
                ref_res = ref_dict[key]
                other_res = other_dict[key]
                
                # Check if both are protonation variants
                ref_in_variants = ref_res.res_name in self.PROTONATION_VARIANTS
                other_in_variants = other_res.res_name in self.PROTONATION_VARIANTS
                
                if ref_in_variants and other_in_variants and ref_res.res_name != other_res.res_name:
                    protonation_diffs.append(
                        f"{ref_res.res_name}{ref_res.res_seq} (chain {ref_res.chain_id}): "
                        f"{ref.name}={ref_res.res_name}, {other.name}={other_res.res_name}"
                    )
        
        if protonation_diffs:
            report.warnings.append(
                f"Detected {len(protonation_diffs)} protonation state differences. "
                f"First few: {', '.join(protonation_diffs[:3])}"
            )


def validate_pdb_compatibility(
    pdb_files: list[str | Path],
    state_names: list[str] | None = None,
    verbose: bool = True
) -> None:
    """
    Convenience function for PDB compatibility validation.
    
    Validates that multiple PDB files are compatible for multiple-basin
topology generation. Raises ValueError with detailed message if incompatible.
    
    Args:
        pdb_files: List of PDB file paths
        state_names: Optional names for each structure (for error messages)
        verbose: Whether to print progress messages
        
    Raises:
        FileNotFoundError: If a PDB file doesn't exist
        ValueError: If structures are incompatible
    """
    if len(pdb_files) < 2:
        return
    
    if verbose:
        print("\n  [Pre-check] Validating PDB file compatibility...")
    
    validator = PDBCompatibilityValidator(verbose=verbose)
    
    try:
        reports = validator.validate(pdb_files, state_names)
        if verbose:
            print("    ✓ All PDB files are compatible")
            # Print summary of findings
            for report in reports:
                if report.warnings:
                    print(f"\n    ⚠️  Warnings for {report.other_name}:")
                    for warning in report.warnings[:3]:  # Show first 3
                        print(f"       - {warning[:70]}...")
    except ValueError:
        raise
