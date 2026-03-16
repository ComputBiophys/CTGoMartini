"""
End-to-end tests for vermouth >= 0.15.0 compatibility.

Tests the new Go contact workflow using martinize2 -go option:
    1. Convert rCSU .map files to contact_map.out format
    2. Run martinize2 with -go flag to generate Go contacts
    3. Verify go_nbparams.itp correctly maps to [contacts] in final topology

This validates the LJ to contacts conversion mechanism in CTGoMartini,
which converts nonbond_params Go interactions to bond-like contacts
for computational efficiency and accuracy.
"""

import os
import subprocess
import tempfile
import shutil

import pytest

import ctgomartini
from ctgomartini.data.ctgomartinize import SBGOMartinize
from ctgomartini.utils import convert_map_format


class TestVermouth015Contacts:
    """
    Test the Go contacts generation and conversion for vermouth >= 0.15.0
    """

    @pytest.fixture
    def test_data_dir(self):
        """Path to test data directory."""
        return os.path.join(
            os.path.dirname(__file__),
            "../fixtures/WriteItp/WriteItp/ref"
        )

    @pytest.fixture
    def temp_work_dir(self):
        """Create a temporary working directory."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir, ignore_errors=True)

    def parse_go_nbparams(self, filepath: str) -> set:
        """
        Parse go_nbparams.itp and return set of contacts.
        
        Args:
            filepath: Path to go_nbparams.itp file
            
        Returns:
            Set of ((atom1, atom2), sigma) tuples
        """
        contacts = set()
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('['):
                    continue
                parts = line.split()
                # Format: Open_X Open_Y 1 sigma epsilon ; comment
                if len(parts) >= 5 and parts[2] == '1':
                    atom1 = int(parts[0].replace('Open_', ''))
                    atom2 = int(parts[1].replace('Open_', ''))
                    sigma = float(parts[3])
                    # Store as sorted tuple to handle bidirectional
                    key = tuple(sorted([atom1, atom2]))
                    contacts.add((key, sigma))
        return contacts

    def parse_protein_contacts(self, filepath: str) -> set:
        """
        Parse Protein.itp [contacts] section and return set of contacts.
        
        Virtual site IDs (starting from 513) are converted back to 
        atom IDs by subtracting 512.
        
        Args:
            filepath: Path to Protein.itp file
            
        Returns:
            Set of ((atom1, atom2), sigma) tuples
        """
        contacts = set()
        in_contacts = False
        
        with open(filepath) as f:
            for line in f:
                line = line.strip()
                if line == '[ contacts ]':
                    in_contacts = True
                    continue
                if in_contacts and line.startswith('['):
                    break
                if in_contacts and line and line[0].isdigit():
                    parts = line.split()
                    # Format: atom1 atom2 1 sigma epsilon
                    if len(parts) >= 5 and parts[2] == '1':
                        atom1 = int(parts[0])
                        atom2 = int(parts[1])
                        sigma = float(parts[3])
                        # Convert virtual_site ID to atom ID
                        # Virtual sites start at 513 (512 offset from atoms)
                        atom1_real = atom1 - 512 if atom1 > 512 else atom1
                        atom2_real = atom2 - 512 if atom2 > 512 else atom2
                        key = tuple(sorted([atom1_real, atom2_real]))
                        contacts.add((key, sigma))
        return contacts

    def test_map_format_conversion(self, test_data_dir, temp_work_dir):
        """
        Test convert_map_format correctly converts .map to contact_map.out
        """
        # Copy input files
        pdb_file = os.path.join(test_data_dir, "1GGG_1_clean.pdb")
        map_file = os.path.join(test_data_dir, "1GGG_1_clean.map")
        
        if not os.path.exists(map_file):
            pytest.skip(f"Test data not found: {map_file}")
        
        work_pdb = os.path.join(temp_work_dir, "1GGG_1_clean.pdb")
        work_map = os.path.join(temp_work_dir, "1GGG_1_clean.map")
        shutil.copy(pdb_file, work_pdb)
        shutil.copy(map_file, work_map)
        
        # Convert map format
        output_file = os.path.join(temp_work_dir, "contact_map.out")
        result = convert_map_format(
            input_file=work_map,
            output_file=output_file,
            pdb_name="1GGG_1_clean.pdb",
            force=True
        )
        
        # Verify output file exists
        assert os.path.exists(result), "contact_map.out was not created"
        
        # Verify file has correct header format
        with open(result) as f:
            content = f.read()
            assert "CONTACT MAPS FROM PDB FILES" in content
            assert "Residue-Residue Contacts" in content
            assert "I1,I2    - serial residue id" in content

    def test_martinize2_go_generation(self, test_data_dir, temp_work_dir):
        """
        Test martinize2 -go generates correct go_nbparams.itp
        """
        # Copy input files
        pdb_file = os.path.join(test_data_dir, "1GGG_1_clean.pdb")
        map_file = os.path.join(test_data_dir, "1GGG_1_clean.map")
        
        if not os.path.exists(map_file):
            pytest.skip(f"Test data not found: {map_file}")
        
        work_pdb = os.path.join(temp_work_dir, "1GGG_1_clean.pdb")
        work_map = os.path.join(temp_work_dir, "1GGG_1_clean.map")
        shutil.copy(pdb_file, work_pdb)
        shutil.copy(map_file, work_map)
        
        # Convert map format first
        contact_map = os.path.join(temp_work_dir, "contact_map.out")
        convert_map_format(
            input_file=work_map,
            output_file=contact_map,
            pdb_name="1GGG_1_clean.pdb",
            force=True
        )
        
        # Run martinize2 with -go
        os.chdir(temp_work_dir)
        result = subprocess.run(
            f"martinize2 -f 1GGG_1_clean.pdb -x protein_cg.pdb -o system.top "
            f"-ff martini3001 -p backbone -dssp "
            f"-go contact_map.out -go-low 0.3 -go-up 1.1 -go-eps 12.0 "
            f"-name Open",
            shell=True,
            capture_output=True,
            encoding='utf-8'
        )
        
        assert result.returncode == 0, f"martinize2 failed: {result.stderr}"
        
        # Verify go files are generated
        assert os.path.exists(os.path.join(temp_work_dir, "go_atomtypes.itp"))
        assert os.path.exists(os.path.join(temp_work_dir, "go_nbparams.itp"))
        
        # Verify go_nbparams has correct format
        nbparams = self.parse_go_nbparams(
            os.path.join(temp_work_dir, "go_nbparams.itp")
        )
        assert len(nbparams) > 0, "go_nbparams.itp is empty"

    def test_sbgomartinize_contacts_conversion(self, test_data_dir, temp_work_dir):
        """
        Test SBGOMartinize correctly converts go_nbparams to [contacts]
        
        This is the core test validating that:
        1. martinize2 -go generates go_nbparams.itp with LJ parameters
        2. GenSBPTop reads the topology including go_nbparams
        3. Final Protein.itp has [contacts] section with converted interactions
        4. All contacts are correctly mapped (virtual_site ID -> atom ID)
        """
        # Copy input files
        pdb_file = os.path.join(test_data_dir, "1GGG_1_clean.pdb")
        map_file = os.path.join(test_data_dir, "1GGG_1_clean.map")
        
        if not os.path.exists(map_file):
            pytest.skip(f"Test data not found: {map_file}")
        
        work_pdb = os.path.join(temp_work_dir, "1GGG_1_clean.pdb")
        work_map = os.path.join(temp_work_dir, "1GGG_1_clean.map")
        shutil.copy(pdb_file, work_pdb)
        shutil.copy(map_file, work_map)
        
        # Copy force field
        ff_source = os.path.join(
            ctgomartini.__path__[0],
            "data/ForceFields/Martini300/martini_v3.0.0.itp"
        )
        ff_dest = os.path.join(temp_work_dir, "martini_v3.0.0.itp")
        shutil.copy(ff_source, ff_dest)
        
        # Run SBGOMartinize
        os.chdir(temp_work_dir)
        SBGOMartinize(
            aa_strfile_list=["1GGG_1_clean.pdb"],
            map_file_list=["1GGG_1_clean.map"],
            state_name_list=["Open"],
            sbmol_name="Protein",
            method="SBP",
            dssp=None,  # Use MDTraj
            ff="martini3001",
            other_params=""
        )
        
        # Parse go_nbparams.itp from Open directory
        go_nbparams_file = os.path.join(temp_work_dir, "Open", "go_nbparams.itp")
        assert os.path.exists(go_nbparams_file), "go_nbparams.itp not found"
        
        nb_contacts = self.parse_go_nbparams(go_nbparams_file)
        assert len(nb_contacts) > 0, "go_nbparams.itp has no contacts"
        
        # Parse [contacts] from Protein.itp
        protein_itp_file = os.path.join(temp_work_dir, "Protein.itp")
        assert os.path.exists(protein_itp_file), "Protein.itp not found"
        
        protein_contacts = self.parse_protein_contacts(protein_itp_file)
        assert len(protein_contacts) > 0, "Protein.itp has no [contacts]"
        
        # Verify all go_nbparams contacts are in Protein.itp
        missing_in_protein = nb_contacts - protein_contacts
        extra_in_protein = protein_contacts - nb_contacts
        
        assert len(missing_in_protein) == 0, (
            f"{len(missing_in_protein)} contacts missing in Protein.itp: "
            f"{list(missing_in_protein)[:5]}"
        )
        
        assert len(extra_in_protein) == 0, (
            f"{len(extra_in_protein)} extra contacts in Protein.itp: "
            f"{list(extra_in_protein)[:5]}"
        )
        
        # Verify counts match
        assert len(nb_contacts) == len(protein_contacts), (
            f"Contact count mismatch: go_nbparams={len(nb_contacts)}, "
            f"Protein.itp={len(protein_contacts)}"
        )

    def test_contact_sigma_values_match(self, test_data_dir, temp_work_dir):
        """
        Test that sigma values are preserved during contacts conversion
        """
        # Copy input files
        pdb_file = os.path.join(test_data_dir, "1GGG_1_clean.pdb")
        map_file = os.path.join(test_data_dir, "1GGG_1_clean.map")
        
        if not os.path.exists(map_file):
            pytest.skip(f"Test data not found: {map_file}")
        
        work_pdb = os.path.join(temp_work_dir, "1GGG_1_clean.pdb")
        work_map = os.path.join(temp_work_dir, "1GGG_1_clean.map")
        shutil.copy(pdb_file, work_pdb)
        shutil.copy(map_file, work_map)
        
        # Copy force field
        ff_source = os.path.join(
            ctgomartini.__path__[0],
            "data/ForceFields/Martini300/martini_v3.0.0.itp"
        )
        ff_dest = os.path.join(temp_work_dir, "martini_v3.0.0.itp")
        shutil.copy(ff_source, ff_dest)
        
        # Run SBGOMartinize
        os.chdir(temp_work_dir)
        SBGOMartinize(
            aa_strfile_list=["1GGG_1_clean.pdb"],
            map_file_list=["1GGG_1_clean.map"],
            state_name_list=["Open"],
            sbmol_name="Protein",
            method="SBP",
            dssp=None,
            ff="martini3001",
            other_params=""
        )
        
        # Parse both files
        go_nbparams_file = os.path.join(temp_work_dir, "Open", "go_nbparams.itp")
        protein_itp_file = os.path.join(temp_work_dir, "Protein.itp")
        
        nb_contacts = self.parse_go_nbparams(go_nbparams_file)
        protein_contacts = self.parse_protein_contacts(protein_itp_file)
        
        # Find common contacts and compare sigma values
        nb_dict = {k: s for k, s in nb_contacts}
        protein_dict = {k: s for k, s in protein_contacts}
        
        sigma_mismatches = []
        for key in nb_dict:
            if key in protein_dict:
                nb_sigma = nb_dict[key]
                protein_sigma = protein_dict[key]
                if abs(nb_sigma - protein_sigma) > 1e-6:
                    sigma_mismatches.append((key, nb_sigma, protein_sigma))
        
        assert len(sigma_mismatches) == 0, (
            f"{len(sigma_mismatches)} sigma value mismatches: "
            f"{sigma_mismatches[:5]}"
        )

    def test_virtual_site_id_mapping(self, test_data_dir, temp_work_dir):
        """
        Test virtual site ID to atom ID mapping (512 offset)
        
        In vermouth >= 0.15.0:
        - Atoms are numbered 1-512 (for 220 residues with Martini 3)
        - Virtual sites for Go contacts start at 513
        - The mapping is: virtual_site_id = atom_id + 512
        """
        # Copy input files
        pdb_file = os.path.join(test_data_dir, "1GGG_1_clean.pdb")
        map_file = os.path.join(test_data_dir, "1GGG_1_clean.map")
        
        if not os.path.exists(map_file):
            pytest.skip(f"Test data not found: {map_file}")
        
        work_pdb = os.path.join(temp_work_dir, "1GGG_1_clean.pdb")
        work_map = os.path.join(temp_work_dir, "1GGG_1_clean.map")
        shutil.copy(pdb_file, work_pdb)
        shutil.copy(map_file, work_map)
        
        # Copy force field
        ff_source = os.path.join(
            ctgomartini.__path__[0],
            "data/ForceFields/Martini300/martini_v3.0.0.itp"
        )
        ff_dest = os.path.join(temp_work_dir, "martini_v3.0.0.itp")
        shutil.copy(ff_source, ff_dest)
        
        # Run SBGOMartinize
        os.chdir(temp_work_dir)
        SBGOMartinize(
            aa_strfile_list=["1GGG_1_clean.pdb"],
            map_file_list=["1GGG_1_clean.map"],
            state_name_list=["Open"],
            sbmol_name="Protein",
            method="SBP",
            dssp=None,
            ff="martini3001",
            other_params=""
        )
        
        # Check that virtual site IDs in contacts are in correct range
        protein_itp_file = os.path.join(temp_work_dir, "Protein.itp")
        
        with open(protein_itp_file) as f:
            content = f.read()
            in_contacts = False
            for line in content.split('\n'):
                line = line.strip()
                if line == '[ contacts ]':
                    in_contacts = True
                    continue
                if in_contacts and line.startswith('['):
                    break
                if in_contacts and line and line[0].isdigit():
                    parts = line.split()
                    if len(parts) >= 2:
                        atom1 = int(parts[0])
                        atom2 = int(parts[1])
                        # Both should be >= 1, virtual sites >= 513
                        assert atom1 >= 1, f"Invalid atom ID: {atom1}"
                        assert atom2 >= 1, f"Invalid atom ID: {atom2}"
                        # At least one should be a virtual site (>= 513)
                        assert atom1 >= 513 or atom2 >= 513, (
                            f"No virtual site in contact: {atom1}-{atom2}"
                        )
