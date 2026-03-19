"""Tests for new contact types: Contacts6_12 and Contacts10_12.

This module tests:
1. Contacts6_12 (functype=2) - 6-12 LJ with adjustable rcut
2. Contacts10_12 (functype=3) - 10-12 LJ with adjustable rcut
3. Corresponding multi-state handlers
"""

import unittest
import numpy as np
import openmm as mm
import openmm.unit as unit

from ctgomartini.topology.interactions.bonds import ContactsLJ, Contacts6_12, Contacts10_12
from ctgomartini.topology.interactions.multi_state import (
    MultiContactsLJ, MultiContacts6_12, MultiContacts10_12
)


class TestContacts6_12(unittest.TestCase):
    """Test Contacts6_12 class."""

    def test_initialization(self):
        """Test that the class initializes correctly."""
        contacts = Contacts6_12(use_sigma_eps=True)
        self.assertEqual(contacts.name, 'contacts_6_12_adjustable')
        self.assertEqual(contacts.category, 'contacts')
        self.assertEqual(contacts.FUNCTYPE, "2")
        self.assertEqual(contacts.EXPECTED_FIELDS, 6)

    def test_mm_force_parameters(self):
        """Test that mm_force has correct per-bond parameters."""
        contacts = Contacts6_12()
        # Should have C12, C6, rcut as per-bond parameters
        self.assertEqual(contacts.mm_force.getNumPerBondParameters(), 3)

    def test_add_interaction_with_c6_c12(self):
        """Test adding interaction with C6/C12 values."""
        contacts = Contacts6_12(use_sigma_eps=False)
        # fields: [ai, aj, functype, C6, C12, rcut]
        fields = ["1", "2", "2", "1.0", "2.0", "1.2"]
        contacts.add_interaction(fields, base_atom_index=0, offset=-1)

        self.assertEqual(contacts.mm_force.getNumBonds(), 1)
        # Check parameters: [C12, C6, rcut]
        _, _, params = contacts.mm_force.getBondParameters(0)
        self.assertAlmostEqual(params[0], 2.0)  # C12
        self.assertAlmostEqual(params[1], 1.0)  # C6
        self.assertAlmostEqual(params[2], 1.2)  # rcut

    def test_add_interaction_with_sigma_eps(self):
        """Test adding interaction with sigma/epsilon values."""
        contacts = Contacts6_12(use_sigma_eps=True)
        # fields: [ai, aj, functype, sigma, epsilon, rcut]
        sigma = 0.5
        eps = 1.0
        fields = ["1", "2", "2", str(sigma), str(eps), "1.1"]
        contacts.add_interaction(fields, base_atom_index=0, offset=-1)

        # Check calculated C6 and C12
        expected_C6 = 4 * eps * sigma ** 6
        expected_C12 = 4 * eps * sigma ** 12

        _, _, params = contacts.mm_force.getBondParameters(0)
        self.assertAlmostEqual(params[0], expected_C12, places=10)
        self.assertAlmostEqual(params[1], expected_C6, places=10)
        self.assertAlmostEqual(params[2], 1.1)


class TestContacts10_12(unittest.TestCase):
    """Test Contacts10_12 class."""

    def test_initialization(self):
        """Test that the class initializes correctly."""
        contacts = Contacts10_12(use_sigma_eps=True)
        self.assertEqual(contacts.name, 'contacts_10_12')
        self.assertEqual(contacts.category, 'contacts')
        self.assertEqual(contacts.FUNCTYPE, "3")
        self.assertEqual(contacts.EXPECTED_FIELDS, 6)

    def test_mm_force_parameters(self):
        """Test that mm_force has correct per-bond parameters."""
        contacts = Contacts10_12()
        # Should have C12, C10, rcut as per-bond parameters
        self.assertEqual(contacts.mm_force.getNumPerBondParameters(), 3)

    def test_add_interaction_with_c10_c12(self):
        """Test adding interaction with C10/C12 values."""
        contacts = Contacts10_12(use_sigma_eps=False)
        # fields: [ai, aj, functype, C10, C12, rcut]
        fields = ["1", "2", "3", "1.0", "2.0", "1.5"]
        contacts.add_interaction(fields, base_atom_index=0, offset=-1)

        self.assertEqual(contacts.mm_force.getNumBonds(), 1)
        # Check parameters: [C12, C10, rcut]
        _, _, params = contacts.mm_force.getBondParameters(0)
        self.assertAlmostEqual(params[0], 2.0)  # C12
        self.assertAlmostEqual(params[1], 1.0)  # C10
        self.assertAlmostEqual(params[2], 1.5)  # rcut

    def test_add_interaction_with_sigma_eps(self):
        """Test adding interaction with sigma/epsilon values."""
        contacts = Contacts10_12(use_sigma_eps=True)
        # fields: [ai, aj, functype, sigma, epsilon, rcut]
        sigma = 0.5
        eps = 1.0
        fields = ["1", "2", "3", str(sigma), str(eps), "1.0"]
        contacts.add_interaction(fields, base_atom_index=0, offset=-1)

        # Check calculated C10 and C12
        expected_C10 = 4 * eps * sigma ** 10
        expected_C12 = 4 * eps * sigma ** 12

        _, _, params = contacts.mm_force.getBondParameters(0)
        self.assertAlmostEqual(params[0], expected_C12, places=10)
        self.assertAlmostEqual(params[1], expected_C10, places=10)
        self.assertAlmostEqual(params[2], 1.0)


class TestMultiContacts6_12(unittest.TestCase):
    """Test MultiContacts6_12 handler."""

    def test_class_attributes(self):
        """Test class-level attributes."""
        self.assertEqual(MultiContacts6_12.name, "contact_lj_adj")
        self.assertEqual(MultiContacts6_12.category, "multi_contacts")
        self.assertEqual(MultiContacts6_12.functype, "2")
        self.assertEqual(MultiContacts6_12.expected_fields, 8)

    def test_can_handle(self):
        """Test can_handle method."""
        # fields: [ai, aj, nstates, stateid, functype, C6, C12, rcut]
        fields = ["1", "2", "2", "1", "2", "0.5", "1.0", "1.2"]
        self.assertTrue(MultiContacts6_12.can_handle("multi_contacts", fields))

        # Wrong functype
        fields_wrong = ["1", "2", "2", "1", "1", "0.5", "1.0", "1.2"]
        self.assertFalse(MultiContacts6_12.can_handle("multi_contacts", fields_wrong))

    def test_energy_expr(self):
        """Test energy expression."""
        handler = MultiContacts6_12()
        expr = handler.energy_expr
        self.assertIn("rcut", expr)  # rcut should be in expression (parameter)
        self.assertIn("C12", expr)
        self.assertIn("C6", expr)

    def test_per_bond_params(self):
        """Test per-bond parameters."""
        handler = MultiContacts6_12()
        params = handler.per_bond_params
        self.assertEqual(params, ["delta_contact_adj", "C12", "C6", "rcut"])

    def test_extract_params(self):
        """Test parameter extraction."""
        # With C6/C12
        fields = ["1", "2", "2", "1", "2", "1.0", "2.0", "1.2"]
        params = MultiContacts6_12.extract_params(fields, use_sigma_eps=False)

        self.assertAlmostEqual(params["delta_contact_adj"], 1.0)
        self.assertAlmostEqual(params["C6"], 1.0)
        self.assertAlmostEqual(params["C12"], 2.0)
        self.assertAlmostEqual(params["rcut"], 1.2)

        # With sigma/epsilon
        sigma = 0.5
        eps = 1.0
        fields = ["1", "2", "2", "1", "2", str(sigma), str(eps), "1.1"]
        params = MultiContacts6_12.extract_params(fields, use_sigma_eps=True)

        expected_C6 = 4 * eps * sigma ** 6
        expected_C12 = 4 * eps * sigma ** 12
        self.assertAlmostEqual(params["C6"], expected_C6)
        self.assertAlmostEqual(params["C12"], expected_C12)


class TestMultiContacts10_12(unittest.TestCase):
    """Test MultiContacts10_12 handler."""

    def test_class_attributes(self):
        """Test class-level attributes."""
        self.assertEqual(MultiContacts10_12.name, "contact_1012")
        self.assertEqual(MultiContacts10_12.category, "multi_contacts")
        self.assertEqual(MultiContacts10_12.functype, "3")
        self.assertEqual(MultiContacts10_12.expected_fields, 8)

    def test_can_handle(self):
        """Test can_handle method."""
        # fields: [ai, aj, nstates, stateid, functype, C10, C12, rcut]
        fields = ["1", "2", "2", "1", "3", "0.5", "1.0", "1.2"]
        self.assertTrue(MultiContacts10_12.can_handle("multi_contacts", fields))

        # Wrong functype
        fields_wrong = ["1", "2", "2", "1", "2", "0.5", "1.0", "1.2"]
        self.assertFalse(MultiContacts10_12.can_handle("multi_contacts", fields_wrong))

    def test_energy_expr(self):
        """Test energy expression."""
        handler = MultiContacts10_12()
        expr = handler.energy_expr
        self.assertIn("rcut", expr)
        self.assertIn("C12", expr)
        self.assertIn("C10", expr)  # Note: C10, not C6
        self.assertNotIn("C6", expr)

    def test_per_bond_params(self):
        """Test per-bond parameters."""
        handler = MultiContacts10_12()
        params = handler.per_bond_params
        self.assertEqual(params, ["delta_contact_1012", "C12", "C10", "rcut"])

    def test_extract_params(self):
        """Test parameter extraction."""
        # With C10/C12
        fields = ["1", "2", "2", "1", "3", "1.0", "2.0", "1.3"]
        params = MultiContacts10_12.extract_params(fields, use_sigma_eps=False)

        self.assertAlmostEqual(params["delta_contact_1012"], 1.0)
        self.assertAlmostEqual(params["C10"], 1.0)
        self.assertAlmostEqual(params["C12"], 2.0)
        self.assertAlmostEqual(params["rcut"], 1.3)

        # With sigma/epsilon
        sigma = 0.5
        eps = 1.0
        fields = ["1", "2", "2", "1", "3", str(sigma), str(eps), "1.1"]
        params = MultiContacts10_12.extract_params(fields, use_sigma_eps=True)

        expected_C10 = 4 * eps * sigma ** 10
        expected_C12 = 4 * eps * sigma ** 12
        self.assertAlmostEqual(params["C10"], expected_C10)
        self.assertAlmostEqual(params["C12"], expected_C12)


class TestContactsComparison(unittest.TestCase):
    """Compare different contact types."""

    def test_different_functypes(self):
        """Verify all three contact types have different functypes."""
        # ContactsLJ uses type_label=[2, "1"] not FUNCTYPE
        self.assertEqual(ContactsLJ().type_label, [2, "1"])
        self.assertEqual(Contacts6_12.FUNCTYPE, "2")
        self.assertEqual(Contacts10_12.FUNCTYPE, "3")

    def test_different_names(self):
        """Verify all contact types have different names."""
        names = [
            ContactsLJ().name,
            Contacts6_12().name,
            Contacts10_12().name,
        ]
        self.assertEqual(len(names), len(set(names)))

    def test_rcut_behavior(self):
        """Test rcut behavior differences."""
        # Original ContactsLJ: fixed rcut in expression
        c1 = ContactsLJ()
        expr1 = c1.mm_force.getEnergyFunction()
        self.assertIn("rcut=1.1", expr1)  # Fixed value

        # Adjustable Contacts: rcut as parameter
        c2 = Contacts6_12()
        expr2 = c2.mm_force.getEnergyFunction()
        self.assertNotIn("rcut=", expr2)  # No fixed value, it's a parameter

        c3 = Contacts10_12()
        expr3 = c3.mm_force.getEnergyFunction()
        self.assertNotIn("rcut=", expr3)


class TestContactsEnergyCalculation(unittest.TestCase):
    """Test energy calculation values match theoretical expectations."""

    def _calculate_energy(self, contacts, positions, box_vectors=None):
        """Helper to calculate energy using OpenMM."""
        system = mm.System()
        system.addParticle(1.0)
        system.addParticle(1.0)
        system.addForce(contacts.mm_force)

        if box_vectors is not None:
            system.setDefaultPeriodicBoxVectors(*box_vectors)

        integrator = mm.VerletIntegrator(0.001)
        platform = mm.Platform.getPlatformByName('Reference')
        context = mm.Context(system, integrator, platform)
        context.setPositions(positions)

        state = context.getState(getEnergy=True)
        return state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

    def test_contacts_6_12_energy_at_half_rcut(self):
        """Test 6-12 LJ energy at r = 0.5 * rcut."""
        contacts = Contacts6_12(use_sigma_eps=False)
        C12, C6, rcut = 1.0, 1.0, 1.1
        contacts.mm_force.addBond(0, 1, [C12, C6, rcut])

        # Test at r = 0.5 nm
        r = 0.5
        positions = [mm.Vec3(0, 0, 0), mm.Vec3(r, 0, 0)]
        energy = self._calculate_energy(contacts, positions)

        # Manual calculation: V = (C12/r^12 - C6/r^6) - (C12/rcut^12 - C6/rcut^6)
        expected = (C12/r**12 - C6/r**6) - (C12/rcut**12 - C6/rcut**6)
        self.assertAlmostEqual(energy, expected, places=5,
                               msg=f"6-12 energy at r={r}: expected {expected}, got {energy}")

    def test_contacts_6_12_energy_at_rcut(self):
        """Test 6-12 LJ energy is zero at r = rcut."""
        contacts = Contacts6_12(use_sigma_eps=False)
        C12, C6, rcut = 1.0, 1.0, 1.1
        contacts.mm_force.addBond(0, 1, [C12, C6, rcut])

        # Test at r = rcut
        r = rcut
        positions = [mm.Vec3(0, 0, 0), mm.Vec3(r, 0, 0)]
        energy = self._calculate_energy(contacts, positions)

        # At r = rcut, V should be exactly 0 (due to shifting)
        self.assertAlmostEqual(energy, 0.0, places=10,
                               msg=f"6-12 energy at r=rcut should be 0, got {energy}")

    def test_contacts_6_12_energy_beyond_rcut(self):
        """Test 6-12 LJ energy is zero when r > rcut."""
        contacts = Contacts6_12(use_sigma_eps=False)
        C12, C6, rcut = 1.0, 1.0, 1.1
        contacts.mm_force.addBond(0, 1, [C12, C6, rcut])

        # Test at r = 1.5 * rcut
        r = 1.5 * rcut
        positions = [mm.Vec3(0, 0, 0), mm.Vec3(r, 0, 0)]
        energy = self._calculate_energy(contacts, positions)

        # Beyond rcut, step function should give V = 0
        self.assertAlmostEqual(energy, 0.0, places=10,
                               msg=f"6-12 energy beyond rcut should be 0, got {energy}")

    def test_contacts_10_12_energy_at_half_rcut(self):
        """Test 10-12 LJ energy at r = 0.5 * rcut."""
        contacts = Contacts10_12(use_sigma_eps=False)
        C12, C10, rcut = 1.0, 1.0, 1.1
        contacts.mm_force.addBond(0, 1, [C12, C10, rcut])

        # Test at r = 0.5 nm
        r = 0.5
        positions = [mm.Vec3(0, 0, 0), mm.Vec3(r, 0, 0)]
        energy = self._calculate_energy(contacts, positions)

        # Manual calculation: V = (C12/r^12 - C10/r^10) - (C12/rcut^12 - C10/rcut^10)
        expected = (C12/r**12 - C10/r**10) - (C12/rcut**12 - C10/rcut**10)
        self.assertAlmostEqual(energy, expected, places=5,
                               msg=f"10-12 energy at r={r}: expected {expected}, got {energy}")

    def test_contacts_10_12_energy_at_rcut(self):
        """Test 10-12 LJ energy is zero at r = rcut."""
        contacts = Contacts10_12(use_sigma_eps=False)
        C12, C10, rcut = 1.0, 1.0, 1.1
        contacts.mm_force.addBond(0, 1, [C12, C10, rcut])

        # Test at r = rcut
        r = rcut
        positions = [mm.Vec3(0, 0, 0), mm.Vec3(r, 0, 0)]
        energy = self._calculate_energy(contacts, positions)

        # At r = rcut, V should be exactly 0 (due to shifting)
        self.assertAlmostEqual(energy, 0.0, places=10,
                               msg=f"10-12 energy at r=rcut should be 0, got {energy}")

    def test_contacts_10_12_energy_beyond_rcut(self):
        """Test 10-12 LJ energy is zero when r > rcut."""
        contacts = Contacts10_12(use_sigma_eps=False)
        C12, C10, rcut = 1.0, 1.0, 1.1
        contacts.mm_force.addBond(0, 1, [C12, C10, rcut])

        # Test at r = 1.5 * rcut
        r = 1.5 * rcut
        positions = [mm.Vec3(0, 0, 0), mm.Vec3(r, 0, 0)]
        energy = self._calculate_energy(contacts, positions)

        # Beyond rcut, step function should give V = 0
        self.assertAlmostEqual(energy, 0.0, places=10,
                               msg=f"10-12 energy beyond rcut should be 0, got {energy}")

    def test_contacts_lj_energy_at_rcut(self):
        """Test original ContactsLJ energy is zero at r = rcut (1.1 nm)."""
        contacts = ContactsLJ(use_sigma_eps=False)
        C12, C6 = 1.0, 1.0
        contacts.mm_force.addBond(0, 1, [C12, C6])

        # Test at r = 1.1 nm (fixed rcut)
        r = 1.1
        positions = [mm.Vec3(0, 0, 0), mm.Vec3(r, 0, 0)]
        energy = self._calculate_energy(contacts, positions)

        # At r = rcut, V should be exactly 0
        self.assertAlmostEqual(energy, 0.0, places=10,
                               msg=f"ContactsLJ energy at r=rcut should be 0, got {energy}")

    def test_contacts_lj_energy_beyond_rcut(self):
        """Test original ContactsLJ energy is zero when r > rcut."""
        contacts = ContactsLJ(use_sigma_eps=False)
        C12, C6 = 1.0, 1.0
        contacts.mm_force.addBond(0, 1, [C12, C6])

        # Test at r = 1.5 nm
        r = 1.5
        positions = [mm.Vec3(0, 0, 0), mm.Vec3(r, 0, 0)]
        energy = self._calculate_energy(contacts, positions)

        # Beyond rcut, energy should be 0
        self.assertAlmostEqual(energy, 0.0, places=10,
                               msg=f"ContactsLJ energy beyond rcut should be 0, got {energy}")


if __name__ == '__main__':
    unittest.main()
