"""
Basic test suite for PyOCN.

These tests verify core functionality with deterministic outputs.
"""
import unittest
import numpy as np
import networkx as nx
import PyOCN as po

class TestBasicOCN(unittest.TestCase):
    """Basic tests for OCN creation and energy computation."""
    
    def test_ocn_energy(self):
        """Test that OCN creation produces expected energy value."""
        ocn = po.OCN.from_net_type(
            net_type="V",
            dims=(64, 64),
            random_state=8472,
        )
        
        expected_energy = 16469.684
        actual_energy = ocn.energy
        
        self.assertAlmostEqual(
            actual_energy, 
            expected_energy, 
            places=3,
            msg=f"Expected energy {expected_energy}, got {actual_energy}"
        )
    
    def test_ocn_copy(self):
        """Test that copying an OCN instance produces an identical copy."""
        ocn = po.OCN.from_net_type(
            net_type="V",
            dims=(32, 32),
            random_state=1234,
        )
        
        ocn_copy = ocn.copy()
        
        # Check that the copy has the same attributes
        self.assertEqual(ocn.dims, ocn_copy.dims, "Dimensions do not match.")
        self.assertEqual(ocn.energy, ocn_copy.energy, "Energies do not match.")
        self.assertEqual(ocn.nroots, ocn_copy.nroots, "Number of roots do not match.")
        self.assertEqual(ocn.wrap, ocn_copy.wrap, "Wrap settings do not match.")
        
        # Check that the internal data arrays are equal
        original_array = ocn.to_numpy()
        copy_array = ocn_copy.to_numpy()
        
        np.testing.assert_array_equal(
            original_array, 
            copy_array, 
            err_msg="Internal data arrays do not match."
        )

        # fit both and check energies are still equal
        ocn.fit(array_reports=0)
        ocn_copy.fit(array_reports=0)
        
        self.assertEqual(ocn.energy, ocn_copy.energy, "Energies do not match after fitting.")

    def test_ocn_custom_cooling(self):
        """Test that custom cooling schedule affects energy as expected."""
        ocn = po.OCN.from_net_type(
            net_type="V",
            dims=(64, 64),
            random_state=143798,
        )
        
        # Fit with a custom cooling schedule
        energy = ocn.energy
        def custom_schedule(iter):
            return energy / (iter + 1)
        
        ocn_copy = ocn.copy()

        ocn_copy.fit_custom_cooling(custom_schedule, max_iterations_per_loop=100, n_iterations=50_000, array_reports=0, iteration_start=0)
        ocn.fit_custom_cooling(custom_schedule, max_iterations_per_loop=100, n_iterations=5_000, array_reports=0, iteration_start=0)
        ocn.fit_custom_cooling(custom_schedule, max_iterations_per_loop=100, n_iterations=45_000, array_reports=0, iteration_start=5_000)

        self.assertAlmostEqual(ocn.energy, ocn_copy.energy, places=6)

    def test_ocn_rng(self):
        """Test that the rng property returns the expected random state."""
        ocn = po.OCN.from_net_type(
            net_type="V",
            dims=(64, 64),
            random_state=8472,
        )
        expected_rng = (510574073, 2087720647, 3836914231, 3781483648)
        
        self.assertEqual(
            ocn.rng, 
            expected_rng, 
            f"Expected RNG state {expected_rng}, got {ocn.rng}"
        )

        ocn.rng = 8470
        self.assertNotEqual(
            ocn.rng, 
            expected_rng, 
            f"After changing seed, RNG state should differ from {expected_rng}"
        )

        # Reset to original seed and check
        ocn.rng = 8472
        self.assertEqual(
            ocn.rng, 
            expected_rng, 
            f"After resetting, expected RNG state {expected_rng}, got {ocn.rng}"
        )

        ocn.rng = 8472
        ocn.rng = 8472
        ocn.rng = 8472
        self.assertEqual(
            ocn.rng, 
            expected_rng, 
            f"Setting rng to the same seed multiple times should not change the state. Expected {expected_rng}, got {ocn.rng}"
        )

        ocn.single_iteration(0)
        self.assertNotEqual(
            ocn.rng, 
            expected_rng, 
            f"After iteration, RNG state should differ from {expected_rng}. Got {ocn.rng}"
        )




    


if __name__ == "__main__":
    unittest.main()