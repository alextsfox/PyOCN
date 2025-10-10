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

    def test_from_digraph(self):
        """Test creating OCN from custom NetworkX digraph."""
        # Create a simple 3x3 cross pattern like in demo
        dag = nx.DiGraph()
        for i in range(9):
            row, col = divmod(i, 3)
            dag.add_node(i, pos=(row, col))
        
        # Simple flow pattern: all flow to center (node 4)
        for i in [0, 1, 2, 3, 5, 6, 7, 8]:
            dag.add_edge(i, 4)
        
        ocn = po.OCN.from_digraph(dag, random_state=1234)
        
        self.assertEqual(ocn.dims, (3, 3))
        self.assertEqual(ocn.nroots, 1)  # Only center node should be root
        self.assertAlmostEqual(ocn.energy, 11.0, places=6)
        
        ocn = po.OCN.from_digraph(dag, random_state=1234, gamma=1)
        self.assertAlmostEqual(ocn.energy, 17.0, places=6)

    def test_export_formats(self):
        """Test different export formats produce consistent data."""
        ocn = po.OCN.from_net_type("V", dims=(10, 10), random_state=9999)
        
        # Test numpy export
        numpy_array = ocn.to_numpy(unwrap=False)
        self.assertEqual(numpy_array.shape, (3, 10, 10))  # 3 channels: energy, area, watershed
        
        # Test digraph export  
        dag = ocn.to_digraph()
        self.assertEqual(len(dag.nodes), 100)  # 10x10 = 100 nodes
        self.assertTrue(nx.is_directed_acyclic_graph(dag))
        
        # Check that node attributes exist
        for node in dag.nodes:
            attrs = dag.nodes[node]
            self.assertIn('pos', attrs)
            self.assertIn('drained_area', attrs)
            self.assertIn('energy', attrs)
            self.assertIn('watershed_id', attrs)

    def test_periodic_boundaries(self):
        """Test OCN with periodic boundary conditions."""
        ocn_wrap = po.OCN.from_net_type("H", dims=(10, 10), wrap=True, random_state=5555)
        ocn_no_wrap = po.OCN.from_net_type("H", dims=(10, 10), wrap=False, random_state=5555)
        ocn_wrap.fit(n_iterations=500, pbar=False)
        ocn_no_wrap.fit(n_iterations=500, pbar=False)
        
        self.assertTrue(ocn_wrap.wrap)
        self.assertFalse(ocn_no_wrap.wrap)
        
        # Both should have same dimensions but potentially different energies
        self.assertEqual(ocn_wrap.dims, ocn_no_wrap.dims)
        # Wrapped version typically has different energy due to edge connections
        self.assertNotEqual(ocn_wrap.energy, ocn_no_wrap.energy)

    def test_fit_convergence(self):
        """Test that early exit works as intended."""
        ocn = po.OCN.from_net_type("E", dims=(20, 20), random_state=7777)
        
        ocn.fit(n_iterations=20*20*500, pbar=False, tol=1e-4)
        final_energy = ocn.energy

        self.assertIsInstance(final_energy, float)
        self.assertGreater(len(ocn.history), 0)
        
        # History should have correct structure: [iteration, energy, temperature]
        self.assertEqual(ocn.history.shape[1], 3)
        self.assertTrue(np.all(ocn.history[:, 0] >= 0))  # iterations >= 0
        self.assertTrue(np.all(ocn.history[:, 1] > 0))   # energy > 0
        self.assertTrue(np.all(ocn.history[:, 2] >= 0))  # temperature >= 0
        
        # final energy should match last history entry
        self.assertEqual(ocn.energy, ocn.history[-1, 1])

        # iteration should have stopped before max iterations and should be less than tol
        self.assertLess(ocn.history[-1, 0], 20*20*500)
        self.assertLessEqual((ocn.history[-1, 1] - ocn.history[-2, 1])/ocn.history[-1, 1], 1e-4)

    def test_single_iteration(self):
        """Test single iteration method."""
        ocn = po.OCN.from_net_type("V", dims=(16, 16), random_state=3333)
        initial_history_len = len(ocn.history)
        
        result = ocn.single_iteration(temperature=ocn.energy, array_report=False)
        
        # Should return None when array_report=False
        self.assertIsNone(result)
        
        # History should have one more entry
        self.assertEqual(len(ocn.history), initial_history_len + 1)

    def test_error_handling(self):
        """Test error handling for invalid inputs."""
        # Invalid gamma
        with self.assertRaises(TypeError):
            po.OCN.from_net_type("V", dims=(10, 10), gamma="invalid")
        
        # Invalid net_type
        with self.assertRaises(ValueError):
            po.OCN.from_net_type("INVALID", dims=(10, 10))
        
        # Invalid dimensions
        with self.assertRaises((ValueError, TypeError)):
            po.OCN.from_net_type("V", dims=(5,))  # Should need 2D dims


if __name__ == "__main__":
    unittest.main()