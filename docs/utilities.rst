Utility Functions
=================

PyOCN includes several utility functions for working with channel networks, optimization schedules, and network analysis.

Network Generation
------------------

.. autofunction:: PyOCN.utils.net_type_to_dag

Optimization
------------

.. autofunction:: PyOCN.utils.simulated_annealing_schedule

Example usage:

.. code-block:: python

   import PyOCN as po
   
   ocn = po.OCN.from_net_type("I", dims=(64, 64))
   n_iterations = 50000
   # Create a custom cooling schedule
   schedule = po.utils.simulated_annealing_schedule(
       dims=(64, 64),
       E0=ocn.energy,
       constant_phase=0.1,
       n_iterations=n_iterations,
       cooling_rate=1.0
   )
   
   # Use in optimization
   ocn.fit_custom_cooling(cooling_func=schedule, n_iterations=n_iterations)

Network Analysis
----------------

.. autofunction:: PyOCN.utils.unwrap_digraph

.. autofunction:: PyOCN.utils.assign_subwatersheds

.. autofunction:: PyOCN.utils.get_subwatersheds

Example watershed analysis:

.. code-block:: python

   import PyOCN as po
   
   # Create and optimize network
   ocn = po.OCN.from_net_type("I", dims=(32, 32))
   ocn.fit(n_iterations=10000)
   
   # Get the network as a graph
   dag = ocn.to_dag()
   
   # Analyze subwatersheds
   subwatersheds = po.utils.get_subwatersheds(dag, node=5)
   print(f"Found {len(subwatersheds)} outlet nodes")