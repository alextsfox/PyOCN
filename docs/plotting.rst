Plotting Functions
==================

PyOCN provides several visualization functions for exploring and analyzing channel networks.

Overview
--------

The plotting module offers three main visualization approaches:

* **Raster plots**: Show the network as an array
* **Graph plots**: Display the network structure as a directed graph
* **Positional plots**: Combine spatial positioning with graph structure

Basic Examples
--------------

.. code-block:: python

   import PyOCN as po
   import matplotlib.pyplot as plt
   
   # Create and optimize a network
   ocn = po.OCN.from_net_type("I", dims=(32, 32), gamma=0.5)
   ocn.fit()
   
   # Raster visualization
   fig, ax = plt.subplots()
   po.plot_ocn_raster(ocn, ax=ax)
   plt.show()
   
   # Graph visualization
   fig, ax = plt.subplots()
   po.plot_ocn_as_dag(ocn, ax=ax)
   plt.show()

Function Reference
------------------

.. autofunction:: PyOCN.plot_ocn_raster

.. autofunction:: PyOCN.plot_ocn_as_dag

.. autofunction:: PyOCN.plot_positional_digraph