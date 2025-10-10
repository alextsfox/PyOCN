OCN Class
=========

The :class:`OCN` class is the main interface for creating and optimizing channel networks.

Basic Usage
-----------

.. code-block:: python

   import PyOCN as po
   
   # Create from a network type
   ocn = po.OCN.from_net_type("I", dims=(64, 64), gamma=0.5)
   
   # Optimize the network
   ocn.fit()
   
   # Access properties
   print(f"Final energy: {ocn.energy}")

Class Reference
---------------

.. autoclass:: PyOCN.OCN
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__