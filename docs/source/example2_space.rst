
Example 2: Geometry-related features
------------------------------------

The code snippet below shows a simple scenario for geometry analysis of the reconstructed cytoskelton.

.. code-block:: python
   :linenos:

   from cytoskeleton_analyser.inout import Paths      # simulation-specific paths
   from cytoskeleton_analyser.position.drivers import process  # driver fuction
   from sim_specs import cell, data_in, data_out, rinds

   # Set general parameters of the input simulation. (a shortlist
   # of config settings).
   params = {
      'edge_len_fine': 0.008,    # (μm) length of filament edge.
      'iscoarse': True,          # use coarse-grained representation
      'coarsegraining': 25,      # coarse-graining factor
      'use_final': False,        # if True, limit to final snapshot only
      'cell': cell,              # cell object
      'slice_limits': {'bottom': 0., 'top': 0.8},  # z-pos. of slicing planes
   }

   # Specify features to process. If this is None,
   # all features are included.
   features = [
      'RadialMass',
   ]

   # Decide if this is a new analysis, or we merely want to import one.
   new = True

   # This dict will become populated only if analysis
   # is imported ('new' = False).
   # Otherwise, the results are stored on a disc.
   result = {}

   if __name__ == '__main__':

      for ci in rinds:
         paths = Paths(data_in, data_out, params['cell'], ci)
         result[str(ci)] = process(new, paths, ci, params, True, features)



