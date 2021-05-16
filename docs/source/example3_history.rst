
Example 3: Examination of microtubule dynamic instability
---------------------------------------------------------

The code snippet below shows a scenario for analysis of the microtubule dynamics where only *duration*,
*spatial_maps* and *event_collections* features are selected.


.. code-block:: python
   :linenos:

   from cytoskeleton_analyser.inout import Paths      # simulation-specific paths
   from cytoskeleton_analyser.history.drivers import process  # driver fuction
   from sim_specs import cell, data_in, data_out, rinds

   params = {
      'edge_len': 0.008,     # (μm) length of filament edge.
      'cell': cell,          # cell object
   }

   # Specify features to process. If this is None,
   # all features are included.
   features = [
      'duration',
      'spatial_maps',
      'event_collections',
   ]

   if __name__ == '__main__':

      for ci in rinds:
         paths = Paths(data_in, data_out, params['cell'], ci)
         process(paths, ci, 1, params, True, features)
