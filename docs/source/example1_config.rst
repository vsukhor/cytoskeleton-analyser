
Example 1: Organizing Simulation Parameters
-------------------------------------------

Two versions of configuration processing are shown. The first one merely collects the simulation parameters
into a dictionary. The Second version stores them additionally in a database.

Without the use of a database:
""""""""""""""""""""""""""""""

.. code-block:: python
   :linenos:

   from cytoskeleton_analyser.inout import Paths     # simulation-specific paths
   from cytoskeleton_analyser.config.drivers import process  # driver fuction
   from sim_specs import cell, data_in, data_out, rinds

   # Decide if this is a new analysis, or we merely want to import one.
   new = True

   if __name__ == '__main__':

      for ci in rinds:
         paths = Paths(data_in, data_out, cell, ci)
         process(new, paths, ci, cell)


Accumulate the settings in a database:
"""""""""""""""""""""""""""""""""""""""
Similar to the above, but the database module is included (line 3).
The data are collected in 'cfs' (lines 11,12) list before sending to the database (line 13).

.. code-block:: python
   :linenos:

   from cytoskeleton_analyser.inout import Paths     # simulation-specific paths
   from cytoskeleton_analyser.config.drivers import process  # driver fuction
   import cytoskeleton_analyser.database.sqlite_alchemy_orm.db as db  # database
   from sim_specs import cell, data_in, data_out, rinds

   # Decide if this is a new analysis, or we merely want to import one.
   new = True

   if __name__ == '__main__':

      cfs = [process(new, Paths(data_in, data_out, cell, ci), ci, cell)
             for ci in rinds]
      db.update(cfs, cell, Paths(data_in, data_out, cell).data_out.parent)


