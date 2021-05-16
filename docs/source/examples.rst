
Examples
========

Below are simple templeates cases of simulation postprocessing in a typical application.

Before proceeding to the proper analysis, feature-independent settings need to be sspecified.

* External input and output paths for the data. For the examples, data input files can be downloaded from `data/in`_ directory of the GitHub repository.
* Index of the simulation run.
* The cell type. It has three attributes packed inside a class: *typename* and *plmind* are used in constructing cell-specific path to input files and correspond to internal directory structure the Explicit Microtubules. *regions* cell field defines borders of cell internal subcompartments.

.. _`data/in`: https://github.com/vsukhor/cytoskeleton-analyser/

Below, we will process the same simulation, so it is convenient to store
the settings in a separate python module for importing them into the process-specific scripts.

.. code-block:: python
   :caption: sim_specs.py
   :linenos:

   from pathlib import PurePath                     # OS-independedt paths
   from numpy import Inf
   from cytoskeleton_analyser.cells import CellType
   from cytoskeleton_analyser.cells import Region, Regions

   # Set data source and destination. !! CHANGE !!
   data_in = PurePath('/Full/path/to/input/directory')
   data_out = PurePath('Full/path/to/output/directory')

   # Specify the simulation runs to analyse:
   rinds = [1]

   # The runs represent the same cell type.
   # Derive the cell class from CellType to ensure
   # that no attributes are missing.
   cell = CellType(
      # Name of the data foder for this cell type.
      typename = 'irreg',
      # Index of the membrane mesh.
      plmind = 1,
      # Radial distance limits (in μm) for the cell subcompartments.
      regions = Regions(
         cytosol=Region(0., Inf),
         soma=Region(0., 6.),
         lamella=Region(6., Inf),
         edge=Region(10., Inf),
      )
   )

Now let us see the use cases:

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   example1_config
   example2_space
   example3_history

