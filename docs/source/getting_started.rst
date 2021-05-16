
Getting Started
===============

Simulation of the microtubule cytoskeleton with **Explicit Microtubules** produces output stored
in text and binary files. The results carry information on spatial structure, composition and
temporal evolution of the reconstructed system.


Organising configuration settings
---------------------------------

It is often desirable to sample the parameter space by generating ensembles of cell reconstructions
differing in one or many parameters of interest. Being able to automatically collect and categorise
the applied configuration settings is often helpful as a starting point in examination
of the microtubule bechavior.

The collected configurations may be stored in a database. Cytoskeleton Analyser uses
*sqlite3* database. Everything necessary for hangling it is included in the package (For more
convenience, additional open source GUI viewers, e.g. `SQLiteStudio`_ are also available for installation).
The database itself is a regulary file on a disc and does not require a dedicated server.

.. _`SQLiteStudio`: https://sqlitestudio.pl/

.. note::
    Use of the database is optional. Run settings are also serialised to a dictionary saved in a JSON file.


Features
--------

Cytoskeleton Analyser can handle both the time-dependent and the geometric features of the
reconstructed system of microtubules.

The software may be either configured to examine features of interest by specifying the specific
names listed below or be set to process them all in batch mode.
Depending on particular feature, Cytoskeleton-Analyser produces output saved on a disc as one or several files:
whenever applicable, these are a JSON reporting the main summary, accompanied by the frequency distribution of the feature
saved in CSV format and its (optianally gzipped) histogram as SVG image.

During the feature examination, the summary is also dumped to a log file.

Temporal
^^^^^^^^

cytoskeleton-analyser examines time-dependent charactersitics of microtubule
polymerisation kinetics using its 'history' subpackage.

One of the most prominent and biologically important aspects of microtubules is their dinamic
instability, which induces ongoing switches between polymerisation and depolymerisation states
referred to as recoveries and catastrophes respectively. Explicit-Microtubules simulator
approximates the instability in real time and generates a history of microtubule states.
Model of dymanic instability applied by Explicit-Microtubules utilizes extended dynamics whereby
microtubule ends adopt stochastically one of three states: shrink, growth and pause.
To represent the bechaviour of microtubules known from experiments in cell-specific
environment (as opposed to in-vitro systems), probabilities of these
states in the generated microtubules are position-dependent This accounts for the
local modulation of the dynamic instability parameters by cell cortex and the plasma membrane.

Histories of all cell microtubules collected over a preset time interval are imported from
Explicit-Microtubules output file 'history_eX_R.dat'.
Letters 'X' and 'R' in the file name are placeholders indicating:

 X :  Microtubule end: *0* for minus-end or *1* for plus-end

 R :  Index of the simulation run: non-negative number

Experiments in living cells show that strucure and behavior of the cytoskeleton in the cell periphery
differ strongly from those in the central regions for many cell types.
Cytoskeleton Analyser accounts for cell compartmentation by examining peripheral and
central cell regions separately, along with the analysis done for the whole-cell.
The above discrimination is reflected in prefixes *soma*, *lamella* and *global* attached to output file names.
Other components of the names is the feature description itself and index of the microtubule end.
Below are the feature keywords.

**Major temporal characteristics:**:

 **duration** : Time of history recording.

 **length** : Microtubule lengths.

 **comets** : Radial distribution of growing plus-ends.

 **frequencies** : Transition frequencies between the dynamic instability states.

 **spatial_maps** : Spatial density maps of in cell xy projection for the transition types and their relations.

 **event_collections** : Elongation, duration and velocity of microtubules in each of growth, shrink and pause states. Reorientation angle of the microtubule end over the states. Similar characteristics for two- and three-state sequencies where the consecutive states immediately follow each other.

Positional
^^^^^^^^^^

cytoskeleton-analyser examines time-independent geometric features of the microtubule array
using specific classes of its 'position' subpackage.

Explicit-Microtubules outputs spatially resolved data as a series of
instantaneous cytoskeleton snapshots (see Note below), and these are used as input to Cytoskeleton-Analyser.
Required files are 'positions_G_K_R' and 'csk_ages_G_K_R' containing 3d positions
and time after polymerization of the nodes composing microtubules respectively. In the above names,
the letters 'G', 'K' and 'R' are palceholders representing:

 G :  Granularity level: *fine* or *coarse*

 K :  Index of the snapshot: non-negative number

 R :  Index of the simulation run: non-negative number

While the imported replicas of the virtual cytoskeleton consists of a full 3d-resolved set of microtubules,
experimental results in the literature are often derived from a flattened view based on optical microscopy images.
To make the recovered characteristics comparable with the empirical knowledge, one may choose to examine,
in addition to the full representation, 2d projections of the reconstructed cytoskeleton filaments
to xy plane as well as narrow slices cut parallel to xy plane, which approximate confocal or Total Internal Reflection Microscopy.
Whole-cell 3d plot of the microtubules superimposed onto the plasma membrane may be also saved as an SVG image.
Below are the class names of geometry-related features, which also serve as the key words used for
configuring the cytoskeleton-analyser operation.

Characteristics relevant for both full-depth and sliced representations:
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

 **Lengths3d** : Lengths of microtubules in 3d space.

 **Lengths2d** : Lengths of filament projections to a xy plane.

 **RadialMass** : Radial distribution of mass in filament projections to xy plane from cell geometrical center to periphery.

 **RadialEnds** : Radial distribution of plus-ends in filament projections to xy plane from cell geometrical center to periphery.

 **Curvature3d** : Curvature of microtubules in 3d space.

 **AnglesToRad** : Angles between filament xy-projections and radial direction outwards in the cell marginal zone.

 **SegmentNumbers** : Apparent total number of filaments in the system.

Characteristics relevant for the sliced representation only:
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

 **Curvature2dConv** : Apparent curvature of filament xy-projections.

 **Curvature2dMboc17** : Apparent xy-curvature using a simplified approach developed for expeimental analysis of optical microscopy results [1]_.


Characteristics relevant for the full-depth representation only:
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

 **AgesByNode** : Age distribution of all filament nodes (including internal ones).

 **AgesByFilament** : Age distribution of filament '+'-end nodes.


.. note::

    Stochastic models like that implemented within Explicit-Microtubules represent fluctuating environments.
    As a rule, one either collects several snapshots of the system separated by time an interval large enough for
    for temporal autocorrelations in the parameters of interest to decay or performs several independently seeded runs.

    For the specific case of reconstructed microtubules, also single snapshot may be sufficient for the particular case
    when the simulation was configured with inter-microtubule interactions switched off. Then,
    the cytoskeleton filaments grow independent of each other and their whole-cell ensemble samples the geometric
    characteristics in a satisfactory way.


.. [1] `Zhang Zh., Nishimura Y., and Kanchanawonga P. (2017) Extracting microtubule networks from superresolution
    single-molecule localization microscopy data MBoC, 28:2. <https://www.molbiolcell.org/doi/full/10.1091/mbc.e16-06-0421>`_

