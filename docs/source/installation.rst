
Installation
=============

Requirements
-------------

Cytoskeleton Analyser requires `python3.9`_ or later.

.. _python3.9: https://www.python.org/downloads/

When installing by `pip`_, the following dependencies are installed *automatically*:

* `NumPy`_
* `SciPy`_
* `Matplotlib`_
* `meshio`_
* `SQLAlchemy`_

.. _NumPy: http://www.numpy.org
.. _SciPy: https://www.scipy.org/
.. _Matplotlib: https://matplotlib.org/
.. _meshio: https://github.com/nschloe/meshio
.. _SQLAlchemy: https://www.sqlalchemy.org/

Installation
------------

Provided that python is present,
in the simplest case the Cytoskeleton Analyser may be installed globallly using `pip`_:

.. _pip: https://pip.pypa.io

.. code-block:: text

    $ pip install cytoskeleton-analyser

However, the recommended way is installing the package into a `virtual environment`_, which adds only a few commans:

1. Create a project directory (here named 'myproject') and cd to it.

2. Create a python virtual environment (here named 'venv3.9').

3. Activate the virtual environment.

4. Install as above.

.. _`virtual environment`: https://virtualenv.pypa.io

.. code-block:: text
   :linenos:

    $ mkdir myproject && cd myproject
    $ python3 -m venv venv3.9
    $ source venv3.9/bin/activate
    $ pip install cytoskeleton-analyser


Alternatively, if you want to install the python package straight from the `GitHub`_ repository, replace the emphasized line above with

.. code-block:: text

    $ pip install git+git://github.com/vsukhor/cytoskeleton-analyser

The Source Code
---------------

Source code of Cytoskeleton Analyser is publicly available at `GitHub`_. You can clone the whole repository with

.. code-block:: text

    $ git clone https://github.com/vsukhor/cytoskeleton-analyser.git

.. _`GitHub`: https://github.com/vsukhor/cytoskeleton-analyser

License
-------

**cytoskeleton-analyser** is available under the terms of the BSD 3 Clause license. See `LICENSE.txt`_ for more information.

.. _`LICENSE.txt`: https://github.com/vsukhor/cytoskeleton-analyser/LICENSE.txt

