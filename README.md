<a href="https://cytoskeleton-analyser.readthedocs.io">
    <img src="https://raw.githubusercontent.com/vsukhor/cytoskeleton-analyser/master/docs/source/_static/artwork/logo_t_large.svg"
         alt="cytoskeleton-analyser logo"
         align="top">
</a>

[![pypi](https://img.shields.io/pypi/v/cytoskeleton-analyser.svg)](https://pypi.org/project/cytoskeleton-analyser/)
![python](https://img.shields.io/badge/python-3.9-blue.svg)
[![Documentation Status](https://readthedocs.org/projects/cytoskeleton-analyser/badge/?version=latest)](https://cytoskeleton-analyser.readthedocs.io/en/latest/?badge=latest)
[![License: BSD3](https://img.shields.io/badge/License-BSD3-blueviolet.svg)](./LICENSE.md)

**cytoskeleton-analyser** is a postprocessor for extraction, analysis and human-friendly presentation 
of the raw binary results generated by Explicit-Microtubules.

**Explicit-Microtubules** Explicit-Microtubules is a C++ simulator of the whole-cell microtubule system. 
It reconstructs computationally the microtubule geometry and real-time dynamics in 3d space and 
time using stochastic algorithms. In the course of the simulation based on a multitude of bioogical 
information, it outputs detailed histories of microtubule dynamics along with three-dimensional 
snapshots of microtubules. However, designed to be focuced on versatile generation procedures, 
Explicit-Microtubules only includes basic data analysis.

Complementing the simulator, **cytoskeleton-analyser** is designed to extract, analyse 
and visualize various characteristics of the stochastically reconstructed microtubule system.


## Requirements

Cytoskeleton Analyser requires [python 3.9](https://www.python.org/downloads/) or later. 

When installing by [pip](https://pip.pypa.io), the following dependencies are installed *automatically*:

* [NumPy](http://www.numpy.org)
* [SciPy](https://www.scipy.org/)
* [Matplotlib](https://matplotlib.org/)
* [meshio](https://github.com/nschloe/meshio)
* [SQLAlchemy](https://www.sqlalchemy.org/)


## Installation


Provided that python is present, in the simplest case the Cytoskeleton Analyser 
may be installed globallly using [pip](https://pip.pypa.io):

```bash
pip install cytoskeleton-analyser
```

However, the recommended way is installing the package into a
[python virtual environment](https://virtualenv.pypa.io), which adds only a few commans:

1. Create a project directory (here named ???myproject???) and cd to it.

2. Create a python virtual environment (here named ???venv3.9???).

3. Activate the virtual environment.

4. Install as above.

```bash
mkdir myproject && cd myproject
python3 -m venv venv3.9
source venv3.9/bin/activate
pip install cytoskeleton-analyser
```

Alternatively, if you want to install the python package straight from the GitHub repository, 
replace the last line above with

```bash
pip install git+git://github.com/vsukhor/cytoskeleton-analyser
```

## Getting Started

A more complete introduction, use examples and API are available in the [documentation](https://cytoskeleton-analyser.readthedocs.io/en/latest/index.html).`

## License
**cytoskeleton-analyser** is available under the terms of the BSD 3 Clause license. 
See LICENSE.txt for more information.

