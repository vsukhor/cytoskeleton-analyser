#!/usr/bin/env python3
"""Cytoskeleton Analyser setup script.
"""

from setuptools import setup

# Metadata goes in setup.cfg.
# These are here for GitHub's dependency graph.
setup(
    name="cytoskeleton-analyser",
    install_requires=[
        "matplotlib==3.3.4",
        "meshio>=4.3.11",
        "scipy>=1.6.0",
        "numpy>=1.20.0",
        "SQLAlchemy>=1.3.23",
    ],
)