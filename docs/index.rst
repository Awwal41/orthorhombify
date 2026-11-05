Orthorhombify documentation
===========================

**Orthorhombify** turns a monoclinic or triclinic crystal into a near-orthorhombic
supercell (all angles close to 90 degrees) without straining it, then writes
LAMMPS files with starting charges for a ReaxFF + QEq equilibration.

.. code-block:: bash

   pip install git+https://github.com/Awwal41/orthorhombify.git

This gives you two commands and a Python package:

* ``orthorhombify`` - CIF in, near-orthorhombic CIF / POSCAR out
* ``orthorhombify-lammps`` - CIF / POSCAR in, LAMMPS data + input script out
* ``import orthorhombify`` - the same tools as Python functions

You can read the same guides as plain Markdown on GitHub under ``docs/``.

.. tip::

   New here? Go straight to :doc:`getting-started`. It takes about five minutes.

Quick links
-----------

* :doc:`getting-started` - install, convert one crystal, make LAMMPS files
* :doc:`install` - pip, venv, conda, HPC, from source
* :doc:`orthorhombify` - how the supercell search works and how to tune it
* :doc:`lammps` - charges, neutralizing, and what the LAMMPS script does
* :doc:`python-api` - use it from your own Python scripts
* :doc:`examples` - results for beta-HMX, MDNT, TNT and FOX-7
* :doc:`cli` - every command-line option
* :doc:`api` - function reference

.. toctree::
   :maxdepth: 2
   :caption: User guides

   getting-started
   install
   orthorhombify
   lammps
   python-api
   examples
   cli

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api

.. toctree::
   :maxdepth: 1
   :caption: Project

   changelog
