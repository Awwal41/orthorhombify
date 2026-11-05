Orthorhombify documentation
===========================

**Orthorhombify** turns any crystal (monoclinic, hexagonal, rhombohedral or triclinic)
into a near-orthorhombic supercell (all angles close to 90 degrees) without straining
it, tells you which old lattice directions the new axes follow, writes LAMMPS files with
starting charges for a ReaxFF + QEq equilibration, and takes the result back to the
original cell afterwards.

.. code-block:: bash

   pip install git+https://github.com/Awwal41/orthorhombify.git

This gives you three commands and a Python package:

* ``orthorhombify`` - CIF in, near-orthorhombic CIF / POSCAR / JSON out
* ``orthorhombify-lammps`` - CIF / POSCAR in, LAMMPS data + input script out
* ``orthorhombify-transform`` - any other cell, or back to the original (e.g. after LAMMPS)
* ``import orthorhombify`` - the same tools as Python functions

You can read the same guides as plain Markdown on GitHub under ``docs/``.

.. tip::

   New here? Go straight to :doc:`getting-started`. It takes about five minutes.

Quick links
-----------

* :doc:`getting-started` - install, convert one crystal, make LAMMPS files
* :doc:`install` - pip, venv, conda, HPC, from source
* :doc:`orthorhombify` - how the supercell search works, for every crystal system
* :doc:`lattice-directions` - which old directions and planes the new axes follow
* :doc:`lammps` - charges, neutralizing, and what the LAMMPS script does
* :doc:`transform` - back to the original cell after MD, or to any other cell
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
   lattice-directions
   lammps
   transform
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
