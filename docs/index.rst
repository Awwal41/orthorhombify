Orthorhombify documentation
===========================

**Orthorhombify** turns a monoclinic or triclinic crystal into a near-orthorhombic
supercell (all angles close to 90 degrees) without straining it, then writes
LAMMPS files with starting charges for a ReaxFF + QEq equilibration.

It is two small Python scripts:

* ``orthorhombify_v2.py`` - CIF in, near-orthorhombic CIF / POSCAR out
* ``make_lammps_inputs.py`` - CIF / POSCAR in, LAMMPS data + input script out

You can read the same guides as plain Markdown on GitHub under ``docs/``.

.. tip::

   New here? Go straight to :doc:`getting-started`. It takes about five minutes.

Quick links
-----------

* :doc:`getting-started` - install, convert one crystal, make LAMMPS files
* :doc:`orthorhombify` - how the supercell search works and how to tune it
* :doc:`lammps` - charges, neutralizing, and what the LAMMPS script does
* :doc:`examples` - results for beta-HMX, MDNT, TNT and FOX-7
* :doc:`cli` - every command-line option

.. toctree::
   :maxdepth: 2
   :caption: User guides

   getting-started
   orthorhombify
   lammps
   examples
   cli

.. toctree::
   :maxdepth: 1
   :caption: Project

   changelog
