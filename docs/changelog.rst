Changelog
=========

2.0.0
-----

* Monoclinic support: exact 2D search perpendicular to the unique axis (any setting:
  unique axis a, b or c)
* Auto-detects monoclinic / triclinic / orthorhombic cells
* Best-angle-per-box-size table, ``--tol`` for the smallest good-enough box
* ``--force-orthogonal`` with printed strain
* Several input files and wildcards in one run
* New ``make_lammps_inputs.py``: LAMMPS data file with per-element starting charges
  (neutralized by default) and a ReaxFF + QEq minimize / NVT / NPT input script
* Sphinx documentation

1.0.0
-----

* Original ``orthorhombify.py``: 3D HNF x unimodular search for triclinic cells
