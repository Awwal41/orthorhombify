Changelog
=========

2.1.0
-----

* All crystal systems: cell shape is reported as cubic, tetragonal, orthorhombic,
  hexagonal, monoclinic, rhombohedral or triclinic. Hexagonal cells become exactly
  orthorhombic (N = 2), rhombohedral cells too (N = 6)
* New default 3D search for triclinic / rhombohedral cells (``--method vectors``):
  better angles and much faster than the v1 search (still available as ``--method hnf``)
* Lattice directions: every result prints and saves (``.json``) how the new axes and
  planes are made of the old ones (and back), including mixed directions such as
  ``a' = a + c``. Python: ``LatticeMapping``, ``Candidate.mapping``
* New ``orthorhombify-transform`` command and ``transform`` / ``back_to_original``:
  any cell change (fractional matrices too), e.g. orthorhombic to hexagonal, and the
  way back to the original cell from a LAMMPS ``*_equil.data`` file, undoing
  ``replicate`` and averaging the copies of each atom
* ``orthorhombify`` reads LAMMPS ``.data`` files directly

2.0.0
-----

* Monoclinic support: exact 2D search perpendicular to the unique axis (any setting:
  unique axis a, b or c)
* Auto-detects monoclinic / triclinic / orthorhombic cells
* Best-angle-per-box-size table, ``--tol`` for the smallest good-enough box
* ``--force-orthogonal`` with printed strain
* Several input files and wildcards in one run
* LAMMPS writer: data file with per-element starting charges (neutralized by default)
  and a ReaxFF + QEq minimize / NVT / NPT input script
* Installable package (``pip install``) with ``orthorhombify`` and
  ``orthorhombify-lammps`` commands and a Python API
* beta-HMX example bundled (``orthorhombify.example_path()``), pytest suite
* Sphinx documentation

1.0.0
-----

* Original ``orthorhombify.py``: 3D HNF x unimodular search for triclinic cells
