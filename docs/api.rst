API reference
=============

Everything below is importable from the top-level package::

   import orthorhombify as orth

Supercells
----------

.. autofunction:: orthorhombify.find_supercells

.. autoclass:: orthorhombify.Candidate
   :members: build

.. autofunction:: orthorhombify.best_angle_by_size

.. autofunction:: orthorhombify.detect_system

.. autofunction:: orthorhombify.build_supercell

.. autofunction:: orthorhombify.min_distance

LAMMPS and charges
------------------

.. py:data:: orthorhombify.DEFAULT_CHARGES
   :value: {"O": -0.3668, "N": 0.3668, "C": -0.4320, "H": 0.2160}

   Starting charge (e) of each element, used when no ``charges`` are given.

.. autofunction:: orthorhombify.write_lammps_inputs

.. autofunction:: orthorhombify.assign_charges

.. autofunction:: orthorhombify.suggest_replicate

.. autofunction:: orthorhombify.make_input_script

.. autoclass:: orthorhombify.ChargeInfo

.. autoclass:: orthorhombify.LammpsFiles

Examples
--------

.. autofunction:: orthorhombify.example_path
