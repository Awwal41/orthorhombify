"""Orthorhombify: near-orthorhombic supercells of monoclinic and triclinic crystals,
and LAMMPS ReaxFF + QEq inputs with starting charges.

Example::

    from ase.io import read
    import orthorhombify as orth

    atoms = read("beta_HMX.cif")
    best = orth.find_supercells(atoms, max_cells=8, top=1)[0]
    supercell = best.build(atoms)
    orth.write_lammps_inputs(supercell, "lammps/beta_HMX", "beta_HMX")
"""
from importlib.resources import files as _files

from .core import (Candidate, best_angle_by_size, build_supercell, detect_system,
                   find_supercells, min_distance)
from .lammps import (DEFAULT_CHARGES, ChargeInfo, LammpsFiles, assign_charges,
                     make_input_script, suggest_replicate, write_lammps_inputs)

__version__ = "2.0.0"

__all__ = [
    "Candidate", "find_supercells", "best_angle_by_size", "build_supercell",
    "detect_system", "min_distance",
    "DEFAULT_CHARGES", "ChargeInfo", "LammpsFiles", "assign_charges",
    "make_input_script", "suggest_replicate", "write_lammps_inputs",
    "example_path", "__version__",
]


def example_path(name: str = "beta_HMX") -> str:
    """Path to a bundled example CIF. Available: ``"beta_HMX"`` (COD 4110146)."""
    path = _files("orthorhombify") / "data" / f"{name}.cif"
    if not path.is_file():
        raise FileNotFoundError(f"No bundled example named {name!r}")
    return str(path)
