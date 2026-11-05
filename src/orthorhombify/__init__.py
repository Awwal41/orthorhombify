"""Orthorhombify: near-orthorhombic supercells of monoclinic, hexagonal, rhombohedral
and triclinic crystals, the lattice directions that connect the old and new cells,
the way back, and LAMMPS ReaxFF + QEq inputs with starting charges.

Example::

    from ase.io import read
    import orthorhombify as orth

    atoms = read("beta_HMX.cif")
    best = orth.find_supercells(atoms, max_cells=8, top=1)[0]
    print(best.mapping.summary(atoms.cell))      # which old directions the new axes follow
    supercell = best.build(atoms)
    orth.write_lammps_inputs(supercell, "lammps/beta_HMX", "beta_HMX")
"""
from importlib.resources import files as _files

__version__ = "2.1.0"

from .core import (Candidate, best_angle_by_size, build_supercell, detect_system,  # noqa: E402
                   find_supercells, min_distance)
from .lammps import (DEFAULT_CHARGES, ChargeInfo, LammpsFiles, assign_charges,  # noqa: E402
                     make_input_script, suggest_replicate, write_lammps_inputs)
from .mapping import LatticeMapping, write_mapping_json  # noqa: E402
from .transform import back_to_original, guess_replicate, parse_matrix, transform  # noqa: E402

__all__ = [
    "Candidate", "find_supercells", "best_angle_by_size", "build_supercell",
    "detect_system", "min_distance",
    "LatticeMapping", "write_mapping_json",
    "transform", "back_to_original", "parse_matrix", "guess_replicate",
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
