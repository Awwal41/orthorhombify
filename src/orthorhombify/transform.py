"""Apply any cell transformation, including fractional ones that make the cell
smaller, and undo an earlier orthorhombify step (for example after an MD run)."""
from __future__ import annotations

import json
from fractions import Fraction

import numpy as np
from ase import Atoms
from ase.build import make_supercell
from ase.cell import Cell
from ase.neighborlist import neighbor_list

from .mapping import LatticeMapping, to_fraction_matrix


def parse_matrix(text: str) -> list[list[Fraction]]:
    """Parse a matrix written as text.

    Rows are separated by ``;``, numbers by spaces or commas, and fractions are
    allowed. Three numbers mean a diagonal matrix.

    Examples: ``"1 0 1; 0 1 0; -4 0 3"``, ``"1/2 -1/2 0; 1/2 1/2 0; 0 0 1"``, ``"2 2 1"``.
    """
    rows = [r.replace(",", " ").split() for r in text.strip().split(";") if r.strip()]
    if len(rows) == 1 and len(rows[0]) == 3:
        d = rows[0]
        rows = [[d[0], "0", "0"], ["0", d[1], "0"], ["0", "0", d[2]]]
    elif len(rows) == 1 and len(rows[0]) == 9:
        rows = [rows[0][0:3], rows[0][3:6], rows[0][6:9]]
    return to_fraction_matrix(rows)


def _standardize(atoms, M):
    atoms.set_cell(Cell.fromcellpar(atoms.cell.cellpar()), scale_atoms=True)
    atoms.wrap()
    atoms.info["orthorhombify_P"] = [[str(x) for x in row] for row in M]
    return atoms


def transform(atoms: Atoms, matrix, *, reduce: str = "auto", tol: float = 0.5) -> Atoms:
    """Build the cell ``matrix @ cell`` of ``atoms`` and fill it with the right atoms.

    Args:
        atoms: The structure.
        matrix: 3x3 transformation (ints, floats, fraction strings), or text in the
            command-line form such as ``"1 0 0; -1/2 1/2 0; 0 0 1"`` or ``"2 2 1"``.
            ``new_cell_rows = matrix @ old_cell_rows``.
        reduce: Only used when the new cell is *smaller* (fractional matrix).
            ``"average"`` averages every copy of each atom, which is the right choice
            for a supercell built by Orthorhombify (copies are stored one after another),
            including after an MD run. ``"dedup"`` folds all atoms into the new cell and
            merges atoms of the same element closer than ``tol``. ``"auto"`` uses
            average if the atom order allows it, otherwise dedup.
        tol: Distance (Å) under which folded atoms count as the same atom (dedup), and
            the largest allowed spread between copies of one atom (average).

    Returns:
        New ``Atoms`` in the standard orientation (a along x, b in the xy plane).
        ``info["orthorhombify_reduce"]`` records which method was used and
        ``info["orthorhombify_P"]`` the matrix (as fraction strings).
    """
    M = parse_matrix(matrix) if isinstance(matrix, str) else to_fraction_matrix(matrix)
    mapping = LatticeMapping(M)
    if mapping.det < 0:
        raise ValueError("The transformation flips handedness (det < 0); negate one row.")
    Mf = np.array([[float(x) for x in row] for row in M])
    n_new = len(atoms) * mapping.det

    if all(x.denominator == 1 for row in M for x in row):
        out = make_supercell(atoms, np.array([[int(x) for x in row] for row in M]))
        out.info["orthorhombify_reduce"] = "supercell"
        return _standardize(out, M)

    if n_new.denominator != 1:
        raise ValueError(f"{len(atoms)} atoms x det {mapping.det} is not a whole number of atoms; "
                         "this matrix does not describe a cell of this lattice.")
    n_new = int(n_new)
    new_cell = Mf @ atoms.cell.array

    if reduce in ("auto", "average"):
        out = _reduce_average(atoms, new_cell, n_new, tol)
        if out is not None:
            out.info["orthorhombify_reduce"] = "average"
            return _standardize(out, M)
        if reduce == "average":
            raise ValueError("Cannot average: atoms are not stored as consecutive identical copies "
                             "(use reduce='dedup').")
    out = _reduce_dedup(atoms, new_cell, tol)
    out.info["orthorhombify_reduce"] = "dedup"
    if len(out) != n_new:
        out.info["orthorhombify_warning"] = (f"expected {n_new} atoms, got {len(out)}; "
                                             "try a different tol")
    return _standardize(out, M)


def _reduce_average(atoms, new_cell, n_new, tol):
    """Average the copies of each atom; copy j of atom i is atom j*n_new + i."""
    if n_new == 0 or len(atoms) % n_new:
        return None
    n_copies = len(atoms) // n_new
    symbols = np.array(atoms.get_chemical_symbols()).reshape(n_copies, n_new)
    if not (symbols == symbols[0]).all():
        return None
    frac = (atoms.positions @ np.linalg.inv(new_cell)).reshape(n_copies, n_new, 3)
    shift = frac - frac[0]
    shift -= np.round(shift)              # copies differ by whole new-cell vectors
    spread = np.linalg.norm(shift @ new_cell, axis=2).max()
    if spread > tol:
        return None
    mean = frac[0] + shift.mean(axis=0)
    out = Atoms(symbols=list(symbols[0]), scaled_positions=mean % 1.0, cell=new_cell, pbc=True)
    if atoms.has("initial_charges"):
        out.set_initial_charges(atoms.get_initial_charges().reshape(n_copies, n_new).mean(axis=0))
    out.info["orthorhombify_spread"] = float(spread)
    return out


def _reduce_dedup(atoms, new_cell, tol):
    folded = Atoms(symbols=atoms.get_chemical_symbols(), cell=new_cell, pbc=True,
                   scaled_positions=(atoms.positions @ np.linalg.inv(new_cell)) % 1.0)
    i, j = neighbor_list("ij", folded, tol)
    sym = np.array(folded.get_chemical_symbols())
    same = (sym[i] == sym[j]) & (i < j)
    drop = np.zeros(len(folded), dtype=bool)
    for a, b in sorted(zip(i[same], j[same])):
        if not drop[a]:
            drop[b] = True
    return folded[~drop]


def back_to_original(atoms: Atoms, mapping_json: str, *, replicate=None, reduce: str = "auto",
                     tol: float = 0.5) -> Atoms:
    """Undo an Orthorhombify step: go from the (equilibrated) supercell back to the
    original cell, using the ``*_orth_N.json`` file written next to it.

    Args:
        atoms: The supercell, for example read from the LAMMPS ``*_equil.data`` file.
        mapping_json: The ``.json`` file written by ``orthorhombify``.
        replicate: ``(nx, ny, nz)`` if the box was replicated in LAMMPS. Worked out
            from the atom count and box lengths if not given.
        reduce, tol: See :func:`transform`.

    Returns:
        The structure in the original cell. ``info["orthorhombify_replicate"]`` records
        the replication that was undone.
    """
    with open(mapping_json) as fd:
        data = json.load(fd)
    P_inv = to_fraction_matrix(data["P_inverse"])
    n_sc = data["new_cell"]["n_atoms"]
    if replicate is None:
        replicate = guess_replicate(atoms, data)
    reps = [int(x) for x in replicate]
    if len(atoms) != n_sc * reps[0] * reps[1] * reps[2]:
        raise ValueError(f"{len(atoms)} atoms does not match {n_sc} x replicate {tuple(reps)}")
    undo_rep = [[Fraction(1, reps[0]), 0, 0], [0, Fraction(1, reps[1]), 0], [0, 0, Fraction(1, reps[2])]]
    M = [[sum(P_inv[i][k] * undo_rep[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    out = transform(atoms, M, reduce=reduce, tol=tol)
    out.info["orthorhombify_replicate"] = tuple(reps)
    return out


def guess_replicate(atoms: Atoms, data: dict):
    """Replication of the saved supercell that best matches ``atoms`` (box lengths)."""
    n_sc = data["new_cell"]["n_atoms"]
    if len(atoms) % n_sc:
        raise ValueError(f"{len(atoms)} atoms is not a multiple of the saved supercell ({n_sc} atoms)")
    ref = np.array([data["new_cell"][k] for k in ("a", "b", "c")])
    reps = np.maximum(1, np.round(atoms.cell.lengths() / ref)).astype(int)
    if reps.prod() != len(atoms) // n_sc:
        raise ValueError("Could not work out the replication from the box size; pass replicate=(nx, ny, nz)")
    return tuple(int(x) for x in reps)
