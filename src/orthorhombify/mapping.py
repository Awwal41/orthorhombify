"""How the axes, directions and planes of two cells of the same lattice map onto each
other, for a transformation ``new_cell_rows = P @ old_cell_rows``.

Directions [uvw] are coordinates along the cell vectors; planes (hkl) are Miller
indices. With rows convention:

* direction, old -> new: ``u' = u @ P^-1``      new -> old: ``u = u' @ P``
* plane,     old -> new: ``h' = P @ h``         new -> old: ``h = P^-1 @ h'``

so the zone law ``h . u`` is the same in both cells.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from functools import reduce

import numpy as np

OLD_NAMES = ("a", "b", "c")
NEW_NAMES = ("a'", "b'", "c'")


# ---------- exact rational 3x3 algebra ----------
def to_fraction_matrix(M) -> list[list[Fraction]]:
    """3x3 matrix of :class:`~fractions.Fraction` from ints, floats or strings like ``"1/2"``."""
    rows = [[_to_fraction(x) for x in row] for row in (M.tolist() if hasattr(M, "tolist") else M)]
    if len(rows) != 3 or any(len(r) != 3 for r in rows):
        raise ValueError("Transformation matrix must be 3x3")
    return rows

def _to_fraction(x) -> Fraction:
    if isinstance(x, Fraction):
        return x
    if isinstance(x, str):
        return Fraction(x.strip())
    return Fraction(x).limit_denominator(10000)

def _det(M):
    return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
            - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
            + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))

def _inv(M):
    d = _det(M)
    if d == 0:
        raise ValueError("Transformation matrix is singular")
    cof = [[(M[(i + 1) % 3][(j + 1) % 3] * M[(i + 2) % 3][(j + 2) % 3]
             - M[(i + 1) % 3][(j + 2) % 3] * M[(i + 2) % 3][(j + 1) % 3]) for j in range(3)]
           for i in range(3)]
    return [[cof[j][i] / d for j in range(3)] for i in range(3)]  # adjugate / det

def _row_times(v, M):
    return [sum(v[i] * M[i][j] for i in range(3)) for j in range(3)]

def _times_col(M, v):
    return [sum(M[i][j] * v[j] for j in range(3)) for i in range(3)]

def smallest_integers(v) -> list[int]:
    """Scale a rational vector to the smallest whole numbers with the same direction."""
    v = [_to_fraction(x) for x in v]
    den = reduce(lambda x, y: x * y // math.gcd(x, y), (x.denominator for x in v), 1)
    ints = [int(x * den) for x in v]
    g = reduce(math.gcd, (abs(x) for x in ints)) or 1
    return [x // g for x in ints]


# ---------- formatting ----------
def _fmt_frac(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"

def expression(coeffs, names) -> str:
    """``[1, 0, 1]`` with names a, b, c -> ``"a + c"``; fractions allowed."""
    parts = []
    for x, n in zip(coeffs, names):
        x = _to_fraction(x)
        if x == 0:
            continue
        mag = abs(x)
        term = n if mag == 1 else (f"{_fmt_frac(mag)}{n}" if mag.denominator == 1 else f"({_fmt_frac(mag)}){n}")
        parts.append(("- " if x < 0 else "+ ") + term)
    if not parts:
        return "0"
    s = " ".join(parts)
    return s[2:] if s.startswith("+ ") else "-" + s[2:]

def bracket(v, kind="[]") -> str:
    """``[1, 0, -1]`` -> ``"[1 0 -1]"`` (or ``"(1 0 -1)"`` with kind ``"()"``)."""
    return kind[0] + " ".join(str(int(x)) for x in v) + kind[1]


class LatticeMapping:
    """Axes, directions and planes of an old cell expressed in a new cell, and back.

    Args:
        P: 3x3 transformation, ``new_cell_rows = P @ old_cell_rows``. Integers,
            floats or fraction strings (``"1/2"``) are all fine.

    Example::

        m = LatticeMapping([[1, 0, 1], [0, 1, 0], [-4, 0, 3]])
        m.new_axis_expressions()      # ["a + c", "b", "-4a + 3c"]
        m.direction_to_new([1, 1, 0]) # old [110] in new-cell coordinates
        m.plane_to_old([0, 0, 1])     # Miller indices of the new (001) plane in the old cell
    """

    def __init__(self, P):
        self.P = to_fraction_matrix(P)
        self.P_inv = _inv(self.P)
        self.det = _det(self.P)

    # --- general conversions ---
    def direction_to_new(self, uvw) -> list[Fraction]:
        """Old-cell direction [uvw] in new-cell coordinates."""
        return _row_times([_to_fraction(x) for x in uvw], self.P_inv)

    def direction_to_old(self, uvw) -> list[Fraction]:
        """New-cell direction [uvw] in old-cell coordinates."""
        return _row_times([_to_fraction(x) for x in uvw], self.P)

    def plane_to_new(self, hkl) -> list[int]:
        """Old-cell plane (hkl) as new-cell Miller indices (smallest integers)."""
        return smallest_integers(_times_col(self.P, [_to_fraction(x) for x in hkl]))

    def plane_to_old(self, hkl) -> list[int]:
        """New-cell plane (hkl) as old-cell Miller indices (smallest integers)."""
        return smallest_integers(_times_col(self.P_inv, [_to_fraction(x) for x in hkl]))

    def inverse(self) -> "LatticeMapping":
        """The mapping from the new cell back to the old one."""
        return LatticeMapping(self.P_inv)

    # --- the three axes ---
    def new_axis_expressions(self) -> list[str]:
        """Each new axis as a combination of old axes, e.g. ``"a + c"``."""
        return [expression(self.P[i], OLD_NAMES) for i in range(3)]

    def old_axis_expressions(self) -> list[str]:
        """Each old axis as a combination of new axes, e.g. ``"(3/7)a' - (1/7)c'"``."""
        return [expression(self.P_inv[i], NEW_NAMES) for i in range(3)]

    def new_axes_as_old_directions(self) -> list[list[int]]:
        """Direction [uvw] (old cell) that each new axis points along."""
        return [smallest_integers(self.P[i]) for i in range(3)]

    def old_axes_as_new_directions(self) -> list[list[int]]:
        """Direction [uvw] (new cell) that each old axis points along."""
        return [smallest_integers(self.P_inv[i]) for i in range(3)]

    def axis_angles(self, old_cell) -> list[float]:
        """Angle (degrees) between each new axis and the old axis of the same name."""
        A = np.asarray(old_cell, dtype=float)
        B = np.array([[float(x) for x in row] for row in self.P]) @ A
        out = []
        for i in range(3):
            c = np.dot(A[i], B[i]) / (np.linalg.norm(A[i]) * np.linalg.norm(B[i]))
            out.append(math.degrees(math.acos(max(-1.0, min(1.0, c)))))
        return out

    # --- output ---
    def summary(self, old_cell=None, indent: str = "  ") -> str:
        """Human-readable table of axes and planes (old = a, b, c; new = a', b', c')."""
        new_e, old_e = self.new_axis_expressions(), self.old_axis_expressions()
        new_d, old_d = self.new_axes_as_old_directions(), self.old_axes_as_new_directions()
        ang = self.axis_angles(old_cell) if old_cell is not None else None
        w1 = max(len(s) for s in new_e + old_e) + 2
        lines = [f"{indent}Lattice directions (old cell a, b, c -> new cell a', b', c'):"]
        for i in range(3):
            extra = f"   {ang[i]:6.2f} deg from old {OLD_NAMES[i]}" if ang else ""
            lines.append(f"{indent}  new {NEW_NAMES[i]:<2} = {new_e[i]:<{w1}} along old {bracket(new_d[i]):<12}{extra}")
        for i in range(3):
            lines.append(f"{indent}  old {OLD_NAMES[i]:<2} = {old_e[i]:<{w1}} along new {bracket(old_d[i])}")
        lines.append(f"{indent}Planes:")
        for i, e in enumerate(np.eye(3, dtype=int)):
            lines.append(f"{indent}  new {bracket(e, '()')} = old {bracket(self.plane_to_old(e), '()'):<12}"
                         f"   old {bracket(e, '()')} = new {bracket(self.plane_to_new(e), '()')}")
        return "\n".join(lines)

    def to_dict(self, old_cell=None) -> dict:
        """JSON-friendly description (fractions written as strings)."""
        eye = np.eye(3, dtype=int)
        d = {
            "P": [[_fmt_frac(x) for x in row] for row in self.P],
            "P_inverse": [[_fmt_frac(x) for x in row] for row in self.P_inv],
            "det_P": _fmt_frac(self.det),
            "new_axes": {NEW_NAMES[i]: {"expression": self.new_axis_expressions()[i],
                                        "old_direction": self.new_axes_as_old_directions()[i]}
                         for i in range(3)},
            "old_axes": {OLD_NAMES[i]: {"expression": self.old_axis_expressions()[i],
                                        "new_direction": self.old_axes_as_new_directions()[i]}
                         for i in range(3)},
            "new_planes_in_old": {bracket(e, "()"): self.plane_to_old(e) for e in eye},
            "old_planes_in_new": {bracket(e, "()"): self.plane_to_new(e) for e in eye},
        }
        if old_cell is not None:
            d["axis_angles_deg"] = dict(zip(NEW_NAMES, [round(x, 4) for x in self.axis_angles(old_cell)]))
        return d


def write_mapping_json(path: str, mapping: LatticeMapping, old_atoms, new_atoms, **extra) -> None:
    """Save a mapping plus both cells, so the transformation can be undone later
    (``orthorhombify-transform --inverse``)."""
    from . import __version__
    data = {
        "orthorhombify_version": __version__,
        **extra,
        "old_cell": _cell_dict(old_atoms),
        "new_cell": _cell_dict(new_atoms),
        **mapping.to_dict(old_atoms.cell.array),
    }
    with open(path, "w") as fd:
        json.dump(data, fd, indent=2)

def _cell_dict(atoms):
    a, b, c, al, be, ga = (round(float(x), 6) for x in atoms.cell.cellpar())
    return {"a": a, "b": b, "c": c, "alpha": al, "beta": be, "gamma": ga,
            "n_atoms": len(atoms), "formula": atoms.get_chemical_formula()}
